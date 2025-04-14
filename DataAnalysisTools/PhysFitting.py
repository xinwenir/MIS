import array
import os
from datetime import datetime

import numpy as np

from DataFactory import CommonUtils, ROOTUtils
from DataFactory import PhysModel
from DataAnalysisTools import dic_tree_branch_name
from ROOT import TMath, TFile, TGraphErrors, TTree, TF1, TCanvas
from ROOT import gROOT

import logging

logger = logging.getLogger(__name__)


class PhysFitting:
    def __init__(self):
        self._tree = None
        self._par_fitting = None
        self._gr_err = None
        self._phys_model = None

    def _read_tree_from_file(self, fin_name, tree_name='t_tree'):
        if ROOTUtils.find_key_in_file(fin_name, tree_name) == '':
            raise Exception('Cannot find %s in file %s.' % (tree_name, fin_name))

        fin = TFile(fin_name, 'read')

        t_tree = fin.Get(ROOTUtils.find_key_in_file(fin_name, tree_name))

        t_tree.SetDirectory(gROOT)

        fin.Close()

        return t_tree

    def fit_parameters(self, fin_name, phys_model, tree_name,
                       eps=0.0001, fout_name='', rob=1., lis_fix_number=None):
        self._phys_model = phys_model

        if ROOTUtils.find_key_in_file(fin_name, tree_name) == '':
            raise Exception('Cannot find %s in file %s.' % (tree_name, fin_name))

        fin = TFile(fin_name, 'read')

        self._tree = fin.Get(ROOTUtils.find_key_in_file(fin_name, tree_name))

        n_entries = self._tree.GetEntries()
        print('fit_parameters: Get %d entries' % n_entries)

        df_tr = ROOTUtils.cvt_tree_to_dataframe(self._tree)

        a_theta = array.array('d', df_tr[dic_tree_branch_name['theta']].to_numpy())
        a_sim_Leff = array.array('d', df_tr[dic_tree_branch_name['sim_leff']].to_numpy())
        a_meas_ratio = array.array('d', df_tr[dic_tree_branch_name['meas_ratio']].to_numpy())
        a_meas_ratio_err = array.array('d', df_tr[dic_tree_branch_name['meas_ratio_err']].to_numpy())

        a_emin = array.array('d', df_tr.shape[0]*[0.])
        a_emin_err = array.array('d', df_tr.shape[0]*[0.])

        for i in range(df_tr.shape[0]):
            print('\r %d / %d = %.2f    ' % (i, n_entries, 1. * i / n_entries), end='')
            a_emin[i] = phys_model.get_Emin_from_ratio(theta=a_theta[i], ratio=a_meas_ratio[i])
            de_dl = phys_model.get_demin_dratio(a_theta[i], a_meas_ratio[i], eps=eps, emin=a_emin[i])
            a_emin_err[i] = de_dl * a_meas_ratio_err[i]
        print('')

        self._gr_err = TGraphErrors(n_entries, a_sim_Leff, a_emin, 0, a_emin_err)

        self._gr_err.SetMinimum(0)

        e_max = np.amax(a_emin)

        if phys_model.get_relation_id() == 1:
            f_fit = TF1('f_fit', 'pol1', 0, e_max)
            f_fit.SetParameters(self._phys_model.get_parameters()[0], self._phys_model.get_parameters()[1])
            print('Set parameters: ', f_fit.GetParameters()[0], f_fit.GetParameters()[1])
        elif phys_model.get_relation_id() == 2:
            f_fit = TF1('f_fit', 'pol2', 0, e_max)
            f_fit.SetParameters(self._phys_model.get_parameters()[0],
                                self._phys_model.get_parameters()[1],
                                self._phys_model.get_parameters()[2])
            print('Set parameters: ', f_fit.GetParameters()[0], f_fit.GetParameters()[1], f_fit.GetParameters()[2])
        else:
            raise Exception('PhysicalModel RalationID must be 1 or 2.')

        if lis_fix_number is not None:
            for i_fix in lis_fix_number:
                f_fit.FixParameter(i_fix, self._phys_model.get_parameters()[i_fix])
                print('Fix para %d as' % i_fix, f_fit.GetParameters()[i_fix])

        opt = ''
        if lis_fix_number is not None:
            opt = opt + 'B'
        if rob < 1:
            opt = opt + "+rob=%f" % rob

        self._gr_err.Fit(f_fit, opt)

        self._gr_err.GetXaxis().SetTitle('L_{eff} [m g/cm^{3}]')
        self._gr_err.GetYaxis().SetTitle('E_{min} [Gev]')

        t_par = f_fit.GetParameters()
        if phys_model.get_relation_id() == 1:
            self._par_fitting = [t_par[0], t_par[1]]
        elif phys_model.get_relation_id() == 2:
            self._par_fitting = [t_par[0], t_par[1], t_par[2]]

        root_fname = os.path.splitext(fout_name)[0] + '.root'
        fout = TFile(root_fname, 'recreate')
        self._gr_err.Write()
        can = TCanvas('can_phys_model_fitting', '', 2000, 1500)
        self._gr_err.Draw('ap')
        can.Write()
        fout.Save()
        fout.Close()

        # fill emin, emin_err and fit_leff, and print as csv
        df_tr['emin'] = a_emin
        df_tr['emin_err'] = a_emin_err
        a_fit_leff = array.array('d', df_tr.shape[0] * [0.])
        for i in range(df_tr.shape[0]):
            a_fit_leff[i] = phys_model.get_Leff_from_Emin(a_emin[i], self._par_fitting)

        df_tr['fit_leff'] = a_fit_leff
        csv_fname = os.path.splitext(fout_name)[0] + '.csv'
        df_tr.to_csv(csv_fname, index=False, sep=' ')

    def write_parameters(self, fout_name, phys_model_name='PhysTmp'):
        print('\n\n', datetime.now(), file=open(fout_name, 'a'))
        print('[%s]' % phys_model_name, file=open(fout_name, 'a'))
        print(self._phys_model.get_flux_info(), file=open(fout_name, 'a'))
        print('Relation = %i' % self._phys_model.get_relation_id(), file=open(fout_name, 'a'))
        print('Parameters = %s' % ','.join([str(i) for i in self._par_fitting]), file=open(fout_name, 'a'))
        print('Ecut = %f' % self._phys_model.get_ecut(), file=open(fout_name, 'a'))
