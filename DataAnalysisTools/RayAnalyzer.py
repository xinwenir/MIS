import os, sys, array
import time
from datetime import datetime

import pandas as pd
import numpy as np
from ROOT import TH2F, TTree, TFile, TVector3, TMath, TList, TH1F, TRandom
from ROOT import TGraph2DErrors, TGraphErrors, TMultiGraph, TCanvas, TLegend

from DataFactory import CommonUtils, ROOTUtils
from DataAnalysisTools import _hthetaxy_name_prefix, _hthetaphi_name_prefix, dic_tree_branch_name

import logging

logger = logging.getLogger(__name__)


class RayAnalyzer:
    def __init__(self):
        self.out_dir = None
        self.dic_branch_name = dic_tree_branch_name

    def _fill_hist_in_tree(self, det, hist, branch_name: str, coord: str):
        if branch_name not in self.dic_branch_name.keys():
            raise Exception('branch_name must be in self.dic_branch_name')

        det_id = det.get_det_id()
        det_x = det.get_positions()[0]
        det_y = det.get_positions()[1]
        det_z = det.get_positions()[2]
        rot_phi = det.get_rotAngle()[2]

        xbins = hist.GetNbinsX()
        ybins = hist.GetNbinsY()

        tree_fname = 'tr_%s_%s.root' % (branch_name, det.get_det_name())
        tree_fname = os.path.join(self.out_dir, tree_fname)
        f_tree = TFile(tree_fname, 'recreate')
        tree = TTree('tr_%s' % branch_name, 'tr_%s' % branch_name)

        arr_det_id = array.array('i', [det_id])
        arr_det_x = array.array('d', [det_x])
        arr_det_y = array.array('d', [det_y])
        arr_det_z = array.array('d', [det_z])

        tree.Branch(self.dic_branch_name['det_id'], arr_det_id, '%s/I' % self.dic_branch_name['det_id'])
        tree.Branch(self.dic_branch_name['det_x'], arr_det_x, '%s/D' % self.dic_branch_name['det_x'])
        tree.Branch(self.dic_branch_name['det_y'], arr_det_y, '%s/D' % self.dic_branch_name['det_y'])
        tree.Branch(self.dic_branch_name['det_z'], arr_det_z, '%s/D' % self.dic_branch_name['det_z'])

        arr_angle_x = array.array('d', [-10.])
        arr_angle_y = array.array('d', [-10.])
        if coord.__contains__('phi'):
            tree.Branch(self.dic_branch_name['phi'], arr_angle_x, '%s/D' % self.dic_branch_name['phi'])
            tree.Branch(self.dic_branch_name['theta'], arr_angle_y, '%s/D' % self.dic_branch_name['theta'])
        else:
            tree.Branch(self.dic_branch_name['thetay'], arr_angle_x, '%s/D' % self.dic_branch_name['thetay'])
            tree.Branch(self.dic_branch_name['thetax'], arr_angle_y, '%s/D' % self.dic_branch_name['thetax'])

        arr_val = array.array('d', [-10.])
        arr_val_err = array.array('d', [-10.])
        tree.Branch(self.dic_branch_name[branch_name], arr_val,
                    '%s/D' % self.dic_branch_name[branch_name])
        if not branch_name.__contains__('count') and branch_name != 'ltopo' and not branch_name.__contains__('sig'):
            tree.Branch(self.dic_branch_name['%s_err' % branch_name], arr_val_err,
                        '%s/D' % self.dic_branch_name['%s_err' % branch_name])

        for i_xbin in range(1, xbins + 1):
            xx = hist.GetXaxis().GetBinCenter(i_xbin)
            for i_ybin in range(1, ybins + 1):
                yy = hist.GetYaxis().GetBinCenter(i_ybin)

                ibin = hist.GetBin(i_xbin, i_ybin)

                if coord.__contains__('phi'):
                    arr_angle_x[0] = CommonUtils.phi_cvt_in_minuspi_to_pi(xx + rot_phi)
                    arr_angle_y[0] = yy
                else:
                    vec = TVector3(TMath.Tan(yy), TMath.Tan(xx), 1)
                    _mag = vec.Mag()
                    _theta = vec.Theta()
                    _phi = CommonUtils.phi_cvt_in_minuspi_to_pi(vec.Phi() + rot_phi)
                    vec.SetMagThetaPhi(_mag, _theta, _phi)

                    _thetax = TMath.ATan(vec.X() / vec.Z())
                    _thetay = TMath.ATan(vec.Y() / vec.Z())

                    arr_angle_x[0] = _thetay
                    arr_angle_y[0] = _thetax

                arr_val[0] = hist.GetBinContent(ibin)
                arr_val_err[0] = hist.GetBinError(ibin)

                tree.Fill()

        tree.Write()
        f_tree.Save()
        f_tree.Close()

        return tree_fname

    def _fill_all_det_hists_into_one_tree(self, lis_det: list, lis_fname: list,
                                          hist_name, branch_name, coord, fout_name):
        if len(lis_det) != len(lis_fname):
            raise Exception('FillHistsIntoTree: number of xxRootNames '
                            'and DetNames must be the same.')

        lis_tree_fname = []

        for i in range(len(lis_det)):
            det = lis_det[i]
            fin_name = lis_fname[i]

            _hist_name = ROOTUtils.find_key_in_file(fin_name, hist_name)
            if branch_name == 'sig' and _hist_name.__contains__('smooth'):
                branch_name = 'smooth_sig'

            fin = TFile(fin_name, 'read')
            hist = fin.Get(_hist_name)

            tree_fname = self._fill_hist_in_tree(hist=hist, det=det, branch_name=branch_name, coord=coord)

            lis_tree_fname.append(tree_fname)

        str_lis_tree_fname = ' '.join(lis_tree_fname)
        branch_fname = os.path.splitext(fout_name)[0] + '_%s.root' % branch_name
        str_cmd = 'hadd -f %s %s' % (branch_fname, str_lis_tree_fname)

        if os.path.exists(branch_fname):
            os.remove(branch_fname)
        os.system(str_cmd)
        os.system('rm %s' % str_lis_tree_fname)

        return branch_fname

    def FillHistsIntoTree(self, dic_input: dict, fout_name: str):
        self.out_dir = os.path.split(fout_name)[0]

        lis_det = dic_input.get('Detectors')
        coord_name = dic_input.get('CoordName', 'thetaphi')

        lis_ltopo_fname = dic_input.get('LtopoRootNames', '')
        hist_name = dic_input.get('LtopoHistName', '')
        if hist_name == '':
            if coord_name.__contains__('phi'):
                hist_name = _hthetaphi_name_prefix
            else:
                hist_name = _hthetaxy_name_prefix
        if lis_ltopo_fname != '':
            branch_fname = self._fill_all_det_hists_into_one_tree(lis_det, lis_ltopo_fname, hist_name, 'ltopo',
                                                                  coord_name, fout_name)

        lis_meas_ratio_fname = dic_input.get('MeasRatioRootNames', '')
        hist_name = dic_input.get('MeasRatioHistName', '')
        if hist_name == '':
            if coord_name.__contains__('phi'):
                hist_name = _hthetaphi_name_prefix
            else:
                hist_name = _hthetaxy_name_prefix
        if lis_meas_ratio_fname != '':
            branch_fname = self._fill_all_det_hists_into_one_tree(lis_det, lis_meas_ratio_fname, hist_name,
                                                                  'meas_ratio', coord_name, fout_name)

        lis_sim_ratio_fname = dic_input.get('SimRatioRootNames', '')
        hist_name = dic_input.get('SimRatioHistName', '')
        if hist_name == '':
            if coord_name.__contains__('phi'):
                hist_name = _hthetaphi_name_prefix
            else:
                hist_name = _hthetaxy_name_prefix
        if lis_sim_ratio_fname != '':
            branch_fname = self._fill_all_det_hists_into_one_tree(lis_det, lis_sim_ratio_fname, hist_name, 'sim_ratio',
                                                                  coord_name, fout_name)

        lis_meas_leff_fname = dic_input.get('MeasLeffRootNames', '')
        hist_name = dic_input.get('MeasLeffHistName', '')
        if hist_name == '':
            if coord_name.__contains__('phi'):
                hist_name = _hthetaphi_name_prefix
            else:
                hist_name = _hthetaxy_name_prefix
        if lis_meas_leff_fname != '':
            branch_fname = self._fill_all_det_hists_into_one_tree(lis_det, lis_meas_leff_fname, hist_name, 'meas_leff',
                                                                  coord_name, fout_name)

        lis_sim_leff_fname = dic_input.get('SimLeffRootNames', '')
        hist_name = dic_input.get('SimLeffHistName', '')
        if hist_name == '':
            if coord_name.__contains__('phi'):
                hist_name = _hthetaphi_name_prefix
            else:
                hist_name = _hthetaxy_name_prefix
        if lis_sim_leff_fname != '':
            branch_fname = self._fill_all_det_hists_into_one_tree(lis_det, lis_sim_leff_fname, hist_name, 'sim_leff',
                                                                  coord_name, fout_name)

        lis_tg_stats_fname = dic_input.get('TargetStatsRootNames', '')
        hist_name = dic_input.get('TargetStatsHistName', '')
        if hist_name == '':
            if coord_name.__contains__('phi'):
                hist_name = _hthetaphi_name_prefix
            else:
                hist_name = _hthetaxy_name_prefix
        if lis_tg_stats_fname != '':
            tree = self._fill_all_det_hists_into_one_tree(lis_det, lis_tg_stats_fname, hist_name, 'meas_target_count',
                                                          coord_name, fout_name)

        lis_bg_stats_fname = dic_input.get('BackgroundStatsRootNames', '')
        hist_name = dic_input.get('BackgroundStatsHistName', '')
        if hist_name == '':
            if coord_name.__contains__('phi'):
                hist_name = _hthetaphi_name_prefix
            else:
                hist_name = _hthetaxy_name_prefix
        if lis_bg_stats_fname != '':
            branch_fname = self._fill_all_det_hists_into_one_tree(lis_det, lis_bg_stats_fname, hist_name,
                                                                  'meas_background_count', coord_name, fout_name)

        lis_tg_sim_count_fname = dic_input.get('TargetSimCountRootNames', '')
        hist_name = dic_input.get('TargetSimCountHistName', '')
        if hist_name == '':
            if coord_name.__contains__('phi'):
                hist_name = _hthetaphi_name_prefix
            else:
                hist_name = _hthetaxy_name_prefix
        if lis_tg_sim_count_fname != '':
            branch_fname = self._fill_all_det_hists_into_one_tree(lis_det, lis_tg_sim_count_fname, hist_name,
                                                                  'sim_target_count', coord_name, fout_name)

        lis_bg_sim_count_fname = dic_input.get('BackgroundSimCountRootNames', '')
        hist_name = dic_input.get('BackgroundSimCountHistName', '')
        if hist_name == '':
            if coord_name.__contains__('phi'):
                hist_name = _hthetaphi_name_prefix
            else:
                hist_name = _hthetaxy_name_prefix
        if lis_bg_sim_count_fname != '':
            branch_fname = self._fill_all_det_hists_into_one_tree(lis_det, lis_bg_sim_count_fname, hist_name,
                                                                  'sim_background_count', coord_name, fout_name)

        lis_sig_fname = dic_input.get('SigRootNames', '')
        hist_name = dic_input.get('SigHistName', '')
        if hist_name == '':
            if coord_name.__contains__('phi'):
                hist_name = _hthetaphi_name_prefix
            else:
                hist_name = _hthetaxy_name_prefix
        if lis_sig_fname != '':
            branch_fname = self._fill_all_det_hists_into_one_tree(lis_det, lis_sig_fname, hist_name,
                                                                  'sig', coord_name, fout_name)

    def _check_tree_can_be_friend(self, tr1, tr2):
        h1 = TH1F('h1', 'h1', 100, 0, 100)
        h2 = TH1F('h2', 'h2', 100, 0, 100)

        tr1.Draw("%s>>h1" % self.dic_branch_name['det_id'], '', 'goff')
        tr2.Draw("%s>>h2" % self.dic_branch_name['det_id'], '', 'goff')

        if h1.GetEntries() != h2.GetEntries():
            return False

        if h1.KolmogorovTest(h2) == 1:
            return True

        return False

    def _check_dataframe_can_be_friend(self, df1, df2):
        det_id1 = df1[self.dic_branch_name['det_id']].values
        h1 = ROOTUtils.fill_array_into_h1f(det_id1, 'h1', 'h1')
        det_id2 = df2[self.dic_branch_name['det_id']].values
        h2 = ROOTUtils.fill_array_into_h1f(det_id2, 'h2', 'h2')

        if h1.GetEntries() != h2.GetEntries():
            return False

        if h1.KolmogorovTest(h2) == 1:
            return True

        return False

    @staticmethod
    def _check_dataframe_can_continue(df1, df2):
        try:
            if_con = df1.columns == df2.columns
        except:
            return False
        if if_con:
            return True

        return False

    def MergeTreeBranch(self, lis_tree_fname: list, fout_name='', log_fname=None):
        lis_df = []
        if_fill_first = False
        first_tree_name = None
        tr_name = 'merge'
        print('MergeTreeBranch, tree Friends are:')
        for tree_name in lis_tree_fname:
            fin = TFile(tree_name, 'read')
            for i in fin.GetListOfKeys():
                val = fin.Get(i.GetName())
                if val.Class_Name() != 'TTree':
                    continue

                df = ROOTUtils.cvt_tree_to_dataframe(val)

                if not if_fill_first:
                    lis_df.append(df)
                    first_tree_name = i.GetName()
                    print(i.GetName())
                    tr_name += '_%s' % i.GetName()
                    if_fill_first = True
                else:
                    if not self._check_dataframe_can_be_friend(lis_df[0], df):
                        logger.warning('%s and %s can not be friend.' % (first_tree_name, i.GetName()))
                    else:
                        lis_df.append(df)
                        print(i.GetName())
                        tr_name += '_%s' % i.GetName()

        # add columns and save as tree
        df_merge = pd.DataFrame()
        for i_df in lis_df:
            for i_col in i_df.columns:
                if i_col in df_merge.columns:
                    continue
                if df_merge.shape[0] == 0:
                    df_merge[i_col] = i_df[i_col]
                else:
                    df_merge = pd.concat([df_merge, i_df[i_col]], axis=1)

        if fout_name == '':
            fout_name = os.path.join(self.out_dir, tr_name + '.root')
        ROOTUtils.cvt_dataframe_to_tree(df_merge, tr_name, fout_name)

        if log_fname is not None:
            print('\n\n', datetime.now(), file=open(log_fname, 'a'))
            print('MergeTreeBranch: %s\n' % ';'.join([os.path.split(i)[1] for i in lis_tree_fname]),
                  file=open(log_fname, 'a'))
            print('Save as %s' % fout_name,
                  file=open(log_fname, 'a'))

        return fout_name

    def ContinueTree(self, lis_tree_fname: list, lis_tree_name: list, fout_name='', log_fname=None):
        lis_df = []
        if_fill_first = False
        first_tree_name = None
        tr_name = 'merge'
        print('ContinueTree, trees are:')
        for i in range(len(lis_tree_fname)):
            tree_fname = lis_tree_fname[i]
            tree_name = lis_tree_name[i]
            tree_name = ROOTUtils.find_key_in_file(tree_fname, part_name=tree_name)

            fin = TFile(tree_fname, 'read')
            val = fin.Get(tree_name)
            print('%s in %s' % (tree_name, tree_fname))

            df = ROOTUtils.cvt_tree_to_dataframe(val)

            if i == 0:
                first_tree_name = tree_name
                lis_df.append(df)
            else:
                if_can_continue = self._check_dataframe_can_continue(lis_df[0], df)
                if not if_can_continue:
                    raise Exception('%s and %s can not link up.' % (first_tree_name, tree_name))
                else:
                    lis_df.append(df)
                    print(tree_name)
                    tr_name += '_%s' % tree_name

        # add columns and save as tree
        df_merge = pd.DataFrame()
        for i_df in lis_df:
            for i_col in i_df.columns:
                if df_merge.shape[0] == 0:
                    df_merge[i_col] = i_df[i_col]
                else:
                    df_merge = pd.concat([df_merge, i_df[i_col]], axis=0)

        if fout_name == '':
            fout_name = os.path.join(self.out_dir, tr_name + '.root')
        ROOTUtils.cvt_dataframe_to_tree(df_merge, tr_name, fout_name)

        if log_fname is not None:
            print('\n\n', datetime.now(), file=open(log_fname, 'a'))
            for i in range(len(lis_tree_fname)):
                tree_fname = lis_tree_fname[i]
                tree_name = lis_tree_name[i]
                print('ContinueTree: %s in %s\n' % (tree_name, tree_fname),
                      file=open(log_fname, 'a'))
            print('Save as %s' % fout_name,
                  file=open(log_fname, 'a'))

        return fout_name

    def CopyTreeWithSelection(self, selection: str, lis_tree_fname: list, fout_name: str, log_fname: str):
        self.out_dir = os.path.split(fout_name)[0]

        merge_fname = self.MergeTreeBranch(lis_tree_fname)
        tree_name = os.path.splitext(os.path.split(merge_fname)[1])[0]

        fin = TFile(merge_fname, 'read')
        tree = fin.Get(tree_name)

        fout = TFile(fout_name, 'recreate')
        tr_copy = tree.CopyTree(selection)
        tr_copy.SetName('tr_selected')
        tr_copy.Write()

        n1 = tree.GetEntries()
        n2 = tr_copy.GetEntries()
        print('%d entries has seleted from all %s entries' % (n2, n1))

        fout.Save()
        fout.Close()

        os.remove(merge_fname)

        print('Select %s with: "%s"' % (tree_name, selection))
        print('Save as %s' % fout_name)

        if log_fname is not None:
            print('\n\n', datetime.now(), file=open(log_fname, 'a'))
            print('CopyTreeWithSelection for %s\n' % ';'.join([os.path.split(i)[1] for i in lis_tree_fname]),
                  file=open(log_fname, 'a'))
            print('Select %s with: "%s"' % (tree_name, selection),
                  file=open(log_fname, 'a'))
            print('%d entries has seleted from all %s entries' % (n2, n1),
                  file=open(log_fname, 'a'))
            print('Save as %s' % fout_name,
                  file=open(log_fname, 'a'))

    @classmethod
    def _sample_poisson(cls, lam):
        rnd = np.random.poisson(lam)
        while rnd < 1:
            rnd = np.random.poisson(lam)
        return rnd

    def SampleFromMeasCount(self, dic_input, fout_name: str, log_fname: str = None):
        # noinspection DuplicatedCode
        self.out_dir = os.path.split(fout_name)[0]

        lis_tree_fname = dic_input.get('TreeRootNames')
        phys_model = dic_input.get('PhysModel')
        eps = float(dic_input.get('Epsilon'))

        merge_fname = self.MergeTreeBranch(lis_tree_fname)
        tree_name = os.path.splitext(os.path.split(merge_fname)[1])[0]

        fin = TFile(merge_fname, 'read')
        tree = fin.Get(tree_name)

        lis_col = [i.GetName() for i in tree.GetListOfBranches()]
        if (not self.dic_branch_name['meas_target_count'] in lis_col) \
                or (not self.dic_branch_name['meas_background_count'] in lis_col) \
                or (not self.dic_branch_name['meas_ratio'] in lis_col):
            raise Exception('At least input tr_meas_ratio, tr_meas_target_count and tr_meas_background_count.')

        df_tr = ROOTUtils.cvt_tree_to_dataframe(tree)

        new_mtx = np.zeros([df_tr.shape[0], 4])
        new_mtx[:, :] = -1.

        for i in range(df_tr.shape[0]):
            if (i > df_tr.shape[0] - 5) or (i % 200 == 0):
                print('\rSampleFromMeasCount: %d / %d = %.2f   ' % (i, df_tr.shape[0], 1. * i / df_tr.shape[0]), end='')

            df_row = df_tr.iloc[i]

            if (not df_row[self.dic_branch_name['meas_target_count']] > 5) \
                    or (not df_row[self.dic_branch_name['meas_background_count']] > 5):
                continue
            rnd_tg_count = self._sample_poisson(df_row[self.dic_branch_name['meas_target_count']])
            rnd_bg_count = self._sample_poisson(df_row[self.dic_branch_name['meas_background_count']])
            time_ratio = df_row[self.dic_branch_name['meas_ratio']] \
                         / df_row[self.dic_branch_name['meas_target_count']] \
                         * df_row[self.dic_branch_name['meas_background_count']]

            r, e = CommonUtils.calc_divide_and_err(rnd_tg_count, np.sqrt(rnd_tg_count), rnd_bg_count,
                                                   np.sqrt(rnd_bg_count))
            ratio = r * time_ratio
            ratio_err = e * time_ratio
            new_mtx[i, 0] = ratio
            new_mtx[i, 1] = ratio_err

            if self.dic_branch_name['phi'] in df_tr.columns:
                theta = df_row[self.dic_branch_name['theta']]
            else:
                thetax = df_row[self.dic_branch_name['thetax']]
                thetay = df_row[self.dic_branch_name['thetay']]
                theta = CommonUtils.cvt_thetaxy_to_thetaphi(thetax=thetax, thetay=thetay).get('theta')

            leff = phys_model.get_Leff_from_ratio(theta, ratio)
            dl_dr = phys_model.get_demin_dratio(theta, ratio, eps, leff)
            err_leff = dl_dr * ratio_err

            new_mtx[i, 2] = leff
            new_mtx[i, 3] = err_leff
        print()

        new_cols = [self.dic_branch_name['rnd_meas_ratio'],
                    self.dic_branch_name['rnd_meas_ratio_err'],
                    self.dic_branch_name['rnd_meas_leff'],
                    self.dic_branch_name['rnd_meas_leff_err']]
        df_tr[new_cols] = new_mtx

        df_tr = df_tr[df_tr[new_cols[0]] > 0]

        ROOTUtils.cvt_dataframe_to_tree(df_tr, 'tr_rnd_meas_leff', fout_name)

        print(df_tr[new_cols])
        os.remove(merge_fname)

        if log_fname is not None:
            print('\n\n', datetime.now(), file=open(log_fname, 'a'))
            print('SampleFromMeasCount: %s\n' % ';'.join([os.path.split(i)[1] for i in lis_tree_fname]),
                  file=open(log_fname, 'a'))
            print('Save as %s' % fout_name,
                  file=open(log_fname, 'a'))

    def SampleFromSimCount(self, dic_input, fout_name: str, log_fname: str = None):
        # noinspection DuplicatedCode
        self.out_dir = os.path.split(fout_name)[0]

        tree_fname = dic_input.get('TreeRootName')
        phys_model = dic_input.get('PhysModel')
        eps = float(dic_input.get('Epsilon'))

        tree_name = dic_input.get('TreeName', '')
        if tree_name == '':
            tree_name = ROOTUtils.find_key_in_file(tree_fname, part_name='tr')
        if tree_name == '':
            raise Exception('Please set a correct TreeName.')

        fin = TFile(tree_fname, 'read')
        tree = fin.Get(tree_name)

        lis_col = [i.GetName() for i in tree.GetListOfBranches()]
        if (not self.dic_branch_name['meas_target_count'] in lis_col) \
                or (not self.dic_branch_name['meas_background_count'] in lis_col) \
                or (not self.dic_branch_name['sim_ratio'] in lis_col) \
                or (not self.dic_branch_name['meas_ratio'] in lis_col):
            raise Exception(
                'At least input tr_meas_ratio, tr_sim_ratio, tr_meas_target_count and tr_meas_background_count.')

        df_tr = ROOTUtils.cvt_tree_to_dataframe(tree)

        new_mtx = np.zeros([df_tr.shape[0], 4])
        new_mtx[:, :] = -1.

        for i in range(df_tr.shape[0]):
            if (i > df_tr.shape[0] - 5) or (i % 200 == 0):
                print('\rSampleFromSimCount: %d / %d = %.2f   ' % (i, df_tr.shape[0], 1. * i / df_tr.shape[0]), end='')

            df_row = df_tr.iloc[i]

            if (not df_row[self.dic_branch_name['meas_target_count']] > 5) \
                    or (not df_row[self.dic_branch_name['meas_background_count']] > 5):
                continue
            sim_ratio = df_row[self.dic_branch_name['sim_ratio']]
            meas_ratio = df_row[self.dic_branch_name['meas_ratio']]
            bg_count = df_row[self.dic_branch_name['meas_background_count']]
            rnd_bg_count = self._sample_poisson(bg_count)
            time_ratio = meas_ratio / df_row[self.dic_branch_name['meas_target_count']] \
                         * df_row[self.dic_branch_name['meas_background_count']]
            rnd_tg_count = self._sample_poisson(bg_count * sim_ratio / time_ratio)

            # noinspection DuplicatedCode
            r, e = CommonUtils.calc_divide_and_err(rnd_tg_count, np.sqrt(rnd_tg_count), rnd_bg_count,
                                                   np.sqrt(rnd_bg_count))
            ratio = r * time_ratio
            ratio_err = e * time_ratio
            new_mtx[i, 0] = ratio
            new_mtx[i, 1] = ratio_err

            if self.dic_branch_name['phi'] in df_tr.columns:
                theta = df_row[self.dic_branch_name['theta']]
            else:
                thetax = df_row[self.dic_branch_name['thetax']]
                thetay = df_row[self.dic_branch_name['thetay']]
                theta = CommonUtils.cvt_thetaxy_to_thetaphi(thetax=thetax, thetay=thetay).get('theta')

            leff = phys_model.get_Leff_from_ratio(theta, ratio)
            dl_dr = phys_model.get_demin_dratio(theta, ratio, eps, leff)
            err_leff = dl_dr * ratio_err

            new_mtx[i, 2] = leff
            new_mtx[i, 3] = err_leff
        print()

        new_cols = [self.dic_branch_name['rnd_sim_ratio'],
                    self.dic_branch_name['rnd_sim_ratio_err'],
                    self.dic_branch_name['rnd_sim_leff'],
                    self.dic_branch_name['rnd_sim_leff_err']]
        df_tr[new_cols] = new_mtx

        df_tr = df_tr[df_tr[new_cols[0]] > 0]

        ROOTUtils.cvt_dataframe_to_tree(df_tr, 'tr_rnd_sim_ratio_leff', fout_name)

        print(df_tr[new_cols])

        if log_fname is not None:
            print('\n\n', datetime.now(), file=open(log_fname, 'a'))
            print('SampleFromSimCount: %s in %s\n' % (tree_name, tree_fname),
                  file=open(log_fname, 'a'))
            print('Save as %s' % fout_name,
                  file=open(log_fname, 'a'))

    def TreeThetaxy2Thetaphi(self, fin_name, tree_name, fout_name, if_keep_thetaxy=True):
        if tree_name == '':
            tree_name = 'tr_'

        tree_name = ROOTUtils.find_key_in_file(fin_name, tree_name)

        fin = TFile(fin_name, 'read')
        tree = fin.Get(tree_name)

        df_tr = ROOTUtils.cvt_tree_to_dataframe(tree)

        if (not self.dic_branch_name['thetax'] in df_tr.columns) \
                or (not self.dic_branch_name['thetay'] in df_tr.columns):
            raise Exception('Cannot find %s or %s in %s'
                            % (self.dic_branch_name['thetax'], self.dic_branch_name['thetay'], tree_name))

        new_mtx = np.zeros([df_tr.shape[0], 2])
        for i in range(df_tr.shape[0]):
            df_row = df_tr.iloc[i]
            thetax = df_row[self.dic_branch_name['thetax']]
            thetay = df_row[self.dic_branch_name['thetay']]

            dic_thetaphi = CommonUtils.cvt_thetaxy_to_thetaphi(thetax=thetax, thetay=thetay)
            theta = dic_thetaphi.get('theta')
            phi = dic_thetaphi.get('phi')

            new_mtx[i, 0] = theta
            new_mtx[i, 1] = phi

        df_tr[[self.dic_branch_name['theta'], self.dic_branch_name['phi']]] = new_mtx
        if not if_keep_thetaxy:
            df_tr = df_tr.drop(columns=[self.dic_branch_name['thetax'], self.dic_branch_name['thetay']])

        ROOTUtils.cvt_dataframe_to_tree(df_tr, tree_name, fout_name)

    @classmethod
    def _get_phys_model_1d_theta(cls, phys_model, theta, lmin=0, lmax=200):
        return phys_model.get_phys_model_1d_theta(theta, lmin, lmax)

    def ComparePhysicalModelAndMeasuredData(self, fin_name, tr_name, lis_det_id: list, lis_phys_model: list, fout_name,
                                            bins_dic):
        # get the number of theta bins
        theta_bins = bins_dic['theta_bins'][0]
        theta_min = bins_dic['theta_bins'][1]
        theta_max = bins_dic['theta_bins'][2]

        tr_name = ROOTUtils.find_key_in_file(fin_name, tr_name)
        fin = TFile(fin_name, 'read')
        tr = fin.Get(tr_name)

        fout = TFile(fout_name, 'recreate')
        gr_3D_multi = TMultiGraph('mg_3D_RLTheta',
                                  'Physical Model for all selected detectors; '
                                  '#theta [rad]; '
                                  'Effective Length [m g#dotc cm^{-3}]; '
                                  'Measured Survival Rate'
                                  )

        # fill 3D graph
        legend1 = TLegend(0.85, 0.7, 0.9, 0.9)
        can1 = TCanvas("can_gr_3D_multi", "can_gr_3D_multi", 2880, 1440)
        if_draw = False
        for i_det_id in lis_det_id:
            if isinstance(i_det_id, str):
                i_det_id = int(i_det_id)
            tr_i_det = tr.CopyTree("%s == %d" % (self.dic_branch_name['det_id'], i_det_id))

            dic_gr_3D = {
                'arr_theta': [],
                'arr_theta_err': [],
                'arr_sim_length': [],
                'arr_sim_length_err': [],
                'arr_meas_ratio': [],
                'arr_meas_ratio_err': []
            }

            for iev in tr_i_det:
                theta = getattr(iev, self.dic_branch_name['theta'])
                theta_err = (theta_max - theta_min) / theta_bins / 2
                meas_ratio = getattr(iev, self.dic_branch_name['meas_ratio'])
                meas_ratio_err = getattr(iev, self.dic_branch_name['meas_ratio_err'])
                sim_length = getattr(iev, self.dic_branch_name['sim_leff'])
                sim_length_err = 1

                dic_gr_3D['arr_theta'].append(theta)
                dic_gr_3D['arr_theta_err'].append(theta_err)
                dic_gr_3D['arr_sim_length'].append(sim_length)
                dic_gr_3D['arr_sim_length_err'].append(sim_length_err)
                dic_gr_3D['arr_meas_ratio'].append(meas_ratio)
                dic_gr_3D['arr_meas_ratio_err'].append(meas_ratio_err)

            # fill arrays into graphs and save them in file
            arr_theta = array.array('d', dic_gr_3D['arr_theta'])
            arr_theta_err = array.array('d', dic_gr_3D['arr_theta_err'])
            arr_sim_length = array.array('d', dic_gr_3D['arr_sim_length'])
            arr_sim_length_err = array.array('d', dic_gr_3D['arr_sim_length_err'])
            arr_meas_ratio = array.array('d', dic_gr_3D['arr_meas_ratio'])
            arr_meas_ratio_err = array.array('d', dic_gr_3D['arr_meas_ratio_err'])
            gr_3D = TGraph2DErrors(len(dic_gr_3D['arr_theta']),
                                   arr_theta, arr_sim_length, arr_meas_ratio,
                                   arr_theta_err, arr_sim_length_err, arr_meas_ratio_err)
            gr_3D.SetNameTitle('gr_3D_RLTheta_%d' % i_det_id,
                               'Physical Model for Detector %d; '
                               '#theta [rad]; '
                               'Effective Length [m g#dotc cm^{-3}]; '
                               'Measured Survival Rate' % i_det_id)
            fout.cd()
            gr_3D.Write()

            # gr_3D draw options
            gr_3D.SetMarkerColor(i_det_id + 1)
            gr_3D.SetMarkerStyle(21)
            if if_draw:
                gr_3D.Draw('p same')
            else:
                gr_3D.Draw('p')
                if_draw = True
            gr_3D.SetMaximum(2)
            # gr_3D_multi.Add(gr_3D)
            # gr_3D_multi.Draw('a')
            legend1.AddEntry(gr_3D_multi, "Detector %d" % i_det_id, "p")

        # gr_3D_multi.Write()
        legend1.Write()

        # gr_3D_multi.Draw('p')
        legend1.Draw()
        can1.Write()

        # fill 2D graphs
        # get 2D data from tree
        lis_gr_2D_det = []
        for i_det_id in lis_det_id:
            if isinstance(i_det_id, str):
                i_det_id = int(i_det_id)
            tr_i_det = tr.CopyTree("%s == %d" % (self.dic_branch_name['det_id'], i_det_id))

            dic_gr_2D = {}
            for itheta in range(theta_bins):
                dic_gr_2D['gr_data_%d' % itheta] = {
                    'arr_theta': [],
                    'arr_theta_err': [],
                    'arr_sim_length': [],
                    'arr_sim_length_err': [],
                    'arr_meas_ratio': [],
                    'arr_meas_ratio_err': []
                }

            for iev in tr_i_det:
                theta = getattr(iev, self.dic_branch_name['theta'])
                theta_err = (theta_max - theta_min) / theta_bins / 2
                meas_ratio = getattr(iev, self.dic_branch_name['meas_ratio'])
                meas_ratio_err = getattr(iev, self.dic_branch_name['meas_ratio_err'])
                sim_length = getattr(iev, self.dic_branch_name['sim_leff'])
                sim_length_err = sim_length / 200  # 0

                itheta = int((theta - theta_min) * theta_bins / (theta_max - theta_min))
                dic_gr_2D['gr_data_%d' % itheta]['arr_theta'].append(theta)
                dic_gr_2D['gr_data_%d' % itheta]['arr_theta_err'].append(theta_err)
                dic_gr_2D['gr_data_%d' % itheta]['arr_sim_length'].append(sim_length)
                dic_gr_2D['gr_data_%d' % itheta]['arr_sim_length_err'].append(sim_length_err)
                dic_gr_2D['gr_data_%d' % itheta]['arr_meas_ratio'].append(meas_ratio)
                dic_gr_2D['gr_data_%d' % itheta]['arr_meas_ratio_err'].append(meas_ratio_err)

            lis_gr_2D_det.append(dic_gr_2D)

        # fill data into graphs and draw them
        for itheta in range(theta_bins):
            theta_center = (theta_max - theta_min) / theta_bins * (itheta + 0.5)
            mg_2D_itheta = TMultiGraph('mg_2D_itheta%d' % itheta,
                                       '#theta ~ %.3f; '
                                       'Effective Length [m g cm^{-3}]; '
                                       'Measured Survival Rate' % theta_center)

            legend2 = TLegend(0.8, 0.7, 0.9, 0.9)
            legend2.SetName('lg_mg_2D_itheta%d' % itheta)
            if_has_data = False
            for i in range(len(lis_det_id)):
                i_det_id = lis_det_id[i]
                det_data = lis_gr_2D_det[i]
                det_itheta_data = det_data['gr_data_%d' % itheta]

                n = len(det_itheta_data['arr_sim_length'])
                if n < 1:
                    continue
                if_has_data = True
                arr_sim_length = array.array('d', det_itheta_data['arr_sim_length'])
                arr_sim_length_err = array.array('d', det_itheta_data['arr_sim_length_err'])
                arr_meas_ratio = array.array('d', det_itheta_data['arr_meas_ratio'])
                arr_meas_ratio_err = array.array('d', det_itheta_data['arr_meas_ratio_err'])

                gr_det_itheta = TGraphErrors(n, arr_sim_length, arr_meas_ratio, arr_sim_length_err, arr_meas_ratio_err)
                gr_det_itheta.SetName('gr_2D_itheta%d_det%d' % (itheta, i_det_id))
                gr_det_itheta.SetTitle('#theta ~ %.3f; '
                                       'Effective Length [m g cm^{-3}]; '
                                       'Measured Survival Rate' % theta_center)

                # gr_2D draw options
                gr_det_itheta.SetMarkerColor(i_det_id + 1)
                gr_det_itheta.SetMarkerStyle(21)
                gr_det_itheta.SetMarkerSize(1.8)
                gr_det_itheta.SetLineColor(i_det_id + 1)
                bar_width = (gr_det_itheta.GetXaxis().GetXmax() - gr_det_itheta.GetXaxis().GetXmin()) / 20
                gr_det_itheta.SetFillColorAlpha(i_det_id + 1, 0.5)
                gr_det_itheta.SetFillStyle(3001)
                fout.cd()
                gr_det_itheta.Write()

                mg_2D_itheta.Add(gr_det_itheta, 'p2')
                legend2.AddEntry(gr_det_itheta, "Detector %d" % i_det_id, "lp")

            if not if_has_data:
                continue
            mg_2D_itheta.Write('a')
            legend2.Write()

            # calculate physical model in the graphs
            lis_f1 = []
            ip = 0
            for iphys in lis_phys_model:
                f1 = self._get_phys_model_1d_theta(iphys, theta_center,
                                                   mg_2D_itheta.GetXaxis().GetXmin(), mg_2D_itheta.GetXaxis().GetXmax())
                f1.SetLineColor(20 + ip)
                # f1.Write()
                lis_f1.append(f1)
                ip += 1

            can2 = TCanvas("can_gr_2D_multi_%d" % itheta, "can_gr_2D_multi_%d" % itheta, 2880, 1440)
            mg_2D_itheta.SetMaximum(1.1)
            mg_2D_itheta.SetMinimum(0)

            mg_2D_itheta.GetXaxis().SetTitleSize(0.05)
            mg_2D_itheta.GetYaxis().SetTitleSize(0.05)
            mg_2D_itheta.GetXaxis().SetTitleOffset(0.7)
            mg_2D_itheta.GetYaxis().SetTitleOffset(0.7)

            mg_2D_itheta.Draw("a")
            # for if1 in lis_f1:
            #     if1.Draw("same")
            legend2.Draw()
            can2.Write()
