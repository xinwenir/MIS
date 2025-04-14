##########################################
# Statistics is used for collect all the events in the given tree
# Developed by Guorui Liu
import os.path

import pandas as pd
from ROOT import TH2F, TH1F, TMath, TTree, TFile, TCanvas, THStack, TLegend, TVectorD, TH3F

from DataFactory import CommonUtils, ROOTUtils

import logging

logger = logging.getLogger(__name__)


class MuonStatistics:

    def __init__(self):
        self.h1f_varnames_1layer = ['xo1', 'zx1', 'yo1', 'zy1', 'xo2', 'zx2', 'yo2', 'zy2']
        self.h1f_binnames_1layer = ['x_bins', 'z_bins', 'y_bins', 'z_bins', 'x_bins', 'z_bins', 'y_bins', 'z_bins']
        self.h1f_cuts_1layer = ['',
                                '_R_NBar#==1',
                                '_R_NBar#==2',
                                '_R_NBar#==3',
                                '_R_NBar#==4']
        self.h1f_cuts_name_1layer = ['TOT', 'S', 'D', 'T', 'Q']  # total, double, single, triple, and quadruple

        self.h1f_varnames_2layer_Z = ['Z1', 'Z2']
        self.h1f_varnames_2layer = ['X1', 'Y1', 'X2', 'Y2', 'thetax_trk', 'thetay_trk']
        self.h1f_binnames_2layer = ['x_bins', 'y_bins', 'x_bins', 'y_bins', 'thetax_bins', 'thetay_bins']
        self.h1f_cuts_2layer = ['',
                                '_R_NBar1>1 && _R_NBar2>1',
                                '_R_NBar1==1 && _R_NBar2==1',
                                '(_R_NBar1==1 && _R_NBar2>1) || (_R_NBar1>1 && _R_NBar2==1)']
        self.h1f_cuts_name_2layer = ['TOT', 'MM', 'SS', 'MS']  # M: multiple

        self.h1f_varnames_4layer = ['theta', 'phi', 'thetax', 'thetay']
        self.h1f_binnames_4layer = ['theta_bins', 'phi_bins', 'thetax_bins', 'thetay_bins']

        self.h2f_varnames = ['X1:Y1', 'X2:Y2', 'theta:phi', 'thetax:thetay', 'thetax_trk:thetay_trk',
                             'xo1:yo1', 'xo2:yo2', 'xo1:zx1', 'yo1:zy1', 'xo2:zx2', 'yo2:zy2']
        self.h2f_binnames = ['x_bins:y_bins', 'x_bins:y_bins',
                             'theta_bins:phi_bins', 'thetax_bins:thetay_bins', 'thetax_bins:thetay_bins',
                             'x_bins:y_bins', 'x_bins:y_bins',
                             'x_bins:z_bins', 'y_bins:z_bins',
                             'x_bins:z_bins', 'y_bins:z_bins']
        self.delim = ':'
        self.h2f_cuts = ['',
                         'xNBar1>1 && xNBar2>1 && yNBar1>1 && yNBar2>1',
                         'xNBar1==1 && xNBar2==1 && yNBar1==1 && yNBar2==1']
        self.h2f_cuts_name = ['TOT', 'MMMM', 'SSSS']

    @classmethod
    def _replace_cut(cls, cut, varname):
        if varname.lower().__contains__('x'):
            cut = cut.replace('_R_', 'x')
        elif varname.lower().__contains__('y'):
            cut = cut.replace('_R_', 'y')

        if varname.__contains__('1'):
            cut = cut.replace('NBar#', 'NBar1')
        elif varname.__contains__('2'):
            cut = cut.replace('NBar#', 'NBar2')

        return cut

    @classmethod
    def _calc_statistic_uncertainty(cls, h):
        h1 = h.Clone()
        bins = h1.GetNbinsX() * h1.GetNbinsY() * h1.GetNbinsZ()
        for i in range(1, bins + 1):
            val = h1.GetBinContent(i)
            err = 0
            if val > 0:
                err = TMath.Sqrt(val)
            h1.SetBinError(i, err)
        return h1

    def _draw_h1f(self, tree, varname, h_name, h_title, bins, cut):
        if bins[1] > bins[2]:
            bins[1] = int(1)
            bins[2] = int(-1)
        h1f = TH1F(h_name, h_title, bins[0], bins[1], bins[2])
        tree.Draw("%s>>%s" % (varname, h_name), cut, "goff")
        self._calc_statistic_uncertainty(h1f)
        return h1f

    @classmethod
    def _draw_h_stack(cls, h_lis, hs_name, hs_title):
        hs = THStack(hs_name, hs_title)
        lg = TLegend(0.72, 0.8, 0.9, 0.9)
        lg.SetName('lg_' + hs_name)
        lg.SetFillStyle(0)

        for i in range(len(h_lis)):
            h_lis[i].SetLineColor(1 + i)
            hs.Add(h_lis[i], "hist e")
            lg.AddEntry(h_lis[i], h_lis[i].GetTitle())

        return hs, lg

    def DoStatistics1D_1Layer(self, bins_dic, root_in_fname, root_out_fname):
        """tree draw h1f using only one layer"""
        logger.info('------------ HIST 1D 1Layer ------------')
        fin = TFile(root_in_fname, 'read')
        meas_time = fin.Get('meas_time')
        tree = fin.Get('t_tree')

        fout = TFile(root_out_fname, 'UPDATE')

        # draw variables only use one layer information
        n_vars = len(self.h1f_varnames_1layer)
        n_cuts = len(self.h1f_cuts_1layer)

        for i in range(n_vars):
            var_name = self.h1f_varnames_1layer[i]
            bin_name = self.h1f_binnames_1layer[i]
            lis_h1f = []
            var_unit = 'mm'

            for j in range(n_cuts):
                tmp_cuts = self.h1f_cuts_1layer[j]
                cut_name = self.h1f_cuts_name_1layer[j]
                tmp_bins = bins_dic[bin_name]
                tmp_cuts = self._replace_cut(tmp_cuts, var_name)

                h_name = 'h{}_{}'.format(var_name, cut_name)
                h_title = h_name + '; {} [{}]; counts'.format(var_name, var_unit)
                logger.debug('draw: %s %s' % (h_name, h_title))

                h1f = self._draw_h1f(tree=tree, varname=var_name, h_name=h_name, h_title=h_title,
                                     bins=tmp_bins, cut=tmp_cuts)
                h1f.Write()
                lis_h1f.append(h1f)
            h_lis = lis_h1f[1:]
            hs_name = 'hs_' + var_name
            hs_title = hs_name + '; {} [{}]; counts'.format(var_name, var_unit)
            hs, lg = self._draw_h_stack(h_lis, hs_name, hs_title)
            hs.Write()
            lg.Write()

        if not 'meas_time' in fout.GetListOfKeys():
            meas_time.Write('meas_time')

        fin.Close()
        fout.Save()
        fout.Close()

    def DoStatistics1D_2Layer(self, bins_dic, root_in_fname, root_out_fname):
        """tree draw h1f using two layers' information"""
        logger.info('------------ HIST 1D 2Layer ------------')
        fin = TFile(root_in_fname, 'read')
        meas_time = fin.Get('meas_time')
        tree = fin.Get('t_tree')

        fout = TFile(root_out_fname, 'UPDATE')

        # draw Z1,Z2
        for var_name in self.h1f_varnames_2layer_Z:
            h_name = 'h{}_{}'.format(var_name, 'TOT')
            h_title = h_name + '; {} [{}]; counts'.format(var_name, 'mm')
            logger.debug('draw: %s %s' % (h_name, h_title))

            h1f = self._draw_h1f(tree, varname=var_name, h_name=h_name, h_title=h_title,
                                 bins=bins_dic['z_bins'], cut='')
            h1f.Write()

        # draw variables only use two layer information
        n_vars = len(self.h1f_varnames_2layer)
        n_cuts = len(self.h1f_cuts_2layer)

        for i in range(n_vars):
            bin_name = self.h1f_binnames_2layer[i]
            var_name = self.h1f_varnames_2layer[i]
            var_unit = 'mm'
            if var_name.lower().__contains__('theta') or var_name.lower().__contains__('phi'):
                var_unit = 'rad'

            lis_h1f = []
            for j in range(n_cuts):
                tmp_cuts = self.h1f_cuts_2layer[j]
                cut_name = self.h1f_cuts_name_2layer[j]
                tmp_bins = bins_dic[bin_name]

                tmp_cuts = self._replace_cut(tmp_cuts, var_name)

                h_name = 'h{}_{}'.format(var_name, cut_name)
                h_title = h_name + '; {} [{}]; counts'.format(var_name, var_unit)
                logger.debug('draw: %s %s' % (h_name, h_title))

                h1f = self._draw_h1f(tree, varname=var_name, h_name=h_name, h_title=h_title,
                                     bins=tmp_bins, cut=tmp_cuts)

                h1f.Write()
                lis_h1f.append(h1f)
            h_lis = lis_h1f[1:]
            hs_name = 'hs_' + var_name
            hs_title = hs_name + '; {} [{}]; counts'.format(var_name, var_unit)
            hs, lg = self._draw_h_stack(h_lis, hs_name, hs_title)
            hs.Write()
            lg.Write()

        if not 'meas_time' in fout.GetListOfKeys():
            meas_time.Write('meas_time')

        fin.Close()
        fout.Save()
        fout.Close()

    def DoStatistics1D_4Layer(self, bins_dic, root_in_fname, root_out_fname):
        """tree draw h1f using four layers' information"""
        logger.info('------------ HIST 1D 4Layer ------------')
        fin = TFile(root_in_fname, 'read')
        meas_time = fin.Get('meas_time')
        tree = fin.Get('t_tree')

        fout = TFile(root_out_fname, 'UPDATE')

        # draw variables only use one layer information
        n_vars = len(self.h1f_varnames_4layer)

        for i in range(n_vars):
            bin_name = self.h1f_binnames_4layer[i]
            var_name = self.h1f_varnames_4layer[i]
            var_unit = 'rad'

            tmp_cuts = ''
            cut_name = 'TOT'
            tmp_bins = bins_dic[bin_name]

            h_name = 'h{}_{}'.format(var_name, cut_name)
            h_title = h_name + '; {} [{}]; counts'.format(var_name, var_unit)
            logger.debug('draw: %s %s' % (h_name, h_title))

            h1f = self._draw_h1f(tree, varname=var_name, h_name=h_name, h_title=h_title,
                                 bins=tmp_bins, cut=tmp_cuts)

            h1f.Write()

        if not 'meas_time' in fout.GetListOfKeys():
            meas_time.Write('meas_time')

        fin.Close()
        fout.Save()
        fout.Close()

    @classmethod
    def _draw_h2f(cls, tree, varname, h_name, h_title, bins1, bins2, cut):
        if bins1[1] > bins1[2]:
            bins1[1] = int(1)
            bins1[2] = int(-1)
        if bins2[1] > bins2[2]:
            bins2[1] = int(1)
            bins2[2] = int(-1)
        h2f = TH2F(h_name, h_title, bins1[0], bins1[1], bins1[2], bins2[0], bins2[1], bins2[2])
        tree.Draw("%s>>%s" % (varname, h_name), cut, "goff")
        return h2f

    @classmethod
    def _hist_divide(cls, h1, h2, h_name):
        if h1 is None or h2 is None:
            raise Exception(__name__ + ' hist_divide, two hists should be passed here.')
        if not ROOTUtils.check_two_hists_the_same_size(h1=h1, h2=h2):
            raise ValueError(__name__ + ' hist_divide, passing two hists in different size.')

        h3 = h1.Clone()
        h3.Reset("ICES")
        h3.SetNameTitle(h_name, h_name)

        bins = h1.GetNbinsX() * h1.GetNbinsY() * h1.GetNbinsZ()

        for i in range(1, bins + 1):
            val1 = h1.GetBinContent(i)
            err1 = h1.GetBinError(i)
            val2 = h2.GetBinContent(i)
            err2 = h2.GetBinError(i)

            if val2 < 1e-10:
                h3.SetBinContent(i, 0.)
                h3.SetBinError(i, 0.)
            else:
                val3 = val1 / val2
                term1 = TMath.Power(err1 / val2, 2)
                term2 = TMath.Power(val1 * err2 / val2 / val2, 2)
                err3 = TMath.Sqrt(term1 + term2)

                h3.SetBinContent(i, val3)
                h3.SetBinError(i, err3)

        return h3

    def DoStatistics2D(self, bins_dic, root_in_fname, root_out_fname,
                       det, h_thetaphi_solid_angle, h_thetaxy_solid_angle):
        """tree draw h2f"""
        logger.info('------------ HIST 2D ------------')
        fin = TFile(root_in_fname, 'read')
        meas_time = fin.Get('meas_time')
        tree = fin.Get('t_tree')

        fout = TFile(root_out_fname, 'UPDATE')

        n_vars = len(self.h2f_varnames)
        n_cuts = len(self.h2f_cuts)

        for i in range(n_vars):
            var_unit = 'mm'
            if self.h2f_varnames[i].__contains__('theta'):
                var_unit = 'rad'
            var_name = self.h2f_varnames[i]
            var_lis = self.h2f_varnames[i].split(self.delim)
            bin_lis = self.h2f_binnames[i].split(self.delim)

            var1 = var_lis[0]
            var2 = var_lis[1]
            bins1 = bins_dic[bin_lis[1]]
            bins2 = bins_dic[bin_lis[0]]

            for j in range(n_cuts):
                tmp_cuts = self.h2f_cuts[j]
                cut_name = self.h2f_cuts_name[j]

                h_name = 'h{}_{}_{}'.format(var1, var2, cut_name)
                h_title = h_name + '; {} [{}]; {} [{}]; counts'.format(var2, var_unit, var1, var_unit)

                logger.debug('draw: %s %s' % (h_name, h_title))

                h2f = self._draw_h2f(tree, varname=var_name, h_name=h_name, h_title=h_title,
                                     bins1=bins1, bins2=bins2, cut=tmp_cuts)
                h2f.Write()

                # calc hist_div_acceptance and _div_solid_angle
                if j == 0:
                    h_name_acc = '%s_div_acceptance' % h_name
                    h_name_acc_solid_angle = '%s_div_acceptance_div_solidAngle' % h_name
                    h_accept_thetaphi, h_accept_thetaxy = det.get_acceptance_plane_det()

                    t_h2f = h2f.Clone()

                    logger.debug('draw: %s %s' % (h_name_acc, h_name_acc_solid_angle))

                    if var_name == 'theta:phi':
                        h2f_acc = self._hist_divide(h1=t_h2f, h2=h_accept_thetaphi, h_name=h_name_acc)
                        h2f_acc.Write()
                        h2f_acc_solid_angle = self._hist_divide(h1=h2f_acc, h2=h_thetaphi_solid_angle,
                                                                h_name=h_name_acc_solid_angle)
                        h2f_acc_solid_angle.Write()

                    elif var_name == 'thetax:thetay':
                        h2f_acc = self._hist_divide(h1=t_h2f, h2=h_accept_thetaxy, h_name=h_name_acc)
                        h2f_acc.Write()
                        h2f_acc_solid_angle = self._hist_divide(h1=h2f_acc, h2=h_thetaxy_solid_angle,
                                                                h_name=h_name_acc_solid_angle)
                        h2f_acc_solid_angle.Write()
                    else:
                        continue

        if not 'meas_time' in fout.GetListOfKeys():
            meas_time.Write('meas_time')

        fin.Close()
        fout.Save()
        fout.Close()

    @classmethod
    def _draw_h3f(cls, tree, varname, h_name, h_title, bins1, bins2, bins3, cut):
        if bins1[1] > bins1[2]:
            bins1[1] = int(1)
            bins1[2] = int(-1)
        if bins2[1] > bins2[2]:
            bins2[1] = int(1)
            bins2[2] = int(-1)
        if bins3[1] > bins3[2]:
            bins3[1] = int(1)
            bins3[2] = int(-1)
        h3f = TH3F(h_name, h_title, bins1[0], bins1[1], bins1[2],
                   bins2[0], bins2[1], bins2[2],
                   bins3[0], bins3[1], bins3[2])
        tree.Draw("%s>>%s" % (varname, h_name), cut, "goff")
        return h3f

    @classmethod
    def _get_tree_in_time(cls, tree, time_expression, itime):
        tr = tree.CopyTree(time_expression)

        tr.SetName("tr_slice_%d" % itime)

        return tr

    def _get_lis_time_cut(self, time_interval, time_fname):
        # get the time slice list
        df_meas_time = pd.read_csv(time_fname)
        lis_start_time, lis_end_tme = CommonUtils.get_time_list(df_meas_time, time_interval)

        if time_interval == '':
            return ""

        lis_time_expression = []
        for i in range(len(lis_start_time)):
            st = lis_start_time[i]
            et = lis_end_tme[i]
            lis_time_expression.append("board_0_timeTag>%d && board_0_timeTag<%d" % (st, et + 1))

        return lis_time_expression

    def _drawTH_with_cut(self, tr, var_name, bins_dic, bin_name, h_name, h_title, cut):

        lis_vars = var_name.split(':')
        dimension = len(lis_vars)

        bin_name = bin_name.lower()

        if dimension == 1:
            bins = bins_dic[bin_name + "_bins"]
            bins_title = bins_dic[bin_name + "_title"]
            h_title += ";%s;counts" % bins_title
            h = self._draw_h1f(tr, var_name, h_name, h_title, bins, cut)
        elif dimension == 2:
            split_bin_name = bin_name.split(':')
            binsy = bins_dic[split_bin_name[0] + "_bins"]
            binsx = bins_dic[split_bin_name[1] + "_bins"]
            binsy_title = bins_dic[split_bin_name[0] + "_title"]
            binsx_title = bins_dic[split_bin_name[1] + "_title"]
            h_title += ";%s;%s;counts" % (binsx_title, binsy_title)
            h = self._draw_h2f(tr, var_name, h_name, h_title, binsx, binsy, cut)
        elif dimension == 3:
            split_bin_name = bin_name.split(':')
            binsz = bins_dic[split_bin_name[0] + "_bins"]
            binsy = bins_dic[split_bin_name[1] + "_bins"]
            binsx = bins_dic[split_bin_name[2] + "_bins"]
            binsz_title = bins_dic[split_bin_name[0] + "_title"]
            binsy_title = bins_dic[split_bin_name[1] + "_title"]
            binsx_title = bins_dic[split_bin_name[2] + "_title"]
            h_title += ";%s;%s;%s;counts" % (binsx_title, binsy_title, binsz_title)
            h = self._draw_h3f(tr, var_name, h_name, h_title, binsx, binsy, binsz, cut)
        else:
            raise Exception("Cannot draw hist with dimension %d" % dimension)

        return h

    def AddSelfDefinedStatistics(self, bins_dic, root_in_fname, root_out_fname,
                                 lis_var_exp, lis_bin_exp, lis_cut_exp, lis_cut_name,
                                 csv_fname='', time_interval=''):
        fin = TFile(root_in_fname, 'read')
        meas_time = fin.Get('meas_time')
        tree = fin.Get('t_tree')

        fout = TFile(root_out_fname, 'recreate')

        n_vars = len(lis_var_exp)
        n_cuts = len(lis_cut_exp)

        for i in range(n_vars):
            bin_name: str = lis_bin_exp[i]
            var_name: str = lis_var_exp[i]

            for j in range(n_cuts):
                cut: str = lis_cut_exp[j]
                if cut.lower() == "null":
                    cut = ""
                cut_name: str = lis_cut_name[j]

                if time_interval == "" or time_interval <0:
                    h_name = 'h{}_{}'.format(var_name.replace(":", "_"), cut_name)
                    h_title = cut_name
                    h = self._drawTH_with_cut(tree, var_name, bins_dic, bin_name, h_name, h_title, cut)
                    h.Write()

                    continue
                else:
                    lis_branch_name = [ib.GetName() for ib in tree.GetListOfBranches()]
                    if "board_0_timeTag" not in lis_branch_name:
                        raise Exception("The TimeInterval was given, but no board_0_timeTag was found")

                lis_time_expression = self._get_lis_time_cut(time_interval, csv_fname)
                for k in range(len(lis_time_expression)):
                    print("\rTimeSlice %d:" % k, end=" ")

                    time_expression = lis_time_expression[k]
                    tr = self._get_tree_in_time(tree, time_expression, k)

                    # get the true start and end time
                    tr.GetEntry(1)
                    start_time = tr.board_0_timeTag
                    tr.GetEntry(tr.GetEntries() - 1)
                    end_time = tr.board_0_timeTag
                    t_title = ("TBegin %s, TEnd %s" %
                               (CommonUtils.convert_timeTag_with_timeZone(start_time),
                                CommonUtils.convert_timeTag_with_timeZone(end_time)))
                    print(" %s" % t_title)

                    time_delta = end_time-start_time

                    h_name = 'h%s_%s_timeSlice_%d_%d' % (var_name.replace(":", "_"), cut_name, start_time, end_time)
                    h_title = '{} {}'.format(cut_name, t_title)

                    # draw hists
                    logger.debug('draw: %s %s' % (h_name, h_title))
                    h = self._drawTH_with_cut(tr, var_name, bins_dic, bin_name, h_name, h_title, cut)
                    h.Write()

                    h_scale = h.Clone()
                    h_scale.Reset("ICES")
                    h_scale.SetName(h.GetName() + "_timeNorm")
                    h_scale.SetTitle(h.GetTitle() + " timeNormalized")
                    nbins = h.GetNbinsX() * h.GetNbinsY() * h.GetNbinsZ()
                    for ibin in range(1, nbins+1):
                        v1 = h.GetBinContent(ibin)

                        if time_delta > 1e-10:
                            h_scale.SetBinContent(ibin, v1 * time_interval/time_delta)
                            h_scale.SetBinError(ibin, TMath.Sqrt(v1) * time_interval/time_delta)

                    h_scale.Write()

        fin.Close()
        fout.Save()
        fout.Close()
