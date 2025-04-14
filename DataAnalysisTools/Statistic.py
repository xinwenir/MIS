##########################################
# Statistic is used for collect all the events in the given tree
# Developed by Guorui Liu
import os.path

from ROOT import TH2F, TH1F, TMath, TTree, TFile, TCanvas, THStack, TLegend, TVectorD
from ROOT import gStyle

from DataFactory import CommonUtils

import logging

logger = logging.getLogger(__name__)


def replace_cut(cut, varname):
    if varname.lower().__contains__('x'):
        cut = cut.replace('_R_', 'x')
    elif varname.lower().__contains__('y'):
        cut = cut.replace('_R_', 'y')

    if varname.__contains__('1'):
        cut = cut.replace('NBar#', 'NBar1')
    elif varname.__contains__('2'):
        cut = cut.replace('NBar#', 'NBar2')

    return cut


def form_h1f_name_title(varname, cut_name, var_unit):
    h_name = 'h{}_{}'.format(varname, cut_name)
    h_title = h_name + '; {} [{}]; counts'.format(varname, var_unit)
    return h_name, h_title


def form_h2f_name_title(varname1, varname2, cut_name, var_unit):
    h_name = 'h{}_{}_{}'.format(varname1, varname2, cut_name)
    h_title = h_name + '; {} [{}]; {} [{}]; counts'.format(varname2, var_unit, varname1, var_unit)
    return h_name, h_title


def hist_divide(h1, h2, h_name):
    if h1 is None or h2 is None:
        logger.fatal('hist_divide, two hists should be passed here.')
        raise Exception(__name__ + ' hist_divide, two hists should be passed here.')
    if not Utils.check_two_hists_the_same_size(h1=h1, h2=h2):
        logger.fatal('hist_divide, passing two hists in different size.')
        raise ValueError(__name__ + ' hist_divide, passing two hists in different size.')

    h3 = h1.Clone()
    h3.Reset("ICES")
    h3.SetNameTitle(h_name, h_name)

    bins = h1.GetNbinsX()*h1.GetNbinsY()*h1.GetNbinsZ()

    for i in range(1, bins+1):
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


def calc_statistic_uncertainty(h):
    bins = h.GetNbinsX()*h.GetNbinsY()*h.GetNbinsZ()
    for i in range(1, bins+1):
        val = h.GetBinContent(i)
        err = 0
        if val > 0:
            err = TMath.Sqrt(val)
        h.SetBinError(i, err)


class Statistic:

    def __init__(self, det, global_dir_tag):
        data_tree_fname = det.get_tree_fname()
        if data_tree_fname == '' or not os.path.isfile(data_tree_fname):
            raise Exception('[Statistic] data_tree_fname dose not exist: %s' % data_tree_fname)
        else:
            logger.info('Read tree from file: ' + data_tree_fname)
            self._tree_file = TFile(data_tree_fname, 'read')
            self._tree = self._tree_file.Get('t_tree')
            self._meas_time = self._tree_file.Get('meas_time')

        self._det = det
        self._h_accept_thetaphi, self._h_accept_thetaxy = self._det.get_acceptance_plane_det()
        self._write_fname = os.path.join(global_dir_tag[0],
                                         global_dir_tag[1] + 'Statistic_' + self._det.get_det_name() + '.root')
        if self._det.get_stats_fname().__contains__('.root'):
            self._write_fname = self._det.get_stats_fname()

        self._print_fname = self._write_fname.replace('.root', '.pdf')

        self._fout = None

        self._h1f = []
        self._h2f = []

        self._can = TCanvas("can_stats", "can_stats", 2000, 1500)

    def _draw_h_stack(self, h_lis):
        hs = THStack("hs", "hs")
        lg = TLegend(0.75, 0.83, 0.95, 0.93)

        for i in range(len(h_lis)):
            h_lis[i].SetLineColor(1+i)
            hs.Add(h_lis[i], "hist e")
            lg.AddEntry(h_lis[i], h_lis[i].GetTitle())

        self._can.cd()
        gStyle.SetOptStat(0)
        hs.Draw("nostack")
        lg.Draw()
        self._can.Print(self._print_fname + '(')

    def _write_print_hist(self, h, dim=''):
        self._fout.cd()
        self._can.cd()

        h.Write()
        h_title = '%s' % (h.GetTitle())
        if dim.__contains__('1'):
            h.Draw("hist e")
        elif dim.__contains__('2'):
            h.Draw("colz")
        else:
            logger.error('Do not know the dimension of this hist: %s, given dimension: %s' % (h_title, dim))
            h.Draw()
        self._can.Print(self._print_fname + '(')

    def _draw_h1f(self, varname, h_name, h_title, bins, cut):
        if bins[1] > bins[2]:
            bins[1] = int(1)
            bins[2] = int(-1)
        h1f = TH1F(h_name, h_title, bins[0], bins[1], bins[2])
        self._tree.Draw("%s>>%s" % (varname, h_name), cut, "goff")
        calc_statistic_uncertainty(h1f)
        return h1f

    def _draw_h2f(self, varname, h_name, h_title, bins1, bins2, cut):
        if bins1[1] > bins1[2]:
            bins1[1] = int(1)
            bins1[2] = int(-1)
        if bins2[1] > bins2[2]:
            bins2[1] = int(1)
            bins1[2] = int(-1)
        h2f = TH2F(h_name, h_title, bins1[0], bins1[1], bins1[2], bins2[0], bins2[1], bins2[2])
        self._tree.Draw("%s>>%s" % (varname, h_name), cut, "goff")
        return h2f

    def write_hists_plane_det(self, bins_dic, h_thetaphi_solid_angle=None, h_thetaxy_solid_angle=None):
        """
        draw all the hists of the muon events
        :param bins_dic:
        :return:
        """
        self._fout = TFile(self._write_fname, "recreate")
        self._meas_time.Write('meas_time')

        if h_thetaphi_solid_angle is None or h_thetaxy_solid_angle is None:
            logger.critical('write_hists: Missing arguments, h_thetaphi_solid_angle, h_thetaxy_solid_angle')
            raise Exception('write_hists: Missing arguments, h_thetaphi_solid_angle, h_thetaxy_solid_angle')

        """define variables and cuts"""
        h1f_varnames_1layer = ['xo1', 'zx1', 'yo1', 'zy1', 'xo2', 'zx2', 'yo2', 'zy2']
        h1f_varnames_2layer = ['X1', 'Y1', 'Z1', 'X2', 'Y2', 'Z2',
                               'theta', 'phi', 'thetax', 'thetay', 'thetax_trk', 'thetay_trk']
        h1f_binnames_1layer = ['x_bins', 'z_bins', 'y_bins', 'z_bins', 'x_bins', 'z_bins', 'y_bins', 'z_bins']
        h1f_binnames_2layer = ['x_bins', 'y_bins', 'z_bins', 'x_bins', 'y_bins', 'z_bins',
                               'theta_bins', 'phi_bins', 'thetax_bins', 'thetay_bins', 'thetax_bins', 'thetay_bins']
        h1f_cuts_1layer = ['',
                           '_R_NBar#==1',
                           '_R_NBar#==2',
                           '_R_NBar#==3',
                           '_R_NBar#==4']
        h1f_cuts_2layer = ['',
                           '_R_NBar1>1 && _R_NBar2>1',
                           '_R_NBar1==1 && _R_NBar2==1',
                           '(_R_NBar1==1 && _R_NBar2>1) || (_R_NBar1>1 && _R_NBar2==1)']
        h1f_cuts_name_1layer = ['TOT', 'S', 'D', 'T', 'Q']  # total, double, single, triple, and quadruple
        h1f_cuts_name_2layer = ['TOT', 'MM', 'SS', 'MS']    # M: multiple

        """draw h1f"""
        logger.info('------------ HIST 1D ------------')

        # draw variables only use one layer information
        n_vars = len(h1f_varnames_1layer)
        n_cuts = len(h1f_cuts_1layer)

        for i in range(n_vars):
            for j in range(n_cuts):
                tmp_cuts = h1f_cuts_1layer[j]
                var_unit = 'mm'
                tmp_bins = bins_dic[h1f_binnames_1layer[i]]

                tmp_cuts = replace_cut(tmp_cuts, h1f_varnames_1layer[i])
                h_name, h_title = form_h1f_name_title(h1f_varnames_1layer[i], h1f_cuts_name_1layer[j], var_unit)
                logger.debug('draw: %s %s' % (h_name, h_title))

                h1f = self._draw_h1f(varname=h1f_varnames_1layer[i], h_name=h_name, h_title=h_title, bins=tmp_bins,
                                     cut=tmp_cuts)
                self._write_print_hist(h1f, '1d')
                self._h1f.append(h1f)
            h_lis = self._h1f[-4:]
            self._draw_h_stack([h_lis[0], h_lis[1], h_lis[2], h_lis[3]])

        # draw variables only use two layer information
        n_vars = len(h1f_varnames_2layer)
        n_cuts = len(h1f_cuts_2layer)

        for i in range(n_vars):
            b_compare = True
            for j in range(n_cuts):
                tmp_cuts = h1f_cuts_2layer[j]
                var_unit = 'mm'
                tmp_bins = bins_dic[h1f_binnames_2layer[i]]

                if h1f_varnames_2layer[i].lower().__contains__('theta') or h1f_varnames_2layer[i].lower().__contains__(
                        'phi'):
                    var_unit = 'rad'
                    if not (h1f_varnames_2layer[i].__contains__('x') or h1f_varnames_2layer[i].__contains__(
                            'y')):
                        tmp_cuts = ''
                        h_name, h_title = form_h1f_name_title(h1f_varnames_2layer[i], 'TOT', var_unit)
                        logger.debug('draw: %s %s' % (h_name, h_title))

                        h1f = self._draw_h1f(varname=h1f_varnames_2layer[i], h_name=h_name, h_title=h_title,
                                             bins=tmp_bins,
                                             cut=tmp_cuts)
                        self._write_print_hist(h1f, '1d')
                        self._h1f.append(h1f)
                        b_compare = False
                        break

                if h1f_varnames_2layer[i].lower().__contains__('z'):
                    tmp_cuts = ''
                    h_name, h_title = form_h1f_name_title(h1f_varnames_2layer[i], 'TOT', var_unit)
                    logger.debug('draw: %s %s' % (h_name, h_title))

                    h1f = self._draw_h1f(varname=h1f_varnames_2layer[i], h_name=h_name, h_title=h_title, bins=tmp_bins,
                                         cut=tmp_cuts)
                    self._write_print_hist(h1f, '1d')
                    self._h1f.append(h1f)
                    b_compare = False
                    break

                tmp_cuts = replace_cut(tmp_cuts, h1f_varnames_2layer[i])
                h_name, h_title = form_h1f_name_title(h1f_varnames_2layer[i], h1f_cuts_name_2layer[j], var_unit)
                logger.debug('draw: %s %s' % (h_name, h_title))

                h1f = self._draw_h1f(varname=h1f_varnames_2layer[i], h_name=h_name, h_title=h_title, bins=tmp_bins,
                                     cut=tmp_cuts)
                self._write_print_hist(h1f, '1d')
                self._h1f.append(h1f)
            if b_compare:
                h_lis = self._h1f[-3:]
                self._draw_h_stack([h_lis[0], h_lis[1], h_lis[2]])

        """draw h2f"""
        logger.info('------------ HIST 2D ------------')
        h2f_varnames = ['X1:Y1', 'X2:Y2', 'theta:phi', 'thetax:thetay', 'thetax_trk:thetay_trk']
        h2f_binnames = ['x_bins:y_bins', 'x_bins:y_bins', 'theta_bins:phi_bins', 'thetax_bins:thetay_bins',
                        'thetax_bins:thetay_bins']
        delim = ':'
        h2f_cuts = ['',
                    'xNBar1>1 && xNBar2>1 && yNBar1>1 && yNBar2>1',
                    'xNBar1==1 && xNBar2==1 && yNBar1==1 && yNBar2==1']
        h2f_cuts_name = ['TOT', 'MMMM', 'SSSS']

        # draw variables only use one layer information
        n_vars = len(h2f_varnames)
        n_cuts = len(h2f_cuts)

        for i in range(n_vars):
            for j in range(n_cuts):
                tmp_cuts = h2f_cuts[j]
                var_unit = 'mm'
                if h2f_varnames[i].__contains__('theta'):
                    var_unit = 'rad'
                var_lis = h2f_varnames[i].split(delim)
                bin_lis = h2f_binnames[i].split(delim)
                var1 = var_lis[0]
                var2 = var_lis[1]
                bins1 = bins_dic[bin_lis[1]]
                bins2 = bins_dic[bin_lis[0]]

                h_name, h_title = form_h2f_name_title(var1, var2, h2f_cuts_name[j], var_unit)
                logger.debug('draw: %s %s' % (h_name, h_title))

                h2f = self._draw_h2f(varname=h2f_varnames[i], h_name=h_name, h_title=h_title,
                                     bins1=bins1, bins2=bins2, cut=h2f_cuts[j])
                self._write_print_hist(h2f, '2d')
                self._h2f.append(h2f)

                if j == 0:
                    h_name_acc = '%s_div_acceptance' % h_name
                    h_name_acc_solid_angle = '%s_div_acceptance_div_solidAngle' % h_name
                    logger.debug('draw: %s %s' % (h_name_acc, h_name_acc_solid_angle))

                    if i == 2:
                        h2f_acc = hist_divide(h1=self._h2f[-1], h2=self._h_accept_thetaphi, h_name=h_name_acc)
                        h2f_acc_solid_angle = hist_divide(h1=h2f_acc, h2=h_thetaphi_solid_angle,
                                                          h_name=h_name_acc_solid_angle)

                    elif i == 3:
                        h2f_acc = hist_divide(h1=self._h2f[-1], h2=self._h_accept_thetaxy, h_name=h_name_acc)
                        h2f_acc_solid_angle = hist_divide(h1=h2f_acc, h2=h_thetaxy_solid_angle,
                                                          h_name=h_name_acc_solid_angle)
                    else:
                        continue

                    self._write_print_hist(h2f_acc, '2d')
                    self._write_print_hist(h2f_acc_solid_angle, '2d')
                    self._h2f.append(h2f)

        can1 = TCanvas("can1", "can1", 2000, 1500)
        can1.Print(self._print_fname + ']')
