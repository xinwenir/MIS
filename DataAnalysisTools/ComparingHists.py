import array
import os

from ROOT import TH2F, TH1F, TMath, TFile, TCanvas, TVector3, THStack, TLegend, TGraphErrors, TMultiGraph
from ROOT import gStyle, gROOT
import ROOT

from DataFactory import CommonUtils, ROOTUtils
from DataAnalysisTools import _hthetaxy_name_prefix, _hthetaphi_name_prefix

import logging

logger = logging.getLogger(__name__)


class ComparingHists:
    def __init__(self):
        pass

    @classmethod
    def _calc_significance(cls, h1, h2, h_name, h_title, **kwargs):

        if_abs = kwargs.get('if_abs', False)
        if_err1_zero = kwargs.get('if_err1_zero', False)
        if_err2_zero = kwargs.get('if_err2_zero', False)

        h_sig = h1.Clone()
        h_sig.Reset("ICES")
        h_sig.SetNameTitle(h_name, h_title)

        xbins = h1.GetNbinsX()
        ybins = h1.GetNbinsY()

        for i_xbin in range(1, xbins + 1):
            for i_ybin in range(1, ybins + 1):
                i_bin = h1.GetBin(i_xbin, i_ybin)

                val1 = h1.GetBinContent(i_bin)
                val2 = h2.GetBinContent(i_bin)

                if val1 > 1e-10 and val2 > 1e-10:
                    err1 = h1.GetBinError(i_bin)
                    err2 = h2.GetBinError(i_bin)
                    if if_err1_zero:
                        err1 = 0.
                    if if_err2_zero:
                        err2 = 0.

                    sig = (val1 - val2) / TMath.Sqrt(err1 * err1 + err2 * err2)
                    if if_abs:
                        sig = abs(sig)

                    h_sig.SetBinContent(i_bin, sig)
                    h_sig.SetBinError(i_bin, 0)
                else:
                    h_sig.SetBinContent(i_bin, 0)
                    h_sig.SetBinError(i_bin, 0)

        return h_sig

    @classmethod
    def _get_whole_keys_in_file(cls, fin_name, lis_keys):
        lis_whole_keys = []
        for i in lis_keys:
            i_name = ROOTUtils.find_key_in_file(fin_name, i)
            if i_name == '':
                raise Exception('Cannot find partial key %s in %s' % (i, fin_name))
            else:
                lis_whole_keys.append(i_name)
        return lis_whole_keys

    @classmethod
    def _split_calcSig_expression(cls, expr: str):
        [a, b] = expr.split('/Err')
        dic_expr = {}
        if a.lower().__contains__('abs'):
            dic_expr['if_abs'] = True
        else:
            dic_expr['if_abs'] = False

        a = a[a.index('(') + 1: a.index(')')]
        [n1, n2] = a.split('-')
        dic_expr['name1'] = n1
        dic_expr['name2'] = n2

        b = b[b.index('(') + 1: b.index(')')]
        lis_err = b.split('&')
        if n1 in lis_err:
            dic_expr['if_err1_zero'] = False
        else:
            dic_expr['if_err1_zero'] = True

        if n2 in lis_err:
            dic_expr['if_err2_zero'] = False
        else:
            dic_expr['if_err2_zero'] = True

        if dic_expr['if_err1_zero'] and dic_expr['if_err2_zero']:
            raise Exception('Input.Expressions: Err() must have at least one valid argument.')

        return dic_expr

    def CalcSignificance(self, expr: str, dic_input: dict, fout_name: str, if_smooth=True):
        dic_expr = self._split_calcSig_expression(expr)

        fin1_name = dic_input.get(dic_expr['name1'] + 'RootName')
        fin2_name = dic_input.get(dic_expr['name2'] + 'RootName')

        lis_fin1_hist_name = dic_input.get(dic_expr['name1'] + 'HistNames')
        lis_fin2_hist_name = dic_input.get(dic_expr['name2'] + 'HistNames')

        if len(lis_fin1_hist_name) != len(lis_fin2_hist_name):
            raise Exception('number of lis_fin1_hist_name and lis_fin2_hist_name must be the same.')

        lis_fin1_hist_name = self._get_whole_keys_in_file(fin1_name, lis_fin1_hist_name)
        lis_fin2_hist_name = self._get_whole_keys_in_file(fin2_name, lis_fin2_hist_name)

        fin1 = TFile(fin1_name, 'read')
        fin2 = TFile(fin2_name, 'read')

        fout = TFile(fout_name, 'update')

        for i in range(len(lis_fin1_hist_name)):
            f1_hist_name = lis_fin1_hist_name[i]
            f2_hist_name = lis_fin2_hist_name[i]

            h1 = fin1.Get(f1_hist_name)
            h2 = fin2.Get(f2_hist_name)

            if f1_hist_name.__contains__('phi'):
                h_name = _hthetaphi_name_prefix + '_sig'
                h_title = h_name + '; phi [rad]; theta [rad]'
            else:
                h_name = _hthetaxy_name_prefix + '_sig'
                h_title = h_name + '; thetay [rad]; thetax [rad]; Sig'

            h_sig = self._calc_significance(h1, h2, h_name, h_title,
                                            if_abs=dic_expr['if_abs'],
                                            if_err1_zero=dic_expr['if_err1_zero'],
                                            if_err2_zero=dic_expr['if_err2_zero'])

            fout.cd()
            h_sig.Write()

            if if_smooth:
                h_sig_smooth = h_sig.Clone()
                h_sig_smooth.SetName('smooth_' + h_sig.GetName())
                h_sig_smooth.SetTitle('smooth_' + h_sig.GetTitle())
                h_sig_smooth.Smooth()
                h_sig_smooth.Write()

        fout.Save()
        fout.Close()

    @staticmethod
    def CutSignificance(fin_name, hist_name, low_thr, high_thr):
        hist_name = ROOTUtils.find_key_in_file(fin_name, hist_name)
        fin = TFile(fin_name, 'update')
        h1 = fin.Get(hist_name)

        h_name = h1.GetName() + '_SigCut%d_%d' % (low_thr, high_thr)
        h_title = h1.GetTitle() + '_SigCut%.1f_%.1f' % (low_thr, high_thr)

        h2 = ROOTUtils.cut_hist_with_threshold(h1, h_name, h_title, low_thr, high_thr)
        h2.Write()

        fin.Save()
        fin.Close()

    def CalcKSTest(self):
        pass

    def CalcQQPlot(self):
        pass

    @classmethod
    def _get_hist_lists(cls, lis_fin_name, lis_hist_name, lis_scale):
        lis_hists = []
        lis_legend = []

        for i in range(len(lis_fin_name)):
            fin_name = lis_fin_name[i]
            hist_name = lis_hist_name[i]
            scale = float(lis_scale[i])

            hist_name = ROOTUtils.find_key_in_file(fin_name, hist_name)

            fin = TFile(fin_name, 'read')
            hist = fin.Get(hist_name)
            hist.SetName(hist.GetName() + '_' + os.path.splitext(os.path.basename(fin_name))[0])
            hist.SetDirectory(gROOT)

            hist.Scale(scale)

            (root_dir, root_fname) = os.path.split(fin_name)
            root_name = os.path.splitext(root_fname)[0]

            lis_hists.append(hist)
            lis_legend.append(hist.GetName() + '_' + root_name)

            fin.Close()

        print('Compare these hists:')
        hs = THStack('hs', ';%s;%s;%s' % (lis_hists[0].GetXaxis().GetTitle(),
                                          lis_hists[0].GetYaxis().GetTitle(),
                                          lis_hists[0].GetZaxis().GetTitle()))
        for i in range(len(lis_hists)):
            lg = lis_legend[i]
            print(lg)

            h = lis_hists[i]
            h.SetLineColor(i + 1)
            h.SetMarkerColor(i + 1)
            hs.Add(h, 'e')

        return hs, lis_legend

    @classmethod
    def _integral1D_hstack(cls, hs_in, lis_lg: list):
        xbin = hs_in.GetHists()[0].GetNbinsX()

        hs_int = THStack('hs_int', ';%s;%s'
                         % (hs_in.GetHists()[0].GetXaxis().GetTitle(), hs_in.GetHists()[0].GetYaxis().GetTitle()))
        lg = TLegend(0.1, 0.9, 0.9, 0.99)
        for ih, i_lg in zip(hs_in.GetHists(), lis_lg):
            _sum = 0
            _err = 0
            h1 = ih.Clone()
            h1.Reset('ICES')
            h1.SetName('%s_integral' % ih.GetName())
            for ix in range(1, xbin + 1):
                _sum += ih.GetBinContent(ix)
                _err += ih.GetBinError(ix)
                h1.SetBinContent(ix, _sum)
                h1.SetBinError(ix, _err)

            if h1.GetBinContent(xbin) > 0:
                h1.Scale(1. / h1.GetBinContent(xbin))
            hs_int.Add(h1, 'e')

            ks_test = ih.KolmogorovTest(hs_in.GetHists()[0])
            i_lg += '  KSTest:' + str(ks_test)

            lg.AddEntry(h1, i_lg)

        return hs_int, lg

    @classmethod
    def _cross_plot_hstack(cls, hs_in, lis_lg: list):
        xbin = hs_in.GetHists()[0].GetNbinsX()

        mg_cross = TMultiGraph('mg_cross', '; first hist')
        lg = TLegend(0.1, 0.9, 0.9, 0.99)
        i = 0
        for ih, i_lg in zip(hs_in.GetHists(), lis_lg):
            arr_x = array.array('d', xbin * [0.])
            arr_y = array.array('d', xbin * [0.])
            arr_x_err = array.array('d', xbin * [0.])
            arr_y_err = array.array('d', xbin * [0.])

            for ix in range(1, xbin + 1):
                arr_x[ix - 1] = hs_in.GetHists()[0].GetBinContent(ix)
                arr_y[ix - 1] = ih.GetBinContent(ix)
                arr_x_err[ix - 1] = hs_in.GetHists()[0].GetBinError(ix)
                arr_y_err[ix - 1] = ih.GetBinError(ix)

            gr = TGraphErrors(xbin, arr_x, arr_y, arr_x_err, arr_y_err)
            gr.SetName('gr_%s' % ih.GetName())
            gr.SetLineColor(i + 1)
            gr.SetMarkerColor(i + 1)

            mg_cross.Add(gr, 'P')
            lg.AddEntry(gr, i_lg)

            i += 1

        return mg_cross, lg

    @classmethod
    def get_hist_ratio(cls, h1, h2):
        h_sig = h1.Clone()
        h_sig.Reset("ICES")
        h_sig.SetNameTitle(h1.GetName() + '_sig', '')

        xbins = h1.GetNbinsX()
        ybins = h1.GetNbinsY()

        for i_xbin in range(1, xbins + 1):
            for i_ybin in range(1, ybins + 1):
                i_bin = h1.GetBin(i_xbin, i_ybin)

                val1 = h1.GetBinContent(i_bin)
                val2 = h2.GetBinContent(i_bin)

                if val1 > 1e-10 and val2 > 1e-10:
                    err1 = h1.GetBinError(i_bin)
                    err2 = h2.GetBinError(i_bin)

                    sig = val1 / val2
                    err = TMath.Sqrt(err1 * err1 / val2 / val2 + val1 * val1 / val2 / val2 / val2 / val2 * err2 * err2)

                    h_sig.SetBinContent(i_bin, sig)
                    h_sig.SetBinError(i_bin, err)
                else:
                    h_sig.SetBinContent(i_bin, 1)
                    h_sig.SetBinError(i_bin, 0)

        return h_sig

    def _ratio_hstack(self, hs_in, lis_lg: list):
        href = hs_in.GetHists()[0]

        hs_int = THStack('hs_int', ';%s;%s'
                         % (hs_in.GetHists()[0].GetXaxis().GetTitle(), hs_in.GetHists()[0].GetYaxis().GetTitle()))
        lg = TLegend(0.1, 0.9, 0.9, 0.99)
        for ih, i_lg in zip(hs_in.GetHists(), lis_lg):
            _sum = 0
            _err = 0
            h1 = self.get_hist_ratio(ih, href)

            hs_int.Add(h1, 'e')
            i_lg += '  Ratio:' + str(ih.GetName())

            lg.AddEntry(h1, i_lg)

        return hs_int, lg

    @classmethod
    def _integral2D_hstack(cls, hs_in, lis_lg: list, x_or_y: str):
        hs_title = ';'
        h_name = 'hint_'
        if x_or_y.lower().__contains__('x'):
            hs_title += hs_in.GetHists()[0].GetXaxis().GetTitle()
            h_name += 'x'
            xbin = hs_in.GetHists()[0].GetNbinsX()
        else:
            hs_title += hs_in.GetHists()[0].GetYaxis().GetTitle()
            h_name += 'y'
            xbin = hs_in.GetHists()[0].GetNbinsY()

        hs_int = THStack(h_name, hs_title)
        lg = TLegend(0.1, 0.9, 0.9, 0.99)
        lg.SetName('lg_' + h_name)

        i = 0
        for ih, i_lg in zip(hs_in.GetHists(), lis_lg):
            if x_or_y.lower().__contains__('x'):
                h1 = ih.ProjectionX().Clone()
            else:
                h1 = ih.ProjectionY().Clone()
            h1.SetName(h_name + '_' + ih.GetName())
            h1.Reset("ICES")
            for ix in range(1, xbin + 1):
                err = array.array('d', [0.])
                if x_or_y.lower().__contains__('x'):
                    val_int = ih.IntegralAndError(ix, ix, 1, -1, err)
                else:
                    val_int = ih.IntegralAndError(1, -1, ix, ix, err)

                h1.SetBinContent(ix, val_int)
                h1.SetBinError(ix, err[0])

            h1.SetLineColor(1 + i)
            h1.SetMarkerColor(1 + i)

            hs_int.Add(h1, 'e')
            lg.AddEntry(h1, h1.GetName())
            i += 1

        return hs_int, lg

    def _integral_comparing(self, fout_name, hs, lis_lg, x_or_y: str):
        if x_or_y.lower().__contains__('x'):
            x_or_y = 'x'
        else:
            x_or_y = 'y'

        fout = TFile(fout_name, 'update')

        hs_int, lg = self._integral2D_hstack(hs, lis_lg, x_or_y)
        hs_int.Write()
        lg.Write()

        if x_or_y == 'x':
            hs_name = 'hs_sig_integral_x'
            hs_title = ';%s;significance' % hs.GetHists()[0].GetXaxis().GetTitle()
        else:
            hs_name = 'hs_sig_integral_y'
            hs_title = ';%s;significance' % hs.GetHists()[0].GetYaxis().GetTitle()

        hs_sig = THStack(hs_name, hs_title)
        lg_sig = TLegend(0.1, 0.9, 0.9, 0.99)
        lg_sig.SetName('lg_' + hs_sig.GetName())

        for i in range(len(lis_lg)):
            h1 = hs_int.GetHists()[0]
            h2 = hs_int.GetHists()[i]
            h_name = '%s_vs_%s' % (lis_lg[0], lis_lg[i])
            h_title = ';%s;Sig' % (h1.GetXaxis().GetTitle())
            hsig = self._calc_significance(h1, h2, h_name, h_title)
            hsig.SetLineColor(i + 1)
            hsig.SetMarkerColor(i + 1)

            hs_sig.Add(hsig)

        can = TCanvas('can_hist_list_integral_' + x_or_y, 'can_hist_list_integral_' + x_or_y, 2400, 1800)
        can.Divide(2, 2)

        can.cd(1)
        hs_int.Draw('nostack')
        lg.Draw()

        # significance
        can.cd(2)
        hs_sig.Draw('nostack')
        lg_sig.Draw()

        # integral
        can.cd(3)
        hs1d_integral, lg_integral = self._integral1D_hstack(hs_int, lis_lg=lis_lg)
        hs1d_integral.Draw('nostack')
        lg_integral.Draw()

        # cross plot
        can.cd(4)
        # hs1d_cross, lg_icross = self._cross_plot_hstack(hs_int, lis_lg=lis_lg)
        # hs1d_cross.Draw('A')
        # lg_icross.Draw()
        hs1d_ratio, lg_ratio = self._ratio_hstack(hs_int, lis_lg)
        hs1d_ratio.Draw('nostack')
        lg_ratio.Draw()

        can.Write()

        fout.Save()
        fout.Close()

    def CompareHistList(self, lis_fin_name, lis_hist_name, lis_scale, fout_name: str):
        if len(lis_fin_name) != len(lis_hist_name) or len(lis_fin_name) != len(lis_scale):
            raise Exception('number of lis_fin_name, lis_hist_name and lis_scale'
                            'must be the same. ')
        if len(lis_fin_name) < 2:
            raise Exception('CompareHistList: need at least two comparing hists.')

        hs, lis_lg = self._get_hist_lists(lis_fin_name, lis_hist_name, lis_scale)

        self._integral_comparing(fout_name, hs, lis_lg, 'x')
        self._integral_comparing(fout_name, hs, lis_lg, 'y')

        # significance
        hs_sig = THStack('hs_sig', ';%s;%s;significance' %
                         (hs.GetHists()[0].GetXaxis().GetTitle(), hs.GetHists()[0].GetYaxis().GetTitle()))
        for i in range(len(lis_lg)):
            h1 = hs.GetHists()[0]
            h2 = hs.GetHists()[i]
            h_name = '%s_vs_%s' % (lis_lg[0], lis_lg[i])
            h_title = ';%s;%s;Sig' % (h1.GetXaxis().GetTitle(), h1.GetYaxis().GetTitle())
            hsig = self._calc_significance(h1, h2, h_name, h_title)
            hsig.SetLineColor(i + 1)
            hsig.SetMarkerColor(i + 1)

            hs_sig.Add(hsig)

        fout = TFile(fout_name, 'update')

        can = TCanvas('can_hist_list', 'can_hist_list', 2400, 1800)
        can.Divide(2, 2)

        xbin = hs.GetHists()[0].GetNbinsX()
        ybin = hs.GetHists()[0].GetNbinsY()

        for i in range(1, xbin + 1):
            if hs.GetHists()[0].Integral(i, i, 1, -1) < 1:
                continue
            can.SetName('can_slice2d_x%d' % i)
            can.cd(1)
            hs1d, lg1d = self._get_hstack_1D(hs, i, 'x', lis_lg=lis_lg)
            hs1d.Draw('nostack')
            lg1d.Draw()

            # significance
            can.cd(2)
            hs1d_sig, lg1d_sig = self._get_hstack_1D(hs_sig, i, 'x', lis_lg=lis_lg)
            hs1d_sig.Draw('nostack')
            lg1d_sig.Draw()

            # integral
            can.cd(3)
            hs1d_integral, lg_integral = self._integral1D_hstack(hs1d, lis_lg=lis_lg)
            hs1d_integral.Draw('nostack')
            lg_integral.Draw()

            # cross plot
            can.cd(4)
            # hs1d_cross, lg_icross = self._cross_plot_hstack(hs1d, lis_lg=lis_lg)
            # hs1d_cross.Draw('A')
            # lis_lg.Draw()
            hs1d_ratio, lg_ratio = self._ratio_hstack(hs1d, lis_lg)
            hs1d_ratio.Draw('nostack')
            lg_ratio.Draw()

            can.Write()

        for i in range(1, ybin + 1):
            if hs.GetHists()[0].Integral(1, -1, i, i) < 1:
                continue
            can.SetName('can_slice2d_y%d' % i)
            can.cd(1)
            hs1d, lg1d = self._get_hstack_1D(hs, i, 'y', lis_lg=lis_lg)
            hs1d.Draw('nostack')
            lg1d.Draw()

            # significance
            can.cd(2)
            hs1d_sig, lg1d_sig = self._get_hstack_1D(hs_sig, i, 'y', lis_lg=lis_lg)
            hs1d_sig.Draw('nostack')
            lg1d_sig.Draw()

            # integral
            can.cd(3)
            hs1d_integral, lg_integral = self._integral1D_hstack(hs1d, lis_lg=lis_lg)
            hs1d_integral.Draw('nostack')
            lg_integral.Draw()

            # cross plot
            can.cd(4)
            # hs1d_cross, lg_icross = self._cross_plot_hstack(hs1d, lis_lg=lis_lg)
            # hs1d_cross.Draw('A')
            # lg_icross.Draw()
            hs1d_ratio, lg_ratio = self._ratio_hstack(hs1d, lis_lg)
            hs1d_ratio.Draw('nostack')
            lg_ratio.Draw()

            can.Write()

        fout.Save()
        fout.Close()

    @classmethod
    def _get_hist_with_prefix(cls, fin_name, h_prefix, lis_key_parts=None):
        if lis_key_parts is None:
            h_name = ROOTUtils.find_key_in_file(fin_name, h_prefix)
        else:
            h_name = ROOTUtils.find_key_in_file(fin_name, h_prefix, lis_key_parts=lis_key_parts)
            # if h_name == '':
            #     h_name = Utils.find_key_in_file(fin_name, h_prefix)
        if h_name == '':
            raise Exception('Cannot find %s in %s' % (h_prefix, fin_name))

        fin = TFile(fin_name, 'read')
        h = fin.Get(h_name)
        h.SetDirectory(gROOT)

        fin.Close()
        return h

    @classmethod
    def _check_hstack_same_size(cls, hs):
        for i in range(1, len(hs.GetHists())):
            if not ROOTUtils.check_two_hists_the_same_size(hs.GetHists()[0], hs.GetHists()[i]):
                raise Exception('%s and %s are not the same size.'
                                % (hs.GetHists()[0].GetName(), hs.GetHists()[i].GetName()))
        return True

    @classmethod
    def _get_hstack_1D(cls, hs_in, i_bin: int, x_or_y: str, lis_lg: list = None, first_no_option=True):
        lg = TLegend(0.1, 0.9, 0.9, 0.99)
        if x_or_y.lower().__contains__('x'):
            x_title = hs_in.GetHists()[0].GetYaxis().GetTitle()
            hs_name = hs_in.GetName() + '_%s%d' % (x_title.split()[0], i_bin)
            lg.SetHeader('%s = %.3f' % (hs_in.GetHists()[0].GetXaxis().GetTitle(),
                                        hs_in.GetHists()[0].GetXaxis().GetBinCenter(i_bin)))
        else:
            x_title = hs_in.GetHists()[0].GetXaxis().GetTitle()
            hs_name = hs_in.GetName() + '_%s%d' % (x_title.split()[0], i_bin)
            lg.SetHeader('%s = %.3f' % (hs_in.GetHists()[0].GetYaxis().GetTitle(),
                                        hs_in.GetHists()[0].GetYaxis().GetBinCenter(i_bin)))
        hs_title = '; %s; %s' % (x_title, hs_in.GetHists()[0].GetZaxis().GetTitle())
        hs = THStack(hs_name, hs_title)
        lg.SetName('lg_%s' % hs_name)

        i = 0
        for ih in hs_in.GetHists():
            if x_or_y.lower().__contains__('x'):
                h1d = ih.ProjectionY(ih.GetName() + '_%s%d' % (x_or_y, i_bin),
                                     i_bin, i_bin)
            else:
                h1d = ih.ProjectionX(ih.GetName() + '_%s%d' % (x_or_y, i_bin),
                                     i_bin, i_bin)
            h1d.SetLineColor(1 + i)
            h1d.SetMarkerColor(1 + i)
            if i == 0 and first_no_option:
                hs.Add(h1d, 'e')
            else:
                hs.Add(h1d, 'hist e')
            if lis_lg is not None and len(hs.GetHists()) == len(lis_lg):
                lg.AddEntry(h1d, lis_lg[i], 'l')
            else:
                lg.AddEntry(h1d, h1d.GetName(), 'l')
            i += 1

        return hs, lg

    def _get_target_background_hstack(self, dic_input, coord_name):
        meas_target_stats_fname = dic_input.get('MeasTargetStatsName')
        meas_background_stats_fname = dic_input.get('MeasBackgroundStatsName')
        h_cut_name = dic_input.get('StatsCutName', 'TOT')

        # find the hist name
        if coord_name.__contains__('phi'):
            h_prefix = _hthetaphi_name_prefix
        else:
            h_prefix = _hthetaxy_name_prefix
        h_meas_target_stats = self._get_hist_with_prefix(meas_target_stats_fname, h_prefix, [h_cut_name])
        fname = os.path.splitext(os.path.split(meas_target_stats_fname)[1])[0]
        h_meas_target_stats.SetName(h_meas_target_stats.GetName() + '_%s' % fname)
        h_meas_background_stats = self._get_hist_with_prefix(meas_background_stats_fname, h_prefix, [h_cut_name])
        fname = os.path.splitext(os.path.split(meas_background_stats_fname)[1])[0]
        h_meas_background_stats.SetName(h_meas_background_stats.GetName() + '_%s' % fname)

        hs = THStack('hs_target_background',
                     ';%s;%s' % (h_meas_target_stats.GetXaxis().GetTitle(), h_meas_target_stats.GetYaxis().GetTitle()))

        hs.Add(h_meas_target_stats)
        hs.Add(h_meas_background_stats)

        lg = TLegend(0.1, 0.9, 0.9, 0.99)
        lg.SetName('lg_hs_target_background')
        lg.SetFillStyle(0)

        lis_lg = []
        lg.AddEntry(h_meas_target_stats, h_meas_target_stats.GetName())
        lis_lg.append(h_meas_target_stats.GetName())
        lg.AddEntry(h_meas_background_stats, h_meas_background_stats.GetName())
        lis_lg.append(h_meas_background_stats.GetName())

        return hs, lg, lis_lg

    def _slice_and_write_in_root(self, hs, lg, fout_name):
        fout = TFile(fout_name, 'update')

        hs.Write()
        lg.Write()

        xbin = hs.GetHists()[0].GetNbinsX()
        ybin = hs.GetHists()[0].GetNbinsY()

        for i in range(1, xbin + 1):
            if hs.GetHists()[0].Integral(i, i, 1, -1) > 1e-15:
                hs1d, lg1d = self._get_hstack_1D(hs, i, 'x')
                hs1d.Write()
                lg1d.Write()
        for i in range(1, ybin + 1):
            if hs.GetHists()[0].Integral(1, -1, i, i) > 1e-15:
                hs1d, lg1d = self._get_hstack_1D(hs, i, 'y')
                hs1d.Write()
                lg1d.Write()

        fout.Save()
        fout.Close()

    def _get_special_hstack(self, dic_input, coord_name, special_name: str):
        if special_name.__contains__('ratio'):
            special_name = 'ratio'
            config_name = 'Ratio'
        elif special_name.lower().__contains__('leff'):
            special_name = 'Leff'
            config_name = 'Leff'
        elif special_name.__contains__('sig'):
            special_name = 'sig'
            config_name = 'Sig'
        else:
            raise Exception('only three special_name are allowed: ratio, leff and sig')

        if special_name != 'sig':
            meas_spc_fname = dic_input.get('Meas%sRootName' % config_name)
            lis_sim_spc_fname = dic_input.get('Sim%sRootNames' % config_name)
            if not isinstance(lis_sim_spc_fname, list):
                raise Exception('Input.Sim%sRootNames must be a list.' % config_name)
        else:
            meas_spc_fname = None
            lis_sim_spc_fname = dic_input.get('SigRootNames')
            if not isinstance(lis_sim_spc_fname, list):
                raise Exception('Input.SigRootNames must be a list.')
        h_cut_name = dic_input.get('StatsCutName', 'TOT')

        # find the hist name
        if coord_name.__contains__('phi'):
            h_prefix = _hthetaphi_name_prefix
        else:
            h_prefix = _hthetaxy_name_prefix

        lis_h = []
        for sim_fname in lis_sim_spc_fname:
            h_sim = self._get_hist_with_prefix(sim_fname, h_prefix, [special_name, h_cut_name])
            fname = os.path.splitext(os.path.split(sim_fname)[1])[0]
            h_sim.SetName(h_sim.GetName() + '_%s' % fname)
            lis_h.append(h_sim)

        if meas_spc_fname is not None:
            h_meas_spc = self._get_hist_with_prefix(meas_spc_fname, h_prefix, [special_name, h_cut_name])
            fname = os.path.splitext(os.path.split(meas_spc_fname)[1])[0]
            h_meas_spc.SetName(h_meas_spc.GetName() + '_%s' % fname)
        else:
            h_meas_spc = lis_h[0].Clone()
            h_meas_spc.Reset("ICES")
            h_meas_spc.SetName('h_useless')

        # add into a hstack
        hs = THStack('hs_meas_sim_%s' % special_name,
                     ';%s;%s;%s' % (h_meas_spc.GetXaxis().GetTitle(),
                                    h_meas_spc.GetYaxis().GetTitle(),
                                    h_meas_spc.GetZaxis().GetTitle()))

        if meas_spc_fname is not None:
            hs.Add(h_meas_spc)

        for ih in lis_h:
            hs.Add(ih)

        lg = TLegend(0.1, 0.9, 0.9, 0.99)
        lg.SetName('lg_hs_meas_sim_%s' % special_name)
        lg.SetFillStyle(0)

        lis_lg = []
        if meas_spc_fname is not None:
            lg.AddEntry(h_meas_spc, h_meas_spc.GetName())
            lis_lg.append(h_meas_spc.GetName())

        for ih in lis_h:
            lg.AddEntry(ih, ih.GetName())
            lis_lg.append(ih.GetName())

        return hs, lg, lis_lg

    def TargetBackgroundSlice(self, dic_input, fout_name, coord_name: str):
        hs, lg, lis_lg = self._get_target_background_hstack(dic_input, coord_name)
        self._check_hstack_same_size(hs)

        self._slice_and_write_in_root(hs, lg, fout_name)

    def MeasSimLeffSlice(self, dic_input, fout_name, coord_name: str):
        hs, lg, lis_lg = self._get_special_hstack(dic_input, coord_name, 'Leff')
        self._check_hstack_same_size(hs)

        self._slice_and_write_in_root(hs, lg, fout_name)

    def MeasSimRatioSlice(self, dic_input, fout_name, coord_name):
        hs, lg, lis_lg = self._get_special_hstack(dic_input, coord_name, 'ratio')
        self._check_hstack_same_size(hs)

        self._slice_and_write_in_root(hs, lg, fout_name)

    def MeasSimSigSlice(self, dic_input, fout_name, coord_name: str):
        hs, lg, lis_lg = self._get_special_hstack(dic_input, coord_name, 'sig')
        self._check_hstack_same_size(hs)

        self._slice_and_write_in_root(hs, lg, fout_name)

    def _check_four_hstack_size(self, hs1, hs2, hs3, hs4):
        hs_check = THStack()
        hs_check.Add(hs1.GetHists()[0])
        hs_check.Add(hs2.GetHists()[0])
        hs_check.Add(hs3.GetHists()[0])
        hs_check.Add(hs4.GetHists()[0])
        self._check_hstack_same_size(hs_check)

    @classmethod
    def _adjust_draw_range(cls, h):
        h_min = h.GetMinimum()
        h_max = h.GetMaximum()

        h.SetMaximum(2 * h_max - h_min)
        h.SetMinimum(2 * h_min - h_max)

    @classmethod
    def _set_h_min_max(cls, h, h_min, h_max):
        if h_max is not None:
            h.SetMaximum(h_max)
        if h_min is not None:
            h.SetMinimum(h_min)

    def Slice2D(self, dic_input, fout_name, coord_name='thetaphi', adjust_draw_range=False):
        self.TargetBackgroundSlice(dic_input, fout_name, coord_name)
        self.MeasSimLeffSlice(dic_input, fout_name, coord_name)
        self.MeasSimRatioSlice(dic_input, fout_name, coord_name)
        self.MeasSimSigSlice(dic_input, fout_name, coord_name)

        hs1, lg1, lis_lg1 = self._get_target_background_hstack(dic_input, coord_name)
        self._check_hstack_same_size(hs1)

        hs2, lg2, lis_lg2 = self._get_special_hstack(dic_input, coord_name, 'Leff')
        self._check_hstack_same_size(hs2)

        hs3, lg3, lis_lg3 = self._get_special_hstack(dic_input, coord_name, 'ratio')
        self._check_hstack_same_size(hs3)

        hs4, lg4, lis_lg4 = self._get_special_hstack(dic_input, coord_name, 'sig')
        self._check_hstack_same_size(hs4)

        # check four h2d same size
        self._check_four_hstack_size(hs1, hs2, hs3, hs4)

        meas_ratio_min = None
        try:
            meas_ratio_min = float(dic_input.get('MeasRatioSetMinMax')[0])
        except:
            pass
        meas_ratio_max = None
        try:
            meas_ratio_max = float(dic_input.get('MeasRatioSetMinMax')[1])
        except:
            pass
        meas_leff_min = None
        try:
            meas_leff_min = float(dic_input.get('MeasLeffSetMinMax')[0])
        except:
            pass
        meas_leff_max = None
        try:
            meas_leff_max = float(dic_input.get('MeasLeffSetMinMax')[1])
        except:
            pass

        fout = TFile(fout_name, 'update')
        can = TCanvas('can_slice2d', 'can_slice2d', 2400, 1800)
        can.Divide(2, 2)

        xbin = hs1.GetHists()[0].GetNbinsX()
        ybin = hs1.GetHists()[0].GetNbinsY()

        for i in range(1, xbin + 1):
            if hs1.GetHists()[0].Integral(i, i, 1, -1) > 1:
                hs1d_1, lg1d_1 = self._get_hstack_1D(hs1, i, 'x')
                hs1d_2, lg1d_2 = self._get_hstack_1D(hs2, i, 'x')
                hs1d_3, lg1d_3 = self._get_hstack_1D(hs3, i, 'x')
                hs1d_4, lg1d_4 = self._get_hstack_1D(hs4, i, 'x')

                can.SetName('can_slice2d_%s_x%d' % (coord_name, i))
                can.cd(1)
                hs1d_1.Draw('nostack')
                lg1d_1.Draw()
                can.cd(2)
                if adjust_draw_range:
                    self._adjust_draw_range(hs1d_2)
                self._set_h_min_max(hs1d_2, meas_leff_min, meas_leff_max)
                hs1d_2.Draw('nostack')
                lg1d_2.Draw()
                can.cd(3)
                if adjust_draw_range:
                    self._adjust_draw_range(hs1d_3)
                self._set_h_min_max(hs1d_3, meas_ratio_min, meas_ratio_max)
                hs1d_3.Draw('nostack')
                lg1d_3.Draw()
                can.cd(4)
                hs1d_4.Draw('nostack')
                lg1d_4.Draw()

                can.Write()

        for i in range(1, ybin + 1):
            if hs1.GetHists()[0].Integral(1, -1, i, i) > 1:
                hs1d_1, lg1d_1 = self._get_hstack_1D(hs1, i, 'y')
                hs1d_2, lg1d_2 = self._get_hstack_1D(hs2, i, 'y')
                hs1d_3, lg1d_3 = self._get_hstack_1D(hs3, i, 'y')
                hs1d_4, lg1d_4 = self._get_hstack_1D(hs4, i, 'y')

                can.SetName('can_slice2d_%s_y%d' % (coord_name, i))
                can.cd(1)
                hs1d_1.Draw('nostack')
                lg1d_1.Draw()
                can.cd(2)
                if adjust_draw_range:
                    self._adjust_draw_range(hs1d_2)
                self._set_h_min_max(hs1d_2, meas_leff_min, meas_leff_max)
                hs1d_2.Draw('nostack')
                lg1d_2.Draw()
                can.cd(3)
                if adjust_draw_range:
                    self._adjust_draw_range(hs1d_3)
                self._set_h_min_max(hs1d_3, meas_ratio_min, meas_ratio_max)
                hs1d_3.Draw('nostack')
                lg1d_3.Draw()
                can.cd(4)
                hs1d_4.Draw('nostack')
                lg1d_4.Draw()

                can.Write()

        fout.Save()
        fout.Close()

    def SliceOneTH2(self, fin_name, hist_name, fout_name):
        hist_name = ROOTUtils.find_key_in_file(fin_name, hist_name)
        fin = TFile(fin_name, 'read')
        hist = fin.Get(hist_name)

        hs = THStack('hs_%s' % hist_name, ';%s;%s' % (hist.GetXaxis().GetTitle(), hist.GetYaxis().GetTitle()))
        hs.Add(hist)
        lg = TLegend(0.75, 0.85, 0.9, 0.9)
        lg.SetName('lg_hs_%s' % hist_name)
        lg.AddEntry(hist, hist.GetName(), 'l')
        self._slice_and_write_in_root(hs, lg, fout_name)

    @classmethod
    def PrintSpectrumSlides(cls, fin_name, lis_slide_number, fout_name, if_smooth=False):
        if not isinstance(lis_slide_number[0], int):
            lis_slide_number = [int(i) for i in lis_slide_number]

        dic_hname = ROOTUtils.classify_by_board_channel_ID(fin_name)

        fin = TFile(fin_name, 'read')
        dic_spec = {}
        for ikey in dic_hname:
            lis_hname = dic_hname[ikey]
            lis_h_draw = []
            for ihname in lis_hname:
                for i_num in lis_slide_number:
                    islide = ihname.split('_')[-1]
                    if int(islide) == i_num:
                        h = fin.Get(ihname)
                        lis_h_draw.append(h)

            if len(lis_h_draw) == 0:
                continue
            dic_spec[ikey] = lis_h_draw

        fout = TFile(fout_name, 'recreate')

        for ikey in dic_spec:
            lis_h = dic_spec[ikey]
            hs = THStack('hs_%s' % ikey,
                         'slides %s;%s;%s' % (lis_slide_number.__str__(), lis_h[0].GetXaxis().GetTitle(), lis_h[0].GetYaxis().GetTitle()))
            lg = TLegend(0.3, 0.8, 0.9, 0.9)
            lg.SetName('lg_hs_%s' % ikey)

            for i in range(len(lis_h)):
                ih = lis_h[i]
                area = ih.Integral()
                meanx = ih.GetMean(1)
                if if_smooth:
                    ih.Smooth()
                    ih.SetName('smooth_'+ih.GetName())
                ih.SetLineColor(1 + i)
                ih.SetMarkerColor(1 + i)
                hs.Add(ih)
                lg.AddEntry(ih,
                            ih.GetName() + '; Area: %.1f' % area
                            + '; MeanX: %.1f' % meanx, 'l')

            hs.Write()
            lg.Write()

        fout.Save()
        fout.Close()
