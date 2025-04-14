import os
import time

from ROOT import TH2F, TH1F, TMath, TFile, TCanvas, TVector3
from ROOT import gStyle, gROOT
from DataAnalysisTools import _hthetaxy_name_prefix, _hthetaphi_name_prefix

from DataFactory import CommonUtils, ROOTUtils

import logging

logger = logging.getLogger(__name__)


class GetRatioLeff:

    def __init__(self):
        pass

    @classmethod
    def _calc_hist_ratio(cls, h1, h2, time_ratio, h_name, h_title):
        h_ratio = h1.Clone()
        h_ratio.Reset("ICES")
        h_ratio.SetNameTitle(h_name, h_title)

        xbins = h1.GetNbinsX()
        ybins = h1.GetNbinsY()

        for i_xbin in range(1, xbins + 1):
            for i_ybin in range(1, ybins + 1):
                i_bin = h1.GetBin(i_xbin, i_ybin)

                val1 = h1.GetBinContent(i_bin)
                val2 = h2.GetBinContent(i_bin)

                if val1 > 1e-10 and val2 > 1e-10:
                    _ratio = (val1 / val2) / time_ratio
                    err1 = TMath.Sqrt(val1)
                    err2 = TMath.Sqrt(val2)
                    _ratio_err = TMath.Sqrt(
                        err1 * err1 / val2 / val2 + pow((val1 / val2 / val2 * err2), 2.0)) / time_ratio

                    h_ratio.SetBinContent(i_bin, _ratio)
                    h_ratio.SetBinError(i_bin, _ratio_err)
                else:
                    h_ratio.SetBinContent(i_bin, 0)
                    h_ratio.SetBinError(i_bin, 0)

        return h_ratio

    def GetMeasRatio(self, tg_stats_fname, bg_stats_fname, ratio_fname,
                     ratio_mode=None, ratio_scale=1., log_fname=None):
        if ratio_mode is None:
            ratio_mode = ['TOT']
        if not os.path.isfile(tg_stats_fname) or not os.path.isfile(bg_stats_fname):
            raise Exception('Missing stats files for getting ratio.')

        f_tg = TFile(tg_stats_fname, 'read')
        f_bg = TFile(bg_stats_fname, 'read')

        meas_time_tg = f_tg.Get('meas_time')[0]
        meas_time_bg = f_bg.Get('meas_time')[0]
        time_ratio = float(ratio_scale) * meas_time_tg / meas_time_bg

        fout = TFile(ratio_fname, 'update')

        for i_mode in ratio_mode:
            hthetaxy_name = _hthetaxy_name_prefix + '_%s' % i_mode
            hthetaphi_name = _hthetaphi_name_prefix + '_%s' % i_mode

            hthetaxy_tg = f_tg.Get(hthetaxy_name)
            hthetaxy_tg.SetDirectory(gROOT)
            hthetaxy_bg = f_bg.Get(hthetaxy_name)
            hthetaxy_bg.SetDirectory(gROOT)

            logger.debug('Get %s,%s from %s and %s' % (hthetaxy_name, hthetaphi_name, tg_stats_fname, bg_stats_fname))

            if not ROOTUtils.check_two_hists_the_same_size(hthetaxy_tg, hthetaxy_bg):
                raise Exception('hthetax_thetay of target and of background have '
                                'different size. Please check input files')

            logger.debug('Begin to calculate ratio thetax_thetay: %s' % ratio_fname)
            h_name = _hthetaxy_name_prefix + '_ratio'
            h_title = h_name + '; thetay [rad]; thetax [rad]; ratio'
            h_ratio_thetaxy = self._calc_hist_ratio(hthetaxy_tg, hthetaxy_bg, time_ratio, h_name, h_title)

            h_ratio_thetaxy.Write()

            hthetaphi_tg = f_tg.Get(hthetaphi_name)
            hthetaphi_tg.SetDirectory(gROOT)
            hthetaphi_bg = f_bg.Get(hthetaphi_name)
            hthetaphi_bg.SetDirectory(gROOT)

            if not ROOTUtils.check_two_hists_the_same_size(hthetaphi_tg, hthetaphi_bg):
                raise Exception('htheta_phi of target and of background have '
                                'different size. Please check input files')

            logger.debug('Begin to calculate ratio theta_phi: %s' % ratio_fname)
            h_name = _hthetaphi_name_prefix + '_ratio'
            h_title = h_name + '; phi [rad]; theta [rad]; ratio'
            h_ratio_thetaphi = self._calc_hist_ratio(hthetaphi_tg, hthetaphi_bg, time_ratio, h_name, h_title)

            h_ratio_thetaphi.Write()

        fout.Save()
        fout.Close()

        if log_fname is not None:
            print('\n', time.time(), file=open(log_fname, 'a'))
            print('GetMeasRatio for Target:%s and Background:%s, output as %s'
                  % (tg_stats_fname, bg_stats_fname, ratio_fname),
                  file=open(log_fname, 'a'))

    @classmethod
    def _ratio_to_leff(cls, h, phys_model, coord: str, eps):
        xbins = h.GetNbinsX()
        ybins = h.GetNbinsY()

        h_leff = h.Clone()
        h_leff.Reset("ICES")
        h_name = h.GetName()
        h_name = h_name.replace('ratio', 'Leff')
        h_title = h.GetName()
        h_title = h_title.replace('ratio', 'Leff') + ';%s;%s;Leff [m g/cm^{3}]' % (h.GetXaxis().GetTitle(),
                                                                    h.GetYaxis().GetTitle())
        h_leff.SetNameTitle(h_name, h_title)

        for i_xbin in range(1, xbins + 1):
            print('\r%s: %d / %d = %.2f    ' % (h.GetName(), i_xbin, xbins, 1. * i_xbin / xbins), end='')
            xx = h.GetXaxis().GetBinCenter(i_xbin)
            for i_ybin in range(1, ybins + 1):
                yy = h.GetYaxis().GetBinCenter(i_ybin)

                if coord.lower().__contains__('phi'):
                    theta = yy
                else:
                    dic_thetaphi = CommonUtils.cvt_thetaxy_to_thetaphi(thetax=yy, thetay=xx)
                    theta = dic_thetaphi.get('theta')

                i_bin = h.GetBin(i_xbin, i_ybin)

                ratio = h.GetBinContent(i_bin)
                err_ratio = h.GetBinError(i_bin)

                if ratio >= 1:
                    leff = 0
                elif ratio < 1e-10:
                    leff = -10
                else:
                    leff = phys_model.get_Leff_from_ratio(theta, ratio)

                if leff > -9:
                    dl_dr = phys_model.get_demin_dratio(theta, ratio, eps, leff)
                    err_leff = dl_dr * err_ratio

                    h_leff.SetBinContent(i_bin, leff)
                    h_leff.SetBinError(i_bin, err_leff)
        print('')
        return h_leff

    def GetMeasLeff(self, ratio_fname, leff_fname, phys_model, eps=0.0005, log_fname=None):
        fin = TFile(ratio_fname, 'read')

        hratio_thetaphi_name = _hthetaphi_name_prefix + '_ratio'
        hratio_thetaphi_name = ROOTUtils.find_key_in_file(ratio_fname, hratio_thetaphi_name)
        hratio_thetaxy_name = _hthetaxy_name_prefix + '_ratio'
        hratio_thetaxy_name = ROOTUtils.find_key_in_file(ratio_fname, hratio_thetaxy_name)

        fout = TFile(leff_fname, 'update')

        hratio_thetaphi = fin.Get(hratio_thetaphi_name)
        hleff_thetaphi = self._ratio_to_leff(hratio_thetaphi, phys_model, 'thetaphi', eps)
        hleff_thetaphi.Write()

        hratio_thetaxy = fin.Get(hratio_thetaxy_name)
        hleff_thetaxy = self._ratio_to_leff(hratio_thetaxy, phys_model, 'thetaxy', eps)
        hleff_thetaxy.Write()

        fout.Save()
        fout.Close()

        if log_fname is not None:
            print('\n', time.time(), file=open(log_fname, 'a'))
            print('GetMeasLeff for %s, output as %s.\nPhysical model:' % (ratio_fname, leff_fname),
                  file=open(log_fname, 'a'))
            phys_model.print_physModel(log_fname)
            print('Epsilon = %f' % eps, file=open(log_fname, 'a'))
