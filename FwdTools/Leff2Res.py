from ROOT import TH2F, TFile, gROOT, TVector3, TMath
from DataFactory import PhysModel, CommonUtils
from FwdTools import dic_hist_name
import logging

logger = logging.getLogger(__name__)


class Leff2Res:
    def __init__(self):
        self._hthetaphi_ratio = None
        self._hthetaxy_ratio = None
        self._hthetaphi_flux = None
        self._hthetaxy_flux = None
        self._hthetaphi_flux_x_SldAgl = None
        self._hthetaxy_flux_x_SldAgl = None
        self._hthetaphi_counting_rate = None
        self._hthetaxy_counting_rate = None
        self._hthetaphi_counts = None
        self._hthetaxy_counts = None
        self._h_thetaphi_solid_angle = None
        self._h_thetaxy_solid_angle = None
        self._h_accept_thetaphi = None
        self._h_accept_thetaxy = None

    def write_pred_res(self, fout_name):
        lis_need_check = [self._hthetaphi_ratio, self._hthetaxy_ratio,
                          self._hthetaphi_flux, self._hthetaxy_flux,
                          self._hthetaphi_flux_x_SldAgl, self._hthetaxy_flux_x_SldAgl,
                          self._hthetaphi_counting_rate, self._hthetaxy_counting_rate,
                          self._hthetaphi_counts, self._hthetaxy_counts,
                          self._h_thetaphi_solid_angle, self._h_thetaxy_solid_angle,
                          self._h_accept_thetaphi, self._h_accept_thetaxy]
        fout = TFile(fout_name, 'recreate')
        for i in lis_need_check:
            if i is not None:
                i.Write()

        fout.Save()
        fout.Close()

    def calc_ratio(self, phys_model: PhysModel, leff_fname, opt='', eps=0.001):
        # read leff hist from file
        _f = TFile(leff_fname, 'read')
        _hthetaphi_leff = _f.Get(dic_hist_name['hthetaphi_leff'])
        _hthetaxy_leff = _f.Get(dic_hist_name['hthetaxy_leff'])
        _hthetaphi_leff.SetDirectory(gROOT)
        _hthetaxy_leff.SetDirectory(gROOT)

        if not opt.__contains__('thetax'):
            self._calc_ratio_given_leff(phys_model, _hthetaphi_leff, 'thetaphi', eps)

        if not opt.__contains__('phi'):
            self._calc_ratio_given_leff(phys_model, _hthetaxy_leff, 'thetaxy', eps)


    def _calc_ratio_given_leff(self, phys_model: PhysModel, h_leff, opt: str, eps):
        if opt.lower().__contains__('phi'):
            mode = 0
        elif opt.lower().__contains__('thetax'):
            mode = 1
        else:
            raise ValueError('opt should be thetaphi or thetaxy; opt = %s' % opt)

        h_ratio = h_leff.Clone()
        h_ratio.SetDirectory(gROOT)
        h_ratio.Reset("ICES")

        xbins = h_leff.GetNbinsX()
        ybins = h_leff.GetNbinsY()

        for i_xbin in range(1, xbins + 1):
            xx = h_leff.GetXaxis().GetBinCenter(i_xbin)
            for i_ybin in range(1, ybins + 1):
                yy = h_leff.GetYaxis().GetBinCenter(i_ybin)

                ibin = h_leff.GetBin(i_xbin, i_ybin)
                leff = h_leff.GetBinContent(ibin)
                leff_err = h_leff.GetBinError(ibin)

                if leff < 1e-15:
                    continue

                if mode == 0:
                    theta = yy
                else:
                    vec = TVector3(TMath.Tan(yy), TMath.Tan(xx), 1)
                    theta = vec.Theta()

                _ratio = phys_model.get_ratio_from_Leff(theta, leff)
                h_ratio.SetBinContent(ibin, _ratio)

                if leff_err > 1e-10:
                    dr_dl = phys_model.get_dflux_dleff(theta, leff, eps, _ratio)
                    ratio_err = dr_dl * leff_err
                else:
                    ratio_err = 0
                h_ratio.SetBinError(ibin, ratio_err)

        if mode == 0:
            h_ratio.SetNameTitle(dic_hist_name['hthetaphi_ratio'],
                                 '%s; %s; ratio' % (dic_hist_name['hthetaphi_ratio'], dic_hist_name['thetaphi_axis_title']))
            self._hthetaphi_ratio = h_ratio
        else:
            h_ratio.SetNameTitle(dic_hist_name['hthetaxy_ratio'],
                                 '%s; %s; ratio' % (dic_hist_name['hthetaxy_ratio'], dic_hist_name['thetaxy_axis_title']))
            self._hthetaxy_ratio = h_ratio

    def calc_flux_free_sky(self, phys_model: PhysModel, _h_accept_thetaphi, _h_accept_thetaxy, opt='', eps=0.001):
        if not opt.__contains__('thetax'):
            self._calc_free_sky_flux_given_bin(phys_model, _h_accept_thetaphi, 'thetaphi', eps)

        if not opt.__contains__('phi'):
            self._calc_free_sky_flux_given_bin(phys_model, _h_accept_thetaxy, 'thetaxy', eps)

    def _calc_free_sky_flux_given_bin(self, phys_model: PhysModel, h_leff, opt: str, eps=0.001, meas_time: int = None):
        mode = CommonUtils.check_angle_mesh_opt(opt)

        h_flux = h_leff.Clone()
        h_flux.Reset("ICES")
        h_flux.SetDirectory(gROOT)

        xbins = h_leff.GetNbinsX()
        ybins = h_leff.GetNbinsY()

        for i_xbin in range(1, xbins + 1):
            xx = h_leff.GetXaxis().GetBinCenter(i_xbin)
            for i_ybin in range(1, ybins + 1):
                yy = h_leff.GetYaxis().GetBinCenter(i_ybin)

                ibin = h_leff.GetBin(i_xbin, i_ybin)

                if mode == 0:
                    theta = yy
                else:
                    vec = TVector3(TMath.Tan(yy), TMath.Tan(xx), 1)
                    theta = vec.Theta()

                _flux = phys_model.get_flux_from_Leff(theta, 0)

                h_flux.SetBinContent(ibin, _flux)
                h_flux.SetBinError(ibin, 0)

        if mode == 0:
            h_flux.SetNameTitle(dic_hist_name['hthetaphi_flux'],
                                '%s; %s; flux' % (dic_hist_name['hthetaphi_flux'], dic_hist_name['thetaphi_axis_title']))
            self._hthetaphi_flux = h_flux
        else:
            h_flux.SetNameTitle(dic_hist_name['hthetaxy_flux'],
                                '%s; %s' % (dic_hist_name['hthetaxy_flux'], dic_hist_name['thetaxy_axis_title']))
            self._hthetaxy_flux = h_flux

    def calc_flux(self, phys_model: PhysModel, leff_fname, opt='', eps=0.001):
        # read leff hist from file
        _f = TFile(leff_fname, 'read')
        _hthetaphi_leff = _f.Get(dic_hist_name['hthetaphi_leff'])
        _hthetaxy_leff = _f.Get(dic_hist_name['hthetaxy_leff'])
        _hthetaphi_leff.SetDirectory(gROOT)
        _hthetaxy_leff.SetDirectory(gROOT)

        if not opt.__contains__('thetax'):
            self._calc_flux_given_leff(phys_model, _hthetaphi_leff, 'thetaphi', eps)

        if not opt.__contains__('phi'):
            self._calc_flux_given_leff(phys_model, _hthetaxy_leff, 'thetaxy', eps)

    def _calc_flux_given_leff(self, phys_model: PhysModel, h_leff, opt: str, eps=0.001):
        mode = CommonUtils.check_angle_mesh_opt(opt)

        h_flux = h_leff.Clone()
        h_flux.Reset("ICES")
        h_flux.SetDirectory(gROOT)

        xbins = h_leff.GetNbinsX()
        ybins = h_leff.GetNbinsY()

        for i_xbin in range(1, xbins + 1):
            xx = h_leff.GetXaxis().GetBinCenter(i_xbin)
            for i_ybin in range(1, ybins + 1):
                yy = h_leff.GetYaxis().GetBinCenter(i_ybin)

                ibin = h_leff.GetBin(i_xbin, i_ybin)
                leff = h_leff.GetBinContent(ibin)
                leff_err = h_leff.GetBinError(ibin)

                if mode == 0:
                    theta = yy
                else:
                    vec = TVector3(TMath.Tan(yy), TMath.Tan(xx), 1)
                    theta = vec.Theta()

                _flux = phys_model.get_flux_from_Leff(theta, leff)
                h_flux.SetBinContent(ibin, _flux)

                if leff_err > 1e-10:
                    df_dl = phys_model.get_dflux_dleff(theta, leff, eps, _flux)
                    flux_err = df_dl * leff_err
                else:
                    flux_err = 0.
                h_flux.SetBinError(ibin, flux_err)

        if mode == 0:
            h_flux.SetNameTitle(dic_hist_name['hthetaphi_flux'],
                                '%s; %s; flux' % (dic_hist_name['hthetaphi_flux'], dic_hist_name['thetaphi_axis_title']))
            self._hthetaphi_flux = h_flux
        else:
            h_flux.SetNameTitle(dic_hist_name['hthetaxy_flux'],
                                '%s; %s' % (dic_hist_name['hthetaxy_flux'], dic_hist_name['thetaxy_axis_title']))
            self._hthetaxy_flux = h_flux

    def calc_flux_x_SldAgl(self, h_solid_angle, opt: str):
        mode = CommonUtils.check_angle_mesh_opt(opt)
        if mode == 0:
            if self._hthetaphi_flux is None:
                raise Exception('self._hthetaphi_flux is None; Run [calc_flux] first.')
            if not h_solid_angle.GetName().__contains__('phi'):
                logger.warning('h_solid_angle\'s name dose\'nt have "phi". It may be wrong.')
            h_flux = self._hthetaphi_flux.Clone()
            self._h_thetaphi_solid_angle = h_solid_angle
        else:
            if self._hthetaxy_flux is None:
                raise Exception('self._hthetaxy_flux is None; Run [calc_flux] first.')
            if not h_solid_angle.GetName().__contains__('thetax'):
                logger.warning('h_solid_angle\'s name dose\'nt have "thetax". It may be wrong.')
            h_flux = self._hthetaxy_flux.Clone()
            self._h_thetaxy_solid_angle = h_solid_angle

        h_flux_x_SldAgl = h_flux.Clone()
        h_flux_x_SldAgl.SetDirectory(gROOT)
        h_flux_x_SldAgl.Multiply(h_solid_angle)

        max_bin = h_flux_x_SldAgl.GetMaximumBin()
        if h_flux.GetBinError(max_bin) > 1e-5:
            nbin = h_flux_x_SldAgl.GetNbinsX() * h_flux_x_SldAgl.GetNbinsY()

            for ibin in range(1, nbin + 1):
                err_flux = h_flux.GetBinError(ibin)
                scale = h_solid_angle.GetBinContent(ibin)
                h_flux_x_SldAgl.SetBinError(ibin, err_flux * scale)

        if mode == 0:
            h_flux_x_SldAgl.SetNameTitle(dic_hist_name['hthetaphi_flux_x_SldAgl'],
                                         '%s; %s' % (dic_hist_name['hthetaphi_flux_x_SldAgl'],
                                                     dic_hist_name['thetaphi_axis_title']))
            self._hthetaphi_flux_x_SldAgl = h_flux_x_SldAgl
        else:
            h_flux_x_SldAgl.SetNameTitle(dic_hist_name['hthetaxy_flux_x_SldAgl'],
                                         '%s; %s' % (dic_hist_name['hthetaxy_flux_x_SldAgl'],
                                                     dic_hist_name['thetaxy_axis_title']))
            self._hthetaxy_flux_x_SldAgl = h_flux_x_SldAgl

    def calc_counting_rate(self, h_accpt, area, opt: str):
        mode = CommonUtils.check_angle_mesh_opt(opt)
        if mode == 0:
            if self._hthetaphi_flux_x_SldAgl is None:
                raise Exception('self._hthetaphi_flux_x_SldAgl is None; Run [calc_flux_x_SldAgl] first.')
            if not h_accpt.GetName().__contains__('phi'):
                logger.warning('h_accpt\'s name does not have "phi". It may be wrong.')
            h_flux_x_SldAgl = self._hthetaphi_flux_x_SldAgl.Clone()
            self._h_accept_thetaphi = h_accpt
        else:
            if self._hthetaxy_flux_x_SldAgl is None:
                raise Exception('self._hthetaxy_flux_x_SldAgl is None; Run [calc_flux_x_SldAgl] first.')
            if not h_accpt.GetName().__contains__('thetax'):
                logger.warning('h_accpt\'s name does not have "thetax". It may be wrong.')
            h_flux_x_SldAgl = self._hthetaxy_flux_x_SldAgl.Clone()
            self._h_accept_thetaxy = h_accpt

        h_counting = h_flux_x_SldAgl.Clone()
        h_counting.SetDirectory(gROOT)
        h_counting.Multiply(h_accpt)
        h_counting.Scale(area)

        # calculate uncertainty
        max_bin = h_counting.GetMaximumBin()
        if h_flux_x_SldAgl.GetBinError(max_bin) > 1e-5:
            nbin = h_counting.GetNbinsX() * h_counting.GetNbinsY()

            for ibin in range(1, nbin + 1):
                err_flux_x_SldAgl = h_flux_x_SldAgl.GetBinError(ibin)
                scale = h_accpt.GetBinContent(ibin)
                h_counting.SetBinError(ibin, err_flux_x_SldAgl * scale)

        if mode == 0:
            h_counting.SetNameTitle(dic_hist_name['hthetaphi_counting_rate'],
                                    '%s; %s' % (dic_hist_name['hthetaphi_counting_rate'],
                                                dic_hist_name['thetaphi_axis_title']))
            self._hthetaphi_counting_rate = h_counting
        else:
            h_counting.SetNameTitle(dic_hist_name['hthetaxy_counting_rate'],
                                    '%s; %s' % (dic_hist_name['hthetaxy_counting_rate'],
                                                dic_hist_name['thetaxy_axis_title']))
            self._hthetaxy_counting_rate = h_counting

    def calc_counts(self, meas_time, opt: str):
        mode = CommonUtils.check_angle_mesh_opt(opt)
        if mode == 0:
            if self._hthetaphi_counting_rate is None:
                raise Exception('self._hthetaphi_flux_x_SldAgl is None; Run [calc_flux_x_SldAgl] first.')
            h_counting_rate = self._hthetaphi_counting_rate.Clone()
        else:
            if self._hthetaxy_counting_rate is None:
                raise Exception('self._hthetaxy_flux_x_SldAgl is None; Run [calc_flux_x_SldAgl] first.')
            h_counting_rate = self._hthetaxy_counting_rate.Clone()

        h_counts = h_counting_rate.Clone()
        h_counts.SetDirectory(gROOT)
        h_counts.Scale(meas_time)

        max_bin = h_counts.GetMaximumBin()
        if h_counting_rate.GetBinError(max_bin) > 1e-5:
            nbin = h_counts.GetNbinsX() * h_counts.GetNbinsY()
            for ibin in range(1, nbin + 1):
                err_counting_rate = h_counting_rate.GetBinError(ibin)
                h_counts.SetBinError(ibin, err_counting_rate * meas_time)
        else:
            nbin = h_counts.GetNbinsX() * h_counts.GetNbinsY()
            for ibin in range(1, nbin + 1):
                val = h_counts.GetBinError(ibin)
                h_counts.SetBinError(ibin, TMath.Sqrt(val))

        if mode == 0:
            h_counts.SetNameTitle(dic_hist_name['hthetaphi_counts'],
                                  '%s; %s' % (dic_hist_name['hthetaphi_counts'],
                                              dic_hist_name['thetaphi_axis_title']))
            self._hthetaphi_counts = h_counts
        else:
            h_counts.SetNameTitle(dic_hist_name['hthetaxy_counts'],
                                  '%s; %s' % (dic_hist_name['hthetaxy_counts'],
                                              dic_hist_name['thetaxy_axis_title']))
            self._hthetaxy_counts = h_counts
