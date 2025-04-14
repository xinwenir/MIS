from ROOT import TH2F, TMath, TFile, TRandom, TVector3, gROOT
from DataFactory import CommonUtils, ROOTUtils
from array import array
import os
import logging

logger = logging.getLogger(__name__)


class AngleMesh:
    """
    create histograms for two angular mesh, and calculate the solid angle in each bin
    """

    def __init__(self):
        self._h_thetaphi_solid_angle = None
        self._h_thetaxy_solid_angle = None

    def write_solid_angle(self, bins_dic, solid_angle_fname='', nrand=int(1e7), enforce=False):
        if solid_angle_fname != '' and os.path.isfile(solid_angle_fname) and not enforce:
            b_get_hists = self._read_solid_angle(solid_angle_fname, bins_dic)
            if not b_get_hists:
                self._calc_solid_angle(bins_dic=bins_dic, nrand=nrand)
            else:
                logger.info('Get hists of solid angle from file: %s', solid_angle_fname)
                return 1
        else:
            self._calc_solid_angle(bins_dic=bins_dic, nrand=nrand)

        fout = TFile(solid_angle_fname, 'recreate')
        fout.cd()
        self._h_thetaphi_solid_angle.Write()
        self._h_thetaxy_solid_angle.Write()
        fout.Save()
        fout.Close()

    def get_thetaphi_solid_angle(self):
        return self._h_thetaphi_solid_angle

    def get_thetaxy_solid_angle(self):
        return self._h_thetaxy_solid_angle

    def _fill_sphere_random(self, nrand):
        my_rand = TRandom()

        _x = array('d', [0])
        _y = array('d', [0])
        _z = array('d', [0])

        for i in range(nrand):
            my_rand.Sphere(_x, _y, _z, 1.)
            _z[0] = abs(_z[0])
            v_sphere = TVector3(_x[0], _y[0], _z[0])

            theta = v_sphere.Theta()
            phi = v_sphere.Phi()
            if v_sphere(2) < 1e-10:
                thetax = TMath.Pi() / 2
                thetay = TMath.Pi() / 2
            else:
                thetax = TMath.ATan(v_sphere(0) / v_sphere(2))
                thetay = TMath.ATan(v_sphere(1) / v_sphere(2))

            self._h_thetaphi_solid_angle.Fill(phi, theta)
            self._h_thetaxy_solid_angle.Fill(thetay, thetax)

        scale = 2. * TMath.Pi() / nrand
        self._h_thetaphi_solid_angle.Scale(scale)
        self._h_thetaxy_solid_angle.Scale(scale)

        # set error as 0, meaning that no uncertainty
        bins_thetaphi = self._h_thetaphi_solid_angle.GetNbinsX() * \
                        self._h_thetaphi_solid_angle.GetNbinsY() * self._h_thetaphi_solid_angle.GetNbinsZ()
        for i in range(1, bins_thetaphi + 1):
            self._h_thetaphi_solid_angle.SetBinError(i, 0)
        bins_thetaxy = self._h_thetaxy_solid_angle.GetNbinsX() * \
                       self._h_thetaxy_solid_angle.GetNbinsY() * self._h_thetaxy_solid_angle.GetNbinsZ()
        for i in range(1, bins_thetaxy + 1):
            self._h_thetaxy_solid_angle.SetBinError(i, 0)

    """calculate solid angle in each bin by random filling"""

    def _calc_solid_angle(self, bins_dic, nrand):

        theta_bins = bins_dic['theta_bins']
        phi_bins = bins_dic['phi_bins']
        thetax_bins = bins_dic['thetax_bins']
        thetay_bins = bins_dic['thetay_bins']

        self._h_thetaphi_solid_angle = TH2F("h_thetaphi_solidAngle", "h_thetaphi_solidAngle",
                                            phi_bins[0], phi_bins[1], phi_bins[2], theta_bins[0], theta_bins[1],
                                            theta_bins[2])
        self._h_thetaxy_solid_angle = TH2F("h_thetaxy_solidAngle", "h_thetaxy_solidAngle",
                                           thetay_bins[0], thetay_bins[1], thetay_bins[2],
                                           thetax_bins[0], thetax_bins[1], thetax_bins[2])
        self._fill_sphere_random(nrand)

    """read hists of solid angle from file"""

    def _read_solid_angle(self, solid_angle_fname, bins_dic):
        _fin = TFile(solid_angle_fname, 'read')
        if not (_fin.GetListOfKeys().Contains('h_thetaphi_solidAngle') and
                _fin.GetListOfKeys().Contains('h_thetaxy_solidAngle')):
            return False
        self._h_thetaphi_solid_angle = _fin.Get("h_thetaphi_solidAngle")
        self._h_thetaphi_solid_angle.SetDirectory(gROOT)
        self._h_thetaxy_solid_angle = _fin.Get("h_thetaxy_solidAngle")
        self._h_thetaxy_solid_angle.SetDirectory(gROOT)
        _fin.Close()

        # judge whether hists in the file are valid
        b_same_size = True
        theta_bins = bins_dic['theta_bins']
        phi_bins = bins_dic['phi_bins']
        thetax_bins = bins_dic['thetax_bins']
        thetay_bins = bins_dic['thetay_bins']
        tmp_h_thetaphi = TH2F("tmp_h_thetaphi", "tmp_h_thetaphi", phi_bins[0], phi_bins[1], phi_bins[2],
                              theta_bins[0], theta_bins[1], theta_bins[2])
        tmp_h_thetaxy = TH2F("tmp_h_thetaxy", "tmp_h_thetaxy", thetay_bins[0], thetay_bins[1], thetay_bins[2],
                             thetax_bins[0], thetax_bins[1], thetax_bins[2])

        if not ROOTUtils.check_two_hists_the_same_size(self._h_thetaphi_solid_angle, tmp_h_thetaphi):
            b_same_size = False
        if not ROOTUtils.check_two_hists_the_same_size(self._h_thetaxy_solid_angle, tmp_h_thetaxy):
            b_same_size = False

        if not b_same_size:
            self._h_thetaphi_solid_angle = None
            self._h_thetaxy_solid_angle = None
            return False

        return True
