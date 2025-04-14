import os.path

from DataFactory import CommonUtils
from ROOT import TMath, TH2F, TVector3, TCanvas, TString, TFile
import logging

logger = logging.getLogger(__name__)


class Detector:

    def __init__(self, det_name: str, dic_det_info: dict, bins_dic: dict):
        self._det_name = det_name
        self._dic_det_info = dic_det_info
        self._bins_dic = bins_dic

        self._ax = self._ay = self._h = None # unit is mm
        self._x = self._y = self._z = None  # unit is m
        self._rot_thetax = self._rot_thetay = self._rot_phi = None
        self._meas_time = None

        self._h_accept_thetaphi = None
        self._h_accept_thetaxy = None

    """
    get acceptance hist if it exists, or calculate it out if not
    """

    def get_acceptance_plane_det(self, prefix=None):
        if self._rot_thetax is None:
            lis = self.get_rotAngle()
        if self._ax is None:
            lis = self.get_size_plane_det()

        if self._h_accept_thetaphi is None:
            self._calc_accept_thetaphi()
        if self._h_accept_thetaxy is None:
            self._calc_accept_thetaxy()

        if prefix is None:
            return self._h_accept_thetaphi, self._h_accept_thetaxy

        fout = TFile('%s_%s.root' % (prefix, self._det_name), 'recreate')
        self._h_accept_thetaphi.Write()
        self._h_accept_thetaxy.Write()
        fout.Save()
        fout.Close()

        return self._h_accept_thetaphi, self._h_accept_thetaxy

    def _fun_acc_cart(self, theta=None, phi=None, thetax=None, thetay=None):
        v_cart = TVector3()

        if theta is not None and phi is not None and thetax is None and thetay is None:
            # figure out thetax_trk and thetay_trk in the tracker coordinate
            v_cart.SetMagThetaPhi(1, theta, phi)
        elif theta is None and phi is None and thetax is not None and thetay is not None:
            # figure out thetax_trk and thetay_trk in the tracker coordinate
            v_cart.SetXYZ(TMath.Tan(thetax), TMath.Tan(thetay), 1)
        else:
            logger.fatal('Detector._fun_acc_cart: Pass in only thetaphi or thetaxy')
            raise ValueError('Detector._fun_acc_cart: Pass in only thetaphi or thetaxy')

        v_tracker = v_cart
        if self._rot_thetax < 1e-10:
            v_tracker.RotateX(-1. * self._rot_thetay)
        if self._rot_thetay < 1e-10:
            v_tracker.RotateY(-1. * self._rot_thetax)
        thetax_trk = TMath.Pi() / 2
        thetay_trk = TMath.Pi() / 2
        if v_tracker(2) > 1e-10:
            thetax_trk = TMath.ATan(v_tracker(0) / v_tracker(2))
            thetay_trk = TMath.ATan(v_tracker(1) / v_tracker(2))

        tmp_acc1 = (1 - self._h / self._ax * abs(TMath.Tan(thetax_trk)))
        tmp_acc2 = (1 - self._h / self._ay * abs(TMath.Tan(thetay_trk)))
        acc = 0.

        if tmp_acc1 > 1e-10 and tmp_acc2 > 1e-10:
            acc = tmp_acc1 * tmp_acc2

        return acc

    @classmethod
    def _calc_projection_given_theta(cls, theta):
        return TMath.Cos(theta)

    def _calc_accept_thetaphi(self):
        # using _fun_acc_thetaphi_cart to fill the acceptance hist thetaphi
        theta_bins = self._bins_dic['theta_bins']
        phi_bins = self._bins_dic['phi_bins']

        self._h_accept_thetaphi = TH2F("h_accept_thetaphi_%s" % self._det_name, "h_accept_thetaphi_%s" % self._det_name,
                                       phi_bins[0], phi_bins[1], phi_bins[2],
                                       theta_bins[0], theta_bins[1], theta_bins[2])

        xbins = self._h_accept_thetaphi.GetNbinsX()
        ybins = self._h_accept_thetaphi.GetNbinsY()

        for i in range(1, xbins + 1):
            phi = self._h_accept_thetaphi.GetXaxis().GetBinCenter(i)
            for j in range(1, ybins + 1):
                theta = self._h_accept_thetaphi.GetYaxis().GetBinCenter(j)
                ibin = self._h_accept_thetaphi.GetBin(i, j)

                acc = self._fun_acc_cart(theta=theta, phi=phi)
                prj = self._calc_projection_given_theta(theta)

                self._h_accept_thetaphi.SetBinContent(ibin, acc * prj)

        # set error as 0, meaning that no uncertainty
        bins_thetaphi = self._h_accept_thetaphi.GetNbinsX() * \
                        self._h_accept_thetaphi.GetNbinsY() * self._h_accept_thetaphi.GetNbinsZ()
        for i in range(1, bins_thetaphi + 1):
            self._h_accept_thetaphi.SetBinError(i, 0)

    def _calc_accept_thetaxy(self):
        # using _fun_acc_thetaxy_cart to fill the acceptance hist thetaxy
        thetax_bins = self._bins_dic['thetax_bins']
        thetay_bins = self._bins_dic['thetay_bins']

        self._h_accept_thetaxy = TH2F("h_accept_thetaxy_%s" % self._det_name, "h_accept_thetaxy_%s" % self._det_name,
                                      thetay_bins[0], thetay_bins[1], thetay_bins[2],
                                      thetax_bins[0], thetax_bins[1], thetax_bins[2])

        xbins = self._h_accept_thetaxy.GetNbinsX()
        ybins = self._h_accept_thetaxy.GetNbinsY()

        for i in range(1, xbins + 1):
            thetay = self._h_accept_thetaxy.GetXaxis().GetBinCenter(i)
            for j in range(1, ybins + 1):
                thetax = self._h_accept_thetaxy.GetYaxis().GetBinCenter(j)
                ibin = self._h_accept_thetaxy.GetBin(i, j)

                vec = TVector3(TMath.ATan(thetax), TMath.ATan(thetay), 1.)
                theta = vec.Theta()

                acc = self._fun_acc_cart(thetax=thetax, thetay=thetay)
                prj = self._calc_projection_given_theta(theta)

                self._h_accept_thetaxy.SetBinContent(ibin, acc * prj)

        # set error as 0, meaning that no uncertainty
        bins_thetaxy = self._h_accept_thetaxy.GetNbinsX() * \
                       self._h_accept_thetaxy.GetNbinsY() * self._h_accept_thetaxy.GetNbinsZ()
        for i in range(1, bins_thetaxy + 1):
            self._h_accept_thetaxy.SetBinError(i, 0)

    """
    get variables
    """

    def get_det_name(self):
        return self._det_name

    def get_det_id(self):
        if self._dic_det_info.get('DetID', '') == '':
            raise ValueError('Please set detector DetID')
        return int(self._dic_det_info.get('DetID'))

    def get_positions(self):
        # positions of the detector
        if self._dic_det_info.get('Position', '') == '':
            raise ValueError('Please set detector Position')
        self._x = float(self._dic_det_info.get('Position')[0])
        self._y = float(self._dic_det_info.get('Position')[1])
        self._z = float(self._dic_det_info.get('Position')[2])
        lis = [self._x, self._y, self._z]
        return lis

    def get_size_plane_det(self):
        # some properties of the detector
        if self._dic_det_info.get('Size', '') == '':
            raise ValueError('Please set detector Size')
        self._ax = float(self._dic_det_info.get('Size')[0])
        self._ay = float(self._dic_det_info.get('Size')[1])
        self._h = float(self._dic_det_info.get('Size')[2])
        return [self._ax, self._ay, self._h]

    def get_rotAngle(self):
        # the rotation angle of the detector
        if self._dic_det_info.get('RotAngle', '') == '':
            raise ValueError('Please set detector RotAngle')
        self._rot_thetax = float(self._dic_det_info.get('RotAngle')[0]) * TMath.Pi() / 180
        self._rot_thetay = float(self._dic_det_info.get('RotAngle')[1]) * TMath.Pi() / 180
        self._rot_phi = float(self._dic_det_info.get('RotAngle')[2]) * TMath.Pi() / 180
        lis = [self._rot_thetax, self._rot_thetay, self._rot_phi]
        return lis

    def get_det_type(self):
        det_type = self._dic_det_info.get('Type', '')
        det_type = det_type.lower()
        if det_type == 'planea' or det_type.__contains__('xujia'):
            det_type_str = 'plane_a'
            det_type_id = 1
        elif det_type == 'boreholea' or det_type.__contains__('zhuodai'):
            det_type_str = 'borehole_a'
            det_type_id = 2
        elif det_type == 'boreholeb' or det_type.__contains__('kaiqiang'):
            det_type_str = 'borehole_b'
            det_type_id = 3
        elif det_type == 'boreholec' or det_type.__contains__('baopeng'):
            det_type_str = 'borehole_c'
            det_type_id = 4
        elif det_type == 'planeb' or det_type.__contains__('wenjing') or det_type.__contains__('jiangkun'):
            det_type_str = 'plane_b'
            det_type_id = 5
        elif det_type == 'plane':
            det_type_str = 'plane_a'
            det_type_id = 1
            print("ATTENTION: type 'plane' is defaulted as 'planeA' --------------.")
        else:
            raise ValueError('Unknown detector type: %s\n'
                             'The known types are: -----------------\n'
                             '* planeA: plane detector assembled by triangular scintillators, '
                             'size: 50cm * 50cm, developed by Xujia Luo in 2019;\n'
                             '* planeB: plane detector with larger and complete scintillator plane, '
                             'planed to be developed by Wenjing Liu and Jiangkun Li\n'
                             '* boreholeA: cylinder detector used in bolehole, '
                             'diameter: 26cm, developed by Zhuodai Li;\n'
                             '* boreholeB: borehole detector used solid scintillators, '
                             'planed to be developed by Kaiqiang Yao;\n'
                             '* boreholeC: borehole detector used liquid scintillators, '
                             'planed to be developed by Baopeng Su.\n'
                             'Thanks a lot for their contributions.\n\n'
                             'Please input a type name or the developer\'s given name above.' % det_type)

        return det_type_id, det_type_str

    def get_meas_time(self):
        # the measurement duration of the detector
        if self._meas_time is None:
            meas_time = self._dic_det_info.get('MeasTime', '')
            if meas_time == '':
                self._meas_time = 0
            else:
                self._meas_time = int(CommonUtils.time_convert(meas_time))
        return self._meas_time

    def get_background_meas_time(self):
        # the background measurement duration of the detector
        bg_meas_time = self._dic_det_info.get('BackgroundMeasTime', '')
        return int(CommonUtils.time_convert(bg_meas_time))

    def get_det_pro_fname(self):
        det_pro_fname = self._dic_det_info.get('DetProperty', '')
        if det_pro_fname == '':
            raise Exception('Cannot fine DetProperty for %s' % self.get_det_name())
        det_pro_fname = os.path.join(self._dic_det_info.get('GlobalDir', ''), det_pro_fname)
        return det_pro_fname

    def get_Z_aver(self):
        zaver = self._dic_det_info.get('ZAver', '')
        if isinstance(zaver, list) and len(zaver) == 2:
            z_aver = [eval(i) for i in zaver]
            return z_aver
        print('Do not get ZAver.')
        return None
