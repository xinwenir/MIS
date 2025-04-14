import os.path

from DataFactory import CommonUtils, ROOTUtils
from ROOT import TFile, gROOT, TF1, TF2

from ROOT import TMath, TString
import ROOT
from array import array
import logging

logger = logging.getLogger(__name__)


class PhysicalModel:
    def __init__(self):
        self._relation_id = 0
        self._paras = None
        self._ecut = None

        (phys_calc_dir, fname) = os.path.split(__file__)
        phys_calc_fname = os.path.join(phys_calc_dir, 'PhysCalc.h')
        ROOT.gInterpreter.ProcessLine('#include ' + '\"' + phys_calc_fname + '\"')

    @classmethod
    def _read_relation(cls, rel):
        if rel == '':
            return 1
        mstr = TString(rel)

        if mstr.IsFloat():
            i = int(mstr.Atof())
            if i < 3:
                return i
            else:
                raise ValueError('Only have two relations:\nlinear: 1\t\tReyna-quadratic: 2')
        else:
            if rel.lower().__contains__('linear'):
                return 1
            elif rel.lower().__contains__('quadratic') or rel == 'Q':
                return 2
            else:
                raise ValueError('Cannot identify the given relation "%s"' % rel)

    def get_ecut(self):
        return self._ecut

    def get_relation_id(self):
        return self._relation_id

    def get_parameters(self):
        return self._paras

    def set_ecut(self, ecut=0.2):
        self._ecut = ecut

    def set_relation(self, relation=1):
        self._relation_id = self._read_relation(relation)

    def get_Emin_from_Leff(self, leff):
        if leff < 1e-10:
            return self._ecut

        x1 = array('d', [leff])
        p1 = array('d', self._paras)
        # x[0] = leff, p = parameters
        if self._relation_id == 1:
            emin = ROOT.leff_to_emin_linear(x1, p1)
        elif self._relation_id == 2:
            emin = ROOT.leff_to_emin_quadratic(x1, p1)
        else:
            raise Exception('Unknown relation_id %d\nInput 1: linear; 2: quadratic' % self._relation_id)

        return emin

    def get_Leff_from_Emin(self, emin, par: list):
        p = array('d', par)
        leff = None
        if len(par) == 2:
            leff = ROOT.emin_to_leff_linear(emin, p)
        elif len(par) == 3:
            leff = ROOT.emin_to_leff_quadratic(emin, p)
        else:
            raise Exception('get_Leff_from_Emin: input 2 or 3 parameters (%d are given).' % len(par))
        return leff

    @classmethod
    def _adjust_eps(cls, val, eps):
        while val / 4 < eps:
            eps = eps / 2
        return eps

    def get_dleff_dratio(self, theta, ratio, eps, leff=None):
        if leff is None:
            leff = self.get_Leff_from_ratio(theta, ratio)

        eps = self._adjust_eps(ratio, eps)
        ratio_minus = ratio - eps
        leff_minus = self.get_Leff_from_ratio(theta, ratio_minus)
        dl_dr = (leff - leff_minus) / eps
        return abs(dl_dr)

    def get_demin_dratio(self, theta, ratio, eps=0.001, emin=None):
        if ratio < 1e-15:
            return 0

        if emin is None:
            emin = self.get_Emin_from_ratio(theta, ratio)

        eps = self._adjust_eps(ratio, eps)

        ratio_minus = ratio - eps
        emin_minus = self.get_Emin_from_ratio(theta, ratio_minus)
        de_dr = (emin - emin_minus) / eps
        return abs(de_dr)

    def get_drario_dleff(self, theta, leff, eps, ratio=None):
        if ratio is None:
            ratio = self.get_ratio_from_Leff(theta, leff)

        eps = self._adjust_eps(leff, eps)

        leff_plus = leff + eps
        ratio_plus = self.get_ratio_from_Leff(theta, leff_plus)
        dr_dl = (ratio - ratio_plus) / eps
        return abs(dr_dl)

    def get_dflux_dleff(self, theta, leff, eps, flux=None):
        if flux is None:
            flux = self.get_flux_from_Leff(theta, leff)

        eps = self._adjust_eps(leff, eps)

        leff_plus = leff + eps
        flux_plus = self.get_flux_from_Leff(theta, leff_plus)
        df_dl = (flux - flux_plus) / eps
        return abs(df_dl)


class PhysModelFromFormula(PhysicalModel):
    def __init__(self, dic_phys_model):
        super().__init__()
        self.dic_phys_model = dic_phys_model
        self._ecut = float(dic_phys_model.get('Ecut', '0.2'))
        self._paras = [eval(i) for i in dic_phys_model.get('Parameters', ['0', '0.2'])]
        self._flux_id = self._read_flux_model(dic_phys_model.get('FluxModelFromFormula', ''))
        self._relation_id = self._read_relation(dic_phys_model.get('Relation', '1'))  # 1: linear; 2: quadratic

        self._dic_flux_id = {1: 'Bogdanova',
                             2: 'Reyna-Bugaev',
                             3: 'Reyna-Hebbeker'}
        self._dic_raltion_id = {1: 'linear',
                                2: 'quadratic'}

    def print_physModel(self, out_fname=None):
        if out_fname is not None:
            print('FluxModelFromFormula = %d(%s)' % (self._flux_id, self._dic_flux_id.get(self._flux_id, '')),
                  file=open(out_fname, 'a'))
            print('Relation = %d(%s)' % (self._relation_id, self._dic_raltion_id.get(self._relation_id, '')),
                  file=open(out_fname, 'a'))
            print('Parameters = %s' % ','.join([str(i) for i in self._paras]),
                  file=open(out_fname, 'a'))
            print('Ecut = %f' % self._ecut,
                  file=open(out_fname, 'a'))

        print('FluxModelFromFormula = %d(%s)' % (self._flux_id, self._dic_flux_id.get(self._flux_id, '')))
        print('Relation = %d(%s)' % (self._relation_id, self._dic_raltion_id.get(self._relation_id, '')))
        print('Parameters = %s' % ','.join([str(i) for i in self._paras]))
        print('Ecut = %f' % self._ecut)

    @classmethod
    def _read_flux_model(cls, model):
        if model == '':
            return 3
        mstr = TString(model)

        if mstr.IsFloat():
            i = int(mstr.Atof())
            if i < 4:
                return i
            else:
                raise ValueError('Only have three flux models:\nBogdanova: 1\t\tReyna-Hebbeker: 2\t\tReyna-Hebbeker: 3')
        else:
            if model.lower().__contains__('bogdanova'):
                return 1
            elif model.lower().__contains__('bugaev') or model == 'R-B':
                return 2
            elif model.lower().__contains__('hebbeker') or model == 'R-H':
                return 3
            else:
                raise ValueError('Cannot identify the given flux model "%s"' % model)

    def get_flux_id(self):
        return self._flux_id

    def get_flux_info(self):
        return 'FluxModelFromFormula = %i' % self.get_flux_id()

    def set_flux(self, flux=2):
        self._flux_id = self._read_flux_model(flux)

    def get_Leff_from_ratio(self, theta, ratio):
        if ratio > 1.:
            return 0.

        emin = ROOT.get_emin_from_ratio(theta, ratio, self._ecut, self._flux_id)
        p = array('d', self._paras)
        leff = ROOT.get_leff_from_emin(emin, self._ecut, self._relation_id, p)
        return leff

    def get_Emin_from_ratio(self, theta, ratio):
        if ratio > 1.:
            return self._ecut

        emin = ROOT.get_emin_from_ratio(theta, ratio, self._ecut, self._flux_id)

        return emin

    def get_ratio_from_Leff(self, theta, leff):
        emin = self.get_Emin_from_Leff(leff)

        x2 = array('d', [emin])
        p2 = array('d', [theta, self._ecut, self._flux_id])
        # x[0] = emin, p[0] = theta, p[1] = ecut, p[2] = flux_id
        ratio = ROOT.calc_ratio_with_emin(x2, p2)
        return ratio

    def get_ratio_from_Emin(self, theta, emin):
        if emin < self._ecut:
            return 1.

        x2 = array('d', [emin])
        p2 = array('d', [theta, self._ecut, self._flux_id])
        # x[0] = emin, p[0] = theta, p[1] = ecut, p[2] = flux_id
        ratio = ROOT.calc_ratio_with_emin(x2, p2)
        return ratio

    def get_flux_from_Leff(self, theta, leff):
        emin = self.get_Emin_from_Leff(leff)

        x2 = array('d', [emin])
        p2 = array('d', [theta, self._ecut, self._flux_id])
        # x[0] = emin, p[0] = theta, p[1] = ecut, p[2] = flux_id
        flux = ROOT.calc_flux_with_function(x2, p2)
        return flux

    def get_phys_model_1d_theta(self, theta, lmin=0, lmax=200):

        def get_ratio_from_length(arr, par):
            # arr[0] = length, par[0] = theta, par[1] = ecut, par[2] = flux_id
            # par[3] = _relation_id, par[3 +] = parameters of relation

            if par[3] == 1:
                emin = par[4] + arr[0] * par[5]
            elif par[3] == 2:
                emin = par[4] + arr[0] * par[5] + arr[0] * arr[0] * par[6]
            else:
                return -10000

            arr[0] = emin
            res = ROOT.calc_ratio_with_emin(arr, par)

            return res

        # // x[0] = length, p[0] = theta, p[1] = ecut, p[2] = flux_id
        #   // p[3] = _relation_id, p[3+] = parameters of relation
        lis_par = [theta, self._ecut, self._flux_id, self._relation_id]
        for i in self._paras:
            lis_par.append(i)
        par = array('d', lis_par)

        f1 = None
        if self._relation_id == 1:
            f1 = TF1("f1", get_ratio_from_length, lmin, lmax, 6)
            f1.SetParameters(par)
        elif self._relation_id == 2:
            f1 = TF1("f1", get_ratio_from_length, lmin, lmax, 7)
            f1.SetParameters(par)

        # check if it works
        if f1.Eval(10) < -10:
            raise Exception('Unknown relation_id %d\nInput 1: linear; 2: quadratic' % self._relation_id)

        return f1

    def get_phys_model_2d(self, theta, xmin=0, xmax=200, ):
        f2 = None
        if self._flux_id == 1:
            f2 = TF2("f1", ROOT.bogdanova_2d, xmin, xmax, 2)
        elif self._flux_id == 2:
            f2 = TF2("f1", ROOT.reyna_bugaev_2d, xmin, xmax, 2)
        elif self._flux_id == 3:
            f2 = TF2("f1", ROOT.reyna_hebbeker_2d, xmin, xmax, 2)

        f2.SetParameters(theta, self._ecut)
        return f2


class PhysModelFromHist(PhysicalModel):
    def __init__(self, dic_phys_model):
        super().__init__()
        self.dic_phys_model = dic_phys_model
        self._ecut = float(dic_phys_model.get('Ecut', '0.2'))
        self._paras = [eval(i) for i in dic_phys_model.get('Parameters', ['0', '0.2'])]
        self._relation_id = self._read_relation(dic_phys_model.get('Relation', '1'))  # 1: linear; 2: quadratic
        self._dic_raltion_id = {1: 'linear',
                                2: 'quadratic'}

        self._flux_fname = dic_phys_model.get('FluxModelRootName')
        self._flux_hist_name = dic_phys_model.get('FluxModelHistName')
        h = self._read_flux_hist(self._flux_fname, self._flux_hist_name)
        h.Smooth()
        self._lis_eup, self._lis_eup_bin = self._read_energy_up(h)
        self._elow = self._read_energy_low(h)
        self._flux_hist = self._hist_complement(h)
        self._ratio_hist = self._calc_ratio_hist(self._flux_hist)

    def print_physModel(self, out_fname=None):
        if out_fname is not None:
            print('FluxModelFromHist =\n%s: %s' % (self._flux_fname, self._flux_hist_name),
                  file=open(out_fname, 'a'))
            print('Relation = %d(%s)' % (self._relation_id, self._dic_raltion_id.get(self._relation_id, '')),
                  file=open(out_fname, 'a'))
            print('Parameters = %s' % ','.join([str(i) for i in self._paras]),
                  file=open(out_fname, 'a'))
            print('Ecut = %f' % self._ecut,
                  file=open(out_fname, 'a'))

        print('FluxModelFromHist =\n%s: %s' % (self._flux_fname, self._flux_hist_name))
        print('Relation = %d(%s)' % (self._relation_id, self._dic_raltion_id.get(self._relation_id, '')))
        print('Parameters = %s' % ','.join([str(i) for i in self._paras]))
        print('Ecut = %f' % self._ecut)

    @classmethod
    def _read_flux_hist(cls, fin_name, hist_name):
        hist_name = ROOTUtils.find_key_in_file(fin_name, hist_name)

        fin = TFile(fin_name, 'read')
        h = fin.Get(hist_name)
        h.SetDirectory(gROOT)

        fin.Close()

        return h

    @classmethod
    def _check_flux_model_xytitle(cls, h):
        if h.GetXaxis().GetTitle().lower().__contains__('theta'):
            raise Exception('The given flux model hist must be like: Energy-Theta(x-y).')
        ybin = h.GetNbinsY()
        if h.GetYaxis().GetBinCenter(ybin) > TMath.Pi() / 2:
            raise Exception('The given flux model hist must be like: Energy-Theta(x-y).')

    def _read_energy_up(self, h):
        self._check_flux_model_xytitle(h)

        xbin = h.GetNbinsX()
        ybin = h.GetNbinsY()

        lis_eup = []
        lis_eup_bin = []
        for iy in range(1, ybin + 1):
            h_eng = h.ProjectionX('h_pro_eng_%d' % iy, iy, iy)
            ibin = h_eng.GetMinimumBin()
            lis_eup_bin.append(ibin)
            lis_eup.append(h.GetXaxis().GetBinCenter(ibin) - h.GetXaxis().GetBinWidth(xbin) / 2)
        return lis_eup, lis_eup_bin

    def _read_energy_low(self, h):
        self._check_flux_model_xytitle(h)

        return h.GetXaxis().GetBinCenter(1) - h.GetXaxis().GetBinWidth(1) / 2

    def _hist_complement(self, h):
        eng_bin = h.GetNbinsX()
        ybin = h.GetNbinsY()

        for i in range(1, ybin + 1):
            theta = h.GetYaxis().GetBinCenter(i)
            x2 = array('d', [self._lis_eup[i - 1]])
            p2 = array('d', [theta, self._ecut, 1])
            # x[0] = emin, p[0] = theta, p[1] = ecut, p[2] = flux_id
            frac = ROOT.calc_ratio_with_emin(x2, p2)

            a = h.Integral(1, -1, i, i)
            # a / cont = (1-frac) / frac ==> cont = a * frac / (1 - frac)
            cont = a * frac / (1 - frac)
            ix = self._lis_eup_bin[i - 1]
            ibin = h.GetBin(ix, i)
            h.SetBinContent(ibin, cont)
            for j in range(ix + 1, eng_bin + 1):
                ibin = h.GetBin(j, i)
                h.SetBinContent(ibin, 0)

        return h

    def _calc_ratio_hist(self, h):
        self._check_flux_model_xytitle(h)

        eng_bin = h.GetNbinsX()
        theta_bin = h.GetNbinsY()

        h_ratio = h.Clone()
        h_ratio.Reset("ICES")

        for itheta in range(1, theta_bin + 1):
            total = h.Integral(1, -1, itheta, itheta)

            for ie in range(1, eng_bin + 1):
                val_int = h.Integral(ie, -1, itheta, itheta)

                ibin = h.GetBin(ie, itheta)
                h_ratio.SetBinContent(ibin, val_int / total)

        return h_ratio

    def get_flux_file_hist_name(self):
        return self._flux_fname, self._flux_hist_name

    def get_flux_info(self):
        return 'FluxModelRootName = %s\nFluxModelHistName = %s' \
            % (self.get_flux_file_hist_name()[0], self.get_flux_file_hist_name()[1])

    def get_flux_hist(self):
        return self._flux_hist

    def get_ratio_hist(self):
        return self._ratio_hist

    def set_flux(self, hist):
        self._flux_hist = hist

    @classmethod
    def _find_theta_bin(cls, h, theta):
        h_pro = h.ProjectionY("h_pro_theta")
        return h_pro.FindBin(theta)

    def _find_energy_bin_from_ratio(self, ratio, theta):
        """find in which bin the given emin is, support function for get_Emin_from_ratio"""
        theta_bin = self._find_theta_bin(self._ratio_hist, theta)
        h_pro_energy = self._ratio_hist.ProjectionX("h_pro_energy", theta_bin, theta_bin)

        return h_pro_energy.FindLastBinAbove(ratio)

    def get_Emin_from_ratio(self, theta, ratio):
        if ratio > 1.:
            return self._ecut

        ie = self._find_energy_bin_from_ratio(ratio, theta)
        if ie < 0:
            return -11
        itheta = self._find_theta_bin(self._ratio_hist, theta)
        if ie >= self._lis_eup_bin[itheta - 1]:
            return -10

        x1 = self._ratio_hist.GetXaxis().GetBinCenter(ie)
        x2 = self._ratio_hist.GetXaxis().GetBinCenter(ie)
        ibin1 = self._ratio_hist.GetBin(ie, itheta)
        y1 = self._ratio_hist.GetBinContent(ibin1)
        ibin2 = self._ratio_hist.GetBin(ie + 1, itheta)
        y2 = self._ratio_hist.GetBinContent(ibin2)

        emin = (ratio - y1) * (x1 - x2) / (y1 - y2) + x1

        return emin

    def get_Leff_from_ratio(self, theta, ratio):
        emin = self.get_Emin_from_ratio(theta, ratio)

        if emin < 0:
            return -10

        p = array('d', self._paras)
        leff = ROOT.get_leff_from_emin(emin, self._ecut, self._relation_id, p)
        return leff

    @classmethod
    def _integral_energy_for_certain_theta(cls, h, theta, emin):
        ibin = h.FindBin(emin, theta)

        a_xbin = array('i', [0])
        a_ybin = array('i', [0])
        a_zbin = array('i', [0])
        h.GetBinXYZ(ibin, a_xbin, a_ybin, a_zbin)

        val_int = h.Integral(a_xbin[0] + 1, -1, a_ybin[0], a_ybin[0])

        # calculate the fraction in the a_xbin[0]
        ebin_up = h.GetXaxis().GetBinCenter(a_xbin[0]) + h.GetXaxis().GetBinWidth(a_xbin[0]) / 2
        ebin_width = h.GetXaxis().GetBinWidth(a_xbin[0])
        val_ebin = h.GetBinContent(ibin)
        val_frac = (ebin_up - emin) / ebin_width * val_ebin

        return val_frac + val_int

    def get_ratio_from_Emin(self, theta, emin):
        itheta = self._find_theta_bin(self._ratio_hist, theta)
        if emin > self._lis_eup[itheta - 1]:
            logger.warning('Emin(Leff) is larger than the maximum energy flux hist gives.')
            return 0
        if emin < self._elow:
            return 1

        surv = self._integral_energy_for_certain_theta(self._flux_hist, theta, emin)
        total = self._integral_energy_for_certain_theta(self._flux_hist, theta, self._elow)

        return surv / total

    def get_ratio_from_Leff(self, theta, leff):
        emin = self.get_Emin_from_Leff(leff)

        return self.get_ratio_from_Emin(theta, emin)

    def get_flux_from_Leff(self, theta, leff):
        emin = self.get_Emin_from_Leff(leff)

        flux = self._integral_energy_for_certain_theta(self._flux_hist, theta, emin)

        return flux
