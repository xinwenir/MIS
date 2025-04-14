import os.path

from DataFactory.PhysModel import PhysModelFromFormula

from ROOT import TMath, TString, TH2F, TFile
import ROOT
from array import array
import logging

logger = logging.getLogger(__name__)


class FluxModelFromFormula(PhysModelFromFormula):
    @classmethod
    def _fill_dic_xy_range(cls, dic_xy_range: dict = None):
        if dic_xy_range is None:
            dic_xy_range = {}

        dic_range = {'e_low': eval(dic_xy_range.get('e_low', '0.106')),
                     'e_up': eval(dic_xy_range.get('e_up', '50')),
                     'theta_low': eval(dic_xy_range.get('theta_low', '0')),
                     'theta_up': eval(dic_xy_range.get('theta_up', 'TMath.Pi()/2')),
                     'e_bins': int(eval(dic_xy_range.get('e_bins', '100'))),
                     'theta_bins': int(eval(dic_xy_range.get('theta_bins', '60')))}

        return dic_range

    def _get_TF2(self, dic_xy_range: dict = None):
        dic_xy_range = self._fill_dic_xy_range(dic_xy_range)

        e_low = dic_xy_range.get('e_low')
        e_up = dic_xy_range.get('e_up')
        theta_low = dic_xy_range.get('theta_low')
        theta_up = dic_xy_range.get('theta_up')

        if self._flux_id == 1:
            f2 = ROOT.func2d_bogdanova(e_low, e_up, theta_low, theta_up)
        elif self._flux_id == 2:
            f2 = ROOT.func2d_reyna_bugaev(e_low, e_up, theta_low, theta_up)
        elif self._flux_id == 3:
            f2 = ROOT.func2d_reyna_hebbeker(e_low, e_up, theta_low, theta_up)
        else:
            raise ValueError('FluxID must be 1, 2 or 3. %d is given' % self._flux_id)

        f2.SetName("f_flux_model_%d" % self._flux_id)
        return f2

    def DrawFluxModel(self, out_fname: str, dic_xy_range: dict = None):
        f2 = self._get_TF2(dic_xy_range)

        fout = TFile(out_fname, 'update')
        f2.Write()

        fout.Save()
        fout.Close()

    @classmethod
    def _integral2D(cls, h2, h_name, h_title, x_or_y: str):
        if x_or_y.lower().__contains__('x'):
            h1 = h2.ProjectionX().Clone()
            xbin = h2.GetNbinsX()
        else:
            h1 = h2.ProjectionY().Clone()
            xbin = h2.GetNbinsY()
        h1.SetNameTitle(h_name, h_title)
        h1.Reset("ICES")

        for ix in range(1, xbin):
            err = array('d', [0.])
            if x_or_y.lower().__contains__('x'):
                val_int = h2.IntegralAndError(ix, ix, 1, -1, err)
            else:
                val_int = h2.IntegralAndError(1, -1, ix, ix, err)

            h1.SetBinContent(ix, val_int)
            h1.SetBinError(ix, err[0])

        return h1

    @classmethod
    def _hist_divide_sintheta(cls, h2, x_or_y: str):
        h2_multi_sintheta = h2.Clone()
        h2_multi_sintheta.SetName(h2.GetName() + '_multi_sintheta')

        if x_or_y.lower().__contains__('x'):
            thetabin = h2.GetNbinsX()
            ebin = h2.GetNbinsY()

            for ix in range(1, thetabin + 1):
                theta = h2.GetXaxis().GetBinCenter(ix)
                for iy in range(1, ebin + 1):
                    ibin = h2.GetBin(ix, iy)
                    h2_multi_sintheta.SetBinContent(ibin, h2.GetBinContent(ibin) * TMath.Sin(theta))
        else:
            thetabin = h2.GetNbinsY()
            ebin = h2.GetNbinsX()

            for ix in range(1, ebin + 1):
                for iy in range(1, thetabin + 1):
                    theta = h2.GetYaxis().GetBinCenter(iy)
                    ibin = h2.GetBin(ix, iy)
                    h2_multi_sintheta.SetBinContent(ibin, h2.GetBinContent(ibin) * TMath.Sin(theta))

        return h2_multi_sintheta

    def SampleFromFluxModel(self, out_fname: str, n_entries: int, dic_xy_range: dict = None, scale: float = 1.):
        f2 = self._get_TF2(dic_xy_range)

        dic_xy_range = self._fill_dic_xy_range(dic_xy_range)

        e_low = dic_xy_range.get('e_low')
        e_up = dic_xy_range.get('e_up')
        theta_low = dic_xy_range.get('theta_low')
        theta_up = dic_xy_range.get('theta_up')
        e_bins = dic_xy_range.get('e_bins')
        theta_bins = dic_xy_range.get('theta_bins')

        h2 = TH2F("h_flux_model_%d" % self._flux_id, ";Energy [GeV];theta [rad]", e_bins, e_low, e_up, theta_bins, theta_low, theta_up)
        h2.FillRandom(f2.GetName(), n_entries)

        h2.Scale(scale)

        h2_multi_sintheta = self._hist_divide_sintheta(h2, 'y')

        h_eng = self._integral2D(h2_multi_sintheta, "%s_energy" % h2_multi_sintheta.GetName(), '; Energy [GeV]', 'x')
        h_theta = self._integral2D(h2_multi_sintheta, "%s_energy" % h2_multi_sintheta.GetName(), '; theta [rad]', 'y')

        fout = TFile(out_fname, 'update')
        h2.Write()
        h2_multi_sintheta.Write()
        h_eng.Write()
        h_theta.Write()

        fout.Save()
        fout.Close()
