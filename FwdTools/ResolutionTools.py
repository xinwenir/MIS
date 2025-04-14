import os.path

from ROOT import TH2F, TFile, gROOT, TMath
import re

from FwdTools import dic_hist_name
from DataFactory import ROOTUtils
from FwdTools import dic_hist_name

import logging

logger = logging.getLogger(__name__)


class ResolutionTools:
    def __init__(self):
        pass

    @classmethod
    def _calc_meastime_with_hcounting_rate(cls, h_without, h_with, hname, htitle, sig):
        h_time = h_without.Clone()
        h_time.SetDirectory(gROOT)
        h_time.Reset("ICES")
        h_time.SetName(hname)
        h_time.SetTitle(htitle)
        h_time.GetZaxis().SetTitle("day")

        xbin = h_without.GetNbinsX()
        ybin = h_without.GetNbinsY()

        for i in range(1, xbin + 1):
            for j in range(1, ybin + 1):
                ibin = h_without.GetBin(i, j)
                v_without = h_without.GetBinContent(ibin)  # counting rate for without deposit
                v_with = h_with.GetBinContent(ibin)  # counting rate for with deposit

                if abs(v_without - v_with) > 1e-100:
                    t = sig * sig * v_without / (v_with - v_without) / (v_with - v_without)
                    t=t/24/60/60 # convert the time unit into day
                    h_time.SetBinContent(ibin, t)

        return h_time

    def _calc_meastime_by_given_file(self, f_without_name, f_with_name, suffix, opt="thetaphi", sig=2.):
        if opt.__contains__("phi"):
            opt = "thetaphi"
        elif opt.__contains__("x"):
            opt = "thetaxy"
        else:
            raise Exception("No such option %s\nAvailable options are thetaphi and thetaxthetay" % opt)

        suffix_with = suffix.split("-")[0]
        suffix_without = suffix.split("-")[1]

        f_without = TFile(f_without_name, "read")
        h_without = f_without.Get(dic_hist_name["h%s_counting_rate" % opt])
        h_without.SetDirectory(gROOT)
        h_without.SetName(dic_hist_name["h%s_counting_rate" % opt] + "_" + suffix_without)

        f_with = TFile(f_with_name, "read")
        h_with = f_with.Get(dic_hist_name["h%s_counting_rate" % opt])
        h_with.SetDirectory(gROOT)
        h_with.SetName(dic_hist_name["h%s_counting_rate" % opt] + "_" + suffix_with)

        if not ROOTUtils.check_two_hists_the_same_size(h_without, h_with):
            raise Exception("Two histograms don't have the same size: \n"
                            "%s in %s, %s in %s" %
                            (h_without.GetName(), f_without_name, h_with.GetName(), f_with_name))

        hname = "h_time_%s_" % opt + suffix
        htitle = ("h_time_%s_" % opt + suffix +
                  ";%s;%s" % (h_with.GetXaxis().GetTitle(), h_with.GetYaxis().GetTitle()))
        h_time = self._calc_meastime_with_hcounting_rate(h_without, h_with, hname, htitle, sig)

        return h_time

    def ExpectMeasTime(self, lis_f_without_name, lis_f_with_name, lis_suffix, fout_name, sig=2.):
        if not len(lis_f_without_name) == len(lis_f_with_name) or \
                not len(lis_f_without_name) == len(lis_suffix):
            raise Exception("The length of lis_f_without_name, lis_f_with_name and lis_suffix must match")

        fout = TFile(fout_name, "update")

        for i in range(len(lis_f_without_name)):
            f_without_name = lis_f_without_name[i]
            f_with_name = lis_f_with_name[i]
            suffix = lis_suffix[i]

            h_time = self._calc_meastime_by_given_file(f_without_name, f_with_name, suffix, "thetaphi", sig)
            fout.cd()
            h_time.Write()

            h_time = self._calc_meastime_by_given_file(f_without_name, f_with_name, suffix, "thetaxy", sig)
            fout.cd()
            h_time.Write()

        fout.Save()
        fout.Close()
