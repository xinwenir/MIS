import re
import sys
import time
from array import array
import os, logging

import numpy as np
import pandas as pd
from scipy.stats import linregress

default_value = -100000

from DatabaseTool import *
from DataFactory import CommonUtils, ROOTUtils
from FwdTools import dic_hist_name

from ROOT import TH1F, TFile, TMath, TTree, RDF, TVectorD, TVector3, TGraph, TMultiGraph, TLegend, TGraphErrors

logger = logging.getLogger(__name__)


class TimeSeriesTools(object):
    def __init__(self):
        pass

    @classmethod
    def _get_total_counts_in_time_series(cls, lis_h_sel, hname, htitle):
        nbins_h_total = len(lis_h_sel)
        h_total = TH1F(hname, htitle, nbins_h_total, 0, nbins_h_total)
        for i in range(len(lis_h_sel)):
            ih = lis_h_sel[i]
            lis_time_reg = re.search(r"\d+_\d+", ih.GetName()).group().split("_")
            start_time = int(lis_time_reg[0])
            end_time = int(lis_time_reg[1])

            dt_start_time = CommonUtils.convert_timeTag_with_timeZone(start_time)
            dt_end_time = CommonUtils.convert_timeTag_with_timeZone(end_time)

            # get total counts and error
            err = array('d', [0.])
            nbinsx_ih = ih.GetNbinsX()
            nbinsy_ih = ih.GetNbinsY()
            total_counts = ih.IntegralAndError(1, nbinsx_ih, 1, nbinsy_ih, err)
            h_total.SetBinContent(i+1, total_counts)
            h_total.SetBinError(i+1, err[0])

            h_total.GetXaxis().SetBinLabel(i + 1, '%s - %s'
                                           % (dt_start_time.strftime('%m.%d %H:%M'),
                                              dt_end_time.strftime('%m.%d %H:%M')))

        return h_total

    @classmethod
    def _get_selected_hist_list(cls, lis_hcounts, hleff_sel, suffix):
        lis_h_sel = []
        for hts in lis_hcounts:
            # check the dimension of two hists
            if not ROOTUtils.check_two_hists_the_same_size(hts, hleff_sel):
                raise Exception("The dimension or size of the hist to be selected and "
                                "the Leff hist are not the same."
                                "hist to be selected: (%d,%d); Leff hist: (%d,%d)" %
                                (hts.GetXbins(), hts.GetYbins(),
                                 hleff_sel.GetXbins(), hleff_sel.GetYbins()))

            h_sel = hts.Clone()
            h_sel.SetName(hts.GetName() + suffix)
            h_sel.Multiply(hleff_sel)

            lis_h_sel.append(h_sel)

        return lis_h_sel

    def GetCountsTimeSeries(self, hts_fname, hts_reg, leff_file, sel_exp, fout_name):
        fleff = TFile(leff_file, "read")

        slice_xy_title = None
        if hts_reg.__contains__("phi"):
            hleff_name = dic_hist_name['hthetaphi_leff']
            slice_xy_title = ["#phi [rad]", "#theta [rad]"]
        elif hts_reg.__contains__("x"):
            hleff_name = dic_hist_name['hthetaxy_leff']
            slice_xy_title = ["#theta_{y} [rad]", "#theta_{x} [rad]"]
        else:
            raise Exception("The HistRegularExpression should be in (theta, phi) or in (thetax, thetay)")

        hleff = fleff.Get(hleff_name)
        hleff_sel = ROOTUtils.convert_hist_to_bool_with_selection(hleff, sel_exp)

        fout = TFile(fout_name, "recreate")

        # get list of selected hist
        fin = TFile(hts_fname, 'read')
        lis_hcounts = []
        for i in fin.GetListOfKeys():
            if not re.search(hts_reg, i.GetName()):
                continue

            hts = fin.Get(i.GetName())
            lis_hcounts.append(hts)

        # get total counts in time series
        lis_h_sel = self._get_selected_hist_list(lis_hcounts, hleff_sel, "_sel")
        h_total = self._get_total_counts_in_time_series(lis_h_sel,
                                                        "htotal_counts_selected",
                                                        "Total Counts In Selection;;counts")
        fout.cd()
        h_total.Write()

        # get sliced counts in time series
        xbins = hleff_sel.GetNbinsX()
        ybins = hleff_sel.GetNbinsY()
        # slice x
        lis_slice_x = []
        lis_slice_x_label = []
        for i in range(1, xbins + 1):
            htemp = hleff_sel.Clone()
            htemp.Reset("ICES")
            htemp.SetName(hleff_sel.GetName() + "_sliceX%d" % i)

            for j in range(1, ybins + 1):
                ibin = hleff_sel.GetBin(i, j)
                htemp.SetBinContent(ibin, hleff_sel.GetBinContent(ibin))

            if htemp.Integral()<1e-10:
                continue

            hname = "hcounts_selected_sliceX%d" % i
            x = htemp.GetXaxis().GetBinCenter(i)
            htitle = ("Counts In Selection %s = %.2f %s;;counts" %
                      (slice_xy_title[0].split(" ")[0], x, slice_xy_title[0].split(" ")[1]))

            lis_h_sel_slice_x = self._get_selected_hist_list(lis_h_sel, htemp, "_sliceX%d" % i)
            h_slice = self._get_total_counts_in_time_series(lis_h_sel_slice_x, hname, htitle)

            if h_slice.Integral()>0:
                fout.cd()
                h_slice.Write()
                lis_slice_x.append(h_slice)
                lis_slice_x_label.append("%s = %.2f" % (slice_xy_title[0].split(" ")[0], x))

        hname = "hcounts_selected_sliceX"
        h2d_slicex, h2d_slicex_norm = ROOTUtils._get_h2d_with_lis_h1d(lis_slice_x, lis_slice_x_label, hname)
        h2d_slicex.Write()
        h2d_slicex_norm.Write()

        # slice y
        lis_slice_y = []
        lis_slice_y_label = []
        for j in range(1, ybins + 1):
            htemp = hleff_sel.Clone()
            htemp.Reset("ICES")
            htemp.SetName(hleff_sel.GetName() + "_sliceY%d" % j)

            for i in range(1, xbins + 1):
                ibin = hleff_sel.GetBin(i, j)
                htemp.SetBinContent(ibin, hleff_sel.GetBinContent(ibin))

            if htemp.Integral() < 1e-10:
                continue

            hname = "hcounts_selected_sliceY%d" % j
            y = htemp.GetYaxis().GetBinCenter(j)
            htitle = ("Counts In Selection %s = %.2f %s;;counts" %
                      (slice_xy_title[1].split(" ")[0], y, slice_xy_title[1].split(" ")[1]))

            lis_h_sel_slice_y = self._get_selected_hist_list(lis_h_sel, htemp, "_sliceY%d" % j)
            h_slice = self._get_total_counts_in_time_series(lis_h_sel_slice_y, hname, htitle)

            if h_slice.Integral()>0:
                fout.cd()
                h_slice.Write()
                lis_slice_y.append(h_slice)
                lis_slice_y_label.append("%s = %.2f" % (slice_xy_title[1].split(" ")[0], y))

        hname = "hcounts_selected_sliceY"
        h2d_slicey, h2d_slicey_norm = ROOTUtils._get_h2d_with_lis_h1d(lis_slice_y, lis_slice_y_label, hname)
        h2d_slicey.Write()
        h2d_slicey_norm.Write()

        fout.Save()
        fout.Close()


    def ScaleHCountsByTotalAirCounts(self, hts_fname, hts_reg, air_fname, air_hname, fout_name):
        fair = TFile(air_fname, "read")

        hair = fair.Get(air_hname)
        n_time_slice = hair.GetNbinsX()

        fout = TFile(fout_name, "recreate")

        # get list of selected hist
        fin = TFile(hts_fname, 'read')
        lis_hcounts = []
        for i in fin.GetListOfKeys():
            if not re.search(hts_reg, i.GetName()):
                continue

            hts = fin.Get(i.GetName())
            lis_hcounts.append(hts)

        if len(lis_hcounts) != n_time_slice:
            raise Exception("The scaling hist has %d time slices, but there are %d hists to be scaled."
                            % (n_time_slice, len(lis_hcounts)))

        fout = TFile(fout_name, "recreate")

        h_all_time = lis_hcounts[0].Clone()
        h_all_time.Reset("ICES")
        h_all_time.SetName("hcounts_all_time")

        for i in range(len(lis_hcounts)):
            ih = lis_hcounts[i]
            ih_scaled = ih.Clone()
            ih_scaled.SetName(ih_scaled.GetName() + "_scaled")

            scale_val = hair.GetBinContent(i+1)
            ih_scaled.Scale(scale_val)

            ih_scaled.Write()
            h_all_time.Add(ih_scaled)

        h_all_time.Write()

        fout.Save()
        fout.Close()
