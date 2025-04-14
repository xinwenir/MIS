# from ROOT import TH1F, TH2F, TCanvas, TFile, THStack, TLegend, TVectorD, gROOT, TMath, TPad
# import ROOT
#
# from array import array
#
# from DataFactory import Utils
# from DataFactory.Utils import calc_significance
# import logging
#
# logger = logging.getLogger(__name__)
#
#
# def cvt_meas_time_to_day(meas_time):
#     return meas_time / 24 / 60 / 60
#
#
# class SliceTools:
#     def __init__(self):
#         self._fout = None
#         self._pdf_name = None
#         self._can = TCanvas("can_slice", "can_slice", 3000, 3000)
#
#         # ROOT.gInterpreter.ProcessLine('#include "./DataTools/PrintTFile.h"')
#
#     def compare_th2s(self, l_comp_name: list, l_fin_name: list, l_hist_name: list, fout_name: str, l_scale=None):
#         if not fout_name.__contains__('.root'):
#             fout_name = 'compare_' + '_'.join([i for i in l_comp_name]) + '.root'
#
#         logger.info('Compare 2d hists from files: %s.' % ', '.join(i for i in l_fin_name))
#
#         if not len(l_comp_name) == len(l_fin_name) or not len(l_comp_name) == len(l_hist_name):
#             raise Exception('[compare_th2s] Given lists do not have the same length.')
#
#         if l_scale is None:
#             l_scale = len(l_comp_name)*[1.]
#
#         for i_h in range(len(l_hist_name[0])):
#             l_h2 = []
#             l_meas_time = []
#             t_hist_names = [hs[i_h] for hs in l_hist_name]
#             for i in range(len(l_comp_name)):
#                 fin = TFile(l_fin_name[i], 'read')
#                 t_h2 = fin.Get(t_hist_names[i])
#                 t_h2.SetDirectory(gROOT)
#                 t_meas_time = fin.Get('meas_time')[0] / l_scale[i]
#                 l_h2.append(t_h2)
#                 l_meas_time.append(t_meas_time)
#
#             for i in range(1, len(l_h2)):
#                 if not Utils.check_two_hists_the_same_size(l_h2[0], l_h2[i]):
#                     logger.error('%s and %s do not have the same size' % (l_h2[0].GetName(), l_h2[i].GetName()))
#                     raise Exception('%s and %s do not have the same size' % (l_h2[0].GetName(), l_h2[i].GetName()))
#
#             xbin = l_h2[0].GetNbinsX()
#             ybin = l_h2[0].GetNbinsY()
#
#             self._fout = TFile(fout_name, 'recreate')
#             self._pdf_name = fout_name.replace('.root', '.pdf')
#
#             for i in range(1, xbin + 1):
#                 if l_h2[0].Integral(i, i, 1, -1) > 1:
#                     self._slice_th2s_1dim(l_h2, l_comp_name, l_meas_time, i, 'x')
#
#             for i in range(1, ybin + 1):
#                 if l_h2[0].Integral(1, -1, i, i) > 1:
#                     self._slice_th2s_1dim(l_h2, l_comp_name, l_meas_time, i, 'y')
#
#         self._can.Print('%s]' % self._pdf_name)
#         self._fout.Save()
#         self._fout.Close()
#
#     def _slice_th2s_1dim(self, l_h2: list, l_comp_name: list, l_meas_time: list, i_bin: int, x_or_y: str):
#         l_h1d = []
#         t_h2d = l_h2[0]
#
#         title_x = t_h2d.GetXaxis().GetTitle()
#         title_y = t_h2d.GetYaxis().GetTitle()
#
#         x_or_y = x_or_y.lower()
#         if (x_or_y.__contains__('x') and x_or_y.__contains__('y')) or \
#                 (not x_or_y.__contains__('x') and not x_or_y.__contains__('y')):
#             raise ValueError('[slice_th2s_1dim] Give function x or y')
#
#         slice_value = -10
#         hs_name = ''
#         hs_title = ''
#         if x_or_y.__contains__('x'):
#             slice_value = t_h2d.GetXaxis().GetBinCenter(i_bin)
#             title = 'x'
#             if not title_x == '':
#                 title = title_x.split()[0]
#             hs_name = 'CompareHists_%s_%d' % (title, i_bin)
#             hs_title = 'CompareHists_%s=%.2f; %s' % (title, slice_value, t_h2d.GetXaxis().GetTitle())
#
#             for i_h in range(len(l_h2)):
#                 t_h1 = l_h2[i_h].ProjectionY(hs_name.replace('CompareHists', '%s_%s' %
#                                                              (l_h2[i_h].GetName(), l_comp_name[i_h])), i_bin, i_bin)
#                 l_h1d.append(t_h1)
#
#         else:
#             slice_value = t_h2d.GetYaxis().GetBinCenter(i_bin)
#             title = 'y'
#             if not title_y == '':
#                 title = title_y.split()[0]
#             hs_name = 'CompareHists_%s_%d' % (title, i_bin)
#             hs_title = 'CompareHists_%s=%.2f' % (title, slice_value)
#
#             for i_h in range(len(l_h2)):
#                 t_h1 = l_h2[i_h].ProjectionX(hs_name.replace('CompareHists', '%s_%s' %
#                                                              (l_h2[i_h].GetName(), l_comp_name[i_h])),i_bin, i_bin)
#                 l_h1d.append(t_h1)
#
#         self._draw_h1stack(l_h1d, l_meas_time, hs_name, hs_title)
#
#     def _draw_h1stack(self, l_h1d, l_meas_time, hs_name, hs_title):
#         hs = THStack(hs_name, hs_title)
#         hs_sig = None
#         lg = TLegend(0.6, 0.8, 0.95, 0.93)
#         self._can.Clear()
#         self._can.cd()
#         self._fout.cd()
#
#         for i_h in range(len(l_h1d)):
#             l_h1d[i_h].Write()
#             meas_time = cvt_meas_time_to_day(l_meas_time[i_h])
#             l_h1d[i_h].Scale(1. / meas_time)
#             peak_area = l_h1d[i_h].Integral()
#             l_h1d[i_h].SetLineColor(i_h + 1)
#
#             hs.Add(l_h1d[i_h], "hist e")
#
#             hs_sig = None
#             if i_h > 0 and l_h1d[0].Integral() > 10 and l_h1d[i_h].Integral() > 10:
#                 chi2 = l_h1d[i_h].Chi2Test(l_h1d[0], "CHI2/NDF")
#                 hs_sig = THStack('%s_sig' % hs_name, '')
#                 h_sig = calc_significance(l_h1d[0], l_h1d[i_h],
#                                           '%s_%s_sig' % (l_h1d[0].GetName(), l_h1d[i_h].GetName()),
#                                           '%s_%s_sig' % (l_h1d[0].GetName(), l_h1d[i_h].GetName()), 'real')
#                 h_sig.SetLineColor(i_h + 1)
#                 hs_sig.Add(h_sig, "hist")
#                 lg.AddEntry(l_h1d[i_h], "%s, PeakArea=%.1f, Chi2Test=%.1f" % (l_h1d[i_h].GetName(), peak_area, chi2))
#             else:
#                 lg.AddEntry(l_h1d[i_h], "%s, PeakArea=%.1f" % (l_h1d[i_h].GetName(), peak_area))
#
#         self._can.cd()
#
#         pad1 = TPad("pad1", "pad1", 0, 0.33, 1, 1.0)
#         pad1.SetBottomMargin(0)
#         pad1.SetGridx()
#         pad1.Draw()
#         pad1.cd()
#         hs.Draw("nostack")
#         lg.Draw()
#
#         self._can.cd()
#         pad2 = TPad("pad2", "pad2", 0, 0.08, 1, 0.33)
#         pad2.SetTopMargin(0)
#         pad2.SetBottomMargin(0.2)
#         pad2.SetGridx()
#         pad2.Draw()
#         pad2.cd()
#
#         if hs_sig is not None:
#             hs_sig.Draw("nostack")
#
#             hs_sig.GetYaxis().SetLabelSize(15)
#             hs_sig.GetYaxis().SetTitle("significance ")
#             hs_sig.GetYaxis().SetNdivisions(505)
#             hs_sig.GetYaxis().SetTitleSize(30)
#             hs_sig.GetYaxis().SetTitleFont(43)
#             hs_sig.GetYaxis().SetTitleOffset(1.55)
#             hs_sig.GetYaxis().SetLabelFont(43)
#             hs_sig.GetYaxis().SetLabelSize(25)
#
#             hs_sig.GetXaxis().SetTitle(l_h1d[0].GetXaxis().GetTitle())
#             hs_sig.GetXaxis().SetTitleSize(30)
#             hs_sig.GetXaxis().SetTitleFont(43)
#             hs_sig.GetXaxis().SetTitleOffset(4.)
#             hs_sig.GetXaxis().SetLabelFont(43)
#             hs_sig.GetXaxis().SetLabelSize(25)
#
#             pad2.Update()
#
#         self._can.Print("%s(" % self._pdf_name)
#         hs.Write()
#
