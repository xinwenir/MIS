import array
import math
import os

from ROOT import TH1F, TH2F, TFile, TVector3, TMath, TCanvas, TLine, TLegend, TList, THStack
from ROOT import gROOT, gStyle, gPad
from DataFactory import CommonUtils, ROOTUtils
import pandas as pd


class PrintTObject:
    def __init__(self):
        pass

    @classmethod
    def PrintLeffHist2Obs(cls, fin_name, det, fout, opt='thetaphi', leff_uncertainty_level=0.05, if_print_zero=True):
        # find the hist name and get it
        if opt.lower().__contains__('phi'):
            mode = 0
            h_name = ROOTUtils.find_key_in_file(fin_name, 'htheta_phi')
        elif opt.lower().__contains__('thetax'):
            mode = 1
            h_name = ROOTUtils.find_key_in_file(fin_name, 'hthetax_thetay')
        else:
            h_name = ''
            raise ValueError('opt should be thetaphi pr thetaxy')

        if h_name == '':
            raise Exception('Cannot find hist name in file %s.' % fin_name)

        fin = TFile(fin_name, 'read')

        h = fin.Get(h_name)

        rot_phi = det.get_rotAngle()[2]

        xbins = h.GetXaxis().GetNbins()
        ybins = h.GetYaxis().GetNbins()

        for i_x in range(1, xbins + 1):
            xx = h.GetXaxis().GetBinCenter(i_x)

            for i_y in range(1, ybins + 1):
                yy = h.GetYaxis().GetBinCenter(i_y)

                if mode == 0:
                    theta = yy
                    phi = xx + rot_phi
                else:
                    vec = TVector3(TMath.Tan(yy), TMath.Tan(xx), 1)
                    theta = vec.Theta()
                    phi = vec.Phi() + rot_phi

                ibin = h.GetBin(i_x, i_y)
                leff = h.GetBinContent(ibin)

                leff_err = leff * leff_uncertainty_level
                if abs(leff_err) < 1e-5:
                    leff_err = 1e-5

                str_print = '%i %.3f %.3f %.3f %.3f %.3f %.3f %.5f %.5f' % (
                    det.get_det_id(), det.get_positions()[0], det.get_positions()[1], det.get_positions()[2],
                    theta, phi, -1, leff, abs(leff_err))
                if not if_print_zero:
                    if abs(leff) > 1e-5:
                        print(str_print, file=fout)
                else:
                    print(str_print, file=fout)

    @classmethod
    def _check_draw_option_dict(cls, dic_draw_option: dict):
        dic_draw_option_new = {'TH1F': dic_draw_option.get('TH1F', 'histe'),
                               'TH1D': dic_draw_option.get('TH1D', 'histe'),
                               'TH2F': dic_draw_option.get('TH2F', 'colz'),
                               'TH2D': dic_draw_option.get('TH2D', 'colz'),
                               'TCanvas': dic_draw_option.get('TCanvas', ''),
                               'THStack': dic_draw_option.get('TCanvas', 'nostack'),
                               'TMultiGraph': dic_draw_option.get('TMultiGraph', 'apl'),
                               'TGraphErrors': dic_draw_option.get('TGraphErrors', 'ap'),
                               'gStyle.SetOptStat': eval(dic_draw_option.get('gStyle.SetOptStat', '0'))}

        return dic_draw_option_new

    @classmethod
    def _check_draw_option_str(cls, draw_option: str):
        dic_option = {1: 'E',
                      2: 'conterr',
                      3: 'texterr',
                      4: 'legoerr'}
        try:
            opt = int(eval(draw_option))
        except:
            if isinstance(draw_option, str) and draw_option.lower() == 'e':
                opt = 1
            elif isinstance(draw_option, str) and draw_option.lower().__contains__('conterr'):
                opt = 2
            elif isinstance(draw_option, str) and draw_option.lower().__contains__('texterr'):
                opt = 3
            elif isinstance(draw_option, str) and draw_option.lower().__contains__('legoerr'):
                opt = 4
            else:
                opt = 1
        print('Draw with option ', opt, dic_option[opt])

        return opt

    def PrintTFile2PDF(self, root_name, pdf_fname, dic_draw_option=None, if_only_canvas=False, if_show=False):
        if not if_show:
            gROOT.SetBatch()
        if dic_draw_option is None:
            dic_draw_option = {}
        dic_draw_option = self._check_draw_option_dict(dic_draw_option)

        fin = TFile(root_name, 'read')

        can = TCanvas('can_PrintTFile2PDF', 'can_PrintTFile2PDF', 2000, 1500)

        if 'gStyle.SetOptStat' in dic_draw_option.keys():
            gStyle.SetOptStat(dic_draw_option.get('gStyle.SetOptStat'))

        for i in fin.GetListOfKeys():
            val = fin.Get(i.GetName())
            class_name = val.Class_Name()

            print('Draw ', i)

            if if_only_canvas:
                if class_name == 'TCanvas':
                    val.Print(pdf_fname + '(')
                continue

            if class_name == 'TCanvas':
                val.Print(pdf_fname + '(')
            elif class_name == 'TLegend':
                continue
            elif class_name == 'THStack' or class_name == 'TMultiGraph':
                can.cd()
                val.Draw(dic_draw_option.get(class_name, ''))
                lg_val = fin.Get('lg_' + i.GetName())
                lg_val.SetFillStyle(0)
                lg_val.Draw()
                can.Print(pdf_fname + '(')
            else:
                can.cd()
                val.Draw(dic_draw_option.get(class_name, ''))
                can.Print(pdf_fname + '(')

        can.Print(pdf_fname + ']')

        fin.Close()

    def PrintTH2WithError(self, lis_root_name, lis_hname, pdf_fname, draw_option=None, if_setoptstat=0):
        opt = self._check_draw_option_str(draw_option)

        gROOT.SetBatch()

        if not if_setoptstat == 1:
            if_setoptstat = 0
            print('set: gStyle.SetOptStat(0)')
        gStyle.SetOptStat(if_setoptstat)

        fout_name = pdf_fname.replace('.pdf', '_canvas.root')
        fout = TFile(fout_name, 'recreate')
        can = TCanvas('can_PrintTH2WithError', 'can_PrintTH2WithError', 2000, 1500)

        for i in range(len(lis_root_name)):
            fin_name = lis_root_name[i]
            hname = lis_hname[i]

            hname = ROOTUtils.find_key_in_file(fin_name, hname)
            fin = TFile(fin_name, 'read')
            h = fin.Get(hname)
            h.SetName(h.GetName()+os.path.splitext(os.path.split(fin_name)[1])[0])

            h.SetDirectory(gROOT)
            herr = h.Clone()
            herr.SetName('err_'+h.GetName())

            nbins = h.GetNbinsX() * h.GetNbinsY() * h.GetNbinsZ()
            for ib in range(1, nbins + 1):
                err = h.GetBinError(ib)
                herr.SetBinContent(ib, err)

            fout.cd()
            h.Write()
            herr.Write()

            if opt == 1:
                h.Draw('E')
            elif opt == 2:
                h.Draw('colz')
                herr.Draw('cont1 same')
            elif opt == 3:
                h.Draw('colz')
                herr.Draw('text same')
            elif opt == 4:
                h1 = h.Clone()
                h1.SetName('low_'+h.GetName())
                h1.Add(herr, -1)
                h2 = herr.Clone()
                h2.SetName('2err_' + h.GetName())
                h2.Scale(2)

                h1.SetFillColor(2)
                h2.SetFillColor(1)

                hs = THStack('hs_' + h.GetName(), 'hs_' + h.GetName())
                hs.Add(h1)
                hs.Add(h2)

                hs.Draw('lego2')
                # h1.Draw('surf1 same')

                h1.Write()
                h2.Write()
                hs.Write()

            can.Write()
            can.Print(pdf_fname + '(')
            fin.Close()

        can.Print(pdf_fname + ']')

        fout.Save()
        fout.Close()

    @staticmethod
    def DrawBaselineThreshold(root_name, pdf_fname, baseline_fname, threshold_fname, if_show=False):
        if not if_show:
            gROOT.SetBatch()

        df_baseline = pd.read_csv(baseline_fname)
        df_threshold = pd.read_csv(threshold_fname)

        # read hists from root file
        fin = TFile(root_name, 'read')

        can = TCanvas('can_DrawBaselineThreshold', 'can_DrawBaselineThreshold', 2000, 1500)
        fout = TFile(pdf_fname.replace(".pdf","_baseline.root"), 'recreate')

        for i in fin.GetListOfKeys():
            if not i.GetName().__contains__('h_energy_spectrum_board'):
                continue

            h_name = i.GetName()
            h1 = fin.Get(h_name)

            t_name = h_name[h_name.index('board') + 5:]
            [i_board, i_channel] = [int(j) for j in t_name.split('_chn')]

            val_baseline = df_baseline.iloc[i_board, i_channel]
            val_threshold = df_threshold.iloc[i_board, i_channel]

            if val_threshold < 1:
                continue

            h1_smooth = h1.Clone()
            h1_smooth.Smooth(2)
            h1_smooth.GetXaxis().SetRangeUser(val_threshold + 10, val_threshold + 300)
            peak = h1_smooth.GetMaximum()
            height = peak * 2 + 10

            h1.GetXaxis().SetRangeUser(val_baseline - 50, val_threshold + 320)
            h1.SetMaximum(height + 1)

            line_baseline = TLine(val_baseline, 0, val_baseline, height)
            line_threshold = TLine(val_threshold, 0, val_threshold, height)

            lg = TLegend(0.78, 0.67, 0.98, 0.77)
            lg.AddEntry(line_baseline, 'baseline', 'l')
            lg.AddEntry(line_threshold, 'threshold', 'l')

            can.cd()
            h1.Draw('hist')
            line_baseline.SetLineColor(2)
            line_threshold.SetLineColor(1)
            line_baseline.Draw()
            line_threshold.Draw()
            lg.Draw()

            can.Print('%s(' % pdf_fname)

            fout.cd()
            can.Write()

        can.Print('%s]' % pdf_fname)

    @staticmethod
    def PrintAllCountsMeasTime(root_dir, fout_name):
        df_root_file = pd.DataFrame(columns=['total counts',
                                             'meas_time / second', 'count rate / per second',
                                             'meas_time / day', 'count rate / per day'])
        for root, dirs, files in os.walk(root_dir):
            df_root_file.loc[root] = ['~~~~~~', '~~~~~~', '~~~~~~', '~~~~~~', '~~~~~~']
            for name in files:
                if os.path.splitext(name)[1] == ".root":
                    try:
                        fin = TFile(os.path.join(root, name), 'read')
                        tree = fin.Get('t_tree')
                        n_rows = tree.GetEntries()
                        meas_time_sec = fin.Get('meas_time')[0]
                        meas_time_day = meas_time_sec / 60. / 60. / 24.
                        df_root_file.loc[os.path.join(root, name)] = [n_rows,
                                                                      meas_time_sec, 1. * n_rows / meas_time_sec,
                                                                      meas_time_day, 1. * n_rows / meas_time_day]
                    except:
                        try:
                            fin = TFile(os.path.join(root, name), 'read')
                            meas_time_sec = fin.Get('meas_time')[0]
                            meas_time_day = meas_time_sec / 60. / 60. / 24.
                            df_root_file.loc[os.path.join(root, name)] = ['-',
                                                                          meas_time_sec, '-',
                                                                          meas_time_day, '-']
                        except:
                            df_root_file.loc[os.path.join(root, name)] = ['-', '-', '-', '-', '-']

        df_root_file.to_csv(fout_name)

    @staticmethod
    def PrintTree2Csv(fin_name, tree_name, lis_var_name: list, fout_name):
        if tree_name == '':
            tree_name = 'tr_'

        tree_name = ROOTUtils.find_key_in_file(fin_name, tree_name)

        fin = TFile(fin_name, 'read')
        tree = fin.Get(tree_name)

        df_tr = ROOTUtils.cvt_tree_to_dataframe(tree)

        try:
            df_var = df_tr[lis_var_name]
        except Exception as e:
            print(e)
            for i_var in lis_var_name:
                if i_var not in df_tr.columns:
                    print('%s is not in given tree' % i_var)
            raise Exception('Insure the input variables are in the given tree %s.' % tree.GetName())

        df_var.to_csv(fout_name, sep=' ', index=False)
