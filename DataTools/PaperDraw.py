import array, os

import ROOT
from ROOT import TH1F, TH2F, TFile, TCanvas, TLegend, THStack, TMath, TList, TGraph, TMultiGraph
from ROOT import gROOT, gStyle, gPad
from DataFactory import CommonUtils, ROOTUtils


class PaperDraw:
    def __init__(self):
        pass
        # (head_dir, fname) = os.path.split(__file__)
        # head_fname = os.path.join(head_dir, 'RootDraw.h')
        # ROOT.gInterpreter.ProcessLine('#include ' + '\"' + head_fname + '\"')

    @classmethod
    def title_convert(cls, title):
        return title.replace("phi", "#varphi").replace("theta", "#theta")

    @classmethod
    def get_hist_transformed(cls, h1, x_offset=0, y_offset=0):
        h2 = h1.Clone()

        xmax = h2.GetXaxis().GetXmax()
        xmin = h2.GetXaxis().GetXmin()
        h2.GetXaxis().SetLimits(xmin + x_offset, xmax + x_offset)

        ymax = h2.GetYaxis().GetXmax()
        ymin = h2.GetYaxis().GetXmin()
        h2.GetYaxis().SetLimits(ymin + y_offset, ymax + y_offset)

        return h2

    @classmethod
    def get_hist_setRangeUser(cls, h1, xmin=None, xmax=None, ymin=None, ymax=None):
        if xmin is None:
            xmin = h1.GetXaxis().GetXmin()
        if xmax is None:
            xmax = h1.GetXaxis().GetXmax()
        if ymin is None:
            ymin = h1.GetYaxis().GetXmin()
        if ymax is None:
            ymax = h1.GetYaxis().GetXmax()

        h2 = h1.Clone()
        h2.GetXaxis().SetRangeUser(xmin, xmax)
        h2.GetYaxis().SetRangeUser(ymin, ymax)

        return h2

    @classmethod
    def set_hist_large_label(cls, h, label_size=0.05, title_size=0.06, title_offset=0.78,
                             xtitle_offset=None, ytitle_offset=None):
        h.GetXaxis().SetLabelSize(label_size)
        h.GetYaxis().SetLabelSize(label_size)
        h.GetXaxis().SetTitleSize(title_size)
        h.GetYaxis().SetTitleSize(title_size)
        try:
            h.GetZaxis().SetLabelSize(label_size)
            h.GetZaxis().SetTitleSize(title_size)
        except:
            pass

        if xtitle_offset is None:
            xtitle_offset = title_offset
        if ytitle_offset is None:
            ytitle_offset = title_offset
        h.GetXaxis().SetTitleOffset(xtitle_offset)
        h.GetYaxis().SetTitleOffset(ytitle_offset)

    @classmethod
    def set_obj_marker(cls, o, marker_size=2., line_width=2.,
                       marker_color=None, line_color=None,
                       marker_style=None, line_style=None):
        o.SetMarkerSize(marker_size)
        o.SetLineWidth(line_width)

        if marker_color is not None:
            o.SetMarkerColor(marker_color)
        if marker_style is not None:
            o.SetMarkerStyle(marker_style)
        if line_color is not None:
            o.SetLineColor(line_color)
        if line_style is not None:
            o.SetLineStyle(line_style)

    def get_obj_title(self, o):
        x_title = self.title_convert(o.GetXaxis().GetTitle())
        y_title = self.title_convert(o.GetYaxis().GetTitle())
        return x_title, y_title

    @classmethod
    def get_hist_ratio(cls, h1, h2):
        h_sig = h1.Clone()
        h_sig.Reset("ICES")
        h_sig.SetNameTitle(h1.GetName() + '_sig', '')

        xbins = h1.GetNbinsX()
        ybins = h1.GetNbinsY()

        for i_xbin in range(1, xbins + 1):
            for i_ybin in range(1, ybins + 1):
                i_bin = h1.GetBin(i_xbin, i_ybin)

                val1 = h1.GetBinContent(i_bin)
                val2 = h2.GetBinContent(i_bin)

                if val1 > 1e-10 and val2 > 1e-10:
                    err1 = h1.GetBinError(i_bin)
                    err2 = h2.GetBinError(i_bin)

                    sig = val1 / val2
                    err = TMath.Sqrt(err1 * err1 / val2 / val2 + val1 * val1 / val2 / val2 / val2 / val2 * err2 * err2)

                    h_sig.SetBinContent(i_bin, sig)
                    h_sig.SetBinError(i_bin, err)
                else:
                    h_sig.SetBinContent(i_bin, 1)
                    h_sig.SetBinError(i_bin, 0)

        return h_sig

    # @classmethod
    # def get_ratio_canvas(cls, hs, hs_ratio):
    #     gStyle.SetOptStat(0)
    #     cr = TCanvas("cr", "cr", 1024, 640)
    #     cr.SetFillStyle(4000)
    #
    #     nx = 1
    #     ny = 2
    #
    #     lmargin = 0.05
    #     rmargin = 0.12
    #     bmargin = 0.12
    #     tmargin = 0.12
    #
    #     ROOT.CanvasPartition(cr, nx, ny, lmargin, rmargin, bmargin, tmargin)
    #     lis_h = TList()
    #     lis_h.Add(hs_ratio)
    #     lis_h.Add(hs)
    #
    #     ROOT.canvas_multiple_pad(cr, nx, ny, lis_h)
    #
    #
    #     return 1

    def DrawSigPNG4Paper(self, fin_name, h_name, png_name, rot_phi,
                         lis_x_range=None, lis_y_range=None,
                         low_thr=-2., high_thr=2.):

        h_name = ROOTUtils.find_key_in_file(fin_name, h_name)
        fin = TFile(fin_name, 'read')
        hsig = fin.Get(h_name)
        x_title, y_title = self.get_obj_title(hsig)
        hsig.SetTitle(";%s;%s;" % (x_title, y_title))

        if lis_x_range is None:
            lis_x_range = [hsig.GetXaxis().GetXmin() + rot_phi, hsig.GetXaxis().GetXmax() + rot_phi]
        if lis_y_range is None:
            lis_y_range = [hsig.GetYaxis().GetXmin(), hsig.GetYaxis().GetXmax()]

        hcut1 = hsig.Clone()
        contours1 = array.array('d', [low_thr])
        hcut1.SetContour(1, contours1)
        hcut2 = hsig.Clone()
        contours2 = array.array('d', [high_thr])
        hcut2.SetContour(1, contours2)

        hsig = self.get_hist_transformed(hsig, x_offset=rot_phi)
        hcut1 = self.get_hist_transformed(hcut1, x_offset=rot_phi)
        hcut2 = self.get_hist_transformed(hcut2, x_offset=rot_phi)
        hsig = self.get_hist_setRangeUser(hsig, lis_x_range[0], lis_x_range[1], lis_y_range[0], lis_y_range[1])
        hcut1 = self.get_hist_setRangeUser(hcut1, lis_x_range[0], lis_x_range[1], lis_y_range[0], lis_y_range[1])
        hcut2 = self.get_hist_setRangeUser(hcut2, lis_x_range[0], lis_x_range[1], lis_y_range[0], lis_y_range[1])

        gStyle.SetOptStat(0)
        can = TCanvas("can", "", 1600, 1200)

        self.set_hist_large_label(hsig)
        hsig.Draw("colz")

        hcut1.SetLineColor(3)
        hcut1.SetLineWidth(5)
        hcut1.Draw("cont3same")
        hcut2.SetLineColor(2)
        hcut2.SetLineWidth(5)
        hcut2.Draw("cont3same")

        gPad.SetTickx()
        gPad.SetTicky()

        can.Print(png_name)

    def DrawRatioLeff4Paper(self, fin_meas_name, h_meas_name,
                            fin_sim_name, h_sim_name, png_name,
                            rot_phi, lis_x_range=None, lis_y_range=None,
                            lis_contours=None):
        h_meas_name = ROOTUtils.find_key_in_file(fin_meas_name, h_meas_name)

        fin_meas = TFile(fin_meas_name, 'read')
        h_meas = fin_meas.Get(h_meas_name)
        fin_pre = TFile(fin_sim_name, 'read')
        h_pre = fin_pre.Get(h_meas_name)

        if lis_x_range is None:
            lis_x_range = [h_meas.GetXaxis().GetXmin() + rot_phi, h_meas.GetXaxis().GetXmax() + rot_phi]
        if lis_y_range is None:
            lis_y_range = [h_meas.GetYaxis().GetXmin(), h_meas.GetYaxis().GetXmax()]
        if lis_contours is None:
            vmin = h_pre.GetMinimum()
            vmax = h_pre.GetMaximum()
            lis_contours = [(i + 0.5) * (vmax - vmin) + vmin for i in range(6)]

        h_meas = self.get_hist_transformed(h_meas, x_offset=rot_phi)
        h_pre = self.get_hist_transformed(h_pre, x_offset=rot_phi)
        h_meas = self.get_hist_setRangeUser(h_meas, lis_x_range[0], lis_x_range[1], lis_y_range[0], lis_y_range[1])
        h_pre = self.get_hist_setRangeUser(h_pre, lis_x_range[0], lis_x_range[1], lis_y_range[0], lis_y_range[1])

        h_meas.SetMaximum(h_pre.GetMaximum()*0.8)

        gStyle.SetOptStat(0)
        can = TCanvas("can", "", 1600, 1200)

        self.set_hist_large_label(h_meas)
        contours1 = array.array('d', lis_contours)
        h_pre.SetContour(len(lis_contours), contours1)
        self.set_hist_large_label(h_pre)
        h_pre.Draw("cont z list")
        can.Update()

        cont = gROOT.GetListOfSpecials().FindObject("contours")
        try:
            n_cont = cont.GetSize()
        except:
            print("*** No Contours Were Extracted!")
            n_cont = 0

        x_title, y_title = self.get_obj_title(h_meas)
        h_meas.SetTitle(";%s;%s;" % (x_title, y_title))

        h_meas.Draw("colz")
        gPad.SetTickx()
        gPad.SetTicky()

        lg = TLegend(0.01, 0.01, 0.99, 0.99)
        lg.SetNColumns(6)

        for i in range(n_cont):
            contLevel = cont.At(i)
            curv = contLevel.First()
            for j in range(contLevel.GetSize()):
                if j == 0:
                    lg.AddEntry(curv, "%.2f" % contours1[i], "l")
                curv.SetLineColor(2)
                curv.SetLineWidth(6)
                curv.SetLineStyle(i + 1)
                curv.Draw("same")

                curv = contLevel.After(curv)

        can.SetRightMargin(0.16)
        can.Print(png_name)

        can_lg = TCanvas("can", "", 1600, 100)
        lg.Draw()
        can_lg.Print(png_name.replace(".png", "_lg.png"))

    @classmethod
    def _get_thstack_or_multigraph(cls, name, title="", opt='h'):
        if opt.__contains__('h'):
            a = THStack(name, title)
        elif opt.__contains__('g'):
            a = TMultiGraph(name, title)
        else:
            raise Exception('input thstack or multigraph.')
        return a

    def TH1sComparison(self, lis_fname, lis_hname, png_name, lis_scale=None, lis_label=None, legend_columns=3,
                       xtitle='', ytitle='', h_or_g=True):
        if h_or_g:
            draw = THStack('draw', '')
            draw_ratio = THStack('draw_ratio', '')
        else:
            draw = TMultiGraph("draw", "")
            draw_ratio = TMultiGraph("draw_ratio", "")
        lg = TLegend(0.68, 0.68, 0.9, 0.9)
        lg.SetNColumns(legend_columns)

        h0 = None
        for i in range(len(lis_fname)):
            fname = lis_fname[i]
            hname = lis_hname[i]
            hname = ROOTUtils.find_key_in_file(fname, hname)

            fin = TFile(fname, 'read')
            h = fin.Get(hname)
            h.SetName('h%d' % i)
            h.SetDirectory(0)

            if i == 0:
                h0 = h

            # get scale and label for the hist
            if lis_scale is None:
                scale = 1. / h.Integral()
            else:
                scale = lis_scale[i]
            h.Scale(scale)

            x_title, y_title = self.get_obj_title(h)
            if lis_scale is not None:
                y_title = y_title + '/day'
            else:
                y_title = 'normalized ' + y_title
            if lis_label is None:
                label = os.path.splitext(os.path.split(fname)[1])[0]
            else:
                label = lis_label[i]

            if xtitle != '':
                x_title = self.title_convert(xtitle)
            if ytitle != '':
                y_title = self.title_convert(ytitle)
            draw.SetTitle(";%s;%s;" % (x_title, y_title))
            self.set_obj_marker(h, line_width=4, marker_color=i + 1, line_color=i + 1,
                                marker_style=1, line_style=i + 1)

            print('%s peak_area: %.5f' % (label, h.Integral()), file=open(png_name.replace('.png', '.log'), "a"))

            h_ratio = self.get_hist_ratio(h, h0)
            h_ratio.SetDirectory(0)
            draw_ratio.SetTitle(";%s;%s;" % (x_title, 'ratio'))
            self.set_obj_marker(h_ratio, line_width=4, marker_color=i + 1, line_color=i + 1,
                                marker_style=1, line_style=i + 1)
            if h_or_g:
                draw.Add(h, "histe")
                draw_ratio.Add(h_ratio, "e")
                lg.AddEntry(h, label, "l")
            else:
                g = TGraph(h)
                draw.Add(g)
                draw_ratio.Add(TGraph(h_ratio))
                lg.AddEntry(g, label, "l")

        gStyle.SetOptStat(0)
        can = TCanvas("can", "", 1600, 1200)
        can.Divide(1,2)
        can.cd(1)
        gPad.SetPad(0,0.35,1,1)
        gPad.SetBottomMargin(0.001)

        if h_or_g:
            draw.Draw("nostack")
        else:
            draw.Draw("ap")
        self.set_hist_large_label(draw,  label_size=0.07, title_size=0.07, xtitle_offset=0.8, ytitle_offset=0.6)

        lg.Draw()

        can.cd(2)
        gPad.SetPad(0,0,1,0.35)
        gPad.SetTopMargin(0.001)
        gPad.SetBottomMargin(0.3)
        if h_or_g:
            draw_ratio.Draw("nostack")
        else:
            draw_ratio.Draw("ap")
        self.set_hist_large_label(draw_ratio, label_size=0.14, title_size=0.14, xtitle_offset=0.8, ytitle_offset=0.3)
        gPad.Update()
        can.Print(png_name)

    def DrawMultiGraphs(self, lis_fname, lis_hname, png_name, lis_label=None, legend_columns=3,
                        xtitle='', ytitle=''):
        draw = TMultiGraph("draw", "")
        lg = TLegend(0.6, 0.7, 0.9, 0.9)
        lg.SetNColumns(legend_columns)

        for i in range(len(lis_fname)):
            fname = lis_fname[i]
            gname = lis_hname[i]
            gname = ROOTUtils.find_key_in_file(fname, gname)

            fin = TFile(fname, 'read')
            g = fin.Get(gname)
            g.SetName('g%d' % i)

            x_title, y_title = self.get_obj_title(g)
            if lis_label is None:
                label = os.path.splitext(os.path.split(fname)[1])[0]
            else:
                label = lis_label[i]

            if xtitle != '':
                x_title = self.title_convert(xtitle)
            if ytitle != '':
                y_title = self.title_convert(ytitle)
            draw.SetTitle(";%s;%s;" % (x_title, y_title))
            self.set_obj_marker(g, marker_size=1, line_width=2,
                                marker_color=int(i/4)+1, line_color=int(i/4)+1,
                                marker_style=7, line_style=1+i%4)

            g = TGraph(g)
            draw.Add(g)
            lg.AddEntry(g, label, "pl")

        can = TCanvas("can", "", 1600, 1200)

        draw.Draw("al")
        # self.set_hist_large_label(draw)
        draw.GetXaxis().SetRangeUser(0, 100)
        draw.GetYaxis().SetRangeUser(0.001, 1)
        can.SetLeftMargin(0.1)
        can.SetBottomMargin(0.11)
        gPad.SetLogy()

        lg.Draw()
        can.Print(png_name)
