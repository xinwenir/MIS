import random
import os.path
from array import array

import numpy as np
import pandas as pd

from DataFactory import CommonUtils

from ROOT import TFile, TMath, TString, TCanvas, TGraph2D, TGraphErrors, TGraph2DErrors, TH1F, TH2F, TVector3, TTree, \
    RDataFrame, RDF, TVectorD
from ROOT import gStyle

import logging

logger = logging.getLogger(__name__)


def find_key_in_file(fin_name, part_name: str, opt='part', lis_key_parts=None):
    t_fin = TFile(fin_name, 'read')
    key_name_list = [i.GetName() for i in t_fin.GetListOfKeys()]
    key_name = ''
    lis_keys = []
    for i_key in key_name_list:
        if opt.__contains__('part') and i_key.__contains__(part_name):
            key_name = i_key
            lis_keys.append(i_key)
        elif opt.__contains__('whole') and i_key == part_name:
            key_name = part_name
    t_fin.Close()

    if len(lis_keys) <= 1:
        return key_name

    if lis_key_parts is not None:
        n_contain = np.zeros(len(lis_keys))
        for i in range(len(lis_keys)):
            i_key = lis_keys[i]
            for j in lis_key_parts:
                if i_key.__contains__(j):
                    n_contain[i] += 1
        contain_max = n_contain.max(axis=0)
        remove_ind = np.where(n_contain == contain_max)[0]
        t_lis_key = []
        for ind in remove_ind:
            t_lis_key.append(lis_keys[ind])
        lis_keys = t_lis_key

    shortest_index = np.array([len(i) for i in lis_keys]).argmin()
    t_fin.Close()

    return lis_keys[shortest_index]


def find_all_keys_in_file(fin_name, part_name: str):
    t_fin = TFile(fin_name, 'read')
    lis_keys = []
    for i_key in t_fin.GetListOfKeys():
        if i_key.GetName().__contains__(part_name):
            lis_keys.append(i_key.GetName())
    t_fin.Close()

    return lis_keys


def get_bin_position(axis, ibin, opt="center"):
    up = axis.GetBinUpEdge(ibin)
    low = axis.GetBinLowEdge(ibin)
    if opt.lower().__contains__("c"):
        pos = axis.GetBinCenter(ibin)
    elif opt.lower().__contains__("up"):
        pos = up
    elif opt.lower().__contains__("low"):
        pos = low
    elif opt.lower().__contains__("nd"):
        pos = random.uniform(low, up)
    else:
        raise Exception("opt = %s. Please give opt = center or up or low or random" % opt)

    return pos


def check_two_hists_the_same_size(h1, h2):
    if h1.GetNbinsX() != h2.GetNbinsX():
        return False
    if h1.GetNbinsY() != h2.GetNbinsY():
        return False
    if abs(h1.GetXaxis().GetBinWidth(1) - h2.GetXaxis().GetBinWidth(1)) > 1e-5:
        return False
    if abs(h1.GetYaxis().GetBinWidth(1) - h2.GetYaxis().GetBinWidth(1)) > 1e-5:
        return False
    return True


def fill_array_into_h1f(arr, h_name, h_title, _bin_width=1, _xmin=1, _xmax=-1):
    tmp_array = array('d', arr)  # for pyroot
    tmp_weight = array('d', np.ones(len(arr)))
    if _xmin >= _xmax:
        _xmin = arr.min()
        _xmax = arr.max()
    nbin = int(float(_xmax - _xmin) / _bin_width + 0.99) + 1
    h1f = TH1F(h_name, h_title, nbin, _xmin, _xmin + nbin * _bin_width)
    h1f.FillN(len(arr), tmp_array, tmp_weight)
    return h1f


def add_array_to_h1f(arr, h1f):
    tmp_array = array('d', arr)  # for pyroot
    tmp_weight = array('d', np.ones(len(arr)))
    h1f.FillN(len(arr), tmp_array, tmp_weight)
    return h1f


def cvt_tree_to_dataframe(tree: TTree):
    rdf = RDataFrame(tree)
    dic_numpy = rdf.AsNumpy()
    df = pd.DataFrame(dic_numpy)
    return df


def cvt_dataframe_to_tree(df: pd.DataFrame, tree_name, fout_name):
    data = {key: np.array(df[key].values) for key in df.columns}
    try:
        rdf = RDF.MakeNumpyDataFrame(data)
    except:
        try:
            rdf = RDF.FromNumpy(data)
        except:
            raise Exception('cvt_dataframe_to_tree: fail to use RDF to convert dataframe to tree.')

    rdf.Snapshot(tree_name, fout_name)


def hist2csv(fin_name: str, hist_name: str, fout_name: str):
    h_name = find_key_in_file(fin_name, hist_name)
    if h_name == '':
        raise Exception('Cannot find %s in %s' % (hist_name, fin_name))

    fin = TFile(fin_name, 'read')
    h = fin.Get(h_name)

    xbins = h.GetNbinsX()
    ybins = h.GetNbinsY()
    zbins = h.GetNbinsZ()

    print('%s shape: (%d,%d,%d)' % (h_name, xbins, ybins, zbins))

    if not os.path.exists(fout_name):
        print('ix,iy,iz,x,y,z,val,err', file=open(fout_name, 'a'))

    for i_xbin in range(1, xbins + 1):
        x = h.GetXaxis().GetBinCenter(i_xbin)
        for i_ybin in range(1, ybins + 1):
            y = h.GetYaxis().GetBinCenter(i_ybin)
            for i_zbin in range(1, zbins + 1):
                z = h.GetZaxis().GetBinCenter(i_zbin)

                ibin = h.GetBin(i_xbin, i_ybin, i_zbin)

                val = h.GetBinContent(ibin)
                err = h.GetBinError(ibin)

                print('%d,%d,%d,%.5f,%.5f,%.5f,%.5f,%.5f' % (i_xbin - 1, i_ybin - 1, i_zbin - 1, x, y, z, val, err),
                      file=open(fout_name, 'a'))

    print('printed %s into %s' % (h_name, fout_name))


def draw_graph_with_csv(fin_name: str, lis_columns: list, fout_name: str, graph_name='', fit_func=None):
    data = pd.read_csv(fin_name, index_col=False)
    if data.shape[1] < 2:
        data = pd.read_csv(fin_name, index_col=False, sep=' ')
    fin_name = os.path.splitext(os.path.split(fin_name)[1])[0]

    data = data[lis_columns]
    lis_val_col = []
    lis_err_col = []
    for i in lis_columns:
        if not i.lower().__contains__('err'):
            lis_val_col.append(i)
        else:
            lis_err_col.append(i)

    if len(lis_val_col) == 2:
        x_title = lis_columns[0]
        y_title = lis_columns[1]
        print('Draw a 2D graph: %s-%s' % (x_title, y_title))

        if graph_name == '':
            graph_name = 'gr_%s_%s_%s' % (x_title, y_title, fin_name)
        graph_title = graph_name + '; %s; %s' % (x_title, y_title)

        n_data = data.shape[0]
        x = array('d', data[x_title].to_numpy())
        y = array('d', data[y_title].to_numpy())
        x_err = array('d', n_data * [0.])
        y_err = array('d', n_data * [0.])

        for i in lis_err_col:
            if i.__contains__(x_title):
                x_err = array('d', data[i].to_numpy())
            elif i.__contains__(y_title):
                y_err = array('d', data[i].to_numpy())
            else:
                logger.warning('Cannot identify of what the error %s is.' % i)

        gr = TGraphErrors(n_data, x, y, x_err, y_err)
    elif len(lis_val_col) == 3:
        x_title = lis_columns[0]
        y_title = lis_columns[1]
        z_title = lis_columns[2]
        print('Draw a 3D graph: %s-%s-%s' % (x_title, y_title, z_title))

        if graph_name == '':
            graph_name = 'gr_%s_%s_%s_%s' % (x_title, y_title, z_title, fin_name)
        graph_title = graph_name + '; %s; %s; %s' % (x_title, y_title, z_title)

        n_data = data.shape[0]
        x = array('d', data[x_title].to_numpy())
        y = array('d', data[y_title].to_numpy())
        z = array('d', data[z_title].to_numpy())
        x_err = array('d', n_data * [0.])
        y_err = array('d', n_data * [0.])
        z_err = array('d', n_data * [0.])

        for i in lis_err_col:
            if i.__contains__(x_title):
                x_err = array('d', data[i].to_numpy())
            elif i.__contains__(y_title):
                y_err = array('d', data[i].to_numpy())
            elif i.__contains__(z_title):
                z_err = array('d', data[i].to_numpy())
            else:
                logger.warning('Cannot identify of what the error %s is.' % i)

        gr = TGraph2DErrors(n_data, x, y, z, x_err, y_err, z_err)
    else:
        raise Exception('Cannot identify what graph to draw.\n'
                        'values are %s, errors are %s' % (lis_val_col.__str__(), lis_err_col.__str__()))

    gr.SetName(graph_name)
    gr.SetTitle(graph_title)
    if fit_func is not None:
        gStyle.SetOptFit(1)
        gr.Fit(fit_func)
        gr.Draw('ap')

    fout = TFile(fout_name, 'update')
    gr.Write()
    fout.Save()
    fout.Close()


def cut_hist_with_threshold(h1, h_name, h_title, low_thr, high_thr):
    h2 = h1.Clone()
    h2.Reset("ICES")
    h2.SetName(h_name)
    h2.SetTitle(h_title)

    xbins = h1.GetNbinsX()
    ybins = h1.GetNbinsY()

    for i in range(1, xbins + 1):
        for j in range(1, ybins + 1):
            ibin = h1.GetBin(i, j)
            c = h1.GetBinContent(ibin)
            if low_thr < c < high_thr:
                h2.SetBinContent(ibin, c)

    return h2


def classify_by_board_channel_ID(fin_name):
    fin = TFile(fin_name, 'read')

    dic_hname = {}

    for i in fin.GetListOfKeys():
        hname = i.GetName()
        lis_hname = hname.split('_')
        bid = -1
        cid = -1
        for j in lis_hname:
            if j.__contains__('board'):
                bid = float(j.replace('board', ''))
            if j.__contains__('chn'):
                cid = float(j.replace('chn', ''))

        if bid > -1 and cid > -1:
            key_name = 'board%d_chn%d' % (bid, cid)
            if dic_hname.get(key_name) is None:
                dic_hname[key_name] = [hname]
            else:
                dic_hname[key_name].append(hname)

    fin.Close()

    return dic_hname


def merge_root_files(lis_root_files, out_root_file, if_delete=True):
    print('Have %d sub-files' % len(lis_root_files))
    # remove existed root file and hadd all split into it
    if os.path.exists(out_root_file):
        os.remove(out_root_file)
    if not out_root_file.__contains__(".root"):
        out_root_file += ".root"

    str_cmd = 'hadd -f %s ' % out_root_file
    split_files = ''
    for i in lis_root_files:
        split_files += i + ' '
    str_cmd += split_files

    os.system(str_cmd)
    if if_delete:
        os.system('rm %s' % split_files)


def convert_hist_to_bool_with_selection(h, sel_exp):
    h_bool = h.Clone()
    h_bool.Reset("ICES")

    lis_sel = sel_exp.split("&&")
    sel_exp_new = " && ".join(["val%s" % i for i in lis_sel])

    nbins = h.GetNbinsX() * h.GetNbinsY() * h.GetNbinsZ()

    for i in range(1, nbins + 1):
        val = h.GetBinContent(i)
        if sel_exp_new == "val":
            h_bool.SetBinContent(i, 1)
            continue
        if eval(sel_exp_new):
            h_bool.SetBinContent(i, 1)
        else:
            h_bool.SetBinContent(i, 0)

    return h_bool


def select_events_with_hist2d_selection(tree, x_name, y_name, h2d, sel_exp):
    h_bool = convert_hist_to_bool_with_selection(h2d, sel_exp)

    tr = tree.CloneTree(0)

    x = array('d', [0.])
    y = array('d', [0.])
    tree.SetBranchAddress(x_name, x)
    tree.SetBranchAddress(y_name, y)

    for i in range(tree.GetEntries()):
        tr.GetEntry(i)

        ibin = h_bool.FindBin(x[0], y[0])
        val = h_bool.GetBinContent(ibin)

        if val > 1e-10:
            tr.Fill()

    return tr


def merge_meas_time_in_tfile(fname):
    fin = TFile(fname, "update")
    fin.ls()

    lis_meas_time = []
    for i in fin.GetListOfKeys():
        if i.GetName() == "meas_time":
            val = fin.Get(i.GetName())
            lis_meas_time.append(val.Norm1())
            i.Delete()

    meas_time_total = sum(lis_meas_time)

    # add the measure time information
    measTime = TVectorD(1)
    measTime[0] = meas_time_total
    measTime.Write('meas_time')
    fin.Save()
    fin.Close()


def _get_h2d_with_lis_h1d(lis_h1d, lis_h1d_label, h2d_name):
    ybins = len(lis_h1d)
    xbins = lis_h1d[0].GetNbinsX()

    htitle = CommonUtils.find_same_prefix([i.GetTitle() for i in lis_h1d])

    h2d = TH2F(h2d_name, htitle, xbins, 0, xbins, ybins, 0, ybins)
    h2d_norm = TH2F(h2d_name + "_normalized", htitle + "(normalized)", xbins, 0, xbins, ybins, 0, ybins)

    for j in range(1, ybins + 1):
        h1d = lis_h1d[j-1]
        h1d_label = lis_h1d_label[j-1]
        for i in range(xbins):
            norm_scale = h1d.GetBinContent(xbins-i)
            if norm_scale > 1e-10:
                break
        for i in range(1, xbins + 1):
            # set content
            ibin = h2d.GetBin(i, j)
            val = h1d.GetBinContent(i)
            h2d.SetBinContent(ibin, val)
            h2d_norm.SetBinContent(ibin, val / norm_scale)

            #set y label
            h2d.GetYaxis().SetBinLabel(j, h1d_label)
            h2d_norm.GetYaxis().SetBinLabel(j, h1d_label)

    # set x label
    for i in range(1, xbins + 1):
        h2d.GetXaxis().SetBinLabel(i, lis_h1d[0].GetXaxis().GetBinLabel(i))
        h2d_norm.GetXaxis().SetBinLabel(i, lis_h1d[0].GetXaxis().GetBinLabel(i))

    return h2d, h2d_norm
