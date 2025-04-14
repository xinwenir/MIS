import os.path

import pandas as pd
import numpy as np

from DataFactory import CommonUtils, ROOTUtils
from DataAnalysisTools import baseCore, _Index, _barIndex
from ROOT import TH1F, TF1, TCanvas, TFile, TSpectrum, TLine, TTree, TMath, TVector3, TList, TVectorD, gPad
import time
import logging

logger = logging.getLogger(__name__)


def get_n_boards(data: pd.DataFrame):
    n = data.drop_duplicates(['boardID'])['boardID'].max() + 1
    return n


def cut_dataframe(data: pd.DataFrame, time_cut=-1, num_cut=-1):
    if num_cut > 0:
        data = data[:num_cut]
    if time_cut > 0:
        data = data[data['timeTag'] < (data['timeTag'].iloc[0] + time_cut)]
    return data


def select_hit_all_boards(data: pd.DataFrame):
    n_board = get_n_boards(data)
    trigger_count = data.groupby(['triggerID'])['triggerID'].count()
    hit_all_filter = trigger_count.index.values[trigger_count.values == n_board]
    selected_data = data.set_index('triggerID').loc[hit_all_filter, :]
    return selected_data


def getData_with_cut(h5_obj: baseCore.h5Data, num_cut=-1, time_cut=-1):
    n_split = h5_obj.get_n_splits()
    if num_cut < 0 and time_cut < 0:
        return h5_obj.getData(-2)
    elif num_cut < 0 <= time_cut:
        num_cut = int(1e15)
    elif num_cut >= 0 > time_cut:
        time_cut = int(1e15)
    _data = h5_obj.getData(0)
    i = 1
    while _data.iloc[-1, :]['timeTag'] < (_data.iloc[0, :]['timeTag'] + time_cut) and len(_data) < num_cut:
        t_data = h5_obj.getData(i)
        t_data['triggerID'] = t_data['triggerID'] + i * 65535
        _data = pd.concat([_data, t_data], ignore_index=True)
        i = i + 1
        if i >= n_split:
            break
    _data = cut_dataframe(_data, time_cut=time_cut, num_cut=num_cut)
    return _data


def get_channel_replaced_data(data, baseline_mtx, threshold_mtx, bar_channel_config, n_board=4, n_channel=36, n_bar=32):
    _data_n_channels = data.loc[:, _Index[:n_channel]].values
    n_events = int(_data_n_channels.shape[0] / n_board)
    _data_n_channels = _data_n_channels.reshape(n_events, n_board, -1)

    # larger than threshold
    np.place(_data_n_channels, _data_n_channels < threshold_mtx, 0)
    # substrate baseline
    _data_n_channels = _data_n_channels - baseline_mtx
    np.place(_data_n_channels, _data_n_channels < 0, 0)

    # connect the channels
    _data_mtx = np.full([n_events, n_board, n_bar], 0)
    for i_board in range(n_board):
        for i_bar in range(n_bar):
            t_chn = bar_channel_config.iloc[i_board, i_bar]
            _data_mtx[:, i_board, i_bar] = _data_n_channels[:, i_board, t_chn]

    return _data_mtx


def calc_one_event_4plane_det(_data, det_info):
    n_board = _data.shape[0]
    if not n_board == 4:
        _data = _data[0:3, :]
        n_board = 4
        logger.warning('Attention: Only use the first four board information.')

    event_mtx = np.full([n_board, 4], -1.)

    # calc for each two board
    lis_x_board = []
    lis_y_board = []
    for i_board in range(n_board):
        if str(det_info['xy_board'][i_board]).__contains__('y') or int(det_info['xy_board'][i_board]) == 0:
            lis_y_board.append(i_board)
        elif str(det_info['xy_board'][i_board]).__contains__('x') or int(det_info['xy_board'][i_board]) == 1:
            lis_x_board.append(i_board)
        else:
            logger.error('In DetInfo, xy_board is not x nor y.')
            raise ValueError('In DetInfo, xy_board is not x nor y.')
    if not len(lis_x_board) == 2 or not len(lis_x_board) == 2:
        logger.error('In DetInfo, you should give two x board and two y board.\nNumber of x_board = %d; number of '
                     'y_board = %d' % (len(lis_x_board) == 2, len(lis_x_board) == 2))
        raise ValueError('In DetInfo, you should give two x board and two y board.')

    i_board1 = lis_x_board[0]
    i_board2 = lis_x_board[1]
    lis1, lis2 = calc_two_board_plane_det(_data, i_board1, i_board2, det_info)
    event_mtx[0, :] = np.array(lis1)
    event_mtx[1, :] = np.array(lis2)

    i_board1 = lis_y_board[0]
    i_board2 = lis_y_board[1]
    lis1, lis2 = calc_two_board_plane_det(_data, i_board1, i_board2, det_info)
    event_mtx[2, :] = np.array(lis1)
    event_mtx[3, :] = np.array(lis2)

    # x1: event_mtx[0,:]; x2: event_mtx[0,:]; y1: event_mtx[0,:]; y2: event_mtx[0,:]
    return event_mtx


class RawProcessing:
    def __init__(self, det, global_dir_tag):
        self._det = det
        self._global_dir_tag = global_dir_tag
        self._can = TCanvas('can_rawpro_%s' % self._det.get_det_name(),
                            'can_rawpro_%s' % self._det.get_det_name(), 2000, 1500)
        self._pdf_fname = os.path.join(global_dir_tag[0], 'rawpro_%s.pdf' % self._det.get_det_name())
        self._root_fname = os.path.join(global_dir_tag[0], 'rawpro_%s.root' % self._det.get_det_name())

        self.pro_data_frame = None

        self._baseline_csv = None
        self._threshold_csv = None
        self._baseline_mtx = None
        self._threshold_mtx = None
        self._fout_occ = None
        self._can_occ = None

        self._tree_list = TList()

    def find_baseline(self):
        if self._baseline_csv is None:
            baseline_csv_fname = self._det.get_baseline_csv_fname()
            if os.path.isfile(baseline_csv_fname):
                logger.info('Get baseline csv from file: %s' % baseline_csv_fname)
                self._baseline_csv = pd.read_csv(baseline_csv_fname)
            else:
                baseline_h5_fname = self._det.get_baseline_h5_fname()
                logger.info('Calculate baseline from h5 data: %s' % baseline_h5_fname)
                self._calc_baseline(baseline_h5_fname=baseline_h5_fname, baseline_csv_fname=baseline_csv_fname)

    def _calc_baseline(self, baseline_h5_fname, baseline_csv_fname):
        # get data from h5 file, cut according to the config file
        h5_obj = baseCore.h5Data(baseline_h5_fname, 'r')
        num_cut = self._det.get_baseline_num_cut()
        time_cut = self._det.get_baseline_time_cut()
        baseline_data = getData_with_cut(h5_obj, num_cut=num_cut, time_cut=time_cut)
        baseline_data = select_hit_all_boards(baseline_data)

        # find the baseline for every board every channel
        n_board = get_n_boards(baseline_data)
        n_channel = self._det.get_n_channels()

        baseline_mtx = np.full([n_board, n_channel], 0)

        fname_prefix = baseline_h5_fname.replace('.h5', '_baseline')
        fout = TFile('%s.root' % fname_prefix, 'recreate')
        fout.cd()

        for i_board in range(n_board):
            board_data = baseline_data[baseline_data['boardID'] == i_board]
            for i_channel in range(n_channel):
                # get the certain data
                channel_data = board_data.iloc[:, i_channel]
                tmp_arr = channel_data[channel_data.values > 0].to_numpy()
                if tmp_arr.shape[0] < 10:
                    continue

                # fill data into hist
                h_name = 'spectrum_board%d_channel%d' % (i_board, i_channel)
                h_title = 'Bl_' + h_name + '; Energy channel'
                tmp_h1f = ROOTUtils.fill_array_into_h1f(tmp_arr, h_name, h_title)

                # fine peak
                peak_pos = int(tmp_h1f.GetBinCenter(tmp_h1f.GetMaximumBin()))
                peak_height = tmp_h1f.GetBinContent(tmp_h1f.GetMaximumBin())

                self._can.cd()
                tmp_h1f.Draw('hist')
                line = TLine(peak_pos, 0, peak_pos, peak_height)
                line.SetLineColor(2)
                line.Draw()
                self._can.SetName('can_spectrum_board%d_channel%d' % (i_board, i_channel))
                self._can.Write()
                self._can.Print('%s.pdf(' % fname_prefix)

                baseline_mtx[i_board, i_channel] = peak_pos

        self._can.Print('%s.pdf]' % fname_prefix)
        self._baseline_csv = pd.DataFrame(baseline_mtx.astype('int'), columns=_Index[:n_channel])
        self._baseline_csv.to_csv(baseline_csv_fname, index=False)

        fout.Save()
        fout.Close()

    def data_time_slice(self):
        # global bin_width
        logger.info('slice time for detector %s' % self._det.get_det_name())
        measData_h5_fname = self._det.get_measData_h5_fname()
        fname_prefix = measData_h5_fname.replace('.h5', '_time_slice')
        pdf_fname = fname_prefix + '.pdf'
        root_fname = fname_prefix + '.root'
        fout = TFile(root_fname, 'recreate')
        fout.cd()
        can = TCanvas('can_time_slice', 'can_time_slice', 2000, 1500)

        h5_obj = baseCore.h5Data(measData_h5_fname, 'r')
        time_interval = self._det.get_measData_time_interval()
        n_split = h5_obj.get_n_splits()
        _data0 = h5_obj.getData(0)
        start_time = _data0['timeTag'][0]
        n_board = get_n_boards(_data0)
        n_channel = self._det.get_n_channels()

        # slice every 20 splits
        i = 0
        while i < n_split:
            from_split = i
            _data = h5_obj.getData(i)
            i = i + 1
            print('data_time_slice ', i)
            for j in range(1, 20):
                t_data = h5_obj.getData(i)
                t_data['triggerID'] = t_data['triggerID'] + j * 65535
                _data = pd.concat([_data, t_data], ignore_index=True)
                i = i + 1
                if i >= n_split:
                    break

            _data = select_hit_all_boards(_data)

            for i_board in range(n_board):
                board_data = _data[_data['boardID'] == i_board]
                for i_channel in range(n_channel):
                    channel_data = board_data.loc[:, ['chn_%s' % i_channel, 'timeTag']]
                    _time_arr = channel_data[channel_data.iloc[:, 0] > 0]['timeTag'].values
                    _time_arr = _time_arr - start_time
                    h_name = 'time_slice_board%d_channel%d_f%dt%d' % (i_board, i_channel, from_split, i - 1)
                    h_title = 'time_slice_board%d_channel%d_startingAt%s_split%d_to_split%d' \
                              % (i_board, i_channel, time.strftime('%D', time.localtime(start_time)),
                                 from_split, i - 1)
                    if len(_time_arr) <= 0:
                        time_min = channel_data['timeTag'].iloc[0] - start_time
                        time_max = channel_data['timeTag'].iloc[-1] - start_time
                        tmp_h1f = TH1F(h_name, h_title, 100, time_min, time_max)
                    else:
                        time_min = channel_data['timeTag'].iloc[0] - start_time
                        time_max = channel_data['timeTag'].iloc[-1] - start_time
                        bin_width = time_interval
                        if (time_max - time_min) < time_interval:
                            bin_width = int(time_max - time_min)
                        tmp_h1f = ROOTUtils.fill_array_into_h1f(_time_arr, h_name, h_title, _bin_width=bin_width,
                                                            _xmin=time_min, _xmax=time_max)

                        if len(_time_arr) > 10:
                            peak_height = tmp_h1f.GetBinContent(tmp_h1f.GetMaximumBin())
                            minima_height = tmp_h1f.GetBinContent(tmp_h1f.GetMinimumBin())
                            if float(peak_height) / minima_height > 2 \
                                    and tmp_h1f.GetMinimumBin() < tmp_h1f.GetNbinsX():
                                logger.warning('Attention: board%d_channel%d may encounters a break at %.1f s'
                                               % (i_board, i_channel, tmp_h1f.GetBinCenter(tmp_h1f.GetMinimumBin())))

                    for i_bin in range(1, tmp_h1f.GetNbinsX() + 1, 3):
                        str_time = time.strftime('%H:%M:%S', time.localtime(int(tmp_h1f.GetBinCenter(i_bin))))
                        tmp_h1f.GetXaxis().SetBinLabel(i_bin, str_time)

                    can.cd()
                    tmp_h1f.Draw('hist')
                    can.Print('%s(' % pdf_fname)
                    can.Write()

        can.Print('%s]' % pdf_fname)
        fout.Save()
        fout.Close()

    def find_threshold(self):
        if self._threshold_csv is None:
            threshold_csv_fname = self._det.get_threshold_csv_fname()
            if os.path.isfile(threshold_csv_fname):
                logger.info('Get threshold csv from file: %s' % threshold_csv_fname)
                self._threshold_csv = pd.read_csv(threshold_csv_fname)
            else:
                threshold_h5_fname = self._det.get_measData_h5_fname()
                logger.info('Calculate threshold from h5 data: %s' % threshold_h5_fname)
                self._calc_threshold(threshold_h5_fname=threshold_h5_fname, threshold_csv_fname=threshold_csv_fname)

    def _calc_threshold(self, threshold_h5_fname, threshold_csv_fname):
        if self._baseline_csv is None:
            self.find_baseline()
        # get data from h5 file, cut according to the config file
        h5_obj = baseCore.h5Data(threshold_h5_fname, 'r')
        num_cut = self._det.get_measData_num_cut()
        if num_cut < 0:
            num_cut = 1000000
        time_cut = self._det.get_measData_time_cut()
        threshold_data = getData_with_cut(h5_obj, num_cut=num_cut, time_cut=time_cut)
        threshold_data = select_hit_all_boards(threshold_data)

        # find the threshold for every board every channel
        n_board = get_n_boards(threshold_data)
        n_channel = self._det.get_n_channels()

        threshold_mtx = np.full([n_board, n_channel], 0)

        fname_prefix = threshold_h5_fname.replace('.h5', '_threshold')
        fout = TFile('%s.root' % fname_prefix, 'recreate')
        fout.cd()

        for i_board in range(n_board):
            board_data = threshold_data[threshold_data['boardID'] == i_board]
            for i_channel in range(n_channel):
                # get the certain data
                channel_data = board_data.iloc[:, i_channel]
                tmp_arr = channel_data[channel_data.values > 0].to_numpy()
                if tmp_arr.shape[0] < 10:
                    continue

                # fill data into hist
                h_name = 'spectrum_board%d_channel%d' % (i_board, i_channel)
                h_title = 'Tr_' + h_name + '; Energy channel'
                tmp_h1f = ROOTUtils.fill_array_into_h1f(tmp_arr, h_name, h_title)
                tmp_h1f.Write()

                # fine the minima
                # tmp_h1f.Smooth()
                t_baseline = self._baseline_csv.iloc[i_board, i_channel]
                baseline_bin = tmp_h1f.FindBin(t_baseline)

                tmp_h1f.GetXaxis().SetRange(baseline_bin, baseline_bin + 20)
                threshold_pos = int(tmp_h1f.GetBinCenter(tmp_h1f.GetMinimumBin()))
                y_axis_max = int(tmp_h1f.GetBinContent(tmp_h1f.GetMinimumBin())) * 15

                self._can.cd()
                tmp_h1f.GetXaxis().SetRange(baseline_bin - 30, baseline_bin + 300)
                tmp_h1f.SetMaximum(int(y_axis_max) + 1)
                tmp_h1f.Draw('hist')
                line_baseline = TLine(t_baseline, 0, t_baseline, y_axis_max)
                line_threshold = TLine(threshold_pos, 0, threshold_pos, y_axis_max)
                line_baseline.SetLineColor(2)
                line_threshold.SetLineColor(1)
                line_baseline.Draw()
                line_threshold.Draw()
                self._can.Write()
                self._can.Print('%s.pdf(' % fname_prefix)

                threshold_mtx[i_board, i_channel] = threshold_pos

        fout.Save()
        fout.Close()

        self._can.Print('%s.pdf]' % fname_prefix)
        self._threshold_csv = pd.DataFrame(threshold_mtx.astype('int'), columns=_Index[:n_channel])
        self._threshold_csv.to_csv(threshold_csv_fname, index=False)

    def _get_channel_replace_information(self, channel_replace_config_fname=''):
        if os.path.isfile(channel_replace_config_fname):
            self.bar_channel_config = pd.read_csv(channel_replace_config_fname)
        else:
            if self._threshold_csv is None:
                self.find_baseline()
                self.find_threshold()
            n_board = self._threshold_csv.shape[0]
            n_channel = self._det.get_n_channels()
            n_bar = self._det.get_n_bars()
            channel_replace_mtx = np.full([n_board, n_bar], -1)

            for i_board in range(n_board):
                n_replace = 0
                for i_bar in range(n_bar):
                    if self._threshold_csv.iloc[i_board, i_bar] <= 0:
                        if (n_bar + n_replace) >= n_channel:
                            logger.error('Error with auto-finding channel replacing information.')
                            raise Exception('Error with auto-finding channel replacing information.\nFollow the '
                                            'replacing rule or write a config file.')
                        elif (n_bar + n_replace) < n_channel and \
                                self._threshold_csv.iloc[i_board, n_bar + n_replace] <= 0:
                            logger.error('Error with auto-finding channel replacing information.')
                            raise Exception('Error with auto-finding channel replacing information.\nFollow the '
                                            'replacing rule or write a config file.')
                        channel_replace_mtx[i_board, i_bar] = n_bar + n_replace
                        n_replace = n_replace + 1
                    else:
                        channel_replace_mtx[i_board, i_bar] = i_bar
            channel_replace_csv = pd.DataFrame(channel_replace_mtx.astype('int'), columns=_barIndex[:n_bar])
            channel_replace_config_outfname = os.path.join(self._global_dir_tag[0],
                                                           'ChannelReplaceConfig_%s.csv' % self._det.get_det_name())
            channel_replace_csv.to_csv(channel_replace_config_outfname, index=False)
            self.bar_channel_config = channel_replace_csv

    def _replace_baseline_threshold(self):
        n_board = self._threshold_csv.shape[0]
        n_bar = self._det.get_n_bars()
        baseline_mtx = np.full([n_board, n_bar], -1)
        threshold_mtx = np.full([n_board, n_bar], -1)

        for i_board in range(n_board):
            for i_bar in range(n_bar):
                i_channel = self.bar_channel_config.iloc[i_board, i_bar]
                baseline_mtx[i_board, i_bar] = self._baseline_csv.iloc[i_board, i_channel]
                threshold_mtx[i_board, i_bar] = self._threshold_csv.iloc[i_board, i_channel]

        return baseline_mtx, threshold_mtx

    def convert_txt_to_tree(self):
        t_tree = TTree('t_tree', 't_tree')
        _tree_fname = self._det.get_tree_fname()

        txt_fname = self._det.get_txt_fname()

        det_info_fname = self._det.get_det_pro_fname()
        det_info = pd.read_csv(det_info_fname)
        z_aver1 = (det_info['Origin_z'][0] + det_info['Origin_z'][1]) / 2
        z_aver2 = (det_info['Origin_z'][2] + det_info['Origin_z'][3]) / 2
        rotAnglex = self._det.get_rotAngle()[0]
        rotAngley = self._det.get_rotAngle()[1]

        meas_time = self._det.get_meas_time()

        str_txt2root = "root -q -b -l './DataAnalysisTools/txt2root.C(\"%s\",\"%s\",%.3f,%.3f,%.3f, %.3f, %d)' " \
                       % (txt_fname, _tree_fname, z_aver1, z_aver2, rotAnglex, rotAngley, meas_time)
        logger.info(str_txt2root)
        os.system(str_txt2root)

    def convert_h5_to_tree(self):
        # get baseline and threshold for the following use
        self.find_baseline()
        self.find_threshold()
        self._get_channel_replace_information(self._det.get_channel_config_fname())
        threshold_mtx = self._threshold_csv.values
        baseline_mtx = self._baseline_csv.values

        # get properties
        _h5_fname = self._det.get_measData_h5_fname()

        h5_obj = baseCore.h5Data(_h5_fname, 'r')
        num_cut = self._det.get_measData_num_cut()
        time_cut = self._det.get_measData_time_cut()
        _data0 = h5_obj.getData(0)
        start_time = _data0['timeTag'][0]
        end_time = -1
        n_split = h5_obj.get_n_splits()
        n_board = get_n_boards(_data0)
        n_bar = self._det.get_n_bars()
        n_channel = self._det.get_n_channels()

        if num_cut < 0:
            num_cut = int(1e15)
        if time_cut < 0:
            time_cut = int(1e15)

        n_split_used = 0
        for i_split in range(n_split):
            print('Split %d ...' % i_split)
            _data = h5_obj.getData(i_split)
            _data['triggerID'] = _data['triggerID'] + i_split * 65535

            time_last_event = end_time - start_time
            num_last_event = _data.index[-1]
            b_cut = False
            if time_last_event > time_cut or num_last_event > num_cut:
                b_cut = True
                _data = cut_dataframe(_data, time_cut=time_cut, num_cut=num_cut)
                left_event = int(len(_data) / n_board)
                logger.info('Cut detector %s split%d, left_events = %d, %.2f%%'
                            % (self._det.get_det_name, i_split, left_event, 100. * left_event / 65535))

            _data = select_hit_all_boards(_data)
            hit_all_events = int(len(_data) / n_board)
            end_time = _data['timeTag'].iloc[-1]
            logger.info('Detector %s split%d: hit_all_boards_events = %d, %.3f%%'
                        % (self._det.get_det_name, i_split, hit_all_events, 100. * hit_all_events / 65535))

            # get channel-replaced data matrix with shape (, n_bar) and calculate muon hit information
            _data_mtx = get_channel_replaced_data(_data, bar_channel_config=self.bar_channel_config,
                                                  threshold_mtx=threshold_mtx, baseline_mtx=baseline_mtx,
                                                  n_board=n_board, n_channel=n_channel, n_bar=n_bar)
            if self.pro_data_frame is None:
                self.pro_data_frame = _data_mtx
            else:
                np.concatenate((self.pro_data_frame, _data_mtx), axis=0)
            self._calc_hit_info_4plane_det(_data_mtx, i_split=i_split)

            if b_cut:
                n_split_used = i_split
                break
            n_split_used = i_split

        logger.info('Detector %s measured time: %s' % (
            self._det.get_det_name(), time.strftime('%D,%H:%M:%S', time.localtime(start_time))))
        logger.info('Detector %s measured time: %d' % (self._det.get_det_name(), int(end_time - start_time)))
        meas_time = TVectorD(1)
        meas_time[0] = end_time - start_time

        _tree_fname = self._det.get_tree_fname()
        str_hadd = 'hadd %s' % _tree_fname
        split_files = ''
        for i in range(n_split_used + 1):
            split_fname = _tree_fname.replace('.root', '_split%d.root' % i)
            split_files = split_files + ' %s' % split_fname
        str_hadd = str_hadd + split_files

        os.system(str_hadd)
        os.system('rm %s' % split_files)

        fout = TFile(_tree_fname, 'UPDATE')
        meas_time.Write('meas_time')
        fout.Save()
        fout.Close()

        # self._draw_occupancy()

    def occupancy_time_slice(self):
        self.find_baseline()
        self.find_threshold()
        # global bin_width
        logger.info('occupancy_time_slice for detector %s' % self._det.get_det_name())
        # fname_prefix = os.path.join(self._global_dir_tag[0], 'Occupancy_slice_%s' % self._det.get_det_name())
        self._can_occ = TCanvas('can_occupancy', 'can_occupancy', 2000, 1500)

        measData_h5_fname = self._det.get_measData_h5_fname()
        pdf_fname = measData_h5_fname.replace('.h5', '_occpancy.pdf')
        root_fname = measData_h5_fname.replace('.h5', '_occpancy.root')
        self._fout_occ = TFile(root_fname, 'recreate')

        h5_obj = baseCore.h5Data(measData_h5_fname, 'r')
        time_interval = self._det.get_measData_time_interval()
        n_split = h5_obj.get_n_splits()
        _data0 = h5_obj.getData(0)
        n_board = get_n_boards(_data0)
        n_channel = self._det.get_n_channels()
        n_bar = self._det.get_n_bars()
        self._get_channel_replace_information(self._det.get_channel_config_fname())
        threshold_mtx = self._threshold_csv.values
        baseline_mtx = self._baseline_csv.values

        # slice every 20 splits
        i = 0
        while i < n_split:
            _data = h5_obj.getData(i)
            i = i + 1
            for j in range(1, 20):
                t_data = h5_obj.getData(i)
                t_data['triggerID'] = t_data['triggerID'] + j * 65535
                _data = pd.concat([_data, t_data], ignore_index=True)
                i = i + 1
                if i >= n_split:
                    break

            _data = select_hit_all_boards(_data)

            i_slice = 0
            start_time = _data['timeTag'].iloc[0]
            end_time = _data['timeTag'].iloc[0] + time_interval
            while start_time < _data['timeTag'].iloc[-1]:
                tmp_data = _data[_data['timeTag'] >= start_time]
                tmp_data = tmp_data[tmp_data['timeTag'] < end_time]

                # get channel-replaced data matrix with shape (, n_bar) and calculate muon hit information
                _data_mtx = get_channel_replaced_data(tmp_data, bar_channel_config=self.bar_channel_config,
                                                      threshold_mtx=threshold_mtx, baseline_mtx=baseline_mtx,
                                                      n_board=n_board, n_channel=n_channel, n_bar=n_bar)
                self._draw_occupancy(_data_mtx, start_time=start_time, end_time=end_time, pdf_fname=pdf_fname)

                start_time = start_time + time_interval
                end_time = end_time + time_interval

        self._can_occ.Print('%s]' % pdf_fname)
        self._fout_occ.Save()
        self._fout_occ.Close()

    def _draw_occupancy(self, pro_data_mtx=None, start_time=0, end_time=0, pdf_fname=''):
        if pro_data_mtx is None:
            pro_data_mtx = self.pro_data_frame
        if pro_data_mtx is not None:
            n_board = self._threshold_csv.shape[0]
            n_bar = self._det.get_n_bars()
            self._can_occ.cd()
            self._fout_occ.cd()
            for i_board in range(n_board):
                h_name = h_title = 'h_occupancy_board%d' % i_board
                if not end_time == 0:
                    h_name = h_name + '_%d_to_%d' % (start_time, end_time)
                    h_title = h_title + '_%s-%s' % (time.strftime('%D,%H:%M:%S', time.localtime(start_time)),
                                                    time.strftime('%D,%H:%M:%S', time.localtime(end_time)))
                h1 = TH1F(h_name, h_title, n_bar, 0, n_bar)
                h1_single = TH1F(h_name + '_single', h_title + ': single', n_bar, 0, n_bar)
                h1_double = TH1F(h_name + '_double', h_title + ': double', n_bar, 0, n_bar)
                for i_bar in range(n_bar):
                    _larger_than_0 = pro_data_mtx[:, i_board, i_bar] > 0

                    _count_single = _count_double = 0
                    if i_bar == 0:
                        _larger_than_0_next = pro_data_mtx[:, i_board, i_bar + 1] > 0
                        _count_double = (_larger_than_0 & _larger_than_0_next).sum()
                        _count_single = (_larger_than_0 & ~_larger_than_0_next).sum()
                    elif i_bar == n_bar - 1:
                        _larger_than_0_prev = pro_data_mtx[:, i_board, i_bar - 1] > 0
                        _count_double = (_larger_than_0 & _larger_than_0_prev).sum()
                        _count_single = (_larger_than_0 & ~_larger_than_0_prev).sum()
                    else:
                        _larger_than_0_next = pro_data_mtx[:, i_board, i_bar + 1] > 0
                        _larger_than_0_prev = pro_data_mtx[:, i_board, i_bar - 1] > 0
                        _count_double = ((_larger_than_0 & _larger_than_0_prev) | (_larger_than_0 & _larger_than_0_next)).sum()
                        _count_single = ((_larger_than_0 & ~_larger_than_0_prev) & (_larger_than_0 & ~_larger_than_0_next)).sum()

                    _count = _larger_than_0.sum()
                    h1.SetBinContent(i_bar + 1, _count)
                    h1_single.SetBinContent(i_bar + 1, _count_single)
                    h1_double.SetBinContent(i_bar + 1, _count_double)
                h1.Draw('hist')
                gPad.SetGridx()
                self._can_occ.Print('%s(' % pdf_fname)
                h1.Write()
                h1_single.Draw('hist')
                gPad.SetGridx()
                self._can_occ.Print('%s(' % pdf_fname)
                h1_single.Write()
                h1_double.Draw('hist')
                gPad.SetGridx()
                self._can_occ.Print('%s(' % pdf_fname)
                h1_double.Write()

    def _calc_hit_info_4plane_det(self, _data_mtx, i_split):
        # get detector properties
        det_info_fname = self._det.get_det_pro_fname()
        det_info = pd.read_csv(det_info_fname)
        t_tree = TTree('t_tree', 't_tree')

        # define the tree by detector's type
        plane_det_vars = ['X1', 'Y1', 'Z1', 'X2', 'Y2', 'Z2',
                          'xo1', 'zx1', 'xo2', 'zx2', 'yo1', 'zy1', 'yo2', 'zy2',
                          'thetax_trk', 'thetay_trk', 'thetax', 'thetay', 'theta', 'phi',
                          'xNBar1', 'xNBar2', 'yNBar1', 'yNBar2']

        if self._det.get_det_type().__contains__('plane'):
            for i_var in plane_det_vars:
                define_var = "%s = array('d', [0])" % i_var
                tree_branch = "t_tree.Branch('%s', %s, '%s/D')" % (i_var, i_var, i_var)
                if i_var.__contains__('NBar'):
                    define_var = define_var.replace("'d'", "'i'")
                    tree_branch = tree_branch.replace('/D', '/I')
                exec(define_var)
                exec(tree_branch)

        axis0 = _data_mtx.shape[0]
        axis1 = _data_mtx.shape[1]
        if not axis1 == 4:
            logger.warning('Attention: Using _calc_one_event_4plane_det, but n_board = %d.' % axis1)
        logger.info('Total events: %d' % axis0)

        count_mode = []
        Z1_aver = (det_info['Origin_z'][0] + det_info['Origin_z'][1]) / 2
        Z2_aver = (det_info['Origin_z'][2] + det_info['Origin_z'][3]) / 2
        rotAnglex = self._det.get_rotAngle()[0]
        rotAngley = self._det.get_rotAngle()[1]

        for i_event in range(axis0):
            event_data = _data_mtx[i_event, :, :]
            event_mtx = calc_one_event_4plane_det(event_data, det_info)
            for i in range(4):
                count_mode.append(event_mtx[i, 2])
            if not -1 in event_mtx[:, -1]:
                lis_var = tree_variables_plane_det(event_mtx, Z1_aver, Z2_aver, rotAnglex, rotAngley)
                for i_var in range(len(plane_det_vars)):
                    str_pass_value = '%s[0] = lis_var[%d]' % (plane_det_vars[i_var], i_var)
                    exec(str_pass_value)

                t_tree.Fill()

        fout = TFile(self._det.get_tree_fname().replace('.root', '_split%d.root' % i_split), 'recreate')
        t_tree.Write()
        fout.Save()
        fout.Close()
        print('After calculate every event hit info, # board info: %d' % len(count_mode))
        for i_mode in list(set(count_mode)):
            i_mode_count = count_mode.count(i_mode)
            logger.info('Bar trigger mode %d: %d, %.3f' % (int(i_mode), i_mode_count, i_mode_count / len(count_mode)))


def tree_variables_plane_det(event_mtx, Z1_aver, Z2_aver, rotAnglex, rotAngley):
    xo1 = event_mtx[0, 0]
    zx1 = event_mtx[0, 1]
    xNBar1 = int(event_mtx[0, -1])
    xo2 = event_mtx[1, 0]
    zx2 = event_mtx[1, 1]
    xNBar2 = int(event_mtx[1, -1])
    yo1 = event_mtx[2, 0]
    zy1 = event_mtx[2, 1]
    yNBar1 = int(event_mtx[2, -1])
    yo2 = event_mtx[3, 0]
    zy2 = event_mtx[3, 1]
    yNBar2 = int(event_mtx[3, -1])
    Z1 = Z1_aver
    Z2 = Z2_aver
    X1 = (Z1 - zx1) * (xo1 - xo2) / (zx1 - zx2) + xo1
    X2 = (Z2 - zx1) * (xo1 - xo2) / (zx1 - zx2) + xo1
    Y1 = (Z1 - zy1) * (yo1 - yo2) / (zy1 - zy2) + yo1
    Y2 = (Z2 - zy1) * (yo1 - yo2) / (zy1 - zy2) + yo1
    v_det_coord1 = TVector3(X1, Y1, Z1)
    v_det_coord2 = TVector3(X2, Y2, Z2)

    thetax_trk = TMath.ATan2((xo1 - xo2), (zx1 - zx2))
    thetay_trk = TMath.ATan2((yo1 - yo2), (zy1 - zy2))

    v_det_angle = TVector3(TMath.Tan(thetax_trk), TMath.Tan(thetay_trk), 1)
    v_det_coord = v_det_coord1 - v_det_coord2
    v_car_coord = TVector3(v_det_coord)

    if rotAnglex > 1e-10 and rotAngley > 1e-10:
        logger.error('rotAnglex>1e-10 && rotAngley>1e-10.')
        raise ValueError('rotAnglex>1e-10 && rotAngley>1e-10.')

    if rotAnglex < 1e-10:
        v_car_coord.RotateX(rotAngley)

    if rotAngley < 1e-10:
        v_car_coord.RotateY(rotAnglex)

    thetax = thetay = -10
    if rotAnglex < 1e-10:
        thetax = TMath.ATan2(v_car_coord(0), v_car_coord(2))
        thetay = thetay_trk + rotAngley

    if rotAngley < 1e-10:
        thetay = TMath.ATan2(v_car_coord(1), v_car_coord(2))
        thetax = thetax_trk + rotAnglex

    theta = v_car_coord.Theta()
    phi = v_car_coord.Phi()

    return [X1, Y1, Z1, X2, Y2, Z2, xo1, zx1, xo2, zx2, yo1, zy1, yo2, zy2, thetax_trk, thetay_trk, thetax, thetay,
            theta, phi, xNBar1, xNBar2, yNBar1, yNBar2]


def hit_position_nbar1(triggered_bar, det_info):
    start_x = get_start_x(det_info)
    start_z = det_info['Origin_z']
    tri_up = get_tri_up(det_info)
    bar_interval = det_info['barInterval']
    bar_height = det_info['barHeight']

    if triggered_bar % 2 == 1:
        tri_up = not tri_up
    if tri_up:
        sign = 1.
    else:
        sign = -1.
    x = start_x + bar_interval / 2 * (triggered_bar + 1)
    z = start_z + sign * bar_height / 2
    return x, z


def hit_position_nbar2(triggered_bar, triggered_energy, det_info):
    start_x = get_start_x(det_info)
    start_z = det_info['Origin_z']
    tri_up = get_tri_up(det_info)
    bar_interval = det_info['barInterval']
    bar_height = det_info['barHeight']

    if triggered_bar[0] % 2 == 1:
        tri_up = not tri_up
    if tri_up:
        sign = 1.
    else:
        sign = -1.

    energy_ratio = 1. * triggered_energy[1] / (triggered_energy[0] + triggered_energy[1])
    x = start_x + bar_interval / 2 * (triggered_bar[0] + 1 + energy_ratio)
    z = start_z + sign * (0.5 - energy_ratio) * bar_height

    return x, z


def hit_position_center2(triggered_bar, triggered_energy, det_info):
    start_x = get_start_x(det_info)
    start_z = det_info['Origin_z']
    tri_up = get_tri_up(det_info)
    bar_interval = det_info['barInterval']
    bar_height = det_info['barHeight']

    if triggered_bar[0] % 2 == 1:
        tri_up = not tri_up
    if tri_up:
        sign = 1.
    else:
        sign = -1.

    energy_ratio = 1. * triggered_energy[0] / (triggered_energy[0] + triggered_energy[1])
    x = start_x + bar_interval / 2 * (triggered_bar[0] + 1 + energy_ratio)
    z = start_z + sign * (0.5 - energy_ratio) * bar_height

    return x, z


def hit_position_nbar3(triggered_bar, triggered_energy, det_info, itself, other_lis):
    if len(triggered_bar) == 4:
        if triggered_bar[-1] - triggered_bar[-2] > 1:
            triggered_bar = triggered_bar[:-1]
            triggered_energy = triggered_energy[:-1]
        if triggered_bar[1] - triggered_bar[0] > 1:
            triggered_bar = triggered_bar[1:]
            triggered_energy = triggered_energy[:-1]

    other_x = other_lis[0]
    other_z = other_lis[1]
    it_x = itself[0]
    it_z = itself[1]

    if (it_z - other_z) * (it_x - other_x) > 0:  # go left
        go_where = 'go left'
        up_down = [False, True]
    else:  # go right
        go_where = 'go right'
        up_down = [True, False]

    tri_up = get_tri_up(det_info)

    trigger_up = triggered_bar % 2 == 0
    trigger_up = trigger_up == tri_up

    if trigger_up[:-1].tolist() == up_down:
        triggered_bar = triggered_bar[:-1]
        triggered_energy = triggered_energy[:-1]
    elif trigger_up[1:].tolist() == up_down:
        triggered_bar = triggered_bar[1:]
        triggered_energy = triggered_energy[1:]
    else:
        logger.error('hit_position_nbar3:\nray direction: %s; triggered_bar: %s' %
                     (go_where, ' '.join([str(i) for i in triggered_bar])))

    x, z = hit_position_center2(triggered_bar, triggered_energy, det_info)

    return x, z


def hit_position_nbar4(triggered_bar, triggered_energy, det_info):
    triggered_bar = triggered_bar[1:2]
    triggered_energy_new = np.full([2], 0)
    triggered_energy_new[0] = triggered_energy[0] + triggered_energy[1]
    triggered_energy_new[1] = triggered_energy[2] + triggered_energy[3]

    x, z = hit_position_center2(triggered_bar, triggered_energy_new, det_info)

    return x, z


def calc_two_board_plane_det(_data, i_board1, i_board2, det_info):
    board1_data = _data[i_board1, :]
    triggered_energy1 = board1_data[board1_data > 0]
    triggered_bar1 = np.where(board1_data > 0)[0]
    board2_data = _data[i_board2, :]
    triggered_energy2 = board2_data[board2_data > 0]
    triggered_bar2 = np.where(board2_data > 0)[0]
    lis1 = calc_one_board_plane_det(triggered_bar=triggered_bar1, triggered_energy=triggered_energy1,
                                    det_info=det_info.iloc[i_board1])
    lis2 = calc_one_board_plane_det(triggered_bar=triggered_bar2, triggered_energy=triggered_energy2,
                                    det_info=det_info.iloc[i_board2])

    invalid_lis = [0, 0, -1, -1]
    if lis1[-1] > 2 or lis2[-1] > 2:
        del_x = abs(lis1[0] - lis2[0])
        board1_z = det_info['Origin_z'][i_board1]
        board2_z = det_info['Origin_z'][i_board2]
        del_z = abs(board1_z - board2_z)
        ray_agl = TMath.ATan2(del_x, del_z)
        agl1 = TMath.ATan2(0.5 * det_info['barInterval'][i_board1], det_info['barHeight'][i_board1])
        agl2 = TMath.ATan2(0.5 * det_info['barInterval'][i_board2], det_info['barHeight'][i_board2])
        bar_agl = min(agl1, agl2)
        if ray_agl < bar_agl:
            if lis1[-2] == 4 and not lis2[-2] == -1:
                lis1 = calc_one_board_plane_det(triggered_bar=triggered_bar1, triggered_energy=triggered_energy1,
                                                det_info=det_info.iloc[i_board1], mode4_NBar=2)
            if lis2[-2] == 4 and not lis1[-2] == -1:
                lis2 = calc_one_board_plane_det(triggered_bar=triggered_bar2, triggered_energy=triggered_energy2,
                                                det_info=det_info.iloc[i_board2], mode4_NBar=2)
        else:
            lis1 = invalid_lis
            lis2 = invalid_lis

    if lis1[-1] == 3 and not lis2[-1] == -1:
        x, z = hit_position_nbar3(triggered_bar1, triggered_energy1, det_info.iloc[i_board1],
                                  itself=lis1, other_lis=lis2)
        lis1 = [x, z, lis1[-2], lis1[-1]]
    if lis2[-1] == 3 and not lis1[-1] == -1:
        x, z = hit_position_nbar3(triggered_bar2, triggered_energy2, det_info.iloc[i_board2],
                                  itself=lis2, other_lis=lis1)
        lis2 = [x, z, lis2[-2], lis2[-1]]
    return lis1, lis2


def calc_one_board_plane_det(triggered_bar, triggered_energy, det_info, mode4_NBar=0):
    # trigger mode for plane detector
    # one triggered: mode 1 ; []; NBar 1
    # two neighbour triggerd: mode 2; [1]; NBar 2
    # two interval triggered: mode 7; [2]; discard
    # two far triggered: mode -1; [a] a>2; discard
    # three neighbour triggered: mode 3; [1,1]; NBar 3
    # three one neighbour one interval: mode 4; [1,2] or [2,1]; NBar 4
    # three two two interval: mode -1; [2,2]; discard
    # three far triggered: mode 4; [a,b] a>2 or b>2; discard
    # four neighbour triggered: mode 6; [1,1,1]; NBar 4
    # four triggered, three neighbour: mode 5; [1,1,2] or [2,1,1]; NBar 3
    # four triggered, two two neighbour: mode -1; [1,2,1] discard
    # four triggered, far: mode -1; [a,b,c] >2 .count()>1; discard

    if len(triggered_bar) == 0:
        return [0, 0, 8, -1]

    if len(triggered_bar) == 1:
        mode = 1
        NBar = 1
        x, z = hit_position_nbar1(triggered_bar[0], det_info)
        return [x, z, mode, NBar]

    if triggered_bar[-1] - triggered_bar[0] > 4:
        return [0, 0, -1, -1]

    minus1 = triggered_bar[1:]
    minus2 = triggered_bar[:-1]
    minus = (minus1 - minus2).tolist()

    if minus == [1]:
        mode = 2
        NBar = 2
        x, z = hit_position_nbar2(triggered_bar, triggered_energy, det_info)
        return [x, z, mode, NBar]

    if minus == [2]:
        mode = 7
        return [0, 0, mode, -1]

    if minus == [1, 1]:
        mode = 3
        NBar = 3
        center_bar = triggered_bar[1]
        x, z = hit_position_nbar1(center_bar, det_info)
        return [x, z, mode, NBar]

    if minus == [1, 2] or minus == [2, 1]:
        if mode4_NBar == 2:
            mode = 4
            NBar = 2
            if minus[0] == 2:
                triggered_bar_new = triggered_bar[1:]
                triggered_energy_new = triggered_energy[1:]
            else:
                triggered_bar_new = triggered_bar[:-1]
                triggered_energy_new = triggered_energy[:-1]
            x, z = hit_position_nbar2(triggered_bar_new, triggered_energy_new, det_info)
            return [x, z, mode, NBar]

        mode = 4
        NBar = 4
        triggered_bar_new = np.array(list(range(triggered_bar[0], triggered_bar[-1] + 1)))
        triggered_energy_new = np.full([4], 0)
        ii = 0
        for i in range(4):
            if triggered_bar_new[i] in triggered_bar:
                triggered_energy_new[i] = triggered_energy[ii]
                ii = ii + 1
        x, z = hit_position_nbar4(triggered_bar_new, triggered_energy_new, det_info)
        return [x, z, mode, NBar]

    if len(minus) == 2 and minus[0] == 1 and minus[1] > 2:
        mode = 4
        NBar = 2
        triggered_bar_new = triggered_bar[:-1]
        triggered_energy_new = triggered_energy[:-1]
        x, z = hit_position_nbar2(triggered_bar_new, triggered_energy_new, det_info)
        return [x, z, mode, NBar]

    if len(minus) == 2 and minus[1] == 1 and minus[0] > 2:
        mode = 4
        NBar = 2
        triggered_bar_new = triggered_bar[1:]
        triggered_energy_new = triggered_energy[1:]
        x, z = hit_position_nbar2(triggered_bar_new, triggered_energy_new, det_info)
        return [x, z, mode, NBar]

    if minus == [1, 1, 1]:
        mode = 6
        NBar = 4
        x, z = hit_position_nbar4(triggered_bar, triggered_energy, det_info)
        return [x, z, mode, NBar]

    if len(minus) == 3 and minus[:-1] == [1, 1] and minus[-1] > 1:
        mode = 5
        NBar = 3
        triggered_bar_new = triggered_bar[:-1]
        triggered_energy_new = triggered_energy[:-1]
        lis = calc_one_board_plane_det(triggered_bar_new, triggered_energy_new, det_info, mode4_NBar)
        return [lis[0], lis[1], mode, NBar]

    if len(minus) == 3 and minus[1:] == [1, 1] and minus[0] > 1:
        mode = 5
        NBar = 3
        triggered_bar_new = triggered_bar[1:]
        triggered_energy_new = triggered_energy[1:]
        lis = calc_one_board_plane_det(triggered_bar_new, triggered_energy_new, det_info, mode4_NBar)
        return [lis[0], lis[1], mode, NBar]

    return [0, 0, -1, -1]


def get_start_x(det_info):
    start_x = 0
    if str(det_info['xy_board']).__contains__('y') or int(det_info['xy_board']) == 0:  # y_board
        start_x = det_info['Origin_y']
    elif str(det_info['xy_board']).__contains__('x') or int(det_info['xy_board']) == 1:
        start_x = det_info['Origin_x']
    else:
        logger.error('In DetInfo, xy_board is not x nor y.')
        raise ValueError('In DetInfo, xy_board is not x nor y.')
    return start_x


def get_tri_up(det_info):
    tri_up = True
    if str(det_info['TriBarUp']).lower().__contains__('up') or int(det_info['TriBarUp']) == 1:  # y_board
        tri_up = True
    elif str(det_info['TriBarUp']).lower().__contains__('down') or int(det_info['TriBarUp']) == 0:
        tri_up = False
    else:
        logger.error('In DetInfo, TriBarUp is not up nor down.')
        raise ValueError('In DetInfo, TriBarUp is not up nor down.')
    return tri_up

# spc = TSpectrum()
# spc.Search(tmp_h1f, 1)
# peak_bin = spc.GetPositionX()
# logger.info('Get %d peaks from board %d channel %d' % (len(peak_bin), i_board, i_channel))
# peak_pos = tmp_h1f.GetBinCenter(1 + int(peak_bin[0] + 0.5))
# peak_height = tmp_h1f.GetBinContent(1 + int(peak_bin[0] + 0.5))
