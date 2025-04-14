import sys
from array import array
import os, logging

import datetime

from DatabaseTool import *
from DataFactory import CommonUtils, ROOTUtils

# from MyTools.Monitor.Process_monitor import process_monitoring
# process_monitoring(interval=5)

from ROOT import TH1F, TFile, TMath, TTree, RDF, TVectorD, TVector3, TGraph, TMultiGraph, TLegend, TGraphErrors
import ROOT

logger = logging.getLogger(__name__)


class DataframeAnalyzer:
    def __init__(self):
        pass

    @staticmethod
    def FindBarChannelRelation(h5_fname, csv_out_fname, directory, n_bars=0):
        # (root_dir, h5_name) = os.path.split(h5_fname)
        h5_name = CommonUtils.splitext_h5_fname(h5_fname)

        df_data = None

        for i in Cormis_DataReader_hdf.getCaching(h5_name, 100000, directory):
            df_data = i.copy()
            break

        n_boards = CommonUtils.count_n_boards(df_data)
        n_channels = CommonUtils.count_n_channels(df_data)

        chn_col_names = df_data.columns[df_data.columns.str.contains(r'chn_\d+')]

        if n_bars != 0:
            n_bars = CommonUtils.count_n_bars(df_data)
            # check data trigger count
            if (df_data[chn_col_names] > 0).sum(axis=1).max() > n_bars:
                raise Exception('%d bars are triggered in %s' % (n_bars, h5_fname))
        else:
            n_bars = CommonUtils.count_n_bars(df_data)

        lis_bar = []  # columns name for building the relation dataframe
        lis_relation = []  # bar channel relation for building the relation dataframe
        for i_board in range(n_boards):
            certain_board_data = df_data[df_data['boardID'] == i_board]
            chn_data = certain_board_data[chn_col_names]

            sum_chn_data = (chn_data > 0).sum()
            _n_bars = sum_chn_data[sum_chn_data > 5].shape[0]
            if _n_bars != n_bars:
                logger.warning(
                    'data of all boards indicates that one plane has %d bars,'
                    ' but board %d indicates %d bars.' % (n_bars, i_board, _n_bars))
                print('data of all boards indicates that one plane has %d bars,'
                      ' but board %d indicates %d bars.' % (n_bars, i_board, _n_bars))
            # fill lis_bar
            if i_board == 0:
                for i_bar in range(n_bars):
                    lis_bar.append('bar_%d' % i_bar)

            lis_channel_number = []

            # fill lis_channel_number
            replace_number = 0
            for i_bar in range(n_bars):
                chn_name = 'chn_%d' % i_bar
                if sum_chn_data[chn_name] > 0:  # this channel connect with a bar
                    lis_channel_number.append(i_bar)
                else:
                    while 1:
                        if (n_bars + replace_number) > n_channels:
                            CommonUtils.pop_up_messages(
                                'Cannot find a replaced channel for board%d bar %d' % (i_board, i_bar))
                            lis_channel_number.append(i_bar)
                        if sum_chn_data['chn_%d' % (n_bars + replace_number)] > 0:
                            lis_channel_number.append(n_bars + replace_number)
                            replace_number += 1
                            break
                        replace_number += 1

            lis_relation.append(lis_channel_number)

        # build the dataframe and save it into csv file
        df_relation = pd.DataFrame(lis_relation, columns=lis_bar)
        df_relation.to_csv(csv_out_fname, index=False)

    @staticmethod
    def FindBaseline(h5_fname, csv_out_fname, directory, root_fname=''):
        if root_fname == '':
            root_fname = csv_out_fname.replace('.csv', '.root')

        # noinspection DuplicatedCode
        h5_name = CommonUtils.splitext_h5_fname(h5_fname)

        df_data = None

        for i in Cormis_DataReader_hdf.getCaching(h5_name, 1000000, directory):
            df_data = i.copy()
            break

        n_boards = CommonUtils.count_n_boards(df_data)
        n_channels = CommonUtils.count_n_channels(df_data)

        fout = TFile(root_fname, 'recreate')

        # find baseline
        lis_channel = []
        lis_baseline = []
        for i_board in range(n_boards):
            certain_board_data = df_data[df_data['boardID'] == i_board]

            lis_baseline_board = []
            for i_channel in range(n_channels):
                chn_name = 'chn_%d' % i_channel

                if i_board == 0:
                    lis_channel.append(chn_name)

                if certain_board_data[chn_name].sum() < 1:
                    lis_baseline_board.append(0)
                    continue

                h_name = 'h_board%d_chn%d' % (i_board, i_channel)
                h_title = h_name + '; Energy; Counts'

                # save the spectrum
                val = certain_board_data[chn_name].values
                h1_chn_data_write = ROOTUtils.fill_array_into_h1f(val, h_name, h_title)
                h1_chn_data_write.Write()

                # find the baseline
                larger_than_zero = val[val > 0]
                h1_chn_data = ROOTUtils.fill_array_into_h1f(larger_than_zero, h_name + '_larger_than_zero', h_title)

                peak_pos = int(h1_chn_data.GetBinCenter(h1_chn_data.GetMaximumBin()))
                lis_baseline_board.append(peak_pos)

            lis_baseline.append(lis_baseline_board)

        fout.Save()
        fout.Close()

        # build the dataframe and save it into csv file
        df_baseline = pd.DataFrame(lis_baseline, columns=lis_channel)
        df_baseline.to_csv(csv_out_fname, index=False)

    @staticmethod
    def DrawBaselineList(lis_csv_fname, out_fname):
        lis_df_baseline = []
        lis_csv_name = []
        n_channel = pd.read_csv(lis_csv_fname[0]).shape[1]
        n_board = pd.read_csv(lis_csv_fname[0]).shape[0]

        for i_csv_fname in lis_csv_fname:
            df_baseline = pd.read_csv(i_csv_fname)
            csv_name = os.path.splitext(os.path.split(i_csv_fname)[1])[0]

            lis_df_baseline.append(df_baseline)
            lis_csv_name.append(csv_name)

        root_fname = os.path.splitext(out_fname)[0] + '.root'
        fout = TFile(root_fname, 'recreate')

        for i_board in range(n_board):
            mg = TMultiGraph('mg_board%d' % i_board, 'baseline compare board%d; channel; baseline' % i_board)
            lg = TLegend(0.6, 0.74, 0.90, 0.90)
            lg.SetFillStyle(0)
            lg.SetName('lg_mg_board%d' % i_board)
            for i in range(len(lis_df_baseline)):
                df_baseline = lis_df_baseline[i]
                csv_name = lis_csv_name[i]

                xx = array('d', range(n_channel))
                yy = array('d', df_baseline.iloc[i_board].values)

                gr = TGraph(n_channel, xx, yy)
                gr.SetName('%s_board%d' % (csv_name, i_board))
                gr.SetLineColor(i + 1)
                gr.SetMarkerColor(i + 1)

                mg.Add(gr, "lp")
                lg.AddEntry(gr, gr.GetName(), 'l')

            mg.Write()
            lg.Write()

        fout.Save()
        fout.Close()

    @staticmethod
    def DrawEnergySpectrum(h5_fname, directory, root_fname, key_name: str = 'data', dic_split=None):
        # noinspection DuplicatedCode
        h5_name = CommonUtils.splitext_h5_fname(h5_fname)
        n_rows = CommonUtils.get_dataframe_total_lines(os.path.join(directory, '%s.h5' % h5_name), key_name=key_name)

        t_df = None
        if key_name == 'data':
            for i in Cormis_DataReader_hdf.getCaching(h5_name, 1000, directory):
                t_df = i.copy()
                break
        else:
            try:
                t_df = pd.read_hdf(os.path.join(directory, '%s.h5' % h5_name), key=key_name, start=0, stop=1000)
                if t_df is None:
                    raise Exception('Cannot find %s in %s.h. Please correct the key name.' % (key_name, h5_name))
            except Exception as e:
                print(e)

        n_boards = CommonUtils.count_n_boards(t_df)
        n_channels = CommonUtils.count_n_channels(t_df)
        chn_col_names = t_df.columns[t_df.columns.str.contains(r'chn_\d+')]

        lis_hists = n_boards * n_channels * [None]

        i_page = 0
        count = 0
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 100000, directory):
            count += i.shape[0]

            for i_board in range(n_boards):
                certain_board_data = i[i['boardID'] == i_board]

                for i_channel in range(n_channels):
                    chn_name = 'chn_%d' % i_channel

                    np_chn_data = certain_board_data[chn_name].values
                    np_chn_data = np_chn_data[np_chn_data > 0]

                    if lis_hists[n_channels * i_board + i_channel] is None:
                        if len(np_chn_data) > 1:
                            h_name = 'h%s_spectrum_board%d_chn%d' % (key_name, i_board, i_channel)
                            h_title = h_name + '; ADC; Counts'

                            if dic_split is not None:
                                h_name += '_%d' % dic_split['i_time']
                                h_title = h_name + ' %s - %s' % (
                                    dic_split['start_time'], dic_split['end_time']) + '; ADC; Counts'

                            h1 = ROOTUtils.fill_array_into_h1f(np_chn_data, h_name, h_title)
                            lis_hists[n_channels * i_board + i_channel] = h1
                        else:
                            lis_hists[n_channels * i_board + i_channel] = -1

                    else:
                        h1 = lis_hists[i_channel + i_board * n_channels]
                        if h1 == -1 or len(np_chn_data) < 1:
                            continue
                        h1 = ROOTUtils.add_array_to_h1f(np_chn_data, h1)
                        lis_hists[i_channel + i_board * n_channels] = h1

            print('\rDrawEnergySpectrum for %s in %s: %d / %d    ' % (h5_name, key_name,
                                                                      count, n_rows), end='')
            i_page += 1

        print('')

        fout = TFile(root_fname, 'recreate')

        for h in lis_hists:
            if h == -1:
                continue
            h.Write()

        fout.Save()
        fout.Close()

    @staticmethod
    def FindThreshold(h5_fname, directory, baseline_csv, csv_out_fname, root_fname='', smooth_level=2, offset=0):
        if root_fname == '':
            root_fname = csv_out_fname.replace('.csv', '.root')

        # noinspection DuplicatedCode
        h5_name = CommonUtils.splitext_h5_fname(h5_fname)
        n_rows = CommonUtils.get_dataframe_total_lines(os.path.join(directory, '%s.h5' % h5_name))

        print('\rFindThreshold for %s' % h5_name, end='')

        df_data = None
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 1000, directory):
            df_data = i.copy()
            break

        n_boards = CommonUtils.count_n_boards(df_data)
        n_channels = CommonUtils.count_n_channels(df_data)

        # find threshold
        lis_channel = []
        lis_threshold = []
        lis_hists = []

        pageIndex = 0
        count = 0
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 100000, directory):
            count += i.shape[0]
            for i_board in range(n_boards):
                certain_board_data = i[i['boardID'] == i_board]

                for i_channel in range(n_channels):
                    chn_name = 'chn_%d' % i_channel

                    np_chn_data = certain_board_data[chn_name].values
                    np_chn_data = np_chn_data[np_chn_data > 0]

                    if pageIndex == 0:
                        if i_board == 0:
                            lis_channel.append(chn_name)

                        if len(np_chn_data) > 1:
                            h_name = 'h_energy_spectrum_board%d_chn%d' % (i_board, i_channel)
                            h_title = h_name + '; ADC; counts'

                            h1 = ROOTUtils.fill_array_into_h1f(np_chn_data, h_name, h_title)
                            lis_hists.append(h1)
                        else:
                            lis_hists.append(-1)

                    else:
                        h1 = lis_hists[i_channel + i_board * n_channels]
                        if h1 == -1 or len(np_chn_data) < 1:
                            continue
                        h1 = ROOTUtils.add_array_to_h1f(np_chn_data, h1)
                        lis_hists[i_channel + i_board * n_channels] = h1

            pageIndex += 1
            print('\rFindThreshold for %s: %d / %d    ' % (h5_name, count, n_rows), end='')
        print('')

        # read baseline csv as dataframe
        df_baseline = pd.read_csv(baseline_csv)

        # find threshold and fill into dataframe
        fout = TFile(root_fname, 'recreate')
        for i_board in range(n_boards):
            lis_threshold_board = []
            for i_channel in range(n_channels):
                i_hist = lis_hists[i_channel + i_board * n_channels]

                if i_hist == -1:
                    lis_threshold_board.append(0)
                    continue

                i_hist.Write()

                h_smooth = i_hist.Clone()
                if smooth_level > 0:
                    h_smooth.Smooth(smooth_level)

                t_baseline = df_baseline.iloc[i_board, i_channel]
                baseline_bin = h_smooth.FindBin(t_baseline)

                h_smooth.GetXaxis().SetRange(baseline_bin + offset, baseline_bin + 20 + int(offset / 3))
                threshold_pos = int(h_smooth.GetBinCenter(h_smooth.GetMinimumBin()))

                lis_threshold_board.append(threshold_pos)

            lis_threshold.append(lis_threshold_board)

        fout.Save()
        fout.Close()

        # build the dataframe and save it into csv file
        df_threshold = pd.DataFrame(lis_threshold, columns=lis_channel)
        df_threshold.to_csv(csv_out_fname, index=False)

    def TimeSlice(self, h5_fname, directory, root_fname, time_interval: int, pageCapacity=10000, checkFactor=3,
                  if_temperature=False):
        h5_name = CommonUtils.splitext_h5_fname(h5_fname)
        log_name = root_fname.replace('.root', '.ATTENTION')
        if os.path.exists(log_name):
            os.remove(log_name)
        dic_count_rate = self._get_time_slice_array(h5_fname, directory, time_interval, pageCapacity)

        df_data = None
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 10000, directory):
            df_data = i.copy()
            break

        n_boards = CommonUtils.count_n_boards(df_data)
        n_channels = CommonUtils.count_n_channels(df_data)
        del df_data

        dic_count_rate = self._add_dict_info(dic_count_rate, n_boards, n_channels,
                                             checkA=False, checkB=True, checkC=False)

        fout = TFile(root_fname, 'recreate')

        # time slice for trigger bar counts
        lis_start_time = dic_count_rate.get('start_time')
        lis_end_time = dic_count_rate.get('end_time')
        lis_temp_mean = dic_count_rate.get('temperature_mean')
        lis_temp_std = dic_count_rate.get('temperature_std')

        for i in range(len(lis_start_time)):
            start_time = lis_start_time[i]
            end_time = lis_end_time[i]

            for i_board in range(n_boards):
                lis_trigger_counts = dic_count_rate.get('board%d' % i_board)[i]

                if lis_trigger_counts.shape[0] < 1:
                    continue

                h_name = 'h_trigger_bar_count_board%d_slice%d' % (i_board, i)
                h_title = 'trigger_bar_count for board%d %s - %s; number of triggered channels; counts' % (
                    i_board, start_time, end_time)
                h1 = ROOTUtils.fill_array_into_h1f(lis_trigger_counts, h_name, h_title)

                h1.Write()

        # time slice for each channel
        lis_start_time = dic_count_rate.get('start_time')
        lis_end_time = dic_count_rate.get('end_time')
        for i_board in range(n_boards):
            if len(lis_start_time) < 1:
                break
            for i_channel in range(n_channels):
                chn_name = 'board%s_channel%d' % (i_board, i_channel)
                lis_count_rate = dic_count_rate.get(chn_name)
                lis_count_err = dic_count_rate.get(chn_name + '_err')
                h_name = 'h_count_rate_' + chn_name
                h_title = 'count_rate for ' + chn_name + '; ; counts'
                h1 = TH1F(h_name, h_title, len(lis_start_time), 0, len(lis_start_time))

                n_slice = len(lis_count_rate)
                arr_count = array('d', n_slice * [0.])
                arr_count_err = array('d', n_slice * [0.])
                arr_temp = array('d', n_slice * [0.])
                arr_temp_err = array('d', n_slice * [0.])

                for i in range(n_slice):
                    start_time = lis_start_time[i]
                    end_time = lis_end_time[i]

                    count_rate = lis_count_rate[i]
                    err_count_rate = lis_count_err[i]

                    h1.SetBinContent(i + 1, count_rate)
                    h1.SetBinError(i + 1, err_count_rate)
                    h1.GetXaxis().SetBinLabel(i + 1, '%s - %s'
                                              % (start_time.strftime('%m.%d %H:%M'), end_time.strftime('%m.%d %H:%M')))

                    arr_count[i] = count_rate
                    arr_count_err[i] = err_count_rate
                    arr_temp[i] = lis_temp_mean[i]
                    arr_temp_err[i] = lis_temp_std[i]

                h1.GetXaxis().SetLabelSize(0.03)
                h1.Write()

                if if_temperature:
                    gr_err = TGraphErrors(n_slice, arr_temp, arr_count, arr_temp_err, arr_count_err)
                    gr_err.SetName('gr_count_temperature_' + chn_name)
                    gr_err.SetTitle('gr_count_temperature_' + chn_name + '; temperature; counts')
                    gr_err.Write()

                # check count rate
                mean_val = np.array(lis_count_rate).mean()
                max_val = h1.GetMaximum()
                min_val = h1.GetMinimum()
                max_err = max(h1.GetBinError(h1.GetMaximumBin()), h1.GetBinError(h1.GetMaximumBin()))
                if ((max_val - mean_val) > checkFactor * max_err) \
                        or ((mean_val - min_val) > checkFactor * max_err):
                    logger.warning('ATTENTION!!! The count rate jumps: %s in %s'
                                   % (chn_name, h5_name))
                    print('ATTENTION!!! The count rate jumps: %s in %s'
                          % (chn_name, h5_name), file=open(log_name, 'a'))

        # temperature
        lis_temp_mean = dic_count_rate.get('temperature_mean')
        lis_temp_std = dic_count_rate.get('temperature_std')

        n_pnts = len(lis_temp_mean)
        h_temp = TH1F('h_temp', 'h_temp;;temperature', n_pnts, 0, n_pnts)

        for i in range(n_pnts):
            start_time = lis_start_time[i]
            end_time = lis_end_time[i]
            h_temp.SetBinContent(i + 1, lis_temp_mean[i])
            h_temp.SetBinError(i + 1, lis_temp_std[i])

            h_temp.GetXaxis().SetBinLabel(i + 1, '%s - %s'
                                          % (start_time.strftime('%m.%d %H:%M'), end_time.strftime('%m.%d %H:%M')))

        h_temp.GetXaxis().SetLabelSize(0.03)
        # gr = TGraphErrors(h_temp)
        # gr.SetName('gr_temperature')
        h_temp.Write()
        # gr.Write()

        # occupancy
        lis_start_time = dic_count_rate.get('start_time')
        lis_end_time = dic_count_rate.get('end_time')
        for i in range(len(lis_start_time)):
            start_time = lis_start_time[i]
            end_time = lis_end_time[i]

            for i_board in range(n_boards):
                h_name = 'h_occupancy_board%d_slice%d' % (i_board, i)
                h_title = 'occupancy for board%d %s - %s; channel; counts' % (i_board, start_time, end_time)
                h1 = TH1F(h_name, h_title, n_channels, 0, n_channels)

                for i_channel in range(n_channels):
                    chn_name = 'board%s_channel%d' % (i_board, i_channel)
                    lis_count_rate = dic_count_rate.get(chn_name)
                    h1.SetBinContent(i_channel + 1, lis_count_rate[i])
                    h1.SetBinError(i_channel + 1, TMath.Sqrt(lis_count_rate[i]))

                h1.Write()

        fout.Save()
        fout.Close()

    @classmethod
    def _calc_effective_count_rate(cls, count_rate, time_interval: int, time_delta: datetime.timedelta):
        if count_rate == 0:
            return 0, 0
        _time_delta = time_delta.total_seconds() + 1
        _count_rate = count_rate * time_interval / _time_delta
        _err_count = _count_rate / TMath.Sqrt(count_rate)
        return _count_rate, _err_count

    @classmethod
    def get_df_meas_time(cls, h5_fname):
        h5_name = CommonUtils.splitext_h5_fname(h5_fname)
        try:
            df_meas_time = pd.read_hdf(h5_fname, key="meas_time")
        except:
            _cr_offline = Cormis_DataReader_offline()
            start_time = _cr_offline.getStartTime(h5_fname)
            end_time = _cr_offline.getStopTime(h5_fname)
            meas_time = end_time - start_time
            df_meas_time = pd.DataFrame(columns=['time_delta', 'start_time'])
            df_meas_time.loc[h5_name] = [meas_time, start_time]
            print('Create key "meas_time" with start_time and end_time.')
        return df_meas_time

    def _get_time_slice_array(self, h5_fname, directory, time_interval: int, pageCapacity):
        # noinspection DuplicatedCode
        h5_name = CommonUtils.splitext_h5_fname(h5_fname)
        h5_fname = os.path.join(directory, '%s.h5' % h5_name)

        n_rows = CommonUtils.get_dataframe_total_lines(os.path.join(directory, '%s.h5' % h5_name))

        df_meas_time = self.get_df_meas_time(h5_fname)

        # fill the time slice list
        lis_start_time, lis_end_tme = CommonUtils.get_time_list(df_meas_time, time_interval)

        # get some properties of the dataframe
        df_data = None
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 10000, directory):
            df_data = i.copy()
            break
        n_boards = CommonUtils.count_n_boards(df_data)
        n_channels = CommonUtils.count_n_channels(df_data)
        n_bars = CommonUtils.count_n_bars(df_data)
        chn_col_names = df_data.columns[df_data.columns.str.contains(r'chn_\d+')]
        del df_data, i

        # initialize a dict for the array of time slice
        dic_count_rate = {'start_time': [], 'end_time': [],
                          'temperature_mean': [], 'temperature_std': []}
        for i_board in range(n_boards):
            dic_count_rate['board%d' % i_board] = []
            dic_count_rate['board%d_stats' % i_board] = []
            for i_channel in range(n_channels):
                dic_count_rate['board%s_channel%d' % (i_board, i_channel)] = []
                dic_count_rate['board%s_channel%d_err' % (i_board, i_channel)] = []

        # begin to slice
        df_remains = None
        df_data = None
        count = 0
        getCaching_iter = Cormis_DataReader_hdf.getCaching(h5_name, pageCapacity, directory).__iter__()

        for i_slice in range(len(lis_start_time)):
            _st = lis_start_time[i_slice]
            _et = lis_end_tme[i_slice]

            if df_remains is None or df_remains.shape[0] == 0:
                try:
                    i = getCaching_iter.__next__()
                    count += i.shape[0]
                    print('\rTimeSlice %d / %d = %.2f    ' % (count, n_rows, 1. * count / n_rows), end='')
                except:
                    print('%s: cannot get __next__(); i_slice = %d' % (h5_name, i_slice))
                    break
            else:
                i = df_remains

            while i['timeTag'].iloc[-1] <= (_et - datetime.timedelta(minutes=5)):
                try:
                    df_data = getCaching_iter.__next__()
                    count += df_data.shape[0]
                    print('\rTimeSlice %d / %d = %.2f    ' % (count, n_rows, 1. * count / n_rows), end='')
                except:
                    print('%s: cannot get __next__(); i_slice = %d' % (h5_name, i_slice))
                    break

                i = pd.concat([i, df_data], axis=0)

                # check memory
                if (sys.getsizeof(i) / 1024 / 1024 / 1024) > 0.8:
                    CommonUtils.pop_up_messages('time interval too large!')

            df_slice = i[(i['timeTag'] >= _st) & (i['timeTag'] < _et)]
            df_remains = i[i['timeTag'] >= _et]

            dic_count_rate = self._fill_time_slice_in_dict(df_slice, dic_count_rate, n_boards,
                                                           n_channels, chn_col_names, _st, _et, time_interval)
            print('\rslice %d    ' % i_slice, end='')

            del df_slice, i

        print('')

        return dic_count_rate

    def _fill_time_slice_in_dict(self, df_slice, dic_count_rate, n_boards, n_channels,
                                 chn_col_names, begin_time, end_time, time_interval):
        dic_count_rate['start_time'].append(begin_time)
        try:
            dic_count_rate['end_time'].append(df_slice['timeTag'].iloc[-1])
            end_time = df_slice['timeTag'].iloc[-1]
        except:
            dic_count_rate['end_time'].append(end_time)

        for i_board in range(n_boards):
            chn_data = df_slice[df_slice['boardID'] == i_board]
            chn_data = chn_data[chn_col_names]
            mtx_data, dic_data = CommonUtils.cvt_dataframe_to_numpy(chn_data)

            mtx_signal = mtx_data > 0
            # for trigger bar count
            mxt_trigger_count = mtx_signal.sum(axis=1)
            dic_count_rate['board%d' % i_board].append(mxt_trigger_count)
            dic_stats = {}
            for num in mxt_trigger_count:
                dic_stats[num] = dic_stats.get(num, 0) + 1
            dic_count_rate['board%d_stats' % i_board].append(dic_stats)

            # for channel slice and occupancy
            mtx_signal = mtx_signal.sum(axis=0)

            for i_channel in range(n_channels):
                ncount = mtx_signal[i_channel]
                ncount, err = self._calc_effective_count_rate(ncount, time_interval, time_delta=(end_time - begin_time))
                dic_count_rate['board%s_channel%d' % (i_board, i_channel)].append(ncount)
                dic_count_rate['board%s_channel%d_err' % (i_board, i_channel)].append(err)

        df_temp = df_slice['temperature']
        temp_mean = df_temp.mean()
        temp_std = df_temp.std()
        dic_count_rate['temperature_mean'].append(temp_mean)
        dic_count_rate['temperature_std'].append(temp_std)

        return dic_count_rate

    def _add_dict_info(self, dic, n_board, n_channel, checkA=False, checkB=False, checkC=False):
        dic['n_board'] = n_board
        dic['n_channel'] = n_channel

        if checkA:
            pass

        if checkB:
            dic['checkB'] = self._checkB_func(dic, n_board, n_channel)

        if checkC:
            pass

        return dic

    @classmethod
    def _checkB_func(cls, dic, n_board, n_channel):
        lis_checkB = []
        for i_board in range(n_board):
            for i_chn in range(n_channel):
                arr_counts = np.array(dic['board%d_channel%d' % (i_board, i_chn)])
                arr_err = np.array(dic['board%d_channel%d_err' % (i_board, i_chn)])

                if len(dic['board%d_channel%d' % (i_board, i_chn)]) != 0 and arr_counts.max() != 0:
                    err_max = arr_err[arr_counts.argmax()]
                    err_min = arr_err[arr_counts.argmin()]
                    alpha = (arr_counts.max() - arr_counts.min()) / np.sqrt(err_max * err_max + err_min * err_min)
                    lis_checkB.append(alpha)
                else:
                    lis_checkB.append(0)

        return lis_checkB

    @classmethod
    def _get_tri_up(cls, tri_up_info):
        tri_up = True
        if str(tri_up_info).lower().__contains__('up') or int(tri_up_info) == 1:
            tri_up = True
        elif str(tri_up_info).lower().__contains__('down') or int(tri_up_info) == 0:
            tri_up = False
        else:
            raise ValueError('In DetInfo, TriBarUp is not up nor down.')
        return tri_up

    @classmethod
    def _get_xy_board(cls, xy_info):
        xy = True
        if str(xy_info).lower().__contains__('x') or int(xy_info) == 1:
            xy = True
        elif str(xy_info).lower().__contains__('y') or int(xy_info) == 0:
            xy = False
        else:
            raise ValueError('In DetInfo, xy_board is not x nor y.')
        return xy

    @classmethod
    def _get_origin(cls, origin_info, if_x_board):
        if if_x_board:
            origin = origin_info['Origin_x']
        else:
            origin = origin_info['Origin_y']
        return origin

    def _calc_trigger_12(self, df, df_det_pro, n_boards=4, n_bars=32):
        column_names = ['chn_%d' % i for i in range(n_bars)]
        mtx_data, dic_data = CommonUtils.cvt_dataframe_to_numpy(df[column_names])
        mtx_data = mtx_data.reshape(-1, n_boards, n_bars)
        lis_mtx_pos = []

        for i_board in range(n_boards):
            bar_interval = df_det_pro['barInterval'].iloc[i_board]
            bar_height = df_det_pro['barHeight'].iloc[i_board]
            if_x_board = self._get_xy_board(df_det_pro['xy_board'].iloc[i_board])
            tri_bar_up = self._get_tri_up(df_det_pro['TriBarUp'].iloc[i_board])
            x_origin = self._get_origin(df_det_pro[['Origin_x', 'Origin_y']].iloc[i_board], if_x_board)
            z_origin = df_det_pro['Origin_z'].iloc[i_board]

            mtx_dac = np.zeros((mtx_data.shape[0], 4))  # e1, e2, trigger_count, sign

            mtx_board = mtx_data[:, i_board, :]
            mtx_board = np.insert(mtx_board, n_bars, values=0, axis=1)
            first_nonzero_index = (mtx_board > 0).argmax(axis=1)

            # fill mtx_dac
            mtx_dac[:, 0] = mtx_board[np.arange(mtx_board.shape[0]), first_nonzero_index]
            mtx_dac[:, 1] = mtx_board[np.arange(mtx_board.shape[0]), first_nonzero_index + 1]
            mtx_dac[mtx_dac[:, 1] > 0, 2] = 2
            mtx_dac[mtx_dac[:, 1] == 0, 2] = 1
            mtx_dac[first_nonzero_index % 2 == 0, 3] = tri_bar_up
            mtx_dac[first_nonzero_index % 2 == 1, 3] = not tri_bar_up
            mtx_dac[mtx_dac[:, 3] == 0, 3] = -1

            # calculate position
            mtx_pos = np.zeros((mtx_data.shape[0], 3))  # x,z,trigger_count
            energy_ratio = 1. * mtx_dac[:, 1] / (mtx_dac[:, 0] + mtx_dac[:, 1])
            mtx_pos[:, 0] = x_origin + bar_interval / 2 * (first_nonzero_index + 1 + energy_ratio)
            mtx_pos[:, 1] = z_origin + mtx_dac[:, 3] * (0.5 - energy_ratio) * bar_height
            mtx_pos[:, 2] = mtx_dac[:, 2]

            lis_mtx_pos.append(mtx_pos)

        return lis_mtx_pos

    def _calc_trigger_34(self, df, df_det_pro: pd.DataFrame, n_boards=4, n_bars=32):
        column_names = ['chn_%d' % i for i in range(n_bars)]
        mtx_data, dic_data = CommonUtils.cvt_dataframe_to_numpy(df[column_names])
        mtx_data = mtx_data.reshape(-1, n_boards, n_bars)
        lis_mtx_pos = []

        for i_board in range(n_boards):
            bar_interval = df_det_pro['barInterval'].iloc[i_board]
            bar_height = df_det_pro['barHeight'].iloc[i_board]
            if_x_board = self._get_xy_board(df_det_pro['xy_board'].iloc[i_board])
            tri_bar_up = self._get_tri_up(df_det_pro['TriBarUp'].iloc[i_board])
            x_origin = self._get_origin(df_det_pro[['Origin_x', 'Origin_y']].iloc[i_board], if_x_board)
            z_origin = df_det_pro['Origin_z'].iloc[i_board]

            mtx_dac = np.zeros((mtx_data.shape[0], 5))  # e1, e2, trigger_count, sign

            mtx_board = mtx_data[:, i_board, :]
            trigger_count = (mtx_board > 0).sum(axis=1)
            mtx_board = np.insert(mtx_board, n_bars, values=0, axis=1)
            first_nonzero_index = (mtx_board > 0).argmax(axis=1)

            # fill 1 2 trigger
            mtx_dac[trigger_count < 3, 0] = \
                mtx_board[np.arange(mtx_board.shape[0])[trigger_count < 3], first_nonzero_index[trigger_count < 3]]
            mtx_dac[trigger_count < 3, 1] = \
                mtx_board[np.arange(mtx_board.shape[0])[trigger_count < 3], first_nonzero_index[trigger_count < 3] + 1]
            mtx_dac[(mtx_dac[:, 1] > 0) & (trigger_count < 3), 2] = 2
            mtx_dac[(mtx_dac[:, 1] == 0) & (trigger_count < 3), 2] = 1
            mtx_dac[(first_nonzero_index % 2 == 0) & (trigger_count < 3), 3] = tri_bar_up
            mtx_dac[(first_nonzero_index % 2 == 1) & (trigger_count < 3), 3] = not tri_bar_up
            mtx_dac[(mtx_dac[:, 3] == 0) & (trigger_count < 3), 3] = -1
            mtx_dac[trigger_count < 3, 4] = first_nonzero_index[trigger_count < 3]

            # fill 3 4 trigger
            mtx_dac_34 = mtx_dac[trigger_count > 2]
            mtx_board_34 = mtx_board[trigger_count > 2]
            max_ind = mtx_board[trigger_count > 2].argmax(axis=1)
            first_nonzero_ind = (mtx_board[trigger_count > 2] > 0).argmax(axis=1)
            last_nonzero_ind = mtx_board.shape[1] - 1 - (mtx_board[trigger_count > 2][:, ::-1] > 0).argmax(axis=1)

            continue_four_trigger = last_nonzero_ind - first_nonzero_ind == 3
            mtx_dac_34[continue_four_trigger, 1] = mtx_board_34[
                                                       np.arange(mtx_board_34.shape[0])[continue_four_trigger],
                                                       first_nonzero_ind[continue_four_trigger]] \
                                                   + mtx_board_34[
                                                       np.arange(mtx_board_34.shape[0])[continue_four_trigger],
                                                       first_nonzero_ind[continue_four_trigger] + 1]
            mtx_dac_34[continue_four_trigger, 0] = mtx_board_34[
                                                       np.arange(mtx_board_34.shape[0])[continue_four_trigger],
                                                       first_nonzero_ind[continue_four_trigger] + 2] \
                                                   + mtx_board_34[
                                                       np.arange(mtx_board_34.shape[0])[continue_four_trigger],
                                                       first_nonzero_ind[continue_four_trigger] + 3]
            mtx_dac_34[continue_four_trigger, 2] = 4
            mtx_dac_34[continue_four_trigger & ((first_nonzero_ind + 1) % 2 == 0), 3] = tri_bar_up
            mtx_dac_34[continue_four_trigger & ((first_nonzero_ind + 1) % 2 == 1), 3] = not tri_bar_up
            mtx_dac_34[continue_four_trigger, 4] = first_nonzero_ind[continue_four_trigger]

            continue_three_trigger = last_nonzero_ind - first_nonzero_ind == 2
            max_first_bar = (max_ind == first_nonzero_ind) & continue_three_trigger
            mtx_dac_34[max_first_bar, 1] = \
                mtx_board_34[np.arange(mtx_board_34.shape[0])[max_first_bar], first_nonzero_ind[max_first_bar]]
            mtx_dac_34[max_first_bar, 0] = \
                mtx_board_34[np.arange(mtx_board_34.shape[0])[max_first_bar], first_nonzero_ind[max_first_bar] + 1] \
                + mtx_board_34[np.arange(mtx_board_34.shape[0])[max_first_bar], first_nonzero_ind[max_first_bar] + 2]
            mtx_dac_34[max_first_bar & (first_nonzero_ind % 2 == 0), 3] = tri_bar_up
            mtx_dac_34[max_first_bar & (first_nonzero_ind % 2 == 1), 3] = not tri_bar_up
            mtx_dac_34[max_first_bar, 4] = first_nonzero_ind[max_first_bar]

            max_second_bar = (max_ind == first_nonzero_ind + 1) & continue_three_trigger
            mtx_dac_34[max_second_bar, 1] = \
                mtx_board_34[np.arange(mtx_board_34.shape[0])[max_second_bar], max_ind[max_second_bar]]
            mtx_dac_34[max_second_bar & (max_ind % 2 == 0), 3] = tri_bar_up
            mtx_dac_34[max_second_bar & (max_ind % 2 == 1), 3] = not tri_bar_up
            mtx_dac_34[max_second_bar, 4] = first_nonzero_ind[max_second_bar]

            max_third_bar = (max_ind == last_nonzero_ind) & continue_three_trigger
            mtx_dac_34[max_third_bar, 1] = \
                mtx_board_34[np.arange(mtx_board_34.shape[0])[max_third_bar], first_nonzero_ind[max_third_bar]] \
                + mtx_board_34[np.arange(mtx_board_34.shape[0])[max_third_bar], first_nonzero_ind[max_third_bar] + 1]
            mtx_dac_34[max_third_bar, 0] = \
                mtx_board_34[np.arange(mtx_board_34.shape[0])[max_third_bar], first_nonzero_ind[max_third_bar] + 2]
            mtx_dac_34[max_third_bar & ((first_nonzero_ind + 1) % 2 == 0), 3] = tri_bar_up
            mtx_dac_34[max_third_bar & ((first_nonzero_ind + 1) % 2 == 1), 3] = not tri_bar_up
            mtx_dac_34[max_third_bar, 4] = first_nonzero_ind[max_third_bar]

            mtx_dac_34[continue_three_trigger, 2] = 3
            mtx_dac_34[mtx_dac_34[:, 3] == 0, 3] = -1

            # add mtx_dac_34 into mtx_dac
            mtx_dac[trigger_count > 2] = mtx_dac_34

            # calculate position
            mtx_pos = np.zeros((mtx_data.shape[0], 3))  # x,z,trigger_count
            energy_ratio = 1. * mtx_dac[:, 1] / (mtx_dac[:, 0] + mtx_dac[:, 1])
            mtx_pos[:, 0] = x_origin + bar_interval / 2 * (first_nonzero_index + 1 + energy_ratio)
            mtx_pos[:, 1] = z_origin + mtx_dac[:, 3] * (0.5 - energy_ratio) * bar_height
            mtx_pos[:, 2] = mtx_dac[:, 2]

            lis_mtx_pos.append(mtx_pos)

        return lis_mtx_pos

    def _get_raw_dataframe(self, lis_mtx_pos: list, df_det_pro):
        n_board = len(lis_mtx_pos)

        layer_x = layer_y = 1
        df_data = pd.DataFrame()
        for i_board in range(n_board):
            mtx_pos = lis_mtx_pos[i_board]
            if_x_board = self._get_xy_board(df_det_pro['xy_board'].iloc[i_board])

            if if_x_board:
                df_data['xo%d' % layer_x] = mtx_pos[:, 0]
                df_data['zx%d' % layer_x] = mtx_pos[:, 1]
                df_data['xNBar%d' % layer_x] = mtx_pos[:, 2].astype('int')
                layer_x += 1
            else:
                df_data['yo%d' % layer_y] = mtx_pos[:, 0]
                df_data['zy%d' % layer_y] = mtx_pos[:, 1]
                df_data['yNBar%d' % layer_y] = mtx_pos[:, 2].astype('int')
                layer_y += 1
        return df_data

    @classmethod
    def _get_thetaxy_trk(cls, df_data):
        thetax_trk = TMath.ATan2((df_data.xo1 - df_data.xo2), (df_data.zx1 - df_data.zx2))
        thetay_trk = TMath.ATan2((df_data.yo1 - df_data.yo2), (df_data.zy1 - df_data.zy2))
        return thetax_trk, thetay_trk

    @classmethod
    def _get_thetaxy_thetaphi(cls, df_data, rotAnglex, rotAngley):
        v_det_coord1 = TVector3(df_data.X1, df_data.Y1, df_data.Z1)
        v_det_coord2 = TVector3(df_data.X2, df_data.Y2, df_data.Z2)

        v_det_angle = TVector3(TMath.Tan(df_data.thetax_trk), TMath.Tan(df_data.thetay_trk), 1)
        v_det_coord = v_det_coord1 - v_det_coord2
        v_car_coord = TVector3(v_det_coord)

        if rotAnglex > 1e-10 and rotAngley > 1e-10:
            raise ValueError('rotAnglex>1e-10 && rotAngley>1e-10.')

        if rotAnglex < 1e-10:
            v_car_coord.RotateX(rotAngley)

        if rotAngley < 1e-10:
            v_car_coord.RotateY(rotAnglex)

        thetax = thetay = -10
        if rotAnglex < 1e-10:
            thetax = TMath.ATan2(v_car_coord(0), v_car_coord(2))
            # thetay = TMath.ATan2(v_car_coord(1), v_car_coord(2))
            # I am not sure is df_data.thetay_trk + rotAngley the same as TMath.ATan2(v_car_coord(1), v_car_coord(2)) here
            # It seems to be the same
            thetay = df_data.thetay_trk + rotAngley

        if rotAngley < 1e-10:
            thetay = TMath.ATan2(v_car_coord(1), v_car_coord(2))
            # thetax = TMath.ATan2(v_car_coord(0), v_car_coord(2))
            thetax = df_data.thetax_trk + rotAnglex

        theta = v_car_coord.Theta()
        phi = v_car_coord.Phi()
        return thetax, thetay, theta, phi

    def _get_whole_dataframe_4plane(self, df_data, Z1_aver, Z2_aver, rotAnglex, rotAngley, check_angle=0):
        df_data['Z1'] = Z1_aver
        df_data['Z2'] = Z2_aver

        df_data.eval("""X1 = (Z1 - zx1) * (xo1 - xo2) / (zx1 - zx2) + xo1
        X2 = (Z2 - zx1) * (xo1 - xo2) / (zx1 - zx2) + xo1
        Y1 = (Z1 - zy1) * (yo1 - yo2) / (zy1 - zy2) + yo1
        Y2 = (Z2 - zy1) * (yo1 - yo2) / (zy1 - zy2) + yo1""", inplace=True)

        df_data[["thetax_trk", "thetay_trk"]] = df_data.apply(self._get_thetaxy_trk, axis=1, result_type="expand")
        df_data[["thetax", "thetay", "theta", "phi"]] = \
            df_data.apply(self._get_thetaxy_thetaphi, rotAnglex=rotAnglex, rotAngley=rotAngley, axis=1,
                          result_type="expand")

        if not check_angle > 0:
            return df_data

        df_data = df_data[abs(df_data['thetay_trk']) < check_angle]
        df_data = df_data[abs(df_data['thetax_trk']) < check_angle]

        return df_data

    @classmethod
    def _get_Z_aver_4plane(cls, det_info, det):
        if det.get_Z_aver() is None:
            print('Do not get ZAver, use the default ones.')
            Z1_aver = det_info['Origin_z'][0] + det_info['Origin_z'][1]
            Z2_aver = det_info['Origin_z'][2] + det_info['Origin_z'][3]
            return [Z1_aver, Z2_aver]
        else:
            [Z1_aver, Z2_aver] = det.get_Z_aver()
            print('Get ZAver1 = %f, ZAver2 = %f' % (Z1_aver, Z2_aver))
            return [Z1_aver, Z2_aver]

    def CalcHitPoints_4Plane(self, h5_fname, directory, root_fname, det, calc_multi_trigger=True,
                             use_det_meas_time=False):
        # noinspection DuplicatedCode
        h5_name = CommonUtils.splitext_h5_fname(h5_fname)

        n_rows = CommonUtils.get_dataframe_total_lines('%s.h5' % os.path.join(directory, h5_name))

        print('\rCalcHitPoints_Plane for %s' % h5_name, end='')

        detpro_fname = det.get_det_pro_fname()
        lis_rot_angle = det.get_rotAngle()
        if use_det_meas_time:
            meas_time = det.get_meas_time()
        else:
            meas_time = 0
        det_pro = pd.read_csv(detpro_fname)
        [Z1_aver, Z2_aver] = self._get_Z_aver_4plane(det_pro, det)
        rotAnglex = lis_rot_angle[0]
        rotAngley = lis_rot_angle[1]

        bar_interval = det_pro['barInterval'].max()
        bar_height = det_pro['barInterval'].min()
        check_angle = TMath.ATan2(bar_interval / 2, bar_height)

        df_data = None
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 100000, directory):
            df_data = i.copy()
            break
        n_boards = CommonUtils.count_n_boards(df_data)
        n_channels = CommonUtils.count_n_channels(df_data)
        n_bars = CommonUtils.count_n_bars(df_data)
        chn_col_names = ['chn_%d' % i for i in range(n_bars)]

        # get measure time
        if meas_time == 0:
            df_meas_time = self.get_df_meas_time(os.path.join(directory, h5_name + '.h5'))
            df_meas_time = df_meas_time.drop_duplicates()
            meas_time = df_meas_time['time_delta'].sum()
            meas_time = meas_time.total_seconds()

        (root_name, file_type) = os.path.splitext(root_fname)

        df_data = None
        pageIndex = 0
        count = 0
        lis_root_splits = []
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 100000, directory):
            count += i.shape[0]
            # split out the three or four trigger event
            chn_data = i[chn_col_names]
            trigger_count = (chn_data > 1).sum(axis=1)
            trigger_34 = trigger_count > 3
            trigger_34_triggerID = i[trigger_34]['triggerID'].values
            i_triggerID_index = i.set_index('triggerID')
            trigger_34_event = i_triggerID_index.loc[trigger_34_triggerID]
            trigger_12_event = i_triggerID_index.drop(index=trigger_34_triggerID)

            # calculate trigger_12 positions and fill into df_data
            lis_trigger_12 = self._calc_trigger_12(trigger_12_event, det_pro, n_boards, n_bars)
            df_data = self._get_raw_dataframe(lis_trigger_12, df_det_pro=det_pro)
            df_data = self._get_whole_dataframe_4plane(df_data, Z1_aver, Z2_aver, rotAnglex, rotAngley)

            del lis_trigger_12, trigger_12_event, trigger_count, trigger_34, trigger_34_triggerID, i_triggerID_index

            # calculate trigger_34 positions and fill into df_data
            if calc_multi_trigger:
                lis_trigger_34 = self._calc_trigger_34(trigger_34_event, det_pro, n_boards, n_bars)
                t_df_data = self._get_raw_dataframe(lis_trigger_34, df_det_pro=det_pro)
                check_angle = TMath.ATan2(det_pro['barInterval'].iloc[0] / 2, det_pro['barHeight'].iloc[0])
                t_df_data = self._get_whole_dataframe_4plane(t_df_data, Z1_aver, Z2_aver,
                                                             rotAnglex, rotAngley, check_angle=check_angle)
                df_data = pd.concat((df_data, t_df_data), axis=0)
                del lis_trigger_34, t_df_data, check_angle

            t_root_fname = '%s_%d.root' % (root_name, pageIndex)
            ROOTUtils.cvt_dataframe_to_tree(df_data, 't_tree', t_root_fname)
            lis_root_splits.append(t_root_fname)
            print(r'Save split %s' % t_root_fname, end='')

            pageIndex += 1
            print('\rCalcHitPoints_Plane for %s: %d / %d    ' % (h5_name, count, n_rows), end='')
            del i

        print('')

        ROOTUtils.merge_root_files(lis_root_splits, root_name)

        # add the measure time information
        measTime = TVectorD(1)
        measTime[0] = meas_time
        fout = TFile(root_fname, 'UPDATE')
        measTime.Write('meas_time')
        fout.Save()
        fout.Close()
