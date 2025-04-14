import os, logging
import shutil
import datetime

from DatabaseTool.dataLayer.baseCore import h5Data
from DatabaseTool import *
from DataFactory import CommonUtils, ROOTUtils

logger = logging.getLogger(__name__)


class DataframeProcessor:
    def __init__(self):
        self._cr_offline = Cormis_DataReader_offline()

    @staticmethod
    def SelectWithCondition(h5_in_fname: str, h5_out_fname: str, directory: str, dict_condition: dict):
        """
        Select events with given condition
        @param h5_in_fname: input hdf5 file name
        @param h5_out_fname: output hdf5 file absolute path
        @param directory: input hdf5 directory
        @param dict_condition: dict {start_time, end_time, start_trigger, end_trigger}
        @return:
        """
        h5_name = CommonUtils.splitext_h5_fname(h5_in_fname)
        n_rows = CommonUtils.get_dataframe_total_lines(os.path.join(directory, '%s.h5' % h5_name))

        start_time = dict_condition.get('start_time', '')
        end_time = dict_condition.get('end_time', '')

        start_trigger = dict_condition.get('start_trigger')
        end_trigger = dict_condition.get('end_trigger')

        if os.path.exists(h5_out_fname):
            os.remove(h5_out_fname)

        pageIdx = 0
        count = 0

        # convert input timedelta to datetime
        df = None
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 1000, directory):
            df = i.copy()
            break

        start_time, end_time = CommonUtils.transform_df_start_end_time(df, start_time, end_time)
        _start_time = None
        _end_time = end_time

        for i in Cormis_DataReader_hdf.getCaching(h5_name, 20000, directory):
            count += i.shape[0]
            if pageIdx == 0:
                _start_time = i['timeTag'].iloc[0]

            if_break = False
            # select with conditions
            if start_time is not None:
                i = i[i['timeTag'] >= start_time]

            if end_time is not None and i.shape[0] > 0:
                if end_time < i['timeTag'].iloc[-1]:
                    if_break = True
                i = i[i['timeTag'] < end_time]

            if start_trigger is not None and i.shape[0] > 0:
                start_trigger = int(start_trigger)
                i = i[i['triggerID'] >= start_trigger]

            if end_trigger is not None and i.shape[0] > 0:
                end_trigger = int(end_trigger)
                if end_trigger < i['triggerID'].iloc[-1]:
                    if_break = True
                i = i[i['triggerID'] < end_trigger]

            i.to_hdf(h5_out_fname, key="data", append=True)

            if i.shape[0] > 0:
                _end_time = i['timeTag'].iloc[-1]
            if if_break:
                break

            pageIdx += 1
            print('\rSelectWithCondition for %s: %d / %d    ' % (h5_name, count, n_rows), end='')
        print('')

        dict_condition['start_time'] = _start_time if start_time is None else max(_start_time, start_time)
        dict_condition['end_time'] = _end_time if end_time is None else min(_end_time, end_time)

        # edit df_meas_time
        df_meas_time = None
        try:
            df_meas_time = pd.read_hdf(os.path.join(directory, '%s.h5' % h5_name), key="meas_time")
        except:
            pass
        if df_meas_time is None or not isinstance(df_meas_time, pd.DataFrame):
            return dict_condition

        df_new = df_meas_time.copy()
        if not dict_condition['start_time'] in df_meas_time['start_time'].tolist():
            row_id = (df_new['start_time'] < dict_condition['start_time']).sum()
            df_new = df_new.iloc[row_id - 1:, :]
            df_new.iloc[0, 0] = df_new.iloc[0, 0] - (dict_condition['start_time'] - df_new.iloc[0, 1])
            df_new.iloc[0, 1] = dict_condition['start_time']
            df_new = df_new.rename(index={df_new.index[0]: df_new.index[0] + '_cut'})

        pd_end_time = df_new['start_time'] + df_new['time_delta']
        if not dict_condition['end_time'] in pd_end_time.to_list():
            keep_num = (pd_end_time < dict_condition['end_time']).sum()
            df_new = df_new.iloc[:keep_num + 1, :]
            df_new.iloc[-1, 0] = dict_condition['end_time'] - df_new.iloc[-1, 1]
            df_new = df_new.rename(index={df_new.index[-1]: df_new.index[-1] + '_cut'})

        print(df_new)
        df_new.to_hdf(os.path.join(directory, '%s.h5' % h5_name), key="meas_time")

        return dict_condition

    @staticmethod
    def PrintHdf5ToFile(h5_fname, out_fname, directory):
        """
        print dataframe in hdf5 into csv file
        @param h5_fname: input hdf5 file name
        @param out_fname: output hdf5 file absolute path
        @param directory: input hdf5 directory
        @return:
        """
        h5_name = CommonUtils.splitext_h5_fname(h5_fname)

        pd.set_option('display.max_columns', None)

        if os.path.exists(out_fname):
            os.remove(out_fname)

        use_head = True
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 100000, directory):
            if use_head:
                print(i)
            i.to_csv(out_fname, mode='a', header=use_head)
            use_head = False

    @staticmethod
    def ReplaceChannel(h5_in_fname, h5_out_fname, directory, relation_info_fname):
        """
        fill dac information of channels with given bar-channel relation
        @param h5_in_fname: input hdf5 file name
        @param h5_out_fname: output hdf5 file absolute path
        @param directory: input hdf5 directory
        @param relation_info_fname: csv file of bar-channel relation
        @return:
        """
        h5_name = CommonUtils.splitext_h5_fname(h5_in_fname)

        n_rows = CommonUtils.get_dataframe_total_lines(os.path.join(directory, '%s.h5' % h5_name))

        # read relation_info as a dataframe
        df_rel_info = pd.read_csv(relation_info_fname)

        df_data = None
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 100000, directory):
            df_data = i.copy()
            break
        n_boards = CommonUtils.count_n_boards(df_data)
        n_channels = CommonUtils.count_n_channels(df_data)
        n_bars = CommonUtils.count_n_bars(df_data)

        if os.path.exists(h5_out_fname):
            os.remove(h5_out_fname)

        pageIdx = 0
        count = 0
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 100000, directory):
            count += i.shape[0]
            chn_col_names = i.columns[i.columns.str.contains(r'chn_\d+')]

            chn_data = i[chn_col_names]
            mtx_data, dic_data = CommonUtils.cvt_dataframe_to_numpy(chn_data)

            mtx_data = mtx_data.reshape(-1, n_boards, len(dic_data))
            mtx_data_replace = np.full([mtx_data.shape[0], mtx_data.shape[1], n_channels], 0)

            for i_board in range(n_boards):
                for i_bar in range(n_bars):
                    t_chn = df_rel_info.iloc[i_board, i_bar]
                    mtx_data_replace[:, i_board, i_bar] = mtx_data[:, i_board, t_chn]

            mtx_data_replace = mtx_data_replace.reshape(-1, n_channels)

            i[chn_col_names] = mtx_data_replace

            i.to_hdf(h5_out_fname, key="data", append=True)
            pageIdx += 1
            print('\rReplaceChannel for %s: %d / %d' % (h5_in_fname, count, n_rows), end='')
        print('')

        try:
            df_meas_time = pd.read_hdf(os.path.join(directory, '%s.h5' % h5_name), key="meas_time")
            df_meas_time.to_hdf(h5_out_fname, key="meas_time", append=True)
        except:
            pass

    @classmethod
    def _select_all_hit(cls, data: pd.DataFrame, n_boards: int):
        """
        select events that hit all boards
        @param data: pandas.DataFrame
        @param n_boards: int; number of boards
        @return:
        """
        trigger_count = data.groupby(['triggerID'])['triggerID'].count()
        hit_all_filter = trigger_count.index.values[trigger_count.values == n_boards]
        selected_data = data.set_index('triggerID').loc[hit_all_filter, :]
        selected_data = selected_data.reset_index()

        return selected_data

    @staticmethod
    def FiltrateAllHit(h5_in_fname, h5_out_fname, directory, log_fname=None, check_continue_trigger=False):
        """
        filtrate events that hit all boards
        @param h5_in_fname: input hdf5 file name
        @param h5_out_fname: output hdf5 file absolute path
        @param directory: input hdf5 directory
        @param log_fname: log file absolute path
        @param check_continue_trigger: boolean; if check triggered bars are adjoining
        @return:
        """
        h5_name = CommonUtils.splitext_h5_fname(h5_in_fname)

        total_hit = CommonUtils.get_dataframe_total_lines('%s.h5' % os.path.join(directory, h5_name))

        n_boards = -1
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 5000, directory):
            n_boards: int = CommonUtils.count_n_boards(i)
            break

        if os.path.exists(h5_out_fname):
            os.remove(h5_out_fname)

        # initialize
        pageIdx = 0
        dic_trigger_count = {}
        for i_board in range(n_boards + 1):
            dic_trigger_count['trigger_%d_board' % i_board] = 0

        count = 0
        remain_count = 0
        start_time = 0
        end_time = 0
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 20000, directory):
            if pageIdx == 0:
                start_time = i['timeTag'].iloc[0]
            end_time = i['timeTag'].iloc[-1]

            count += i.shape[0]
            # delete all 0 rows
            chn_col_names = i.columns[i.columns.str.contains(r'chn_\d+')]
            chn_data = i[chn_col_names]
            if_trigger = chn_data.sum(axis=1) > 0
            i = i[if_trigger]

            # count trigger number
            trigger_count = i.groupby(['triggerID'])['triggerID'].count()
            hit_all_filter = trigger_count.index.values[trigger_count.values == n_boards]
            selected_data = i.set_index('triggerID').loc[hit_all_filter, :]

            # select continue trigger event
            if check_continue_trigger:
                mtx_select, dic_select = CommonUtils.cvt_dataframe_to_numpy(selected_data[chn_col_names])
                first_nonzero_ind = (mtx_select > 0).argmax(axis=1)
                last_nonzero_ind = mtx_select.shape[1] - 1 - (mtx_select[:, ::-1] > 0).argmax(axis=1)
                trigger_interval = last_nonzero_ind - first_nonzero_ind
                trigger_bar_count = (mtx_select > 0).sum(axis=1)
                if_continue = (trigger_interval + 1) == trigger_bar_count
                if_continue = if_continue.reshape(-1, n_boards)
                continue_trigger = if_continue.sum(axis=1) == n_boards
                for i_board in range(n_boards):
                    if_continue[:, i_board] = continue_trigger
                if_continue = if_continue.reshape(-1)

                selected_data = selected_data[if_continue]

            selected_data = selected_data.reset_index()

            dic_trigger_count['trigger_%d_board' % n_boards] += selected_data.shape[0]
            remain_count += selected_data.shape[0]

            selected_data.to_hdf(h5_out_fname, key="data", append=True)

            # count trigger board
            for i_board in range(n_boards):
                hit_all_filter = trigger_count.index.values[trigger_count.values == i_board]
                t_data = i.set_index('triggerID').loc[hit_all_filter, :]
                dic_trigger_count['trigger_%d_board' % i_board] += t_data.shape[0]

            pageIdx += 1
            print('\rFiltrateAllHit for %s: %d / %d     ' % (h5_name, count, total_hit), end='')
        print('')

        try:
            df_meas_time = pd.read_hdf(os.path.join(directory, '%s.h5' % h5_name), key="meas_time")
            df_meas_time.to_hdf(h5_out_fname, key="meas_time", append=True)
        except:
            pass

        time_interval = end_time - start_time
        logger.info('%s\nhit_all_board / total = %d / %d = %.3f\n'
                    'trigger rate (dataframe number of rows / time interval, unit: per hour): %d / %s = %.2f' %
                    (h5_name, dic_trigger_count['trigger_%d_board' % n_boards], total_hit,
                     1. * dic_trigger_count['trigger_%d_board' % n_boards] / total_hit,
                     remain_count, time_interval, 3600. * remain_count / time_interval.total_seconds()))
        if log_fname is not None:
            print('\n', datetime.datetime.now(), file=open(log_fname, "a"))
            print('trigger rate (dataframe number of rows / time interval, unit: per hour): %d / %s = %.2f'
                  % (remain_count, time_interval, 3600. * remain_count / time_interval.total_seconds()),
                  file=open(log_fname, "a"))
            print('\nFiltrateAllHit for %s:' % h5_in_fname + '.h5', file=open(log_fname, "a"))
            for i_board in range(n_boards + 1):
                i_key = 'trigger_%d_board' % i_board
                n_trigger = dic_trigger_count[i_key]
                print('\n%s / total = %d / %d = %.3f\n' %
                      (i_key, n_trigger, total_hit,
                       1. * n_trigger / total_hit), file=open(log_fname, "a"))
            print('Save as %s' % h5_out_fname, file=open(log_fname, "a"))

    @staticmethod
    def DeleteGivenChannels(h5_in_fname, h5_out_fname, directory, lis_del_chn_number, log_fname=None):
        """
        set dac information of certain channels to zero
        @param h5_in_fname: input hdf5 file name
        @param h5_out_fname: output hdf5 file absolute path
        @param directory: input hdf5 directory
        @param lis_del_chn_number: list of channel numbers to be deleted
        @param log_fname: log file absolute path
        @return:
        """
        h5_name = CommonUtils.splitext_h5_fname(h5_in_fname)

        n_rows = CommonUtils.get_dataframe_total_lines('%s.h5' % os.path.join(directory, h5_name))

        if os.path.exists(h5_out_fname):
            os.remove(h5_out_fname)

        pageIdx = 0
        count = 0
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 1000000, directory):
            count += i.shape[0]
            for i_num in lis_del_chn_number:
                [i_board, i_chn] = i_num.split(':')
                i_chn_name = 'chn_' + str(i_chn)

                i = i.set_index('boardID')
                i.loc[int(i_board), i_chn_name] = 0
                i = i.reset_index()

            i.to_hdf(h5_out_fname, key="data", append=True)

            pageIdx += 1
            print('\rDeleteGivenChannels for %s: %d / %d    ' % (h5_name, count, n_rows), end='')
        print('')

        try:
            df_meas_time = pd.read_hdf(os.path.join(directory, '%s.h5' % h5_name), key="meas_time")
            df_meas_time.to_hdf(h5_out_fname, key="meas_time", append=True)
        except:
            pass

        logger.info('delete channel %s in %s' %
                    (' '.join([str(j) for j in lis_del_chn_number]), h5_in_fname + '.h5'))
        if log_fname is not None:
            print(datetime.datetime.now(), 'DeleteGivenChannel channel %s from %s\nSave as %s' %
                  (' '.join([str(j) for j in lis_del_chn_number]), h5_in_fname + '.h5', h5_out_fname),
                  file=open(log_fname, "a"))

    @staticmethod
    def MergeHdf5(lis_h5_fname, h5_out_fname: str, log_fname=None):
        """
        merge dataframes from several hdf files, record the measured time and save in one hdf5 file
        @param lis_h5_fname: list of input hdf5 files' absolute path
        @param h5_out_fname: output hdf5 file absolute path
        @param log_fname: log file absolute path
        @return:
        """
        lis_h5_fname = [CommonUtils.splitext_h5_fname(i) for i in lis_h5_fname]

        if os.path.exists(h5_out_fname):
            os.remove(h5_out_fname)

        pageIdx = 0
        all_count = 0
        df_meas_time = pd.DataFrame(columns=['time_delta', 'start_time'])
        trigger_index = 0
        for h5_fname in lis_h5_fname:
            n_rows = CommonUtils.get_dataframe_total_lines('%s.h5' % h5_fname)

            (h5_dir, h5_name) = os.path.split(h5_fname)
            if log_fname is not None:
                print(datetime.datetime.now(),
                      '\nMerge %s, all counts = %d + %d' % (h5_fname, all_count, n_rows),
                      file=open(log_fname, "a"))

            count = 0
            start_time = end_time = 0
            last_trigger = 0
            for i in Cormis_DataReader_hdf.getCaching(h5_name, 100000, h5_dir):
                if start_time == 0:
                    start_time = i['timeTag'].iloc[0]
                if 'id' in i.columns:
                    i = i.drop(columns=['id'])
                i.triggerID = i['triggerID'] + trigger_index + 1
                i.to_hdf(h5_out_fname, key="data", append=True)

                count += i.shape[0]
                print('\rMerge %s, all counts = %d +  %d / %d    ' % (h5_name, all_count, count, n_rows), end='')

                pageIdx += 1
                end_time = i['timeTag'].iloc[-1]
                last_trigger = i.triggerID.iloc[-1]

            try:
                df_meas_time_from_file = pd.read_hdf(os.path.join(h5_dir, h5_name + '.h5'), key="meas_time")
                df_meas_time = pd.concat([df_meas_time, df_meas_time_from_file], axis=0)
            except:
                df_meas_time.loc[h5_name] = [end_time - start_time, start_time]

            all_count += n_rows
            trigger_index = last_trigger

        # df_meas_time.convert_dtypes() can avoid the error "Cannot serialize the column [time_delta]
        # because its data contents are not [string] but [timedelta] object dtype". However, it deprecated. So sad
        df_meas_time['start_time'] = pd.to_datetime(df_meas_time['start_time'])
        df_meas_time['time_delta'] = pd.to_timedelta(df_meas_time['time_delta'])
        df_meas_time.to_hdf(h5_out_fname, key="meas_time", append=True)

    @classmethod
    def _checkIndex(cls, rawData: h5Data):
        for i in range(rawData.index.shape[0]):
            if rawData.index[i] == 0:
                rawData.index.resize((i,))
                return i

    @classmethod
    def _checkDataLength(cls, rawData: h5Data):
        startIdx = rawData.index[rawData.index.shape[0] - 1]
        stopIdx = rawData.index[0]
        for i in range(startIdx, stopIdx):
            lineData = rawData.data[i]
            lineTime = rawData.timeTag[i]
            lineSum = np.sum(lineData)
            if (lineSum == 0) or (lineTime == 0):
                print("\n line data:{}\n line time:{}\n idx:{}; break here.".format(lineData, lineTime, i))
                rawData.index[0] = i - 1
                return i
            else:
                tLen = stopIdx - startIdx
                per = (i + 1) / tLen
                print("\r \033[31m{}\033[0m/\033[32m{}\033[0m\t{}\tsum:{};timeTage:{}".format(i, tLen, per, lineSum,
                                                                                              lineTime), end="")

    def _fixH5Data(self, path: str):
        rawData = h5Data(path, "r+")
        print(rawData.data.shape)
        print(rawData.timeTag.shape)
        print(rawData.index[:])
        print("resize rawData index length from {} to {};".format(rawData.index.shape[0], self._checkIndex(rawData)))
        print(rawData.index[:])
        print(
            "resize rawData data length from {} to {};".format(rawData.index.shape[0], self._checkDataLength(rawData)))
        # for i in range(rawData.index.shape[0]):
        #     print(rawData.getData(i))

    def FixHdf5Data(self, h5_in_fname, directory):
        """
        fix broken hdf5
        @param h5_in_fname: list of input hdf5 files' name
        @param directory: input hdf5 directory
        @return:
        """
        h5_name = CommonUtils.splitext_h5_fname(h5_in_fname)
        h5_fname = os.path.join(directory, h5_name + '.h5')

        if os.path.exists(h5_fname):
            shutil.copy(h5_fname, os.path.join(directory, h5_name + '_bkup.h5'))

        self._fixH5Data(h5_fname)

    def AddKeyMeasTime(self, h5_in_fname, directory):
        """
        add key 'meas_time' (measurement starting time and duration) into a hdf file
        @param h5_in_fname: input hdf5 file name
        @param directory: input hdf5 directory
        @return:
        """
        h5_name = CommonUtils.splitext_h5_fname(h5_in_fname)
        h5_fname = os.path.join(directory, h5_name + '.h5')

        try:
            df_meas_time = pd.read_hdf(h5_fname, key="meas_time")
            logger.info('%s already has meas_time:\n%s' % (h5_fname, df_meas_time))
        except:
            start_time = self._cr_offline.getStartTime(h5_fname)
            end_time = self._cr_offline.getStopTime(h5_fname)

            meas_time = end_time - start_time
            df_meas_time = pd.DataFrame(columns=['time_delta', 'start_time'])
            df_meas_time.loc[h5_name] = [meas_time, start_time]

            print(df_meas_time)

            try:
                df_meas_time.to_hdf(h5_fname, key="meas_time", mode='a', append=True)
            except:
                h5_fname = os.path.join(directory, 'addKey_' + h5_name + '.h5')
                if os.path.exists(h5_fname):
                    os.remove(h5_fname)
                print('Convert old hdf5 file to new format', h5_fname)

                for i in Cormis_DataReader_hdf.getCaching(h5_name, 1000000, directory):
                    i.to_hdf(h5_fname, key="data", append=True)
                df_meas_time.to_hdf(h5_fname, key="meas_time", mode='a', append=True)

    @staticmethod
    def CutWithThreshold(h5_in_fname, h5_out_fname, directory, threshold_csv_fname):
        """
        set dac lower than threshold to zero
        @param h5_in_fname: input hdf5 file name
        @param h5_out_fname: output hdf5 file absolute path
        @param directory: input hdf5 directory
        @param threshold_csv_fname: csv file absolute path of threshold
        @return:
        """
        h5_name = CommonUtils.splitext_h5_fname(h5_in_fname)

        n_rows = CommonUtils.get_dataframe_total_lines(os.path.join(directory, '%s.h5' % h5_name))

        df_threshold = pd.read_csv(threshold_csv_fname)
        mtx_threshold = df_threshold.to_numpy()

        n_boards = df_threshold.shape[0]
        n_channels = df_threshold.shape[1]

        if os.path.exists(h5_out_fname):
            os.remove(h5_out_fname)

        pageIdx = 0
        count = 0
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 1000000, directory):
            this_count = i.shape[0]

            chn_col_names = i.columns[i.columns.str.contains(r'chn_\d+')]

            chn_data = i[chn_col_names]
            mtx_data, dic_data = CommonUtils.cvt_dataframe_to_numpy(chn_data)

            mtx_data = mtx_data.reshape(-1, n_boards, len(dic_data))
            mtx_data_cut = mtx_data
            mtx_data_cut[mtx_data < mtx_threshold] = 0

            mtx_data_cut = mtx_data_cut.reshape(-1, n_channels)
            i[chn_col_names] = mtx_data_cut
            i.to_hdf(h5_out_fname, key="data", append=True)

            pageIdx += 1
            count += this_count
            print('\rCut %d / %d for %s    ' % (count, n_rows, h5_name), end='')
        print('')

        try:
            df_meas_time = pd.read_hdf(os.path.join(directory, '%s.h5' % h5_name), key="meas_time")
            df_meas_time.to_hdf(h5_out_fname, key="meas_time", append=True)
        except:
            pass

    @staticmethod
    def SubtractBaseline(h5_in_fname, h5_out_fname, directory, baseline_csv_fname):
        """
        substrate dac with baseline
        @param h5_in_fname: input hdf5 file name
        @param h5_out_fname: output hdf5 file absolute path
        @param directory: input hdf5 directory
        @param baseline_csv_fname: csv file absolute path of baseline
        @return:
        """
        h5_name = CommonUtils.splitext_h5_fname(h5_in_fname)

        n_rows = CommonUtils.get_dataframe_total_lines(os.path.join(directory, '%s.h5' % h5_name))

        df_baseline = pd.read_csv(baseline_csv_fname)
        mtx_baseline = df_baseline.to_numpy()

        n_boards = df_baseline.shape[0]
        n_channels = df_baseline.shape[1]

        if os.path.exists(h5_out_fname):
            os.remove(h5_out_fname)

        pageIdx = 0
        count = 0
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 1000000, directory):
            count += i.shape[0]
            chn_col_names = i.columns[i.columns.str.contains(r'chn_\d+')]

            chn_data = i[chn_col_names]
            mtx_data, dic_data = CommonUtils.cvt_dataframe_to_numpy(chn_data)

            mtx_data = mtx_data.reshape(-1, n_boards, len(dic_data))
            mtx_data_sub = mtx_data
            mtx_data_sub = mtx_data_sub - mtx_baseline

            mtx_data_sub = mtx_data_sub.reshape(-1, n_channels)
            mtx_data_sub[mtx_data_sub < 0] = 0
            i[chn_col_names] = mtx_data_sub
            i.to_hdf(h5_out_fname, key="data", append=True)

            pageIdx += 1
            print('\rSubtractBaseline for %s: %d / %d     ' % (h5_name, count, n_rows), end='')
        print('')

        try:
            df_meas_time = pd.read_hdf(os.path.join(directory, '%s.h5' % h5_name), key="meas_time")
            df_meas_time.to_hdf(h5_out_fname, key="meas_time", append=True)
        except:
            pass

    @staticmethod
    def SeparateFalseTrigger_Plane(h5_in_fname, h5_out_fname, directory, log_fname=None):
        """
        separate false trigger events to key 'false_trigger'
        @param h5_in_fname: input hdf5 file name
        @param h5_out_fname: output hdf5 file absolute path
        @param directory: input hdf5 directory
        @param log_fname: log file absolute path
        @return:
        """
        h5_name = CommonUtils.splitext_h5_fname(h5_in_fname)

        all_count = CommonUtils.get_dataframe_total_lines('%s.h5' % os.path.join(directory, h5_name))

        if os.path.exists(h5_out_fname):
            os.remove(h5_out_fname)

        pageIdx_data = pageIdx_false = 0
        data_count = 0
        false_count = 0
        start_time = 0
        end_time = 0
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 1000000, directory):
            print('\rSeparateFalseTrigger_Plane for %s: %d / %d     '
                  % (h5_name, data_count + false_count + i.shape[0], all_count), end='  ')
            if pageIdx_data == 0:
                start_time = i['timeTag'].iloc[0]
            end_time = i['timeTag'].iloc[-1]

            chn_col_names = i.columns[i.columns.str.contains(r'chn_\d+')]

            chn_data = i[chn_col_names]
            sum_trigger_bar = (chn_data > 0).sum(axis=1)

            more_than_four = i[sum_trigger_bar > 4]
            one2four_bars = i[sum_trigger_bar <= 4]

            if one2four_bars.shape[0] > 0:
                one2four_bars.to_hdf(h5_out_fname, key="data", append=True)
                pageIdx_data += 1
                data_count += one2four_bars.shape[0]

            if more_than_four.shape[0] > 0:
                more_than_four.to_hdf(h5_out_fname, key="false_trigger", append=True)
                pageIdx_false += 1
                false_count += more_than_four.shape[0]

            print('Has separated data_count %d + false_count %d / total_count %d for %s'
                  % (data_count, false_count, all_count, h5_name), end='')
        print('')

        try:
            df_meas_time = pd.read_hdf(os.path.join(directory, '%s.h5' % h5_name), key="meas_time")
            df_meas_time.to_hdf(h5_out_fname, key="meas_time", append=True)
        except:
            pass

        time_interval = end_time - start_time

        logger.info('SeparateFalseTrigger_Plane for %s\n'
                    '1-4 bar trigger: %d / %d = %.3f\n'
                    'more than 4 bars trigger: %d / %d = %.3f\n'
                    'trigger rate (dataframe number of rows / time interval, unit: per hour): %d / %s = %.2f'
                    % (h5_in_fname, data_count, all_count, 1. * data_count / all_count,
                       false_count, all_count, 1. * false_count / all_count,
                       data_count, time_interval, 3600. * data_count / time_interval.total_seconds()))
        if log_fname is not None:
            print('\n\n', datetime.datetime.now(), file=open(log_fname, 'a'))
            print('SeparateFalseTrigger_Plane for %s\n1-4 bar trigger: %d / %d = %.3f\n'
                  'more than 4 bars trigger: %d / %d = %.3f\n'
                  'trigger rate (dataframe number of rows / time interval, unit: per hour): %d / %s = %.2f'
                  % (os.path.join(directory, h5_in_fname), data_count, all_count, 1. * data_count / all_count,
                     false_count, all_count, 1. * false_count / all_count,
                     data_count, time_interval, 3600. * data_count / time_interval.total_seconds()),
                  file=open(log_fname, 'a'))

    @staticmethod
    def AddUniformTimeTag(h5_in_fname, h5_out_fname, directory, time_delta: int, begin_time=0):
        h5_name = CommonUtils.splitext_h5_fname(h5_in_fname)

        n_rows = CommonUtils.get_dataframe_total_lines(os.path.join(directory, '%s.h5' % h5_name))
        delta = datetime.timedelta(seconds=time_delta / n_rows)

        count = 0
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 1000000, directory):
            i = i.reset_index()
            multi = i.index.to_numpy() + count
            i['timeTag'] = multi * delta
            i['timeTag'] += pd.Timestamp(begin_time, unit='s')

            i.to_hdf(h5_out_fname, key="data", append=True)

            count += i.shape[0]
            print('\rAddUniformTimeTag for %s: %d / %d    ' % (h5_name, count, n_rows), end='')

        print('')

    @staticmethod
    def Set0Randomly(h5_in_fname, h5_out_fname, directory, lis_del_chn_number, opt=0, probability=0.5):
        """
        set adc as 0 randomly for certain channels.
        @param h5_in_fname: input hdf5 file name
        @param h5_out_fname: output hdf5 file absolute path
        @param directory: input hdf5 directory
        @param lis_del_chn_number: ['board_id:channel_id', ...]
        @param opt: opt=0: act on the channel itself; opt=1: act on neighbour channels
        @param probability: the probability to set as 0
        @return:
        """
        h5_name = CommonUtils.splitext_h5_fname(h5_in_fname)

        n_rows = CommonUtils.get_dataframe_total_lines('%s.h5' % os.path.join(directory, h5_name))

        if os.path.exists(h5_out_fname):
            os.remove(h5_out_fname)

        pageIdx = 0
        count = 0
        for i in Cormis_DataReader_hdf.getCaching(h5_name, 1000000, directory):
            count += i.shape[0]

            df_rand = pd.DataFrame(np.random.random(i.shape[0]))
            for i_num in lis_del_chn_number:
                [i_board, i_chn] = i_num.split(':')
                i_chn_name = 'chn_' + str(i_chn)

                i.reset_index(drop=True, inplace=True)
                ind = i[(df_rand[0] < probability) & (i['boardID'] == eval(i_board))].index
                if opt == 0:
                    i.loc[ind, i_chn_name] = 0
                elif opt == 1:
                    nei_chn_1 = eval(i_chn) - 1
                    nei_chn_name1 = 'chn_' + str(nei_chn_1)
                    if nei_chn_1 >= 0:
                        i.loc[ind, nei_chn_name1] = 0
                    nei_chn_2 = eval(i_chn) + 1
                    nei_chn_name2 = 'chn_' + str(nei_chn_2)
                    if nei_chn_2 < CommonUtils.count_n_channels(i):
                        i.loc[ind, nei_chn_name2] = 0

            i.to_hdf(h5_out_fname, key="data", append=True)

            pageIdx += 1
            print('\rDeleteGivenChannels for %s: %d / %d    ' % (h5_name, count, n_rows), end='')
        print('')

        try:
            df_meas_time = pd.read_hdf(os.path.join(directory, '%s.h5' % h5_name), key="meas_time")
            df_meas_time.to_hdf(h5_out_fname, key="meas_time", append=True)
        except:
            pass
