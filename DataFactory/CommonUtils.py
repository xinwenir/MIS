import os.path
from array import array

import numpy as np
import pandas as pd
import re
from datetime import datetime, timedelta
import tkinter
import tkinter.messagebox

from DatabaseTool.dataLayer.baseCore import h5Data
from FwdTools import points_columns
from DataFactory.Setting import my_tz

from ROOT import TMath, TString, TGraph2D, TVector3
import logging

logger = logging.getLogger(__name__)


def path_isabs(file_path: str):
    if os.path.isabs(file_path):
        return True
    if len(file_path) > 0 and file_path[0] == '~':
        return True
    return False


def phi_cvt_in_minuspi_to_pi(phi: float):
    while phi > TMath.Pi():
        phi -= 2 * TMath.Pi()
    while phi < -1. * TMath.Pi():
        phi += 2 * TMath.Pi()
    return phi


def check_angle_mesh_opt(opt: str):
    if opt.lower().__contains__('phi'):
        mode = 0
    elif opt.lower().__contains__('thetax'):
        mode = 1
    else:
        raise ValueError('opt should be thetaphi or thetaxy; opt = %s' % opt)
    return mode


def calc_divide_and_err(num, err_num, deno, err_deno):
    div = num / deno
    err = TMath.Sqrt(pow(err_num / deno, 2) + pow(num / deno / deno * err_deno, 2))
    return div, err


def dimension_convert(dim: str):
    if dim == '':
        return 10
    tstr = TString(dim)

    if tstr.IsFloat():
        return tstr.Atof()
    else:
        if tstr.Index('m') >= 0:
            string = TString(tstr.Data(), tstr.Index('m'))
            return string.Atof()
        elif tstr.Index('km') >= 0:
            string = TString(tstr.Data(), tstr.Index('km'))
            return string.Atof() * 1000.
        elif tstr.Index('dm') >= 0:
            string = TString(tstr.Data(), tstr.Index('dm'))
            return string.Atof() / 10.
        elif tstr.Index('cm') >= 0:
            string = TString(tstr.Data(), tstr.Index('cm'))
            return string.Atof() / 100.
        else:
            raise ValueError('Do not know what the dimension is: %s. The default unit is meter.' % dim)


def draw_tri_points(df):
    x = array('d',
              [df[points_columns[0]], df[points_columns[3]], df[points_columns[6]]])
    y = array('d',
              [df[points_columns[1]], df[points_columns[4]], df[points_columns[7]]])
    z = array('d',
              [df[points_columns[2]], df[points_columns[5]], df[points_columns[8]]])
    gr2d = TGraph2D(3, x, y, z)
    gr2d.SetMarkerStyle(20)
    gr2d.SetMarkerColor(2)

    return gr2d


def time_convert(time: str):
    if time == '':
        return -1

    try:
        _time = datetime.strptime(time, '%Y-%m-%d %H:%M:%S')
        return _time
    except:
        pass

    tstr = TString(time)

    if tstr.IsFloat():
        return tstr.Atof()
    else:
        if tstr.Index('s') >= 0:
            string = TString(tstr.Data(), tstr.Index('s'))
            return string.Atof()
        elif tstr.Index('min') >= 0:
            string = TString(tstr.Data(), tstr.Index('min'))
            return string.Atof() * 60.
        elif tstr.Index('h') >= 0:
            string = TString(tstr.Data(), tstr.Index('h'))
            return string.Atof() * 60. * 60
        elif tstr.Index('d') >= 0:
            string = TString(tstr.Data(), tstr.Index('d'))
            return string.Atof() * 60. * 60 * 24.
        else:
            return -1


def check_timeTag_with_timeZone(df):
    try:
        a = pd.Timestamp(1, tz='Asia/Shanghai')
        b = df['timeTag'].iloc[0] >= a
        return True
    except:
        return False


def transform_df_start_end_time(df, start_time, end_time):
    start_time = start_time.split('+')[0]
    end_time = end_time.split('+')[0]
    start_second = time_convert(start_time)
    end_second = time_convert(end_time)
    if isinstance(start_second, float):
        start_second = int(start_second + 0.5)
        start_time = df['timeTag'].iloc[0] + timedelta(seconds=start_second)
    elif isinstance(start_second, datetime):
        start_time = start_second
    else:
        start_time = None
    if isinstance(end_second, float):
        end_second = int(end_second + 0.5)
        end_time = df['timeTag'].iloc[0] + timedelta(seconds=end_second)
    elif isinstance(end_second, datetime):
        end_time = end_second
    else:
        end_time = None
    # transform start_time and end_time to pandas datetime and my time zone
    if start_time is not None:
        if check_timeTag_with_timeZone(df):
            start_time = pd.to_datetime(start_time).tz_localize(my_tz)
    if end_time is not None:
        if check_timeTag_with_timeZone(df):
            end_time = pd.to_datetime(end_time).tz_localize(my_tz)
    return start_time, end_time


def convert_timeTag_with_timeZone(timetag):
    if isinstance(timetag, datetime):
        tmp = int(timetag.timestamp())
        return tmp
    if isinstance(timetag, int) or TString(timetag).IsFloat():
        return pd.to_datetime(timetag * 10 ** 9).tz_localize(my_tz)


def count_n_boards(data: pd.DataFrame):
    n = data.drop_duplicates(['boardID'])['boardID'].max() + 1
    return n


def count_n_channels(data: pd.DataFrame):
    lis_column = data.columns.to_list()
    lis_chn = list(filter(lambda x: re.match('chn_\d+', x) != None, lis_column))
    n_channels = len(lis_chn)
    return n_channels


def count_n_bars(data: pd.DataFrame):
    n_board = count_n_boards(data)
    lis_n_bars = []
    for i_board in range(n_board):
        chn_col_names = data.columns[data.columns.str.contains(r'chn_\d+')]
        chn_data = data[data['boardID'] == i_board][chn_col_names]
        sum_chn_data = (chn_data > 0).sum()
        n_bar = sum_chn_data[sum_chn_data > 5].shape[0]
        lis_n_bars.append(n_bar)

    return max(lis_n_bars)


def splitext_h5_fname(h5_fname):
    if h5_fname.__contains__('.h5'):
        h5_name = h5_fname[:h5_fname.index('.h5')]
    else:
        h5_name = h5_fname
    return h5_name


def cvt_dataframe_to_numpy(data):
    _np_df = data.to_numpy()
    _dic_df_columns = dict(zip(data.columns, list(range(0, len(data.columns)))))
    return _np_df, _dic_df_columns


def cvt_numpy_to_dataframe(mtx: np.matrix, dic_df: dict):
    df = pd.DataFrame(mtx, columns=dic_df.keys())
    return df


def replace_csv_channel(csv_in_fname, csv_out_fname, bar_shn_rel_fname):
    df_csv = pd.read_csv(csv_in_fname)
    df_bar_chn_rel = pd.read_csv(bar_shn_rel_fname)

    n_board = df_bar_chn_rel.shape[0]
    n_bar = df_bar_chn_rel.shape[1]
    n_channel = df_csv.shape[1]

    df_csv_out = pd.DataFrame(np.zeros((n_board, n_channel), dtype=int), columns=df_csv.columns)
    for i_board in range(n_board):
        for i_bar in range(n_bar):
            i_chn = df_bar_chn_rel.iloc[i_board, i_bar]
            df_csv_out.iloc[i_board, i_bar] = df_csv.iloc[i_board, i_chn]

    df_csv_out.to_csv(csv_out_fname, index=False)


def get_dataframe_total_lines(tableName: str, key_name: str = 'data'):
    try:
        data = h5Data(tableName)
        return data.index[0]
    except:
        pass
    try:
        try:
            store = pd.HDFStore(tableName)
            n_rows = store.get_storer(key_name).nrows
            store.close()
            return n_rows
        except:
            pass
        data = pd.read_hdf(tableName, key_name, start=-2, stop=-1)
        try:
            n_rows = data['id'].values[0]
            return n_rows
        except:
            pass
        if data.shape[0] > 0:
            return data.index.values[0]
    except:
        print(tableName)
    raise IOError


def pop_up_messages(message, mode='Warning'):
    tkinter.messagebox.showwarning(mode + ': ', message)


def cvt_thetaxy_to_thetaphi(**kwargs):
    thetax = kwargs.get('thetax')
    thetay = kwargs.get('thetay')
    if thetax is None or thetay is None:
        raise ValueError('Utils.thetaxy_to_thetaphi: Please use kwargs here.\nthetax=  thetay=')

    vec = TVector3(TMath.ATan(thetax), TMath.ATan(thetay), 1.)
    dic_thetaphi = {'theta': vec.Theta(),
                    'phi': vec.Phi()}
    return dic_thetaphi


def get_time_list(df_meas_time, time_interval):
    lis_start_time = []
    lis_end_time = []
    _time_interval = timedelta(seconds=time_interval)
    for i in range(df_meas_time.shape[0]):
        df_one = df_meas_time.iloc[i]
        _start_time = pd.to_datetime(df_one['start_time'])
        _time_delta = pd.to_timedelta(df_one['time_delta'])
        n_slice = int(_time_delta / _time_interval) + 1

        for i_slice in range(n_slice):
            _st = i_slice * _time_interval + _start_time
            if i_slice < (n_slice - 1):
                _et = (i_slice + 1) * _time_interval + _start_time
            else:
                _et = _start_time + _time_delta
            if (_et - _st) > timedelta(minutes=5):
                lis_start_time.append(_st)
                lis_end_time.append(_et)

    if not isinstance(lis_start_time[0], pd._libs.tslibs.timestamps.Timestamp):
        lis_start_time = [convert_timeTag_with_timeZone(i) for i in lis_start_time]
        lis_end_time = [convert_timeTag_with_timeZone(i) for i in lis_end_time]

    return lis_start_time, lis_end_time


def cut_dataframe_with_selection(df: pd.DataFrame, selection: str):
    try:
        df = df[eval(selection)]
    except Exception as e:
        print(e)
        raise Exception('validate the selection %s.')
    return df


def search_files_with_re(path, pattern):
    lis_f = os.listdir(path)
    re_pattern = re.compile(pattern)

    lis_f = [f for f in lis_f if re.search(re_pattern, f)]
    print("Files founded by %s:\n" % pattern, lis_f)

    lis_f = [os.path.join(path, i) for i in lis_f]

    return lis_f


def find_same_prefix(lis_strs):
    prefix = ''
    for x in zip(*lis_strs):  # 每个单词依次取一个字符组合
        if len(set(x)) == 1:  # 若去重后长度为1则为共有前缀字符
            prefix = prefix + x[0]  # 将共有字符拼接到字符串上
        else:  # 若出现不同字符结束循环
            break
    if prefix:  # 若字符串非空
        return prefix  # 返回共有前缀
    else:
        return 'NOT FOUND'
