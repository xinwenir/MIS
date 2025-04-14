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

# from MyTools.Monitor.Process_monitor import process_monitoring
# process_monitoring(interval=5)
from memory_profiler import profile

from ROOT import TH1F, TFile, TMath, TTree, RDF, TVectorD, TVector3, TGraph, TMultiGraph, TLegend, TGraphErrors

logger = logging.getLogger(__name__)


class CalculatePositionsDF:
    def __init__(self, det):
        self._det_type_id, self._det_type_str = det.get_det_type()

    @classmethod
    def get_df_meas_time(cls, h5_fname):
        # noinspection DuplicatedCode
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

    def get_meas_duration(self, abs_h5_fname, use_det_meas_time=False, det=None):
        if use_det_meas_time:
            meas_time = det.get_meas_time()
        else:
            meas_time = 0
        # get measure time
        if meas_time == 0:
            df_meas_time = self.get_df_meas_time(abs_h5_fname)
            df_meas_time = df_meas_time.drop_duplicates()
            meas_time = df_meas_time['time_delta'].sum()
            meas_time = meas_time.total_seconds()

        return meas_time


class CalculatePositionsDF_PlaneA(CalculatePositionsDF):
    def __init__(self, det):
        super().__init__(det)

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

    @classmethod
    def _get_mtx_data_and_timeTag(cls, df, n_boards, n_bars):
        column_names = ['chn_%d' % i for i in range(n_bars)]
        mtx_data, dic_data = CommonUtils.cvt_dataframe_to_numpy(df[column_names])
        mtx_data = mtx_data.reshape(-1, n_boards, n_bars)
        mtx_timeTag = df['timeTag'].to_numpy().reshape(-1, n_boards)

        return mtx_data, mtx_timeTag

    def _get_energy_matrix_trigger12(self, df, df_det_pro, n_boards, n_bars):
        mtx_data, mtx_timeTag = self._get_mtx_data_and_timeTag(df, n_boards, n_bars)
        lis_mtx_dac = []

        for i_board in range(n_boards):
            tri_bar_up = self._get_tri_up(df_det_pro['TriBarUp'].iloc[i_board])

            mtx_dac = np.zeros((mtx_data.shape[0], 6))  # e1, e2, trigger_count, sign, first_nonzero_index, timeTag

            mtx_board = mtx_data[:, i_board, :]
            timeTag_board = mtx_timeTag[:, i_board]
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

            mtx_dac[:, 4] = first_nonzero_index
            mtx_dac[:, 5] = timeTag_board

            lis_mtx_dac.append(mtx_dac)

        return lis_mtx_dac

    def _get_energy_matrix_trigger34(self, df, df_det_pro: pd.DataFrame, n_boards, n_bars):
        mtx_data, mtx_timeTag = self._get_mtx_data_and_timeTag(df, n_boards, n_bars)
        lis_mtx_dac = []

        for i_board in range(n_boards):
            tri_bar_up = self._get_tri_up(df_det_pro['TriBarUp'].iloc[i_board])

            mtx_dac = np.zeros((mtx_data.shape[0], 6))  # e1, e2, trigger_count, sign, first_nonzero_index, timeTag

            mtx_board = mtx_data[:, i_board, :]
            trigger_count = (mtx_board > 0).sum(axis=1)
            mtx_board = np.insert(mtx_board, n_bars, values=0, axis=1)
            timeTag_board = mtx_timeTag[:, i_board]
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

            mtx_dac[:, 4] = first_nonzero_index
            mtx_dac[:, 5] = timeTag_board

            lis_mtx_dac.append(mtx_dac)

        return lis_mtx_dac

    def _calc_mtx_pos(self, lis_mtx_dac, df_det_pro, n_boards):
        # calculate hit position with energy_matrix
        lis_mtx_pos = []
        for i_board in range(n_boards):
            bar_interval = df_det_pro['barInterval'].iloc[i_board]
            bar_height = df_det_pro['barHeight'].iloc[i_board]
            if_x_board = self._get_xy_board(df_det_pro['xy_board'].iloc[i_board])
            x_origin = self._get_origin(df_det_pro[['Origin_x', 'Origin_y']].iloc[i_board], if_x_board)
            z_origin = df_det_pro['Origin_z'].iloc[i_board]

            mtx_dac = lis_mtx_dac[i_board]

            mtx_pos = np.zeros((mtx_dac.shape[0], 5))  # x,z,trigger_count, first_nonzero_index, timeTag
            energy_ratio = 1. * mtx_dac[:, 1] / (mtx_dac[:, 0] + mtx_dac[:, 1])
            mtx_pos[:, 0] = x_origin + bar_interval / 2 * (mtx_dac[:, 4] + 1 + energy_ratio)
            mtx_pos[:, 1] = z_origin + mtx_dac[:, 3] * (0.5 - energy_ratio) * bar_height
            mtx_pos[:, 2] = mtx_dac[:, 2]

            mtx_pos[:, 3] = mtx_dac[:, 4]
            mtx_pos[:, 4] = mtx_dac[:, 5]

            lis_mtx_pos.append(mtx_pos)

        return lis_mtx_pos

    def _get_raw_dataframe(self, lis_mtx_pos: list, df_det_pro, det_info):
        n_board = len(lis_mtx_pos)

        layer_x = layer_y = 1
        df_data = pd.DataFrame()
        for i_board in range(n_board):
            mtx_pos = lis_mtx_pos[i_board]
            if_x_board = self._get_xy_board(df_det_pro['xy_board'].iloc[i_board])

            if if_x_board:
                df_data['xo%d' % layer_x] = mtx_pos[:, 0]
                df_data['zx%d' % layer_x] = mtx_pos[:, 1]
                df_data['NBar_x%d' % layer_x] = mtx_pos[:, 2].astype('int')
                df_data['board_%d_timeTag' % i_board] = mtx_pos[:, 4].astype('int')
                df_data['ZPlaneX%d' % layer_x] = det_info['Origin_z'][layer_x + layer_y - 2]
                layer_x += 1
            else:
                df_data['yo%d' % layer_y] = mtx_pos[:, 0]
                df_data['zy%d' % layer_y] = mtx_pos[:, 1]
                df_data['NBar_y%d' % layer_y] = mtx_pos[:, 2].astype('int')
                df_data['board_%d_timeTag' % i_board] = mtx_pos[:, 4].astype('int')
                df_data['ZPlaneY%d' % layer_y] = det_info['Origin_z'][layer_x + layer_y - 2]
                layer_y += 1
        return df_data

    @classmethod
    def _get_thetaxy_trk(cls, df_data):
        thetax_trk = TMath.ATan2((df_data.xo1 - df_data.xo2), (df_data.zx1 - df_data.zx2))
        thetay_trk = TMath.ATan2((df_data.yo1 - df_data.yo2), (df_data.zy1 - df_data.zy2))
        return thetax_trk, thetay_trk

    @classmethod
    def _2d_fitting(cls, lis_xo, lis_zx):
        # fit all points
        arr_x = np.array(lis_xo)
        arr_z = np.array(lis_zx)
        slope, intercept, r_value, p_value, str_err = linregress(arr_x, arr_z)
        theta_trk = np.arctan(1 / slope)

        return theta_trk, slope, intercept, r_value, p_value, str_err

    def _get_thetaxy_trk_one_row_numpy(self, arr_xo_zx, n_xo):
        # calculate direction 1d
        lis_xo = arr_xo_zx[:n_xo]
        lis_zx = arr_xo_zx[n_xo:]

        if len(set(lis_xo)) == 1:
            slope = 1e10
            [intercept, r_value, p_value, str_err] = [default_value] * 4
            theta_trk = 0
        elif len(lis_xo) > 1 and len(set(lis_xo)) > 1:
            theta_trk, slope, intercept, r_value, p_value, str_err = self._2d_fitting(lis_xo, lis_zx)
        else:
            [theta_trk, slope, intercept, r_value, p_value, str_err] = [default_value] * 6

        return [theta_trk, slope, intercept, r_value, p_value, str_err]

    @classmethod
    def _get_thetaxy_thetaphi_numpy(cls, arr_data, rotAnglex, rotAngley):
        thetax_trk = arr_data[0]
        thetay_trk = arr_data[1]

        v_det_coord = TVector3(TMath.Tan(thetax_trk), TMath.Tan(thetay_trk), 1)
        v_car_coord = TVector3(v_det_coord)

        if rotAnglex > 1e-10 and rotAngley > 1e-10:
            raise ValueError('rotAnglex>1e-10 && rotAngley>1e-10.')

        thetax = thetay = -10
        if rotAnglex < 1e-10:
            v_car_coord.RotateX(rotAngley)
            thetax = TMath.ATan2(v_car_coord(0), v_car_coord(2))
            thetay = TMath.ATan2(v_car_coord(1), v_car_coord(2))
            # thetay = thetay_trk + rotAngley

        if rotAngley < 1e-10:
            v_car_coord.RotateY(rotAnglex)
            thetay = TMath.ATan2(v_car_coord(1), v_car_coord(2))
            thetax = TMath.ATan2(v_car_coord(0), v_car_coord(2))
            # thetax = thetax_trk + rotAnglex

        theta = v_car_coord.Theta()
        phi = v_car_coord.Phi()
        return [thetax, thetay, theta, phi]

    def _get_whole_dataframe_nplane(self, df_data, rotAnglex, rotAngley, check_angle=0):
        if df_data.empty == True:
            return df_data

        df_data = df_data.sort_index()

        # calculate thetax thetay
        lis_xo_columns = [icol for icol in df_data.columns.values if icol.__contains__('xo')]
        lis_zx_columns = [icol for icol in df_data.columns.values if icol.__contains__('zx')]
        lis_yo_columns = [icol for icol in df_data.columns.values if icol.__contains__('yo')]
        lis_zy_columns = [icol for icol in df_data.columns.values if icol.__contains__('zy')]

        mtx_xo, dic_xo = CommonUtils.cvt_dataframe_to_numpy(df_data[lis_xo_columns])
        mtx_zx, dic_zx = CommonUtils.cvt_dataframe_to_numpy(df_data[lis_zx_columns])
        mtx_x = np.concatenate((mtx_xo, mtx_zx), axis=1)
        mtx_yo, dic_yo = CommonUtils.cvt_dataframe_to_numpy(df_data[lis_yo_columns])
        mtx_zy, dic_zy = CommonUtils.cvt_dataframe_to_numpy(df_data[lis_zy_columns])
        mtx_y = np.concatenate((mtx_yo, mtx_zy), axis=1)

        arr_trk_x = np.apply_along_axis(self._get_thetaxy_trk_one_row_numpy, 1, mtx_x, mtx_xo.shape[1])
        arr_trk_y = np.apply_along_axis(self._get_thetaxy_trk_one_row_numpy, 1, mtx_y, mtx_yo.shape[1])

        df_track_x = pd.DataFrame(arr_trk_x,
                                  columns=['thetax_trk', 'x_slope', 'x_intercept', 'x_r_value', 'x_p_value',
                                           'x_str_err'])
        df_track_y = pd.DataFrame(arr_trk_y,
                                  columns=['thetay_trk', 'y_slope', 'y_intercept', 'y_r_value', 'y_p_value',
                                           'y_str_err'])

        # calculate theta phi
        mtx_thetaxy = np.vstack([arr_trk_x[:,0], arr_trk_y[:,0]])
        arr_thetaphi = np.apply_along_axis(self._get_thetaxy_thetaphi_numpy, 1,
                                           mtx_thetaxy.T,
                                           rotAnglex=rotAnglex, rotAngley=rotAngley)
        df_thetaphi = pd.DataFrame(arr_thetaphi,
                                   columns=["thetax", "thetay", "theta", "phi"])

        df_new = pd.concat([df_data, df_track_x, df_track_y, df_thetaphi],axis=1)

        if not check_angle > 0:
            return df_new

        df_new = df_new[abs(df_new['thetay_trk']) < check_angle]
        df_new = df_new[abs(df_new['thetax_trk']) < check_angle]

        return df_new

    @classmethod
    def _get_Z_XY(cls, df_data):
        if len([icol for icol in df_data.columns.values if icol.__contains__('xo')]) > 0:
            lis_xp_columns = [icol for icol in df_data.columns.values if icol.__contains__('ZPlaneX')]
            ser_xslope = df_data['x_slope']

            for i in lis_xp_columns:
                ser_xo = df_data[i.replace('ZPlaneX', 'xo')]
                ser_zx = df_data[i.replace('ZPlaneX', 'zx')]
                # z-zo = k (x-xo) ==> x=(z-zo)/k+xo
                df_data[i.replace('ZPlaneX', 'PX')] = (df_data[i] - ser_zx) / ser_xslope + ser_xo

        if len([icol for icol in df_data.columns.values if icol.__contains__('yo')]) > 0:
            lis_yp_columns = [icol for icol in df_data.columns.values if icol.__contains__('ZPlaneY')]
            ser_yslope = df_data['y_slope']

            for i in lis_yp_columns:
                ser_yo = df_data[i.replace('ZPlaneY', 'yo')]
                ser_zy = df_data[i.replace('ZPlaneY', 'zy')]
                # z-zo = k (y-yo) ==> y=(z-zo)/k+yo
                df_data[i.replace('ZPlaneY', 'PY')] = (df_data[i] - ser_zy) / ser_yslope + ser_yo

        return df_data

    def _print_calc_matrix(self):
        pass

    @profile(precision=4, stream=open('/home/zxw/data/mis_log/log_mem.txt', 'w+', encoding="utf-8"))
    def CalcHitPoints_NPlane(self, h5_fname, root_fname, det, calc_multi_trigger=True,
                             use_det_meas_time=False):
        h5_file_name = os.path.split(CommonUtils.splitext_h5_fname(h5_fname))[1]
        directory = os.path.split(h5_fname)[0]
        abs_h5_fname = os.path.join(directory, h5_file_name + '.h5')
        n_rows = CommonUtils.get_dataframe_total_lines(abs_h5_fname)

        print('\rCalcHitPoints_NPlane for %s' % h5_file_name, end='')

        detpro_fname = det.get_det_pro_fname()
        lis_rot_angle = det.get_rotAngle()

        df_det_pro = pd.read_csv(detpro_fname)
        rotAnglex = lis_rot_angle[0]
        rotAngley = lis_rot_angle[1]

        df_data = None
        for i in Cormis_DataReader_hdf.getCaching(h5_file_name, 10000, directory):
            df_data = i.copy()
            break
        n_boards = CommonUtils.count_n_boards(df_data)
        n_bars = CommonUtils.count_n_bars(df_data)
        del i,df_data

        meas_time = self.get_meas_duration(abs_h5_fname, use_det_meas_time, det)

        (root_name, file_type) = os.path.splitext(root_fname)

        df_data = None
        pageIndex = 0
        count = 0
        lis_root_splits = []
        for i in Cormis_DataReader_hdf.getCaching(h5_file_name, 100000, directory):
            count += i.shape[0]
            tmp = pd.to_datetime(i['timeTag']).astype('int64') / 10 ** 9
            i['timeTag'] = tmp.astype('int')

            # split out the three or four trigger event
            trigger_12_event, trigger_34_event = self._split_12_34_trigger_events(i)

            # calculate trigger_12 positions and fill into df_data
            lis_mtx_dac = self._get_energy_matrix_trigger12(trigger_12_event, df_det_pro, n_boards, n_bars)
            lis_trigger_12 = self._calc_mtx_pos(lis_mtx_dac, df_det_pro, n_boards)
            df_data = self._get_raw_dataframe(lis_trigger_12, df_det_pro=df_det_pro, det_info=df_det_pro)
            df_data = self._get_whole_dataframe_nplane(df_data, rotAnglex, rotAngley)
            df_data = self._get_Z_XY(df_data)

            del trigger_12_event, lis_mtx_dac, lis_trigger_12

            # calculate trigger_34 positions and fill into df_data
            if calc_multi_trigger:
                lis_mtx_dac = self._get_energy_matrix_trigger34(trigger_34_event, df_det_pro, n_boards, n_bars)
                lis_trigger_34 = self._calc_mtx_pos(lis_mtx_dac, df_det_pro, n_boards)
                t_df_data = self._get_raw_dataframe(lis_trigger_34, df_det_pro=df_det_pro, det_info=df_det_pro)
                check_angle = TMath.ATan2(df_det_pro['barInterval'].iloc[0] / 2, df_det_pro['barHeight'].iloc[0])
                t_df_data = self._get_whole_dataframe_nplane(t_df_data, rotAnglex, rotAngley, check_angle=check_angle)
                df_data = pd.concat((df_data, t_df_data), axis=0)

                del t_df_data, lis_trigger_34, trigger_34_event

            t_root_fname = '%s_%d.root' % (root_name, pageIndex)
            ROOTUtils.cvt_dataframe_to_tree(df_data, 't_tree', t_root_fname)
            lis_root_splits.append(t_root_fname)
            print(r'Save split %s' % t_root_fname, end='')

            del df_data, tmp     # try to avoid large memory
            del i

            pageIndex += 1
            print('\rCalcHitPoints_Plane for %s: %d / %d    ' % (h5_file_name, count, n_rows), end='')
        print('')

        ROOTUtils.merge_root_files(lis_root_splits, root_name)

        # add the measure time information
        measTime = TVectorD(1)
        measTime[0] = meas_time
        fout = TFile(root_fname, 'UPDATE')
        measTime.Write('meas_time')
        fout.Save()
        fout.Close()

    @classmethod
    def _split_12_34_trigger_events(cls, df_data):
        n_bars = CommonUtils.count_n_bars(df_data)
        chn_col_names = ['chn_%d' % i for i in range(n_bars)]
        chn_data = df_data[chn_col_names]
        trigger_count = (chn_data > 1).sum(axis=1)
        trigger_34 = trigger_count > 3
        trigger_34_triggerID = df_data[trigger_34]['triggerID'].values
        i_triggerID_index = df_data.set_index('triggerID')
        trigger_34_event = i_triggerID_index.loc[trigger_34_triggerID]
        trigger_12_event = i_triggerID_index.drop(index=trigger_34_triggerID)

        del i_triggerID_index, trigger_34_triggerID, trigger_34, trigger_count, chn_data

        return trigger_12_event, trigger_34_event
