import os, logging

import pandas as pd
import numpy as np

from DatabaseTool.dataLayer.baseCore import h5Data
from DatabaseTool import *
from DataFactory import CommonUtils, ROOTUtils

logger = logging.getLogger(__name__)


class DataFormatTransformer:
    def __init__(self):
        self._cr_offline = Cormis_DataReader_offline()

    def _reshape_data_by_numpy(self, df, n_boards, lis_col_names):
        # reshape n boards to 1 event
        if len(df) % n_boards != 0:
            raise ValueError('Convert2OneEventPerRow requires the data to be complete with %d boards.\n'
                             'Try to use FiltrateAllHit first' % n_boards)

        np_df, dic_col = CommonUtils.cvt_dataframe_to_numpy(df)
        np_df = np_df.reshape(int(np_df.shape[0] / 4), np_df.shape[1] * 4)

        df_new = pd.DataFrame(np_df, columns=lis_col_names)

        return df_new

    def Convert2OneEventPerRow(self, h5_fname, root_fname, tree_name='tr_data'):
        h5_file_name = os.path.split(CommonUtils.splitext_h5_fname(h5_fname))[1]
        directory = os.path.split(h5_fname)[0]
        abs_h5_fname = os.path.join(directory, h5_file_name + '.h5')
        n_rows = CommonUtils.get_dataframe_total_lines(abs_h5_fname)

        out_h5_fname = root_fname.replace('.root', '.h5')

        df_data = None
        for i in Cormis_DataReader_hdf.getCaching(h5_file_name, 100000, directory):
            df_data = i.copy()
            break
        n_boards = CommonUtils.count_n_boards(df_data)
        n_channels = CommonUtils.count_n_channels(df_data)
        n_bars = CommonUtils.count_n_bars(df_data)

        # initialize a list of column names for new dataframe
        lis_col_names = []
        lis_colID = []
        for iboard in range(n_boards):
            min_colID = len(lis_col_names)
            for i in df_data.columns:
                lis_col_names.append('board_%d_%s' % (iboard, i))
            max_colID = len(lis_col_names)
            lis_colID.append([min_colID, max_colID])

        count = 0
        pageIndex = 0
        lis_root_splits = []
        for i in Cormis_DataReader_hdf.getCaching(h5_file_name, 100000, directory):
            count += i.shape[0]

            tmp = pd.to_datetime(i['timeTag']).astype('int64')/ 10**9
            i['timeTag'] = tmp.astype('int')

            df_recon = self._reshape_data_by_numpy(i, n_boards, lis_col_names)

            t_root_fname = '%s_%d.root' % (root_fname, pageIndex)
            ROOTUtils.cvt_dataframe_to_tree(df_recon, tree_name, t_root_fname)
            lis_root_splits.append(t_root_fname)

            df_recon.to_hdf(out_h5_fname, key="data", append=True)

            pageIndex += 1
            print('\rConvertHdf5ToTree for %s: %d / %d    ' % (h5_file_name, count, n_rows), end='')

        ROOTUtils.merge_root_files(lis_root_splits, root_fname + '.root')
