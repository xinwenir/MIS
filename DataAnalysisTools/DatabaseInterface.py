import os

from DatabaseTool import *
from DatabaseTool import *
import pandas as pd
import time

import logging

logger = logging.getLogger(__name__)


class DatabaseInterface:
    def __init__(self):
        self._cr_online = Cormis_DataReader_online()
        self._cr_offline = Cormis_DataReader_offline()

    def PrintDevicesInfo(self, outf_name):
        deviceInfo = self._cr_online.deviceStatus()
        deviceInfo.to_csv(outf_name)

    def PrintOnlineTablesInfo(self, serial_number, outf_name, serial_name):
        datatableInfos = self._cr_online.deviceDatatablesInfo(serial_number)
        datatableInfos['SerialName'] = serial_name

        datatableInfos = self._add_column_in_table_info(datatableInfos)

        datatableInfos.to_csv(outf_name, index=False)

    def DownloadHdf5(self, df_table_infos, directory, pageCapacity=100000):
        for i in range(df_table_infos.shape[0]):
            table_name = df_table_infos['tableName'].values[i]

            if table_name.index('_tmpData') > 0 and df_table_infos['SerialName'].values[i] != '':
                givenName = table_name.replace(table_name[:table_name.index('_tmpData')],
                                           df_table_infos['SerialName'].values[i])
            else:
                givenName = table_name

            self._cr_online.saveCaching(table_name, pageCapacity, givenName, directory)

            logger.info('Download %s to %s' % (table_name, os.path.join(directory, givenName + ".h5")))

    def PrintOfflineTablesInfo(self, h5_dir, outf_name):
        datatableInfos = self._cr_offline.deviceDatatablesInfo(h5_dir)
        datatableInfos = self._add_column_in_table_info(datatableInfos)

        datatableInfos.to_csv(outf_name, index=False)

    @classmethod
    def _add_column_in_table_info(cls, datatableInfos):
        time_interval = (datatableInfos['stopTime']-datatableInfos['startTime']).astype('int64') \
                        / 10 ** 9 / 60. / 60. / 24
        count_rate = datatableInfos['NRows'] / time_interval / 24
        datatableInfos['timeInterval (day)'] = time_interval.values
        datatableInfos['NRow Rate (per hour)'] = count_rate.values

        return datatableInfos

    # def DrawPie4Tables(self, csv_fname, threshold_low=0, threshold_up=100):
    #     datatableInfos = pd.read_csv(csv_fname)
    #
    #     datatableInfos['timeInterval'].plot.pie()
