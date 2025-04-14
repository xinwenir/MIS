# -*- coding: utf-8 -*-
# @Author  : jimbook
# @email   : tianh2021@lzu.edu.cn
# @Time    : 2022/9/15
# @Note    : an offline data reader from local data directory
# @modify  : None
from DatabaseTool import *
from DatabaseTool.Interface import *
from DatabaseTool.onlineReader import Cormis_CachingHDFIterator
from sqlalchemy import create_engine
import re
import os
import time
from .dataLayer.baseCore import h5Data
from .dataLayer import _Index
from DataFactory.Setting import my_tz
from DataFactory import CommonUtils


class Cormis_DataReader_offline(Cormis_DataReader_Interface):
    def dataTableList(self, rootPath):
        tableList = []
        for root, dirs, files in os.walk(rootPath):
            for name in files:
                if os.path.splitext(name)[1] == ".h5":
                    tableList.append(os.path.join(root, name))
        return tableList

    def getStartTime(self, tableName: str):
        from DataFactory.Setting import my_tz
        try:
            data = h5Data(tableName)
            first = data.timeTag[0]
            tm = pd.Timestamp(first, unit="s", tz=my_tz)
            data.close()
            return tm
        except:
            pass
        try:
            data = pd.read_hdf(tableName, "data", start=0, stop=1)
            time_int = data["timeTag"].values[0].astype('datetime64[s]').astype('int')
            return pd.Timestamp(time_int, unit="s", tz=my_tz)
        except:
            pass
        raise IOError

    def getStopTime(self, tableName: str):
        from DataFactory.Setting import my_tz
        try:
            data = h5Data(tableName)
            idx = data.index[0]
            last = data.timeTag[idx - 1]
            tm = pd.Timestamp(last, unit="s", tz=my_tz)
            data.close()
            return tm
        except:
            pass
        try:
            data = pd.read_hdf(tableName, "data", start=-2, stop=-1)
            time_int = data["timeTag"].values[0].astype('datetime64[s]').astype('int')
            return pd.Timestamp(time_int, unit="s", tz=my_tz)
        except:
            pass
        raise IOError

    def getCount(self, tableName: str):
        try:
            data = h5Data(tableName)
            return data.index[0]
        except:
            pass
        try:
            data = pd.read_hdf(tableName, "data", start=-2, stop=-1)
            if data.shape[0] > 0:
                return data.index.values[0]
            else:
                return 0
        except:
            pass
        raise IOError

    def deviceDatatablesInfo(self, rootPath: str = "./"):
        tableNames = []
        startTimes = []
        stopTimes = []
        counts = []
        for name in self.dataTableList(rootPath):
            (name_dir, name_file) = os.path.split(name)
            if os.path.getsize(name) < 1025:
                continue

            print('\rGet deviceDataTablesInfo %s                ' % name_file, end='')
            count = self.getCount(name)
            if count > 0:
                tableNames.append(name)
                startTimes.append(self.getStartTime(name))
                stopTimes.append(self.getStopTime(name))
                counts.append(count)
        print('')
        return pd.DataFrame({"tableName": tableNames, "startTime": startTimes, "stopTime": stopTimes, "NRows": counts})

    # 获取对应数据表的数据
    def deviceData(self, tableName: str, startID: int = None, stopID: int = None) -> pd.DataFrame:
        try:
            result = self.getH5UsingOldAPI(tableName, startID, stopID)
        except:
            result = pd.read_hdf(tableName, key="data", start=startID, stop=stopID)
            if len(result) == 0:
                return result

            if CommonUtils.check_timeTag_with_timeZone(result):
                return result

            def func(timeInt: int):
                tm = pd.Timestamp(timeInt / 10 ** 9, unit='s', tz=my_tz)
                return tm

            result['timeTag'] = result['timeTag'].astype('int64').apply(func)
        return result

    # 以迭代器形式返回数据表的数据，每次迭代的数据量由用户指定
    def deviceDataSplitIter(self, tableName: str, pageCapacity: int):
        try:
            Iterator = Cormis_DataIterator_Interface(self, tableName, pageCapacity)
        except:
            Iterator = Cormis_DataIterator_Interface(self, tableName, pageCapacity)
        return Iterator

    def getH5UsingOldAPI(self, tableName: str, start: int = None, stop: int = None):
        from DataFactory.Setting import my_tz

        data = h5Data(tableName)
        result = pd.DataFrame(columns=_Index).astype(int)
        if start is None:
            start = 0
        if stop is None:
            stop = data.index[0]
        count = stop - start
        if count <= 0:
            return result
        for i in range(data.index.shape[0]):
            _s = data.getData(i, startIndex=start)
            _s['triggerID'] = _s['triggerID'] + i * 65535
            # result = result.append(_s, ignore_index=True)
            result = pd.concat([result, _s], ignore_index=True)
            if result.shape[0] >= count:
                break
        result = result.iloc[:count]
        if len(result) == 0:
            return result

        def func(timeInt: int):
            tm = pd.Timestamp(timeInt, unit='s', tz=my_tz)
            return tm

        result['timeTag'] = result['timeTag'].apply(func)
        return result

    @staticmethod
    def getCaching(tableName: str, pageCapacity: int = 100000, rootPath: str = ""):
        cr = Cormis_DataReader_offline()
        path = os.path.join(rootPath, tableName + ".h5")
        return cr.deviceDataSplitIter(path, pageCapacity)

#
