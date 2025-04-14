from DatabaseTool import *
from sqlalchemy import create_engine
import re
import os
import time
import shutil

import logging
logger = logging.getLogger(__name__)

class Cormis_DataReader():
    def __init__(self, useCaching: bool = False, mergeCaching: bool = False):
        self.DB = create_engine(f'mysql+pymysql://{UserName}:{PassWord}@{Host}:{Port}/{ProjectName}')
        self.useCaching = useCaching
        self.mergeCaching = mergeCaching

    def exeCommand(self, command: str):
        return pd.read_sql(command, self.DB)

    def deviceStatus(self) -> pd.DataFrame:
        return self.exeCommand("SELECT * FROM device;")

    def deviceDatatableList(self, serialNumber: str):
        result = []
        allTablesList = self.exeCommand("SHOW TABLES;")
        for i in allTablesList.iloc[:, 0]:
            d = re.match(serialNumber + "_.+", i)
            if d is not None:
                result.append(i)
        return result

    def getStartTime(self, tableName: str):
        cmd = "SELECT timeTag from {} LIMIT 0, 1;".format(tableName)
        line = self.exeCommand(cmd)
        return line.iloc[0, 0]

    def getStopTime(self, tableName: str):
        cmd = "SELECT timeTag from {0} where id = (SELECT max(id) FROM {0});".format(tableName)
        line = self.exeCommand(cmd)
        return line.iloc[0, 0]

    def getCount(self, tableName: str):
        cmd = "SELECT count(id) FROM {};".format(tableName)
        line = self.exeCommand(cmd)
        return line.iloc[0, 0]

    def deviceDatatablesInfo(self, serialNumber: str):
        tableNames = []
        startTimes = []
        stopTimes = []
        counts = []
        allTablesList = self.exeCommand("SHOW TABLES;")
        for i in allTablesList.iloc[:, 0]:
            d = re.match(serialNumber + "_.+", i)
            if d is not None:
                count = self.getCount(i)
                if count > 0:
                    tableNames.append(i)
                    startTimes.append(self.getStartTime(i))
                    stopTimes.append(self.getStopTime(i))
                    counts.append(count)
        return pd.DataFrame({"tableName": tableNames, "startTime": startTimes, "stopTime": stopTimes, "count": counts})

    def deviceData(self, tableName: str, stopID: int = None, startID: int = None):
        s = startID is not None
        e = stopID is not None
        scmd = None
        if s and e:
            scmd = "LIMIT {}, {}".format(startID, stopID - startID)  # offset, count
        elif e:
            scmd = "LIMIT {}".format(stopID)
        elif s:
            return None
        if scmd is None:
            cmd = "SELECT * FROM {};".format(tableName)
        else:
            cmd = "SELECT * FROM {} {};".format(tableName, scmd)
        return self.exeCommand(cmd)

    def deviceDataSplitIter(self, tableName: str, pageCapacity: int = 100000):
        """
        拿到对应数据表的数据分片迭代器，数据分片会根据triggerID调整每片的行数，来防止将一个event切开
        :param tableName: 数据表名
        :param pageCapacity: 每个分片期望的行数，最小为100
        :return:
        """
        iterator = Cormis_DataIterator(self, tableName, pageCapacity)
        return iterator

    @staticmethod
    def getCaching(tableName: str):
        return Cormis_CachingIterator(tableName)

    def saveCaching(self, tableName: str, pageCapacity: int = 100000):
        cachingPath = os.path.join(os.path.split(__file__)[0], "dataCaching")
        dirPath = os.path.join(cachingPath, tableName)
        allStartTime = time.time()
        allCount = self.getCount(tableName)
        count = 0
        thisTime = time.time()
        if os.path.exists(dirPath):
            shutil.rmtree(dirPath)
        pageIdx = 0
        os.mkdir(dirPath)
        filePath = os.path.join(dirPath, tableName + "_{}.txt")
        for i in self.deviceDataSplitIter(tableName, pageCapacity):
            i.to_csv(filePath.format(pageIdx))
            thisCount = i.shape[0]
            count += thisCount
            pageIdx += 1
            print("\rcount: {}/{} + {}, time use:{}; all time use: {}".format(count, allCount, thisCount,
                                                                              time.time() - thisTime,
                                                                              time.time() - allStartTime), end="")
        print()

    @classmethod
    def isCached(cls, cachingPath: str = None):
        if cachingPath is None:
            cachingPath = os.path.join(os.path.split(__file__)[0], "dataCaching")
        result = []
        for root, dirs, files in os.walk(cachingPath):
            for file in files:
                print(file)
                result.append(file)
        return result


class Cormis_DataIterator():
    def __init__(self, parent: Cormis_DataReader, tableName: str, pageCapacity: int = 10000):
        self.parent = parent
        self.tableName = tableName
        self.pageCapacity = pageCapacity if pageCapacity > 100 else 100
        self.count = parent.getCount(tableName)

    def getTriggerIDOffset(self, startID: int):
        cmd = "SELECT triggerID FROM {} LIMIT {}, 1;"
        d = self.parent.exeCommand(cmd.format(self.tableName, startID))
        if d.shape[0] == 0:
            return 0
        tID_ref = d.iloc[0, 0]
        for i in range(1, 30):
            d = self.parent.exeCommand(cmd.format(self.tableName, startID + i))
            tID = d.iloc[0, 0]
            if d.shape[0] == 0:
                return i
            if tID != tID_ref:
                return i
        return -1

    def __iter__(self):
        self.endTag = False
        self.endID = 0
        self.offset = 0
        return self

    def __next__(self):
        if self.endTag:
            raise StopIteration
        else:
            nextEND = self.endID + self.pageCapacity
            offset = self.getTriggerIDOffset(nextEND)
            d = self.parent.deviceData(self.tableName, nextEND + offset, self.endID + self.offset)
            self.offset = offset
            self.endID = nextEND
            if self.offset == -1:
                self.endTag = True
            elif nextEND >= self.count:
                self.endTag = True
            return d


class Cormis_DataReader_hdf(Cormis_DataReader):
    def __init__(self, useCaching: bool = False, mergeCaching: bool = False):
        super(Cormis_DataReader_hdf, self).__init__(useCaching, mergeCaching)

    def saveCaching(self, tableName: str, pageCapacity: int = 100000, saveAsName: str = None, rootDir: str = None):
        if rootDir is not None:
            cachingPath = rootDir
        else:
            cachingPath = os.path.join(os.path.split(__file__)[0], "dataCaching")
        # dirPath = os.path.join(cachingPath, tableName)
        dirPath = cachingPath
        allStartTime = time.time()
        allCount = self.getCount(tableName)

        # if os.path.exists(dirPath):
        #     shutil.rmtree(dirPath)
        # os.mkdir(dirPath)
        thisTime = time.time()

        pageIdx = 0
        if saveAsName is None:
            filePath = os.path.join(dirPath, tableName + ".h5")
        else:
            filePath = os.path.join(dirPath, saveAsName + ".h5")

        if os.path.exists(filePath):
            os.remove(filePath)

        count = 0
        for i in self.deviceDataSplitIter(tableName, pageCapacity):
            i.to_hdf(filePath.format(pageIdx), key="data", append=True)
            # print('print to', filePath.format(pageIdx))
            thisCount = i.shape[0]
            count += thisCount
            pageIdx += 1
            print("\rcount: {}/{} + {}, time use:{}; all time use: {}".format(count, allCount, thisCount,
                                                                              time.time() - thisTime,
                                                                              time.time() - allStartTime), end="")
            thisTime = time.time()

    @staticmethod
    def getCaching(tableName: str, pageCapacity: int = 100000, rootDir: str = None):
        return Cormis_CachingHDFIterator(tableName, pageCapacity, rootDir)


class Cormis_CachingIterator():
    def __init__(self, tableName: str):
        self.tableName = tableName
        cachingPath = os.path.join(os.path.split(__file__)[0], "dataCaching")
        # self.path = os.path.join(cachingPath,tableName)
        self.path = cachingPath
        self.exit = os.path.exists(self.path)

    def __iter__(self):
        self.subPath = os.path.join(self.path, self.tableName + "_{}.txt")
        self.pageIdx = 0
        return self

    def __next__(self):
        if not self.exit:
            raise StopIteration
        while True:
            p = self.subPath.format(self.pageIdx)
            self.pageIdx += 1
            if os.path.exists(p):
                return pd.read_csv(p, index_col=0)
            else:
                raise StopIteration


class Cormis_CachingHDFIterator(Cormis_CachingIterator):
    def __init__(self, tableName: str, pageCapacity, cachingPath: str = None):
        super(Cormis_CachingHDFIterator, self).__init__(tableName)
        if cachingPath is not None:
            # self.path = os.path.join(cachingPath, tableName)
            self.path = cachingPath
            self.exit = os.path.exists(self.path)
        self.pageCapacity = pageCapacity

    def __iter__(self):
        self.subPath = os.path.join(self.path, self.tableName + ".h5")
        if os.path.exists(self.subPath):
            self.countIdx = 0
            return self
        else:
            logger.error('ERROR: %s dose not exist.' % self.subPath)
            return super(Cormis_CachingHDFIterator, self).__iter__()

    def __next__(self):
        if not self.exit:
            raise StopIteration
        while True:
            data: pd.DataFrame = pd.read_hdf(self.subPath, key="data", mode="r", start=self.countIdx,
                                             stop=self.countIdx + self.pageCapacity)
            if data.shape[0] == 0:
                raise StopIteration
            elif data.shape[0] < self.pageCapacity:
                self.countIdx += self.pageCapacity
                return data.reset_index().iloc[:, 1:]
            else:
                idx = data.shape[0] - 1
                endTID = data['triggerID'].iloc[idx]
                for i in range(100):
                    thisTID = data['triggerID'].iloc[idx - i]
                    if thisTID != endTID:
                        self.countIdx += idx - i + 1
                        return data.iloc[:idx - i + 1]
