# -*- coding: utf-8 -*-
# @Author  : jimbook
# @email   : tianh2021@lzu.edu.cn
# @Time    : 2022/9/15
# @Note    : an online data reader from MySQL DATABASE
# @modify  : None
from DatabaseTool import *
from DatabaseTool.Interface import *
from sqlalchemy import create_engine
import re
import os
import time


# 从在先的MySQL服务器中获取数据的类，MySQL服务器需要在校园网内才能访问
class Cormis_DataReader_online(Cormis_DataReader_Interface):
    def __init__(self):
        '''
        :param useCaching: 是否使用本地缓存，ture则会先检索本地缓存，如果没有则从数据访问数据
        '''
        self.DB = create_engine(f'mysql+pymysql://{UserName}:{PassWord}@{Host}:{Port}/{ProjectName}')

    # 在MySQL上执行对应命令
    def exeCommand(self, command: str):
        return pd.read_sql(command, self.DB)

    # 获取当前探测器连接状态
    def deviceStatus(self) -> pd.DataFrame:
        return self.exeCommand("SELECT * FROM device;")

    # 获取对应探测器的数据表名列表
    def deviceDatatableList(self, serialNumber: str):
        result = []
        allTablesList = self.exeCommand("SHOW TABLES;")
        for i in allTablesList.iloc[:, 0]:
            d = re.match(serialNumber + "_.+", i)
            if (d is not None):
                result.append(i)
        return result

    # 获取对应数据表的开始测量时间
    def getStartTime(self, tableName: str):
        cmd = "SELECT timeTag from {} LIMIT 0, 1;".format(tableName)
        line = self.exeCommand(cmd)
        return line.iloc[0, 0]

    # 获取对应数据表的结束测量时间
    def getStopTime(self, tableName: str):
        cmd = "SELECT timeTag from {0} where id = (SELECT max(id) FROM {0});".format(tableName)
        line = self.exeCommand(cmd)
        return line.iloc[0, 0]

    # 获取对应数据表的行数（数据量）
    def getCount(self, tableName: str):
        cmd = "SELECT count(id) FROM {};".format(tableName)
        line = self.exeCommand(cmd)
        return line.iloc[0, 0]

    # 获取对应数据集合内的数据信息，返回包括数据表名、开始时间、结束时间、数据量
    def deviceDatatablesInfo(self, serialNumber: str = None):
        '''
        :param serialNumber: 如果给出了对应的探测器序列号则只返回对应探测器的数据表，若给出None，则将所有数据数据表均返回
        :return:数据表名、开始时间、结束时间、数据量
        '''
        tableNames = []
        startTimes = []
        stopTimes = []
        counts = []
        allTablesList = self.exeCommand("SHOW TABLES;")
        for i in allTablesList.iloc[:, 0]:
            if serialNumber is None:
                d = re.match(r"\d+?_.+", i)
            else:
                d = re.match(serialNumber + "_.+", i)
            if (d is not None):
                count = self.getCount(i)
                if count > 0:
                    tableNames.append(i)
                    startTimes.append(self.getStartTime(i))
                    stopTimes.append(self.getStopTime(i))
                    counts.append(count)
        return pd.DataFrame({"tableName": tableNames, "startTime": startTimes, "stopTime": stopTimes, "NRows": counts})

    # 获取对应数据表的数据(注意这里更改了start 和 stop 的顺序）
    def deviceData(self, tableName: str, startID: int = None, stopID: int = None, ):
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

    # 以迭代器形式返回数据表的数据，每次迭代的数据量由用户指定
    def deviceDataSplitIter(self, tableName: str, pageCapacity: int = 100000):
        """
        拿到对应数据表的数据分片迭代器，数据分片会根据triggerID调整每片的行数，来防止将一个event切开
        :param tableName: 数据表名
        :param pageCapacity: 每个分片期望的行数，最小为100
        :return:
        """
        iterator = Cormis_DataIterator_Interface(self, tableName, pageCapacity)
        return iterator

    # 从数据库保存缓存数据到本地
    def saveCaching(self, tableName: str, pageCapacity: int = 100000, saveAsName: str = None, rootDir: str = None,
                    skipIfExist=False):
        '''
        :param tableName: 需要下载的数据表名
        :param pageCapacity: 下载时，按对应的分片不断下载
        :param saveAsName: 以特定的数据名下载，默认会以数据表名存储
        :param rootDir: 下载到指定文件夹下，默认为模块中的dataCaching文件夹
        :param skipIfExist: 如果监测到存在同名文件，是否跳过下载，默认为否，表示将会覆盖本地数据从新下载
        :return:
        '''
        if rootDir is not None:
            cachingPath = rootDir
        else:
            cachingPath = os.path.join(os.path.split(__file__)[0], "dataCaching")
        dirPath = cachingPath
        allStartTime = time.time()
        allCount = self.getCount(tableName)
        thisTime = time.time()
        pageIdx = 0
        if not os.path.exists(dirPath):
            os.mkdir(dirPath)
        if saveAsName is None:
            filePath = os.path.join(dirPath, tableName + ".h5")
        else:
            filePath = os.path.join(dirPath, saveAsName + ".h5")
        if os.path.exists(filePath):
            if skipIfExist:
                print("{} has cached, skip to caching it.".format(filePath))
            os.remove(filePath)
        count = 0
        for i in self.deviceDataSplitIter(tableName, pageCapacity):
            i = i.set_index('id')
            i.to_hdf(filePath.format(pageIdx), key="data", append=True)
            thisCount = i.shape[0]
            count += thisCount
            pageIdx += 1
            print("\rcount: {}/{} + {}, time use:{}; all time use: {}".format(count, allCount, thisCount,
                                                                              time.time() - thisTime,
                                                                              time.time() - allStartTime), end="")
            thisTime = time.time()

    @classmethod
    # changed by green.Liu
    def isCached(cls, cachingPath: str = None):
        if cachingPath is None:
            cachingPath = os.path.join(os.path.split(__file__)[0], "dataCaching")
        result = []
        for root, dirs, files in os.walk(cachingPath):
            for file in files:
                print(file)
                result.append(file)
        return result

    @staticmethod
    # 获取缓存数据的切片迭代器
    def getCaching(tableName: str, pageCapacity: int = 100000, rootDir: str = None):
        return Cormis_CachingHDFIterator(None, tableName, pageCapacity, rootDir)

    # 缓存对应探测器数据的所有数据表
    def saveAllAboutDeivce(self, serialNumber: str, pageCapacity: int = 100000, rootDir: str = None, skipIfExist=False):
        tableList = self.deviceDatatableList(serialNumber)
        for name in tableList:
            self.saveCaching(name, pageCapacity=pageCapacity, rootDir=rootDir, skipIfExist=skipIfExist)


# 获取缓存数据切片的迭代器
class Cormis_CachingHDFIterator(Cormis_DataIterator_Interface):
    def __init__(self, parent: Cormis_DataReader_Interface, tableName: str, pageCapacity: int, cachingPath: str = None):
        super(Cormis_CachingHDFIterator, self).__init__(parent, tableName, pageCapacity)
        if cachingPath is not None:
            self.path = cachingPath
        else:
            self.path = os.path.join(os.path.split(__file__)[0], "dataCaching")
        self.exit = os.path.exists(self.path)
        self.pageCapacity = pageCapacity

    def getDataSPlited(self, start: int, stop: int) -> pd.DataFrame:
        return pd.read_hdf(self.subPath, key="data", mode="r", start=start, stop=stop)

    def __iter__(self):
        self.subPath = os.path.join(self.path, self.tableName + ".h5")
        if (os.path.exists(self.subPath)):
            return super(Cormis_CachingHDFIterator, self).__iter__()
        else:
            raise IOError
