'''
包含一些基础使用的数据对象
'''
import numpy as np
import pandas as pd
from datetime import datetime
import h5py
from DataAnalysisTools import _Index
import time

__all__ = ['h5Data', 'LUMISData']


# 操作h5文件
class h5Data(object):
    '''├│└
    =========FORMAT OF H5 FILE========
    root
     ├info(group)
     │ ├time(dataset)
     │ │ ├start time
     │ │ ├stop time
     │ │ └total time
     │ └device_status(dataset)
     │   ├trigger mode
     │   ├board 0 threshold
     │   ├~~~~~~~
     │   ├board 7 threshold
     │   ├board 0 bias voltage
     │   ├~~~~~
     │   └board 7 bias voltage
     └dataGroup(group)
        ├data(dataSet)
        │ └data
        ├timeTag(dataSet)
        │ └timeTag
        └index(dataSet)
          └index
    '''

    def __init__(self, path: str, mode: str = 'r'):
        '''
        :param path: the path of h5 file
        :param mode: "w" or "r"
        内部参数：
            file:h5文件
            readMode:是否是读模式
            SetsIndex:数据集的索引
            dataIndex:数据集的数据索引
        ==========dataset==========
            timeInfo:'time'
            statusInfo:'device_status'

            data:'data',数据
            timeTag:'timeTag', 时间标签
            index:'index'-----数据存储结构：为线性数组，第一个数据(索引为0)为当前数据总行数，
                        之后的每一个都是triggerID置零时的索引(第一个triggerID为0时的索引)
                        index的shape和重置次数对应
        '''
        if mode[0] == 'w':
            self.file = h5py.File(path, mode=mode, libver='latest')
            self.readMode = False
            # 写模式单独参数
            self.dataIndex = 0  # 当前写到第几行
            self.setsIndex = 0  # 当前存在几次triggerID置零
            # 存储信息
            infoGroup = self.file.create_group('info')  # 信息组
            self.timeInfo = infoGroup.create_dataset('time', shape=(3,), maxshape=(None,),
                                                     dtype='float64')  # 开始时间、结束时间、各分片运行时间
            self.statusInfo = infoGroup.create_dataset('device_status', shape=(17,), dtype='uint16')  # 各个板的状态信息
            # 数据信息
            dataGroup = self.file.create_group('dataGroup')  # 数据组
            self.data = dataGroup.create_dataset('data', shape=(5000, 39), maxshape=(None, 39), dtype='uint16')
            self.timeTag = dataGroup.create_dataset('timeTag', shape=(5000,), maxshape=(None,), dtype='uint64')
            self.index = dataGroup.create_dataset('index', shape=(1,), maxshape=(None,), dtype='int64')
            # 设置为single writer multiply reader 模式
            self.file.swmr_mode = True
        elif mode[0] == 'r':
            self.file = h5py.File(path, mode=mode, libver='latest', swmr=True)
            self.readMode = True
            # 存储信息
            infoGroup = self.file['info']
            self.timeInfo = infoGroup['time']
            self.statusInfo = infoGroup['device_status']
            # 数据信息
            dataGroup = self.file['dataGroup']
            self.data = dataGroup['data']
            self.timeTag = dataGroup['timeTag']
            self.index = dataGroup['index']

    def get_n_splits(self):
        return self.index.shape[0]

    # ===========信息==========
    # 添加开始时间
    def startTime(self) -> datetime:
        '''
        "w":add start time (now time) to h5 file,and return it
        "r":just return start time (from h5)
        :return: start time
        '''
        if self.readMode:
            start = datetime.fromtimestamp(self.timeInfo[0])
            return start
        else:
            self.start = datetime.now()
            self.timeInfo[0] = self.start.timestamp()
            return self.start

    # 添加结束时间和总运行时间,并关闭h5文件
    def stopTime(self) -> datetime:
        '''
        "w":add stop time (now time) and total running time to h5 file,and return it
        "r":just return stop time (from h5)
        :return: stop time
        '''
        if self.readMode:
            stop = datetime.fromtimestamp(self.timeInfo[1])
            return stop
        else:
            stop = datetime.now()
            t = stop.timestamp() - self.start.timestamp()
            self.timeInfo[1] = stop.timestamp()
            self.timeInfo[2] = t
            self.close()
            return stop

    # 返回总运行时间
    def totalTime(self) -> float:
        '''
        "w":raise Exception
        "r": get total running time (from h5)
        :return: total running time (unit:s)
        '''
        return self.timeInfo[1] - self.timeInfo[0]

    # 录入设备状态信息
    def putDeviceStatus(self, status: dict):
        '''
        "w": write device configuration state to h5 file
        "r": raise Exception
        :param status:dict:key:boardID;value:tuple(threshold, biasVoltage, triggerMode)
        :return:
        '''
        # print(status)
        triggerMode = status[0][2]
        self.statusInfo[0] = triggerMode
        for i in range(8):
            t = status.get(i, (0, 0, 0))
            self.statusInfo[i + 1] = t[0]
            self.statusInfo[i + 9] = t[1]

    # 获取h5文件中存储的设备状态信息
    def getDeviceStatus(self) -> dict:
        '''
        :return: dict:key:boardID;value:list(threshold, bias voltage, trigger mode)
        '''
        status = {}
        triggerMode = self.statusInfo[0]
        for i in range(8):
            threshold = self.statusInfo[i + 1]
            biasVoltage = self.statusInfo[i + 9]
            if threshold == 0 and biasVoltage == 0:
                break
            else:
                status[i] = [threshold, biasVoltage, triggerMode]
        return status

    # ===============数据=================
    # 将数据添加到h5数据集中
    def addToDataSet(self, value: np.array, timeTag=None):
        '''
        "w": add data to h5 file
        "r": raise Exception
        :param value: data,type:np.array,dtype:'uint16',shape:(n,39)
                    column 0-38 will be pushed into group(data), column 39 will be pushed into group(timeTag)
        :param NewSets: if it set as True,h5 file will create a new data set to contain data
        :return:
        '''
        # 如果数据集的长度不够，则扩展5000行
        if self.dataIndex + value.shape[0] <= self.data.shape[0]:
            self.data.resize(self.data.shape[0] + 5000, axis=0)
            self.timeTag.resize(self.data.shape[0] + 5000, axis=0)
        # 数据存储
        self.data.write_direct(value, source_sel=np.s_[0:value.shape[0]],
                               dest_sel=np.s_[self.dataIndex:value.shape[0] + self.dataIndex])
        if timeTag is not None:
            self.timeTag.write_direct(timeTag, source_sel=np.s_[0:value.shape[0]],
                                      dest_sel=np.s_[self.dataIndex:value.shape[0] + self.dataIndex])
        self.dataIndex += value.shape[0]
        self.index[0] = self.dataIndex

    # 读取h5文件中的数据
    def getData(self, index: int = -1, startIndex: int = 0, dtype='int64') -> pd.DataFrame:
        '''
            read data in h5 file
            :param index: When the index greater than or equal to 0, the corresponding set of data will be returned.
            when the index equal to -1,all the data will be returned.
            when the index equal to -2,all the data will be returned and the triggerID in returned data has been accumulated.
            :param startIndex: data will be returned from the startIndex
            :return:
        '''
        if index + 1 > self.index.shape[0] or index < -2:
            raise ValueError('index({}) out of range (-2-{})'.format(index, self.index.shape[0] - 1))
        if startIndex > self.index[0]:
            raise ValueError('startIndex({}) out of range (0-{})'.format(startIndex, self.index[0]))
        # ========特殊索引========
        if index == -1:
            set = self.data[startIndex:self.index[0]]
            timetage = self.timeTag[startIndex:self.index[0]].reshape((-1, 1))
            # 合并能量信息和时间信息-shape——>(-1,40)
            set = np.append(set, timetage, axis=1)
            # print("getData:{}".format(set))
            _d = pd.DataFrame(set, columns=_Index, dtype=dtype)
            return _d
        elif index == -2:
            _d = pd.DataFrame(columns=_Index, dtype=dtype)
            for i in range(self.index.shape[0]):
                _s = self.getData(i, startIndex=startIndex)
                _s['triggerID'] = _s['triggerID'] + i * 65535
                _d = pd.concat([_d, _s], ignore_index=True)
            return _d
        # ========普通索引========
        if index == 0:  # 选择第一个数据集时
            if self.index.shape[0] == 1:
                stop = self.index[0]
            else:
                stop = self.index[1]
            if startIndex >= stop:
                _d = pd.DataFrame(columns=_Index, dtype=dtype)
            else:
                set = np.empty((stop - startIndex, 39), dtype='uint16')
                self.data.read_direct(set, source_sel=np.s_[startIndex:stop], dest_sel=np.s_[0:set.shape[0]])
                timetag = np.empty((stop - startIndex,), dtype='int64')
                self.timeTag.read_direct(timetag, source_sel=np.s_[startIndex:stop], dest_sel=np.s_[0:set.shape[0]])
                timetage = timetag.reshape((-1, 1))
                # 合并能量信息和时间信息-shape——>(-1,40)
                set = np.append(set, timetage, axis=1)
                _d = pd.DataFrame(set, columns=_Index, dtype=dtype)
        elif index + 1 == self.index.shape[0]:  # 选择最后一个数据集时
            stop = self.index[0]
            start = self.index[index]
            if startIndex >= stop:
                _d = pd.DataFrame(columns=_Index, dtype=dtype)
            elif startIndex < start:
                print(self.data.shape)
                if stop - start != 0:
                    set = np.empty((stop - start, 39), dtype='uint16')
                    self.data.read_direct(set, source_sel=np.s_[start:stop], dest_sel=np.s_[0:set.shape[0]])
                    timetag = np.empty((stop - start,), dtype='int64')
                    self.timeTag.read_direct(timetag, source_sel=np.s_[start:stop], dest_sel=np.s_[0:set.shape[0]])
                    timetage = timetag.reshape((-1, 1))
                    # 合并能量信息和时间信息-shape——>(-1,40)
                    set = np.append(set, timetage, axis=1)
                    _d = pd.DataFrame(set, columns=_Index, dtype=dtype)
                else:
                    set = np.empty((stop - start, 40), dtype='int64')
                    _d = pd.DataFrame(set, columns=_Index, dtype=dtype)
            else:
                print(self.data.shape)
                if stop - start != 0:
                    set = np.empty((stop - startIndex, 39), dtype='uint16')
                    self.data.read_direct(set, source_sel=np.s_[startIndex:stop], dest_sel=np.s_[0:set.shape[0]])
                    timetag = np.empty((stop - startIndex,), dtype='int64')
                    self.timeTag.read_direct(timetag, source_sel=np.s_[startIndex:stop], dest_sel=np.s_[0:set.shape[0]])
                    timetage = timetag.reshape((-1, 1))
                    # 合并能量信息和时间信息-shape——>(-1,40)
                    set = np.append(set, timetage, axis=1)
                    _d = pd.DataFrame(set, columns=_Index, dtype=dtype)
                else:
                    set = np.empty((stop - start, 40), dtype='int64')
                    _d = pd.DataFrame(set, columns=_Index, dtype=dtype)
        else:  # 选择中间数据集时
            # 当判断到此时，index比self.index.shape[0]至少小2，比最大索引小1，index+1不会超出边界
            start = self.index[index]
            stop = self.index[index + 1]
            if startIndex >= stop:
                _d = pd.DataFrame(columns=_Index, dtype=dtype)
            elif startIndex < start:
                set = np.empty((stop - start, 39), dtype='uint16')
                self.data.read_direct(set, source_sel=np.s_[start:stop], dest_sel=np.s_[0:set.shape[0]])
                timetag = np.empty((stop - start,), dtype='int64')
                self.timeTag.read_direct(timetag, source_sel=np.s_[start:stop], dest_sel=np.s_[0:set.shape[0]])
                timetage = timetag.reshape((-1, 1))
                # 合并能量信息和时间信息-shape——>(-1,40)
                set = np.append(set, timetage, axis=1)
                _d = pd.DataFrame(set, columns=_Index, dtype=dtype)
            else:
                set = np.empty((stop - startIndex, 39), dtype='uint16')
                self.data.read_direct(set, source_sel=np.s_[startIndex:stop], dest_sel=np.s_[0:set.shape[0]])
                timetag = np.empty((stop - startIndex,), dtype='int64')
                self.timeTag.read_direct(timetag, source_sel=np.s_[startIndex:stop], dest_sel=np.s_[0:set.shape[0]])
                timetage = timetag.reshape((-1, 1))
                # 合并能量信息和时间信息-shape——>(-1,40)
                set = np.append(set, timetage, axis=1)
                _d = pd.DataFrame(set, columns=_Index, dtype=dtype)
        return _d

    # ==============内部=============
    # 获取triggerID重置次数
    def getSetsIndex(self) -> int:
        return self.setsIndex

    # 获取读取模式
    def getReadMode(self) -> bool:
        return self.readMode

    # 获取文件路径/文件名
    def getFileName(self) -> str:
        return self.file.filename

    # =============辅助===============
    # 将存储在缓冲区的数据输出到磁盘中去
    def flush(self):
        if self.readMode:
            self.data.refresh()
            self.timeTag.refresh()
            self.index.refresh()
            self.timeInfo.refresh()
            self.statusInfo.refresh()
        else:
            self.file.flush()

    # 表示triggerID重置，添加重置处索引并添加分片时间
    def newSets(self):
        '''
        set reset triggerID in index
        '''
        self.setsIndex += 1
        self.index.resize(self.index.shape[0] + 1, axis=0)
        self.timeInfo.resize(self.timeInfo.shape[0] + 1, axis=0)
        # 添加分片时间
        self.index[self.setsIndex] = self.dataIndex
        self.timeInfo[self.timeInfo.shape[0] - 1] = datetime.now().timestamp()

    # 关闭h5文件流
    def close(self):
        if self.file:
            if not self.readMode:
                print('h5 index:', self.dataIndex)
                self.data.resize(self.dataIndex, axis=0)
            self.file.close()

    # def __del__(self):
    #     self.close()
    # super(h5Data, self).__del__()


# 一个存储文件路径
class filePath(object):
    def __init__(self, **kwargs):
        self.p = kwargs.get("path", '')

    def get(self) -> str:
        return self.p

    def set(self, value: str):
        self.p = value


# server操作h5文件
class LUMISData():
    def __init__(self, path: str, mode: str = 'r', **kwargs):
        if mode[0] == 'w':
            self.file = h5py.File(path, mode=mode, libver='latest')
            self.readMode = False
            # 文件信息
            self.file.attrs["version"] = "0.2"
            self.file.attrs["project"] = kwargs.get("project", "unknown")
            self.file.attrs["create time"] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            self.file.attrs["note"] = ""
        if mode[0] == 'r':
            self.file = h5py.File(path, mode=mode, libver='latest', swmr=True)
            # 检查读取文件的版本
            try:
                self.file.attrs["version"]
            except KeyError:
                raise ValueError(
                    "LUMISData can only read Lumis Data version = 0.2. if version = 0.1, please try h5Data to load it.")
            if self.file.attrs["version"] != "0.2":
                raise ValueError(
                    "LUMISData can only read Lumis Data version = 0.2. if version = 0.1, please try h5Data to load it.")

            self.readMode = True

    def version(self):
        return self.file.attrs["version"]

    def showNote(self):
        return self.file.attrs["note"]

    def project(self):
        return self.file.attrs["project"]

    def createTime(self):
        return self.file.attrs["create time"]

    def deviceList(self):
        devList = set()
        for dataPath in self.file.keys():
            deviceName = self.file[dataPath].attrs["name"]
            devList.add(deviceName)
        return list(devList)

    def dataPathList(self, deviceName=None):
        if deviceName is None:
            return list(self.file.keys())
        else:
            result = []
            for dataPath in self.file.keys():
                name = self.file[dataPath].attrs["name"]
                if name == deviceName:
                    result.append(dataPath)
            return result

    def getDevData(self, path: str):
        return LUMISDevData(self.file[path], self)

    def newDevData(self, path, deviceName: str):
        dataGroup = self.file.create_group(path)
        data = dataGroup.create_dataset("data", shape=(5000, 39), maxshape=(None, 39), dtype='uint16')
        index = dataGroup.create_dataset("index", shape=(1,), maxshape=(None,), dtype='int64')
        timeTag = dataGroup.create_dataset("timeTag", shape=(5000,), maxshape=(None,), dtype='uint64')
        device_status = dataGroup.create_dataset("device_status", shape=(17,), dtype='uint16')
        t = time.time()
        dataGroup.attrs["start time"] = t
        dataGroup.attrs["stop time"] = t
        index[0] = 0
        dataGroup.attrs["name"] = deviceName
        return LUMISDevData(dataGroup, self)

    # 获取读取模式
    def getReadMode(self) -> bool:
        return self.readMode


# 一个设备的数据
class LUMISDevData():
    def __init__(self, dataGroup: h5py.Group, parent: LUMISData):
        self.parent = parent
        self.data = dataGroup['data']
        self.timeTag = dataGroup['timeTag']
        self.index = dataGroup['index']
        self.DevName = dataGroup.attrs["name"]
        self.statusInfo = dataGroup['device_status']
        self.dataGroup = dataGroup
        self.dataIndex = self.index[0]  # 当前写到第几行
        self.setsIndex = self.index.shape[0] - 1  # 当前存在几次triggerID置零

    def setStartTime(self, t=None):
        if not self.parent.readMode:
            self.dataGroup.attrs["start time"] = time.time() if t is None else t

    def setStopTime(self, t=None):
        if not self.parent.readMode:
            self.dataGroup.attrs["stop time"] = time.time() if t is None else t

    def startTime(self):
        return self.dataGroup.attrs["start time"]

    def stopTime(self):
        return self.dataGroup.attrs["stop time"]

    def totalTime(self):
        return self.dataGroup.attrs["stop time"] - self.dataGroup.attrs["start time"]

    # 录入设备状态信息
    def putDeviceStatus(self, status: dict):
        '''
        "w": write device configuration state to h5 file
        "r": raise Exception
        :param status:dict:key:boardID;value:tuple(threshold, biasVoltage, triggerMode)
        :return:
        '''
        print(status)
        triggerMode = status[0][2]
        self.statusInfo[0] = triggerMode
        for i in range(8):
            t = status.get(i, (0, 0, 0))
            self.statusInfo[i + 1] = t[0]
            self.statusInfo[i + 9] = t[1]

    # 获取h5文件中存储的设备状态信息
    def getDeviceStatus(self) -> dict:
        '''
        :return: dict:key:boardID;value:tuple(threshold, bias voltage, trigger mode)
        '''
        status = {}
        triggerMode = self.statusInfo[0]
        for i in range(8):
            threshold = self.statusInfo[i + 1]
            biasVoltage = self.statusInfo[i + 9]
            if threshold == 0 and biasVoltage == 0:
                break
            else:
                status[i] = (threshold, biasVoltage, triggerMode)
        return status

    # ===============数据=================
    # 将数据添加到h5数据集中
    def addToDataSet(self, value: np.array, timeTag=None):
        '''
        "w": add data to h5 file
        "r": raise Exception
        :param value: data,type:np.array,dtype:'uint16',shape:(n,39)
                    column 0-38 will be pushed into group(data), column 39 will be pushed into group(timeTag)
        :param NewSets: if it set as True,h5 file will create a new data set to contain data
        :return:
        '''
        # 如果数据集的长度不够，则扩展5000行
        if self.dataIndex + value.shape[0] <= self.data.shape[0]:
            self.data.resize(self.data.shape[0] + 5000, axis=0)
            self.timeTag.resize(self.data.shape[0] + 5000, axis=0)
        # 数据存储
        self.data.write_direct(value, source_sel=np.s_[0:value.shape[0]],
                               dest_sel=np.s_[self.dataIndex:value.shape[0] + self.dataIndex])
        if timeTag is not None:
            self.timeTag.write_direct(timeTag, source_sel=np.s_[0:value.shape[0]],
                                      dest_sel=np.s_[self.dataIndex:value.shape[0] + self.dataIndex])
        self.dataIndex += value.shape[0]
        self.index[0] = self.dataIndex

        # 表示triggerID重置，添加重置处索引并添加分片时间

    # 读取h5文件中的数据
    def getData(self, index: int = -1, startIndex: int = 0, dtype='int64') -> pd.DataFrame:
        '''
        read data in h5 file
        :param index: When the index greater than or equal to 0, the corresponding set of data will be returned.
        when the index equal to -1,all the data will be returned.
        when the index equal to -2,all the data will be returned and the triggerID in returned data has been accumulated.
        :param startIndex: data will be returned from the startIndex
        :return:
        '''
        if index + 1 > self.index.shape[0] or index < -2:
            raise ValueError('index({}) out of range (-2-{})'.format(index, self.index.shape[0] - 1))
        if startIndex > self.index[0]:
            raise ValueError('startIndex({}) out of range (0-{})'.format(startIndex, self.index[0]))
        # ========特殊索引========
        if index == -1:
            set = self.data[startIndex:self.index[0]]
            timetage = self.timeTag[startIndex:self.index[0]].reshape((-1, 1))
            # 合并能量信息和时间信息-shape——>(-1,40)
            set = np.append(set, timetage, axis=1)
            # print("getData:{}".format(set))
            _d = pd.DataFrame(set, columns=_Index, dtype=dtype)
            return _d
        elif index == -2:
            _d = pd.DataFrame(columns=_Index, dtype=dtype)
            for i in range(self.index.shape[0]):
                _s = self.getData(i, startIndex=startIndex)
                _s['triggerID'] = _s['triggerID'] + i * 65535
                _d = _d.append(_s, ignore_index=True)
            return _d
        # ========普通索引========
        if index == 0:  # 选择第一个数据集时
            if self.index.shape[0] == 1:
                stop = self.index[0]
            else:
                stop = self.index[1]
            if startIndex >= stop:
                _d = pd.DataFrame(columns=_Index, dtype=dtype)
            else:
                set = np.empty((stop - startIndex, 39), dtype='uint16')
                self.data.read_direct(set, source_sel=np.s_[startIndex:stop], dest_sel=np.s_[0:set.shape[0]])
                timetag = np.empty((stop - startIndex,), dtype='int64')
                self.timeTag.read_direct(timetag, source_sel=np.s_[startIndex:stop], dest_sel=np.s_[0:set.shape[0]])
                timetage = timetag.reshape((-1, 1))
                # 合并能量信息和时间信息-shape——>(-1,40)
                set = np.append(set, timetage, axis=1)
                _d = pd.DataFrame(set, columns=_Index, dtype=dtype)
        elif index + 1 == self.index.shape[0]:  # 选择最后一个数据集时
            stop = self.index[0]
            start = self.index[index]
            if startIndex >= stop:
                _d = pd.DataFrame(columns=_Index, dtype=dtype)
            elif startIndex < start:
                set = np.empty((stop - start, 39), dtype='uint16')
                self.data.read_direct(set, source_sel=np.s_[start:stop], dest_sel=np.s_[0:set.shape[0]])
                timetag = np.empty((stop - start,), dtype='int64')
                self.timeTag.read_direct(timetag, source_sel=np.s_[start:stop], dest_sel=np.s_[0:set.shape[0]])
                timetage = timetag.reshape((-1, 1))
                # 合并能量信息和时间信息-shape——>(-1,40)
                set = np.append(set, timetage, axis=1)
                _d = pd.DataFrame(set, columns=_Index, dtype=dtype)
            else:
                set = np.empty((stop - startIndex, 39), dtype='uint16')
                self.data.read_direct(set, source_sel=np.s_[startIndex:stop], dest_sel=np.s_[0:set.shape[0]])
                timetag = np.empty((stop - startIndex,), dtype='int64')
                self.timeTag.read_direct(timetag, source_sel=np.s_[startIndex:stop], dest_sel=np.s_[0:set.shape[0]])
                timetage = timetag.reshape((-1, 1))
                # 合并能量信息和时间信息-shape——>(-1,40)
                set = np.append(set, timetage, axis=1)
                _d = pd.DataFrame(set, columns=_Index, dtype=dtype)
        else:  # 选择中间数据集时
            # 当判断到此时，index比self.index.shape[0]至少小2，比最大索引小1，index+1不会超出边界
            start = self.index[index]
            stop = self.index[index + 1]
            if startIndex >= stop:
                _d = pd.DataFrame(columns=_Index, dtype=dtype)
            elif startIndex < start:
                set = np.empty((stop - start, 39), dtype='uint16')
                self.data.read_direct(set, source_sel=np.s_[start:stop], dest_sel=np.s_[0:set.shape[0]])
                timetag = np.empty((stop - start,), dtype='int64')
                self.timeTag.read_direct(timetag, source_sel=np.s_[start:stop], dest_sel=np.s_[0:set.shape[0]])
                timetage = timetag.reshape((-1, 1))
                # 合并能量信息和时间信息-shape——>(-1,40)
                set = np.append(set, timetage, axis=1)
                _d = pd.DataFrame(set, columns=_Index, dtype=dtype)
            else:
                set = np.empty((stop - startIndex, 39), dtype='uint16')
                self.data.read_direct(set, source_sel=np.s_[startIndex:stop], dest_sel=np.s_[0:set.shape[0]])
                timetag = np.empty((stop - startIndex,), dtype='int64')
                self.timeTag.read_direct(timetag, source_sel=np.s_[startIndex:stop], dest_sel=np.s_[0:set.shape[0]])
                timetage = timetag.reshape((-1, 1))
                # 合并能量信息和时间信息-shape——>(-1,40)
                set = np.append(set, timetage, axis=1)
                _d = pd.DataFrame(set, columns=_Index, dtype=dtype)
        return _d

    def newSets(self):
        '''
        set reset triggerID in index
        '''
        self.setsIndex += 1
        self.index.resize(self.index.shape[0] + 1, axis=0)
        # 添加分片时间
        self.index[self.setsIndex] = self.dataIndex

    def refreshIndex(self, indexList: list):
        if self.index.shape[0] - 1 < len(indexList):
            lNum = len(indexList)
            for i in range(len(indexList)):
                if self.dataIndex < indexList[i]:
                    lNum = i
            self.index.resize(lNum + 1, axis=0)
            self.index[1:lNum + 1] = np.array(indexList[:lNum])

    def __del__(self):
        if self.parent.file:
            if not self.parent.readMode:
                print('h5 index:', self.dataIndex)
                self.data.resize(self.dataIndex, axis=0)
