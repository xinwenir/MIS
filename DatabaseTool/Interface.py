# -*- coding: utf-8 -*-
# @Author  : jimbook
# @email   : tianh2021@lzu.edu.cn
# @Time    : 2022/9/15
# @Note    : an interface for all cormis dataReader class
# @modify  : None


from DatabaseTool import *
#数据读取器
class Cormis_DataReader_Interface():
    # 获取对应数据表的开始测量时间
    def getStartTime(self, tableName: str) -> pd.Timestamp:
        pass

    # 获取对应数据表的结束测量时间
    def getStopTime(self, tableName: str) -> pd.Timestamp:
        pass

    # 获取对应数据表的行数（数据量）
    def getCount(self, tableName: str) -> int:
        pass

    # 获取对应数据集合内的数据信息，返回包括数据表名、开始时间、结束时间、数据量
    def deviceDatatablesInfo(self) -> pd.DataFrame:
        pass

    # 获取对应数据表的数据
    def deviceData(self, tableName:str, start:int, stop:int) -> pd.DataFrame:
        pass

    # 以迭代器形式返回数据表的数据，每次迭代的数据量由用户指定
    def deviceDataSplitIter(self,tableName:str, pageCapacity:int):
        pass

# 数据切片迭代器-接口定义
class Cormis_DataIterator_Interface():
    def __init__(self, parent:Cormis_DataReader_Interface, tableName:str, pageCapacity:int, reindex:bool = True):
        '''
        :param parent: 隶属的数据读取器
        :param tableName: 读取的数据表名
        :param pageCapacity: 每次迭代期望拿到的数据量，为了融合
        '''
        self.parent = parent
        self.tableName = tableName
        self.pageCapacity = pageCapacity if pageCapacity > 100 else 100

    # 获取最后被切开的Event行数
    def getTriggerIDOffset(self,triggerID:np.array):
        '''
        :param triggerID: 需要检查的triggerID一维数组
        :return:
        '''
        count = triggerID.shape[0]-1
        endTriggerID = triggerID[count]
        for i in range(count):
            if endTriggerID !=triggerID[count-i]:
                return i
        return 0

    #获取数据逻辑
    def getDataSPlited(self,start:int, stop:int):
        data = self.parent.deviceData(self.tableName, start, stop)
        return data

    # 初始化迭代器
    def __iter__(self):
        self.endTag = False
        self.startID = 0
        return self

    # 迭代逻辑
    def __next__(self):
        if self.endTag:
            raise StopIteration
        else:
            endID = self.startID + self.pageCapacity
            d = self.getDataSPlited(self.startID,endID)#self.parent.deviceData(self.tableName, self.startID, endID)
            if d.shape[0] == 0:
                raise StopIteration
            elif d.shape[0] < self.pageCapacity:
                offset = 0
                self.endID = True
            else:
                offset = self.getTriggerIDOffset(d["triggerID"].values)
            data = d.iloc[:d.shape[0]-offset]
            self.startID += data.shape[0]
            return data