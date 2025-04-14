# author:高金磊
# datetime:2021/11/13 11:20

'''
 这种查找方式快很多但是效果不好
class topo_tools:
    def __init__(self,file="data/58Mamian_surface_1m_withoutDuoZi_52814.topo"):
        self._file=file
        data = []
        self._data_tree:dict={}
        with open(file=file, encoding='utf-8') as file_obj:
            line = file_obj.readline()#舍弃第一行
            while 1:
                line = file_obj.readline()
                if not line:
                    break
                middle=[]
                for datum in line.replace("\n","").split(" "):
                    if datum !="":
                        middle.append(float(datum))
                data.append(middle)
        file_obj.close()
        for xyz in data:
            x,y,z = xyz
            if not self._data_tree.keys().__contains__(x):
                self._data_tree.__setitem__(x,{})
            xs:dict=self._data_tree.get(x)
            if not xs.__contains__(y):
                xs.__setitem__(y,[])
            xs.get(y).append(z)
        self._xs=sorted([i for i in self._data_tree.keys()])
        del data
        pass
    def  check_point(self,x,y,z):
        middle=1000000
        min_x=0
        for i in self._xs:
            if middle<abs(x-i):
                break
            middle=(x-i)
            min_x=i
        ys=self._data_tree.get(min_x).keys()
        middle=1000000
        min_y=0
        for i in ys:
            if middle>abs(i-y):
                middle=abs(i-y)
                min_y=i
            else:
                break
        return [z<=self._data_tree.get(min_x).get(min_y)[0],[min_x,min_y,self._data_tree.get(min_x).get(min_y)[0]]]

'''


def sort_key(element):
    """
    排序的依据

    :param element: 需要排序的数据
    :return:
    """
    return element[0], element[1]


class topo_tools:
    """
    读取”topo“文件 ,初始化 Topo工具
    """
    def __init__(self, file):
        """
        通过读取文件来构建
        :param file: topo文件
        """
        self._file = file
        self._data = []
        self._x_cache = {}
        self.__point_cache_x = {}
        self.cache_count = 0
        with open(file=file, encoding='utf-8') as file_obj:
            line = file_obj.readline()  # 舍弃第一行
            while 1:
                line = file_obj.readline()
                if not line:
                    break
                middle = []
                for datum in line.replace("\n", "").split(" "):
                    if datum != "":
                        middle.append(float(datum))
                self._data.append(middle)
        file_obj.close()
        # 排序来加快check_point模式
        self._data.sort(key=sort_key)
        self.init_x_cache()
        pass

    def init_x_cache(self):
        """
        缓存x加快索引,x按照0.25一个单位进行离散化

        :return:
        """
        for _datum in self._data:
            x = int(_datum[0] * 4)
            middle = self._x_cache.get(x)
            if not middle:
                middle = []
            middle.append(_datum)
            self._x_cache[x] = middle
        pass

    def get_x_range(self, x, y, loose=1):
        """
        获取x0属于[x-loose/4,x-loose/4]且|y-y0|<1的所有数据

        :param x: x 坐标
        :param y: y 坐标
        :param loose: 宽松度
        :return: (x,y)附近的所有的topo信息
        """
        x = int(x * 4)
        x_range = [x + i for i in range(-loose, loose)]
        res = []
        for target in x_range:
            data = self._x_cache.get(target)
            if not data:
                break
            for datum in data:
                if abs(datum[1] - y) <= 1:
                    res.append(datum)
        return res
        pass

    def check_point(self, x, y, z):
        """
        判断格子是不是在topo文件上方_离散化toppo
        :param x:
        :param y:
        :param z:
        :return: list : [是不是在topo文件上方,用来判定的坐标,需要判定的坐标]
        """
        target = [0, 0, 0]
        middle = 100000
        cache = self.get_x_range(x, y, loose=1)
        if len(cache) > 0:
            # 从x缓存中进一步找到匹配项
            for re in cache:
                distinc = (x - re[0]) ** 2 + (y - re[1]) ** 2
                if distinc < middle:
                    target = [*re]
                    middle = distinc
            return [z > target[2], target, [x, y, z]]
        else:
            for data in self._data:
                # if abs(x - data[0]) > 0.5*2 or (y - data[1]) > 0.25*2:
                #     continue
                distinc = (x - data[0]) ** 2 + (y - data[1]) ** 2
                if distinc < middle:
                    target = [*data]
                    middle = distinc
            return [z > target[2], target, [x, y, z]]

    def check_point_cache(self, x, y, z, cache=True):
        """
        判断格子是不是在topo文件上方--缓存每次判断的结果中不为空气的最大值

        :param x:
        :param y:
        :param z:
        :param cache: 是否使用缓存,使用将降低一些精度,但是可以大幅度提升速度(约1000倍)
        :return: [是不是在topo文件上方,用来判定的坐标,需要判定的坐标]
        """
        if not cache:
            return self.check_point(x, y, z)
        if self.__point_cache_x.keys().__contains__(x):
            point_cache_y = self.__point_cache_x.get(x)
            if point_cache_y.keys().__contains__(y):
                point_cache_z = point_cache_y.get(y)
                # 由于坐标系因素这里z大是空气
                if point_cache_z[1][2] > z:
                    self.cache_count += 1
                    return point_cache_z
                else:
                    self.cache_count += 1
                    middle = [*point_cache_z]
                    middle[0] = True
                    return middle
                    # print("请将z从大到小排，否则效果不佳")
        point_cache_z = self.check_point(x, y, z)
        if point_cache_z[0]:
            # 非空气，缓存真实的数据
            if self.__point_cache_x.keys().__contains__(x):
                point_cache_y = self.__point_cache_x.get(x)
            else:
                point_cache_y = {}
                self.__point_cache_x[x] = point_cache_y
            self.__point_cache_x[x] = point_cache_y
            # 要求z从大到小排列
            # point_cache_y[y] = point_cache_z
            point_cache_y[y] = [False, point_cache_z[1], [x, y, z]]
        else:
            # 空气，缓存topo的数据
            if self.__point_cache_x.keys().__contains__(x):
                point_cache_y = self.__point_cache_x.get(x)
            else:
                point_cache_y = {}
                self.__point_cache_x[x] = point_cache_y
            self.__point_cache_x[x] = point_cache_y
            point_cache_y[y] = [False, point_cache_z[1], [x, y, point_cache_z[1][2]]]
        return point_cache_z

    def get_cache_count(self):
        """
        返回自从上次调用以来有多少数据被加速

        :return:
        """
        middle = self.cache_count
        self.cache_count = 0
        return middle


if __name__ == '__main__':
    obj = topo_tools(r"E:\vscode\Muon_Imaging_Algorithm\dataTools\data\58mamian_survey_821486.topo")
    obj.check_point(-105, 14.7, 9.9)
    obj.check_point(-7.13, -12.58, 11.72)
    pass
