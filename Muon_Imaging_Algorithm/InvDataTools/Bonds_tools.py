# author:高金磊
# datetime:2022/3/23 10:26
import InvDataFactory
from InvDataTools.res_tools import res_tools


class Bonds_tool:
    """
    读取“Bounds”文件并存储相应的数据

    """

    def __init__(self, boods_file,
                 bound_min_value=None, bound_max_value=None):
        self.boods_file=boods_file
        file = open(boods_file, mode='r')
        lines = file.readlines()
        data = []
        for line in lines:
            min_max = line.split()
            min = float(min_max[0])
            if bound_min_value is not None:
                if min < bound_min_value:
                    min = bound_min_value
            max = float(min_max[1])
            if bound_max_value is not None:
                if max > bound_max_value:
                    max = bound_max_value
            if max < min:
                max = min
            data.append([min, max])
            # if max>2:
            #     data.append([-0.4, 0.41])
            # else:
            #     data.append([min, max])
        self.__data = data
        file.close()

    def get_bonds_min_max(self):
        """
        获取格子密度约束的左边界和右边界

        :return: 格子密度约束的左边界和右边界
        """
        return self.__data

    def get_bonds_min(self):
        """
        获取格子密度约束的左边界

        :return: 格子密度约束的左边界
        """
        return [i[0] for i in self.__data]

    def get_bonds_max(self):
        """
        获取格子密度约束的右边界

        :return: 格子密度约束的右边界
        """
        return [i[1] for i in self.__data]

    def show_min(self, res_file):
        """
        将格子密度约束的左边界存放到文件中

        :param res_file: 存放格子密度约束左边界的文件路径
        """
        file = open(file=res_file, mode='w')
        for i in self.get_bonds_min():
            file.write(str(i))
            file.write('\n')
        file.close()
        print("已经写入到:", res_file)

    def show_max(self, res_file):
        """
        将格子密度约束的右边界存放到文件中

        :param res_file: 存放格子密度约束右边界的文件路径
        """
        file = open(file=res_file, mode='w')
        for i in self.get_bonds_max():
            file.write(str(i))
            file.write('\n')
        file.close()
        print("已经写入到:", res_file)
    def update_data(self,min_max):
        """
        将当前结果更新到文件和对象中
        :param min_max: list(list)
        :return:
        """
        self.__data=min_max
        file=open(self.boods_file,'w')
        for data in min_max:
            file.write(str(data[0]))
            file.write(" ")
            file.write(str(data[1]))
            file.write("\n")

        file.close()


class Make_Bounds_Tools:
    """
    生成"bounds"文件

    """

    def __init__(self, airj, shape):
        # 此处可以兼容数组,文件等
        self.airj = airj
        self.shape = shape

    def save_bounds(self, path, air_bound, wall_bound, inner_bound):  #####写对应的设置文件
        """
        确定每个格子的类型(空气、墙皮、内部城墙)并对格子的密度约束赋值，所有格子处理完后，将密度约束存储在文件中

        :param path: 存放结果的文件路径
        :param air_bound: 空气格子的密度约束
        :param wall_bound: 墙皮格子的密度约束
        :param inner_bound: 内部的城墙格子的密度约束
        :return:
        """
        from InvDataTools.Jxyz_Tools import getj_from_xyz
        shape = self.shape
        airj = self.airj
        if type(self.airj) is not set:
            airj = set(self.airj)
        import numpy
        res = numpy.zeros((shape[0] * shape[1] * shape[2],), dtype=list)
        for x in range(shape[0]):
            x += 1
            for y in range(shape[1]):
                y += 1
                # count=0
                for z in range(shape[2], 0, -1):
                    j = getj_from_xyz(shape, (x, y, z))
                    if j in airj:
                        res[j - 1] = air_bound
                    else:
                        if z == shape[2]:
                            res[j - 1] = wall_bound
                        else:
                            flag = []
                            if x + 1 <= shape[0]:
                                flag.append(getj_from_xyz(shape, (x + 1, y, z)) in airj)
                            if x - 1 > 0:
                                flag.append(getj_from_xyz(shape, (x - 1, y, z)) in airj)
                            if y + 1 <= shape[1]:
                                flag.append(getj_from_xyz(shape, (x, y + 1, z)) in airj)
                            if y - 1 > 0:
                                flag.append(getj_from_xyz(shape, (x, y - 1, z)) in airj)
                            if z - 1 > 0:
                                flag.append(getj_from_xyz(shape, (x, y, z - 1)) in airj)  # 已经通过count简化
                            if z + 1 <= shape[2]:
                                flag.append(getj_from_xyz(shape, (x, y, z + 1)) in airj)
                            if flag.count(True) > 0:
                                res[j - 1] = wall_bound
                            else:
                                res[j - 1] = inner_bound
                        # count+=1
        file = open(path, 'w')
        for value in res:
            file.write("%s %s" % (str(value[0]), str(value[1])))
            file.write('\n')
        return True


if __name__ == '__main__':
    # import Setting
    #
    # setting = Setting.Setting.get_instance()
    #
    # tool = Bonds_tool(setting.get_bnd_file())
    # tool.show_max(r'M:\pycharm\Inversion\InvDataTools\data\bonds_max')
    # rt = res_tools()
    # rt.Conversion_2(shape=(280, 144, 63), file_res=r'M:\pycharm\Inversion\InvDataTools\data\bonds_max',
    #                 file_out=r'M:\pycharm\Inversion\InvDataTools\data\bonds_max')
    import InvDataFactory.DataManage

    datamanage = InvDataFactory.DataManage.DataManager()
    # setting=Setting.Setting()
    tool = Make_Bounds_Tools(datamanage.get_unneed_j(), datamanage.mesh.get_shape())
    tool.save_bounds(r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\bounds", (-0.001, 0.001), (2.39, 2.41), (0, 2.1))
