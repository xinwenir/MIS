# author:高金磊
# datetime:2021/11/13 10:09
from copy import copy


class obs_tools:
    """
    读取并处理”obs“文件

    """
    def __init__(self, file):
        self._file = file
        self._data = []
        with open(file=file, encoding='utf-8') as file_obj:
            file_obj.readline()  # 去掉第一个其他数据
            while 1:
                line = file_obj.readline()
                if not line:
                    break
                middle = []
                for datum in line.replace("\n", "").split(" "):
                    if datum != "":
                        middle.append(datum)


                self._data.append(middle)

    def make_d_from_obs(self, d_file):
        """
        将obs中关于d的数据(射线的等效长度d、d的相对误差d_err)存放到文件中

        :param d_file: 存放d数据的文件路径(写入)
        """
        file = open(d_file, 'w')
        for _datum in self._data:
            file.write(_datum[-2])
            file.write(' ')
            file.write(_datum[-1])
            file.write('\n')
        file.close()

    def get_obs_id_count(self):
        """
        获取每个探测器的射线总数量

        :return: 每个探测器的射线总数量
        """
        res = []
        count = 0
        id_middle = int(self._data[0][0])
        for _datum in self._data:
            id = int(_datum[0])
            if id_middle == id:
                count += 1
            else:
                id_middle = id
                res.append(count)
                count = 1
        if count != 0:
            res.append(count)
        return res

    def get_data(self):
        """
        获取所有的射线信息

        :return: 所有的射线信息
        """
        return copy(self._data)

    def shape(self):
        """
        获取探测信息的总行数以及每个探测信息包含数据的个数

        :return: 探测信息的总行数以及每个探测信息包含数据的个数
        """
        return [len(self._data) - 1, len(self.get_data()[0])]

    def get_d_form_obs(self):
        """
        获取d的值(射线的等效长度)

        :return: d的值(射线的等效长度)
        """
        return [i[-2] for i in self.get_data()]

    def get_d_absolute_err_form_obs(self):
        """
        获取d的相对误差

        :return: d的相对误差
        """
        return [i[-1] for i in self.get_data()]

    def get_d_relative_error_form_obs(self):
        """
        获取d的相对误差

        :return: d的相对误差
        """
        return [float(i[-1]) / float(i[-2]) for i in self.get_data()]

    def get_receiver_list(self):
        """
        获取所有探测器的位置坐标

        :return: 所有探测器的位置坐标
        """
        res = []
        for i in self._data:
            res.append([float(i[1]), float(i[2]), float(i[3])])
        return res
    def screening_data(self,Threshold=1.5,new_obs_file=None):
        middle=[]
        for data in self._data:
            if abs(float(data[7]) - float(data[6]) * 1.89) / float(data[8]) < 1.5:
                middle.append(data)
        if new_obs_file is not None:
            new_obs=open(new_obs_file,'w')
            for i in middle:
                line=""
                for union in i:
                    line+=union+" "
                new_obs.write(line)
                new_obs.write('\n')
        return middle


if __name__ == '__main__':
    # obj = obs_tools(r"M:\pycharm\Inversion\DataTools\data\cube01_obs.dat")
    # obj.make_d_from_obs(r'M:\pycharm\Inversion\DataTools\data\d')
    # print(obj.get_obs_id_count())
    # print(obj.get_data()[0])
    # print(obj.shape())
    # print(obj.get_d_relative_error_form_obs())
    obj = obs_tools(r"M:\pycharm\Inversion\data\Input\real_data\27_obs.dat")
    obj.screening_data(Threshold=2,new_obs_file=r"M:\pycharm\Inversion\data\Input\real_data\27gjl_obs.dat")
