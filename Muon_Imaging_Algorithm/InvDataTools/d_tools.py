# author:高金磊
# datetime:2021/11/18 16:05
"""
#d不在从单独的文件中获取，而是从观测文件中获取d和误差d——err
 def __init__(self, file=r"E:\vscode\Muon_Imaging_Algorithm\dataTools\data\d"):
    self._file = file
    self._data = []
    with open(file=file, encoding='utf-8') as file_obj:
        while 1:
            line = file_obj.readline()
            if not line:
                break
            self._data.append(float(line))
"""


class d_tools:
    """
    读取并处理“d"文件

    """
    def __init__(self, file):
        self._data = []
        self._data_err = []
        with open(file=file, mode="r") as file_data:
            # file_data.readline()
            while 1:
                line = file_data.readline()
                if not line:
                    break
                lines = line.split()
                self._data.append(float(lines[0]))
                self._data_err.append(float(lines[1]))

    def get_d(self):
        """
        获取观测的d

        :return: 观测的d
        """
        return self._data

    def get_d_err(self):
        """
        获取d的绝对误差

        :return: d的绝对误差
        """
        return [i for i in self._data_err]

    def get_d_err_range(self):
        """
        获取d的取值范围，上下限，其中上限与下限互为相反数

        :return:[下限,上限]，也即[-d的绝对误差，d的绝对误差]
        """
        d_lt = []
        d_lg = []
        for i in range(len(self._data)):
            # 相对误差-取值区间
            # d_lt.append(self._data[i]*(1+self._data_err[i]))
            # d_lg.append(self._data[i]*(1-self._data_err[i]))

            # 绝对误差-取值区间
            # d_lt.append(self._data[i] + self._data_err[i])
            # d_lg.append(self._data[i] - self._data_err[i])
            # 绝对误差
            d_lt.append(self._data_err[i])
            d_lg.append(-self._data_err[i])
        return [d_lg, d_lt]


if __name__ == '__main__':
    d = d_tools()
    d.get_d()
    d.get_d_err()
    d.get_d_err_range()
    pass
