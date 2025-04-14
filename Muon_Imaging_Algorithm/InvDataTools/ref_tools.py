# author:高金磊
# datetime:2022/3/24 14:09
import numpy


class Ref_tools:
    """
    读取并处理”ref“文件

    """
    def __init__(self, path):
        self.ref_file=path
        file = open(path, 'r')
        self.data = []
        lines = file.readlines()
        for line in lines:
            self.data.append(float(line))
            # self.data.append(0.0)

        file.close()

    def get_data(self):
        """
        获取格子密度的参考值

        :return: 格子密度的参考值
        """
        return self.data

    def recover_resj_by_refj(self, res, refs=None):
        """
        将计算结果中的0值用对应的参考值代替

        :param res: 格子密度的计算值(当refs=None时,此处的值需为恢复后的值)
        :param refs: 格子密度的参考值
        :return:
        """
        if refs is None:
            refs = self.get_data()
        for i in range(len(res)):
            if res[i] == 0:  # 空气也是0吗??会不会对平滑度有影响
                res[i] = refs[i]
        return res

    def make_refs_ps_err_list(self, refs, res):
        """
        获取格子密度的差值(=密度计算值-密度参考值)

        :param refs: 格子密度的参考值
        :param res: 格子密度的计算值
        :return: 格子密度的差值
        """
        data = []
        for i in range(len(refs)):
            data.append(res[i] - refs[i])
        return data

    def make_refs_ps_err_file(self, refs_file, ps, out_file):
        """
        此方法已删除

        将refs和ps之间的差值存放到文件中

        :param refs_file: 格子密度参考值的文件路径(读取)
        :param ps: 格子密度的计算值
        :param out_file: 存放结果的文件路径(写入)
        """
        import warnings
        warnings.warn("make_refs_ps_err_file方法已经过时,删除原因:灵活度太低", DeprecationWarning)
        """
        显示refs和ps之间的差值
        :param refs:
        :param ps:
        :param out_file:
        :return:
        """
        try:
            file = open(out_file, 'w')
            refs = Ref_tools(refs_file).get_data()
            for i in range(len(refs)):
                file.write(str(abs(refs[i] - ps[i])))
                file.write('\n')
            file.close()
        except Exception as e:
            print(e)
    def update_data(self,refs):
        """
        将refs更新到文件和对象中
        :param min_max: list(float)
        :return:
        """
        self.data=refs
        
        file=open(self.ref_file,'w')
        for data in refs:
            file.write(str(data))
            file.write("\n")

        file.close()


class Make_refcence_tools():
    """
    生成”ref“文件

    """
    def __init__(self, airj, shape):
        # 此处可以兼容数组,文件等
        self.airj = airj
        self.shape = shape

    def save_refs(self, path, air_value, wall_value, inner_value):  #####写对应的设置文件
        """
        保存密度参考值

        :param path: 存放结果的文件路径(写入)
        :param air_value: 空气密度的参考值
        :param wall_value: 墙皮密度的参考值
        :param inner_value: 城墙内部密度的参考值
        :return: 定值True
        """
        from InvDataTools.Jxyz_Tools import getj_from_xyz
        shape = self.shape
        airj = self.airj
        if type(self.airj) is not set:
            airj = set(self.airj)
        res = numpy.zeros((shape[0] * shape[1] * shape[2],), dtype=numpy.float)
        for x in range(shape[0]):
            x += 1
            for y in range(shape[1]):
                y += 1
                # count=0
                for z in range(shape[2], 0, -1):
                    j = getj_from_xyz(shape, (x, y, z))
                    if j in airj:
                        res[j - 1] = air_value
                    else:
                        if z == shape[2]:
                            res[j - 1] = wall_value
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
                                res[j - 1] = wall_value
                            else:
                                res[j - 1] = inner_value
                        # count+=1
        file = open(path, 'w')
        for value in res:
            file.write(str(value))
            file.write('\n')
        return True


if __name__ == '__main__':
    import InvDataFactory.DataManage

    datamanage = InvDataFactory.DataManage.DataManager()
    # setting=Setting.Setting()
    tool = Make_refcence_tools(datamanage.get_unneed_j(), datamanage.mesh.get_shape())
    tool.save_refs(r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\refs", -10, 2.4, 1.89)
