# author:高金磊
# datetime:2022/2/24 15:59
from copy import copy

import numpy

from DataTools.Jxyz_Tools import getj_from_xyz, getxyz_from_shape


class new_res_tools():
    """
    对结果进行一些处理

    """

    def __init__(self, res_list=None, res_file=None):
        if res_list is not None:
            self.data = copy(res_list)
            return
        if res_file is None:
            raise Exception("两个参数不可都为空")
        with open(file=res_file, mode='r') as file:
            line = file.readline().replace('\n', '')
            self.data = []
            while line:
                self.data.append(float(line))
                line = file.readline().replace('\n', '')
        file.close()

    def mode0(self, target_file):
        """
        将结果存到文件中

        :param target_file: 存放结果的文件路径
        """
        with open(target_file, 'w') as file:
            for re in self.data:
                file.write(str(re))
                file.write('\n')
            file.close()

    def mode1(self, min, max, ignore_value, target_file=r"M:\pycharm\Inversion\DataTools\data\new_mode1"):
        """
        检查结果，将结果中不在[min，max]的数据调整为min、max数据

        :param min: 最小值
        :param max: 最大值
        :param ignore_value: 忽略的值，即当ignore_value包含值value时，值value不需要判断是否位于[min,max]
        :param target_file: 存放处理结果的文件路径
        """
        res = []
        count = 0
        for re in self.data:
            if re in ignore_value:
                res.append(re)
                continue
            if re < min:
                res.append(min)
                count += 1
            elif re > max:
                res.append(max)
                count += 1
            else:
                res.append(re)
        with open(target_file, 'w') as file:
            for re in res:
                file.write(str(re))
                file.write('\n')
            file.close()
        print("以将结果调整为最小值为%f最大值%f其中共有%d个数据被调整，新文件路径为%s" % (min, max, count, target_file))

    def mode2(self, ignore_value, need_min_value=0.00001, target_file=r"M:\pycharm\Inversion\DataTools\data\new_mode2"):
        """
        检查结果，使结果整体偏移保证最小值为min，此处的min分两种情况，当数据的最小值大于等于0时，不进行任何处理，当数据的最小值小于0时，处理后的结果中的最小值应为need_min_value

        :param ignore_value: 忽略的值，即当ignore_value包含值value时，值value不需要进行额外处理
        :param need_min_value: 结果中的最小值应为need_min_value，当数据的最小值大于等于0时，此参数失效
        :param target_file: 存放处理结果的文件路径
        :return:
        """
        res = []
        a = self.data
        res_min = min(self.data)
        if res_min >= 0:
            res = self.data
            print("最小值已经大于等于%f,不需要偏移" % (res_min))
        else:
            need_min_value -= res_min
            for re in self.data:
                if re in ignore_value:
                    res.append(re)
                    continue
                res.append(re + need_min_value)
            print("已将结果数据偏移，最小值为%f，偏移量为%f，新文件路径为%s" % (need_min_value + res_min, -res_min, target_file))
        with open(target_file, 'w') as file:
            for re in res:
                file.write(str(re))
                file.write('\n')
            file.close()

    def smooth_cells(self, cells_ref_res, cell_refs, air_js, oldj, shape, bounds):
        """
        对结果进行平滑处理

        :param cells_ref_res: 格子的密度值(结果值-参考值)--注意这个必须已经经过 restore_res后的结果,即可以直接显示的结果
        :param cell_refs: 格子密度的参考值
        :param air_js: 空气格子的编号
        :param oldj: 不是空气的格子编号(此处为未压缩前的格子编号、且格子有无被射线穿过取决于输入)
        :param shape: 模型在x、y、z方向上的格子数目
        :param bounds: 格子密度的约束值
        :return: 平滑处理的结果
        """

        for j in range(len(cells_ref_res)):
            if cell_refs[j]==0:
                air_js.append(j+1)

        threshold_value = 0.1
        # weight = 48 / 1  # 当前的格子和周围6个格子的权重--一般格子
        weight = 15 / 1  # 当前的格子和周围6个格子的权重--一般格子
        # coefficient = (3, 3, 1.5, 1.5, 6, 6)
        coefficient = (1, 1, 1, 1, 1, 1)
        oldj = set(oldj)  # 不是空气的格子
        air_js = set(air_js)  # 空气
        all_j = set([j + 1 for j in range(shape[0] * shape[1] * shape[2])])  # 所有格子
        all_j -= air_js  # 除去空气的所有格子

        all_j_middle = set()
        un_smooth_j = []
        for j in all_j:
            if bounds[j - 1][1] - bounds[j - 1][0] < threshold_value * 2:
                # 保护外皮
                un_smooth_j.append(j)
                continue
            all_j_middle.add(j)

        print("%d个格子不需要处理" % (len(all_j) - len(all_j_middle)))
        all_j = all_j_middle
        for j in all_j:
            # 纠正
            if abs(cells_ref_res[j - 1]) < 0.001:
                # if cells_ref_res[j - 1] + cell_refs[j - 1] < 0.001:
                cells_ref_res[j - 1] = 0
        import tqdm
        # oldj_size=0
        count = 0
        for i in range(3):
            count += 1

            middle = list(all_j)
            if count % 2 == 0:
                middle.reverse()
            for j in tqdm.tqdm(middle):
                if j - 1 not in oldj:
                    continue
                x, y, z = getxyz_from_shape(shape, j)
                j -= 1
                # if j in air_js:
                ##不平滑空气
                # continue
                # if abs(cells_ref_res[j]) <= threshold_value:
                # 忽略掉过小的差异
                # cells_ref_res[j] = 0  # ???
                # continue
                # if j  in oldj:
                # 当前的格子属于一般格子(有射线穿过,并且不是空气)
                values = self._get_all_neighbor_node_values(cells_ref_res, shape, oldj, (x, y, z), j, threshold_value,
                                                            coefficient)
                cells_ref_res[j] = (cells_ref_res[j] * weight + sum(values[1])) / (weight + values[0])
                # else:
                #     # pass
                #     values = self._get_all_neighbor_node_values(cells_ref_res, shape, oldj, (x, y, z),j,2)
                #     count=0
                #     sum_values=0
                #     cells_ref_res[j] = (cells_ref_res[j] * 0.1 + 12*sum(values[1])) / (0.1 + 12*values[0])
                #     oldj.add(j + 1)
        data = []
        for i in range(len(cells_ref_res)):
            data.append(cells_ref_res[i] + cell_refs[i])

        un_raycells = all_j - set([j + 1 for j in oldj])
        un_raycells_middle = []
        # 填充没有射线穿过的格子(快速版)
        for i in range(20):
            un_raycells_middle = []
            for j in un_raycells:
                x, y, z = getxyz_from_shape(shape, j)
                j -= 1
                values = self._get_all_neighbor_node_values(data, shape, None, (x, y, z), j, 100)
                count = 0
                sum_value = 0
                for value in values[1]:
                    if value > 0 and value <= bounds[j][1]:
                        count += 1
                        sum_value += value
                if count > 0:
                    data[j] = (sum_value + cell_refs[j] * count * 0.02) / (count * 1.02)
                else:
                    un_raycells_middle.append(j + 1)
                # data[j]=10
            if len(un_raycells_middle) == len(un_raycells):
                break
        un_raycells = list(all_j - set([j + 1 for j in oldj]))
        # 对新的填充数据进行平滑处理
        for i in range(1):
            un_raycells.reverse()
            for j in un_raycells:
                x, y, z = getxyz_from_shape(shape, j)
                j -= 1
                values = self._get_all_neighbor_node_values(data, shape, None, (x, y, z), j, 100)
                count = 0
                sum_value = 0
                for value in values[1]:
                    if value > 0 and value <= bounds[j][1]:
                        count += 1
                        sum_value += value
                if count > 0:
                    data[j] = (sum_value + cell_refs[j] * count * 0.02) / (count * 1.02)

        for j in un_smooth_j:
            data[j - 1] = cell_refs[j - 1]
        for j in un_raycells_middle:
            data[j - 1] = cell_refs[j - 1]
        ###针对于不是空气密度还低
        # 检查是不是应该是城墙表面:空气是非0的定值 非空气也不是0,so 0是一个特殊的群体
        wall_need_fix = True
        if wall_need_fix:
            for i in range(len(data)):
                if data[i] <= 0:
                    if cell_refs[i] > 0:
                        data[i] = cell_refs[i]
                    else:
                        if i + 1 in air_js:
                            continue
                        x, y, z = getxyz_from_shape(shape, i + 1)
                        flag = []
                        if x + 1 <= shape[0]:
                            flag.append(getj_from_xyz(shape, (x + 1, y, z)) in air_js)
                        if x - 1 > 0:
                            flag.append(getj_from_xyz(shape, (x - 1, y, z)) in air_js)
                        if y + 1 <= shape[1]:
                            flag.append(getj_from_xyz(shape, (x, y + 1, z)) in air_js)
                        if y - 1 > 0:
                            flag.append(getj_from_xyz(shape, (x, y - 1, z)) in air_js)
                        if z - 1 > 0:
                            flag.append(getj_from_xyz(shape, (x, y, z - 1)) in air_js)  # 已经通过count简化
                        if z + 1 <= shape[2]:
                            flag.append(getj_from_xyz(shape, (x, y, z + 1)) in air_js)
                        if len(flag) == 6 and flag.count(True) < 6:
                            data[i] = 2.4
        return data

    def _get_all_neighbor_node_values(self, data, shape, oldjs, xyz, j, threshold_value,
                                      coefficient=(1, 1, 1, 1, 1, 1)):
        """
        得到所有的邻居格子的编号

        :param data: 格子的密度差值(=密度结果值-密度参考值)
        :param shape: 模型在x、y、z方向上的格子数目
        :param oldjs: 未压缩前的格子编号
        :param xyz: 格子的离散坐标
        :param coefficient: 邻居格子的权重，顺序依次为左、右、前、后、下、上
        :return: 满足条件的邻居格子的总权重、所有邻居格子的密度差值(不满足条件的赋值为0)
        """

        x, y, z = xyz
        nodes = []
        count = 0
        nodes.append(self._get_neighbor_node_value(shape, -1, oldjs, (x - 1, y, z)))
        nodes.append(self._get_neighbor_node_value(shape, -1, oldjs, (x + 1, y, z)))
        nodes.append(self._get_neighbor_node_value(shape, -1, oldjs, (x, y - 1, z)))
        nodes.append(self._get_neighbor_node_value(shape, -1, oldjs, (x, y + 1, z)))
        nodes.append(self._get_neighbor_node_value(shape, -1, oldjs, (x, y, z - 1)))
        nodes.append(self._get_neighbor_node_value(shape, -1, oldjs, (x, y, z + 1)))
        for i in range(len(nodes)):
            if nodes[i] == -1:
                nodes[i] = 0
                continue
            # count += 1
            nodes[i] = data[nodes[i]] * coefficient[i]
            # if abs(nodes[i]-data[j])>threshold_value:
            #     nodes[i]=0
            #     continue
            count += coefficient[i]
        return (count, nodes)

    # neighbor_node_cache=None
    def _get_neighbor_node_value(self, shape, default_j, oldjs, xyz):
        """
        得到单个格子的邻居格子的编号，邻居格子指格子的上、下、左、右、前、后共六个格子

        :param shape: 模型在x、y、z方向上的格子数目
        :param default_j: 邻居格子编号的默认值，如果邻居格子的离散坐标和编号不满足一定的条件，那么返回默认值
        :param oldjs: 非空气格子的编号
        :param xyz: 邻居格子的离散坐标
        :return: 邻居格子的编号
        """
        # if self.neighbor_node_cache is None:
        #     self.neighbor_node_cache=numpy.zeros(shape,dtype=int)
        x, y, z = xyz
        if x < 1 or y < 1 or z < 1 or x > shape[0] or y > shape[1] or z > shape[2]:
            return default_j
        # cache_j=self.neighbor_node_cache[x-1][y-1][z-1]
        # if cache_j !=0:
        #     return cache_j

        res = getj_from_xyz(shape, xyz) - 1
        if oldjs is not None and res not in oldjs:
            # self.neighbor_node_cache[x - 1][y - 1][z - 1]=default_j
            return default_j
        # self.neighbor_node_cache[x - 1][y - 1][z - 1]=res
        return res


class Merge_res:
    """
    此类未完成

    合并结果

    """
    def __init__(self):
        pass

    def Merge_res_from_file(self, res_file1, res_file2):
        """
        此方法未完成

        合并从文件得到的数据

        :param res_file1: 第一个数据文件路径
        :param res_file2: 第二个数据文件路径
        """
        file1 = open(res_file1)
        file2 = open(res_file2)


if __name__ == '__main__':
    # tool = new_res_tools()
    # tool.mode1(min=0.01, ignore_value=[0, -0.1234], max=4.0001)
    # tool.mode2(ignore_value=[0, -0.1234])
    import DataFactory.DataManage

    # datamanager = DataFactory.DataManage.DataManager()
    datamanager=DataFactory.DataManage.DataManager.get_instance()
    tool = new_res_tools(res_file=r"M:\pycharm\Inversion\data\output\rel\ref_ps")
    refs = datamanager.refs_tool.get_data()
    res_refs = tool.data
    others = datamanager.Make_A()
    tool = new_res_tools(
        tool.smooth_cells(res_refs, refs, datamanager.get_unneed_j(), others[-1].keys(), datamanager.mesh.get_shape(),
                          datamanager.bonds_tool.get_bonds_min_max()))
    tool.mode0(target_file=r"M:\pycharm\Inversion\data\output\rel\res_smooth")
