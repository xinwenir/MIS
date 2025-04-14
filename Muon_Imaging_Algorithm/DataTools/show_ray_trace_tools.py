# author:高金磊
# datetime:2022/2/25 17:59
import warnings
from math import floor, ceil

import numpy
import tqdm

import SysTools
from DataTools.Air_j import Air_j
from DataTools.Bonds_tools import Bonds_tool
from DataTools.Jxyz_Tools import getj_from_xyz, getxyz_from_shape
from DataTools.MeshTools import MeshTools
from DataTools.ref_tools import Ref_tools
from SysTools.MyTools import myPrint
from DataTools.obs_tools import obs_tools
from DataTools.res_tools import res_tools
from SysTools.tools import Desc_Text_alignment_tools, Points_In_Coners_Tools


class show_ray_trace_tools:
    """
    展示射线轨迹

    """
    def __init__(self, shape):
        self.data = [0.0] * (shape[0] * shape[1] * shape[2])
        self.shape = shape

    def mark_ray(self, ids: list, file_res,file_ij, show=None, group_method: list = None, id_values: dict = None):
        """
        根据所给的射线编号来标记格子,支持使用区间编号,指定编号等方式
        
        :param ids: 射线编号
        :param show: 是否显示数据
        :param group_method: 编号区间[a,b,c]  (a,b]标记为1 (b,c]标记为2 其余为0  区间数量最好小于10个 def:一个射线一个编号
        :param id_values: 固定编号,没有为-1 优先级最高
        :param file_ij: ij文件路径(读取)
        :param file_res: 存放结果的文件路径(写入)
        """
        if show is None:
            if len(ids) > 20 and group_method is None:
                print("数据太多默认不显示")
                show = False
            else:
                show = True
        id_value = {}
        id_j = {}
        if show:
            print("将按照以下规则标记:")
        if not id_values:
            # 无固定编号
            if not group_method:
                # 无区间编号
                for i in range(len(ids)):
                    id_value[ids[i]] = i + 1
                    id_j[ids[i]] = []
                    if show:
                        print("射线%s使用值%s" % (str(ids[i]), str(i + 1)))
            else:
                for id in ids:
                    value = 0
                    for j in range(len(group_method)):
                        if id < group_method[j]:
                            value = j+1
                            break
                    id_value[id] = value
                    id_j[id] = []
                    if show:
                        print("射线%s使用值%s" % (str(id), str(value)))
        else:
            for id in ids:
                if id in id_values:
                    id_value[id] = id_values.get(id, -1)
                else:
                    id_value[id] = -1
                id_j[id] = []
                if show:
                    print("射线%s使用值%s" % (str(id), str(id_value[id])))
        ij_data = open(file_ij, mode='r')
        ij_data.readline()
        res_data = open(file_res, mode='w')
        ijs = ij_data.readlines()
        for ij in tqdm.tqdm(ijs, Desc_Text_alignment_tools("正在对部分射线标记")):
            if not ij:
                break
            middle = ij.split()
            if int(middle[0]) in id_j.keys():
                if int(middle[1]) > len(self.data):
                    myPrint.myPrint_Err("出现越界值:", middle[1])
                    continue
                if id_values is None:
                    # 不使用固定编号,尽可能显示所有射线的编号
                    self.data[int(middle[1]) - 1] /= 10
                    self.data[int(middle[1]) - 1] += id_value[int(middle[0])]
                else:
                    if self.data[int(middle[1]) - 1] == 0:
                        # 有固定编号,显示最小值
                        self.data[int(middle[1]) - 1] = id_value[int(middle[0])]
                    else:
                        if self.data[int(middle[1]) - 1] > id_value[int(middle[0])]:
                            self.data[int(middle[1]) - 1] = id_value[int(middle[0])]
                id_j[int(middle[0])].append(int(middle[1]))
        for datum in tqdm.tqdm(self.data, desc=Desc_Text_alignment_tools("正在写入中间结果")):
            res_data.write(str(float(datum)))
            res_data.write("\n")
        res_data.close()

        # 由于代码更新,此处不需要再对结果进行转化
        # res_tools().Conversion_2(self.shape, file_res=file_res, file_out=file_res)
        if show:
            for id in ids:
                print("射线%s穿过的格子编号" % (id), id_j[id])

    def mark_ray_all(self, file_ij,
                     file_res=r"M:\pycharm\Inversion\DataTools\data\ray_trace"):
        """
        显示所有格子穿过射线的条数，其中格子的密度值为射线穿过的条数

        :param file_ij: ij文件路径(读取)
        :param file_res: 存放结果的文件路径(写入)
        """
        print("将显示所有格子穿过射线的条数,其中格子的密度值为射线穿过的个数")
        ij_data = open(file_ij, mode='r')
        res_data = open(file_res, mode='w')
        ijs = ij_data.readlines()
        for ij in tqdm.tqdm(ijs, desc=Desc_Text_alignment_tools("正在统计所有格子穿过射线的条数")):
            if not ij:
                break
            middle = ij.split()
            self.data[int(middle[1]) - 1] += 1

        for datum in tqdm.tqdm(self.data, desc=Desc_Text_alignment_tools("正在写入中间结果")):
            res_data.write(str(float(datum)))
            res_data.write("\n")
        res_data.close()


        #由于代码更新,此处不需要再对结果进行转化
        # res_tools().Conversion_2(self.shape, file_res=file_res, file_out=file_res)

    def find_ray_by_cells(self, ids,file_ij, show=None):
        """
        找出穿过ids包含的所有格子的射线

        :param ids: 格子的id
        :param show: 是否打印详情
        :param file_ij: ij文件路径(读取)
        :return: 射线的id
        """
        if show is None:
            if len(ids) > 20:
                show = False
            else:
                show = True
        if type(ids) is not set:
            ids = set(ids)
        ij_data = open(file_ij, mode='r')
        ijs = ij_data.readlines()
        cell_rays = {}
        res = set()
        for id in ids:
            cell_rays[id] = []
        for ij in tqdm.tqdm(ijs, desc=SysTools.tools.Desc_Text_alignment_tools("正在寻找穿过格子射线的id")):
            if not ij:
                continue
            middle = ij.split()
            if int(middle[1]) in ids:
                cell_rays[int(middle[1])].append(int(middle[0]))
                res.add(int(middle[0]))
        if show:
            for id in ids:
                print("穿过格子", id, "的射线有", cell_rays[id])
        return res

    def get_j_from_scope_general(self,mesh_tool,corners):
        """找出四面体的区域内所有格子的编号"""
        pass

    def get_j_from_scope(self, corners, xcmin, ycmin, zcmin, x_step, y_step, z_step,
                         jxyz_file, show=None) -> set:
        """
        找出四面体的区域内所有格子的编号
        要求格子是均匀划分的
        :param corners: 8*3,表示每个顶点的xyz
        :param xcmin: x起始坐标
        :param ycmin: y起始坐标
        :param zcmin: z起始坐标
        :param x_step: 格子x方向的宽度
        :param y_step: 格子y方向的宽度
        :param z_step: 格子z方向的宽度
        :param jxyz_file: jxyz文件路径(读取)
        :param show: 是否显示详情
        :return: 所有符合的格子
        """
        warnings.warn("请选用支持不均匀划分的版本", DeprecationWarning)
        tool = Points_In_Coners_Tools(coners=corners)
        self.jxyz_file = jxyz_file
        xl, yl, zl = [], [], []
        j_set = set()
        for i in corners:
            xl.append(i[0])
            yl.append(i[1])
            zl.append(i[2])
        xmin, xmax = floor((min(xl) - xcmin) / x_step), ceil((max(xl) - xcmin) / x_step)
        ymin, ymax = floor((min(yl) - ycmin) / y_step), ceil((max(yl) - ycmin) / y_step)
        zmin, zmax = floor((min(zl) - zcmin) / z_step), ceil((max(zl) - zcmin) / z_step)
        middle = []
        with open(jxyz_file, 'r') as f:
            for i_line in tqdm.tqdm(f.readlines(), desc=SysTools.tools.Desc_Text_alignment_tools("正在筛选在区间内的格子")):
                j, x, y, z = [int(ii) for ii in i_line.strip().split()]
                if xmin <= x <= xmax and ymin <= y <= ymax and zmin <= z <= zmax:
                    points = [[xcmin + (x - 1) * x_step, ycmin + (y - 1) * y_step, zcmin + (z - 1) * z_step],
                              [xcmin + x * x_step, ycmin + (y - 1) * y_step, zcmin + (z - 1) * z_step],
                              [xcmin + (x - 1) * x_step, ycmin + (y - 1) * y_step, zcmin + z * z_step],
                              [xcmin + x * x_step, ycmin + (y - 1) * y_step, zcmin + z * z_step],
                              [xcmin + (x - 1) * x_step, ycmin + y * y_step, zcmin + (z - 1) * z_step],
                              [xcmin + x * x_step, ycmin + y * y_step, zcmin + (z - 1) * z_step],
                              [xcmin + (x - 1) * x_step, ycmin + y * y_step, zcmin + z * z_step],
                              [xcmin + x * x_step, ycmin + y * y_step, zcmin + z * z_step]]
                    res = tool.is_in(points)

                    for i in res:
                        if i[0]:
                            j_set.add(j)
                            middle.append(res)
                            break
        if show is None:
            if len(middle) > 20:
                print("数据太多,默认不打印")
                show = False
            else:
                show = True
        if show:
            for m in middle:
                print(m)
        return j_set

    def get_j_from_scope_mesh(self, corners, mesh:MeshTools,
                         jxyz_file=None, show=None) -> set:
        """
        找出四面体的区域内所有格子的编号
        格子可以是不均匀的
        :param corners: 8*3,表示每个顶点的xyz
        :param mesh : Mesh_Tools对象
        :param jxyz_file: jxyz文件路径(读取)
        :param show: 是否显示详情
        :return: 所有符合的格子
        """
        tool = Points_In_Coners_Tools(coners=corners)
        self.jxyz_file = jxyz_file
        j_set = set()

        middle = []
        if jxyz_file is not None:
            with open(jxyz_file, 'r') as f:
                for i_line in tqdm.tqdm(f.readlines(), desc=SysTools.tools.Desc_Text_alignment_tools("正在筛选在区间内的格子")):
                    j, x, y, z = [int(ii) for ii in i_line.strip().split()]
                        #体素的8个点有一个在就认为在
                        # points = [[xcmin + (x - 1) * x_step, ycmin + (y - 1) * y_step, zcmin + (z - 1) * z_step],
                        #           [xcmin + x * x_step, ycmin + (y - 1) * y_step, zcmin + (z - 1) * z_step],
                        #           [xcmin + (x - 1) * x_step, ycmin + (y - 1) * y_step, zcmin + z * z_step],
                        #           [xcmin + x * x_step, ycmin + (y - 1) * y_step, zcmin + z * z_step],
                        #           [xcmin + (x - 1) * x_step, ycmin + y * y_step, zcmin + (z - 1) * z_step],
                        #           [xcmin + x * x_step, ycmin + y * y_step, zcmin + (z - 1) * z_step],
                        #           [xcmin + (x - 1) * x_step, ycmin + y * y_step, zcmin + z * z_step],
                        #           [xcmin + x * x_step, ycmin + y * y_step, zcmin + z * z_step]]
                    #使用中点来代表这个格子
                    points=[[(mesh.get_coordinates_form_xyz(x,y,z,strategy=(0,0,0)))]]
                    res = tool.is_in(points)
                    for i in res:
                        if i[0]:
                            j_set.add(j)
                            middle.append(res)
                            break
        else:
            for j in tqdm.trange(mesh.cells_count(), desc=SysTools.tools.Desc_Text_alignment_tools("正在筛选在区间内的格子")):
                j+=1
                x,y,z=getxyz_from_shape(shape=mesh.get_shape(),j=j)
                points = [[(mesh.get_coordinates_form_xyz(x, y, z, strategy=(0, 0, 0)))]]
                res = tool.is_in(points)
                for i in res:
                    if i[0]:
                        j_set.add(j)
                        middle.append(res)
                        break
        if show is None:
            if len(middle) > 20:
                print("数据太多,默认不打印")
                show = False
            else:
                show = True
        if show:
            for m in middle:
                print(m)
        return j_set

    def draw_points(self, cells: set, res_file, jxyz_file):
        """
        将cells包含的格子的密度值置为1

        :param cells: 要显示的格子编号
        :param res_file: 存放结果的文件路径(写入)
        :param jxyz_file: jxyz文件路径(读取)
        """
        res = numpy.zeros((self.shape[0] * self.shape[1] * self.shape[2]))
        for cell in cells:
            res[int(cell) - 1] = 1
        file = open(res_file, "w")
        for re in res:
            file.write(str(re))
            file.write('\n')
        file.close()
        res_tools().Conversion_2(self.shape, file_res=res_file, file_out=res_file, file_xyz=jxyz_file)
def show_points_3D(data):
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    x = []
    y = []
    z = []
    for datum in data:
        x.append(datum[0])
        y.append(datum[1])
        z.append(datum[2])

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, c='r', marker='o')
    for i in range(len(x)):
        ax.text(x[i],y[i],z[i],str(i+1))
    plt.show()

if __name__ == '__main__':
    model=3
    if model==1:
        tool = show_ray_trace_tools([140, 72, 42])
        # target=[i for i in range(26612,26621)]+[i  for i in range(1240,1249) ]
        # target=[i for i in range(6714)]
        # target=[100,300,555]
        # target=[i for i in range(25770,40000)]
        # target=[i for i in range(26483,40000)]
        # tool.mark_ray(target)
        # tool.mark_ray([1,2,3])
        # tool.mark_ray_all(r"M:\pycharm\Inversion\data\Temp\Gij")
        # tool.mark_ray_all(file_res=r"M:\pycharm\Inversion\DataTools\data\ray_trace")
        # tool.mark_ray([33450, 33451, 33452, 31151, 31152, 31153, 31154, 31155, 31156, 31157, 31158, 31159, 31160, 31161, 31162, 31163, 31164, 31165, 31166, 31191, 31192, 31193, 31194, 31195, 31196, 31197, 31198, 31199, 31200, 31201, 31202, 31203, 31204, 31205, 31206, 31231, 31232, 31233, 31234, 31235, 31236, 31237, 31238, 31239, 31240, 31241, 31242, 31243, 31244, 31245, 31246, 31247, 33482, 33483, 33484, 31272, 31273, 31274, 31275, 31276, 31277, 31278, 31279, 31280, 31281, 31282, 31283, 31284, 31285, 31286, 31287, 31288, 33481, 33496, 31314, 31315, 31316, 31317, 31318, 31319, 31320, 31321, 31322, 31323, 31324, 31325, 31326, 31327, 31328, 31329, 33498, 31358, 31359, 31360, 31361, 31362, 31363, 31364, 31365, 31366, 31367, 31368, 31369, 31370, 31371, 33425, 33426, 33427, 33428, 33429, 33430, 33431, 33432, 33433, 33434, 33435, 33436, 33437, 33438, 33439, 33440, 33441, 33442, 33443, 33444, 33445, 33446, 33447, 33448, 33449, 31402, 31403, 31404, 31405, 31406, 31407, 31408, 31409, 31410, 31411, 31412, 31413, 33453, 33454, 33455, 33456, 33457, 33458, 33459, 33460, 33461, 33462, 33463, 33464, 33465, 33466, 33467, 33468, 33469, 33470, 33472, 33473, 33474, 33475, 33476, 33477, 33478, 33479, 33488, 33480, 33490, 33491, 33492, 33493, 33494, 33495, 33487, 31446, 31447, 31448, 31449, 31450, 31451, 31452, 31453, 31454, 31455, 31456, 33499, 33504, 33505, 33506, 33507, 33508, 33509, 33510, 33511, 33512, 33513, 33514, 33515, 33516, 33522, 33523, 33524, 33525, 33526, 33527, 33528, 33529, 33530, 33531, 33532, 33489, 33533, 33534, 33535, 33542, 31495, 31496, 31497, 33543, 33544, 33545, 33546, 33547, 33548, 33549, 33550, 33551, 33552, 33553, 33554, 33555, 33556, 33563, 33564, 33565, 33566, 33567, 33568, 33569, 33570, 33571, 33572, 33573, 33574, 33575, 33576, 33577, 33587, 33588, 33589, 33590, 33591, 33592, 33593, 33594, 33595, 33596, 33597, 33598, 33599, 33600, 33601, 33602, 33613, 33614, 33615, 33616, 33617, 33618, 33619, 33620, 33621, 33622, 33623, 33624, 33625, 33626, 33627, 33628, 33640, 33641, 33642, 33643, 33644, 33645, 33646, 33647, 33648, 33649, 33650, 33651, 33652, 33653, 33654, 33655, 33669, 33670, 33671, 33672, 33673, 33674, 33675, 33676, 33677, 33678, 33679, 33680, 33681, 33682, 33683, 33684, 33685, 33686, 33687, 33702, 33703, 33704, 33705, 33706, 33707, 33708, 33709, 33710, 33711, 33712, 33713, 33714, 33715, 33716, 33717, 33718, 33719, 33720, 33497, 33735, 33736, 33737, 33738, 25547, 33739, 33740, 33741, 33742, 33743, 33744, 33745, 33746, 33747, 33748, 33749, 33750, 33751, 33752, 33753, 33754, 25574, 25575, 33771, 33772, 33773, 33774, 33775, 33776, 33777, 33778, 33779, 33780, 33781, 33782, 33783, 33784, 33785, 33786, 33787, 33788, 33789, 33790, 33791, 25605, 25606, 33811, 33812, 33813, 33814, 33815, 33816, 33817, 33818, 33819, 33820, 33821, 33822, 33823, 33824, 33825, 33826, 33827, 33828, 25637, 33829, 33830, 33852, 33853, 33854, 33855, 33856, 33857, 33858, 33859, 33860, 33861, 25670, 33862, 33863, 33864, 33865, 33866, 33867, 33868, 33895, 33896, 33897, 33898, 33899, 33900, 33901, 33902, 33903, 33904, 33905, 33906, 33907, 33908, 33941, 33942, 33943, 33944, 33945, 33946, 33947, 33948, 33949, 33950, 33951, 33952, 33988, 33989, 33990, 33991, 33992, 33993, 33994, 33995, 33996, 34034, 34035, 34036, 34037, 34038, 34039, 34040, 34083, 34084, 34085, 34086, 34087, 34133, 34134, 34135])

        id_values = {}
        obs_tool = obs_tools(r"M:\pycharm\Inversion\data\Input\real_data\27_obs.dat")
        # d_err = obs_tool.get_d_relative_error_form_obs()
        # for i in range(len(d_err)):
        #     id_values[i] = d_err[i]

        # corners = [[4.34, -13.5, 12.7463], [2.5025, -13.49, 12.7596], [2.5025, -13.6682, 10.7596],
        #            [4.1874, -13.6682, 10.7724],
        #            [4.37, -11, 12.7582], [2.5319, -10.98, 12.7715], [2.5319, -11.0258, 10.7715],
        #            [4.2111, -11.0258, 10.8021]]

        # corners=[[12,8,7],[9,8,7],[9,8,5],
        #            [12,8,5],
        #            [12,12,7],[9,12,7],[9,12,5],
        #            [12,12,5]
        #            ]
        x_range=[15,22.5]
        y_range=[-12.7,-10.5]
        z_range=[0.5,3.5]
        corners = [[max(x_range), min(y_range), max(z_range)], [min(x_range), min(y_range), max(z_range)], [min(x_range), min(y_range), min(z_range)],
                   [max(x_range), min(y_range), min(z_range)],
                   [max(x_range), max(y_range), max(z_range)], [min(x_range), max(y_range), max(z_range)], [min(x_range), max(y_range), min(z_range)],
                   [max(x_range), max(y_range), min(z_range)]
                   ]

        # for corner in corners:
        #     corner[2] -= 2
        # if corner[2]<12:
        #     corner[2]-=2
        xcmin = -25
        ycmin = -32
        zcmin = 0.1
        x_step = 0.5
        y_step = 0.5
        z_step = 0.3
        # tool.get_j_from_scope()
        cells = tool.get_j_from_scope(corners, xcmin, ycmin, zcmin, x_step, y_step, z_step,
                                      r"M:\pycharm\Inversion\data\Temp\jxyz", show=False)
        myPrint.myPrint_Hint("区域格子的ids:", cells)

        tool.draw_points(cells, res_file=r"M:\pycharm\Inversion\data\output\rel\scop_cells",
                         jxyz_file=r"M:\pycharm\Inversion\data\Temp\jxyz")
        rays = tool.find_ray_by_cells(cells, show=True,file_ij=r"M:\pycharm\Inversion\data\Temp\Gij")

        all_rays=obs_tool.get_data()
        # middle=[]
        # for i in rays:
        #     ray_data=all_rays[i-1]
        #     if abs(float(ray_data[7])-float(ray_data[6])*1.89)/float(ray_data[8])>=2:
        #         middle.append(i)
        # rays=middle


        myPrint.myPrint_Hint("所有的射线编号:", rays)
        # middle=set()
        # for ray in rays:
        #     if ray>=26483:
        #         middle.add(ray)

        # group_method=obs_tool.get_obs_id_count()
        # for i in range(1,len(group_method)):
        #     group_method[i]+=group_method[i-1]
        # tool.mark_ray(list(rays),show=True,group_method=[0,6715,13167,19158,25526,31132,37328])
        # tool.mark_ray(list(rays),show=True,group_method=[0,2591,5087,7037,9676,11769,14061])

        tool.mark_ray(list(rays), show=True, group_method=None, id_values=id_values,file_res=r"M:\pycharm\Inversion\data\output\rel\ray_way",file_ij=r"M:\pycharm\Inversion\data\Temp\Gij")
    elif model==2:
        tool=show_ray_trace_tools([82,82,82])
        obs_tool=obs_tools(r"M:\pycharm\Inversion\data\seed_study\CBC1.dat")
        detectoers=obs_tool.get_obs_id_count()
        for i in range(1,len(detectoers)):
            detectoers[i]+=detectoers[i-1]
        tool.mark_ray(file_ij=r"M:\pycharm\Inversion\data\Temp\Gij",file_res=r"M:\pycharm\Inversion\data\output\rel\ray_way",group_method=detectoers,ids=[i for i in range(0,detectoers[-1])])
    elif model==3:
        tool = show_ray_trace_tools([82, 82, 82])
        tool.mark_ray([i for i in range(8000)], show=True, group_method=None,
                      file_res=r"M:\pycharm\Inversion\data\output\rel\ray_way",
                      file_ij=r"M:\pycharm\Inversion\data\Temp\Gij")
