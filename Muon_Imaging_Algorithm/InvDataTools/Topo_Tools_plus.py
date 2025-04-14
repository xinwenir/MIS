# author:高金磊
# datetime:2022/8/10 14:45

"""
通过正演的射线与topo的交点信息来标记各个部分(城墙/空气)
"""
from ast import If
from logging import warning
from InvDataTools.Jxyz_Tools import getxyz_from_shape, getj_from_xyz
from InvDataTools.MeshTools import MeshTools
from InvDataTools.calculate.Calcu_sensitivity import Calcsensitivity


class Topo_Tools_plus:
    """
    需要ryj修改计算GIJ的程序拆分模块来计算一个射线的路径J
    """

    def __init__(self, points_file, obs_file, mesh_file):
        """
        :param data_file: lgr定义的文件格式
        :param temp_obs_file:
        :param temp_points_info:
        """
        points = open(points_file, 'r').readlines()
        mesh_tool = MeshTools(mesh_file)
        points_info = []
        for line in points:
            elements = line.split()
            if len(elements) % 3 != 0:
                print("数据不完整或者有误,请检查此行数据的位数")
                print(line)
            point_info = []
            point = []
            for i in range(len(elements)):
                if len(point) == 3:
                    point_info.append(point)
                    point = []
                point.append(float(elements[i]))
            if len(point) != 0:
                point_info.append(point)
            points_info.append(point_info)
        # 调用求解GIJ的算法
        tool = Calcsensitivity(mesh_tool,loger=None,obs_file=obs_file)
        import os
        if not os.path.exists(r"Temp/Topo_Tools_plus"):
            os.makedirs(r"Temp/Topo_Tools_plus")
        tool.calc_all_rays_from_obs_file(r"Temp/Topo_Tools_plus/Gij", r"Temp/Topo_Tools_plus/points_info")

        # 加载射线穿过的j编号
        points_info_discrete_datas = open(r"Temp/Topo_Tools_plus/points_info", 'r').readlines()
        points_info_discretes = []
        for points_info_discrete_data in points_info_discrete_datas:
            middle = []
            for points_info_discrete_datum in points_info_discrete_data.split():
                middle.append(int(points_info_discrete_datum))
            middle.pop(0)
            middle.sort()
            points_info_discretes.append(middle)
        self.mesh_file = mesh_file
        self.points_info_discretes = points_info_discretes
        self.points_info = points_info

    def make_topo_j_by_points_info(self):
        """根据交点信息将格子分为两类

        Returns:
            list: 分类后结果list[2][]
        """
        # j_coordinates={}
        mesh_tool = MeshTools(self.mesh_file)
        # 离散交点坐标
        points_j = []
        for point_infos in self.points_info:
            middle = []
            flag = True
            for point_info in point_infos:
                xyz = mesh_tool.discretize_Physical_coordinates(point_info[0], point_info[1], point_info[2], False)
                j = getj_from_xyz(mesh_tool.get_shape(), xyz)
                if j in middle:
                    flag = False
                if j not in middle:
                    middle.append(j)
                else:
                    middle.remove(j)
            # if len(middle) % 2 != 0 or not flag:
            #     print(middle)
            #     middle = []
            points_j.append(middle)

        all_j = set()
        ##通过射线编号找区间---如果有必要需要放宽搜索范围
        for i in range(len(self.points_info_discretes)):
            points_info_discrete = self.points_info_discretes[i]
            points = points_j[i]

            # if len(points)%2!=0 or len(points)==0:
            #     continue
            for p in range(len(points)):
                point=points[p]
                if point not in points_info_discrete:
                    #修正模式
                    new_point=-10000
                    for j in points_info_discrete:
                        if abs(j-point)<abs(point-new_point):
                            new_point=j
                    points[p]=new_point
                all_j.add(points[p])
        self.all_j=all_j
        ######通过射线的格子编号和交点格子编号二分区域
        area_js = [set(), set()]
        for i in range(len(self.points_info_discretes)):
            flag = 1
            points_info_discrete = self.points_info_discretes[i]
            points = points_j[i]
            # if len(points) % 2 != 0:
            #     continue
            for j in points_info_discrete:
                if j in points:
                    flag += 1
                area_js[flag % len(area_js)].add(j)



        # area_js[1] -= area_js[0]
        self.area_js = area_js
        return area_js
    def update_refs_bounds(self,mesh_file,refs_file,bounds_file,refs_value,bounds_value:tuple,default_refs_value,default_bounds_value:tuple):
        """根据当前的分类信息更新refs和bounds文件,
        !如果这两个文件不存在将自动新建(包括目录)
        Args:
        mesh_file (String): mesh文件路径 必须存在
        refs_file (String): refs文件路径 存在将进行原地修改
        bounds_file (String): refs文件路径 存在将进行原地修改
        refs_value (float): 该类别refs的值
        bounds_value (tuple): 该类别bounds的值
        default_refs_value (int): 新建文件时使用的默认值. Defaults to 0.
        default_bounds_value (tuple): 新建文件时使用的默认值. Defaults to (-0.001,0.001).

        Raises:
            Exception: 数据类型错误
            Exception: default_bounds_value不是元组
            Exception: mesh文件路径不正确
            Exception: 数据数量不匹配
        """
        #检查数据是否合法
        if type(bounds_value) is not tuple:
            print(bounds_value,"不是元组")
            raise Exception("数据类型错误")
        if type(default_bounds_value) is not tuple:
            print(default_bounds_value,"不是元组")
            raise Exception("数据类型错误")
        #判断文件是否存在
        import os
        if not os.path.exists(mesh_file):
            raise Exception("mesh文件路径不正确:%s"%(mesh_file))
        mesh_tool=MeshTools(mesh_file)
        if not os.path.exists(refs_file):
            print("文件不存在将会新建文件",refs_file)
            if not os.path.exists(os.path.split(refs_file)[0]):#检查文件夹
                os.makedirs(os.path.split(refs_file)[0])
            file=open(refs_file,'w')
            for i in range(mesh_tool.cells_count()):
                file.write(str(default_refs_value))
                file.write("\n")
            file.close()
        if not os.path.exists(bounds_file):
            print("文件不存在将会新建文件",bounds_file)
            if not os.path.exists(os.path.split(bounds_file)[0]):#检查文件夹
                os.makedirs(os.path.split(bounds_file)[0])
            file=open(bounds_file,'w')
            for i in range(mesh_tool.cells_count()):
                file.write(str(default_bounds_value[0]))
                file.write(" ")
                file.write(str(default_bounds_value[1]))
                file.write("\n")
            file.close()
        #加载数据
        refs_data=open(refs_file,'r').readlines()
        bounds_data=open(bounds_file,'r').readlines()
        if len(refs_data) !=len(bounds_data) or len(refs_data) != mesh_tool.cells_count():
            print("refs文件体素数量:%s \n bounds文件体素数量:%s \n mesh文件要求数量:%s"%(str(len(refs_data)),str(len(bounds_data)),str(mesh_tool.cells_count())))
            raise Exception("数据数量不匹配")        
        #对处于交点1和交点2之间地格子进行赋值
        for j in self.area_js[0]:
            refs_data[j-1]=str(refs_value)+"\n"
            bounds_data[j-1]=str(bounds_value[0])+" "+str(bounds_value[1])+"\n"
        #回写数据
        file=open(refs_file,'w')
        for s in refs_data:
            file.write(s)
            # file.write("\n")
        file.close()
        file=open(bounds_file,"w")
        for s in bounds_data:
            file.write(s)
            # file.write('\n')
        file.close()
        
        
        
    def show_distribution(self, res_file):
        import warnings
        warnings.warn("此方法已经被弃用")
        """停止使用

        Args:
            res_file (_type_): _description_
        """
        res = open(res_file, 'w')
        mesh_tool = MeshTools(self.mesh_file)
        for j in range(mesh_tool.cells_count()):
            j = j + 1
            flag = True
            for i in range(len(self.area_js)):
                if j in self.all_j:
                    res.write('1.5')
                    flag = False
                    break
                if j in self.area_js[i]:
                    res.write(str(i+1))
                    flag = False
                    break
            if flag:
                res.write('0')

            res.write('\n')


# if __name__ == '__main__':
#     tool = Topo_Tools_plus(
#         r"E:\vscode\Muon_Imaging_Algorithm\data\Looking_for_area_through_the_ray_intersection\points_bottom_mesh_SF.txt",
#         r"E:\vscode\Muon_Imaging_Algorithm\data\Looking_for_area_through_the_ray_intersection\points_bottom_mesh_SF.dat",
#         r"E:\vscode\Muon_Imaging_Algorithm\data\Input\real_data\17_58MaMian.msh"
#     )
#     tool.make_topo_j_by_points_info()
#     tool.show_distribution(
#         r"E:\vscode\Muon_Imaging_Algorithm\data\Looking_for_area_through_the_ray_intersection\points_bottom_mesh_SF_res.den")

    # tool = Topo_Tools_plus(
    #         r"E:\vscode\Muon_Imaging_Algorithm\data\Looking_for_area_through_the_ray_intersection\test_points",
    #         r"E:\vscode\Muon_Imaging_Algorithm\data\Looking_for_area_through_the_ray_intersection\test_obs",
    #         r"E:\vscode\Muon_Imaging_Algorithm\data\Input\real_data\17_58MaMian.msh"
    #         )
    # tool.make_topo_j_by_points_info()
    # tool.show_distribution(
    #         r"E:\vscode\Muon_Imaging_Algorithm\data\Looking_for_area_through_the_ray_intersection\58mamian_mesh_D1_res.den")
