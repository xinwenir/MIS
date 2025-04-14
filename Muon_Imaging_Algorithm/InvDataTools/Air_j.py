# author:高金磊
# datetime:2021/12/1 14:28
import time

import tqdm

from InvDataTools.Jxyz_Tools import getj_from_xyz
from InvDataTools.MeshTools import MeshTools
from InvDataTools.topo_tools import topo_tools
from InvSysTools.MyTools.myPrint import myPrint_Err


class Air_j:
    """
    判断空气
    """
    def get_air_j_to_file(self, mesh: MeshTools, topo_file, unneed_j_file, strategy):
        """
        判断格子是否是空气，并将空气格子的编号存放在文件中

        :param mesh: mesh文件路径(读取)
        :param topo_file: topo文件路径(读取)
        :param unneed_j_file: 存放空气格子编号的文件路径(写入)
        :param strategy: 取哪个点作为格子的离散坐标
        :return: 是空气的格子编号
        """
        un_need_j = []
        count = 0
        shape = mesh.get_shape()
        topo = topo_tools(topo_file)
        j_max = shape[0] * shape[1] * shape[2]
        print("共有", j_max, "个数据要判断是否是空气")
        # for j in tqdm.trange(1,j_max+1):#速度太慢了
        #     x,y,z=getxyz_from_shape(shape,j)
        #     rx,ry,rz=mesh.get_coordinates_form_xyz(x,y,z,(1,1,1))
        #     check_info=topo.check_point_cache(rx,ry,rz)
        #     if check_info[0]:
        #         un_need_j.append(j)
        #         count += 1
        for x in tqdm.trange(1, shape[0] + 1):
            for y in range(1, shape[1] + 1):
                for z in range(1, shape[2] + 1):
                    rx, ry, rz = mesh.get_coordinates_form_xyz(x, y, z, strategy)
                    check_info = topo.check_point_cache(rx, ry, rz)
                    if check_info[0]:
                        un_need_j.append(getj_from_xyz(shape, (x, y, z)))
                        count += 1

        un_cache = j_max - topo.get_cache_count()
        if un_cache != 0:
            myPrint_Err("共有", un_cache, "次缓存没有被命中,这些将要遍历整个topo文件,请确认topo数据的精度,调整缓存宽松度")
        self.un_need_j = un_need_j
        file = open(unneed_j_file, 'w')
        for j in un_need_j:
            file.write(str(j))
            file.write('\n')
        file.close()
        print("共有", len(un_need_j), "个数据被判断为空气")
        return un_need_j

    def get_air_j_to_file_Jxyzfile_cache(self, mesh: MeshTools, topo_file, unneed_j_file, jxyz_file, strategy):
        """
        功能等同 get_air_j_to_file
        因为使用缓存可能会大幅度提升性能,减少查找xyz和j关系的时间

        :param mesh: mesh文件路径
        :param topo_file: topo文件路径
        :param unneed_j_file: 存放空气格子编号的文件路径
        :param jxyz_file: jxyz文件路径
        :param strategy: 取哪个点作为格子的离散坐标
        :return: 空气格子的编号
        """

        un_need_j = []
        shape = mesh.get_shape()
        topo = topo_tools(topo_file)
        j_max = shape[0] * shape[1] * shape[2]
        data_jxyz = []
        data_j = set()
        count = 0
        file = open(file=jxyz_file, mode="r")
        while 1:
            line = file.readline()
            if not line:
                break
            middle = [int(i) for i in line.split()]
            if len(middle) != 4:
                continue
            j, x, y, z = middle
            if not data_j.__contains__(j):
                data_j.add(j)
                data_jxyz.append([j, x, y, z])
        for i in tqdm.tqdm(data_jxyz):
            j, x, y, z = i
            rx, ry, rz = mesh.get_coordinates_form_xyz(x, y, z, strategy=strategy)
            res = topo.check_point_cache(rx, ry, rz)
            if res[0]:
                un_need_j.append(j)
                count += 1

        un_cache = j_max - topo.get_cache_count()
        if un_cache != 0:
            myPrint_Err("共有", un_cache, "次缓存没有被命中,这些将要遍历整个topo文件,请确认topo数据的精度,调整缓存宽松度")
        self.un_need_j = un_need_j
        file = open(unneed_j_file, 'w')
        for j in un_need_j:
            file.write(str(j))
            file.write('\n')
        file.close()
        print("共有", count, "个数据被判断为空气")
        return un_need_j

    def get_air_j(self, xyzcells, step_x, step_y, step_z, topo_file,
                  jxyz=r"E:\vscode\Muon_Imaging_Algorithm\dataTools\data\jxyz",
                  un_need_file=r"E:\vscode\Muon_Imaging_Algorithm\dataTools\data\unneed_j"):
        """
        弃用方案

        功能等同 get_air_j_to_file

        :param xyzcells: x、y、z方向的起始点
        :param step_x: 离散x的步长
        :param step_y: 离散y的步长
        :param step_z: 离散z的步长
        :param topo_file:
        :param jxyz: jxyz文件路径
        :param un_need_file: 存放空气格子编号的文件路径
        :return: 空气格子的编号
        """
        import warnings
        warnings.warn("这个方法已经过时,可参考处理方案,不再使用", DeprecationWarning)
        stx, sty, stz = xyzcells
        # jxyz将要不依赖于文件
        global num
        file = open(file=jxyz, mode="r")
        self.un_need_j = []
        topo = topo_tools(topo_file)
        data_jxyz = []
        data_j = set()
        count = 0
        while 1:
            line = file.readline()
            if not line:
                break
            middle = [int(i) for i in line.split()]
            if len(middle) != 4:
                continue
            j, x, y, z = middle
            if not data_j.__contains__(j):
                data_j.add(j)
                data_jxyz.append([j, x, y, z])
        num = len(data_j)
        print("共有", num, "个数据要判断是否是空气")
        # #排序来适应check_point_cache模式
        # def sort_key(element):
        #     return element[3]
        # data_jxyz.sort(key=sort_key,reverse=True)
        start_time = time.time()
        all_time = 1
        res_interval = 5000
        for i in tqdm.tqdm(data_jxyz):
            j, x, y, z = i
            x = stx + (x - 1) * step_x + step_x / 2
            y = sty + (y - 1) * step_y + step_y / 2
            z = stz + (z - 1) * step_z  # + step_z / (2 + 1)  # ?????????????????此处不同写法会导致边缘部分出问题

            # if count == 85:
            #     print(1)
            # print(count)
            #
            res = topo.check_point_cache(x, y, z)

            if res[0]:
                self.un_need_j.append(j)
            count += 1

            # 以下代码为手动实现进度条功能
            """if count % res_interval == 0:
                end_time=time.time()
                all_time+=end_time-start_time
                print("已完成", count,"剩余",num-count,"本轮速度每%d耗时%.3f s"%(res_interval,end_time-start_time),"本轮加速格子数量:",topo.get_cache_count()
                      ,"预计剩余时间:%.3f 分钟"%((num-count)/count*all_time/60),"总耗时:%.3f"%(all_time/60),"分钟")
                start_time = end_time"""
        un_cache = len(data_jxyz) - topo.get_cache_count()
        if un_cache != 0:
            myPrint_Err("共有", un_cache, "次缓存没有被命中,这些将要遍历整个topo文件,请确认topo数据的精度,调整缓存宽松度")
        file.close()
        file = open(un_need_file, mode="w")
        for j in self.un_need_j:
            file.write(str(j))
            file.write('\n')
        file.close()
        print("共有", len(self.un_need_j), "个数据被判断为空气")
        return self.un_need_j

    def get_air_j_from_file(self, unneed_j=r"E:\vscode\Muon_Imaging_Algorithm\dataTools\data\unneed_j"):
        """
        从文件中获取空气格子编号

        :param unneed_j: 存放客气格子编号的文件路径
        :return: 空气格子编号
        :rtype: 列表
        """
        file = open(file=unneed_j, mode="r")
        self.un_need_j = set()
        while 1:
            line = file.readline()
            if not line:
                break
            self.un_need_j.add(int(line))

        return self.un_need_j

    def recover_resfile_by_airj(self, file_res, air_j=None, default_value=0):
        """
        通过空气格子的编号将结果中的空气格子密度值设为默认值。

        :param file_res: 存放密度值结果的文件路径
        :param air_j: 空气格子的编号
        :param default_value: 空气的默认密度值
        :return: 将空气格子密度设为默认值后的密度值结果
        """
        if air_j is None:
            air_j = set(self.un_need_j)
        if type(air_j) is not set:
            air_j = set(air_j)
        old = []
        with open(file=file_res, encoding='utf-8') as file_obj:
            while 1:
                line = file_obj.readline()
                if not line:
                    break
                old.append(line.replace('\n', ''))
            file_obj.close()
        for j in air_j:
            old[j - 1] = default_value
            # try:
            #     old[j-1]=default_value
            # except :
            #     old[j - 1] = default_value#调试用非正常代码
            #     pass
        fo = open(file_res, "w")
        for datum in old:
            fo.write(str(datum))
            fo.write('\n')
        fo.close()
        return old

    def recover_resj_by_airj(self, res, air_j=None, default=-0.1234):
        """
        通过空气格子的编号将结果中的空气格子密度值设为默认值。

        :param res: 计算得到的密度值列表
        :param air_j: 空气格子的编号
        :param default: 空气格子的默认密度值
        :return: 将空气密度设置为默认值后的密度值结果
        """
        if air_j is None:
            air_j = set(self.un_need_j)
        if type(air_j) is not set:
            air_j = set(air_j)
        for i in range(len(res)):
            if i+1 in air_j:
                res[i] = default
        return res


if __name__ == '__main__':
    a = Air_j()
    start_time = time.time()
    mesh_tool = MeshTools(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\real_data\45_58MaMian.msh")
    # a.get_air_j_to_file(mesh_tool,r"E:\vscode\Muon_Imaging_Algorithm\data\Input\real_data\58mamian_survey_821486.topo",r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\unneed_j1"
    #                     ,(1,1,1)
    #                     )
    a.get_air_j_to_file_Jxyzfile_cache(mesh_tool,
                                       r"E:\vscode\Muon_Imaging_Algorithm\data\Input\real_data\58mamian_survey_821486.topo",
                                       r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\unneed_j2",
                                       r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\jxyz"
                                       , (1, 1, 1)
                                       )
    end_time = time.time()
    print("标记空气花费时间为", end_time - start_time, "s")
    # a.get_air_j_from_file()
