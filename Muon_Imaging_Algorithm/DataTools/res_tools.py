# author:高金磊
# datetime:2021/11/22 20:00
from datetime import datetime

import numpy
import tqdm

from DataTools.Air_j import Air_j
from SysTools.tools import Desc_Text_alignment_tools


class res_tools():
    """
    对结果进行处理

    """
    def __init__(self):
        pass

    # def Conversion_1(self, file_res="M:\pycharm\Inversion\DataTools\data\res",file_out="M:\pycharm\Inversion\DataTools\data\new"):
    #     old = []
    #     with open(file=file_res, encoding='utf-8') as file_obj:
    #         while 1:
    #             line=file_obj.readline()
    #             if not line:
    #                 break
    #             old.append(line.replace('\n',''))
    #     x,y,z=48,48,45
    #     res=[]
    #     for i in range(x):
    #         i+=1
    #         middle1=[]
    #         for j in range(z):
    #             j+=1
    #             middle2=[]
    #             for k in range(y):
    #                 k+=1
    #                 middle2.append(old[i*j*k-1])
    #             middle1.append(middle2)
    #         res.append(middle1)
    #     fo = open("M:\pycharm\Inversion\DataTools\data\middle", "w")
    #     for re in res:
    #         for r in re:
    #             for data in r:
    #
    #                 fo.write(str(data))
    #                 fo.write(' ')
    #         fo.write('\n')
    #     for i in range(z):
    #         i+=1
    #         for j in range(x):
    #             j+=1
    #             for k in range(y):
    #                 k+=1
    #                 old[i*j*k-1]=res[k-1][i-1][j-1]
    #     new_data=old
    #     fo = open("../../../DataTools/data/new", "w")
    #     for datum in new_data:
    #         fo.write(datum)
    #         fo.write('\n')
    def Conversion_2(self, shape=(14, 21, 7), file_res=r"M:\pycharm\Inversion\DataTools\data\res",
                     file_xyz=r"M:\pycharm\Inversion\data\Temp\jxyz",
                     file_out=r"M:\pycharm\Inversion\DataTools\data\new"):
        """
        在以后再支持文件和数据混合使用

        jxyz需要被替代

        对格子密度的计算值进行重新排序(可以看作是坐标系的转换)

        :param shape: 模型在x、y、z方向上的格子数目
        :param file_res: 密度计算值的文件路径(读取)
        :param file_xyz: 格子的编号以及对应的离散坐标的文件路径(读取)
        :param file_out: 存放结果的文件路径(读取)
        """
        old = []
        process = tqdm.tqdm(total=4, desc=Desc_Text_alignment_tools("正在对结果转换"))
        with open(file=file_res, encoding='utf-8') as file_obj:
            while 1:
                line = file_obj.readline()
                if not line:
                    break
                old.append(line.replace('\n', ''))
            file_obj.close()
        jxyzp = {}
        process.update(1)
        with open(file=file_xyz, encoding='utf-8') as file_obj:
            while 1:
                line = file_obj.readline()
                if not line:
                    break
                middle = line.split()
                key = int(middle[0])
                value = [int(i) for i in middle[1:]]
                value.append(float(old[key - 1]))
                jxyzp[key] = value
            file_obj.close()
        data = numpy.zeros(shape)
        process.update(1)
        x, y, z = data.shape
        for xyzp in jxyzp.values():
            xi = xyzp[0] - 1
            yi = xyzp[1] - 1
            zi = xyzp[2] - 1
            if xi < 0 or yi < 0 or zi < 0 or xi >= x or yi >= y or zi >= z:
                continue
            data[xi][yi][zi] = xyzp[-1]
        process.update(1)
        fo = open(file_out, "w")
        for k in range(y):
            for j in range(x):
                for i in range(z - 1, -1, -1):
                    fo.write(str(data[j][k][i]))
                    fo.write('\n')

        fo.close()
        process.update(1)
        # data.sort(key=lambda x:(x.sort(),x[2],x[0],x[1]))
        # fo = open("new", "w")
        # for datum in data:
        #     fo.write(str(datum[3]))
        #     fo.write('\n')


def restore_res(res, oldj_newj,shape):
    """
    恢复结果(所有格子的密度值，oldj_newj.keys()不包含的格子密度都为0)

    :param res: 格子密度的计算值
    :param oldj_newj: 字典，键为压缩前的格子编号，值为压缩后的格子编号
    :param shape: 模型在x、y、z方向上的格子数目
    :return: 恢复后的密度值
    """

    res_new = [0] * (shape[0] * shape[1] * shape[2])
    for key in oldj_newj.keys():
        try:
            res_new[key] = res[oldj_newj[key]]
        except Exception as e:
            print(key, "发生越界")
            print(e)
    return res_new


def format_res(res, unneed_j, air_j, def_value, log):
    """
    重组结果(将空气标记为定值)

    :param res: 格子密度的计算值
    :param unneed_j: 空气格子编号
    :param air_j: 空气类(Air_j)实例
    :param def_value: 空气格子的默认密度值
    :param log: 日志
    """
    import warnings
    warnings.warn("因为res有序此处被简化,仅供参考禁止使用", DeprecationWarning)

    print("计算完毕,正在读写文件...%s" % (datetime.now()))
    log.write("计算完毕,正在读写文件...%s\n" % (datetime.now()))
    fo = open(r"M:\pycharm\Inversion\DataTools\data\res", "w")
    for re in res:
        fo.write(str(re))
        fo.write('\n')
    fo.close()
    print("正在重组结果...%s" % (datetime.now()))
    log.write("正在重组结果.....%s\n" % (datetime.now()))
    try:
        rt = res_tools()
        # 为了方便观察,将空气标记为定值
        air_j.remark_point(default_value=def_value, js=unneed_j)
        rt.Conversion_2(shape=shape)
    except Exception as e:
        print(e)
        print("结果重组失败,请检查原因并尝试使用res_tools工具重组结果")
        log.write("结果重组失败,请检查原因并尝试使用res_tools工具重组结果%s\n" % (datetime.now()))


if __name__ == '__main__':
    # obj=res_tools()
    # obj.Conversion_2(file_res=r"M:\pycharm\Inversion\DataTools\data\res",shape=(182, 69, 50))
    rt = res_tools()
    # 为了方便观察,将空气标记为定值
    air_j = Air_j()
    def_value = -0.12345
    unneed_j = air_j.get_air_j_from_file()
    shape = (46, 17, 7)
    air_j.remark_point(default_value=def_value, js=unneed_j)
    rt.Conversion_2(shape=shape)
