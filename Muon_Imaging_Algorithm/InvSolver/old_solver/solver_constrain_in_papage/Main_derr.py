# author:高金磊
# datetime:2021/11/25 17:11
import gc
import random
from datetime import datetime

import numpy.linalg.linalg
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import lsqr

from DataTools.Air_j import Air_j
from DataTools.Gij_tools import G_Data
from DataTools.d_tools import d_tools
from DataTools.new_res_tools import new_res_tools
from DataTools.res_tools import res_tools

accuracy = 100  # 计算精度
count = 0  # 计数器
log = open(r"E:\vscode\Muon_Imaging_Algorithm\dataTools\data\log.txt", "a")


def make_k(value, vi, vj, m, n):
    value += [-1.0 for i in range(n)]
    vi += [m + i for i in range(n)]
    vj += [i for i in range(n)]
    value += [1.0 for i in range(n)]
    vi += [m + n + i for i in range(n)]
    vj += [i for i in range(n)]
    value += [x0[i + n] for i in range(2 * n)]
    vi += [m + i for i in range(2 * n)]
    vj += [n + i for i in range(2 * n)]
    return csc_matrix((value, (vi, vj)), (m + 2 * n, 3 * n))


def solver_constrain(K, L, x0, shape, max_iter=100, show=False, min_p=0.0, max_p=7.0):
    print("开始计算.....%s" % (datetime.now()))
    log.write("开始计算.....%s\n" % (datetime.now()))
    count = 0

    res = []
    err_middle = 0

    m, n = shape
    while count < max_iter:
        count += 1
        # gc.collect()
        K = K.tocsr()
        for i in range(n):
            if x0[i] < min_p:
                x0[i] = min_p
            elif x0[i] > max_p:
                x0[i] = max_p
        for i in range(n + 1, n + m):
            x0[i] = 0
        res = lsqr(K, L, show=show, x0=x0, iter_lim=500 + count * 2)
        x0 = res[0]
        xs = x0[:n]
        err_count = 0
        for x in xs[:n - m]:
            if x < min_p or x > max_p:
                err_count += 1
        print(err_count, min(xs[:n - m]), max(xs[:n - m]))
        err = K * x0.T - L
        log.write("当前计算的数值不在约束内的有%s个,其中最小值为%s,最大值为%s" % (err_count, min(xs), max(xs)))
        err_p = numpy.linalg.linalg.norm(err[:n - m])
        print(count, res[0][:4], "密度误差", err_p, "d的误差", numpy.linalg.linalg.norm(err[n - m:n]))
        log.write("%s外部迭代次数%s本次迭代%s次,误差%s,部分值:%s\n" % (datetime.now(), count, res[2], res[3], res[0][:3]))
        # print(res[0],res[3])
        K = K.tolil()
        middle = 0
        for re in res[0][n:]:
            K[m + middle, n + middle] = re
            middle += 1
        if abs(err_middle - err_p) < 0.01 and count > max_iter / 10:
            print("迭代的效果已经不明显,将要跳过计算...%s\n" % (datetime.now()))
            log.write("迭代的效果已经不明显,将要跳过计算...%s\n" % (datetime.now()))
            count = max(max_iter - 2, count)
        else:
            # err_count_middle = err_count
            err_middle = err_p
    return res


def format_res(res, unneed_j, air_j, def_value):
    print("计算完毕,正在读写文件...%s" % (datetime.now()))
    log.write("计算完毕,正在读写文件...%s\n" % (datetime.now()))
    fo = open(r"E:\vscode\Muon_Imaging_Algorithm\dataTools\data\res", "w")
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


def filter_res(res):
    ##去掉自由变量和规范数据范围
    print("去掉无效数据.....%s" % (datetime.now()))
    log.write("去掉无效数据.....%s\n" % (datetime.now()))
    count = 0

    for i in range(len(res)):
        if i not in j_effective:
            res[i] = min_p
            count += 1
        # elif res[i] < min_p:
        #     res[i] = min_p
        #     count += 1
        # elif res[i] > max_p:
        #     res[i] = max_p
        #     count += 1
    print("共有%s个无效数据.....%s" % (count, datetime.now()))
    log.write("共有%s个无效数据.....%s" % (count, datetime.now()))
    return res


def simplify_value(vi, vj, unneed_j, value):
    """
    删除vj中不需要的元素，压缩vj
    :param vi:
    :param vj:
    :param unneed_j:
    :return: 旧的j原来的位置
    """
    oldj_newj = {}
    count = 0
    new_value = []
    # 为每一个j分配新的位置
    for i in range(len(vi)):
        if vj[i] not in unneed_j and vj[i] not in oldj_newj.keys():
            oldj_newj[vj[i]] = count
            count += 1
    j_new = []
    i_new = []
    for i in range(len(vi)):
        if vj[i] in oldj_newj.keys():
            j_new.append(oldj_newj[vj[i]])
            i_new.append(vi[i])
            new_value.append(value[i])
    return i_new, j_new, oldj_newj, new_value


def restore_res(res, oldj_newj, n):
    res_new = [0] * n
    for key in oldj_newj.keys():
        res_new[key] = res[oldj_newj[key]]
    return res_new


if __name__ == '__main__':
    min_p = -0.1
    max_p = 4.0
    print("读取文件.....%s" % (datetime.now()))
    log.write("读取文件.....%s\n" % (datetime.now()))
    shape = (280,144,63)

    g = G_Data(m=37328, n=shape[0] * shape[1] * shape[2])
    d = d_tools()
    m = g.shape[0]
    air_j = Air_j()
    unneed_j = air_j.get_air_j_from_file()
    # unneed_j.clear()

    print("数据预处理.....%s" % (datetime.now()))
    log.write("数据预处理.....%s\n" % (datetime.now()))
    value = g.get_GV()
    ij = g.get_ij()

    middle_1 = g.shape  # 删除g来节省内存,这里拷贝g中有用的信息
    del g

    # vi=[0,0,0,1,1]
    # vi = [i - 1 for i in ij[0]]
    vi = ij[0]
    for i in range(len(vi)):
        vi[i] = vi[i] - 1

    # vj=[0,1,2,0,2]
    # vj = [i - 1 for i in ij[1]]
    vj = ij[1]
    for i in range(len(vj)):
        vj[i] = vj[i] - 1
    del ij

    # 方程组精简
    vi, vj, oldj_newj, value = simplify_value(vi, vj, set(unneed_j), value)
    # 添加d的误差到方程组
    n = max(vj) + 1
    value += [1] * m
    vi += [i for i in range(m)]
    vj += [i + n for i in range(m)]
    n = n + m

    j_effective = {i for i in vj}
    d_data = d.get_d()
    # x0 = [random.random() for i in range(n)] + [0] * (2 * n)
    x0 = [random.uniform(min_p, max_p) for i in range(n - m)] + [0.0 for i in range(m)] + [random.random() * 2 - 1 for i
                                                                                           #                                                                                        in range(2 * n)]
                                                                                           # x0 = [random.uniform(min_p, max_p) for i in range(n - m)] + d_data+ [random.random() * 2 - 1 for i
                                                                                           in range(2 * n)]
    # value=[1.0,2.0,1.0,3.0,1.0]
    K: lil_matrix = make_k(value, vi, vj, m, n).tolil()
    # 删除对象节约内存
    del vi, vj, value
    gc.collect()

    lg = [-min_p] * (n - m)
    # lg=[1,1,100]
    lt = [max_p] * (n - m)
    # lt=[-0.1,-0.1,100]
    # L=[1.0,1.0]+lg+lt

    d_err_range = d.get_d_err_range()
    # L =  d_data+ lg+[-(i*-0.2) for i in d_data] + lt+[i * 0.2 for i in d_data]
    # L = [0.0 for i in range(len(d_data))] + lg + [-i for i in d_err_range[0]] + lt + d_err_range[1]
    L = d_data + lg + [-i for i in d_err_range[0]] + lt + d_err_range[1]

    del d_err_range, d
    gc.collect()
    res = solver_constrain(K, L, x0, (m, n), max_iter=100, show=False, min_p=min_p, max_p=max_p)[0]

    # 删除对象节约内存
    del d_data, K, L
    gc.collect()
    res = restore_res(res, oldj_newj, middle_1[1])

    # filter_res(res)#因为已经没有无效数据，此处无需再删除无效数据
    # 重新排列结果,并对结果进行处理(标记空气)
    format_res(res, unneed_j, air_j, -0.1234)

    tool = new_res_tools()
    tool.mode1(min=0.01, ignore_value=[0, -0.1234], max=4.0001)
    tool.mode2(ignore_value=[0, -0.1234])

    print("程序运行完毕...%s" % (datetime.now()))
    log.write("程序运行完毕...%s\n" % (datetime.now()))
    log.close()
