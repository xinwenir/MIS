# author:高金磊
# datetime:2022/3/27 18:40

from copy import copy

import numpy as np

from InvSysTools.MyTools import myPrint
from InvDataTools.res_tools import restore_res
from InvSolver.Jacobi import get_neighbor_node_diff, calculate_dxyz
from InvDataFactory import Setting,DataManage

count = 0
min_o3 = 0
min_o1 = 0
o1_last = 0
o3_last = 0
min_res = 0


def reset():
    global count, min_o3, min_o1, o1_last, o3_last, min_res
    count = 0
    min_o3 = 0
    min_o1 = 0
    o1_last = 0
    o3_last = 0
    min_res = 0


def constraint(x, bounds):
    count = 0
    for i in range(len(x)):
        if x[i] < bounds[i][0] or x[i] > bounds[i][1]:
            count += 1
    return count * 100


def get_refs(x, m0, newj_oldj):
    dxyz = calculate_dxyz(newj_oldj, newj=-1)
    refs = x - np.array([m0]).T
    res = 0
    for i in range(len(refs)):
        res += refs[i] * refs[i] * dxyz[i]
    return res

x_good=None
def constr_f(x, A, b, shape, d_err, newj_oldj, oldj_newj, m0, bounds, stp, args):
    """
    目标函数
    :param x:
    :param A: Ax=b的A
    :param b: Ax=b的b
    :param col: 列总数
    :param a_col: 除去添加的d误差列以外的列总数
    :param d: 实际测得的d
    :param x_sum: 二维数组，内部的每个一维数组的长度为x轴切得的格子数
    :param y_sum: 二维数组，内部的每个一维数组的长度为y轴切得的格子数
    :param z_sum: 二维数组，内部的每个一维数组的长度为z轴切得的格子数
    :param m: 计算得到的密度
    :param m0: 参考的密度
    :return: x的一维数组
    """
    global x_good
    setting = Setting.Setting.get_instance()
    dataManager=DataManage.DataManager.get_instance()
    if stp.is_STP():
        setting.get_loger().write("捕获指令STP,正在尝试停止求解器", printer=myPrint.myPrint_Wran)
        if x_good is None:
            x_good=x
        raise stp.stp_Exception((x_good,))
    col = len(x)
    # a_col = col - len(b)#考虑b的误差
    a_col = col
    alpha_s, alpha_x, alpha_y, alpha_z = args

    # from InvDataTools.middle.refs_bound_simple_build import Cov_array
    # print(Cov_array(x))

    global count, min_o3, min_o1, min_res, o1_last, o3_last
    x = np.array([x]).T
    b = np.array([b]).T
    m0T = np.array([m0]).T
    pre_d_err = (A * x - b)
    for i in range(len(d_err)):
        #     pre_d[i][0] *= abs(1-d_err[i]/b[i])
        pre_d_err[i][0] /= d_err[i]

    o1 = pow(np.linalg.norm(pre_d_err), 2)
    o2 = pow(np.linalg.norm(x[a_col:]), 2)

    count += 1
    o3s = get_refs(x, m0, newj_oldj)
    # o3s = sum(pow((x[:a_col] - m0T), 2)) * 0.25 * 0.25 * 0.2
    # oo3 = get_xyzp_m0(x_sum, y_sum, z_sum,m0,shape,newj_oldj)?????????
    # x_sum, y_sum, z_sum = getx_y_z_continuous(shape, x[:a_col]
    #                                           - [i for i in m0]
    #                                           , newj_oldj)
    # oo3 = get_xyzp(x_sum, y_sum, z_sum, shape)

    # o3x = alpha_x * oo3[0]
    # o3y = alpha_y * oo3[1]
    # o3z = alpha_z * oo3[2]
    # o3 = o3s + o3x+o3y+o3z
    # oo3 = o3x + o3y + o3z

    ord=2
    if ord==2:
        oo3 = sum(get_smooth(x, m0, oldj_newj, newj_oldj, shape, p=2,coefficients=(alpha_x, alpha_y, alpha_z)))
        # oo3 = pow(oo3, 0.5)
    elif ord==1:
        oo3 = np.linalg.norm(get_smooth(x, m0, oldj_newj, newj_oldj, shape,p=1, coefficients=(alpha_x, alpha_y, alpha_z)),1)
    else:
        smmoth = get_smooth(x, m0, oldj_newj, newj_oldj, shape,p=2,coefficients=(alpha_x, alpha_y, alpha_z))
        oo3 = pow(sum(smmoth), ord / 2)
        # oo3 = pow(sum(smmoth), ord / 2)

    o3 = alpha_s * o3s + oo3
    o3_last = o3
    o1_last = o1
    res = o1 + dataManager.choose_beta().get_beta()[0] * o3  # +constraint(x,bounds)
    info = "%d目标函数误差:%e方程组误差:%e偏离ref:%e平滑度误差: %ed的误差:%e" % (count, res, o1, o3s, oo3, o2)
    setting.get_loger().write(info, printer=myPrint.myPrint_Success)
    dataManager.collect_misfit_refs_smooth(res[0], o1, o3s, oo3)
    # o3=pow(beta/len(b)*o3,0.5)
    # return [i[0]+(o3)[0] for i in pre_d]
    # return [i[0] for i in pre_d]
    if count == 1:  # 避免取到初始值
        pass
        # res=10e10
        min_o1 = res
        min_o3 = res
        min_res = res
    else:
        if res < min_res:
            min_res = res
            setting.get_loger().write("当前的最优解", printer=myPrint.myPrint_Hint)
            x_good=x
        if o1 < min_o1:
            min_o1 = o1
            setting.get_loger().write("当前方程组误差最小", printer=myPrint.myPrint_Hint)
        if o3 < min_o3:
            min_o3 = o3
            setting.get_loger().write("当前平滑度最小", printer=myPrint.myPrint_Hint)

    # jac_refs=get_jac_refs(x,m0)
    # equations_jac=get_jac_equations(x,pre_d,A)
    # jac+=[-i[0]/100 for i in pre_d]
    # CG的jac要求是(一维列表)
    # return res,(20*equations_jac+0*jac_refs).T[0]
    # 标准
    # return res, (200 * equations_jac + 1 * jac_refs)
    # 最小二乘
    if count > setting.get_max_iter():
        stp.commands(1)
    return res


def getxyz_from_shape(shape, j):
    """
    代替jxyz文件的一种补充方案，之后可能代替jxyz文件的读写
    起始坐标是（1，1，1）与jxyz保持一致
    :param shape:
    :param j: range(1,shape[0]*shape[1]*shape[2])
    :return:
    """
    j -= 1
    y = int(j / (shape[0] * shape[2]))
    j = j % (shape[0] * shape[2])
    x = int(j / shape[2])
    z = j % shape[2]
    return [x + 1, y + 1, shape[2] - z]


def getj_from_xyz(shape, xyz):
    """
    getj_from_xyz
    :param shape:
    :param xyz:
    :return: j range(1,shape[0]*shape[1]*shape[2])
    """
    x, y, z = xyz
    return shape[0] * shape[2] * (y - 1) + shape[2] * (x - 1) + (shape[2] - z + 1)


def get_smooth(x0, m0, oldj_newj, newj_oldj, shape, p=2.0,coefficients=None):
    """
    为了保持一致,此处的获取邻居结点的j复用雅可比矩阵中的方法
    :param x0:
    :param m0:
    :param oldj_newj:
    :param newj_oldj:
    :param shape:
    :param p:
    :param coefficients: (ax, ay, az)
    :return:
    """

    mesh = DataManage.DataManager.get_instance().mesh
    if coefficients is None:
        ax, ay, az = 1, 1, 1
    else:
        ax, ay, az = coefficients
    x0 = x0 - np.array([m0]).T
    x0 = x0.T.tolist()[0]
    xx = copy(x0)  # debug用
    x0 = restore_res(res=x0, oldj_newj=oldj_newj, shape=shape)
    jac = [0] * len(m0)
    cell_xs = mesh.get_xs()
    cell_ys = mesh.get_ys()
    cell_zs = mesh.get_zs()
    for i in range(len(jac)):
        j = newj_oldj[i] + 1
        x, y, z = getxyz_from_shape(shape, j)
        j -= 1
        assert xx[i] == x0[j], "转换关系不正确"
        x0.append(0)
        default_j = (j, -1)[1]
        dx = cell_xs[x - 1]
        dy = cell_ys[y - 1]
        dz = cell_zs[z - 1]
        if x - 1 <= 0:
            partial_derivatives = cell_xs[x - 1]
        else:
            partial_derivatives = (cell_xs[x - 1] + cell_xs[x - 2]) / 2
        jac[i] += pow(
            (x0[j] - x0[
                get_neighbor_node_diff(shape, default_j, oldj_newj.keys(), [x - 1, y, z])]) / partial_derivatives,
            p) * ax * dx * dy * dz

        if y - 1 <= 0:
            partial_derivatives = cell_ys[y - 1]
        else:
            partial_derivatives = (cell_ys[y - 1] + cell_ys[y - 2]) / 2
        jac[i] += pow(
            (x0[j] - x0[
                get_neighbor_node_diff(shape, default_j, oldj_newj.keys(), [x, y - 1, z])]) / partial_derivatives,
            p) * ay * dx * dy * dz

        if z - 1 <= 0:
            partial_derivatives = cell_zs[z - 1]
        else:
            partial_derivatives = (cell_zs[z - 1] + cell_zs[z - 2]) / 2
        jac[i] += pow((x0[j] - x0[
            get_neighbor_node_diff(shape, default_j, oldj_newj.keys(), [x, y, z - 1])]) / partial_derivatives,
                      p) * az * dx * dy * dz
    return np.array(jac)


"""
以下代码已经被计算公示替代
"""
j_xyz_cache = {}


def getx_y_z_continuous(shape: list, old_ps: list, newj_oldj: dict):
    """
    通过old_ps，获取x，y，z方向的每个格子的密度值（包含空气）
    :param shape: x,y,z方向的格子数目
    :param old_ps: 不完全的密度值
    :param newj_oldj: 当前密度下标对应的原始密度链表的下标，字典
    :return: 每个都是二维链表，[i]返回yi，zi相同，x不同的密度构成的链表
    """
    global j_xyz_cache
    # new_ps = restore_res(old_ps, newj_oldj, shape)#恢复res
    x_max, y_max, z_max = shape
    x_continuous = np.zeros((y_max, z_max, x_max))
    y_continuous = np.zeros((x_max, z_max, y_max))
    z_continuous = np.zeros((x_max, y_max, z_max))
    # data=numpy.zeros(shape)
    for i in range(len(old_ps)):
        p = old_ps[i]
        if p == 0:
            continue

        j = newj_oldj[i]

        if j in j_xyz_cache.keys():
            x, y, z = j_xyz_cache[j]
        else:
            x, y, z = getxyz_from_shape(shape, j)
            j_xyz_cache[j] = [x, y, z]
        # data[x-1][y-1][z-1]=p
        # if x==1:
        #     print(p)
        x_continuous[y - 1][z - 2][x - 1] = p
        y_continuous[x - 1][z - 2][y - 1] = p
        z_continuous[x - 1][y - 1][z - 2] = p
    middle = []
    for xs in x_continuous:
        for x in xs:
            middle.append(x)
    x_continuous = middle
    middle = []
    for ys in y_continuous:
        for y in ys:
            middle.append(y)
    y_continuous = middle
    middle = []
    for zs in z_continuous:
        for z in zs:
            middle.append(z)
    z_continuous = middle
    return np.array(x_continuous), np.array(y_continuous), np.array(z_continuous)


x_p, y_p, z_p = None, None, None


def get_xyzp_m0(x_sum, y_sum, z_sum, m0, shape, newj_oldj):
    global x_p, y_p, z_p
    if x_p is None:
        x_p, y_p, z_p = getx_y_z_continuous(old_ps=m0, shape=shape, newj_oldj=newj_oldj)
    sum_xp = sum_yp = sum_zp = 0
    for i in range(len(x_p)):
        for j in range(len(x_p[i])):
            sum_xp += pow((x_p[i][j] - x_sum[i][j]) / 0.25, 2)
    for i in range(len(y_p)):
        for j in range(len(y_p[i])):
            sum_yp += pow((y_p[i][j] - y_sum[i][j]) / 0.25, 2)
    for i in range(len(z_p)):
        for j in range(len(z_p[i])):
            sum_zp += pow((z_p[i][j] - z_sum[i][j]) / 0.2, 2)

    return sum_xp, sum_yp, sum_zp

    pass


def get_interval(path_meshf):
    with open(path_meshf) as f:
        f.readline()
        f.readline()
        x = np.array(f.readline().strip().split(), dtype=np.float64)
        y = np.array(f.readline().strip().split(), dtype=np.float64)
        z = np.array(f.readline().strip().split(), dtype=np.float64)
        x_interval = x[1:] + x[:-1]
        y_interval = y[1:] + y[:-1]
        z_interval = z[1:] + z[:-1]
        return x_interval, y_interval, z_interval


def get_xyzp(x_sum: np.array, y_sum, z_sum, shape):
    sum_xp = sum_yp = sum_zp = 0
    # res = get_interval("D:\\Projects\\fortran\\TestData\\Data_2_23\\45_58MaMian.msh")
    res = [np.array([0.25] * (shape[0] - 1)), np.array([0.25] * (shape[1] - 1)), np.array([0.2] * (shape[2] - 1))]
    dv = 0.25 * 0.25 * 0.2
    for x in x_sum:
        # xp = x[1:] - x[:-1]
        # yp = y[1:] - y[:-1]
        # zp = z[1:] - z[:-1]
        xp = []
        flag = True
        for i in range(1, len(x)):
            if x[i] == 0:
                xp.append(0)
                flag = True
                continue
            else:
                if flag:
                    xp.append(x[i] - x[i - 1])
                else:
                    xp.append(0)
                    flag = True
        # sum_xp += sum(pow(np.array(xp), 2) / res[0])
        sum_xp += sum(pow(np.array(xp) / 0.25, 2)) * dv

    for y in y_sum:
        yp = []
        flag = True
        for i in range(1, len(y)):
            if y[i] == 0:
                flag = True
                yp.append(0)
                continue
            else:
                if flag:
                    yp.append(y[i] - y[i - 1])
                else:
                    yp.append(0)
                    flag = True
        # sum_yp += sum(pow(np.array(yp), 2) / res[1])
        sum_yp += sum(pow(np.array(yp) / 0.25, 2)) * dv

    for z in z_sum:
        zp = []
        flag = True
        for i in range(1, len(z)):
            if z[i] == 0:
                flag = True
                zp.append(0)
                continue
            else:
                if flag:
                    zp.append(z[i] - z[i - 1])
                else:
                    zp.append(0)
                    flag = True
        # sum_zp += sum(pow(np.array(zp), 2) / res[2])
        sum_zp += sum(pow(np.array(zp) / 0.2, 2)) * dv
    return sum_xp, sum_yp, sum_zp


if __name__ == '__main__':
    shape = (280, 144, 63)
    xyz = getxyz_from_shape(shape, 1)
    print(xyz)
    print(getj_from_xyz(shape, xyz))
    print(getj_from_xyz(shape, xyz))
    # x = np.array([[random.random() for _ in range(131)]])
    # y = np.array([[random.random() for _ in range(92)]])
    # z = np.array([[random.random() for _ in range(55)]])
    # bounds=[[1,1]]
    # differential_evolution(constr_f,bounds,args=())

# row = [0, 1]
# col = [0, 2]
# data = [1, 2]
# a = sparse.csc_array((data, (row, col)), shape=(2, 3))
# a = [[1, 1, 1], [2, 2, 2]]
# b = [[1], [2]]
# bounds = [[-0.1, 0.1],[-0.1, 0.1],[-0.1, 0.1]]
# result = differential_evolution(constr_f, bounds,disp=True)
# print(result.x, result.fun)
