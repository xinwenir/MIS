# author:高金磊
# datetime:2022/4/27 8:45
from copy import copy

import numpy as np

from InvSysTools.MyTools import myPrint
from InvDataTools.Jxyz_Tools import getxyz_from_shape, getj_from_xyz
from InvDataTools.res_tools import restore_res
from InvDataFactory import Setting,DataManage
from InvDataTools.calculate_Tools import calculate_dxyz


def get_jac_smooth(x0, m0, oldj_newj, newj_oldj, shape, coefficients=None):
    """

    :param x0:
    :param m0:
    :param oldj_newj:
    :param newj_oldj:
    :param shape:
    :param coefficients: (ax, ay, az)
    :return:
    """
    dataManager=DataManage.DataManager.get_instance()
    mesh = dataManager.mesh
    if coefficients is None:
        ax, ay, az = 1, 1, 1
    else:
        ax, ay, az = coefficients
    x0 = x0 - np.array([m0]).T
    x0 = x0.T.tolist()[0]
    xx = copy(x0)  # debug用
    x0 = restore_res(res=x0, oldj_newj=oldj_newj, shape=shape)
    jac = np.zeros((len(m0)),dtype=np.float64)
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
        ord = 2

        if x - 1 <= 0:
            partial_derivatives = cell_xs[x - 1]
        else:
            partial_derivatives = (cell_xs[x - 1] + cell_xs[x - 2]) / 2
        jac[i] += pow(x0[j] - x0[get_neighbor_node_diff(shape, default_j, oldj_newj.keys(), [x - 1, y, z])],ord-1) * ax * (
                ord / partial_derivatives) * dy * dz
        if x >= shape[0]:
            partial_derivatives = cell_xs[x - 1]
        else:
            partial_derivatives = (cell_xs[x - 1] + cell_xs[x]) / 2
        jac[i] += pow(x0[j] - x0[get_neighbor_node_diff(shape, default_j, oldj_newj.keys(), [x + 1, y, z])],ord-1) * ax * (
                ord / partial_derivatives) * dy * dz
        if y - 1 <= 0:
            partial_derivatives = cell_ys[y - 1]
        else:
            partial_derivatives = (cell_ys[y - 1] + cell_ys[y - 2]) / 2
        jac[i] += pow(x0[j] - x0[get_neighbor_node_diff(shape, default_j, oldj_newj.keys(), [x, y - 1, z])],ord-1) \
                  * ay * (ord / partial_derivatives) * dx * dz
        if y >= shape[1]:
            partial_derivatives = cell_ys[y - 1]
        else:
            partial_derivatives = (cell_ys[y - 1] + cell_ys[y]) / 2
        jac[i] += pow(x0[j] - x0[get_neighbor_node_diff(shape, default_j, oldj_newj.keys(), [x, y + 1, z])],ord-1) * ay * (
                ord / partial_derivatives) * dx * dz
        if z - 1 <= 0:
            partial_derivatives = cell_zs[z - 1]
        else:
            partial_derivatives = (cell_zs[z - 1] + cell_zs[z - 2]) / 2
        jac[i] += pow(x0[j] - x0[get_neighbor_node_diff(shape, default_j, oldj_newj.keys(), [x, y, z - 1])],ord-1) * az * (
                ord / partial_derivatives) * dx * dy
        if z >= shape[2]:
            partial_derivatives = cell_zs[z - 1]
        else:
            partial_derivatives = (cell_zs[z - 1] + cell_zs[z]) / 2
        jac[i] += pow(x0[j] - x0[get_neighbor_node_diff(shape, default_j, oldj_newj.keys(), [x, y, z + 1])],ord-1) * az * (
                ord / partial_derivatives) * dx * dy
    return np.array(jac)


def get_neighbor_node_diff(shape, default_j, oldjs, xyz):
    """
    可以改成一次获取一个格子周围6个格子的值,进行缓存能加快速度
    不过因为场景很多,会导致方法并不通用
    :param shape:
    :param default_j:
    :param oldjs:
    :param xyz:
    :return:
    """
    x, y, z = xyz
    if x < 1 or y < 1 or z < 1 or x > shape[0] or y > shape[1] or z > shape[2]:
        return default_j
    res = getj_from_xyz(shape, xyz) - 1
    if res not in oldjs:
        return default_j
    return res


def get_jac_equations(x, pred_obsd, A):
    jac = A.T * pred_obsd

    return jac


def correct_jac(x, jac, bounds):
    """
    纠正jac,类似罚函数
    :param x:
    :param jac:
    :param bounds:
    :return:
    """
    limt_max = max(jac)
    limt_min = min(jac)
    for i in range(len(jac)):
        if x[i] < bounds[i][0] - 0.1 and jac[i] > 0:
            jac[i] = limt_min*2
        if x[i] > bounds[i][1] + 0.1 and jac[i] < 0:
            jac[i] = limt_max*2
    return jac


def get_jac_refs(x, m0, newj_oldj):
    dxyz = calculate_dxyz(newj_oldj, newj=-1)
    jac_refs = x - np.array([m0]).T
    for i in range(len(jac_refs)):
        jac_refs[i] *= dxyz[i]
    return 2 * jac_refs


def get_jac_equation_for_least_squares(x, A, b, shape, d_err, newj_oldj, oldj_newj, m0: np.array, bounds, stp, args):
    setting = Setting.Setting.get_instance()
    dataManager=DataManage.DataManager.get_instance()
    if stp.is_STP():
        setting.get_loger().write("捕获指令STP,正在尝试停止求解器", printer=myPrint.myPrint_Wran)
        raise stp.stp_Exception((x,))
    pred_obsd = (A * x - b)
    for i in range(len(d_err)):
        #     pre_d[i] *= abs(1-d_err[i]/b[i])
        pred_obsd[i] /= d_err[i]  # 文献的处理方法
    jac_equations = get_jac_equations(x, pred_obsd, A)  # 等价于python的包: scipy.optimize._lsq.common.compute_grad(A,pre_d)
    x = np.array([x]).T
    jac_refs = get_jac_refs(x, m0, newj_oldj)
    jac_smooth = get_jac_smooth(x, m0, oldj_newj, newj_oldj, shape, args[-3:])
    jac = jac_equations + (args[0] * jac_refs.T[0] + 1 * jac_smooth.T) * 1 * dataManager.choose_beta().get_beta()[0]

    # # 一种情况,边界已经到了,但是还要求他越过边界,今后考虑---额.....原来这就是罚函数的做法
    # if setting.get_solve_method() == "CG":
    #     jac = correct_jac(x, jac, bounds)

    return jac
class Jac_line_search(object):
    def __init__(self,jac_function,args):
        self.jac_function=jac_function
        self.args=args
        pass
    def Jac_line_search(self,x,*others):
        if len(others) == 0:
            return self.jac_function(x,*self.args)
        else:
            return self.jac_function(x,*others)