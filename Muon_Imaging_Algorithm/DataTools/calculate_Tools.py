# author:高金磊
# datetime:2022/4/27 11:42

from DataTools.Jxyz_Tools import getxyz_from_shape
from DataFactory.DataManage import DataManager
dxyz = None


def calculate_dxyz(newj_oldj: dict, newj=-1):
    """
    计算压缩后的第newj个格子的体积，当newj为-1时，返回压缩后的全部格子的体积

    :param newj_oldj: 格子的新编号与旧编号的字典
    :param newj: 要计算的格子的编号
    :return: 第newj个格子的体积(newj!=-1)或所有格子的体积(newj==-1)
    """
    dataManager=DataManager.get_instance()
    global dxyz
    if dxyz is None:
        dxyz = [0] * len(newj_oldj.keys())
        mesh = dataManager.mesh
        xcells = mesh.get_xs()
        ycells = mesh.get_ys()
        zcells = mesh.get_zs()
        for j in newj_oldj.keys():
            oldj = newj_oldj[j]
            x, y, z = getxyz_from_shape(mesh.get_shape(), oldj)
            dxyz[j] = xcells[x - 1] * ycells[y - 1] * zcells[z - 1]
    if newj == -1:
        # return copy(dxyz)
        return dxyz
    else:
        return dxyz[newj]
