# author:高金磊
# datetime:2021/11/17 10:45
import sys

import numpy
from Cython import profile

from InvDataTools import d_tools
import numpy as np
from scipy import linalg
from scipy.sparse import csc_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
import scipy as sc
from scipy import integrate


class G_Data():
    """
    读取并处理“G”文件和“Gij"文件

    """
    def __init__(self, fileG, file_ij):
        self._GV = []
        # self._ij=[]#此方法内存过高
        self._ij = [[], []]

        # self._G=[[0 for i in range(m)] for i in range(n)]
        with open(file=fileG, encoding='utf-8') as file_obj:
            while 1:
                line = file_obj.readline()
                if not line:
                    break
                if line.replace(' ', '') != '':
                    self._GV.append(float(line))
            file_obj.close()
        with open(file_ij, encoding="utf-8") as file_obj:
            file_obj.readline()  # 需要去掉第一行
            while 1:
                line = file_obj.readline()
                if not line:
                    break
                middle = line.split()
                self._ij[0].append(int(middle[0]))
                self._ij[1].append(int(middle[1]))
            file_obj.close()

    def get_ij(self):
        """
        获取射线及对应的格子编号

        :return: 射线及对应的格子编号
        """
        return self._ij

    def get_GV(self):
        """
        获取射线穿过某个格子的路径长度

        :return: 射线穿过某个格子的路径长度
        """
        return self._GV

    def get_G_csc_G(self, shape):
        """
        此方法禁止使用,灵活性太低不利于内存的优化

        生成G(Gρ=d)

        :param shape:
        :return: 返回使用稀疏矩阵表示方法的矩阵
        """
        import warnings
        warnings.warn("不推荐直接使用此方法", DeprecationWarning)
        return csr_matrix((self._GV, ([i - 1 for i in self._ij[0]], [i - 1 for i in self._ij[1]])), shape, dtype=float)

    def get_GArray(self, shape):
        """
        此方法禁止使用,反演问题中方程组维度太大

        获取G(Gρ=d)的numpy.ndarray形式

        :param shape: 模型的x、y、z方向的格子数目
        :return: 矩阵G的numpy.ndarray形式
        """
        import warnings
        warnings.warn("禁止使用此方法", DeprecationWarning)

        return self.get_G_csc_G(shape).toarray()
        # for i in range(len(self._GV)):
        # print(self._ij[i])
        # self._G[self._ij[i][0]-1][self._ij[i][1]-1]=self._GV[i]


def get_g(Gij, G):
    """
    根据ijg文件得到每个格子的真实长度，并存放在文件中

    :param Gij: ijg文件路径
    :param G: 存放g的文件路径
    """
    ijp = open(Gij, 'r')
    ijps = ijp.readlines()
    res = open(G, 'w')
    for p in ijps[1:]:
        res.write(p.strip().split()[2])
        res.write('\n')
    res.flush()
    res.close()
    ijp.close()
# if __name__ == '__main__':
# g=G_Data(21228,40222)
# print(g)
# d=d_tools.d_tools()
# gM=np.matrix(g._G)
# res=gM.T*(np.linalg.inv(gM*gM.T))
# dM=np.matrix(d.get_data())
# gamma = 100
# lmd = 10
# x0 = [0]*g.shape[1]
# ts = np.linspace(0, 10, 1000)
# sp_A = csc_matrix((g._GV, ([i-1 for i in g._ij[0]], [i-1 for i in g._ij[1]])), shape=g.shape, dtype=float)

# x=spsolve(sp_A.T*sp_A,sp_A.T*d.get_data())
# data = sc.integrate.odeint(OdeRightFunc, x0, ts, (gamma, lmd))
#
# data=np.linalg.pinv(sp_A)*d
# print(x)
