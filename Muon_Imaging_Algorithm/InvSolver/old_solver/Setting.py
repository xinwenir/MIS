# author:高金磊
# datetime:2021/10/26 19:22
from datetime import datetime

import scipy.sparse.linalg
# import sympy as sp
# import numpy as np
# from DataTools.G import G_Data
# from DataTools.d_tools import d_tools
# accuracy=100#计算精度
# count=0#计数器
#
# t = sp.symbols('t')
# # M=sp.Matrix([[sp.sin(0.6*t),sp.cos(0.6*t), -sp.sin(0.6*t)],[-sp.cos(0.6*t),sp.sin(0.6*t),sp.cos(0.6*t)]])
# # M=sp.Matrix([[2+sp.sin(t),3+sp.cos(t),4-sp.sin(2*t)+sp.cos(t)]])
# # g=G_Data(m=7600)
# print("读取文件.....%s"%(datetime.now()))
# shape = (91,34,13)
# g=G_Data(n=21228,m=shape[0]*shape[1]*shape[2])
# # M=np.mat(g.get_GArray())
# M=g.get_G_csc_G()
# Mp=np.mat(M.toarray())
# Mp=np.linalg.pinv(Mp)
# Mp=M.T*scipy.sparse.linalg.inv(M*M.T)
# d=d_tools(r"E:\vscode\Muon_Imaging_Algorithm\dataTools\data\d").get_data()
# B=sp.Matrix([[sp.sin(t*2)+sp.cos(t)]])
# B=np.mat(d).T
# L=None
# L=sp.Matrix(B).T
# a:sp.Matrix=sp.Matrix([[0,1,0,0.25,0.25,0.25]]).T
# L=L.col_join(a)
gamma = 10
lmd = 1

desc_size_globle=15