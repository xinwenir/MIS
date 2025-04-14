# # author:高金磊
# # datetime:2021/11/14 9:43
# from random import random
#
# from solver.Constrains.OdeRight import OdeRightFunc
# from solver.equation.A import *
# from solver.equation.L import *
# from solver.equation.K import *
# from solver.equation.V import *
# from solver.AFM import *
# import matplotlib.pyplot as plot
# import scipy as sc
# from scipy import integrate
# class Solver():
#     """
#     不允许多线程
#     """
#     def __init__(self,M,V,const):
#         self.M=M
#         self.V=V
#         self.const=const
#         self.AFM1=Bound
#         self.AFM2=linear
#         global m,n
#         m,n=M.shape()
#     def set_AFM1(self,AFM1):
#         self.AFM1=AFM1
#         return self
#     def set_AFM2(self,AFM2):
#         self.AFM2=AFM2
#         return self
#     def solver(self,gamma = 10,lmd = 1):
#         # x0 = [0.1 for i in range(3 * n)] + [0] * (m + 2 * n)
#         x0=[random() for i in range(5*n+m)]
#         # x0=[0]*(5*n+m)
#         ts = np.linspace(0, 10, 100)
#         data = sc.integrate.odeint(OdeRightFunc, x0, ts, (gamma, lmd,self.AFM1,self.AFM2))
#         xs = []
#         norms = []
#
#         for i in range(len(data)):
#             xs.append(data[i][:3 * n])
#             ct = data[i][n:3 * n]
#             norms.append(sp.Matrix(At_value(ts[i], ct) * sp.Matrix(data[i][:3 * n]) - lt_value(ts[i])).norm())
#         self.ts=ts
#         self.norms=norms
#         self.data=data
#         self.xs=xs
#
#         for i in range(len(xs)):
#             # print(ts[i],*xs[i][:3],norms[i])
#             print(ts[i], *xs[i][:3])
#             # print(xs[i])
#         pass
#     def show_norms(self):
#         ts=self.ts
#         norms=self.norms
#         xs=self.xs
#         plot.plot(ts, norms)
#         plot.show()
#         plot.figure(1)
#         plot.subplot(221)
#         plot.axhline(-a[0])
#         plot.axhline(a[3])
#         plot.plot(ts, [xs[i][0] for i in range(len(xs))])
#         plot.subplot(222)
#         plot.axhline(-a[1])
#         plot.axhline(a[4])
#         plot.plot(ts, [xs[i][1] for i in range(len(xs))])
#         plot.subplot(223)
#         plot.axhline(-a[2])
#         plot.axhline(a[5])
#         plot.plot(ts, [xs[i][2] for i in range(len(xs))])
#         plot.subplot(224)
#         plot.plot(ts, norms)
#         plot.show()
