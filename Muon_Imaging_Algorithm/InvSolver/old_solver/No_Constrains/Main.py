# author:高金磊
# datetime:2021/10/26 20:55

import random

import matplotlib.pyplot as plot
import scipy as sc
from scipy import integrate
# from solver.equation.MA import *

# from solver.equation.VB import *
from solver.old_solver.Setting import lmd,gamma,sp,np,M
from solver.old_solver.No_Constrains.OdeRight import OdeRightFunc
t = sp.symbols('t')
# x0=[random.random() for i in range(M.cols)]+[0]*M.rows
x0=[random.random() for i in range(M.shape[1])]+[0]*M.shape[0]
if __name__ == '__main__':
    ts = np.linspace(0, 10, 1000)
    data = sc.integrate.odeint(OdeRightFunc, x0, ts,(gamma,lmd))
    xs=[]
    norms=[]
    for i in range(len(data)):
        xs.append(data[i][0:M.cols])
        # norms.append(sp.Matrix(Mt_value(ts[i])*sp.Matrix(xs[i])-Bt_value(ts[i])).norm())
    plot.plot(ts,norms)
    plot.show()
    plot.figure(1)
    plot.subplot(221)
    plot.plot(ts,[xs[i][0] for i in range(len(xs))])
    plot.subplot(222)
    plot.plot(ts, [xs[i][1] for i in range(len(xs))])
    plot.subplot(223)
    # plot.plot(ts, [xs[i][2] for i in range(len(xs))])
    # plot.subplot(224)
    plot.plot(ts,norms)
    plot.show()
    for i in range(len(xs)):
        print(ts[i],xs[i],norms[i])
    pass








#第一种求解器
    # # yjl=-mA(t)'*inv(mA(t)*mA(t)')*dA*x(1:3)-gamma*mA(t)'*inv(mA(t)*mA(t)')*(mA(t)*x(1:3)-vB(t))+mA(t)'*inv(mA(t)*mA(t)')*dB+mA(t)'*inv(mA(t)*mA(t)')*r;%OZNN（不带积分，不带激励函数）
    #
    # tspan=(0,10)
    # ts,xs= ode.backwardeuler(OdeRightFunc,x0,tspan,(tspan[-1]-tspan[0])/1000)
    #
    # # for i in range(len(ts)):
    # #     print(ts[i])
    # #     print(xs[0][i],"  ",xs[1][i],"  ",xs[2][i],"  ")
    #
    # norms=[sp.Matrix(Mt_value(ts[i])*sp.Matrix([xs[0][i],xs[1][i],xs[2][i]])-Bt_value(ts[i])).norm() for i in  range(len(ts))]
    # # for i in range(len(ts)):
    # #     print(ts[i])
    # #     print(sp.Matrix(Mt_value(ts[i])*sp.Matrix([xs[0][i],xs[1][i],xs[2][i]])-Bt_value(ts[i])).norm())
    # plot.figure(1)
    # # plot.subplot(221)
    # # plot.plot(ts,xs[0])
    # # plot.subplot(222)
    # # plot.plot(ts, xs[1])
    # # plot.subplot(223)
    # # plot.plot(ts, xs[2])
    # # plot.subplot(224)
    # plot.plot(ts,norms)
    # plot.show()
    # for i in range(len(ts)):
    #     print(ts[i],"    "*4,norms[i])



    #
    # """解决方法2"""