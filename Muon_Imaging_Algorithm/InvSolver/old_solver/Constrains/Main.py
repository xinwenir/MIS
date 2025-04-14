# author:高金磊
# datetime:2021/10/26 20:55

import matplotlib.pyplot as plot
import scipy as sc
from scipy import integrate

from solver.old_solver.Constrains.OdeRight import OdeRightFunc
from solver.old_solver.equation.A import *
from solver.old_solver.equation.L import *
from solver.old_solver.AFM import *

# x0=[(0.1*random.random()) for i in range(3*n)]+[0]*(m+2*n)
x0=[0.1 for i in range(3*n)]+[0]*(m+2*n)
# x0=[random.random() for i in range(5*n+m)]
# x0=[0]*(5*n+m)
if __name__ == '__main__':
    ts = np.linspace(0, 10, 100)
    data = sc.integrate.odeint(OdeRightFunc, x0, ts,(gamma,lmd,Bound,linear))
    xs=[]
    norms=[]
    for i in range(len(data)):
        xs.append(data[i][:3*n])
        ct=data[i][n:3*n]
        norms.append(sp.Matrix(At_value(ts[i],ct)*sp.Matrix(data[i][:3*n])-lt_value(ts[i])).norm())
    plot.plot(ts,norms)
    plot.show()
    plot.figure(1)
    plot.subplot(221)
    plot.axhline(-a[0])
    plot.axhline(a[3])
    plot.plot(ts,[xs[i][0] for i in range(len(xs))])
    plot.subplot(222)
    plot.axhline(-a[1])
    plot.axhline(a[4])
    plot.plot(ts, [xs[i][1] for i in range(len(xs))])
    plot.subplot(223)
    plot.axhline(-a[2])
    plot.axhline(a[5])
    plot.plot(ts, [xs[i][2] for i in range(len(xs))])
    plot.subplot(224)
    plot.plot(ts,norms)
    plot.show()
    for i in range(len(xs)):
        # print(ts[i],*xs[i][:3],norms[i])
        print(ts[i], *xs[i][:3])
        # print(xs[i])
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