# author:高金磊
# datetme:2021/10/28 20:22
# from solver.equation.MA import *
# from solver.equation.VB import *
from solver.old_solver.Setting import M,B,Mp
from solver.old_solver.AFM import *


def OdeRightFunc(x,t,gamma,lmd):
    print(t)
    # om=sp.Matrix(x[M.cols:])
    # om=sp.Matrix(x[M.shape[1]:])
    # x=sp.Matrix(x[0:M.shape[1]])
    om=[[i] for i in x[M.shape[1]:]]
    x=[[i] for i in x[0:M.shape[1]]]
    # AFM1=powersum
    # AFM2=powersigmoid
    # n = sp.Matrix([0, 1])
    # nt = sp.Matrix([4*t,4*t])
    # f = (-Mt_value(t).pinv() * Mt_diff_value(t) * x+ Mt_value(t).pinv() * Bt_diff_value(t)- lmd * Mt_value(t).pinv() *
    #      (Mt_value(t) * x - Bt_value(t))-gamma*Mt_value(t).pinv() * (Mt_value(t)*x-Bt_value(t)+lmd*om))

    # ##L L 时变噪声
    # f = Mt_value(t).pinv() * (
    #             -Mt_diff_value(t) * x + Bt_diff_value(t) - lmd * (Mt_value(t) * x - Bt_value(t)) - gamma * (
    #                 Mt_value(t) * x - Bt_value(t) + lmd * om) + n)
    # f = Mt_value(t).pinv() * (
    #             -Mt_diff_value(t) * x + Bt_diff_value(t) - lmd * AFM1(Mt_value(t) * x - Bt_value(t)) - gamma * AFM2(
    #                 Mt_value(t) * x - Bt_value(t) + lmd *AFM1(om)))
    f = Mp* (
                -M * x + B - lmd * (M * x - B) - gamma * (
                   M * x - B + lmd *(om)))

    #先计算表达式再计算速度很慢
    # ti = sp.symbols('t')
    # f = (M.pinv() * (-Mt_diff() * x+ Bt_diff()- lmd * (Mt() * x - Bt())-gamma*(Mt()*x-Bt())+lmd*om)).subs(ti,t).evalf(accuracy,chop=True)

    # om=Mt_value(t)*x-Bt_value(t)
    om=M*x-B
    print(np.linalg.norm(om))
    res=[f[i,0] for i in range(len(f))]+[om[i,0] for i in range(len(om))]
    return res

#适用于第一种ode求解器的
# def OdeRightFunc(t,x):
#     print(t)
#     om=sp.Matrix(x[3:5])
#     x=sp.Matrix(x[0:3])
#     # f = (-Mt_value(t).pinv() * Mt_diff_value(t) * x+ Mt_value(t).pinv() * Bt_diff_value(t)- lmd * Mt_value(t).pinv() *
#     #      (Mt_value(t) * x - Bt_value(t))-gamma*Mt_value(t).pinv() * (Mt_value(t)*x-Bt_value(t)+lmd*om))
#     f = Mt_value(t).pinv() * (-Mt_diff_value(t) * x+ Bt_diff_value(t)- lmd * (Mt_value(t) * x - Bt_value(t))-gamma*(Mt_value(t)*x-Bt_value(t)+lmd*om))
#
#
#     #先计算表达式再计算速度很慢
#     # ti = sp.symbols('t')
#     # f = (M.pinv() * (-Mt_diff() * x+ Bt_diff()- lmd * (Mt() * x - Bt())-gamma*(Mt()*x-Bt())+lmd*om)).subs(ti,t).evalf(accuracy,chop=True)
#
#     om=Mt_value(t)*x-Bt_value(t)
#     res=[f[i] for i in range(len(f))]+[om[i] for i in range(len(om))]
#     return res
