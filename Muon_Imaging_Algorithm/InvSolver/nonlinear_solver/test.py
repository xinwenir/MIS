# author:高金磊
# datetime:2022/7/29 21:48
import datetime

# 函数表达式fun
import time
from copy import copy

from numpy import empty
from urllib3.connectionpool import xrange
import numpy as np

from InvSysTools.tools import dot

# fun = lambda x: 100 * (x[0] ** 2 - x[1]) ** 2 + (x[0] - 1) ** 2
#
# # 梯度向量 gfun
# gfun = lambda x: np.array([400 * x[0] * (x[0] ** 2 - x[1]) + 2 * (x[0] - 1), -200 * (x[0] ** 2 - x[1])])
#
# # 海森矩阵 hess
# hess = lambda x: np.array([[1200 * x[0] ** 2 - 400 * x[1] + 2, -400 * x[0]], [-400 * x[0], 200]])


def gradient(fun, gfun, x0):
    # 最速下降法求解无约束问题
    # x0是初始点，fun和gfun分别是目标函数和梯度
    maxk = 5000
    rho = 0.5
    sigma = 0.4
    k = 0
    epsilon = 1e-5

    while k < maxk:
        gk = gfun(x0)
        dk = -gk
        if np.linalg.norm(dk) < epsilon:
            break
        m = 0
        mk = 0
        while m < 20:
            if fun(x0 + rho ** m * dk) < fun(x0) + sigma * rho ** m * np.dot(gk, dk):
                mk = m
                break
            m += 1
        x0 += rho ** mk * dk
        k += 1

    return x0, fun(x0), k


def dampnm(fun, gfun, hess, x0):
    # 用阻尼牛顿法求解无约束问题
    # x0是初始点，fun，gfun和hess分别是目标函数值，梯度，海森矩阵的函数
    maxk = 500
    rho = 0.55
    sigma = 0.4
    k = 0
    epsilon = 1e-5
    while k < maxk:
        gk = gfun(x0)
        Gk = hess(x0)
        dk = -1.0 * np.linalg.solve(Gk, gk)
        if np.linalg.norm(dk) < epsilon:
            break
        m = 0
        mk = 0
        while m < 20:
            if fun(x0 + rho ** m * dk) < fun(x0) + sigma * rho ** m * np.dot(gk, dk):
                mk = m
                break
            m += 1
        x0 += rho ** mk * dk
        k += 1

    return x0, fun(x0), k


def dfp(fun, gfun, hess, x0):
    # 功能：用DFP族算法求解无约束问题：min fun(x)
    # 输入：x0是初始点，fun,gfun分别是目标函数和梯度
    # 输出：x,val分别是近似最优点和最优解,k是迭代次数
    maxk = 1e5
    rho = 0.55
    sigma = 0.4
    epsilon = 1e-5
    k = 0
    n = np.shape(x0)[0]
    # 海森矩阵可以初始化为单位矩阵
    Hk = np.linalg.inv(hess(x0))  # 或者单位矩阵np.eye(n)

    while k < maxk:
        gk = gfun(x0)
        if np.linalg.norm(gk) < epsilon:
            break
        dk = -1.0 * np.dot(Hk, gk)

        m = 0
        mk = 0
        while m < 20:  # 用Armijo搜索求步长
            if fun(x0 + rho ** m * dk) < fun(x0) + sigma * rho ** m * np.dot(gk, dk):
                mk = m
                break
            m += 1
        # print mk
        # DFP校正
        x = x0 + rho ** mk * dk
        sk = x - x0
        yk = gfun(x) - gk

        if np.dot(sk, yk) > 0:
            Hy = np.dot(Hk, yk)
            print(Hy)
            sy = np.dot(sk, yk)  # 向量的点积
            yHy = np.dot(np.dot(yk, Hk), yk)  # yHy是标量
            # 表达式Hy.reshape((n,1))*Hy 中Hy是向量，生成矩阵
            Hk = Hk - 1.0 * Hy.reshape((n, 1)) * Hy / yHy + 1.0 * sk.reshape((n, 1)) * sk / sy

        k += 1
        x0 = x

    return x0, fun(x0), k  # 分别是最优点坐标，最优值，迭代次数


def bfgs(fun, gfun, hess, x0):
    # 功能：用BFGS族算法求解无约束问题：min fun(x)
    # 输入：x0是初始点，fun,gfun分别是目标函数和梯度
    # 输出：x,val分别是近似最优点和最优解,k是迭代次数
    maxk = 1e5
    rho = 0.55
    sigma = 0.4
    epsilon = 1e-5
    k = 0
    n = np.shape(x0)[0]
    # 海森矩阵可以初始化为单位矩阵
    Bk = eye(n)  # np.linalg.inv(hess(x0)) #或者单位矩阵np.eye(n)

    while k < maxk:
        gk = gfun(x0)
        if np.linalg.norm(gk) < epsilon:
            break
        dk = -1.0 * np.linalg.solve(Bk, gk)
        m = 0
        mk = 0
        while m < 20:  # 用Armijo搜索求步长
            if fun(x0 + rho ** m * dk) < fun(x0) + sigma * rho ** m * np.dot(gk, dk):
                mk = m
                break
            m += 1

        # BFGS校正
        x = x0 + rho ** mk * dk
        sk = x - x0
        yk = gfun(x) - gk

        if np.dot(sk, yk) > 0:
            Bs = np.dot(Bk, sk)
            ys = np.dot(yk, sk)
            sBs = np.dot(np.dot(sk, Bk), sk)

            Bk = Bk - 1.0 * Bs.reshape((n, 1)) * Bs / sBs + 1.0 * yk.reshape((n, 1)) * yk / ys

        k += 1
        x0 = x

    return x0, fun(x0), k  # 分别是最优点坐标，最优值，迭代次数


def broyden(fun, gfun, hess, x0):
    # 功能：用Broyden族算法求解无约束问题：min fun(x)
    # 输入：x0是初始点，fun,gfun分别是目标函数和梯度
    # 输出：x,val分别是近似最优点和最优解,k是迭代次数
    x0 = np.array(x0)

    maxk = 1e5
    rho = 0.55
    sigma = 0.4
    epsilon = 1e-5
    phi = 0.5
    k = 0
    n = np.shape(x0)[0]

    Hk = np.linalg.inv(hess(x0))

    while k < maxk:
        gk = gfun(x0)
        if np.linalg.norm(gk) < epsilon:
            break
        dk = -1 * np.dot(Hk, gk)

        m = 0
        mk = 0
        while m < 20:  # 用Armijo搜索求步长
            if fun(x0 + rho ** m * dk) < fun(x0) + sigma * rho ** m * np.dot(gk, dk):
                mk = m
                break
            m += 1
        # Broyden族校正
        x = x0 + rho ** mk * dk
        sk = x - x0
        yk = gfun(x) - gk

        Hy = np.dot(Hk, yk)
        sy = np.dot(sk, yk)
        yHy = np.dot(np.dot(yk, Hk), yk)

        if (sy < 0.2 * yHy):
            theta = 0.8 * yHy / (yHy - sy)
            sk = theta * sk + (1 - theta) * Hy
            sy = 0.2 * yHy

        vk = np.sqrt(yHy) * (sk / sy - Hy / yHy)
        Hk = Hk - Hy.reshape((n, 1)) * Hy / yHy + sk.reshape((n, 1)) * sk / sy + phi * vk.reshape((n, 1)) * vk
        k += 1
        x0 = x
    return x0, fun(x0), k  # 分别是最优点坐标，最优值，迭代次数


"""LBFGS"""


def twoloop(s, y, rho, gk):
    n = len(s)  # 向量序列的长度

    if np.shape(s)[0] >= 1:
        # h0是标量，而非矩阵
        h0 = 1.0 * np.dot(s[-1], y[-1]) / np.dot(y[-1], y[-1])
    else:
        h0 = 1

    a = empty((n,))

    q = gk.copy()
    for i in range(n - 1, -1, -1):
        a[i] = rho[i] * dot(s[i], q)
        q -= a[i] * y[i]
    z = h0 * q

    for i in range(n):
        b = rho[i] * dot(y[i], z)
        z += s[i] * (a[i] - b)

    return z


def lbfgs(fun, gfun, x0, m=7,args=None,step_strategy="Armijo",details=False):
    # fun和gfun分别是目标函数及其一阶导数,x0是初值,m为储存的序列的大小
    step_strategy_dict={"Armijo":Armijo,"wolfe":wolfe}
    if args is None:
        args=()
    else:
        if type(args) is not tuple:
            args=(args,)
    maxk = 1000
    epsilon = 1e-10
    k = 0
    # n = np.shape(x0)[0]  # 自变量的维度

    s, y, rho = [], [], []

    while k < maxk:
        gk = gfun(x0,*args)
        if np.linalg.norm(gk) < epsilon:
            break

        dk = -1.0 * twoloop(s, y, rho, gk)

        step_len = step_strategy_dict[step_strategy](fun, x0, dk, gk,details, args)


        x = x0 +  step_len* dk
        sk = x - x0
        yk = gfun(x,*args) - gk

        if np.dot(sk, yk) > 0:  # 增加新的向量
            rho.append(1.0 / np.dot(sk, yk))
            s.append(sk)
            y.append(yk)
        if np.shape(rho)[0] > m:  # 弃掉最旧向量
            rho.pop(0)
            s.pop(0)
            y.pop(0)

        k += 1
        x0 = x
        # print(x)
        # print(fun(x))

    return x0, fun(x0,*args), k  # 分别是最优点坐标，最优值，迭代次数

def Armijo(fun,x0,dk,gk,details,args):
    start_time=time.time()
    rou = 0.55
    sigma = 0.4
    m0 = 0  # 不明白为啥___0的话while循环不就没变化
    mk = m0
    old_value = fun(x0, *args)
    new_value = old_value
    flag = True
    while m0 < 5:  # 用Armijo搜索求步长
        if flag:
            flag = False
        else:
            new_value = fun(x0 + rou ** m0 * dk, *args)
        if new_value < old_value + sigma * rou ** m0 * np.dot(gk, dk):
            mk = m0
            break
        m0 += 1
    if details:
        print("""""""Armijio执行结果""""""")
        print("用时:",time.time()-start_time)
        print("步长:",rou**mk)
    return rou ** mk
def wolfe(fun,x,dk,gk,gfun,details,args):
    start_time=time.time()
    rho = 0.25
    sigma = 0.75
    alpha = 1
    a = 0
    b = 1
    old_value = fun(x, *args)
    while a < b:
        if fun(x + alpha*dk,*args) <= old_value + rho*alpha*np.dot(dk,gk):
            if np.dot(gfun(x + alpha * dk,*args),gk) >= sigma*np.dot(dk,gk):
                break
            else:
                a = alpha
                alpha = 0.5 * (a + b)
        else:
            b = alpha
            alpha = (alpha + a) / 2
    if details:
        print("""""""wolfe执行结果""""""")
        print("用时:",time.time()-start_time)
        print("步长:",alpha)
    return alpha
def BCD_line_search(x0,num,fun, gfun,args=None):
    if args is None:
        args = ()
    else:
        if type(args) is not tuple:
            args = (args,)
    steps=[1,1e-1,1e-3,1e-5,1e-7]

    while len(steps)>0:
        step=steps.pop(0)
        gun = gfun(x0)
        old_fun = fun(x0)
        while True:

            x0[num]-=gun[num]*step
            new_fun=fun(x0)
            if new_fun<old_fun:
                old_fun=new_fun
            else:
                x0[num]+=gun[num]*step
                break
    return x0[num]
def BCD_algorithm(x0,fun, gfun,args=None,optimization_method=BCD_line_search):
    if args is None:
        args = ()
    else:
        if type(args) is not tuple:
            args = (args,)
    old_fun=fun(x0,*args)
    count=0
    while True:
        for i in range(len(x0)):
            old_value=x0[i]
            new_value=optimization_method(x0,i,fun,gfun,args)
            x0[i]=new_value
        new_fun=fun(x0)
        # if count==40:
        #     break
        # count+=1
        # old_fun = new_fun
        print(count,x0,new_fun)
        if new_fun>=old_fun-0.000001:
            break
        else:
            old_fun=new_fun
    return x0,new_fun

def Rastrigin_fun(xy):
    X = xy[0]
    Y = xy[1]
    return 20 + X * X - 10 * np.cos(X * np.pi) + Y ** 2 - 10 * np.cos(2 * np.pi * Y)


def Rastrigin_gfun(xy):
    X = xy[0]
    Y = xy[1]
    return np.array([2 * X + 10 * np.pi * np.sin(X * np.pi), 2 * Y - 20 * np.pi * np.sin(2 * np.pi * Y)])


def Sphere_fun(xy):
    X = xy[0]
    Y = xy[1]
    return X ** 2 + Y ** 2


def Sphere_gfun(xy):
    X = xy[0]
    Y = xy[1]
    return np.array([2 * X ,  2 * Y])


def my_fun(xy):
    X = xy[0]
    Y = xy[1]
    return (X - 1) ** 4 + Y ** 2


def my_gfun(xy):
    X = xy[0]
    Y = xy[1]
    return np.array([4 * (X - 1) ** 3, 2 * Y])
def my_hessfun(xy):
    X = xy[0]
    Y = xy[1]
    return np.array([[12 * (X - 1) ** 2, 0],
                     [0,2]])


if __name__ == '__main__':
    x0 = [-200.5, 10.8]
    fun=Sphere_fun
    gun=Sphere_gfun
    print("....................python-lbfgs.............................")

    from scipy.optimize import minimize, dual_annealing

    res = minimize(method="l-bfgs-b", fun=fun, x0=np.array(x0), jac=gun)
    print(res)
    print("...................双重退火..............................")
    res = dual_annealing(fun, x0=np.array(x0), bounds=[[-100, 100], [-1000, 1000]])
    print(res)
    print("....................my-lbfgs.............................")
    a = lbfgs(fun, gun, x0)
    res=copy(res)
    res.x=a[0]
    res.fun=a[1]
    res.nfev=a[2]
    print(res)
    print("....................分块坐标下降.............................")
    res=BCD_algorithm(x0,fun=fun,gfun=gun)
    print(res)
