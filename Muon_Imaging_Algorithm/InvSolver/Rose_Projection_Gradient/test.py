# author:高金磊
# datetime:2022/6/15 10:56
import math

import numpy as np
from scipy.optimize import minimize

min=2
max=100
ykxn=0.00000000001
deta=0.45
def fun(x):
    print("x:",x)
    x1 = x[0]
    x2 = x[1]
    if x1<=min:
        x1=min+ykxn
    if x1>=max:
        x1=max-ykxn
    if x2<=min:
        x2=min+ykxn
    if x2>=max:
        x2=max-ykxn

    # x[0]=0.1
    ##罚函数
    B=deta*math.log2(100-x1)+deta*math.log2(x1-2)+deta*math.log2(100-x2)+deta*math.log2(x2-2)
    x1 = x[0]
    x2 = x[1]
    res= x1 * x1 + \
             x2 * x2 \
           +2 * x1-B
    print("fun:",res)
    return res


def jac(x):
    x1 = x[0]
    x2 = x[1]
    if x1<=min:
        x1=min+ykxn
    if x1>=max:
        x1=max-ykxn
    if x2<=min:
        x2=min+ykxn
    if x2>=max:
        x2=max-ykxn
    # x[0] = 0.001
    B_jac=[deta*(-1/(100-x1)+1/(x1-2)),deta*(-1/(100-x2)+1/(x2-2))]

    x1 = x[0]
    x2 = x[1]

    res=[
        2 * x1 +2-B_jac[0],
        2*x2-B_jac[1]
    ]
    print("grad:",res)
    return res


if __name__ == '__main__':
    x0 = np.array([[50], [50]])
    # bounds=[[2,100],[2,100]]
    result=minimize(fun=fun, jac=jac, x0=x0, method="CG")
    print(result.x)
    print(result.fun)
