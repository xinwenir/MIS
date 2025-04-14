# author:高金磊
# datetme:2021/10/28 20:22

from solver.old_solver.equation.A import *
from solver.old_solver.equation.L import *
from solver.old_solver.equation.K import *

import sympy as sp
from solver.old_solver.AFM import *


def OdeRightFunc(x0,t,gamma,lmd,AFM1=Bound,AFM2=linear):
    print(t)
    # print(x)
    om = sp.Matrix(x0[3 * n:])
    x=sp.Matrix(x0[:3*n])
    ns=sp.Matrix([1]*(m+2*n))
    nt=sp.Matrix([t]*(m+2*n))
    # ct_value=sp.Matrix(np.power(a-H*sp.Matrix(x[:n]),0.5))
    # ct_value=[1]*(2*n)
    ct_value=sp.Matrix(x0[n:3*n])

    # Kt_V=sp.Matrix([[2 * t, t + 3, 0, 0, 0, 0], [-1, 0, 2*x[2], 0, 0, 0], [0, -1, 0, 2*x[3], 0, 0], [1, 0, 0, 0, 2*x[4], 0], [0, 1, 0, 0, 0, 2*x[5]]])
    # Kt_V=sp.Matrix([[2 * t, t + 3, 0, 0, 0, 0], [-1, 0, 2, 0, 0, 0], [0, -1, 0, 2, 0, 0], [1, 0, 0, 0, 2, 0], [0, 1, 0, 0, 0, 2]])
    # At_V=sp.Matrix([[2 * t, t + 3, 0, 0, 0, 0], [-1, 0, x[2], 0, 0, 0], [0, -1, 0, x[3], 0, 0], [1, 0, 0, 0, x[4], 0], [0, 1, 0, 0, 0, x[5]]])
    # At_V=sp.Matrix([[2 * t, t + 3, 0, 0, 0, 0], [-1, 0, 1, 0, 0, 0], [0, -1, 0, 1, 0, 0], [1, 0, 0, 0, 1, 0], [0, 1, 0, 0, 0, 1]])
    # Vt_V=sp.Matrix([[-2, -1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]])
    # Lt_V=sp.Matrix([[t], [-0.250000000000000], [-0.250000000000000], [0.250000000000000], [0.250000000000000]])
    # dLt_V=sp.Matrix([[1], [0], [0], [0], [0]])
    # ns = sp.Matrix([0, 1])
    # nt = sp.Matrix([4*t,4*t])

    # f=Kt_V.pinv()*(
    #     Vt_V*x+dLt_V-lmd*AFM1(At_V*x-Lt_V)
    # )
    # f=Kt_V.pinv()*(
    #     Vt_V*x+dLt_V-lmd*AFM1(At_V*x-Lt_V)
    #     -gamma*AFM2(At_V*x-Lt_V+lmd*om)
    # )
    # Kt_V=sp.Matrix([[sin(t) + 2, cos(t) + 3, -sin(2*t) + cos(t) + 4, 0, 0, 0, 0, 0, 0],
    #                 [-1, 0, 0, ct_value[0]*2, 0, 0, 0, 0, 0],
    #                 [0, -1, 0, 0, ct_value[1]*2, 0, 0, 0, 0],
    #                 [0, 0, -1, 0, 0, ct_value[2]*2, 0, 0, 0],
    #                 [1, 0, 0, 0, 0, 0, ct_value[3]*2, 0, 0],
    #                 [0, 1, 0, 0, 0, 0, 0, ct_value[4]*2, 0],
    #                 [0, 0, 1, 0, 0, 0, 0, 0,ct_value[5]*2]])
    # At_V=sp.Matrix([[sin(t) + 2, cos(t) + 3, -sin(2*t) + cos(t) + 4, 0, 0, 0, 0, 0, 0],
    #                 [-1, 0, 0, ct_value[0], 0, 0, 0, 0, 0],
    #                 [0, -1, 0, 0, ct_value[1], 0, 0, 0, 0],
    #                 [0, 0, -1, 0, 0, ct_value[2], 0, 0, 0],
    #                 [1, 0, 0, 0, 0, 0, ct_value[3], 0, 0],
    #                 [0, 1, 0, 0, 0, 0, 0, ct_value[4], 0],
    #                 [0, 0, 1, 0, 0, 0, 0, 0,ct_value[5]]])
    # Vt_V=sp.Matrix([[-cos(t), sin(t), sin(t) + 2*cos(2*t), 0, 0, 0, 0, 0, 0],
    #                 [0, 0, 0, 0, 0, 0, 0, 0, 0],
    #                 [0, 0, 0, 0, 0, 0, 0, 0, 0],
    #                 [0, 0, 0, 0, 0, 0, 0, 0, 0],
    #                 [0, 0, 0, 0, 0, 0, 0, 0, 0],
    #                 [0, 0, 0, 0, 0, 0, 0, 0, 0],
    #                 [0, 0, 0, 0, 0, 0, 0, 0, 0]])
    # Lt_V=sp.Matrix([[sin(2*t) + cos(t)], [1], [0], [0], [1.0000000000000], [0.10000000000000], [0.10000000000000]])
    # dLt_V=sp.Matrix([[-sin(t) + 2*cos(2*t)], [0], [0], [0], [0], [0], [0]])
    Vt_V=Vt_value(t)
    At_V=At_value(t, ct=ct_value)
    Lt_V=lt_value(t)
    dLt_V=lt_diff_value(t)
    Kt_V=Kt_value(t,ct_value)
    f=Kt_V.pinv()*(
        Vt_V*x+dLt_V
        -lmd*AFM1(At_V*x-Lt_V)
        -gamma*AFM2(At_V*x-Lt_V+lmd*AFM1(om))

    )
    # f = Kt_value(t,ct_value).pinv() * (
    #         Vt_value(t) * x + lt_diff_value(t)
    #         - lmd * AFM1(At_value(t,ct_value) * x - lt_value(t))
    #         - gamma * AFM2(At_value(t,ct_value) * x - lt_value(t) + lmd * AFM1(om))
    # )


    om=At_V*x-Lt_V
    # om=At_value(t,ct=ct_value)*x-lt_value(t)
    res=[f[i] for i in range(len(f))]+[om[i] for i in range(len(om))]
    # print(*f)
    return res
if __name__ == '__main__':
    ct=sp.Matrix([1,1,1,1,1,1]).T
    middle=At_value(1,ct).pinv()*lt_value(1)
    print(middle)
