# author:高金磊
# datetime:2021/11/1 16:36

from solver.old_solver.equation.MA import M,sp,accuracy,t,Mt_value,m,n
from solver.old_solver.equation.L import a,B
import numpy as np

H = sp.Matrix([[-sp.eye(n)], [sp.eye(n)]])
_K=None
def Kt(ct=None):
    global _K
    if _K is not None:
        return _K#懒加载

    _K = sp.Matrix(M)
    _K = _K.col_join(H)
    O = sp.Matrix(sp.zeros(m, 2 * n))
    if ct is None:
        ct: sp.Matrix = np.power(a - H * M.pinv() * B, 0.5)
    Dt=sp.eye(len(ct))
    for i in range(len(ct)):
        Dt[i,i]=2*ct[i]
    OD = O.col_join(Dt)
    _K = _K.row_join(OD)

    return _K


def Kt_value(t_value,ct=None):

    K = sp.Matrix(Mt_value(t_value))
    K = K.col_join(H)
    O = sp.Matrix(sp.zeros(m, 2 * n))
    #为了加快速度,提前带入计算---ct的来源不明确-可以直接通过x带入
    if ct is None:
        ct=np.power(a - H * Mt_value(t_value).pinv() * B, 0.5)
    Dt = sp.eye(len(ct))
    for i in range(len(ct)):
        Dt[i, i] = 2*ct[i]
    OD = O.col_join(Dt)
    K= K.row_join(OD)

    return K.subs(t, t_value).evalf(accuracy, chop=True)

_Kt_diff=None
def Kt_diff(ct=None):
    global _Kt_diff
    if _Kt_diff is None:
        _Kt_diff=Kt(ct).diff()
    return _Kt_diff

def Kt_diff_value(t_value,ct=None):
    return Kt_diff(ct).subs(t, t_value).evalf(accuracy, chop=True)

if __name__ == '__main__':
    print(Kt(ct=[1,2,3,4,5,6]))
