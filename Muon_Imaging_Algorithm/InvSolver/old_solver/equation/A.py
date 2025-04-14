# author:高金磊
# datetime:2021/11/1 10:36

from solver.old_solver.equation.MA import M,sp,accuracy,t,Mt_value,m,n
from solver.old_solver.equation.L import a,B
import numpy as np

H = sp.Matrix([[-sp.eye(n)], [sp.eye(n)]])

_A=None
def At(ct=None):
    global _A
    if _A is not None:
        return _A#懒加载
    _A = sp.Matrix(M)
    _A = _A.col_join(H)
    O = sp.Matrix(sp.zeros(m, 2 * n))
    if ct is None:
        ct: sp.Matrix = np.power(a - H * M.pinv() * B, 0.5)
    Dt=sp.eye(len(ct))
    for i in range(len(ct)):
        Dt[i,i]=ct[i]
    OD = O.col_join(Dt)
    _A = _A.row_join(OD)

    return _A


def At_value(t_value,ct=None):


    A = sp.Matrix(Mt_value(t_value))
    A = A.col_join(H)
    O = sp.Matrix(sp.zeros(m, 2 * n))
    #为了加快速度,提前带入计算---ct的来源不明确-可以直接通过x带入
    if ct is None:
        ct=np.power(a - H * Mt_value(t_value).pinv() * B, 0.5)
    Dt = sp.eye(len(ct))
    for i in range(len(ct)):
        Dt[i, i] = ct[i]
    OD = O.col_join(Dt)
    A = A.row_join(OD)

    return A.subs(t, t_value).evalf(accuracy, chop=True)

_At_diff=None
def At_diff(ct=None):
    global _At_diff
    if _At_diff is None:
        _At_diff=At(ct).diff()
    return _At_diff

def At_diff_value(t_value,ct=None):
    return At_diff(ct).subs(t, t_value).evalf(accuracy, chop=True)


if __name__ == '__main__':
    print(At(ct=[1,2,3,4,5,6]))
    print(At().shape)
