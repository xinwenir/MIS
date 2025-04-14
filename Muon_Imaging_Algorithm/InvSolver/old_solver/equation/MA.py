# author:高金磊
# datetime:2021/10/26 13:43

from solver.old_solver.Setting import accuracy,M,t


# M=sp.Matrix([[2*t,3+t,t]])
m,n=M.shape
def Mt():
    return M
def Mt_value(t_value):
    return M.subs(t,t_value).evalf(accuracy,chop=True)
    # return sp.Matrix(Matrix_eval(M, t, t_value, accuracy))

_Mt_diff=None
def Mt_diff():
    global _Mt_diff
    if _Mt_diff is None:
        _Mt_diff=M.diff()
    return _Mt_diff
    # return Matrix_diff(M,t)
def Mt_diff_value(t_value):
    # return sp.Matrix(Matrix_eval(Mt_diff(), t, t_value, accuracy))
    return Mt_diff().subs(t,t_value).evalf(accuracy,chop=True)

if __name__ == '__main__':
    pass