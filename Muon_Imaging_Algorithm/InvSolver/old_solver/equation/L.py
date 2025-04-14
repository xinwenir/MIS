# author:高金磊
# datetime:2021/11/1 10:36
from solver.old_solver.equation.VB import accuracy,t
from solver.old_solver.Setting import L


def lt():
    return L


def lt_value(t_value):
    return L.subs(t, t_value).evalf(accuracy, chop=True)


def lt_diff():
    return L.diff()


def lt_diff_value(t_value):
    return lt_diff().subs(t, t_value).evalf(accuracy, chop=True)

if __name__ == '__main__':
    print(lt_diff_value(1))
    print(L.shape)
    print(L)
    print(lt_diff())
    print(lt_value(0))