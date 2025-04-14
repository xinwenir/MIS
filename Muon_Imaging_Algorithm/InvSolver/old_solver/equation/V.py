# author:高金磊
# datetime:2021/11/1 20:59
import sympy as sp

_middle=sp.Matrix(-Mt_diff())
_middle=_middle.col_join(sp.zeros(2*n,n))
V=_middle.row_join(sp.zeros((m+2*n),2*n))
def Vt_value(t_value):
    return V.subs(t, t_value).evalf(accuracy, chop=True)

if __name__ == '__main__':
    print(V.shape)
    print(V)