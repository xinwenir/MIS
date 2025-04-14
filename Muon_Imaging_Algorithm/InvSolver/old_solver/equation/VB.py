# author:高金磊
# datetime:2021/10/26 13:43


def Bt():
    return B
def Bt_value(t_value):
    return B.subs(t,t_value).evalf(accuracy,chop=True)
def Bt_diff():
    return B.diff()
def Bt_diff_value(t_value):
    return Bt_diff().subs(t,t_value).evalf(accuracy,chop=True)


# def Bt():
#     return B
# def Bt_value(t_value):
#     return sp.Matrix(Matrix_eval(B,t,t_value))
# def Bt_diff():
#     return Matrix_diff(B,t)
# def Bt_diff_value(t_value):
#     return sp.Matrix(Matrix_eval(Matrix_diff(B,t),t,t_value))