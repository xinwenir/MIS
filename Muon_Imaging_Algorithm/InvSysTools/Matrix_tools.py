# author:高金磊
# datetime:2021/10/26 19:29
import numpy as np
import sympy as sp



#以下方法效率低于官方这里不再使用
def Matrix_eval(matrix,symbol_name,symbol_value:int,accuracy=10):
    """
     计算矩阵表达式的值
    :param: matrix: 待求解的矩阵
    :param:symbol_name: 变量名--目前本函数仅支持一个
    :type: sympy.core.symbol.Symbol
    :param: symbol_value: 变量值
    :param: accuracy: 求解精确度
    :return: 返回求解后的矩阵的
    :rtype:
    """
    if len(matrix)==0:
        return matrix
    res=[]
    for i in range(0,len(matrix)):
        arr=[]
        for j in range(0,len(matrix[0])):
            # arr.append(float(sp.latex(matrix[i][j].subs([(symbol_name,symbol_value)]).evalf(accuracy,chop=True))))
            arr.append(matrix[i][j].subs([(symbol_name,symbol_value)]).evalf(accuracy,chop=True))
        res.append(arr)
    return res
def Matrix_diff(matrix,symbol_name):
    """
    对矩阵进行求导
    :param matrix:
    :param symbol_name:
    :return:
    """
    if len(matrix)==0:
        return matrix
    res = []
    for i in range(len(matrix)):
        arr=[]
        for j in range(len(matrix[0])):
            arr.append(sp.diff(matrix[i][j],symbol_name))
        res.append(arr)
    return res