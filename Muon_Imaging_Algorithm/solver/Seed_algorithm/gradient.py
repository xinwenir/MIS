# author:高金磊
# datetime:2022/9/1 9:48
import numpy


def jacbi(x, A):
    jac = A.T * numpy.dot(A,x)
    return jac