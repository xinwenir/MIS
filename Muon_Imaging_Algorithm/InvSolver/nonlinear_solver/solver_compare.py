def fun(x):
    return (x-2)**2+(x+1)**2
def grad(x):
    return ((x-2)*2)+(x+1)*2

if __name__ == '__main__':
    from scipy.optimize import *
    # x0=100
    # middle=0
    # while(abs(x0-middle)>0.15):
    #     middle=x0
    #     res=line_search(fun,grad,x0,pk=-0.2)
    #     x0-=0.2*res[0]
    # print(x0)
    x0=100
    res=minimize(fun,x0,method="BFGS",jac=grad)
    print(res)
