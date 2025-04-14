#由于考虑平滑性效率低下而且beta不好调节，这里提供了一种只求解决misfit的优化方案
def solver(A,x0,b,bounds):
    import numpy as np
    from scipy.optimize import minimize
    # 定义目标函数
    def objective_function(x):
        Ax = A.dot(x)
        return np.linalg.norm(Ax - b)**2  # 目标函数：||Ax - b||^2

    # 定义雅可比矩阵函数
    def jac_function(x):
        return 2 * A.transpose().dot(A.dot(x) - b)  # 雅可比矩阵：2A^T(Ax - b)
    
    # for i in range(len(bounds)):
    #     bounds[i]=[0,2.7]

    # 调用LBFGSB算法进行优化
    result = minimize(objective_function, x0, jac=jac_function, bounds=bounds,method='L-BFGS-B',options={"maxiter": 300,'disp': True})
    return result