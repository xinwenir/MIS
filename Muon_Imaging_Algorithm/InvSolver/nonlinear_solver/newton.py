import numpy as np
import matplotlib.pyplot as plt
import time

N = 1

# f = lambda x: np.power(x[:, 0], 2) + np.sin(np.pi * x[:, 1])
#
#
# cost = lambda x: np.sum(f(x))
#
#
# def grad(x):
#     gx0 = 2 * x[:, 0]
#     gx1 = np.pi * np.cos(np.pi * x[:, 1])
#     return np.array([np.sum(gx0), np.sum(gx1)]).reshape((2, 1))

#
# def hessian(x):
#     h00 = 2
#     h01 = 0
#     h10 = h01
#     h11 = -np.power(np.pi, 2) * np.sin(np.pi * x[:, 1])
#     return np.array([[np.sum(h00), np.sum(h01)], [np.sum(h10), np.sum(h11)]]).reshape((2, 2))
#
#
# def newton(x_init, epsilon=1e-2, m_lambda=0.01):
#     start_time = time.time()
#     xt = x_init
#     costs = []
#     while True:
#         g = grad(xt)
#         delta_x0, delta_x1 = g
#         if abs(delta_x0) < epsilon and abs(delta_x1) < epsilon:
#             break
#         h = hessian(xt)
#         xt = xt.T - m_lambda * np.linalg.inv(h).dot(g)  # 为了使其收敛，步长放小为0.01
#         xt = xt.T
#         current_cost = cost(xt)
#         costs.append(current_cost)
#     print("Newton: xt^* =\n", xt)
#     print("Newton: time usage =", time.time()-start_time)
#     print("Newton: iter_times =", len(costs))
#     return costs, xt
#
#
# def DFP(x_init, epsilon=1e-2, m_lambda=0.01):
#     start_time = time.time()
#     costs = []
#     xt = x_init
#     gt = grad(xt)
#     Bt = np.array([[1, 0], [0, 1]]).reshape((2, 2))  # B_0
#     while True:
#         delta_x0, delta_x1 = gt
#         if abs(delta_x0) < epsilon and abs(delta_x1) < epsilon:
#             break
#
#         xt1 = xt.T - m_lambda * Bt.dot(gt)  # 为了使其收敛，步长放小为0.01
#         xt1 = xt1.T  # 得到新的迭代点
#
#         gamma_t = (xt1 - xt).T
#         gt1 = grad(xt1)  # 得到新的梯度
#         yt = gt1 - gt
#
#         # 向量化计算Bt时，有点计算技巧，即先算标量，后点乘求和
#         # (gamma_t / gamma_t.T.dot(yt).T).dot(gamma_t.T)
#         Bt = Bt + (gamma_t / gamma_t.T.dot(yt).T).dot(gamma_t.T) - \
#              Bt.dot(yt).dot(yt.T).dot(Bt) / yt.T.dot(Bt).dot(yt)  # 计算新的替代矩阵
#         gt = gt1
#         xt = xt1
#
#         current_cost = cost(xt)
#         costs.append(current_cost)
#     print("DFP: xt^* =\n", xt)
#     print("DFP: time usage =", time.time() - start_time)
#     print("DFP: iter_times =", len(costs))
#     return costs, xt
#
#
# def BFGS(x_init, epsilon=1e-2, m_lambda=0.01):
#     start_time = time.time()
#     costs = []
#     xt = x_init
#     gt = grad(xt)
#     I = np.array([[1, 0], [0, 1]]).reshape((2, 2))
#     Bt = I  # B_0
#     while True:
#         assert N == 1
#         delta_x0, delta_x1 = gt
#         if abs(delta_x0) < epsilon and abs(delta_x1) < epsilon:
#             break
#
#         xt1 = xt.T - m_lambda * Bt.dot(gt)  # 为了使其收敛，步长放小为0.01
#         xt1 = xt1.T  # 得到新的迭代点
#
#         gamma_t = (xt1 - xt).T
#         gt1 = grad(xt1)  # 得到新的梯度
#         yt = gt1 - gt
#
#         Vt = I - (gamma_t / yt.T.dot(gamma_t)).dot(yt.T)
#         Bt = Vt.dot(Bt).dot(Vt.T) + (gamma_t / yt.T.dot(gamma_t)).dot(gamma_t.T)
#         gt = gt1
#         xt = xt1
#
#         current_cost = cost(xt)
#         costs.append(current_cost)
#     print("BFGS: xt^* =\n", xt)
#     print("BFGS: time usage =", time.time() - start_time)
#     print("BFGS: iter_times =", len(costs))
#     return costs, xt
#

def LBFGS(x_init, epsilon=1e-2, m_lambda=0.01, m=10, grad=None, cost=None):
    start_time = time.time()
    costs = []
    rho_gamma_list = []
    V_list = []
    xt = x_init
    gt = grad(xt)
    I = np.array([[1, 0], [0, 1]]).reshape((2, 2))
    Bt = I  # B_0
    while True:
        assert N == 1 and m > 0
        delta_x0, delta_x1 = gt
        if abs(delta_x0) < epsilon and abs(delta_x1) < epsilon:
            break

        xt1 = xt.T - m_lambda * Bt.dot(gt)  # 为了使其收敛，步长放小为0.01
        xt1 = xt1.T  # 得到新的迭代点

        gamma_t = (xt1 - xt).T
        gt1 = grad(xt1)  # 得到新的梯度
        yt = gt1 - gt

        rho_t = 1 / yt.T.dot(gamma_t)
        Vt = I - rho_t * yt.dot(gamma_t.T)
        rho_gamma_t = rho_t * gamma_t.dot(gamma_t.T)

        rho_gamma_list.append(rho_gamma_t)
        V_list.append(Vt)

        V = Vt
        temp_v_list = [V]
        for v in V_list[-2:-min(m-1, len(V_list))-1:-1]:
            V = V.dot(v)
            temp_v_list.append(V)
        Bt_part1 = V.T.dot(I).dot(V)

        Bt_part2 = rho_gamma_t
        for rg in rho_gamma_list[-2:-min(m-1, len(rho_gamma_list))-1:-1]:
            Bt_part2 += V.T.dot(rg).dot(V)

        Bt = Bt_part1 + Bt_part2

        gt = gt1
        xt = xt1

        current_cost = cost(xt)
        costs.append(current_cost)
    print("LBFGS: xt^* =\n", xt)
    print("LBFGS: time usage =", time.time() - start_time)
    print("LBFGS: iter_times =", len(costs))
    return costs, xt

def my_fun(xy):
    X = xy[0]
    Y = xy[1]
    return (X - 1) ** 4 + Y ** 2


def my_gfun(xy):
    X = xy[0]
    Y = xy[1]
    return np.array([4 * (X - 1) ** 3, 2 * Y])

if __name__ == '__main__':



    # np.random.seed(121)
    # X0 = np.random.randn(N, 2)
    X0 = [-20, 100]

    # costs_newton, xt_newton = newton(x_init=X0)
    # costs_dfp, xt_dfp = DFP(x_init=X0)
    # costs_bfgs, xt_bfgs = BFGS(x_init=X0)
    costs_lbfgs, xt_lbfgs = LBFGS(x_init=np.array(X0),cost=my_fun,grad=my_gfun)

    ax = plt.figure().add_subplot(1, 1, 1)
    # ax.plot(range(len(costs_newton)), costs_newton, '-',
    #         range(len(costs_dfp)), costs_dfp, '--',
    #         range(len(costs_bfgs)), costs_bfgs, '-.',
    #         range(len(costs_lbfgs)), costs_lbfgs, '-x')
    ax.plot(range(len(costs_lbfgs)), costs_lbfgs, '-x')
    # ax.set(xlabel='times', ylabel='cost', title="Newton's method(N=%d)" % N)
    plt.legend(('LBFGS'), loc='upper right')
    plt.show()



