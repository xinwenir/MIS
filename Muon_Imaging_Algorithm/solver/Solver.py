# author:高金磊
# datetime:2022/3/21 16:51
import time
from copy import copy
from datetime import datetime

import numpy as np
import scipy.sparse
import optimparallel
from scipy.optimize import differential_evolution, least_squares, OptimizeResult, dual_annealing
from scipy.optimize.linesearch import line_search

from DataTools.Sparse_Matrix_Tools import Find_similar_is
from SysTools.MyTools import myPrint
from solver import Objective_function
from solver.Jacobi import get_jac_equation_for_least_squares, Jac_line_search
from solver.Objective_function import constr_f
from DataFactory import Setting,DataManage


class Solver:
    """
    支持多线程,一次数据存取后通过向run方法中传入新参数(参数,求解器类型等)
    """

    class STP:
        def __init__(self):
            self._stp = 0

        def commands(self, new_stp):
            setting = Setting.Setting.get_instance()
            # 无优先级中断,旧中断将要屏蔽新中断
            if new_stp not in setting.STP_dic.keys():
                setting.get_loger().write("错误的STP指令!!!")
                return False
            if self._stp == new_stp:
                pass
            else:
                if self._stp == 0:
                    self._stp = new_stp
                    setting.get_loger().write(setting.STP_dic[new_stp], printer=myPrint.myPrint_Err)
                else:
                    setting.get_loger().write("指令失败,上一个中断请求没有撤销")

        def stp_Exception(self, env):
            setting = Setting.Setting.get_instance()
            if type(env) is not tuple:
                env = (env,)
                setting.get_loger().write("中断环境应该传入元组!!!!!!", printer=myPrint.myPrint_Hint)
            self._env = env
            return Exception(setting.STP_dic[self._stp])

        def is_STP(self):
            self.timing()
            return self._stp != 0

        def get_env(self):
            return self._env

        def cancel(self):
            self._stp = 0

        _start_time = None

        def timing(self):
            setting = Setting.Setting.get_instance()
            if self._start_time is None:
                setting.get_loger().write("开始计时,限制运行:%ds" % (setting.get_max_time()))
                self._start_time = time.time()
            if time.time() - self._start_time > setting.get_max_time():
                self.commands(2)

    def __init__(self):
        """
        根据设置获取不可变的数据
        """

        global A, d, d_err, bounds, refs, oldi_newi, oldj_newj, args, shape, log, newj_oldj, newi_oldi, x0, m, n,dataManager
        dataManager = DataManage.DataManager.get_instance()
        setting = Setting.Setting.get_instance()
        # beta = setting.get_beta()
        alpha_s, alpha_x, alpha_y, alpha_z = setting.get_alpha()
        args = (alpha_s, alpha_x, alpha_y, alpha_z)  # 今后改成字典的形式
        shape = dataManager.mesh.get_shape()
        log = setting.get_loger()
        log.write("读取并处理数据%s" % (datetime.now()), printer=myPrint.myPrint_Hint)
        A, d, d_err, bounds, refs, oldi_newi, oldj_newj = dataManager.Make_A()
        dataManager.record_data_detail()
        newj_oldj = {}
        for oldj in oldj_newj.keys():
            newj_oldj[oldj_newj[oldj]] = oldj
        newi_oldi = {}
        for oldi in oldi_newi.keys():
            newi_oldi[oldi_newi[oldi]] = oldi
        m, n = A.shape

        x0 = copy(refs)

        # 收紧约束   使得退火算法等算法变成局部算法
        # scope = 0.25
        # for i in range(n):
        #     bounds[i][0]=max(bounds[i][0],middle[i]-scope)
        #     bounds[i][1]=min(bounds[i][1],middle[i]+scope)

        # 初始值选择问题
        # for i in range(n):
        # x0[i]=random.uniform(bounds[i][0],bounds[i][1])
        # x0[i]=(bounds[i][1])
        # x0[i] = (bounds[i][0])
        # pass

    @staticmethod
    def show_env():
        """
        显示当前的基本配置
        :return:
        """
        dataManager.record_data_detail()

    all_time=0
    def run(self, beta=None, def_x=None):
        """
        执行反演
        :return: 本次反演的结果;下次运行的建议(当前和下次运行的beta值,是否需要继续运行等);当前的目标函数
        """
        setting = Setting.Setting.get_instance()
        if beta is not None:
            setting.set_beta(beta)
        stp = self.STP()
        global result, end_time, start_time
        log.write("开始运算.....%s\n" % (datetime.now()))

        # SHGO算法
        # result=scipy.optimize.shgo(constr_f,bounds=bounds,iters=1,options={"maxfev":2,"maxiter":3,"maxev":4,"maxtime":100},
        #                            args=(G, d_values, shape, d_err, newj_oldj, refs))
        method = setting.get_solve_method()
        maxiter = setting.get_max_iter()
        disp = setting.get_show_detailed()
        start_time = time.time()
        result = None

        if def_x is not None:
            x0 = np.array(Solver.in_bounds(def_x,bounds=bounds))
        else:
            x0 = np.array(Solver.in_bounds(refs,bounds=bounds))

        try:
            if method == 'L-BFGS-B':
                result = scipy.optimize.minimize(constr_f,
                                                 method='L-BFGS-B',
                                                 bounds=bounds,
                                                 x0=x0, jac=get_jac_equation_for_least_squares,
                                                 options={"maxfun": maxiter, 'disp': disp},
                                                 args=(
                                                     A, d, shape, d_err, newj_oldj, oldj_newj, np.array(refs), bounds,
                                                     stp,
                                                     args))
            elif method=="L-BFGS-B_optimparallel":
                result = optimparallel.minimize_parallel(constr_f,
                                                 bounds=bounds,
                                                 x0=x0, jac=get_jac_equation_for_least_squares,
                                                 options={"maxfun": maxiter, 'disp': disp},
                                                 args=(
                                                     A, d, shape, d_err, newj_oldj, oldj_newj, np.array(refs), bounds,
                                                     stp,
                                                     args))
            elif method == 'CG':
                myPrint.myPrint_Wran("CG", "不支持增加约束!!!!!罚函数将被添加在jac中")
                result = scipy.optimize.minimize(constr_f,
                                                 method='CG',
                                                 x0=x0, jac=get_jac_equation_for_least_squares,
                                                 options={"maxiter": 10, 'disp': disp},
                                                 args=(
                                                     A, d, shape, d_err, newj_oldj, oldj_newj, np.array(refs), bounds,
                                                     stp,
                                                     args))
            elif method == 'TNC':
                result = scipy.optimize.minimize(constr_f,
                                                 method='TNC',
                                                 bounds=bounds,
                                                 x0=x0, jac=get_jac_equation_for_least_squares,
                                                 options={'disp': disp, 'maxfun': maxiter},
                                                 args=(
                                                     A, d, shape, d_err, newj_oldj, oldj_newj, np.array(refs), bounds,
                                                     stp,
                                                     args))
            elif method == 'least_squares':
                # 非线性最小二乘
                result = least_squares(constr_f,
                                       bounds=[[i[0] for i in bounds], [i[1] for i in bounds]],
                                       jac=get_jac_equation_for_least_squares, verbose=2, x0=x0, max_nfev=maxiter,
                                       method='trf',
                                       # method='dogbox',#不可行
                                       tr_solver='lsmr',
                                       loss="soft_l1",
                                       args=(
                                       A, d, shape, d_err, newj_oldj, oldj_newj, np.array(refs), bounds, stp, args))
            elif method == 'SLSQP':
                result = scipy.optimize.minimize(constr_f,
                                                 method='SLSQP',
                                                 hess='2-point',
                                                 bounds=bounds,
                                                 x0=x0, jac=get_jac_equation_for_least_squares,
                                                 options={"maxiter": maxiter, 'disp': disp},
                                                 args=(
                                                     A, d, shape, d_err, newj_oldj, oldj_newj, np.array(refs), bounds,
                                                     stp,
                                                     args))
            elif method == 'trust-constr':
                result = scipy.optimize.minimize(constr_f,
                                                 method='trust-constr',
                                                 bounds=bounds,
                                                 hess='2-point',
                                                 x0=x0, jac=get_jac_equation_for_least_squares,
                                                 options={"maxiter": maxiter, 'disp': disp},
                                                 args=(
                                                     A, d, shape, d_err, newj_oldj, oldj_newj, np.array(refs), bounds,
                                                     stp, args))
            elif method == 'Powell':
                result = scipy.optimize.minimize(constr_f,
                                                 method='Powell',
                                                 bounds=bounds,
                                                 hess='2-point',
                                                 x0=x0, jac=get_jac_equation_for_least_squares,
                                                 options={"maxiter": maxiter, 'disp': disp},
                                                 args=(
                                                     A, d, shape, d_err, newj_oldj, oldj_newj, np.array(refs), bounds,
                                                     stp,
                                                     args))
            elif method=="dual_annealing":
                # 使用双重退火求函数的全局最小值。
                result = dual_annealing(constr_f, bounds=bounds,
                                        x0=x0, maxiter=300,maxfun=100,initial_temp=300,visit=1.5,
                                        args=(
                                            A, d, shape, d_err, newj_oldj, oldj_newj, np.array(refs), bounds,
                                            stp,
                                            args))

            elif method=="line_search":
                search_gradient=np.array([-1 for i in x0])
                jac_fun=Jac_line_search(get_jac_equation_for_least_squares,(A, d, shape, d_err, newj_oldj, oldj_newj, np.array(refs), bounds,
                                            stp,
                                            args)).Jac_line_search
                xxx=line_search(constr_f,jac_fun,x0,search_gradient,args=(
                                            A, d, shape, d_err, newj_oldj, oldj_newj, np.array(refs), bounds,stp,args
                ))
                result = OptimizeResult()
                result.x=xxx[-1]
                result.success=True
                end_time = time.time()
                self.all_time += end_time - start_time
                misfit = Objective_function.o1_last
                dataManager.collect_beta_misfit_norms(dataManager.choose_beta().get_beta()[0], misfit,
                                                      Objective_function.o3_last)
                # 根据策略更新beta
                beta_info = dataManager.choose_beta().get_beta(A * np.array([result.x]).T, d)
                self.reset()
                return result.x, beta_info, misfit
            elif method=="my_lbfgs":
                from solver.nonlinear_solver.test import lbfgs as my_lbfgs

                res=my_lbfgs(fun=constr_f, gfun=get_jac_equation_for_least_squares, x0=x0, args=(
                    A, d, shape, d_err, newj_oldj, oldj_newj, np.array(refs), bounds,
                    stp,
                    args)
                         )
                result = OptimizeResult()
                result.x=res[0]
                result.success=False
        except Exception as e:
            setting.get_loger().write(str(e))
            if result is None:
                result = OptimizeResult()
            result.x = stp.get_env()[0]
            result.success = False

        end_time = time.time()
        self.all_time+=end_time-start_time
        misfit = Objective_function.o1_last
        dataManager.collect_beta_misfit_norms(dataManager.choose_beta().get_beta()[0], misfit,
                                              Objective_function.o3_last)
        # 根据策略更新beta
        beta_info = dataManager.choose_beta().get_beta(A * np.array([result.x]).T, d)
        self.reset()
        return result.x, beta_info, misfit

    def show_state(self):
        if result is None:
            log.write("求解指令还未执行或者执行未完成", printer=myPrint.myPrint_Err)
            return False

        log.write("本次求解花费时间:%ds,累计用时%d s" % (end_time - start_time,self.all_time))
        log.write("是否成功退出:" + str(result.success))

    def save_res(self):
        if result is None:
            log.write("求解指令还未执行或者执行未完成", printer=myPrint.myPrint_Err)
            return False
        dataManager.Postprocessor(result.x)

    def reset(self):
        global start_time, end_time
        Objective_function.reset()
        # start_time=time.time()

    def close(self):
        log.write("程序运行完毕...%s\n" % (datetime.now()))
        log.flush()
    @staticmethod
    def in_bounds(x,bounds):
        setting = Setting.Setting.get_instance()
        x0=copy(x)
        for i in range(len(x0)):
            if x0[i]<bounds[i][0]:
                setting.get_loger().err("%d 对应的初始值%f不在约束区间内%s"%(i,x0[i],str(bounds[i])))
                x0[i]=bounds[i][0]
            if x0[i]>bounds[i][1]:
                setting.get_loger().err("%d 对应的初始值%f不在约束区间内%s"%(i,x0[i],str(bounds[i])))
                x0[i] = bounds[i][1]
        return x0

# result = least_squares(constr_f, jac_sparsity=jac_sparsity, method='trf', tr_solver='lsmr',
#                        bounds=[[i[0] for i in bounds], [i[1] for i in bounds]], x0=x0,
#                        args=(G, d_values, m, n, d_err,))
# 求解线性最小二乘问题,不能加目标函数
# result=scipy.optimize.lsq_linear(G,d_values,bounds=[[i[0] for i in bounds], [i[1] for i in bounds]],method='trf',verbose=2,max_iter=100)
# def constr_f(x, A, b, shape, d_err, j_oldj, m0: np.array) -> np.array:


# 线性最小二乘
# result=lsq_linear( A=G,b=d_values,
#                    bounds=[[i[0] for i in bounds], [i[1] for i in bounds]],
#                        verbose=2 ,method='trf',max_iter=10,
#                        lsq_solver='lsmr')
# 稀疏矩阵最小二乘
# result = linalg.lsmr(A=G, b=d_values, show=True, maxiter=20)

# 结果和ref文件差的绝对值
# make_refs_ps_err(r'M:\pycharm\Inversion\DataTools\data\04_58MaMian_ref.den',res,r"M:\pycharm\Inversion\DataTools\data\refs_ps")
# filter_res(res)#因为已经没有无效数据，此处无需再删除无效数据
# 重新排列结果,并对结果进行处理(标记空气)
# format_res(res, unneed_j, air_j, -0.1234)

# tool.mode1(min=0.01, ignore_value=[0, -0.1234], max=4.0001)
# tool.mode2(ignore_value=[0, -0.1234])
