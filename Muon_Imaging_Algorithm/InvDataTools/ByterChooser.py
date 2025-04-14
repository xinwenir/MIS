# author:高金磊
# datetime:2022/4/9 9:33
import random

from InvDataFactory import Setting


class ByterChooser:
    """
    beta选择器

    """
    def __init__(self, dobs, derr, beta, beta_coefficient):
        self.dobs = dobs
        self.derr = derr
        self.beta = beta
        self.beta_coefficient = beta_coefficient
        self.beta_min = 0  # 设置beta的下限
        self.beta_max = self.beta * 2

    _count = 0

    def _calculate(self, dpred):
        res = 0
        for i in range(len(dpred)):
            res += pow((dpred[i]) / self.derr[i], 2)
        return res

    def get_beta(self, pred=None, obsd=None):
        """
        请求beta,根据设置来采取beta的选取

        :param pred: 预测的d
        :param obsd: 观测的d
        :return: 经过选择得到的beta
        """
        if pred is None:
            return self.beta, False
        # self._count += 1
        # if self._count <= 3 or self._count % 3 != 0:
        #     return self.beta
        # if self._count < 1:
        #     return self.beta
        # 执行beta自动化选择
        setting = Setting.Setting.get_instance()
        mode = setting.sys_setting.beta_choose_mode
        if mode == 1:
            return setting.get_beta(), False
        elif mode == 2:
            return self._get_beta_by_X2(pred, obsd)
        elif mode == 3:
            return self._get_beta_by_L(pred, obsd)

        return self.beta, False

    _get_beta_by_X2_count=0
    def _get_beta_by_X2(self, pred, obsd):
        """
        使用基于卡方的beta值的寻找

        :param pred: 预测的d
        :param obsd: 观测的d
        :return: beta及本次beta是否作为最终的beta
        """
        self._get_beta_by_X2_count+=1
        if self._get_beta_by_X2_count >= 11:
            print("已经达到最大迭代深度")
            return self.beta, False
        setting = Setting.Setting.get_instance()
        beta_accuracy = setting.sys_setting.beta_accuracy
        N = len(self.derr)
        stat = 0
        for i in range(len(self.derr)):
            stat += pow((pred[i] - obsd[i]) / self.derr[i], 2)
        search_mode = setting.sys_setting.search_mode
        # stat=0
        # data=[]
        # data.append(pred.T.tolist()[0])
        # data.append(obsd.T.tolist()[0])
        # for i in range(len(data[0])):
        #     stat+=pow(data[0][i]-data[1][i],2)/data[1][i]
        # stat,p,dof,expected = chi2_contingency(data)
        coefficient = stat / N
        old_beta = self.beta
        if search_mode == 1:
            # 线性搜索
            if abs(self.beta_coefficient - coefficient) / self.beta_coefficient < 0.005:
                return self.beta, False
            if self.beta_coefficient < coefficient:
                # self.beta-=0.1*self.beta
                self.beta /= random.uniform(1, 2)
            else:
                # self.beta += 0.1*self.beta
                self.beta *= random.uniform(1, 1.5)
            if abs(self.beta - old_beta) < beta_accuracy:
                print("beta搜索满足精度")
                return self.beta, False
            setting.get_loger().write("old_beta:%f new_beta:%f target_coefficient:%f coefficient:%f" % (
            old_beta, self.beta, self.beta_coefficient, coefficient))
            setting.get_loger().write("misfit:%f" % (stat))
            return self.beta, True
        # 二分搜索
        elif search_mode == 2:
            ##beta搜索范围过小
            # if self.beta_max - self.beta_min < beta_accuracy * 2:
            #     print("beta搜索满足精度")
            #     return self.beta, False
            ##不再下降
            if abs(self.beta_coefficient - coefficient) / self.beta_coefficient < 0.01:
                return self.beta, False
            if abs(self.beta_coefficient-coefficient)< 0.02:
                return self.beta, False
            if self.beta_coefficient < coefficient:
                # beta太大
                self.beta_max = self.beta
                self.beta_min-=self.beta/10
                self.beta += self.beta_min
                self.beta /= 2
                setting.get_loger().write("old_beta:%f new_beta:%f 搜索范围:(%f , %f)" % (
                    old_beta, self.beta, self.beta_min, self.beta_max))
                setting.get_loger().write("target_coefficient:%f coefficient:%f" % (self.beta_coefficient, coefficient))
                return self.beta, True
            else:
                self.beta_min = self.beta
                self.beta_max+=self.beta/10 ##略微放大区间
                self.beta += self.beta_max
                self.beta /= 2
                setting.get_loger().write("old_beta:%f new_beta:%f 搜索范围:(%f , %f)" % (
                    old_beta, self.beta, self.beta_min, self.beta_max))
                setting.get_loger().write("target_coefficient:%f coefficient:%f" % (self.beta_coefficient, coefficient))
                return self.beta, True
        # 二分搜索增加不确定因子
        elif search_mode == 3:
            if self.beta_max - self.beta_min < beta_accuracy * 2:
                print("beta搜索满足精度")
                return self.beta, False
            if self.beta_coefficient < coefficient:
                # beta太大
                self.beta_max = self.beta
                self.beta = random.uniform(self.beta_min, self.beta_max)
                setting.get_loger().write("old_beta:%f new_beta:%f 搜索范围:(%f , %f)" % (
                    old_beta, self.beta, self.beta_min, self.beta_max))
                setting.get_loger().write("target_coefficient:%f coefficient:%f" % (self.beta_coefficient, coefficient))
                return self.beta, True
            else:
                self.beta_min = self.beta
                self.beta = random.uniform(self.beta_min, self.beta_max)
                setting.get_loger().write("old_beta:%f new_beta:%f 搜索范围:(%f , %f)" % (
                    old_beta, self.beta, self.beta_min, self.beta_max))
                setting.get_loger().write("target_coefficient:%f coefficient:%f" % (self.beta_coefficient, coefficient))
                return self.beta, True
        elif search_mode == 4:
            # ubc
            target_misfit = self.beta_coefficient * N
            misfit = stat
            print("abs(misfit-target_misfit)/target_misfit:%f" % (abs(misfit - target_misfit) / target_misfit))
            setting.get_loger().write("target_coefficient:%f coefficient:%f" % (self.beta_coefficient, coefficient))
            # if misfit <= target_misfit or abs(misfit - target_misfit) / target_misfit < 0.05:
            if (misfit - target_misfit) / target_misfit < 0.02:
                return self.beta, False
            # 接近的时候

            if abs(misfit - target_misfit) / target_misfit < 0.2:
                self.beta = self.beta * 0.5  # ???????????0.5
                setting.get_loger().write("old_beta:%f new_beta:%f " % (
                    old_beta, self.beta))

                return self.beta, True
            else:
                self.beta = self.beta * 0.25  # ??????0.25
                setting.get_loger().write("old_beta:%f new_beta:%f " % (
                    old_beta, self.beta))
                return self.beta, True

        return self.beta, False

    def _get_beta_by_L(self, pred, obsd):
        pass
