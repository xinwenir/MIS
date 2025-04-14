# author:高金磊
# datetime:2022/4/6 14:56
import math

import numpy as np

from InvSysTools.MyTools import myPrint

from matplotlib import pyplot as plt


class Show_pred_obsd_derr:
    """
    显示/存储反演结果的一些参数
    """

    def __init__(self, pred_obsd_derr_file=None, data=None):
        if data is not None:
            if len(data[0]) != 3:
                raise Exception("请输入n*3的数组,当前维度不满足要求")
            self.data = data
            return
        if pred_obsd_derr_file is None:
            raise Exception("请至少传入一个参数")
        file = open(pred_obsd_derr_file, 'r')
        lines = file.readlines()
        file.close()
        self.data = []
        for line in lines:
            middle = line.split()
            self.data.append([float(middle[0]), float(middle[1]), float(middle[2])])

    def store_res(self, file):
        """
        存储结果到file中
        :param file:文件
        """
        file_obj = open(file, 'w')
        for datum in self.data:
            file_obj.write(str(datum[0]))
            file_obj.write(" ")
            file_obj.write(str(datum[1]))
            file_obj.write(" ")
            file_obj.write(str(datum[2]))
            file_obj.write('\n')
        file_obj.close()

    def show_pred_obsd(self):
        """
        显示Ax和d的分布
        """
        pred = []
        obsd = []
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot()
        for i in range(len(self.data)):
            datum = self.data[i]
            pred.append(datum[0])
            obsd.append(datum[1])
        pred_1 = []
        pred_2 = []
        pred_3 = []
        obsd_1 = []
        obsd_2 = []
        obsd_3 = []
        for i in range(len(pred)):
            if pred[i] < obsd[i]:
                pred_1.append(pred[i])
                obsd_1.append(obsd[i])
            elif pred[i] >= obsd[i]:
                pred_2.append(pred[i])
                obsd_2.append(obsd[i])
            else:
                pred_2.append(pred[i])
                obsd_2.append(obsd[i])

        plt.scatter(pred_1, obsd_1, marker='.', c='red')
        plt.scatter(pred_2, obsd_2, marker='.', c='green')
        plt.scatter(pred_3, obsd_3, marker='o', c='yellow')
        plt.plot([min(pred) - 1, max(pred) + 1], [min(pred) - 1, max(pred) + 1], c='black', ls='-')
        ax.set_xlabel("pred")
        ax.set_ylabel("obsd")
        plt.show()

    def show_pred_obsd_diff(self, interval=0.04, max_ignore=5):
        """
        显示Ax-d的结果分布 使用interval进行离散
        :param interval: 离散的范围如{a,a+interval)被视为a
        :param max_ignore:为了便于显示,去掉数据最多的几个值
        :return:
        """
        diff = []
        for i in range(len(self.data)):
            datum = self.data[i]
            diff.append(datum[0] - datum[1])
        max_value = max(diff)
        min_value = min(diff)
        Statistics = [0] * int((max_value - min_value) / interval + 1)
        for i in range(len(diff)):
            Statistics[int((diff[i] - min_value) / interval)] += 1
        x = [i * interval + min_value for i in range(len(Statistics))]

        Statistics_new = []
        x_new = []
        for i in range(len(Statistics)):
            if Statistics[i] >= max_ignore:
                Statistics_new.append(Statistics[i])
                x_new.append(x[i])
        fig, ax = plt.subplots(dpi=400, figsize=(12, 8))

        plt.plot(x_new, Statistics_new)
        plt.plot([0, 0], [0, max(Statistics) + 1], c='red', ls='dotted')
        # plt.bar(x_new,Statistics_new)
        # plt.hist(Statistics_new,x_new)
        # plt.bar(x,Statistics)
        # for a, b in zip(x, Statistics):
        # if b==0:
        #     continue
        # plt.text(a, b + 1, b, ha='center', va='bottom')

        plt.show()

    def show_pred_obs_obsderr(self):
        """
        """
        import warnings
        warnings.warn("数据太多,不适合这样处理", DeprecationWarning)
        pred = []
        obsd = []
        obsd_up = []
        obsd_down = []
        for i in range(len(self.data)):
            datum = self.data[i]
            if abs(datum[0] - datum[1]) / datum[2] < 10:
                continue

            pred.append(datum[0])
            obsd.append(datum[1])
            obsd_up.append(datum[1] + datum[2])
            # if obsd_up[-1]>60:
            #     print(obsd_up[-1])
            obsd_down.append(datum[1] - datum[2])
        x = [i for i in range(len(obsd))]
        # plt.scatter(x,obsd,marker='.')
        plt.scatter(x, pred, marker='.')
        plt.fill_between(x, obsd_up, obsd_down,  # 上限，下限
                         facecolor='green',  # 填充颜色
                         edgecolor='red',  # 边界颜色
                         alpha=0.1)  # 透明度
        plt.show()


class Norm_tools:
    """
    平滑度,misfit之间的关系处理
    """
    def __init__(self, all_misfit_ms_smooth_file=None, all_misfit_ms_smooth_list=None):
        """
        通过文件初始化Norm_tools的一个对象
        :param all_misfit_ms_smooth_file:
        :param all_misfit_ms_smooth_list:
        """
        if all_misfit_ms_smooth_list is not None:
            if len(all_misfit_ms_smooth_list[0]) != 4:
                raise Exception("请输入n*4的数组,当前维度不满足要求")
            self.data = all_misfit_ms_smooth_list
            return
        if all_misfit_ms_smooth_file is None:
            raise Exception("请至少传入一个参数")
        file = open(all_misfit_ms_smooth_file, 'r')
        lines = file.readlines()
        file.close()
        self.data = []
        for line in lines:
            middle = line.split()
            self.data.append([float(middle[0]), float(middle[1]), float(middle[2]), float(middle[3])])

    def store_res(self, file):
        """
        保存结果到文件中
        :param file:
        """
        file_obj = open(file, 'w')
        for datum in self.data:
            file_obj.write(str(datum[0]))
            file_obj.write(" ")
            file_obj.write(str(datum[1]))
            file_obj.write(" ")
            file_obj.write(str(datum[2]))
            file_obj.write(" ")
            file_obj.write(str(datum[3]))
            file_obj.write('\n')
        file_obj.close()

    def show_norm(self):
        """
        显示norm和迭代次数的关系
        """
        data = self.data
        x = [i for i in range(len(data))]
        all = [i[0] for i in data]
        misfit = [i[1] for i in data]
        ms = [i[2] for i in data]
        smooth = [i[3] for i in data]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        my_x_ticks = np.arange(1, len(misfit), 1)
        plt.xticks(my_x_ticks)
        ax.plot(x, all, '-', c='red', label='all')
        # ax.plot(x, misfit, '-', c='yellow', label='misfit')
        ax2 = ax.twinx()
        # ax2.plot(x, ms, c='green', label='ms')
        ax2.plot(x, smooth, c='black', label='smooth')
        ax.legend(loc='upper center')
        ax2.legend()
        ax.set_ylabel(r"misfit and all")
        ax2.set_ylabel(r"ms and smooth")
        ax.set_xlabel(r"The number of iterations")
        # ax.show()
        plt.show()


class Beta_log:
    """
    Beta和其他数据的关系(misfit,norm)
    """
    def __init__(self, file_path):
        self.file = open(file_path, 'a')
        self.file_path = file_path

    def record(self, beta, misfit, norm):
        """
        记录数据
        :param beta: 要记录的beta
        :param misfit: 需要记录的misfit
        :param norm: 需要记录的norm
        :return:
        """
        try:
            self.file.write("%f %f %f" % (beta, misfit, norm))
            self.file.write('\n')
        except Exception as e:
            myPrint.myPrint_Err(e)
            return False
        return True

    def show_L(self):
        """
        显示L曲线
        """
        self.file.flush()
        plt.rcParams['font.sans-serif'] = ['SimHei']
        plt.rcParams['axes.unicode_minus'] = False
        # matplotlib画图中中文显示会有问题，需要这两行设置默认字体
        file = open(self.file_path, 'r')
        lines = file.readlines()
        betas = []
        misfits = []
        norms = []
        for line in lines:
            middle = line.split()
            betas.append(float(middle[0]))
            misfits.append(float(middle[1]))
            # misfits.append(math.log2(float(middle[1])))
            norms.append(float(middle[2]))
            # norms.append(math.log2(float(middle[2])))
        fig = plt.figure(figsize=(12, 8))
        # plt.xticks(norms)
        # plt.yticks(misfits)
        plt.scatter(norms, misfits)
        plt.xlabel("norm")
        plt.ylabel("misfit")
        # plt.plot(norms,misfits,'-', c='red', label='L')
        plt.show()

    def show_similarL(self):
        """
        显示L曲线的相关信息
        :return:
        """
        self.file.flush()
        plt.rcParams['font.sans-serif'] = ['SimHei']
        plt.rcParams['axes.unicode_minus'] = False
        # matplotlib画图中中文显示会有问题，需要这两行设置默认字体
        file = open(self.file_path, 'r')
        lines = file.readlines()
        middle = []
        for line in lines:
            mm = line.split()
            middle.append([float(mm[0]), float(mm[1]), float(mm[2])])

        def sortkey(num):
            return float(num[0])

        middle.sort(key=sortkey)
        middle.reverse()
        data = []
        # 计算（(dnorm/dmisfit)/dbeta）,这个值越大表示变坏的越快
        for i in range(1, len(middle)):
            if middle[i][0] == middle[i - 1][0] or middle[i][1] == middle[i - 1][1] or middle[i][2] == middle[i - 1][2]:
                myPrint.myPrint_Wran(middle[i], "被舍弃")
                continue
            dbeta = middle[i][0] - middle[i - 1][0]
            dmisfit = middle[i][1] - middle[i - 1][1]
            dnorm = middle[i][2] - middle[i - 1][2]
            data.append(dmisfit / dnorm)
            print(data[-1], middle[i])
            # data.append((float(line[2])-float(pre_line[2]))/(float(line[1])-float(pre_line[1]))/(float(line[0])-float(pre_line[0])))

        # plt.scatter(norms, misfits)
        x = [i for i in range(len(data))]
        plt.xlabel("iter number")
        plt.ylabel("dmisfit/dnorm")
        plt.plot(x, data, '-', c='black', label='L')
        second_derivative = []
        for i in range(1, len(data)):
            second_derivative.append(data[i] - data[i - 1])
        x2 = [i + 1 for i in range(len(second_derivative))]
        plt.plot(x2, second_derivative, '-', c='green', label='second_derivative')
        plt.plot(x, [0 for i in range(len(x))], '.', c='red')
        plt.show()

    def show_Beta_search(self):
        """
        显示迭代次数和beta.norm,misfit的关系
        """
        self.file.flush()
        file = open(self.file_path, 'r')
        lines = file.readlines()
        betas = []
        misfits = []
        norms = []
        for line in lines:
            middle = line.split()
            betas.append(float(middle[0]))
            misfits.append(float(middle[1]))
            norms.append(float(middle[2]))
        fig = plt.figure(figsize=(12, 8))
        xs = [i for i in range(1, len(betas) + 1)]
        ax = fig.add_subplot(111)
        plt.xticks(xs)
        # ax.scatter(xs, misfits, marker='.', c='red', label='all')
        ax.plot(xs, misfits, '-', c='red', label='misfit')
        ax2 = ax.twinx()
        # ax2.scatter(xs, betas, marker='.', c='green', label='ms')
        ax2.plot(xs, betas, '-', c='green', label='beta')
        ax.legend(loc='upper center')
        ax2.legend()
        ax.set_ylabel(r"misfit")
        ax2.set_ylabel(r"beta")
        ax.set_xlabel(r"iter")
        plt.show()

    def show(self):
        """
        显示beta和norm,misfit的关系
        """
        self.file.flush()
        file = open(self.file_path, 'r')
        lines = file.readlines()
        betas = []
        misfits = []
        norms = []
        for line in lines:
            middle = line.split()
            betas.append(float(middle[0]))
            misfits.append(float(middle[1]))
            norms.append(float(middle[2]))
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111)
        plt.xticks(betas)
        ax.scatter(betas, misfits, marker='.', c='red', label='all')
        ax2 = ax.twinx()
        ax2.scatter(betas, norms, marker='.', c='green', label='ms')
        ax.legend(loc='upper center')
        ax2.legend()
        ax.set_ylabel(r"misfit")
        ax2.set_ylabel(r"norm")
        ax.set_xlabel(r"beta")
        plt.show()


class Show_data_distribution:
    """
    x和对应的数量   (x与refx的差)和对应的数量

    """

    def __init__(self):
        pass

    # def _show_res_file_num(self, xfile, split=None):
    #     obsx = []
    #     with open(xfile, 'r') as f:
    #         while (fvalue := f.readline()) != "":
    #             obsx.append(float(fvalue))
    #     self.show_obsx_num(obsx, split)

    # def _show_x_refx_diff_file_num(self, xfile, refxfile, split=None):
    #     diffx = []
    #     with open(xfile, 'r') as x, open(refxfile, 'r') as refx:
    #         while (xfile := x.readline()) != "":
    #             diffx.append(float(xfile) - float(refx.readline()))
    #     self.show_x_refx_diff_num(diffx, split)

    def show_obsx_num(self, obsx, split=None):
        """
        计算值和所对应的数量

        :param obsx:
        :param split:
        """
        self.show_factory_similarx_num(obsx, "x", "count", split)

    def show_x_refx_diff_num(self, diffx, split=None):
        """
        [计算值与参考值的差]和所对应的数量

        :param diffx:
        :param split:
        """
        self.show_factory_similarx_num(diffx, "x-refs", "count", split)

    # 最大值归到最后一个格子内
    def show_factory_similarx_num(self, sx, xlabel, ylabel, split=None):
        """
        展示图形

        :param sx: 数据来源
        :param xlabel: 横坐标标签
        :param ylabel: 纵坐标标签
        :param split: 将sx分割的区间数量
        """
        if split is None:
            split = 100
        xmin = min(sx)
        xmax = max(sx)
        # 每一份的大小
        size = (xmax - xmin) / split
        pos = [0.0] * split
        sum = [0] * split
        for i in range(split):
            pos[i] = xmin + i * size + size / 2
        for i in sx:
            if i == xmax:
                sum[split - 1] += 1
            else:
                sum[math.floor((i - xmin) / size)] += 1
        # plt.plot(pos, [np.log10(i) if i>0 else 0 for i in sum])
        plt.figure(figsize=(12, 8))
        plt.plot(pos, sum)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid()
        plt.show()

class show_2D_the_density_of_res:
    """
    !借助于其他项目包:InvSolver\Seed_algorithm
    """
    def __init__(self,res_file,mesh_file) -> None:
        from InvSolver.Seed_algorithm.Tools import data_tools
        data_tool=data_tools(res_file,mesh_file)
        self.res_file=res_file
        self.data3D=data_tool.get_data_3D()
        
        # tool=Show_data_distribution()
        # tool._show_res_file_num(self.res_file)
        self.statistical_interval=[2.75,3]
        self.mode=1
    def get_z_projection(self):
        z_projection=[]
        for x in range(len(self.data3D)):
            middle=[]
            for y in range(len(self.data3D[0])):
                sum_value=0
                count=0
                min_value=10000
                max_value=-10000
                for z in range(len(self.data3D[0][0])):
                    if self.mode==1:
                        if self.data3D[x][y][z]<2.67 or self.data3D[x][y][z]>100:
                            continue
                        # 平均值策略
                        sum_value+=self.data3D[x][y][z]-2.65
                        count+=1
                    elif self.mode==2:
                        if self.data3D[x][y][z]<1 or self.data3D[x][y][z]>100:
                            continue
                        #最小值策略
                        min_value=min(self.data3D[x][y][z],min_value)
                        count=1
                    elif self.mode==3:
                        if self.data3D[x][y][z]<2.7 or self.data3D[x][y][z]>100:
                            continue
                        #最大值策略
                        max_value=max(self.data3D[x][y][z],max_value)
                        count=1
                    elif self.mode==4:
                        #计数
                        if self.data3D[x][y][z]<self.statistical_interval[1] and self.data3D[x][y][z]>self.statistical_interval[0]:
                            sum_value+=1
                            count=1
                if count==0:
                    if self.mode==1 or self.mode==4:
                        middle.append(0)
                    elif self.mode==2:
                        middle.append(2.65)
                    elif self.mode==3:
                        middle.append(2.65)
                else:
                    if self.mode==1 or self.mode==4:
                        middle.append((sum_value)/count)
                    elif self.mode==2:
                        middle.append((min_value)/count)
                    elif self.mode==3:
                        middle.append((max_value)/count)
            z_projection.append(middle)
        return z_projection
    def get_x_projection(self):
        x_projection=[]
        for z in range(len(self.data3D[0][0])):
            middle=[]
            for y in range(len(self.data3D[0])):
                sum_value=0
                count=0
                min_value=10000
                max_value=-10000
                for x in range(len(self.data3D)):
                    if self.mode==1:
                        if self.data3D[x][y][z]<2.7 or self.data3D[x][y][z]>100:
                            continue
                        # 平均值策略
                        sum_value+=self.data3D[x][y][z]-2.65
                        count+=1
                    elif self.mode==2:
                        if self.data3D[x][y][z]<1 or self.data3D[x][y][z]>100:
                            continue
                        #最小值策略
                        min_value=min(self.data3D[x][y][z],min_value)
                        count=1
                    elif self.mode==3:
                        if self.data3D[x][y][z]<2.7 or self.data3D[x][y][z]>100:
                            continue
                        #最大值策略
                        max_value=max(self.data3D[x][y][z],max_value)
                        count=1
                    elif self.mode==4:
                        #计数
                        if self.data3D[x][y][z]<self.statistical_interval[1] and self.data3D[x][y][z]>self.statistical_interval[0]:
                            sum_value+=1
                            count=1
                if count==0:
                    if self.mode==1 or self.mode==4:
                        middle.append(0)
                    elif self.mode==2:
                        middle.append(2.65)
                    elif self.mode==3:
                        middle.append(2.65)
                else:
                    if self.mode==1 or self.mode==4:
                        middle.append((sum_value)/count)
                    elif self.mode==2:
                        middle.append((min_value)/count)
                    elif self.mode==3:
                        middle.append((max_value)/count)
            x_projection.append(middle)
        return x_projection
    def get_y_projection(self):
        y_projection=[]
        for z in range(len(self.data3D[0][0])):
            middle=[]
            for x in range(len(self.data3D)):
                sum_value=0
                count=0
                min_value=10000
                max_value=-10000
                for y in range(len(self.data3D[0])):
                    if self.mode==1:
                        if self.data3D[x][y][z]<2.7 or self.data3D[x][y][z]>100:
                            continue
                        # 平均值策略
                        sum_value+=self.data3D[x][y][z]-2.65
                        count+=1
                    elif self.mode==2:
                        if self.data3D[x][y][z]<1 or self.data3D[x][y][z]>100:
                            continue
                        #最小值策略
                        min_value=min(self.data3D[x][y][z],min_value)
                        count=1
                    elif self.mode==3:
                        if self.data3D[x][y][z]<2.7 or self.data3D[x][y][z]>100:
                            continue
                        #最大值策略
                        max_value=max(self.data3D[x][y][z],max_value)
                        count=1
                    elif self.mode==4:
                        #计数
                        if self.data3D[x][y][z]<self.statistical_interval[1] and self.data3D[x][y][z]>self.statistical_interval[0]:
                            sum_value+=1
                            count=1
                if count==0:
                    if self.mode==1 or self.mode==4:
                        middle.append(0)
                    elif self.mode==2:
                        middle.append(2.65)
                    elif self.mode==3:
                        middle.append(2.65)
                else:
                    if self.mode==1 or self.mode==4:
                        middle.append((sum_value)/count)
                    elif self.mode==2:
                        middle.append((min_value)/count)
                    elif self.mode==3:
                        middle.append((max_value)/count)
            y_projection.append(middle)
        return y_projection
    
    def show_2D_data(self,data,amplification_factor=10):
        from PIL import Image
        import cv2
        #投影到0-255的空间
        min_value=10000
        max_value=-100000
        for da in data:
            min_value=min(min_value,min(da))
            max_value=max(max_value,max(da))
        print(min_value,max_value)
        diff=(max_value-min_value+0.00001)
        for i in range(len(data)):
            for j in range(len(data[i])):
              data[i][j]=(data[i][j]-min_value)/(diff)*255
         
        data_extend=np.zeros(shape=(len(data)*amplification_factor,len(data[0])*amplification_factor),dtype=float)
        for i in range(data_extend.shape[0]):
            for j in range(data_extend.shape[1]):
                data_extend[i][j]=data[int(i/amplification_factor)][int(j/amplification_factor)]
        data_extend=data_extend.astype(np.uint8)
        cv2.imshow("win", data_extend)
        cv2.waitKey(-1)
        
        
if __name__ == '__main__':
    base_path = r"E:\vscode\Muon_Imaging_Algorithm\data\Temp"
    # tool = Show_pred_obsd_derr(pred_obsd_derr_file=base_path + r"\pred_obsd_derr")
    # tool.show_pred_obs_obsderr()
    # tool.show_pred_obsd()
    # tool.show_pred_obsd_diff()
    norm_tool = Norm_tools(all_misfit_ms_smooth_file=base_path + r"\all_misfit_ms_smooth")
    norm_tool.show_norm()
    # beta_tool = Beta_log(base_path + r"\beta_misfit_norms")
    # beta_tool.show_Beta_search()
    # beta_tool.show()
    # beta_tool.show_similarL()
    # beta_tool.show_L()
    # tool = Show_data_distribution()
    # tool._show_res_file_num(r"E:\vscode\Muon_Imaging_Algorithm\data\output\rel\res")
