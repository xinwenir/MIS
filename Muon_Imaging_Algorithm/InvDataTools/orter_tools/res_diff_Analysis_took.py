# author:高金磊
# datetime:2022/7/29 10:19
import math

from InvDataTools.MeshTools import MeshTools
import matplotlib.pyplot as plt
import matplotlib
# import scienceplots
def res_diff_Analysis(res_file,moudle_file,refs2_file,mesh_file):
    ignore=0.15

    res_data=open(res_file,'r').readlines()
    res_data1=open(refs2_file,'r').readlines()
    moudle_data=open(moudle_file,'r').readlines()
    middle=[]
    for i in res_data:
        middle.append(float(i))
    res_data=middle
    middle=[]
    for i in res_data1:
        middle.append(float(i))
    res_data1=middle
    middle=[]
    for i in moudle_data:
        middle.append(float(i))
    moudle_data=middle

    mesh=MeshTools(mesh_file)
    shape=mesh.get_shape()
    from InvDataTools.Jxyz_Tools import getj_from_xyz
    x_cells_sum=[]
    x_cells_sum1=[]
    x_cells_sum2=[]
    for x in range(shape[0]):
        middle=0
        middle1=0
        middle2=0
        for z in range(shape[2]):
            for y in range(shape[1]):
                j=getj_from_xyz(shape,(x+1,y+1,z+1))-1
                middle += res_data[j] if abs(res_data[j]-2.65)>ignore else 2.65
                middle1+=moudle_data[j] if abs(moudle_data[j]-2.65)>ignore else 2.65
                middle2 += res_data1[j] if abs(res_data1[j]-2.65)>ignore else 2.65
        if x not in [0,1,2,3,shape[0]-1,shape[0]-2,shape[0]-3,shape[0]-4]:
            x_cells_sum.append(middle/((shape[2]-6)*(shape[1]-6)))
            x_cells_sum1.append(middle1/((shape[2]-6)*(shape[1]-6)))
            x_cells_sum2.append(middle2/((shape[2]-6)*(shape[1]-6)))
    print(Cov_array(x_cells_sum,x_cells_sum1))
    y_cells_sum = []
    y_cells_sum1 = []
    y_cells_sum2 = []
    for y in range(shape[1]):
        middle = 0
        middle1 = 0
        middle2=0
        for z in range(shape[2]):
            for x in range(shape[0]):
                j = getj_from_xyz(shape, (x + 1, y + 1, z + 1)) - 1
                middle += res_data[j] if abs(res_data[j]-2.65)>ignore else 2.65
                middle1 += moudle_data[j]
                middle2 += res_data1[j] if abs(res_data1[j]-2.65)>ignore else 2.65
        if y not in [0, 1,2,3, shape[1] - 1, shape[1] - 2,shape[1] - 3,shape[1] - 4]:
            y_cells_sum.append(middle / ((shape[2] - 6) * (shape[0] - 6)))
            y_cells_sum1.append(middle1 / ((shape[2] - 6) * (shape[0] - 6)))
            y_cells_sum2.append(middle2 / ((shape[2] - 6) * (shape[0] - 6)))
    print(Cov_array(y_cells_sum, y_cells_sum1))
    z_cells_sum = []
    z_cells_sum1 = []
    z_cells_sum2 = []
    for z in range(shape[2]):
        middle = 0
        middle1 = 0
        middle2 = 0
        for y in range(shape[1]):
            for x in range(shape[0]):
                j = getj_from_xyz(shape, (x + 1, y + 1, z + 1)) - 1
                middle += res_data[j] if abs(res_data[j]-2.65)>(ignore) else 2.65
                middle1 += moudle_data[j]
                middle2 += res_data1[j] if abs(res_data1[j]-2.65)>(ignore) else 2.65
        if z not in [0, 1,2,3, shape[2] - 1, shape[2] - 2,shape[2]-3,shape[2]-4]:
            z_cells_sum.append(middle / ((shape[1] - 6) * (shape[0] - 6)))
            z_cells_sum1.append(middle1 / ((shape[1] - 6) * (shape[0] - 6)))
            z_cells_sum2.append(middle2 / ((shape[1] - 6) * (shape[0] - 6)))
    print(Cov_array(z_cells_sum, z_cells_sum1))


    #计算理论模型的

    plt.style.use(['science','ieee', 'no-latex'])
    matplotlib.rc("font", family='FangSong')  # 使用代码帮助matplotlib识别中文字体仿宋
    plt.rcParams['savefig.dpi'] = 500  # 图片像素
    # plt.rcParams['figure.dpi'] = 150  # 分辨率
    plt.rcParams['figure.figsize']=(2.9, 1.6)
    ymin=2.4
    ymax=2.68
    fig, ax = plt.subplots()
    # ax.set_ylim(ymin=ymin,ymax=ymax)
    ax.plot([i for i in range(len(x_cells_sum))],x_cells_sum ,label="预测模型")
    ax.plot([i for i in range(len(x_cells_sum1))],x_cells_sum1,label="理论模型")
    # ax.plot([i for i in range(len(x_cells_sum2))],x_cells_sum2 ,label="预测模型_LBFGS")

    # ax.plot([0,len(x_cells_sum)-1],[2.65,2.65],c='y',label="参考模型")
    ax.legend(
        # bbox_to_anchor=(1, 1),  # 图例边界框起始位置
               loc="lower left",  # 图例的位置
               ncol=1,  # 列数
               mode="None",  # 当值设置为“expend”时，图例会水平扩展至整个坐标轴区域
               borderaxespad=0,  # 坐标轴和图例边界之间的间距
            #    title="模型",  # 图例标题
               shadow=False,  # 是否为线框添加阴影
               fancybox=True)  # 线框圆角处理参数
    # ax.grid()
    ax.set_xlim(left=0)
    # ax.set_ylim(bottom=2.4,top=2.7)
    # 设置图表标题并给坐标轴加上标签。
    # ax.set_title("预测模型、理论模型、参考模型x方向的切片密度对比", fontsize=14)
    ax.set_xlabel("东西方向的切片编号")
    ax.set_ylabel("平均密度(x10$^3$kg/m$^3$)")#平均密度
    # 设置刻度标记的大小
    ax.tick_params(axis='both')

    fig, ax = plt.subplots()
    # ax.set_ylim(ymin=ymin,ymax=ymax)

    ax.plot([i for i in range(len(y_cells_sum))], y_cells_sum, label="预测模型")
    ax.plot([i for i in range(len(y_cells_sum1))], y_cells_sum1,label="理论模型")
    # ax.plot([i for i in range(len(y_cells_sum))], y_cells_sum, linewidth=1, c='r',label="结合使用分块坐标下降法")
    # ax.plot([0, len(y_cells_sum)-1], [2.65, 2.65], c='y',label="参考模型")
    ax.legend(
        # bbox_to_anchor=(1, 1),  # 图例边界框起始位置
              loc="lower left",  # 图例的位置
              ncol=1,  # 列数
              mode="None",  # 当值设置为“expend”时，图例会水平扩展至整个坐标轴区域
              borderaxespad=0,  # 坐标轴和图例边界之间的间距
            #   title="模型",  # 图例标题
              shadow=False,  # 是否为线框添加阴影
              fancybox=True)  # 线框圆角处理参数
    # ax.grid()
    ax.set_xlim(left=0)
    # ax.set_ylim(bottom=2.55, top=2.7)
    # 设置图表标题并给坐标轴加上标签。
    # ax.set_title("预测模型、理论模型、参考模型y方向的切片密度对比", fontsize=14)
    ax.set_xlabel("南北方向的切片编号")
    ax.set_ylabel("平均密度(x10$^3$kg/m$^3$)")
    # ax.set_ylabel("g/cm$^3$", fontsize=14)#平均密度
    # 设置刻度标记的大小
    ax.tick_params(axis='both')

    fig, ax = plt.subplots()
    # ax.set_ylim(ymin=ymin,ymax=ymax)
    ax.plot([i for i in range(len(z_cells_sum))], z_cells_sum,label="预测模型")
    ax.plot([i for i in range(len(z_cells_sum1))], z_cells_sum1,label="理论模型")
    # ax.plot([i for i in range(len(z_cells_sum2))], z_cells_sum2, linewidth=1, c='b',label="预测模型_LBFGS")
    # ax.plot([0, len(z_cells_sum)-1], [2.65, 2.65], c='y',label="参考模型")
    ax.legend(bbox_to_anchor=(1, 1),  # 图例边界框起始位置
              loc="lower right",  # 图例的位置
              ncol=1,  # 列数
              mode="None",  # 当值设置为“expend”时，图例会水平扩展至整个坐标轴区域
              borderaxespad=0,  # 坐标轴和图例边界之间的间距
            #   title="模型",  # 图例标题
              shadow=False,  # 是否为线框添加阴影
              fancybox=True)  # 线框圆角处理参数
    # ax.grid()
    ax.set_xlim(left=0)
    # ax.set_ylim(bottom=2.55, top=2.7)
    # 设置图表标题并给坐标轴加上标签。
    # ax.set_title("预测模型、理论模型、参考模型z方向的切片密度对比", fontsize=14)
    ax.set_xlabel("垂直方向的切片编号")
    ax.set_ylabel("平均密度(x10$^3$kg/m$^3$)")
    # ax.set_ylabel("x10$^3$kg/m$^3$", fontsize=14)#平均密度
    # 设置刻度标记的大小
    ax.tick_params(axis='both')

    plt.show()

def Cov(res_file,moudle_file,mesh_file):
    res_data = open(res_file, 'r').readlines()
    moudle_data = open(moudle_file, 'r').readlines()
    middle = []
    for i in res_data:
        middle.append(float(i))
    res_data = middle
    middle = []
    for i in moudle_data:
        middle.append(float(i))
    moudle_data = middle
    return Cov_array(res_data,moudle_data)

# moudle_data = open(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\simulation\moudles", 'r').readlines()
# middle = []
# for i in moudle_data:
#     middle.append(float(i))
# moudle_data = middle
def Cov_array(res_data,moudle_data):
    import numpy as np
    X = np.hstack((np.array(res_data)[:,np.newaxis], np.array(moudle_data)[:,np.newaxis]))
    return np.cov(X.T)[0][1]/math.pow(np.var(np.array(res_data))*np.var(np.array(moudle_data)),0.5)

def res_diff_Analysis_density_anomaly(res_file,moudle_file,mesh_file):
    """
    统计结果的准确率__针对seed

    Args:
        res_file (_type_): _description_
        moudle_file (_type_): _description_
        mesh_file (_type_): _description_
    """
    
    ignore=0.5

    res_data=open(res_file,'r').readlines()
    moudle_data=open(moudle_file,'r').readlines()
    middle=[]
    for i in res_data:
        middle.append(float(i))
    res_data=middle
    middle=[]
    for i in moudle_data:
        middle.append(float(i))
    moudle_data=middle

    mesh=MeshTools(mesh_file)
    shape=mesh.get_shape()
    tt=0
    tf=0
    ff=0
    ft=0
    for i in range(mesh.cells_count()):
        if moudle_data[i]>=2.65 or moudle_data[i]==0:#不是密度异常区的密度
            if res_data[i]==0 or res_data[i]>2.3:
                ff+=1
            else:
                ft+=1
        else:
            if res_data[i]==0 or res_data[i]>2.3:
                tf+=1
            else:
                tt+=1
    print("是密度异常_识别正确:%s,%s"%(str(tt),str(tt/(tt+tf))))
    print("是密度异常_识别错误:%s,%s"%(str(tf),str(tf/(tt+tf))))
    print("不是密度异常_识别正确:%s,%s"%(str(ff),str(ff/(ff+ft))))
    print("不是密度异常_识别错误:%s,%s"%(str(ft),str(ft/(ff+ft))))
            
    
if __name__ == '__main__':
    # res_diff_Analysis(r"E:\vscode\Muon_Imaging_Algorithm\data\output\rel\res_smooth_diff",r"E:\vscode\Muon_Imaging_Algorithm\data\output\rel\res_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\real_data\17_58MaMian.msh")
    # res_diff_Analysis(r"E:\vscode\Muon_Imaging_Algorithm\data\output\simulation\Seed_res\seed_res22",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\simulation\moudles",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\simulation\mesh.txt")
    res_diff_Analysis_density_anomaly(r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\result\wall_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\simulation\moudles",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\simulation\mesh.txt")

