# author:高金磊
# datetime:2022/9/19 15:53

from InvSolver.Seed_algorithm.Tools import data_tools
from InvSolver.Seed_algorithm.Tools import data_tools
from InvSolver.Seed_algorithm.paper.Visualization_tools import Record_Tools
modle="dev" #paper--论文环境（方法较新）；dev--普通开发环境,显示结果 -- 显示迭代指标曲线
if modle=="paper":
    from InvSolver.Seed_algorithm.paper.Solver import Solver
    solver=Solver()
    solver.BCD_solver()

    # data_file_list=[]
    # path=r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\Paper1.6\res\seed_res"
    # # path=r"G:\重要\论文和材料\第二篇\实验数据\无噪声成像结果"
    # for i in range(1,8):
    #     data_file_list.append(path+r"\seed_res%s"%(str(i)))
    # for i in range(0,20):
    #     data_file_list.append(path+ r"\seed_res_smmoth_%s"%(str(i)))
    # # data_file_list.append(r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\Paper1.6\res\seed_res\seed_result_smooth")
    # for data_file in data_file_list:
    #     tool.add_data(data_file)
    # tool.calculate_indicators(1.6,Abnormal_body_rho_interval=0.6)
    # tool.show_res_data()

    # #手动平滑结果
    # data_tool=data_tools(r"E:\vscode\Muon_Imaging_Algorithm\data\output\lgr_new_moudle\seed_res\part1\seed_result", r"E:\vscode\Muon_Imaging_Algorithm\data\Input\lgr_new_moudle\toy3.msh",
    #                      r"E:\vscode\Muon_Imaging_Algorithm\data\Input\lgr_new_moudle\toy3_ref1.den")
    # data_tool.output_smooth_res(r"E:\vscode\Muon_Imaging_Algorithm\data\output\lgr_new_moudle\seed_res\part1\seed_result_smooth")

    
elif modle=="dev": #注意，建议将setting.py中的punishment_diff_value_multiple设置为0 其功能用threshold代替。将punishment_Search_distancevalue_multiple设置为0其功能用max_count代替
    from InvSolver.Seed_algorithm.Solver import Solver
    from InvSolver.Seed_algorithm.Tools import data_tools
    solver=Solver()
    solver.BCD_solver()
elif modle=="显示结果":
    tool=Record_Tools(r"G:\重要\论文和材料\第二篇\实验数据\BCD-NES\input16\moudles",
                r"G:\重要\论文和材料\第二篇\实验数据\mesh.txt",
                r"G:\重要\论文和材料\第二篇\实验数据\BCD-NES\无噪声成像结果\loger.txt")
    #计算模型的查准率查全率误差
    data_file_list=[]
    # pre_moudle=r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\Paper1.6\res\seed_res\seed_result"
    # data_file_list.append(pre_moudle)
    
    path=r"G:\重要\论文和材料\第二篇\实验数据\BCD-NES\无噪声成像结果"
    for i in range(1,8):
        data_file_list.append(path+r"\seed_res%s"%(str(i)))
    for i in range(0,34):
        data_file_list.append(path+ r"\seed_res_smmoth_%s"%(str(i)))
    for data_file in data_file_list:
        tool.add_data(data_file)
    tool.add_data(r"G:\重要\论文和材料\第二篇\实验数据\BCD-NES\无噪声成像结果\seed_result")
    tool.calculate_indicators([1.6],Abnormal_body_rho_interval=0.6)
    tool.show_res_data(show=True)
else:
    data_tool=data_tools(r"E:\vscode\Muon_Imaging_Algorithm\data\output\lgr\paper_dafu\sim1\seed_result",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\lgr_new_moudle\toy4.msh",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\lgr_new_moudle\toy4_ref1.den")
    data_tool.output_smooth_res(r"E:\vscode\Muon_Imaging_Algorithm\data\output\lgr\paper_dafu\test.den")
    print("error modle")
