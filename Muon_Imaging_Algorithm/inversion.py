import sys
import os
# import scienceplots
# import tensorflow
rootpath=str(os.path.split(__file__)[0])
syspath=sys.path
sys.path=[]
sys.path.append(rootpath) #将工程根目录加入到python搜索路径中
 
sys.path.extend([rootpath+i for i in os.listdir(rootpath) if i[0]!="."]) #将工程目录下的一级目录添加到python搜索路径中
 
sys.path.extend(syspath)
# sys.path.append(os.environ.get('NOTEBOOK_ROOT'))

from InvDataFactory.DataManage import DataManager
from InvDataFactory.Setting import Setting
from InvDataTools.MeshTools import MeshTools
from InvDataTools.Topo_Tools_plus import Topo_Tools_plus
from InvDataTools.orter_tools.res_diff_Analysis_took import res_diff_Analysis, res_diff_Analysis_density_anomaly
from InvDataTools.paper_code.refs_bound_simple_build import refs_bounds_build,refs_bounds_build_three_moudle,refs_bounds_build_three_moudle_for_lgr_test


def start_inversion_by_setting(yaml_file):
    """
    根据配置文件yaml_file启动一个新的反演程序
    !将会覆盖内存中的Setting和DataManage对象,此方法是可重复调用的
    !不支持多线程,稍作修改可以支持
    Args:
        yaml_file (string): 配置文件路径
    """
    """更新反演中用到的单例对象"""
    Setting.get_instance(need_new=True,new_setting_file_path=yaml_file)
    DataManager.get_instance(new=True)
    """创建求解器并执行反演指令"""
    from InvSolver import Solver
    solver = Solver.Solver()
    # beta选择模式
    old_beta = 0
    old_misfit = 0
    count = 0
    while True:
        count+=1
        x, beta_info, misfit = solver.run()#根据配置文件执行反演
        solver.show_state()#显示反演程序的状态
        old_beta = beta_info[0]
        old_misfit = misfit
        print("本次使用的beta:",old_beta)
        print("old_misfit:%d, misfit:%d"%(old_misfit, misfit))
        print(count)
        if not beta_info[1]:
            break

    solver.show_state()#显示反演的状态
    solver.save_res()#存储反演结果
    solver.close()#结束反演
    
def build_refs_bounds_file_from_intersection_information(intersection_points_file,obs_file,mesh_file,refs_file,bounds_file,refs_value,bounds_value:tuple,default_refs_value=0,default_bounds_value:tuple=(-0.001,0.001)):
    """
    根据交点信息生成refs和bounds文件
    目前只支持将第1到2,第3到4....之间的数据修改为refs_value和bounds_value
    交点信息可以不是成对出现的
    !如果给定的文件不为空,将在此文件的基础上更新文件

    Args:
        intersection_points_file (String): 交点信息文件路径 必须存在
        obs_file (String): obs文件路径 必须存在
        mesh_file (String): mesh文件路径 必须存在
        refs_file (String): refs文件路径 存在将进行原地修改
        bounds_file (String): refs文件路径 存在将进行原地修改
        refs_value (float): 该类别refs的值
        bounds_value (tuple): 该类别bounds的值
        default_refs_value (int, optional): 新建文件时使用的默认值. Defaults to 0.
        default_bounds_value (tuple, optional): 新建文件时使用的默认值. Defaults to (-0.001,0.001).
    Raises:
        Exception: 数据类型错误
        Exception: default_bounds_value不是元组
        Exception: mesh文件路径不正确
        Exception: 数据数量不匹配
    """
    tool = Topo_Tools_plus(intersection_points_file,obs_file,mesh_file)
    tool.make_topo_j_by_points_info()#触发基于交点信息地分类计算
    tool.update_refs_bounds(mesh_file,refs_file,bounds_file,refs_value,bounds_value,default_refs_value,default_bounds_value)
    
    
    
def get_bottom_cells_xyz_coordinate(mesh_file,cells_xyz_coorfinate_file,strategy=(0,0,1),cells_3d_file=None):
    """获取成像区域最下边格子的成像

    Args:
        mesh_file (String): mesh文件的路径
        cells_xyz_coorfinate_file (String): 结果中地面格子坐标文件保存地址
        strategy (tuple, optional): 0 中点,1 左边/下边/前边 ,2 右边/上边/后边. Defaults to (0,0,1).
        cells_3d_file (String, optional): 是否选择生成三维演示图形,不为None该结果将保存在此位置. Defaults to None.
    """
    mesh_tool=MeshTools(mesh_file)
    res_file=open(cells_xyz_coorfinate_file,'w')
    xyz=mesh_tool.get_shape()
    z = xyz[2]
    j_s=set()
    for x in range(xyz[0]):
        for y in range(xyz[1]):
            # j_s.add(getj_from_xyz(shape=mesh_tool.get_shape(),xyz=(x+1,y+1,1)))
            p_xyz=mesh_tool.get_coordinates_form_xyz(x,y,1,strategy=strategy)
            res_file.write(str(p_xyz[0]))
            res_file.write(' ')
            res_file.write(str(p_xyz[1]))
            res_file.write(' ')
            res_file.write(str(p_xyz[2]))
            res_file.write('\n')
    res_file.close()
    """在三维软件中显示地面格子(值为1)"""
    if cells_3d_file is not None:
        res_3d_file=open(cells_3d_file,'w')
        for j in range(mesh_tool.cells_count()):
            j=j+1
            if j in j_s:
                res_3d_file.write('1')
                res_3d_file.write('\n')
            else:
                res_3d_file.write('-1')
                res_3d_file.write("\n")
        res_3d_file.close()
    
    
def Mesh_Coordinate_transformation(yaml_file):
    """三维坐标转换
    输入物理坐标(顺序依照meshTools.exe软件)显示对应点所在格子的详细信息包括编号j
    Args:
        yaml_file (String): 配置文件地址
    """
    import InvDataFactory.DataManage
    from InvDataTools.Jxyz_Tools import getj_from_xyz,getxyz_from_shape
    from InvSysTools.MyTools.myPrint import myPrint_Err,myPrint_Hint,myPrint_Success,myPrint_Wran
    """更新反演中用到的单例对象"""
    Setting.get_instance(need_new=True,new_setting_file_path=yaml_file)
    datamange=DataManager.get_instance(new=True)
    mesh = datamange.mesh
    while True:
        instruction=input("物理坐标转离散坐标输入xyz,离散坐标转体素中心坐标输入j,终止程序请输入exit:")
        if instruction=="xyz":
            while True:
                xyz = input("输入物理坐标,使用空格隔开,终止程序请输入exit:")
                if xyz == "exit":
                    break
                try:
                    xyz = [float(i) for i in xyz.split()]
                    res = mesh.discretize_Physical_coordinates(xyz[0], xyz[1], xyz[2], True)
                    j = getj_from_xyz(mesh.get_shape(), res)
                    myPrint_Hint("需要判断的数据", xyz)
                    myPrint_Hint("判断结果:", res)
                    myPrint_Hint("所处格子对应的代表坐标(0,0,0):" ,mesh.get_coordinates_form_xyz(res[0],res[1],res[2]))
                    myPrint_Hint("j:", j)
                    myPrint_Wran("对应的参考模型的密度、最小值、最大值",datamange.refs_tool.get_data()[j-1],datamange.bonds_tool.get_bonds_min()[j-1],datamange.bonds_tool.get_bonds_max()[j-1])
                except Exception as e:
                    myPrint_Err(xyz, "查询失败")
                    myPrint_Err(e)
        elif instruction=="j":
            while True:
                j_value = input("输入编号j,终止程序请输入exit:")
                if j_value == "exit":
                    break
                try:
                    myPrint_Hint("需要判断的数据", j_value)
                    j_value=int(j_value)
                    myPrint_Wran("对应的参考模型的密度、最小值、最大值",datamange.refs_tool.get_data()[j_value-1],datamange.bonds_tool.get_bonds_min()[j_value-1],datamange.bonds_tool.get_bonds_max()[j_value-1])
                    x,y,z=getxyz_from_shape(mesh.get_shape(),j=int(j_value))
                    myPrint_Hint("离散结果坐标:", x-1,y-1,z-1)
                    x_up=mesh.x_start_values[x-1]
                    x_low=mesh.x_end_values[x-1]
                    y_up=mesh.y_start_values[y-1]
                    y_low=mesh.y_end_values[y-1]
                    z_up=mesh.z_start_values[z-1]
                    z_low=mesh.z_end_values[z-1]
                    myPrint_Hint("取值范围(x,y,z): ",[x_low,x_up],[y_low,y_up],[z_low,z_up])
                    myPrint_Hint("中点坐标(x,y,z): ",(x_low+x_up)/2,(y_up+y_low)/2,(z_low+z_up)/2)
                except Exception as e:
                    myPrint_Err(j_value, "查询失败")
                    myPrint_Err(e)
        elif instruction=='exit':
            break
        else:
            myPrint_Wran('无效指令')
  
if __name__ == '__main__':
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\configuration.yml")
    # Mesh_Coordinate_transformation(r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\configuration.yml")
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\seed_simulation.yml")
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\data\backup\paper_simulation.yml")
    
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\configuration.yml")
    
#     """测试build_refs_bounds_file_from_intersection_information"""
    # intersection_points_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\inter1.txt"
    # obs_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\direction.dat"
    # mesh_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\ZZG5.msh"
    # refs_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\bkup.ZZG5_ref_1.den"
    # bounds_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\bkup.ZZG5_bnd_1.den"
    # refs_value=10
    # bounds_value=(9,11)
    # build_refs_bounds_file_from_intersection_information(intersection_points_file,obs_file,mesh_file,refs_file,bounds_file,refs_value,bounds_value)

    # res_diff_Analysis(r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\result\wall_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\moudles",
    #                   r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\result\wall_smooth"
    #                     ,r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\seed_opt.msh")
    # res_diff_Analysis_density_anomaly(r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\optimize_res\CRBOX5_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\moudles"
                        # ,r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\seed_opt.msh")
    # from InvDataTools.MeshTools import MeshTools
    # origin_of_coordinates=[0,1,-10]
    # x_directions=[[1,10],[1.1,10],[1,10],[1.0,10]]
    # y_directions=[[1,10],[1.1,10],[1,10],[1.0,10]]
    # z_directions=[[1,10],[1.1,10],[1,10],[1.0,10]]
    # mesh_file_path=r"M:\pycharm\Inversion\Temp\mesh"
    
    # res=MeshTools.generate_mesh_file(origin_of_coordinates,x_directions,y_directions,z_directions,mesh_file_path)
    
    # from InvDataTools.show_ray_trace_tools import show_ray_trace_tools
    # from InvDataTools.obs_tools import obs_tools
    # tool=show_ray_trace_tools([52,34,14])
    # obs_tool=obs_tools(r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\3_obs.dat")
    # detectoers=obs_tool.get_obs_id_count()
    # for i in range(1,len(detectoers)):
    #     detectoers[i]+=detectoers[i-1]
    # tool.mark_ray(file_ij=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\Gij",file_res=r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\output\ray_way",group_method=detectoers,ids=[i for i in range(0,detectoers[-1])])

    # from InvDataTools.Visibility_res_tools import show_2D_the_density_of_res
    # tool=show_2D_the_density_of_res(r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\result\seed_res7",r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\ZZG5.msh")
    # data=tool.get_z_projection()
    # tool.show_2D_data(data)
    # data=tool.get_x_projection()
    # tool.show_2D_data(data)
    # data=tool.get_y_projection()
    # tool.show_2D_data(data)
    
    
    # res_diff_Analysis(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\simulation\res_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\simulation\paper1alt\moudles",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\simulation\res_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\simulation\mesh.txt")
    
    #第二篇论文
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\paper_simulation.yml")
    # refs_bounds_build(r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\mesh.txt",r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\Paper1.6")
    ##多密度异常体
    # refs_bounds_build_three_moudle(r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\mesh.txt",r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\Paper3_box")
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\paper_simulation.yml")
    # Mesh_Coordinate_transformation(r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\paper_simulation.yml")
    
    # res_diff_Analysis(r"G:\重要\论文和材料\第二篇\实验数据\BCD-NES\2.4\seed_result",r"G:\重要\论文和材料\第二篇\实验数据\BCD-NES\2.4\input\moudles",r"G:\重要\论文和材料\第二篇\实验数据\3box\seed_res\seed_result_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\mesh.txt")

    #wall
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\configuration.yml")
    # Mesh_Coordinate_transformation(r"E:\vscode\Muon_Imaging_Algorithm\configuration.yml")
    
    #zzg_topo
    # Mesh_Coordinate_transformation(r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\configuration.yml")
    
    #SGDTest
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\paper_simulation.yml")
    
    #刘国瑞Test--2023年6月20日
    ##直接求解
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\lgr_new_moudle\lgr_new_test.yml")
    ##建立理论模型
    # refs_bounds_build_three_moudle_for_lgr_test(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\lgr_new_moudle\toy4.msh",r"E:\vscode\Muon_Imaging_Algorithm\data\output\lgr_new_moudle\res")
    ##zzgfz
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\lgr\paper_dafu\zzgfz\configuration.yml")
    # #坐标转换和查询工具
    # Mesh_Coordinate_transformation(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\lgr_new_moudle\lgr_new_test.yml")
    # res_diff_Analysis(r"E:\vscode\Muon_Imaging_Algorithm\data\output\lgr_new_moudle\seed_res\seed_result_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\output\lgr_new_moudle\res\moudles",r"E:\vscode\Muon_Imaging_Algorithm\data\output\lgr_new_moudle\seed_res\seed_result_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\lgr_new_moudle\toy2.msh")
    
    ####!第三篇论文
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\simulation\three_paper.yml")
    # refs_bounds_build(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\simulation\mesh.txt",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\simulation")
    # res_diff_Analysis(r"E:\vscode\Muon_Imaging_Algorithm\data\output\simulation\res_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\Paper1.6\moudles",r"E:\vscode\Muon_Imaging_Algorithm\data\output\simulation\res_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\paper2\Differences_in_contrast\data\mesh.txt")
    
    ##todo 三个异常体
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\data\paper3\paper3_box\paper3_simulation_3Box.yml")
    # res_diff_Analysis(r"E:\vscode\Muon_Imaging_Algorithm\data\paper3\paper3_box\res\res_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\paper3\paper3_box\moudles",r"E:\vscode\Muon_Imaging_Algorithm\data\paper3\paper3_box\res\res_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\paper3\paper3_box\mesh.txt")
    
    ##todo 城墙
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\data\paper3\wall\paper3_wall.yml")
    
    
    ##endtodo 城墙2
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\wall2\wall2.yml")
    
    ##endtodo ysy
    # start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\configuration.yml")
    # Mesh_Coordinate_transformation(r"E:\vscode\Muon_Imaging_Algorithm\configuration.yml")
    
    ##todo lgr 论文答复
    # Mesh_Coordinate_transformation(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\lgr\paper_dafu\configuration.yml")
    
    ##毕业论文
    # res_diff_Analysis(r"E:\vscode\Muon_Imaging_Algorithm\data\output\simulation\2.4\res_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\simulation\paper1alt\moudles",r"E:\vscode\Muon_Imaging_Algorithm\data\output\simulation\2.4\res_smooth",r"E:\vscode\Muon_Imaging_Algorithm\data\Input\simulation\mesh.txt")
    
    ##ykq
    start_inversion_by_setting(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\ykq\paper_simulation.yml")