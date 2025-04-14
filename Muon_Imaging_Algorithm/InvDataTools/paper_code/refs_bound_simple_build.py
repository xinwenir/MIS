# author:高金磊
# datetime:2022/6/29 17:33
from json import tool
import math


def refs_bounds_build(mesh_file,res_file_parent_directory):
    from InvDataTools import MeshTools
    mesh=MeshTools.MeshTools(mesh_file=mesh_file)
    from InvDataTools.Jxyz_Tools import getj_from_xyz
    refs_data=[0. for i in range(mesh.cells_count())]
    moudles = [0. for i in range(mesh.cells_count())]#模型的真实情况
    bounds_data=[[-0.01,0.01] for i in range(mesh.cells_count())]

    #paper1
    # x=[7,9]
    # y=[9,12]
    # z=[5,7]
    # corners = [[9, 9, 7], [7, 9, 7], [7, 9, 5],
    #            [9, 9, 5],
    #            [9, 12, 7], [7, 12, 7], [7, 12, 5],
    #            [9, 12, 5]
    #            ]

    #paper2
    x_min=6
    x_max=9
    y_min=8
    y_max=12
    z_min=4
    z_max=6
    corners = [[x_max, y_min, z_max], [x_min, y_min, z_max], [x_min, y_min, z_min],
               [x_max, y_min, z_min],
               [x_max, y_max, z_max], [x_min, y_max, z_max], [x_min, y_max, z_min],
               [x_max, y_max, z_min]
               ]
    #CRBOX
    # x_min=27.5
    # x_max=32.5
    # y_min=17
    # y_max=22
    # z_min=17.5
    # z_max=22.5
    # corners = [[x_max, y_min, z_max], [x_min, y_min, z_max], [x_min, y_min, z_min],
    #            [x_max, y_min, z_min],
    #            [x_max, y_max, z_max], [x_min, y_max, z_max], [x_min, y_max, z_min],
    #            [x_max, y_max, z_min]
    #            ]
    if True:#通过约束构建简单的反演问题
        from InvDataTools.show_ray_trace_tools import show_ray_trace_tools
        tool=show_ray_trace_tools(mesh.get_shape())
        inner_moudle_js=tool.get_j_from_scope_mesh(corners,mesh,r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\jxyz",show=False)
    else:
        inner_moudle_js=set()
    x1=set()
    x1.add(0)
    x1.add(1)
    x1.add(2)
    # x1.add(3)
    x1.add(len(mesh.get_xs())-1)
    x1.add(len(mesh.get_xs())-2)
    x1.add(len(mesh.get_xs()) - 3)
    # x1.add(len(mesh.get_xs()) - 4)
    y1 = set()
    y1.add(0)
    y1.add(1)
    y1.add(2)
    # y1.add(3)
    y1.add(len(mesh.get_ys()) - 1)
    y1.add(len(mesh.get_ys()) - 2)
    y1.add(len(mesh.get_ys()) - 3)
    # y1.add(len(mesh.get_ys()) - 4)

    z1 = set()
    z1.add(0)
    z1.add(1)
    z1.add(2)
    # z1.add(3)
    z1.add(len(mesh.get_zs())-1)
    z1.add(len(mesh.get_zs()) - 2)
    z1.add(len(mesh.get_zs()) - 3)
    # z1.add(len(mesh.get_zs()) - 4)
    for x in range(len(mesh.get_xs())):
        for y in range(len(mesh.get_ys())):
            for z in range(len(mesh.get_zs())):
                j = getj_from_xyz(mesh.get_shape(), (x+1, y+1, z+1))-1
                if x in x1 or y in y1 or z in z1:
                    refs_data[j]=0
                    moudles[j]=0
                    bounds_data[j]=[-0.1,0.1]
                else:
                    if z<15 or z>35: #探测器附近的格子正交性较差,为了避免射线痕迹lbfgs反演中对这种格子增加约束
                        moudles[j]=2.65
                        refs_data[j]=2.65
                        bounds_data[j] = [2.5, 2.7]
                        continue
                    if j+1 in inner_moudle_js:
                        moudles[j]=2.2
                        bounds_data[j] = [0, 2.7]
                    else:
                        moudles[j]=2.65
                        bounds_data[j] = [0, 2.7]
                    refs_data[j] = 2.65

                    # bounds_data[j] = [0,2.75]
    refs_file=open(res_file_parent_directory+r"\refs",'w')
    bounds_file=open(res_file_parent_directory+r"\bounds",'w')
    moudles_file=open(res_file_parent_directory+r"\moudles",'w')
    for i in refs_data:
        refs_file.write(str(i))
        refs_file.write("\n")
    for i in bounds_data:
        bounds_file.write(str(i[0]))
        bounds_file.write(" ")
        bounds_file.write(str(i[1]))
        bounds_file.write("\n")
    for i in moudles:
        moudles_file.write(str(i))
        moudles_file.write('\n')
    refs_file.close()
    bounds_file.close()
    moudles_file.close()
def refs_bounds_build_three_moudle(mesh_file,res_file_parent_directory):
    """论文中三个密度异常体的图形
    #! get_j_from_scope_mesh方法固有问题，处于体素表面的体素会被误判为不在体素内部
    Args:
        mesh_file (_type_): _description_
        res_file_parent_directory (_type_): _description_
    """
    from InvDataTools import MeshTools
    mesh=MeshTools.MeshTools(mesh_file=mesh_file)
    from InvDataTools.Jxyz_Tools import getj_from_xyz
    refs_data=[0. for i in range(mesh.cells_count())]
    moudles = [0. for i in range(mesh.cells_count())]#模型的真实情况
    bounds_data=[[-0.01,0.01] for i in range(mesh.cells_count())]

    #paper2——1.2
    y_min=13.5
    y_max=16.5
    x_min=5.5
    x_max=9.5
    z_min=6
    z_max=8
    corners1_2 = [[x_max, y_min, z_max], [x_min, y_min, z_max], [x_min, y_min, z_min],
               [x_max, y_min, z_min],
               [x_max, y_max, z_max], [x_min, y_max, z_max], [x_min, y_max, z_min],
               [x_max, y_max, z_min]
               ]

    if True:#通过约束构建简单的反演问题
        from InvDataTools.show_ray_trace_tools import show_ray_trace_tools
        tool=show_ray_trace_tools(mesh.get_shape())
        inner_moudle_js1_2=tool.get_j_from_scope_mesh(corners1_2,mesh,r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\jxyz",show=False)
    else:
        inner_moudle_js1_2=set()
        
        #paper2——1.5
    y_min=8.5
    y_max=11.5
    x_min=5.5
    x_max=9.5
    z_min=6
    z_max=8
    corners1_5 = [[x_max, y_min, z_max], [x_min, y_min, z_max], [x_min, y_min, z_min],
               [x_max, y_min, z_min],
               [x_max, y_max, z_max], [x_min, y_max, z_max], [x_min, y_max, z_min],
               [x_max, y_max, z_min]
               ]

    if True:#通过约束构建简单的反演问题
        from InvDataTools.show_ray_trace_tools import show_ray_trace_tools
        tool=show_ray_trace_tools(mesh.get_shape())
        inner_moudle_js1_5=tool.get_j_from_scope_mesh(corners1_5,mesh,r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\jxyz",show=False)
    else:
        inner_moudle_js1_5=set()
        #paper2——1.2
    y_min=3.5
    y_max=6.5
    x_min=5.5
    x_max=9.5
    z_min=6
    z_max=8
    corners1_8 = [[x_max, y_min, z_max], [x_min, y_min, z_max], [x_min, y_min, z_min],
               [x_max, y_min, z_min],
               [x_max, y_max, z_max], [x_min, y_max, z_max], [x_min, y_max, z_min],
               [x_max, y_max, z_min]
               ]

    if True:#通过约束构建简单的反演问题
        from InvDataTools.show_ray_trace_tools import show_ray_trace_tools
        tool=show_ray_trace_tools(mesh.get_shape())
        inner_moudle_js1_8=tool.get_j_from_scope_mesh(corners1_8,mesh,r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\jxyz",show=False)
    else:
        inner_moudle_js1_8=set()
        
    x1=set()
    x1.add(0)
    x1.add(1)
    x1.add(2)
    # x1.add(3)
    x1.add(len(mesh.get_xs())-1)
    x1.add(len(mesh.get_xs())-2)
    x1.add(len(mesh.get_xs()) - 3)
    # x1.add(len(mesh.get_xs()) - 4)
    y1 = set()
    y1.add(0)
    y1.add(1)
    y1.add(2)
    # y1.add(3)
    y1.add(len(mesh.get_ys()) - 1)
    y1.add(len(mesh.get_ys()) - 2)
    y1.add(len(mesh.get_ys()) - 3)
    # y1.add(len(mesh.get_ys()) - 4)

    z1 = set()
    z1.add(0)
    z1.add(1)
    z1.add(2)
    # z1.add(3)
    z1.add(len(mesh.get_zs())-1)
    z1.add(len(mesh.get_zs()) - 2)
    z1.add(len(mesh.get_zs()) - 3)
    # z1.add(len(mesh.get_zs()) - 4)
    for x in range(len(mesh.get_xs())):
        for y in range(len(mesh.get_ys())):
            for z in range(len(mesh.get_zs())):
                j = getj_from_xyz(mesh.get_shape(), (x+1, y+1, z+1))-1
                if x in x1 or y in y1 or z in z1:
                    refs_data[j]=0
                    moudles[j]=0
                    bounds_data[j]=[-0.1,0.1]
                else:
                    if z<20 or z>40: #探测器附近的格子正交性较差,为了避免射线痕迹lbfgs反演中对这种格子增加约束
                        moudles[j]=2.65
                        refs_data[j]=2.65
                        bounds_data[j] = [2.5, 2.7]
                        continue
                    if j+1 in inner_moudle_js1_2:
                        moudles[j]=1.2
                        bounds_data[j] = [0, 2.7]
                    elif j+1 in inner_moudle_js1_5:
                        moudles[j]=1.5
                        bounds_data[j] = [0, 2.7]
                    elif j+1 in inner_moudle_js1_8:
                        moudles[j]=1.8
                        bounds_data[j] = [0, 2.7]
                    else:
                        moudles[j]=2.65
                        bounds_data[j] = [0, 2.7]
                    refs_data[j] = 2.65

                    # bounds_data[j] = [0,2.75]
    refs_file=open(res_file_parent_directory+r"\refs",'w')
    bounds_file=open(res_file_parent_directory+r"\bounds",'w')
    moudles_file=open(res_file_parent_directory+r"\moudles",'w')
    for i in refs_data:
        refs_file.write(str(i))
        refs_file.write("\n")
    for i in bounds_data:
        bounds_file.write(str(i[0]))
        bounds_file.write(" ")
        bounds_file.write(str(i[1]))
        bounds_file.write("\n")
    for i in moudles:
        moudles_file.write(str(i))
        moudles_file.write('\n')
    refs_file.close()
    bounds_file.close()
    moudles_file.close()
    
    
def refs_bounds_build_three_moudle_for_lgr_test(mesh_file,res_file_parent_directory):
    """刘国睿测试中三个密度异常体的模型

    Args:
        mesh_file (_type_): _description_
        res_file_parent_directory (_type_): _description_
    """
    from InvDataTools import MeshTools
    mesh=MeshTools.MeshTools(mesh_file=mesh_file)
    from InvDataTools.Jxyz_Tools import getj_from_xyz
    refs_data=[0. for i in range(mesh.cells_count())]
    moudles = [0. for i in range(mesh.cells_count())]#模型的真实情况
    bounds_data=[[-0.01,0.01] for i in range(mesh.cells_count())]

    #one
    x_min=-35.0001
    x_max=-25.0001
    y_min=-0.0001
    y_max=20.0001
    z_min=40.0001
    z_max=55.0001
    corners1_2 = [[x_max, y_min, z_max], [x_min, y_min, z_max], [x_min, y_min, z_min],
               [x_max, y_min, z_min],
               [x_max, y_max, z_max], [x_min, y_max, z_max], [x_min, y_max, z_min],
               [x_max, y_max, z_min]
               ]

    if True:#通过约束构建简单的反演问题
        from InvDataTools.show_ray_trace_tools import show_ray_trace_tools
        tool=show_ray_trace_tools(mesh.get_shape())
        inner_moudle_js1_2=tool.get_j_from_scope_mesh(corners1_2,mesh,r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\jxyz",show=False)
    else:
        inner_moudle_js1_2=set()
        
    x_min=-5
    x_max=5
    y_min=0.0001
    y_max=20
    z_min=40
    z_max=55
    corners1_5 = [[x_max, y_min, z_max], [x_min, y_min, z_max], [x_min, y_min, z_min],
               [x_max, y_min, z_min],
               [x_max, y_max, z_max], [x_min, y_max, z_max], [x_min, y_max, z_min],
               [x_max, y_max, z_min]
               ]

    if True:#通过约束构建简单的反演问题
        from InvDataTools.show_ray_trace_tools import show_ray_trace_tools
        tool=show_ray_trace_tools(mesh.get_shape())
        inner_moudle_js1_5=tool.get_j_from_scope_mesh(corners1_5,mesh,r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\jxyz",show=False)
    else:
        inner_moudle_js1_5=set()
        #paper2——1.2
    x_min=25
    x_max=35
    y_min=0.0001
    y_max=20
    z_min=40
    z_max=55
    corners1_8 = [[x_max, y_min, z_max], [x_min, y_min, z_max], [x_min, y_min, z_min],
               [x_max, y_min, z_min],
               [x_max, y_max, z_max], [x_min, y_max, z_max], [x_min, y_max, z_min],
               [x_max, y_max, z_min]
               ]

    if True:#通过约束构建简单的反演问题
        from InvDataTools.show_ray_trace_tools import show_ray_trace_tools
        tool=show_ray_trace_tools(mesh.get_shape())
        inner_moudle_js1_8=tool.get_j_from_scope_mesh(corners1_8,mesh,r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\jxyz",show=False)
    else:
        inner_moudle_js1_8=set()
        
    x1=set()
    x1.add(0)
    # x1.add(1)
    # x1.add(2)
    # # x1.add(3)
    x1.add(len(mesh.get_xs())-1)
    # x1.add(len(mesh.get_xs())-2)
    # x1.add(len(mesh.get_xs()) - 3)
    # # x1.add(len(mesh.get_xs()) - 4)
    y1 = set()
    y1.add(0)
    # y1.add(1)
    # y1.add(2)
    # # y1.add(3)
    y1.add(len(mesh.get_ys()) - 1)
    # y1.add(len(mesh.get_ys()) - 2)
    # y1.add(len(mesh.get_ys()) - 3)
    # # y1.add(len(mesh.get_ys()) - 4)

    z1 = set()
    z1.add(0)
    z1.add(1)
    # z1.add(2)
    # z1.add(3)
    z1.add(len(mesh.get_zs())-1)
    # z1.add(len(mesh.get_zs()) - 2)
    # z1.add(len(mesh.get_zs()) - 3)
    # # z1.add(len(mesh.get_zs()) - 4)
    for x in range(len(mesh.get_xs())):
        for y in range(len(mesh.get_ys())):
            for z in range(len(mesh.get_zs())):
                j = getj_from_xyz(mesh.get_shape(), (x+1, y+1, z+1))-1
                if x in x1 or y in y1 or z in z1:
                    refs_data[j]=0
                    moudles[j]=0
                    bounds_data[j]=[-0.1,0.1]
                else:
                    # if z<20 or z>40: #探测器附近的格子正交性较差,为了避免射线痕迹lbfgs反演中对这种格子增加约束
                    #     moudles[j]=2.65
                    #     refs_data[j]=2.65
                    #     bounds_data[j] = [2.5, 2.7]
                    #     continue
                    if j+1 in inner_moudle_js1_2:
                        moudles[j]=3.3
                        bounds_data[j] = [2.6, 4.1]
                    elif j+1 in inner_moudle_js1_5:
                        moudles[j]=2.0
                        bounds_data[j] = [2.6, 4.1]
                    elif j+1 in inner_moudle_js1_8:
                        moudles[j]=1.4
                        bounds_data[j] = [2.6, 4.1]
                    else:
                        moudles[j]=2.65
                        bounds_data[j] = [2.6, 4.1]
                    refs_data[j] = 2.65

                    # bounds_data[j] = [0,2.75]
    refs_file=open(res_file_parent_directory+r"\refs",'w')
    bounds_file=open(res_file_parent_directory+r"\bounds",'w')
    moudles_file=open(res_file_parent_directory+r"\moudles",'w')
    for i in refs_data:
        refs_file.write(str(i))
        refs_file.write("\n")
    for i in bounds_data:
        bounds_file.write(str(i[0]))
        bounds_file.write(" ")
        bounds_file.write(str(i[1]))
        bounds_file.write("\n")
    for i in moudles:
        moudles_file.write(str(i))
        moudles_file.write('\n')
    refs_file.close()
    bounds_file.close()
    moudles_file.close()
    
from InvDataTools.MeshTools import MeshTools
def simulation_moudles(inner_moudle_points,mesh:MeshTools):
    """由于判断格式是不是这个区域的方法支持不均匀划分故此方法失去价值
    建议直接使用get_j_from_scope_mesh
    此方法对于均匀划分的仍然适用且具有便捷性

    Args:
        inner_moudle_points (_type_): _description_
        mesh (MeshTools): _description_

    Returns:
        _type_: _description_
    """
    from InvDataTools.show_ray_trace_tools import show_ray_trace_tools
    tool=show_ray_trace_tools(mesh.get_shape())

    xcmin ,ycmin , zcmin = mesh.get_xyz_start()
    x_step = mesh.get_xs()[0]
    y_step = mesh.get_ys()[0]
    z_step = mesh.get_zs()[0]
    # tool.get_j_from_scope()
    cells = tool.get_j_from_scope(inner_moudle_points, xcmin, ycmin, zcmin, x_step, y_step, z_step,
                                  r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\jxyz", show=False)
    return cells





# import random
# import math

# a = [random.randint(0, 10) for t in range(20)]
# b = [random.randint(0, 10) for t in range(20)]

