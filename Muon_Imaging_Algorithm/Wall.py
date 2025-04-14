# author:高金磊
# datetime:2022/8/5 16:28
#城墙相关的小脚本
import tqdm

from InvDataTools.Air_j import Air_j
from InvDataTools.Bonds_tools import Bonds_tool
from InvDataTools.Jxyz_Tools import getxyz_from_shape, getj_from_xyz
from InvDataTools.MeshTools import MeshTools
from InvDataTools.obs_tools import obs_tools
from InvDataTools.ref_tools import Ref_tools
from InvDataTools.show_ray_trace_tools import show_ray_trace_tools, show_points_3D
from InvSysTools.MyTools import myPrint

if __name__ == '__main__':
    model=8
    if model==1:
        obs_tool1=obs_tools(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\real_data\27_obs.dat")
        all_rays=obs_tool1.get_data()
        obs_tool2=obs_tools(r"E:\vscode\Muon_Imaging_Algorithm\data\no_name\misfit_big100_det5.dat")
        rays_set=set()
        count=0
        # for data in obs_tool2.get_data():
        #     middle=""
        #     for datum in data:
        #         middle+=str(datum)+" "
        #     rays_set.add(middle)
        for i in range(len(all_rays)):
            middle=""
            for datum in all_rays[i]:
                middle+=str(datum)+" "
            if float(all_rays[i][5])>-0.02 and int(all_rays[i][0])==5:
                rays_set.add(middle)
            # if middle in rays_set:
            #     count+=1
            #     all_rays[i][5]=str(float(all_rays[i][5])+0.314)
        file_obj=open(r"E:\vscode\Muon_Imaging_Algorithm\data\no_name\obs_del_100.dat",'w')
        file_obj.write(str(len(all_rays)))
        file_obj.write('\n')
        for all_ray in all_rays:
            middle = ""
            for datum in all_ray:
                middle += str(datum)+" "
            if middle in rays_set:
                continue
            file_obj.write(middle)
            file_obj.write('\n')
        file_obj.close()
        print(count)
        print(len(rays_set))
    elif model == 3:
        """以下代码用来标记加厚的城墙2022-07-30"""
        mesh = MeshTools(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\real_data\17_58MaMian.msh")
        tool = show_ray_trace_tools(mesh.get_shape())

        # x_range=[15,22.5]
        # y_range=[-12.7,-10.5]
        # z_range=[0.5,3.5]
        # corners = [[max(x_range), min(y_range), max(z_range)], [min(x_range), min(y_range), max(z_range)], [min(x_range), min(y_range), min(z_range)],
        #            [max(x_range), min(y_range), min(z_range)],
        #            [max(x_range), max(y_range), max(z_range)], [min(x_range), max(y_range), max(z_range)], [min(x_range), max(y_range), min(z_range)],
        #            [max(x_range), max(y_range), min(z_range)]
        #            ]
        point2 = [21.4972, -1.4525, 10.57]
        point3 = [23.25, 0.31, 0]
        point6 = [21.4972, -2.43, 10.569]
        point7 = [23.2, -1.6, 0]

        point1 = [1.81, -1.69, 10.656]
        point4 = [0, 0, 0]
        point5 = [1.935, -2.43, 10.57]
        point8 = [0.25, -1.9, 0]

        corners = [point6, point5, point8, point7,
                   point2, point1, point4, point3]
        show_points_3D(data=corners)

        cells = tool.get_j_from_scope_mesh(corners, mesh, show=False)
        myPrint.myPrint_Hint("区域格子的ids:", cells)

        point2 = [21.75, -12.74, 10.58]
        point3 = [23.2, -12.74, 0.0001]
        point6 = [21.5, -12.74, 10.58]
        point7 = [21.5, -12.74, 0.0001]

        point1 = [21.75, -2.43, 10.58]
        point4 = [23.2, -2.0, 0.0001]
        point5 = [21.5, -2.43, 10.58]
        point8 = [21.5, -2.0, 0.0001]

        corners = [point1, point5, point8, point4,
                   point2, point6, point7, point3]
        show_points_3D(data=corners)

        cells1 = tool.get_j_from_scope_mesh(corners, mesh, show=False)
        myPrint.myPrint_Hint("区域格子的ids:", cells1)

        cells = cells.union(cells1)

        point2, point3, point6, point7, point1, point4, point5, point8 = [1.95, -2.43, 10.4], \
                                                                         [0.25, -2.0, 0.0001], \
                                                                         [2.2, -2.0, 10.4], \
                                                                         [2.2, -2.43, 0.0001], \
                                                                         [1.953, -12.99, 10.3], \
                                                                         [0.25, -12.99, 0.0001], \
                                                                         [2.2, -12.99, 10.4], \
                                                                         [2.2, -12.99, 0.0001]

        corners = [point5, point1, point4, point8,
                   point6, point2, point3, point7]
        show_points_3D(data=corners)
        cells2 = tool.get_j_from_scope_mesh(corners, mesh, show=False)
        myPrint.myPrint_Hint("区域格子的ids:", cells2)

        cells = cells.union(cells2)

        #####女墙,上表面,小房子等所有10.6m之上的
        x_range = [-100, 100]
        y_range = [-100, 100]
        z_range = [10.6, 13]
        corners = [[max(x_range), min(y_range), max(z_range)], [min(x_range), min(y_range), max(z_range)],
                   [min(x_range), min(y_range), min(z_range)],
                   [max(x_range), min(y_range), min(z_range)],
                   [max(x_range), max(y_range), max(z_range)], [min(x_range), max(y_range), max(z_range)],
                   [min(x_range), max(y_range), min(z_range)],
                   [max(x_range), max(y_range), min(z_range)]
                   ]
        cells3 = tool.get_j_from_scope_mesh(corners, mesh, show=False)
        myPrint.myPrint_Hint("区域格子的ids:", cells3)
        cells = cells.union(cells3)
        #####小房子
        x_range = [2.454, 4.2]
        y_range = [-11.48, -13.2]
        z_range = [12.1, 10.58]
        corners = [[max(x_range), min(y_range), max(z_range)], [min(x_range), min(y_range), max(z_range)],
                   [min(x_range), min(y_range), min(z_range)],
                   [max(x_range), min(y_range), min(z_range)],
                   [max(x_range), max(y_range), max(z_range)], [min(x_range), max(y_range), max(z_range)],
                   [min(x_range), max(y_range), min(z_range)],
                   [max(x_range), max(y_range), min(z_range)]
                   ]
        show_points_3D(data=corners)
        cells4 = tool.get_j_from_scope_mesh(corners, mesh, show=False)
        myPrint.myPrint_Hint("区域格子的ids:", cells4)
        ####后城墙---范围有点大,先复制再用精细的覆盖掉具体细节
        x_range = [-100, 100]
        y_range = [-9.73, -40]
        z_range = [0, 10.7]
        corners = [[max(x_range), min(y_range), max(z_range)], [min(x_range), min(y_range), max(z_range)],
                   [min(x_range), min(y_range), min(z_range)],
                   [max(x_range), min(y_range), min(z_range)],
                   [max(x_range), max(y_range), max(z_range)], [min(x_range), max(y_range), max(z_range)],
                   [min(x_range), max(y_range), min(z_range)],
                   [max(x_range), max(y_range), min(z_range)]
                   ]
        show_points_3D(data=corners)
        middle = tool.get_j_from_scope_mesh(corners, mesh, show=False)
        myPrint.myPrint_Hint("区域格子的ids:", middle)

        cells5 = set()
        xlen = mesh.get_xs()
        ylen = mesh.get_ys()
        for j in tqdm.tqdm(middle):
            x, y, z = getxyz_from_shape(mesh.get_shape(), j)
            if xlen[x - 1] > 0.3 or ylen[y - 1] > 0.3:
                cells5.add(j)

        # 剔除空气
        air_j = Air_j().get_air_j_from_file(r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\unneed_j")
        cells -= set(air_j)
        cells5 -= set(air_j)
        tool.draw_points(cells, res_file=r"E:\vscode\Muon_Imaging_Algorithm\data\output\rel\scop_cells",
                         jxyz_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\jxyz")

        bound_tool = Bonds_tool(r"E:\vscode\Muon_Imaging_Algorithm\data\no_name\17_58MaMian_bnd.den")
        min_max_data = bound_tool.get_bonds_min_max()

        refs_tool = Ref_tools(r"E:\vscode\Muon_Imaging_Algorithm\data\no_name\17_58MaMian_ref.den")
        refs = refs_tool.get_data()

        for j in cells5:
            j -= 1
            refs[j] = 1.89
            min_max_data[j] = [1.86, 1.92]
        for j in cells:
            j -= 1
            refs[j] = 1.7
            min_max_data[j] = [1.68, 1.73]
        for j in cells4:
            j -= 1
            refs[j] = 1.89
            min_max_data[j] = [0, 2]
        for j in tqdm.trange(mesh.cells_count()):
            j -= 1
            if refs[j] > 2.1:
                refs[j] = 1.7
                min_max_data[j] = [1.68, 1.73]

        bound_tool.update_data(min_max_data)
        refs_tool.update_data(refs)
    elif model == 4:
        # 根据射线编号从0开始显示误差最大的100条射线---
        tool = show_ray_trace_tools([131, 92, 55])
        file = open(r"E:\vscode\Muon_Imaging_Algorithm\data\output\rel\d_obs_pred_diff")
        ray_id_errs = file.readlines()
        ids = []
        for i in range(99, -1, -1):
            ids.append(int(ray_id_errs[i].split()[0]))
        obs_tool = obs_tools(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\real_data\27_obs.dat")
        rays_infos = obs_tool.get_data()
        file = open(r"E:\vscode\Muon_Imaging_Algorithm\data\output\rel\d_big_diff_obs", 'w')
        middle = []
        for id in ids:
            ray_info = rays_infos[id]
            # if int(ray_info[0]) == 5 : #and abs(float(ray_info[5])) < 0.1:
            middle.append(id)
            for info in ray_info:
                file.write(str(info))
                file.write(" ")
            file.write('\n')
        file.close()
        tool.mark_ray(ids=middle, file_ij=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\Gij",
                      file_res=r"E:\vscode\Muon_Imaging_Algorithm\data\output\rel\ray_way")
    elif model == 5:
        # 标记夯土,再反推出城墙墙砖--2022-8-4
        air_j = Air_j().get_air_j_from_file(r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\unneed_j copy")
        mesh = MeshTools(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\real_data\17_58MaMian.msh")
        tool = show_ray_trace_tools(mesh.get_shape())
        point1, point2, point3, point4, point5, point6, point7, point8 = [-63, -15.227, 10.325], \
                                                                         [86.678, -12.393, 10.813], \
                                                                         [88.347, -11.177, 0.53938], \
                                                                         [-64.848, -13.837, -0.2743], \
                                                                         [-66.389, -30.252, 10.108], \
                                                                         [91.968, -26.097, 9.8065], \
                                                                         [-68.831, -32.64, 0.2844], \
                                                                         [93.546, -28.995, -0.0243]
        corners = [point6, point5, point7, point8, point2, point1, point4, point3]
        # show_points_3D(corners)
        wall_cells = tool.get_j_from_scope_mesh(corners, mesh, show=False)
        # tool.draw_points(wall_cells, res_file=r"E:\vscode\Muon_Imaging_Algorithm\data\output\rel\scop_cells",
        #                  jxyz_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\jxyz")
        myPrint.myPrint_Hint(wall_cells)

        point1, point2, point3, point4, point5, point6, point7, point8 = [2.7311, -2.7171, 10.66], \
                                                                         [20.658, -2.3212, 10.768], \
                                                                         [21.87, -1.4941, 0.1683], \
                                                                         [1.3744, -1.7812, 0.010574], \
                                                                         [2.9297, -14.979, 10.538], \
                                                                         [20.874, -14.639, 10.599], \
                                                                         [22.025, -13.343, 0.30293], \
                                                                         [1.528, -13.684, -0.0783]
        corners = [point6, point5, point8, point7, point2, point1, point4, point3]
        # show_points_3D(corners)
        mamian_cells = tool.get_j_from_scope_mesh(corners, mesh, show=False)
        # tool.draw_points(mamian_cells, res_file=r"E:\vscode\Muon_Imaging_Algorithm\data\output\rel\scop_cells",
        #                  jxyz_file=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\jxyz")
        myPrint.myPrint_Hint(mamian_cells)

        #####小房子
        x_range = [2.454, 4.2]
        y_range = [-11.48, -13.2]
        z_range = [12.4, 10.58]
        corners = [[max(x_range), min(y_range), max(z_range)], [min(x_range), min(y_range), max(z_range)],
                   [min(x_range), min(y_range), min(z_range)],
                   [max(x_range), min(y_range), min(z_range)],
                   [max(x_range), max(y_range), max(z_range)], [min(x_range), max(y_range), max(z_range)],
                   [min(x_range), max(y_range), min(z_range)],
                   [max(x_range), max(y_range), min(z_range)]
                   ]
        # show_points_3D(data=corners)
        room_cells = tool.get_j_from_scope_mesh(corners, mesh, show=False)
        
        print(room_cells)
        refs_file = open(r"E:\vscode\Muon_Imaging_Algorithm\data\new_methods_ref.den", "w")
        bounds_file = open(r"E:\vscode\Muon_Imaging_Algorithm\data\new_methods_bounds.den", 'w')
        x_cells_len = mesh.get_xs()
        x_cells_len_min = min(x_cells_len)
        y_cells_len = mesh.get_ys()
        y_cells_len_min = min(y_cells_len)
        z_cells_len = mesh.get_zs()
        z_cells_len_min = min(z_cells_len)
        air_j = set(air_j)
        for i in range(mesh.cells_count()):
            if i + 1 in air_j:
                refs_file.write("0")
                refs_file.write('\n')
                bounds_file.write('-0.01')
                bounds_file.write(' ')
                bounds_file.write('0.02')
                bounds_file.write('\n')
            elif i in mamian_cells:
                refs_file.write("1.89")
                refs_file.write('\n')
                bounds_file.write('0')
                bounds_file.write(' ')
                bounds_file.write('2.0')
                bounds_file.write('\n')
            elif i in wall_cells:
                refs_file.write("1.89")
                refs_file.write('\n')
                bounds_file.write('1.87')
                bounds_file.write(' ')
                bounds_file.write('1.91')
                bounds_file.write('\n')
            elif i in room_cells:
                refs_file.write("0")
                refs_file.write('\n')
                bounds_file.write('-0.01')
                bounds_file.write(' ')
                bounds_file.write('0.01')
                bounds_file.write('\n')
            else:
                # 墙皮也分为马面墙皮和城墙墙皮
                x, y, z = getxyz_from_shape(mesh.get_shape(), i + 1)
                if x_cells_len[x - 1] < 0.46 and y_cells_len[y - 1] < 0.46 and z_cells_len[
                    z - 1] <0.3:
                    # 马面墙皮
                    refs_file.write("1.95")
                    refs_file.write('\n')
                    bounds_file.write('1.9')
                    bounds_file.write(' ')
                    bounds_file.write('2.0')
                    bounds_file.write('\n')
                else:
                    # 城墙墙皮
                    refs_file.write("1.7")
                    refs_file.write('\n')
                    bounds_file.write('1.69')
                    bounds_file.write(' ')
                    bounds_file.write('1.71')
                    bounds_file.write('\n')
        refs_file.close()
        bounds_file.close()
    elif model == 6:
        """删除某条探测器数据观察结果情况"""
        obs_tool1 = obs_tools(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\real_data\27_obs.dat")
        all_rays = obs_tool1.get_data()
        rays_set = set()
        count = 0
        flag=False
        for i in range(len(all_rays)):
            middle = ""
            for datum in all_rays[i]:
                middle += str(datum)+" "
            if  int(all_rays[i][0]) == 6:
                if flag:
                    rays_set.add(middle)
                else:
                    flag=True
        file_obj = open(r"E:\vscode\Muon_Imaging_Algorithm\data\no_name\obs_del_detector_6.dat", 'w')
        file_obj.write(str(len(all_rays)-len(rays_set)))
        file_obj.write('\n')
        for all_ray in all_rays:
            middle = ""
            for datum in all_ray:
                middle += str(datum) + " "
            if middle in rays_set:
                continue
            file_obj.write(middle)
            file_obj.write('\n')
            count+=1
        file_obj.close()
        print(count)
        print(len(rays_set))
    elif model==7:
        #城墙底面格子物理坐标生成器
        mesh_tool=MeshTools(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\real_data\17_58MaMian.msh")
        res_file=open(r"E:\vscode\Muon_Imaging_Algorithm\data\Looking_for_area_through_the_ray_intersection\points_bottom.den",'w')
        xyz=mesh_tool.get_shape()
        z = xyz[2]
        j_s=set()
        for x in range(xyz[0]):
            for y in range(xyz[1]):
                # j_s.add(getj_from_xyz(shape=mesh_tool.get_shape(),xyz=(x+1,y+1,1)))
                p_xyz=mesh_tool.get_coordinates_form_xyz(x,y,1)
                res_file.write(str(p_xyz[0]))
                res_file.write(' ')
                res_file.write(str(p_xyz[1]))
                res_file.write(' ')
                res_file.write(str(p_xyz[2]))
                res_file.write('\n')
        # for j in range(mesh_tool.cells_count()):
        #     j=j+1
        #     if j in j_s:
        #         res_file.write('1')
        #         res_file.write('\n')
        #     else:
        #         res_file.write('-1')
        #         res_file.write("\n")
        res_file.close()
    elif model==8:
        refs_tool=Ref_tools(r"E:\vscode\Muon_Imaging_Algorithm\data\58MaMian_ref.den")
        bound_tool=Bonds_tool(r"E:\vscode\Muon_Imaging_Algorithm\data\58MaMian_bnd.den")
        data=refs_tool.get_data()
        bound_data=bound_tool.get_bonds_min_max()
        room_js=[268294, 268295, 268296, 268297, 268298, 268299, 268300, 268301, 268302, 239639, 239640, 239641, 239642, 239643, 239644, 239645, 239646, 239647, 253994, 253995, 253996, 253997, 253998, 253999, 254000, 254001, 254002, 268349, 268350, 268351, 268352, 268353, 268354, 268355, 268356, 268357, 239694, 239695, 239696, 239697, 239698, 239699, 239700, 239701, 239702, 254049, 254050, 254051, 254052, 254053, 254054, 254055, 254056, 254057, 268404, 268405, 268406, 268407, 268408, 268409, 268410, 268411, 268412, 254104, 254105, 254106, 254107, 254108, 254109, 254110, 254111, 254112, 268459, 268460, 268461, 268462, 268463, 268464, 268465, 268466, 268467, 268514, 268515, 268516, 268517, 268518, 268519, 268520, 268521, 268522, 232159, 232160, 232161, 232162, 232163, 232164, 232165, 232166, 232167, 232214, 232215, 232216, 232217, 232218, 232219, 232220, 232221, 232222, 246569, 246570, 246571, 246572, 246573, 246574, 246575, 246576, 246577, 232269, 232270, 232271, 232272, 232273, 232274, 232275, 232276, 232277, 246624, 246625, 246626, 246627, 246628, 246629, 246630, 246631, 246632, 260979, 260980, 260981, 260982, 260983, 260984, 260985, 260986, 260987, 232324, 232325, 232326, 232327, 232328, 232329, 232330, 232331, 232332, 246679, 246680, 246681, 246682, 246683, 246684, 246685, 246686, 246687, 261034, 261035, 261036, 261037, 261038, 261039, 261040, 261041, 261042, 232379, 232380, 232381, 232382, 232383, 232384, 232385, 232386, 232387, 275389, 275390, 275391, 275392, 275393, 275394, 275395, 275396, 275397, 246734, 246735, 246736, 246737, 246738, 246739, 246740, 246741, 246742, 261089, 261090, 261091, 261092, 261093, 261094, 261095, 261096, 261097, 232434, 232435, 232436, 232437, 232438, 232439, 232440, 232441, 232442, 275444, 275445, 275446, 275447, 275448, 275449, 275450, 275451, 275452, 246789, 246790, 246791, 246792, 246793, 246794, 246795, 246796, 246797, 261144, 261145, 261146, 261147, 261148, 261149, 261150, 261151, 261152, 232489, 232490, 232491, 232492, 232493, 232494, 232495, 232496, 232497, 275499, 275500, 275501, 275502, 275503, 275504, 275505, 275506, 275507, 246844, 246845, 246846, 246847, 246848, 246849, 246850, 246851, 246852, 261199, 261200, 261201, 261202, 261203, 261204, 261205, 261206, 261207, 275554, 275555, 275556, 275557, 275558, 275559, 275560, 275561, 275562, 246899, 246900, 246901, 246902, 246903, 246904, 246905, 246906, 246907, 261254, 261255, 261256, 261257, 261258, 261259, 261260, 261261, 261262, 275609, 275610, 275611, 275612, 275613, 275614, 275615, 275616, 275617, 261309, 261310, 261311, 261312, 261313, 261314, 261315, 261316, 261317, 275664, 275665, 275666, 275667, 275668, 275669, 275670, 275671, 275672, 275719, 275720, 275721, 275722, 275723, 275724, 275725, 275726, 275727, 239364, 239365, 239366, 239367, 239368, 239369, 239370, 239371, 239372, 239419, 239420, 239421, 239422, 239423, 239424, 239425, 239426, 239427, 253774, 253775, 253776, 253777, 253778, 253779, 253780, 253781, 253782, 239474, 239475, 239476, 239477, 239478, 239479, 239480, 239481, 239482, 253829, 253830, 253831, 253832, 253833, 253834, 253835, 253836, 253837, 268184, 268185, 268186, 268187, 268188, 268189, 268190, 268191, 268192, 239529, 239530, 239531, 239532, 239533, 239534, 239535, 239536, 239537, 253884, 253885, 253886, 253887, 253888, 253889, 253890, 253891, 253892, 268239, 268240, 268241, 268242, 268243, 268244, 268245, 268246, 268247, 239584, 239585, 239586, 239587, 239588, 239589, 239590, 239591, 239592, 253939, 253940, 253941, 253942, 253943, 253944, 253945, 253946, 253947]
        for i in room_js:
            data[i-1]=0
            bound_data[i-1]=[-0.001,0.001]
        refs_tool.update_data(data)
        bound_tool.update_data(bound_data)
        
        
