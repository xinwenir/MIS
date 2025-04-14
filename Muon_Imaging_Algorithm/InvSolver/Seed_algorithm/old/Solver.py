# author:高金磊
# datetime:2022/9/1 10:05
import dis
from operator import ne
import time
from copy import copy

from scipy.sparse import csc_matrix

from InvDataTools.Gij_tools import G_Data
from InvDataTools.Jxyz_Tools import getj_from_xyz, getxyz_from_shape
from InvDataTools.d_tools import d_tools
from InvDataTools.Bonds_tools import Bonds_tool
from InvSolver.Seed_algorithm.Tools import data_tools
# from InvSysTools.MyTools.Monitor.Process_monitor import process_monitoring
# process_monitoring(interval=5)
from InvSolver.Seed_algorithm.Setting import *
def run_seed():
    start_time=time.time()
    data_tool = data_tools(start_data_file,mesh_file)
    refs_x0=copy(data_tool.get_data())
    data_tool_step = None
    bounds_tool =Bonds_tool(boods_file=bounds_file)
    bounds=bounds_tool.get_bonds_min_max()
    
    moudles=[]
    for s in open(moudles_file,'r').readlines():
        if s!="":
            moudles.append(float(s))
    
    # 加载稀疏矩阵
    g_tool = G_Data(matrix_file%("G"),
                    matrix_file%("Gij"))
    ijs = g_tool.get_ij()
    ij_value = g_tool.get_GV()
    G = csc_matrix((ij_value, ([i - 1 for i in ijs[0]], [i - 1 for i in ijs[1]])),
                   (max(ijs[0]), data_tool.mesh.cells_count()), dtype=float)
    d_tool = d_tools(matrix_file%("d"))
    d=d_tool.get_d()
    d_err=d_tool.get_d_err()
    seed_js = set()
    # seed_js.add(start_seed)#CRB5_中心值
    
    #!取出矿中的一些点作为初始值
    data_tool_middle = data_tools(moudles_file,mesh_file)
    middle_5=data_tool_middle.get_data()
    Interval_of_sampling=50
    interval_count=0
    for i in range(len(middle_5)):
        if middle_5[i]==0.2:
            interval_count+=1
            if interval_count%Interval_of_sampling==0:
                seed_js.add(i+1)
    """"""""""""""""""""""""""""""""""""

    old_seeds = set()
    from InvSolver.Seed_algorithm.objective_function import obj_fun, Obj_fun_Tools
    obj_fun_tool=Obj_fun_Tools(G,d,d_err)

    obj_fun_value = obj_fun(G, data_tool.get_data(), d,d_err)
    count = 0

    middle = 0
    cold_map={}
    if model == 1:
        # 一个一个下降---禁止种子不倍增
        print("测量误差:",obj_fun(G, moudles, d,d_err)) 
        while count<max_count:
            print("误差:",data_tool.calculation_results_gap(moudles))
            print(len(seed_js),len(old_seeds), obj_fun_value,"下降:",obj_fun_value-middle,"已经用时:",time.time()-start_time)
            middle=obj_fun_value
            # 数据备份
            # old_data=copy(data_tool.get_data())
            # 获取邻居
            if count==0:
                new_seeds=seed_js
            else:
                new_seeds = set()
                for seed_j in seed_js:
                    new_seeds = new_seeds.union(data_tool.get_neighbor_js(seed_j))
                seed_js.clear()
            for new_seed in new_seeds:
                if new_seed in old_seeds:
                    continue
                old_value =data_tool.get_value(new_seed)
                # obj_fun_tool.get_j_optimal_value(data_tool,new_seed,bounds[new_seed-1],max_iter=5,accuracy=0.2)
                obj_fun_tool.get_j_line_search_optimal_value(data_tool,new_seed,bounds[new_seed-1],issmooth=False,distance=punishment_Search_distancevalue_multiple*(count)**2,start_seed=start_seed,old_seeds=old_seeds)
                # obj_fun_tool.mult_obj_solver(data_tool,new_seed,bounds,issmooth=len(old_seeds)>30)
                if  abs(refs_x0[new_seed-1]-data_tool.get_value(new_seed))>=seed_Minimum_density_variation:#简单判断算法
                # 接受
                    seed_js.add(new_seed)
                    old_seeds.add(new_seed)
                else:  # 回滚并抛弃
                    old_seeds.add(new_seed)

                # 总目标函数下降
                # fun_value = obj_fun(G, data_tool.get_data(), d,d_err)
                # if fun_value < obj_fun_value:
                #     # 接受
                #     obj_fun_value = fun_value
                #     seed_js.add(new_seed)
                #     old_seeds.add(new_seed)
                # else:  # 回滚并抛弃
                #     old_seeds.add(new_seed)
                #     data_tool.alt_data(new_seed, old_value)
            #小循环纠正密度
            # for j in old_seeds:
            #     obj_fun_tool.get_j_line_search_optimal_value(data_tool,j,bounds[j-1],issmooth=True,distance=0,start_seed=start_seed)
                
                # obj_fun_tool.mult_obj_solver(data_tool,j,bounds,issmooth=(len(old_seeds)>10))
            # cold=2
            # cold_count=len(old_seeds)
            # for old_seed in old_seeds:
            #     seed_cold=cold_map.get(old_seed,0)
            #     if seed_cold==0:
            #         old_j_value=data_tool.get_value(old_seed)
            #         # new_j_value=obj_fun_tool.get_j_optimal_value(data_tool, old_seed,bounds[old_seed-1] ,max_iter=5,accuracy=0.1)
            #         new_j_value=obj_fun_tool.get_j_line_search_optimal_value(data_tool,old_seed,bounds[old_seed-1])
            #         cold_count-=1
            #         if abs(old_j_value-new_j_value) <5*0.1:#单边下降没有往复
            #             seed_cold+=1
            #     elif seed_cold>=cold:
            #         seed_cold=-seed_cold-cold
            #         seed_cold+=1
            #     else:
            #         seed_cold+=1
            #         # pass
            #     cold_map[old_seed]=seed_cold
            # print("当前在冷却:",cold_count)


            """""""不使用小循环"""""""
            # for old_seed in old_seeds:
            #     obj_fun_tool.get_j_optimal_value(data_tool, old_seed, (0, 2.7),max_iter=8,accuracy=0.1)
            # old_seeds.clear()

            obj_fun_value = obj_fun(G, data_tool.get_data(), d,d_err)

            if len(seed_js) == 0:
                break
            count += 1
            if count % 1 == 0:
                # if count %3==0:
                #     for j in old_seeds:
                #         obj_fun_tool.get_j_line_search_optimal_value(data_tool,j,bounds[j-1],issmooth=True,distance=0.0,start_seed=start_seed,old_seeds=old_seeds)
                    # obj_fun_tool.mult_obj_solver(data_tool,j,bounds,issmooth=(len(old_seeds)>10))
                # 快照
                print("快照", count)
                data_tool.output_res(out_put_dir+"\seed_res%s" % (str(count)))
        old_seeds_list=[*old_seeds]
        for i in range(10):
            print("误差:",data_tool.calculation_results_gap(moudles),"  misfit:",obj_fun(G, data_tool.get_data(), d,d_err))
            old_seeds_list.reverse()
            for j in old_seeds_list:
                obj_fun_tool.get_j_line_search_optimal_value(data_tool,j,bounds[j-1],issmooth=i>8,distance=0.0,start_seed=start_seed,old_seeds=old_seeds)
                # obj_fun_tool.mult_obj_solver(data_tool,j,bounds,issmooth=(len(old_seeds)>10))

    elif model == 2:
        while 1:
            print(len(seed_js),len(old_seeds), obj_fun_value,"下降:",obj_fun_value-middle,"已经用时:",time.time()-start_time)

            # 数据备份
            # old_data=copy(data_tool.get_data())
            # 获取邻居
            new_seeds = []
            for seed_j in seed_js:
                new_seeds.append(data_tool.get_neighbor_js(seed_j))
            seed_js.clear()
            for new_seed_set in new_seeds:
                old_seed_value = {}
                sum=0
                for new_seed in new_seed_set:
                    if new_seed in old_seeds:
                        continue
                    old_value = data_tool.get_value(new_seed)
                    obj_fun_tool.get_j_optimal_value(data_tool, new_seed, bounds[new_seed-1], max_iter=5, accuracy=0.2)
                    # old_value = data_tool.alt_data(new_seed, 0.4)
                    # data_tool_step.alt_data(new_seed, -count)
                    old_seed_value[new_seed] = old_value
                    if abs(old_value-data_tool.get_value(new_seed))>=0.6:
                        sum+=1
                new_fun_value = obj_fun(G, data_tool.get_data(), d)
                if new_fun_value < obj_fun_value and sum>=2:
                    obj_fun_value = new_fun_value
                    for new_seed in new_seed_set:
                        seed_js.add(new_seed)
                        old_seeds.add(new_seed)
                # else:
                #     for new_seed in new_seed_set:
                #         if new_seed in old_seeds:
                #             continue
                        # data_tool.alt_data(new_seed, old_seed_value[new_seed])
            count += 1

            cold=5
            cold_count=len(old_seeds)
            for old_seed in old_seeds:
                seed_cold=cold_map.get(old_seed,0)
                if seed_cold==0:
                    old_j_value=data_tool.get_value(old_seed)
                    new_j_value=obj_fun_tool.get_j_optimal_value(data_tool, old_seed,bounds[old_seed-1] ,max_iter=5,accuracy=0.1)
                    cold_count-=1
                    if abs(old_j_value-new_j_value) <5*0.1:#单边下降没有往复
                        seed_cold+=1
                elif seed_cold>=cold:
                    seed_cold=-seed_cold-cold
                    seed_cold+=1
                else:
                    seed_cold+=1
                    # pass
                cold_map[old_seed]=seed_cold
            print("当前在冷却:",cold_count)


            """""""不使用小循环"""""""
            # for old_seed in old_seeds:
            #     obj_fun_tool.get_j_optimal_value(data_tool, old_seed, (0, 2.7),max_iter=8,accuracy=0.1)
            # old_seeds.clear()

            obj_fun_value = obj_fun(G, data_tool.get_data(), d)

            if len(seed_js) == 0:
                break
            count += 1
            if count % 3 == 0:
                # 快照
                print("快照", count)
                data_tool.output_res(out_put_dir+"\seed_res%s" % (str(count)))
    elif model ==3 :
        #理解错误,误认为多个种子增长一个种子
        iteration_count=0#外层迭代次数
        while iteration_count<5:
            new_seed = 268408
            neighbors=set()
            old_seeds.clear()
            middle=obj_fun_value
            print(len(seed_js),len(old_seeds), obj_fun_value,"下降:",obj_fun_value-middle,"已经用时:",time.time()-start_time)
            while True:
                print(middle)
                #更新当前体素的值
                obj_fun_tool.get_j_optimal_value(data_tool,new_seed,(0,3))
                new_neighbors= data_tool.get_neighbor_js(new_seed)
                neighbors=neighbors.union(new_neighbors)
                neighbors=neighbors.difference(old_seeds)
                old_seeds.add(new_seed)
                middle = obj_fun(G, data_tool.get_data(), d)
                #找出邻居中下降效果最佳的格子
                info_data=[]
                for j in neighbors:
                    
                    old_j_value=obj_fun_tool.get_j_optimal_value(data_tool,j,(0,3))
                    obj_fun_value_j=obj_fun(G, data_tool.get_data(), d)
                    new_j_value=data_tool.alt_data(j,old_j_value)
                    if obj_fun_value<=obj_fun_value_j and abs(old_j_value-new_j_value)<0.4:
                        #对目标函数下降丝毫没有贡献
                        old_seeds.add(j)
                        continue
                    
                    info_data.append([j,obj_fun_value_j,old_j_value,new_j_value])
                #选择最佳的邻居格子作为下次的种子
                min_i=0
                min_value=obj_fun_value
                for i in range(len(info_data)):
                    if min_value>info_data[i][1]:
                        min_i=i
                        min_value=info_data[i][1]
                if min_value==obj_fun_value:
                    break
                new_seed=info_data[min_i][0]
                middle=info_data[min_i][1]
                data_tool.alt_data(new_seed,info_data[min_i][3])
                # if len(seed_js) == 0:
                #     break
                
                
                count += 1
                if count % 3 == 0:
                    # 快照
                    print("快照", count)
                    data_tool.output_res(out_put_dir+"\seed_res%s_%s" % (str(iteration_count),str(count)))
            iteration_count+=1
    elif model == 4:
        
        new_seeds=set()
        while True:
            seed_js=seed_js.union(new_seeds)
            seed_js=seed_js.difference(old_seeds)
            new_seeds.clear()
            print(len(seed_js),"  ",middle)
            #遍历所有种子
            for seed_j in seed_js:
                #获取当前种子的邻居
                neighbors=data_tool.get_neighbor_js(seed_j)
                #调整当前种子的密度
                obj_fun_tool.get_j_optimal_value(data_tool,seed_j,bounds[seed_j-1])
                #更新目标函数
                middle = obj_fun(G, data_tool.get_data(), d)
                #遍历其邻居获取新的种子
                info_data=[]
                for j in neighbors:
                    if j in old_seeds or j in seed_js or j in new_seeds:
                        continue
                    old_j_value=data_tool.get_value(j)
                    obj_fun_tool.get_j_optimal_value(data_tool,j,bounds[j-1])
                    obj_fun_value_j=obj_fun(G, data_tool.get_data(), d)
                    new_j_value=data_tool.alt_data(j,old_j_value)
                    if middle<=obj_fun_value_j and abs(old_j_value-new_j_value)<0.4:
                        #对目标函数下降丝毫没有贡献
                        old_seeds.add(j)
                        continue
                        
                    info_data.append([j,obj_fun_value_j,old_j_value,new_j_value])
                #选择最佳的邻居格子作为下次的种子
                min_i=0
                min_value=middle
                for i in range(len(info_data)):
                    if min_value>info_data[i][1]:
                        min_i=i
                        min_value=info_data[i][1]
                if min_value==middle:
                    old_seeds.add(seed_j)
                    break
                new_seed=info_data[min_i][0]
                middle=info_data[min_i][1]
                data_tool.alt_data(new_seed,info_data[min_i][3])
                new_seeds.add(new_seed)
            count += 1
            if count % 5 == 0:
                # 快照
                print("快照", count)
                data_tool.output_res(out_put_dir+"\seed_res%s_%s" % (str(0),str(count)))
    data_tool.output_res(out_put_dir+r"\wall")
    
    # center_of_inertia=[0,0,0]
    weight=0
    center_x=0
    center_y=0
    center_z=0
    #计算中心点
    for i in old_seeds:
        xyz=getxyz_from_shape(data_tool.mesh.get_shape(),i)
        cell_weight=abs(data_tool.get_refs_value(i)-data_tool.get_value(i))
        weight+=cell_weight
        center_x+=xyz[0]*cell_weight
        center_y+=xyz[1]*cell_weight
        center_z+=xyz[2]*cell_weight
    center_x/=weight
    center_y/=weight
    center_z/=weight
    center_xyz=data_tool.mesh.get_coordinates_form_xyz(int(center_x),int(center_y),int(center_z))
    print("建议起始点的中点坐标为:",center_xyz)
    print("建议起始点的编号:",getj_from_xyz(data_tool.mesh.get_shape(),(int(center_x),int(center_y),int(center_z))))
    
    data_tool.output_smooth_res(output_file=out_put_dir+"\wall_smooth")
    print("平滑后的误差:",data_tool.calculation_results_gap(moudles))
    # data_tool_step.output_res(r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\result\CBC1%s" % ("step"))


if __name__ == '__main__':
    run_seed()