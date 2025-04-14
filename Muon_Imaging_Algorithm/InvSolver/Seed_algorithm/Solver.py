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
from InvSolver.Seed_algorithm.BCD_Opt_Tools import BCD_Opt_Tools

from InvSolver.Seed_algorithm.objective_function import obj_fun
from InvSolver.Seed_algorithm.Object_fun import Obj_fun_Tools

class myresult():
    def __init__(self,start_misfit,new_misfit,pre_seed_js,determine_seed_js,opt_result,new_seed_js) -> None:
        self.new_misfit=new_misfit
        self.pre_seed_js=pre_seed_js
        self.determine_seed_js=determine_seed_js
        self.opt_result=opt_result
        self.new_seed_js=new_seed_js
        self.start_misfit=start_misfit
        

class Solver:
    def __init__(self) -> None:
        from InvSolver.Seed_algorithm.Setting import out_put_dir,max_count,start_data_file,mesh_file,matrix_file,start_seed,threshold,seed_Minimum_density_variation,bounds_file,moudles_file,refs_file
        seed_js_history=[]#记录每一次seed的组构成
        data_tool = data_tools(start_data_file,mesh_file)#用来管理结果
        refs=data_tools(refs_file,mesh_file).get_data()#最初状态_或者说是参考值
        bounds_tool =Bonds_tool(boods_file=bounds_file)
        bounds=bounds_tool.get_bonds_min_max()
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
        
        seed_js = set(start_seed)#当前的种子值
        # seed_js.add(47672)
        # seed_js.add(88888)
        # seed_js.add(130104)
        #wall
        # seed_js.add(556462)
        # seed_js.add(556752)
        # seed_js.add(556427)
        
        #zzg
        # seed_js.add(118707)


        #取出矿中的一些点作为初始值
        # data_tool_middle = data_tools(r"E:\vscode\Muon_Imaging_Algorithm\data\zzg\input\20_muto3d.den",mesh_file)
        # middle_5=data_tool_middle.get_data()
        # Interval_of_sampling=40
        # interval_count=0
        # for i in range(len(middle_5)):
        #     if middle_5[i]>0.2:
        #         interval_count+=1
        #         if interval_count%Interval_of_sampling==0:
        #             seed_js.add(i+1)
        
        # print(len(seed_js))
        
        
        
        all_js=set()#记录所有参与者(种子和pre种子)
        self.old_js=set()#记录所有过去的种子
        # seed_js_history.append(seed_js)
        self.seed_js_history=seed_js_history
        self.data_tool=data_tool
        self.seed_js=seed_js
        self.all_js=all_js
        self.refs=refs
        self.bounds=bounds
        self.G=G
        self.d=d
        self.derr=d_err
        self.threshold=threshold
        self.seed_Minimum_density_variation=seed_Minimum_density_variation
        self.mesh_file=mesh_file
        self.moudles_file=moudles_file
        self.max_count=max_count
        self.out_put_dir=out_put_dir
        
        self.Ancestors_seed=[0 for i in range(len(data_tool.get_data()))]
        


        
    def BCD_single_iteration(self,data_tool:data_tools=None,seed_js:set=None,moudles:list=None,threshold=None,distance=0,issmooth=True,find_determine_seeds=True):
        if data_tool==None:
            data_tool=self.data_tool
        if seed_js==None:
            seed_js=self.seed_js
        if threshold ==None:
            from InvSolver.Seed_algorithm.Setting import threshold
            threshold=threshold
        if moudles!=None:
            print("误差:",data_tool.calculation_results_gap(moudles))
        start_misfit=obj_fun(self.G,data_tool.get_data(),self.d,self.derr)
        pre_seed_js=seed_js#为了增加复用性,这里不去掉seed中已经出现过的seed(传什么就计算什么),只在最后返回中去除.
        # for seed in seed_js:
        #     if seed not in self.all_js:
        #         pre_seed_js.append(seed)
        
        obj_fun_tool=Obj_fun_Tools()
        #将pre_seed_js作为一块进行优化
        #构造小的优化问题
        bcd_opt_tools=BCD_Opt_Tools(self.G,self.d,self.derr,self.refs,self.bounds)
        new_G,new_x0,new_d,new_derrs,new_refs,new_bounds,seeds_neighbor_js=bcd_opt_tools.init_BCD_opt_question(data_tool,pre_seed_js,threshold=threshold)
        ancestors_seeds_js=[]
        for j in seed_js:
            if self.Ancestors_seed[j-1]==0:
                self.Ancestors_seed[j-1]=j
            ancestors_seeds_js.append(self.Ancestors_seed[j-1])
        #对该优化问题进性优化
        result=bcd_opt_tools.optimize_single_Block(obj_fun_tool,seed_js, new_G,new_x0,new_d,new_derrs,new_refs,new_bounds,distance,ancestors_seeds_js,issmooth,seeds_neighbor_js,data_tool)
        opt_x=result.x.tolist()
        #将优化后的值写入结果
        determine_seed_js=set()#记录那些种子被确定为合格的种子
        for j in seed_js:
            self.all_js.add(j)
            middle_rho=data_tool.alt_data(j,opt_x.pop(0)[0])#!更改没有验证
            if abs(data_tool.get_refs_value(j)-data_tool.get_value(j))>=self.seed_Minimum_density_variation:
                self.old_js.add(j)
                determine_seed_js.add(j)
            else:
                data_tool.alt_data(j, middle_rho)
                pass
        new_misfit=obj_fun(self.G,data_tool.get_data(),self.d,self.derr)
        
        #获取确定种子附近的seed
        new_seed_js=set()
        if find_determine_seeds:
            for j in determine_seed_js:
                neighbor_js=data_tool.get_neighbor_js(j)
                parent_seed=self.Ancestors_seed[j-1]
                if parent_seed==0:
                    self.Ancestors_seed[j-1]=j
                    parent_seed=j
                for neighbor_j in neighbor_js:
                    if neighbor_j not in self.all_js:
                        new_seed_js.add(neighbor_j) 
                        #确定该种子的祖先节点
                        self.Ancestors_seed[neighbor_j-1]=parent_seed
        
                
        return myresult(start_misfit,new_misfit,pre_seed_js,determine_seed_js,result,new_seed_js)
        
    def BCD_solver(self,seed_js=None):
        
        #!计算时间
        start_time=time.time()
        complete_seed_js=set()
        if seed_js is None:
            seed_js=self.seed_js
        #获取moudles数据
        if self.moudles_file!=None:
            moudles=data_tools(self.moudles_file,self.mesh_file).get_data()
        else:
            moudles=None
        count=0
        while count<self.max_count and len(self.seed_js)>0:
            if moudles is not None:
                print("误差:",self.data_tool.calculation_results_gap(moudles))
            
            count+=1
            #执行当前块的优化算法
            result=self.BCD_single_iteration(distance=count-1,issmooth=False)
            all_time=int(time.time()-start_time)
            print("第%s次初始misfit:%s 优化后的misfit:%s 下降%s 本次存活的种子数量:%s pre新种子数量:%s 现有种子总数:%s 当前用时%s分%s秒"%(
                str(count),
                str(result.start_misfit),str(result.new_misfit),str(result.start_misfit-result.new_misfit),str(len(result.determine_seed_js)),str(len(result.new_seed_js)),
                str(len(self.old_js)),
                str(int(all_time/60)),
                str(all_time%60)
            ))
            self.alt_seed_js(result.new_seed_js)
            
            
            # if count %10==0:
            #     print("执行小循环")
            #     self.threshold*=2
            #     for i in range(len(self.seed_js_history)):
                    
            #         if i in complete_seed_js:#之前优化没有效果,认为已经优化完成
            #             continue
            #         history_js=self.seed_js_history[i]
            #          #执行当前块的优化算法
            #         middle_result=self.BCD_single_iteration(distance=0,seed_js=history_js,issmooth=True,threshold=self.threshold)
            #         all_time=int(time.time()-start_time)
            #         print("初始misfit:%s 优化后的misfit:%s 下降%s 本次存活的种子数量:%s pre新种子数量:%s 现有种子总数:%s 当前用时%s分%s秒"%(
            #             str(middle_result.start_misfit),str(middle_result.new_misfit),str(middle_result.start_misfit-middle_result.new_misfit),str(len(middle_result.determine_seed_js)),str(len(middle_result.new_seed_js)),
            #             str(len(self.old_js)),
            #             str(int(all_time/60)),
            #             str(all_time%60)
            #         ))
            #         self.data_tool.output_res(output_file=self.out_put_dir+"\seed_res%s_%s"%(str(count),str(i)))
            #         if abs(middle_result.start_misfit-middle_result.new_misfit)<0.1 and len(history_js)>40:
            #             complete_seed_js.add(i)
            #     self.threshold/=2
            output_file=self.out_put_dir+"\seed_res%s"%(str(count))
            self.data_tool.output_res(output_file)
            print("当前结果路径:",output_file)
            if len(result.determine_seed_js)==0:
                break
            self.seed_js_history.append(result.determine_seed_js)
        middle_count=0
        history=copy(self.seed_js_history)
        for i in range(2):
            print("做最后的平滑处理")
            self.threshold=0.4
            # history.reverse()
            for history_js in history:
                if len(history_js)==0:
                    continue
                #执行当前块的优化算法
                middle_result=self.BCD_single_iteration(distance=0,seed_js=history_js,issmooth=True,threshold=self.threshold,find_determine_seeds=False)
                all_time=int(time.time()-start_time)
                print("初始misfit:%s 优化后的misfit:%s 下降%s 本次存活的种子数量:%s pre新种子数量:%s 现有种子总数:%s 当前用时%s分%s秒"%(
                    str(middle_result.start_misfit),str(middle_result.new_misfit),str(middle_result.start_misfit-middle_result.new_misfit),str(len(middle_result.determine_seed_js)),str(len(middle_result.new_seed_js)),
                    str(len(self.old_js)),
                    str(int(all_time/60)),
                    str(all_time%60)
                ))
                output_file=self.out_put_dir+"\seed_res_%s_%s"%("smmoth",str(middle_count))
                self.data_tool.output_res(output_file)
                print("当前结果路径:",output_file)
                middle_count+=1
                if moudles is not None:
                    print("误差:",self.data_tool.calculation_results_gap(moudles))
            self.threshold=0.4
        if moudles is not None:
            print("误差:",self.data_tool.calculation_results_gap(moudles))
        self.data_tool.output_res(output_file=self.out_put_dir+"\%s"%("seed_result"))
        self.data_tool.output_smooth_res(output_file=self.out_put_dir+"\seed_result_smooth")
                
            
            
             
        
    
    def alt_seed_js(self, seed_js):
        # self.seed_js_history.append(copy(seed_js))
        self.seed_js=seed_js
    
    def output_Gx_file(self,moudles_file,Gx_file):
        """输出结果的d-Gx"""
        moudles=data_tools(moudles_file,self.mesh_file).get_data()
        d_GX=self.G*moudles
        d_GX=d_GX.tolist()
        d_Gx_file=open(Gx_file,'w')
        for i in d_GX:
            d_Gx_file.write(str(i))
            d_Gx_file.write("\n")
        d_Gx_file.close()
        
        
        
        

