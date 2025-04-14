import copy
from InvSolver.Seed_algorithm.paper.Tools import data_tools
import scipy.sparse
import numpy as np
from scipy.optimize import minimize
from InvSolver.Seed_algorithm.paper.MNI import MNI_solver
from InvSolver.Seed_algorithm.paper.Object_fun import Obj_fun_Tools
class BCD_Opt_Tools:
    def __init__(self,G:scipy.sparse.csc_matrix,d,d_err,refs,bounds):
        self.G=G
        self.j_rays={}
        self.d=d
        self.d_err=d_err
        self.refs=refs
        self.bounds=bounds
    def init_BCD_opt_question(self,data_tool:data_tools,seed_js,threshold=0):
        G=self.G
        G_middle=[]
        new_G=[]
        new_x0=[]
        new_d=[]
        new_refs=[]
        new_bounds=[]
        new_derrs=[]
        rays_id_set=set()
        seeds_neighbor_js=[]
        
        #将待求解的项目设置为0，这样可以求出新问题中的d
        middle_x=copy.copy(data_tool.get_data())
        for j in seed_js:
            middle_x[j-1]=0
        
        #处理射线构造新的G
        for j in seed_js:
            if j in self.j_rays.keys():
                rays_ids=self.j_rays[j]
            else:  
                rays_ids=self.G.getcol(j-1).nonzero()[0]
                for i in rays_ids:
                    rays_id_set.add(i)
                self.j_rays[j]=rays_ids
        #G相关的射线信息
        for rays_id in rays_id_set:
            G_middle.append(self.G.getrow(rays_id))
            new_d.append(self.d[rays_id])
            new_derrs.append(self.d_err[rays_id])
        #固定其他体素密度值得到，非稀疏矩阵G和对应的d
        for i in range(len(G_middle)):
            ray=G_middle[i].toarray()[0]
            new_d[i]-=np.dot(ray,middle_x)
            ray_middle=[]
            for j in seed_js:
                ray_middle.append(ray[j-1])
            new_G.append(ray_middle)
        #创建其他相关的项
        for j in seed_js:
            #混合上下限约束和每次最大下降步长约束
            bound=copy.copy(self.bounds[j-1])
            bound[0]=max(bound[0],data_tool.get_value(j)-threshold)
            # bound[0]=0
            bound[1]=min(bound[1],data_tool.get_value(j)+threshold)
            new_bounds.append(bound)
            new_refs.append(self.refs[j-1])
            new_x0.append(data_tool.get_value(j))
            #每个格子的邻居的下标
            neighbor_js=[]
            for j in data_tool.get_neighbor_js(j,modle=0):
                neighbor_js.append(j)
                # data_tool.get_value(j)-data_tool.get_refs_value(j)
            seeds_neighbor_js.append(neighbor_js)
        return np.array(new_G),np.array([new_x0]).T,np.array([new_d]).T,np.array([new_derrs]).T,new_refs,new_bounds,seeds_neighbor_js
    
    def  optimize_single_Block(self,obj_fun_Tool:Obj_fun_Tools,seed_js,G,x0,d,derrs,refs,bounds,distance,ancestors_seeds_js,issmooth,seeds_neighbor_js,data_tool):
        """
        执行优化算法对指定块进性优化

        Args:
            G (_type_): _description_
            x0 (_type_): _description_
            d (_type_): _description_
            derrs (_type_): _description_
            refs (_type_): _description_
            bounds (_type_): _description_
        """
        obj_fun_Tool.clear_history()
        
        opt_tool=1
        if opt_tool==0:
            res=minimize(fun=obj_fun_Tool.objective_fun,
                        # jac=obj_fun_Tool.objective_gun,# 由于罚函数的出现导致提供的jac效果反而不好
                        method="L-BFGS-B", 
                    x0=x0,
                    bounds=bounds,
                    #  options={"maxfun": max(10,len(x0)),'maxiter':max(10,len(x0)), 'disp':True},
                    args=(seed_js,G,d,derrs,refs,bounds,distance,ancestors_seeds_js,issmooth,seeds_neighbor_js,data_tool,False)
                    )
            #由于精度损失，不一定达到期望的误差
            res.x=np.array([obj_fun_Tool.min_fun_x]).T
            res.fun=obj_fun_Tool.min_fun
            return res
        elif opt_tool==1:
            MNI= MNI_solver(obj_fun_Tool,seed_js,G,x0,d,derrs,refs,bounds,distance,ancestors_seeds_js,issmooth,seeds_neighbor_js,data_tool)
            eta=0.5
            err_integral=None
            for i in range(10):
                if err_integral is None:
                    G,d,err_integral=MNI.convert_problems_to_matrix(x0,eta=eta)
                else:
                    G,d,err_integral_buyao=MNI.convert_problems_to_matrix(x0,eta=eta)
                """MNI"""
                x0,err_integral=MNI_solver.solver(G,d,x0,eta,err_integral)
                
            res=x0
            for i in range(len(res)):
                if res[i]<bounds[i][0]:
                    res[i]=bounds[i][0]
                if res[i]>bounds[i][1]:
                    res[i]=bounds[i][1]
            from scipy.optimize import OptimizeResult
            result=OptimizeResult()
            result.x=res
            return result
        elif opt_tool==2:
            ##和论文保持一致
            res=minimize(fun=obj_fun_Tool.objective_fun,
                        # jac=obj_fun_Tool.objective_gun,# 由于罚函数的出现导致提供的jac效果反而不好
                        method="BFGS", 
                        x0=x0,
                        # bounds=bounds,
                         options={'disp':True},
                        args=(seed_js,G,d,derrs,refs,bounds,distance,ancestors_seeds_js,issmooth,seeds_neighbor_js,data_tool,True,True)
                        )
            res.x=np.array([obj_fun_Tool.min_fun_x]).T
            res.fun=obj_fun_Tool.min_fun
            return res
        elif opt_tool==3:
            from scipy.optimize import lsq_linear
            import random
            MNI= MNI_solver(obj_fun_Tool,seed_js,G,x0,d+0.8,derrs,refs,bounds,distance,ancestors_seeds_js,issmooth,seeds_neighbor_js,data_tool)
            eta=0.0
            # MNI.d-=1.0 
            for i in range(30):
                
                G,d,err_integral=MNI.convert_problems_to_matrix(x0,eta=eta)
                
                # d-=1.0 #恒定噪声--数据
                d-=1*random.random() 
                # res = lsq_linear(G, d.T.tolist()[0], bounds=([i[0] for i in bounds], [i[1] for i in bounds]), lsmr_tol='auto', verbose=1)
                res = lsq_linear(G, d.T.tolist()[0], lsmr_tol='auto', verbose=1)
                x0=np.array([res.x]).T
                # print(res.x)
            return res
            
            