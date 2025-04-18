﻿# author:高金磊
# datetime:2022/9/1 9:48


#简易目标函数
from ast import If
import imp
import re
from this import d
from tkinter.messagebox import NO
from turtle import pu
from unittest import result
import numpy.linalg.linalg
import scipy.sparse
from scipy.optimize import line_search
from solver.Seed_algorithm.Tools import data_tools


def obj_fun(G,x,d):
    return numpy.linalg.linalg.norm(G*numpy.array(x,float)-d)

class Obj_fun_Tools:

    #计算下降步长
    def __init__(self,G:scipy.sparse.csc_matrix,d):
        self.G=G
        self.j_rays={}
        self.d=d

    def update_G(self,G):
        self.G=G

    def get_j_optimal_value(self,data_tool:data_tools,j,bounds,max_iter=5,accuracy=0.2,tol=0.05):
        # min_value,max_value=bounds
        old_value=self.get_j_fun_value(data_tool,j)
        if old_value>tol:
            step=-accuracy
        elif old_value<-tol:
            step=accuracy
        else:
            #不做任何更改
            return data_tool.get_value(j)

        old_value=abs(old_value)
        old_j=data_tool.get_value(j)
        new_j=old_j+step
        while max_iter>0 and new_j>=bounds[0] and new_j<=bounds[1]:
            max_iter-=1
            data_tool.alt_data(j, new_j)
            new_value=self.get_j_fun_value(data_tool,j)
            if abs(new_value)<old_value:
                old_value=new_value
                old_j=new_j
                if new_value > tol:
                    step = -accuracy
                elif new_value < -tol:
                    step = accuracy
                else:
                    # 不做任何更改
                    return data_tool.get_value(j)
                new_j+=step
            else:
                data_tool.alt_data(j,old_j)
                break
        # flag=True
        # while max_iter>0 and new_j >= bounds[0]:
        #     max_iter-=1
        #     data_tool.alt_data(j,new_j)
        #     new_value=abs(self.get_j_fun_value(data_tool,j))
        #     if new_value<old_value:
        #         old_value=new_value
        #         old_j=new_j
        #         new_j+=step
        #         flag=False
        #     else:
        #         data_tool.alt_data(j,old_j)
        #         break
        # if flag:
        #     #反向寻找__
        #     step=accuracy
        #     new_j=old_j+step
        #     while max_iter > 0 and new_j <= bounds[1]:
        #         max_iter -= 1
        #         data_tool.alt_data(j, new_j)
        #         new_value = self.get_j_fun_value(data_tool, j)
        #         if new_value < old_value:
        #             old_value = new_value
        #             old_j = new_j
        #             new_j += step
        #         else:
        #             data_tool.alt_data(j, old_j)
        #             break
        return float(old_j)





    def get_j_fun_value(self,data_tool:data_tools,j):

        if j in self.j_rays.keys():
            rays_d=self.j_rays[j]
        else:
            rays_d = []
            rays_ids=self.G.getcol(j-1).nonzero()[0]
            for rays_id in rays_ids:
                rays_d.append([self.G.getrow(rays_id),self.d[rays_id]])
            self.j_rays[j]=rays_d
        res=0
        res_2=0
        x=data_tool.get_data()
        for data in rays_d:
            ray, d=data
            misfit_ray=(ray*x-d)#能体现梯度的算法
            res+=misfit_ray
            res_2+=misfit_ray**2
            # res += pow(ray * numpy.array(x, float) - d,2)#标准misfit算法
        #平滑性
        smooth=0
        neighbor=data_tool.get_neighbor_js(j)
        for i in neighbor:
            diff=(data_tool.get_value(j)-data_tool.get_value(i))
            smooth+=diff**2
            res+=diff/10
        if len(neighbor)>0:
            smooth/=len(neighbor)
        res_2+=smooth/10
        if res>0:
            return res_2
        else:
            return -res_2
    def get_j_line_search_optimal_value(self,data_tool:data_tools,j,bounds,start_seed,issmooth=False,distance=0.0,old_seeds=set()):  
        old_value=data_tool.get_value(j)
        result= line_search(f=self.line_search_fun,myfprime=self.line_search_grad,xk=data_tool.get_value(j),pk=-0.2,args=(data_tool,j,issmooth,distance,start_seed,old_seeds))

        threshold=0.8
        if result[0] is None:
            data_tool.alt_data(j,old_value)
            return old_value
        res=old_value-0.2*result[0]
        if res<max(bounds[0],old_value-threshold):
            data_tool.alt_data(j,max(bounds[0],old_value-threshold)+0.001)
            return max(bounds[0],old_value-threshold)+0.001
        if res>min(bounds[1],old_value+threshold):
            data_tool.alt_data(j,min(bounds[1],old_value+threshold)+0.001)
            return min(bounds[1],old_value+threshold)-0.001
        data_tool.alt_data(j,res)
        return res
    def line_search_fun(self,x,data_tool:data_tools, j,issmooth,distance,start_seed,old_seeds):
        data_tool.alt_data(j,x)
        if j in self.j_rays.keys():
            rays_d=self.j_rays[j]
        else:
            rays_d = []
            rays_ids=self.G.getcol(j-1).nonzero()[0]
            for rays_id in rays_ids:
                rays_d.append([self.G.getrow(rays_id),self.d[rays_id]])
            self.j_rays[j]=rays_d
        start_seed_value=data_tool.get_value(start_seed)
        punishment_value=data_tool.get_refs_value(j)-data_tool.get_value(j)
        res_2=abs((x-start_seed_value)/(0.5+x-start_seed_value))+distance*abs(punishment_value)
        x=data_tool.get_data()
        for data in rays_d:
            ray, d=data
            misfit_ray=(ray*x-d)#能体现梯度的算法
            res_2+=misfit_ray**2
            # res += pow(ray * numpy.array(x, float) - d,2)#标准misfit算法
        #平滑性
        if issmooth:
            smooth=0
            neighbor=data_tool.get_neighbor_js(j)
            for i in neighbor:
                if i  in old_seeds:
                    smooth+=(data_tool.get_value(j)-data_tool.get_refs_value(j)-data_tool.get_value(i)+data_tool.get_refs_value(i))**2
                else:
                    smooth+=0.05*(data_tool.get_value(j)-data_tool.get_refs_value(j)-data_tool.get_value(i)+data_tool.get_refs_value(i))**2
            res_2+=smooth/10
        return res_2
    
    def line_search_grad(self,x,data_tool:data_tools, j,issmooth,distance,start_seed,old_seeds):
        data_tool.alt_data(j,x)
        if j in self.j_rays.keys():
            rays_d=self.j_rays[j]
        else:
            rays_d = []
            rays_ids=self.G.getcol(j-1).nonzero()[0]
            for rays_id in rays_ids:
                rays_d.append([self.G.getrow(rays_id),self.d[rays_id]])
            self.j_rays[j]=rays_d
        res=0
        punishment_value=data_tool.get_value(j)-data_tool.get_refs_value(j)
        if punishment_value<0:
            res-=distance*0.02*punishment_value
        else:#值大于起始值,惩罚避免过大
            res+=distance*0.02*punishment_value
        x=data_tool.get_data()
        for data in rays_d:
            ray, d=data
            misfit_ray=(ray*x-d)*2
            res+=misfit_ray
        #平滑性
        if issmooth:
            neighbor=data_tool.get_neighbor_js(j)
            for i in neighbor:
                if  i in old_seeds: 
                    diff=(data_tool.get_value(j)-data_tool.get_refs_value(j)+data_tool.get_refs_value(i)-data_tool.get_value(i))
                    res+=diff/5
                else:
                    diff=0.05*(data_tool.get_value(j)-data_tool.get_refs_value(j)-data_tool.get_value(i)+data_tool.get_refs_value(i))**2
                    res+=diff/5
        
        return res
    ##多元目标函数
    def mult_obj_solver(self,data_tool:data_tools,j,bounds,issmooth=False): 
        # old_value=data_tool.get_value(j)
        # result= line_search(f=self.line_search_fun,myfprime=self.line_search_grad,xk=data_tool.get_value(j),pk=-0.2,args=(data_tool,j,issmooth))
        js=list(data_tool.get_neighbor_js(j))
        x0=[]
        js_bounds=[]
        pks=[-0.2]
        for jj in js:
            x0.append(data_tool.get_value(jj))
            js_bounds.append(bounds[jj-1])
            pks.append(-0.2)
        js.append(j)
        x0.append(data_tool.get_value(j))
        js_bounds.append(bounds[j-1])
        
        
        from scipy.optimize import minimize
        # result=minimize(self.mult_obj_fun,x0,method="L-BFGS-B",jac=self.mult_obj_grad,bounds=js_bounds,args=(js,data_tool))
        # return result.x[-1]
        import numpy
        alphs=line_search(self.mult_obj_fun,self.mult_obj_grad,x0,numpy.array(pks),args=(js,data_tool))
        result=[]
        for i in range(len(x0)):
            res=x0[i]
            if alphs[0] is not None:
                res-=0.2*alphs[0]
            if res<bounds[i][0]:
                data_tool.alt_data(j,bounds[i][0]+0.001)
                continue
            if res>bounds[i][1]:
                data_tool.alt_data(j,bounds[i][1]-0.001)
                continue
            data_tool.alt_data(j,res)
            result.append(res)
        return result
        
    def mult_obj_fun(self,js_value,js,data_tool:data_tools):
        #更新数据
        for i in range(len(js)):
            j=js[i]
            value=js_value[i]
            data_tool.alt_data(j,value)
        #获取所有相关的射线信息
        rays_d_set=[]
        for j in js:
            for ray_data in self.get_j_rays(j):
                rays_d_set.append(ray_data)    
        #计算目标函数
        res=0
        x=data_tool.get_data()
        for data in rays_d_set:
            ray, d=data
            misfit_ray=(ray*x-d)
            res+=misfit_ray**2
        return res
    def mult_obj_grad(self,js_value,js,data_tool:data_tools):
        #更新数据
        for i in range(len(js)):
            j=js[i]
            value=js_value[i]
            data_tool.alt_data(j,value)
        #获取所有相关的射线信息
        rays_d_sets=[]
        for j in js:
            rays_d_sets.append(self.get_j_rays(j))    
        #计算目标函数的梯度
        result=[]
        x=data_tool.get_data()
        for rays_d_set in rays_d_sets:
            res=0
            for data in rays_d_set:
                ray, d=data
                misfit_ray=(ray*x-d)
                res+=(misfit_ray*2)
            result.append(res)
        return result
        
    def get_j_rays(self,j):
        if j in self.j_rays.keys():
            rays_d=self.j_rays[j]
        else:
            rays_d = []
            rays_ids=self.G.getcol(j-1).nonzero()[0]
            for rays_id in rays_ids:
                rays_d.append([self.G.getrow(rays_id),self.d[rays_id]])
            self.j_rays[j]=rays_d
        return rays_d
        
