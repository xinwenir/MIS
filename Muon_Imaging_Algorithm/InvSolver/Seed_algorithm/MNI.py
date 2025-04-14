# """增加约束"""
# G=np.array([[1.0,2.0,3.0],
#    [4.0,5.0,6.0]])

# x0=np.array([[1.0,2.0,3.0]]).T

# refs=np.array([[0.1,0.1,0.1]]).T
# ref_c=10.0

# d=np.array([[2.0],[3.0]])

# uper=0.2
# low=-0.3
# phi_c=10.0

# constraint=np.zeros((len(x0),len(x0)))
# constraint_d=np.zeros((len(x0),1))
# refs_constraint=np.diag([phi_c for i in x0])




# err_integral=np.zeros(shape=(len(d)+len(x0)+len(refs),1))
# count=0

# while count<20:
#     for i in range(len(x0)):
#         x_i=x0[i]
#         if x_i<=uper and x_i>=low:
#             constraint[i][i]=0
#             constraint_d[i][0]=0
#         else:
#             if x_i<low:
#                 constraint[i][i]=-phi_c
#                 constraint_d[i][0]=-low*phi_c
#             else:
#                 constraint[i][i]=phi_c
#                 constraint_d[i][0]=uper*phi_c
#     ref_d_constraint=np.array([[i[0]*ref_c] for i in refs])
            
#     G_final=np.vstack((G,constraint,refs_constraint))
#     d_final=np.vstack((d,constraint_d,ref_d_constraint))
      
    
#     err=np.dot(G_final,x0)-d_final
#     print(np.linalg.norm(err))
    
#     G_inv=np.linalg.pinv(G_final)
#     x0=x0-np.dot(G_inv,(err
#                          + 0.0001*err_integral
#                         ))
#     err_integral+=err
#     count+=1
# print(x0)

import copy
from InvSolver.Seed_algorithm.objective_function import Obj_fun_Tools
import numpy as np
from InvSolver.Seed_algorithm.Setting import *
import numpy
import random
class MNI_solver():
    def __init__(self,obj_fun_Tool:Obj_fun_Tools,seed_js,G,x0,d,derrs,refs,bounds,distance,ancestors_seeds_js,issmooth,seeds_neighbor_js,data_tool) -> None:
        self.obj_fun_Tool=obj_fun_Tool
        self.seed_js=seed_js
        self.G=G
        self.x0=x0
        self.d=d
        self.derrs=derrs
        self.refs=refs
        self.bounds=bounds
        self.distance=distance
        self.seeds_neighbor_js=seeds_neighbor_js
        self.issmooth=issmooth
        self.ancestors_seeds_js=ancestors_seeds_js
        self.data_tool=data_tool
    
    def convert_problems_to_matrix(self,x0,eta=0.5):
        deta=0.001
        min_err=10000000
        min_x0=None
        self.x0=x0
        m=len(self.d)
        n=len(x0)
        print("p:%s  q:%s"%(str(m),str(n)))
        #构造优化问题
        ##misfit--就是矩阵本身
        G_final=copy.copy(self.G)
        d_final=copy.copy(self.d)
        ##初始化完整的矩阵
        if self.issmooth:
            G_final=np.vstack((G_final,np.zeros((4*n,n))))
            d_final=np.vstack((d_final,np.zeros((4*n,1))))
            err_integral=np.zeros(shape=(4*n+m,1))
            
        else:
            
            G_final=np.vstack((G_final,np.zeros((3*n,n))))
            d_final=np.vstack((d_final,np.zeros((3*n,1))))
            err_integral=np.zeros(shape=(3*n+m,1))
            
        
        count=0
        
        while count<5:
            ##罚函数限制约束和单次下降距离
            middle_i=m
            for i in range(len(x0)):
                if x0[i][0]<=self.bounds[i][1] and x0[i][0]>=self.bounds[i][0]:
                    G_final[i+middle_i][i]=0
                    d_final[i+middle_i][0]=0
                else:
                    if x0[i][0]<self.bounds[i][0]:
                        G_final[i+middle_i][i]=punishment_factor
                        d_final[i+middle_i][0]=self.bounds[i][0]*punishment_factor
                    else:
                        G_final[i+middle_i][i]=punishment_factor
                        d_final[i+middle_i][0]=self.bounds[i][1]*punishment_factor
            ##搜索距离限制
            middle_i+=n
            for i in range(len(x0)):
                G_final[i+middle_i][i]=self.distance*punishment_Search_distancevalue_multiple
                d_final[i+middle_i][0]=self.refs[i]*(self.distance)*punishment_Search_distancevalue_multiple
            ##refs约束_取消该项
            # for i in range(n):
            #     G_final[i+middle_i][i]=0
            #     d_final[i+middle_i][0]=self.refs[i]*0
            
            ##与祖先种子密度差
            middle_i+=n  
            for i in range(n):
                G_final[i+middle_i][i]=prompt_value_multiple
                d_final[i+middle_i][0]=self.data_tool.get_value(self.ancestors_seeds_js[i])*prompt_value_multiple
                
            
            ##平滑性
            if self.issmooth:
                middle_i+=n 
                for i in range(n):
                    G_final[i+middle_i][i]=smooth_all_multiple
                    for j in self.seeds_neighbor_js[i]:
                        d_final[i+middle_i][0]+=self.data_tool.get_value(j)*smooth_all_multiple
                    d_final[i+middle_i][0]/=len(self.seeds_neighbor_js[i])
                    # d_final[i+middle_i][0]+=x0[i]*smooth_all_multiple*len(self.seeds_neighbor_js[i])
                    
            count+=1
            err_integral*=0
            """MNI"""
            for t in range(10):
                random_noise=2*random.random()
                err=np.dot(G_final,x0)-d_final
                
                # err_value=numpy.linalg.norm(err)
                # # print(err_value)
                # if err_value<min_err:#获取迭代中的最优解
                #     min_x0=copy.copy(x0)
                #     min_err=err_value
                    
                
                G_inv=np.linalg.pinv(G_final)
                err_integral+=np.dot(G_final,x0)-d_final
                x0=x0-np.dot(G_inv,err +
                                    #   (t*deta*100+10) + #线性噪声
                                    1+
                                    # random_noise+ #随机噪声
                                    eta*err_integral
                                    )
                
                # if self.issmooth:
                #     print(numpy.linalg.norm(err))
                
                ##更新data_tool
                # ii=0
                # for seed_j in self.seed_js:
                #     self.data_tool.alt_data(seed_j,x0[ii])
                #     ii+=1
            
        return x0
        
        
        
        
        