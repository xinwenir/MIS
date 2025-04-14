# author:高金磊
# datetime:2022/8/28 18:04
from copy import copy

import numpy
import scipy.sparse.csc
import tqdm


# def Find_similar_is(M:scipy.sparse.csc.csc_matrix):
#     ijs=M.nonzero()
#     nozero_js_list=[]
#     for i in range(ijs[0][-1]+1):
#         nozero_js_list.append(set())
#     for i in range(len(ijs[0])):
#         nozero_js_list[ijs[0][i]].add(ijs[1][i])
#     similar_is=[set() for i in range(len(nozero_js_list))]
#     for i in tqdm.trange(len(ijs[0])):
#         for j in range(i+1,len(ijs[0])):
#             if len(nozero_js_list[ijs[0][i]]) ==len(nozero_js_list[ijs[0][j]]) and nozero_js_list[ijs[0][i]] ==nozero_js_list[ijs[0][j]]:
#                 similar_is[ijs[0][i]].add(ijs[0][j])
#                 similar_is[ijs[0][j]].add(ijs[0][i])
#     return similar_is

def Find_similar_is(ijs):
    nozero_js_list=[]
    for i in range(ijs[0][-1]+1):
        nozero_js_list.append(0)
    for i in range(len(ijs[0])):
        nozero_js_list[ijs[0][i]]+=ijs[1][i]#统计每条射线的key
    nozero_js_is_map={}
    for i in range(len(nozero_js_list)):
        middle=nozero_js_is_map.get(nozero_js_list[i],set())
        middle.add(i)
        nozero_js_is_map[nozero_js_list[i]] = middle
    similar_is=[set() for i in range(len(nozero_js_list))]
    for key in tqdm.tqdm(nozero_js_is_map.keys()):
        for i in nozero_js_is_map[key]:
            similar_is[i]=copy(nozero_js_is_map[key])
            # similar_is[i].remove(i)
    return similar_is,nozero_js_is_map

# def Reducing_dimensions(vi,vj,value,d,d_err,nozero_js_is_map:dict):
#     m=len(nozero_js_is_map.keys())
#     nozero_oldis_newis={}
#     serial_number=0
#     for old_is in nozero_js_is_map.values():
#         #给旧的射线编号重新赋值
#         for old_i in old_is:
#            nozero_oldis_newis[old_i]=serial_number
#         serial_number+=1
#     new_d=[0 for i in range(len(nozero_oldis_newis))]
#     new_d_err = [0 for i in range(len(nozero_oldis_newis))]
#     from scipy.sparse import csc_matrix
#     M_matrix=csc_matrix(([], ([], [])), shape=(m, max(vj)+1))
#     for i in tqdm.trange(len(vi)):#将数据添加到稀疏矩阵中
#         x=nozero_oldis_newis[vi[i]]
#         y=vj[i]
#         v=value[i]
#         M_matrix._set_intXint(x,y,v+M_matrix._get_intXint(x,y))
#         new_d[x]+=d[vi[i]]
#         new_d_err[x]+=d_err[vi[i]]
#     new_i,new_j=M_matrix.nonzero()
#     new_value=M_matrix.data
#     return new_i,new_j,new_value,new_d,new_d_err


def Reducing_dimensions(vi,vj,value,d,d_err,nozero_js_is_map:dict,similar_is:list):
    m=len(nozero_js_is_map.keys())
    nozero_oldis_newis={}
    serial_number=0
    remove_duplicate=set()
    for old_is in similar_is:
        #给旧的射线编号重新赋值
        if not old_is.isdisjoint(remove_duplicate):
            continue
        for old_i in old_is:
           nozero_oldis_newis[old_i]=serial_number
           remove_duplicate.add(old_i)
        serial_number+=1
    new_d=[0 for i in range(len(nozero_oldis_newis))]
    new_d_err = [0 for i in range(len(nozero_oldis_newis))]

    i_js_value=[{} for i in range(m)]

    for i in tqdm.trange(len(vi)):#将数据添加到稀疏矩阵中
        x=nozero_oldis_newis[vi[i]]
        y=vj[i]
        v=value[i]
        i_js_value[x][y]=i_js_value[x].get(y,0)+v
    for i in range(len(d)):
        new_d[nozero_oldis_newis[i]]+=d[i]
        new_d_err[nozero_oldis_newis[i]]+=d_err[i]
    new_i, new_j, new_value=[],[],[]
    for i in range(len(i_js_value)):
        js_value=i_js_value[i]
        for j in js_value.keys():
            new_i.append(i)
            new_j.append(j)
            new_value.append(js_value[j])


    return new_i,new_j,new_value,new_d,new_d_err




