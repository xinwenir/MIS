
#演示seed算法

from InvDataTools.Gij_tools import G_Data
from InvDataTools.Jxyz_Tools import getj_from_xyz, getxyz_from_shape
from InvDataTools.d_tools import d_tools
from InvSolver.Seed_algorithm.Tools import data_tools

from InvDataTools.MeshTools import MeshTools

mesh_tool=MeshTools(r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\seed_test.msh")
shape=mesh_tool.get_shape()
new_seeds=set()
old_seeds=set()

new_seeds.add(int(180253))

new_seeds.add(int(180256))

data_tool=data_tools(r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\refs",r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\seed_test.msh")
data=data_tool.get_data()
for i in range(len(data)):
    data_tool.alt_data(i+1,0)
for i in range(1,6):
    middle=set()
    for j in new_seeds:
        if j in old_seeds:
            continue
        old_seeds.add(j)
        data_tool.alt_data(j,i)
        neighbor_js=data_tool.get_neighbor_js(j)
        for neighbor_j in neighbor_js:
           middle.add(neighbor_j)
    new_seeds=middle
    data_tool.output_res(r"E:\vscode\Muon_Imaging_Algorithm\data\seed_study\seed_step\step%s"%(str(i)))