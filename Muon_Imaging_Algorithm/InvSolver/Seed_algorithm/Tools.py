# author:高金磊
# datetime:2022/9/1 9:48
from copy import copy

import numpy
import tqdm as td
from InvDataTools.Jxyz_Tools import getj_from_xyz, getxyz_from_shape
from InvDataTools.MeshTools import MeshTools


class data_tools:
    def __init__(self, data_file,mesh_file,refs_file=None):
        data = []
        for i in open(data_file, 'r').readlines():
            if i != "":
                rho = float(i)
                # if rho <= 0.2:
                #     rho = -0.1234
                # else:
                #     rho = 2.65
                data.append(rho)
        self.data=numpy.array(data,float)
        self.mesh = MeshTools(mesh_file)
        self._check_()
        # 将数据还原到三维空间
        self._data_3D = None
        self.j_xyz_cache = {}
        self._xyz_j=None
        if refs_file==None:
            self.refs=copy(self.data)
        else:
            refs_data = []
            for i in open(refs_file, 'r').readlines():
                if i != "":
                    rho = float(i)
                    # if rho <= 0.2:
                    #     rho = -0.1234
                    # else:
                    #     rho = 2.65
                    refs_data.append(rho)
            self.refs=numpy.array(refs_data,float)
            if len(self.refs) !=len(self.data):
                raise Exception("数据文件和参考文件不对应")
    def _check_(self):
        if len(self.data) != self.mesh.cells_count():
            raise Exception("数据文件和mesh文件不对应")

    def get_data_3D(self):
        if self._data_3D is None:
            self._data_3D = numpy.zeros(self.mesh.get_shape())
            self._xyz_j = numpy.zeros(self.mesh.get_shape(),dtype=int)
            # data=copy(self.data)

            for y in range(self.mesh.get_shape()[1]):
                for z in range(self.mesh.get_shape()[2]):
                    for x in range(self.mesh.get_shape()[0]):
                        j=getj_from_xyz(self.mesh.get_shape(), (x+1, y+1, z+1))
                        self._data_3D[x][y][z] = self.data[j-1]
                        self._xyz_j[x][y][z] = j
                        # self._data_3D[x][y][z] = data[j]
        return self._data_3D

    def get_xyz_j(self):
        if self._xyz_j is None:
            self.get_data_3D()
        return self._xyz_j

    def alt_data(self, j, value):
        """修改数据

        Args:
            j (int): 格子编号--start by 1
            value (float): 格子的新值

        Returns:
            float: 格子被修改之前的值
        """        
        old_value=self.get_value(j)
        self.data[j - 1] = value
        return old_value

    def update_data(self, new_data):
        self.data = new_data

    def get_data(self):
        return self.data

    def get_value(self, j):
        return self.data[j - 1]

    def get_refs_value(self, j):
        return self.refs[j - 1]

    def get_j_from_xyz(self, xyz):
        # try:
        j = self.get_xyz_j()[xyz[0]-1][xyz[1]-1][xyz[2]-1]
        self.j_xyz_cache[j] = xyz
        # except Exception as e:
        #     print(e)
        return j

    def get_xyz_from_j(self, j):
        if j in self.j_xyz_cache.keys():
            return self.j_xyz_cache[j]
        else:
            xyz = getxyz_from_shape(self.mesh.get_shape(), j)
            self.j_xyz_cache[j] = xyz
            return xyz
    def calculation_results_gap(self,standard_answer_values):
        res=0
        for i in range(len(standard_answer_values)):
            res+=(self.data[i]-standard_answer_values[i])**2
        return pow(res, 1)
        
    def get_neighbor_js(self, j,modle=1):
        """
        获取j附近的格子的j
        :param j:
        :return:
        """
        x, y, z = self.get_xyz_from_j(j)
        res = set()
        # x-=1
        # y-=1
        # z-=1
        
        if modle==0:
            if x-1>=0:
                res.add(self.get_j_from_xyz((x-1, y, z)))
            if x+1<=self.mesh.get_shape()[0]:
                res.add(self.get_j_from_xyz((x+1, y, z)))
            if y-1>=0:
                res.add(self.get_j_from_xyz((x, y-1, z)))
            if y+1<=self.mesh.get_shape()[1]:
                res.add(self.get_j_from_xyz((x, y+1, z)))
            if z-1>=0:
                res.add(self.get_j_from_xyz((x, y, z-1)))
            if z+1<=self.mesh.get_shape()[2]:
                res.add(self.get_j_from_xyz((x, y, z+1)))
        elif modle==1:
            #周围26个
            for xx in range(x - 1, x + 2):
                if xx < 0 or xx >= self.mesh.get_shape()[0]:
                    continue
                for yy in range(y - 1, y + 2):
                    if yy < 0 or yy >= self.mesh.get_shape()[1]:
                        continue
                    for zz in range(z - 1, z+2):
                        if zz < 0 or zz >= self.mesh.get_shape()[2]:
                            continue
                        # if abs(xx-x)+abs(yy-y)+abs(zz-z)==1:
                        res.add(self.get_j_from_xyz((xx, yy, zz)))
        return res

    def output_res(self, output_file):
        output_obj = open(output_file, 'w')
        for i,j in zip(self.data,self.refs):
            if False:####todo临时使用，排除某个值
                # if i==2.65:
                    output_obj.write(str(-0.2))
                # else:
                #     output_obj.write("1.0")
            else:
                output_obj.write(str(i))
            output_obj.write('\n')
        # for i in self.data:#todo 临时关闭
        #     output_obj.write(str(i))
        #     output_obj.write('\n')
        output_obj.close()

    def output_smooth_res(self,output_file):
        for i in range(3):
            for j in td.trange(len(self.data)):
                old_value=self.get_value(j)
                if 0.1<old_value and abs(old_value-2.65)>=0.0:#注意异常体是低密度和高密度，一般通用的不好用
                    neighbor_cells=self.get_neighbor_js(j,modle=0)
                    count=0
                    smooth=0
                    for neighbor_j in neighbor_cells:
                        neighbor_cell_value=self.get_value(neighbor_j)
                        if neighbor_cell_value>0.1: #and neighbor_cell_value!=2.65:  #<2.1:
                            if abs(neighbor_cell_value-2.65)<0.06:
                                smooth+=1*neighbor_cell_value
                                count+=1
                            else:
                                smooth+=2*neighbor_cell_value
                                count+=2
                    if count<=8:
                        count=0
                    
                    smooth+=1*count*self.get_value(j)*1
                    count*=2
                    if count!=0:
                        self.alt_data(j,smooth/count)
        self.output_res(output_file=output_file)

    def output_smooth_res_for_wall(self,output_file):
        for i in range(2):
            for j in td.trange(len(self.data),-1,-1):
                old_value=self.get_value(j)
                ref_value=self.get_refs_value(j)
                # if 0<old_value and old_value<=2.65:#仿真
                if ref_value==1.89:
                    neighbor_cells=self.get_neighbor_js(j)
                    count=0
                    smooth=0
                    for neighbor_j in neighbor_cells:
                        neighbor_cell_value=self.get_value(neighbor_j)
                        if neighbor_cell_value>=0 and self.get_refs_value(neighbor_j)==1.89 and neighbor_cell_value<1.8:
                            smooth+=neighbor_cell_value
                            count+=1
                    smooth*=1
                    smooth+=count*self.get_value(j)*1.2
                    if count!=0:
                        self.alt_data(j,smooth/count/2.2)
            for j in td.trange(len(self.data)):
                old_value=self.get_value(j)
                ref_value=self.get_refs_value(j)
                # if 0<old_value and old_value<=2.65:#仿真
                if ref_value==1.89:
                    neighbor_cells=self.get_neighbor_js(j)
                    count=0
                    smooth=0
                    for neighbor_j in neighbor_cells:
                        neighbor_cell_value=self.get_value(neighbor_j)
                        if neighbor_cell_value>=0 and self.get_refs_value(neighbor_j)==1.89 and neighbor_cell_value<1.9:
                            smooth+=neighbor_cell_value
                            count+=1
                    smooth*=1
                    smooth+=count*self.get_value(j)*3
                    if count!=0:
                        self.alt_data(j,smooth/count/4)
        for j in td.trange(len(self.data),-1,-1):
            ref_value=self.get_refs_value(j)
            if ref_value==1.89:
                if abs(self.get_value(j)-1.89)<0.8:
                        self.alt_data(j, 1.89)
        self.output_res(output_file=output_file)
        


