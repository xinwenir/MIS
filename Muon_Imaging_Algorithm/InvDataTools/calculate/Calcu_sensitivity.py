# 如果对数据要求精确的话，使用Decimal
# from decimal import Decimal
from math import pi
from tqdm import tqdm
from InvDataTools.MeshTools import MeshTools
from InvDataTools.calculate.Rays import Ray
from InvDataTools.calculate.Detector import Detector

class Calcsensitivity:
    def __init__(self, mesh_tool:MeshTools,loger=None,obs_file=None):
        Detector.distory()
        self.mesh_tool=mesh_tool
        self.loger=loger
        self.ray_data=[]
        self.obs_file_info=''
        if obs_file is not None:
            #加载obs数据
            file=open(obs_file,'r')
            data=file.readlines()
            #第一行为文件信息
            self.obs_file_info=data[0]
            for s in data[1:]:
                s=s.replace("\n","")
                
                middle=s.replace("\n","").split(" ")
                if len(middle)!=9:
                    if loger==None:
                        print("射线数据缺失:%s"%(s))
                    else:
                        loger.err("射线数据缺失:%s"%(s))
                    continue
                det_des, y0, x0, z0, theta0, phi0,real_path_len, d, d_err = middle[0], float(middle[1]),float(middle[2]),self.mesh_tool.end_z-float(middle[3]), pi / 2 - float(middle[4]),float(middle[5]),float(middle[6]), float(middle[7]), float(middle[8])
                # for i in range(1,len(middle)):
                #     #除了第一列为不定类型使用字符串保存,其他使用float
                #     middle[i]=float(middle[i])
                ray=Ray(det_des, y0, x0, z0, theta0, phi0,real_path_len, d, d_err,mesh_tool,loger)
                self.ray_data.append(ray)
            file.close()
        
        self.rays_number=len(self.ray_data)
    
    def calc_single_ray_path(self,ray):
        """计算单条射线穿过的格子和穿过格子的路径长度
        本方法仅仅是为了组织ray对象中的两个方法
        #!原则上建议直接对ray对象进行操作
        Args:
            ray (Ray): 射线类对象,包含射线相关的所有操作和数据

        Returns:
            _type_: 两个list,分别为穿过格子的编号和在每个格子的路径长度
        """
        return ray.get_through_cell_js(),ray.get_cross_cell_length()
    
    def calc_all_rays_from_obs_file(self,Gij_file,ray_way_j_file,rm_xy_through_ray=False):
        """根据obs文件生成反演中需要的两个文件
        Gij和ray_way_j

        Args:
            Gij_file (String): Gij要被保存的文件地址
            ray_way_j_file (String): ray_way_j要被保存的文件地址
            rm_xy_through_ray (bool, optional): 是否移除从x和y面穿过的射线. Defaults to False.
        """
        
        if rm_xy_through_ray:
            if self.loger == None:
                print("注意,已经开启了去除xy面穿过的射线")
            else:
                self.loger.waring("注意,已经开启了去除xy面穿过的射线")
        count=0 #可用的射线数量
        all_count=0 #所有射线
        ray_way_j_data=[]
        Gij_data=[]
        for ray in tqdm(self.ray_data):
            cell_js,cell_length=self.calc_single_ray_path(ray)
            ray_way_j_data.append(cell_js)
            all_count+=1
            if not (ray.get_is_xy_through_ray() and rm_xy_through_ray):
                count+=1
                for i in range(len(cell_js)):
                    Gij_data.append([all_count,cell_js[i],cell_length[i]])
                
        
        #写文件
        Gij_file_obj=open(Gij_file,'w')
        Gij_file_obj.write("%s %s"%(str(all_count),str(self.mesh_tool.cells_count())))#记录数量
        Gij_file_obj.write("\n")
        for Gij in Gij_data:
            Gij_file_obj.write("%s %s %s"%(Gij[0],Gij[1],Gij[2]))
            Gij_file_obj.write("\n")
        Gij_file_obj.close()
        ray_way_j_file_obj=open(ray_way_j_file,'w')
        for js in ray_way_j_data:
            for j in js:
                ray_way_j_file_obj.write(str(j))
                ray_way_j_file_obj.write(" ")
            ray_way_j_file_obj.write("\n")
        ray_way_j_file_obj.close()
        if count!=len(self.ray_data):
            if self.loger==None:
                print("存在因为从xy面穿出被舍弃的射线,共计%s个:"%(len(self.ray_data)-count))
            else:
                self.loger.waring("存在因为从xy面穿出被舍弃的射线,共计%s个:"%(len(self.ray_data)-count))
            
               