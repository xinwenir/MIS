from math import sin ,cos
from InvDataTools.MeshTools import MeshTools
from InvDataTools.calculate.Detector import Detector
class Ray(object):
    def __init__(self,detector_id, detector_y,detector_x, detector_z, ray_theta, ray_phi,real_path_len,ray_d, ray_d_err,mesh_tool,loger):
        self._detector=Detector.get_instance(self,detector_id,detector_y, detector_x, detector_z,mesh_tool,loger)
        self.ray_theta=ray_theta
        self.ray_phi=ray_phi
        self.ray_d=ray_d
        self.real_path_len=real_path_len
        self.ray_d_err=ray_d_err
        self.mesh_tool=mesh_tool
        self._cell_js, self._cross_cell_length = None, None
        self._is_xy_through_ray=None
        
    @property
    def detector(self):
        """获得本射线所在的探测器信息
        #!用法:可以根据此探测器获取探测器的所有射线,已经此探测器的信息

        Returns:
            _type_: _description_
        """
        return self._detector 
    
    
    def _calc_through_j_values(self):
        """计算本射线穿过的体素编号和在每个体素中的真实路径长度
        """
        #计算这条射线穿过的格子的编号j以及每个格子的路径长度
        from InvDataTools.MeshTools import MeshTools
        self.mesh_tool:MeshTools=self.mesh_tool
        #获取边界的坐标；n个格子x方向有n+1个边界。
        ynode =self.mesh_tool.x_start_values
        ynode.append(self.mesh_tool.x_end_values[-1])
        xnode = self.mesh_tool.y_start_values
        xnode.append(self.mesh_tool.y_end_values[-1])
        znode =self.mesh_tool.z_start_values
        znode.append(self.mesh_tool.z_end_values[-1])
        znode.reverse()
        for i in range(len(znode)):
            znode[i]=self.mesh_tool.end_z-znode[i]

        my = self.mesh_tool.get_shape()[0]
        mx = self.mesh_tool.get_shape()[1]
        mz = self.mesh_tool.get_shape()[2]
        
        y0=self.detector.detector_x
        x0=self.detector.detector_y
        z0=self.detector.detector_z
        # cell_js存放射线穿过哪些格子，即j，sensrow存放射线穿过格子时，在格子内的路径长度
        cell_js, cross_cell_length = [], []
        # ux是x轴的，uy是y轴的，uz是z轴的，单位长度(可以看作是1)投影到各个轴上的长度，只需要乘以实际长度就可以得到实际的投影长度
        ux = abs(sin(self.ray_phi) * cos(self.ray_theta))
        uy = abs(cos(self.ray_phi) * cos(self.ray_theta))
        uz = sin(self.ray_theta)

        dirr = [0] * 2
        # dirr[0]是x轴的，dirr[1]是y轴的
        dirr[0] = -1 if sin(self.ray_phi) < 0 else 1
        dirr[1] = -1 if cos(self.ray_phi) < 0 else 1

        # 得到探测器的位置离散坐标
        yc=self.detector.get_detector_coordinates_x()
        xc=self.detector.get_detector_coordinates_y()
        zc=self.detector.get_detector_coordinates_z()

        # 射线在每个格子停留时，必定是停留在格子的某个面上
        # 如果射线与z轴的夹角很小
        if abs(uz - 1) < 1e-5:
            # 如果探测器几乎紧挨着格子壁，x0和y0就减小一点
            if abs(x0 - xnode[xc + 1] < 1e-5): x0 -= (xnode[xc + 1] - xnode[xc]) / 10
            if abs(y0 - ynode[yc + 1] < 1e-5): y0 -= (ynode[yc + 1] - ynode[yc]) / 10

            while True:
                zbound = znode[zc]
                dr = z0 - zbound
                z1 = z0 - dr * uz
                xcel = xc * my * mz + yc * mz + zc + 1
                cell_js.append(xcel)
                cross_cell_length.append(dr)

                zc -= 1
                z0 = z1
                if zc == -1:
                    break
        # 射线与z轴的夹角既不大也不小
        else:
            while True:
                # 找到射线此时所在的格子的z方向的上边界值
                try:
                    zbound = znode[zc]
                except Exception as e:
                    print(e)
                # 求射线距离上边界值的距离
                z1 = z0 - zbound
                # dz是斜边，z1是对边，dx,dy,dz的方向都是相同的，只是长短不一
                # 射线要到达z格子上边界所要走的路程
                dz = z1 / uz
                # 根据射线与y轴的角度给dx赋值
                # 与y轴的角度不是很小
                if abs(ux) > 1e-5:  # match考虑到python版本的问题，改用if
                    if dirr[0] == 1:
                        x1 = xnode[xc + 1] - x0
                    else:
                        x1 = x0 - xnode[xc]  # dirr[0]==-1说明与y轴的夹角在180°和360°之间，在x轴的负半轴
                    # 射线要到达x格子边界的所要走的路程
                    dx = x1 / ux
                # 如果射线与y轴的夹角很小，就说明射线接下来必定不会到达x格子边界，所以就把dx设的很大，选取dr的时候肯定不会选择dx，同时可以避免计算dx
                else:
                    dx = 1e32
                # 根据射线与y轴的角度给dy赋值
                # 与y轴的角度不是很大
                if abs(uy) > 1e-5:
                    if dirr[1] == 1:
                        y1 = ynode[yc + 1] - y0
                    else:
                        y1 = y0 - ynode[yc]
                    # 射线要到达y格子边界的所要走的路程
                    dy = y1 / uy
                # 如果射线与y轴的夹角很大(也即射线与x轴的夹角很小)，就说明射线接下来必定不会到达y格子边界，所以就把dy设的很大，选取dr的时候肯定不会选择dy，同时可以避免计算dy
                else:
                    dy = 1e32

                # dx,dy,dz哪个值赋值给dr，射线接下来就穿过哪个格子，e.g. dx最小，那么射线接下来就到达x格子边界
                dr = min(dx, dy, dz)
                # 计算接下来射线走dr路程到达的x1,y1位置,z1是射线走dr路程到达的位置到z方向最大值的差
                x1 = (x0 + dirr[0] * dr * ux)
                y1 = (y0 + dirr[1] * dr * uy)
                z1 = (z0 - dr * uz)
                # 计算射线此时所在格子的j值，因为xc,yc,zc是从0开始的，所以j也是从0开始的，但是需要j从1开始，所以要+1
                # 除去射线的第一个点，其余的点都在格子的表面上，每次到达的下一个点所在的面都是距离当前点最近的面
                # 当一个点在一个面上时，它归属于射线沿着穿出的方向所处的格子
                j = xc * my * mz + yc * mz + zc + 1

                # 存放射线穿过的格子ind以及在格子内的长度sensrow,大小都是numz,即射线穿过的格子数量
                cell_js.append(j)
                cross_cell_length.append(dr)

                # dx最小，说明射线到达了x格子边界，xc改变(根据射线的角度来决定加减)
                # dy最小，说明射线到达了y格子边界，yc改变(根据射线的角度来决定加减)
                # dz最小，说明射线到达了z格子边界，zc改变(只能减，即一直往上方靠)
                # 可能当前的点
                if dx < dy and dx < dz:
                    xc += dirr[0]
                elif dy < dx and dy < dz:
                    yc += dirr[1]
                elif dz < dx and dz < dy:
                    zc -= 1
                else:
                    if dx < dy:
                        xc += dirr[0]
                        zc -= 1
                    else:
                        yc += dirr[1]
                        zc -= 1

                x0, y0, z0 = x1, y1, z1
                
                if -1 < xc < mx and -1 < yc < my:
                    if zc == -1:
                        break
                else:
                    self._is_xy_through_ray=True
                    break
                    
        self._cell_js=cell_js
        self._cross_cell_length=cross_cell_length
            
    def get_through_cell_js(self):
        """
        获取射线穿过格子的编号

        Returns:
            list: 
        """
        if self._cell_js is None:
            self._calc_through_j_values()
        return self._cell_js
    
    def get_cross_cell_length(self):
        """
        射线在每个体素中穿过的长度

        Returns:
            list:
        """
        if self._cross_cell_length is None:
            self._calc_through_j_values()
        return self._cross_cell_length

    def get_is_xy_through_ray(self):
        """
        获取此射线是否是从xy面穿出的

        Returns:
            bool: 
        """
        if self._is_xy_through_ray is None:
            self._calc_through_j_values()
        return self._is_xy_through_ray
        
    
        