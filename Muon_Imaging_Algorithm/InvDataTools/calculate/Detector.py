import time
from InvDataTools.MeshTools import MeshTools
detectors_map={}
class Detector(object):
    _identification=21792071980
    def __init__(self,detector_id,detector_x, detector_y, detector_z,mesh_tool,identification_key,loger) -> None:
        self.detector_id=detector_id
        self.detector_x= detector_x
        self.detector_y=detector_y
        self.detector_z=detector_z
        self.rays=set()
        if identification_key != Detector._identification:
            if loger == None:
                print("Detector不建议直接生成,请使用get_instance来避免程序出错;该对象已经被舍弃")
            else:
                loger.waring("Detector不建议直接生成,请使用get_instance来避免程序出错;该对象已经被舍弃")
        self._check_(mesh_tool)
            
        
    @staticmethod
    def get_instance(ray,detector_id,detector_y, detector_x, detector_z,mesh_tool,loger):
        """
        获取detector对象,通过detector_id来保证复用,降低不必要的计算
        Args:
            detector_id (_type_): 探测器编号 任意类型
            detector_y (_type_): _description_
            detector_x (_type_): _description_
            detector_z (_type_): _description_
            mesh_tool (_type_): _description_
            loger (_type_): _description_

        Returns:
            _type_: _description_
        """
        
        if detector_id in detectors_map.keys():
            detector=detectors_map[detector_id]
            if detector.detector_z!=detector_z or detector.detector_x!=detector_y or detector.detector_y!=detector_x:
                raise Exception("id是探测器的唯一标识,当前id对应的探测器编号已经被占用,请检查该编号下的射线,%s %s %s"%(str(detector_y),str(detector_x),str(detector_z)))
            detector.add(ray)
            return detector
        detectors_map[detector_id]=Detector(detector_id,detector_y, detector_x, detector_z,mesh_tool,Detector._identification,loger)
        detectors_map[detector_id].add(ray)
        return detectors_map[detector_id]
    @staticmethod
    def distory():
        print("已经销毁探测器信息缓存")
        detectors_map.clear()
        
    def _check_(self,mesh_tool:MeshTools):
        """检查探测器是否在成像空间内

        Args:
            mesh_tool (MeshTools): 

        Raises:
            ValueError: 探测器不在mesh边界之内
        """
        #获取三个方向是第几个格子
        coordinates_x=mesh_tool.discretize_Physical_coordinates_x(self.detector_x)[-1]
        coordinates_y=mesh_tool.discretize_Physical_coordinates_y(self.detector_y)[-1]
        coordinates_z=mesh_tool.get_shape()[2]-mesh_tool.discretize_Physical_coordinates_z(mesh_tool.end_z-self.detector_z)[-1]
        if coordinates_x==-1 or coordinates_y==-1 or coordinates_z==-1 or coordinates_z>mesh_tool.get_shape()[2]:
            raise ValueError("探测器不在mesh边界之内:(%s,%s,%s)"%(str(self.detector_x),str(self.detector_y),str(self.detector_z)))
        #得到三个方向格子的编号
        coordinates_x-=1
        coordinates_y-=1
        # coordinates_z-=1
        self._detector_coordinates_x=coordinates_x
        self._detector_coordinates_y=coordinates_y
        self._detector_coordinates_z=coordinates_z
    
    def get_detector_coordinates_x(self):
        return self._detector_coordinates_x
    
    def get_detector_coordinates_z(self):
        return self._detector_coordinates_z
    
    def get_detector_coordinates_y(self):
        return self._detector_coordinates_y
    
    def add(self,ray):
        """记录那些射线属于这个探测器

        Args:
            ray (_type_): _description_
        """
        self.rays.add(ray)

    def get_rays(self):
        return self.rays
    
    @staticmethod
    def get_all_detector():
        """获取程序运行期间所有生成的探测器对象

        Returns:
            list:
        """
        return detectors_map.values()

        
        
        
        
        

        