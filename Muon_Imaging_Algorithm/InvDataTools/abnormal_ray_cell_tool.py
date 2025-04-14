from InvDataFactory import DataManage
from InvDataTools.Cell_Weight_Tools import Cell_Weight_Tools
from InvDataTools.ref_tools import Ref_tools
from InvDataTools.obs_tools import obs_tools
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import numpy as np
# from InvDataTools.bigCellsTool import BigCells


class AbnormalRayCellTool:
    def __init__(self) -> None:
        obsTool = DataManage.DataManager.get_instance().obs_tool
        phyLengths = [float(i) for i in obsTool.get_phyd_from_obs()]
        obsLengths = [float(i) for i in obsTool.get_d_form_obs()] 
        self.rayDensitys = [[index +1, obsLengths[index] / phyLengths[index]] for index in range(len(phyLengths)) if phyLengths[index] >0]        
        self.abnormalRays = None
        self.cellDeensities = None

    def getRayDensityMinMedianMax(self,minPercent=0.1, maxPercent=0.1):
        rayDensitys = sorted(self.rayDensitys, key= lambda x:x[1])
        rayNum = len(self.rayDensitys)
        return [rayDensitys[int(minPercent * rayNum)], 
                rayDensitys[int(rayNum/2)], 
                rayDensitys[int((1 - maxPercent) * rayNum) - 1]]
    
    def getAbnormalRaysCells(self, thresholdRayNum=1, minDensity = 0.8, maxDensity=2.5):
        if self.abnormalRays is None:
            dataManager = DataManage.DataManager.get_instance()
            obsData = dataManager.obs_tool.get_data()
            rayWayJs = dataManager.cell_weight_tool.data
            # refData = dataManager.refs_tool.get_data()
            cellCrossedRaysInfo = [[0, set(), set()] for _ in range(dataManager.mesh.cells_count())]
            # 格子数，射线的探测器编号，射线编号   
            abnormalRays = []
            for rayDensity in self.rayDensitys:
                if not minDensity < rayDensity[1] < maxDensity:
                    abnormalRays.append(obsData[rayDensity[0] - 1])
                    detectorId = obsData[rayDensity[0] - 1][0]
                    for cellId in rayWayJs[rayDensity[0] - 1]:                         
                        if detectorId in cellCrossedRaysInfo[cellId - 1][1]:
                            # pass
                            cellCrossedRaysInfo[cellId - 1][0] += 0.1
                        else:
                             cellCrossedRaysInfo[cellId - 1][1].add(detectorId)
                             cellCrossedRaysInfo[cellId - 1][0] += 1
                        cellCrossedRaysInfo[cellId - 1][2].add(rayDensity[1])
            self.cellDensities = [sum(i[2]) / len(i[2]) if (i[0] > thresholdRayNum and len(i[1]) > 1) else 3 for i in  cellCrossedRaysInfo]
            self.abnormalRays = abnormalRays
        return self.abnormalRays,self.cellDensities
    
    def getAbnormalRaysByPercent(self,minPercent=0.1, maxPercent=0.1):
        minDensity,medianDensity,maxDensity = self.getRayDensityMinMedianMax(minPercent, maxPercent)
        self.getAbnormalRaysCells(minDensity=minDensity[1],maxDensity=maxDensity[1])
        
    def saveAbnormalRaysCells(self):
        abnormalRays, cellDensities = self.getAbnormalRaysCells()
        path = r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\abnormal_rays.dat"
        file = open(path,'w')
        for abnormalRay in abnormalRays:
            for datum in abnormalRay:
                file.write(f"{datum} ")
            file.write("\n")
        file.close()
        file = open(r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\abnormal_res.den",'w')
        for cellDensity in cellDensities:
            file.write(f"{cellDensity}\n")
        file.close()
