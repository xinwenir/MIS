# author:高金磊
# datetime:2022/6/12 10:13

class Cell_Weight_Tools():
    def __init__(self,path):
        file=open(path,'r')
        lines=file.readlines()
        data=[]
        for line in lines:
            if line =="":
                break
            datum=[]
            for num in line.split():
                datum.append(int(num))
            data.append(datum)
        file.close()
        self.data=data
    def get_cell_weight_by_rays(self,amount=5):
        """
        收集射线附近的amount个格子的编号
        :param amount:
        :return:格子编号,从1开始
        """
        cell_ids=set()
        for datum in self.data:
            for i in range(1,min(len(datum),amount+1)):
                cell_ids.add(datum[i])
        return cell_ids

    def get_cell_weight_by_detector(self, amount=30):
        pass


if __name__ == '__main__':
    tool=Cell_Weight_Tools(path=r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\Ray_way_j")
    res=tool.get_cell_weight_by_rays(5)
