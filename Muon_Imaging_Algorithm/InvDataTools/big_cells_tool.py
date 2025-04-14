class BigCellsTool:
    def __init__(self, path) -> None:
        file = open(path, 'r')
        file.readline()
        self.location = file.readline()
        data = []
        for line in file.readlines():
            data.append([float(i) for i in line.strip().split(" ")])
        self.old_data = data
        file.close()
        self.new_data = []
        self.new_shape = []

    def simple_merge(self, x_num, y_num, z_num):
        """
        xyz三个方向几个为一组合并
        """
        xyz_num = [x_num, y_num, z_num]
        for i in range(len(self.old_data)):
            cur_sum = 0
            cur_distances = []
            for j in range(len(self.old_data[i])):
                cur_sum += self.old_data[i][j]
                if (j+1) % xyz_num[i] == 0:
                    cur_distances.append(round(cur_sum, 6))
                    cur_sum = 0
            if cur_sum != 0:
                cur_distances.append(round(cur_sum, 6))
            self.new_data.append(cur_distances)
            self.new_shape.append(len(cur_distances))
                   
    def merge_by_percent(self,xyz_percent:list):
        """
        排序后根据百分之n的数为最小长度合并
        假如传1就是根据每行最大的长度为界限
        """
        max_lens = [sorted(self.old_data[i])[max(0, int(len(self.old_data[i]) * xyz_percent[i]) - 1)] for i in range(len(self.old_data))]
        for i in range(len(self.old_data)):
            cur_sum = 0
            cur_distances = []
            for length in self.old_data[i]:
                cur_sum += length
                if cur_sum >= max_lens[i]:
                    cur_distances.append(round(cur_sum, 6))
                    cur_sum = 0
            if cur_sum != 0:
                cur_distances.append(round(cur_sum, 6))
            self.new_data.append(cur_distances)
            self.new_shape.append(len(cur_distances))    

    def save_big_cells_msh(self):
        """
        存储新的bigmsh文件
        """
        if len(self.new_data) == 0:  # 假如没有merge就用maxMerge
            self.simple_merge([3, 3, 3])
        file = open(r"E:\vscode\Muon_Imaging_Algorithm\data\Temp\Merge.msh",'w')
        for number in self.new_shape:
            file.write(f"{number} ")
        file.write("\n")
        file.write(self.location)
        for distances in self.new_data:
            for distance in distances:
                file.write(f"{distance} ")
            file.write("\n")    
        file.close()


if __name__ == '__main__':
    path = r"E:\vscode\Muon_Imaging_Algorithm\data\paper3\wall\17_58MaMian.msh"
    bigCellsTool = BigCellsTool(path)
    # bigCellsTool.simple_merge([5, 3, 3])
    bigCellsTool.merge_by_percent([1, 1, 1])
    bigCellsTool.save_big_cells_msh()
    



