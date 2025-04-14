# author:高金磊
# datetime:2022/3/29 9:28
"""不可依赖任何存在的包"""
# import sys
# sys.path.append(r"m:/pycharm/Inversion")

from copy import copy


class MeshTools:
    """
    读取并处理”mesh“文件

    """

    def __init__(self, mesh_file):
        file = open(mesh_file, 'r')
        # 第一行
        line = file.readline()
        middle = line.split()
        self._shape = [int(middle[0]),int(middle[1]),  int(middle[2])]
        # 第二行
        line = file.readline()
        middle = line.split()
        self.Start_x = float(middle[0])
        self.Start_y = float(middle[1])
        self.end_z = float(middle[2])
        # 其他行
        line = file.readline()
        middle = line.split()
        self._cell_xs = [float(i) for i in middle]
        line = file.readline()
        middle = line.split()
        self._cell_ys = [float(i) for i in middle]
        line = file.readline()
        middle = line.split()
        self._cell_zs = [float(i) for i in middle]
        self._cell_zs.reverse()
        file.close()

    def get_shape(self):
        """
        获取x、y、z方向上的格子数目

        :return: x、y、z方向上的格子数目
        """
        return copy(self._shape)

    def cells_count(self):
        """
        获取模型被分割后的格子数目

        :return: 模型被分割后的格子数目
        """
        return self._shape[0] * self._shape[1] * self._shape[2]

    def get_xs(self):
        """
        获取x方向上的每个格子的长度

        :return: x方向上的每个格子的长度
        """
        return copy(self._cell_xs)

    _x_start_values = None

    @property
    def x_start_values(self):
        """
        获取每个格子在x方向上的起始值

        :return: 每个格子在x方向上的起始值
        """
        if self._x_start_values is None:
            self._x_start_values = []
            self._x_start_values.append(0 + self.get_xyz_start()[0])
            for i in range(0, len(self._cell_xs) - 1):
                self._x_start_values.append(self._x_start_values[- 1] + self._cell_xs[i])
        return copy(self._x_start_values)

    _x_end_values = None

    @property
    def x_end_values(self):
        """
        获取每个格子在x方向上的结束值

        :return: 每个格子在x方向上的结束值
        """
        if self._x_end_values is None:
            self._x_end_values = []
            self._x_end_values.append(self._cell_xs[0] + self.get_xyz_start()[0])
            for i in range(1, len(self._cell_xs)):
                self._x_end_values.append(self._x_end_values[-1] + self._cell_xs[i])
        return copy(self._x_end_values)

    _y_start_values = None

    @property
    def y_start_values(self):
        """
        获取每个格子在y方向上的起始值

        :return: 每个格子在y方向上的起始值
        """
        if self._y_start_values is None:
            self._y_start_values = []
            self._y_start_values.append(0 + self.get_xyz_start()[1])
            for i in range(0, len(self._cell_ys) - 1):
                self._y_start_values.append(self._y_start_values[- 1] + self._cell_ys[i])
        return copy(self._y_start_values)

    _y_end_values = None

    @property
    def y_end_values(self):
        """
        获取每个格子在y方向上的结束值

        :return: 每个格子在y方向上的结束值
        """
        if self._y_end_values is None:
            self._y_end_values = []
            self._y_end_values.append(self._cell_ys[0] + self.get_xyz_start()[1])
            for i in range(1, len(self._cell_ys)):
                self._y_end_values.append(self._y_end_values[- 1] + self._cell_ys[i])
        return copy(self._y_end_values)

    _z_start_values = None

    @property
    def z_start_values(self):
        """
        获取每个格子在z方向上的起始值

        :return: 每个格子在z方向上的起始值
        """
        if self._z_start_values is None:
            self._z_start_values = []
            self._z_start_values.append(0 + self.get_xyz_start()[2])
            for i in range(0, len(self._cell_zs) - 1):
                self._z_start_values.append(self._z_start_values[- 1] + self._cell_zs[i])
            # self.z_start_values.reverse()  # 注意Z是从上到下排的
        return copy(self._z_start_values)

    _z_end_values = None

    @property
    def z_end_values(self):
        """
        获取每个格子在z方向上的结束值

        :return: 每个格子在z方向上的结束值
        """
        if self._z_end_values is None:
            self._z_end_values = []
            self._z_end_values.append(self._cell_zs[0] + self.get_xyz_start()[2])
            for i in range(1, len(self._cell_zs)):
                self._z_end_values.append(self._z_end_values[- 1] + self._cell_zs[i])
            # self.z_end_values.reverse()#注意Z是从上到下排的
        return copy(self._z_end_values)

    def get_ys(self):
        """
        获取y方向上的每个格子的长度

        :return: y方向上的每个格子的长度
        """
        return copy(self._cell_ys)

    def get_zs(self):
        """
        获取z方向上的每个格子的长度

        :return: z方向上的每个格子的长度
        """
        return copy(self._cell_zs)

    def get_xyz_start(self):
        """
        获取模型在x、y、z方向的起始值

        :return: 模型在x、y、z方向的起始值
        """
        return [self.Start_x, self.Start_y, self.end_z - sum(self.get_zs())]

    def get_xyz_end(self):
        """
        获取模型在x、y、z方向的结束值

        :return: 模型在x、y、z方向的结束值
        """
        return [self.Start_x + sum(self.get_xs()), self.Start_y + sum(self.get_ys()), self.end_z]
    

    _strategy_set = {0, 1, 2}

    def get_coordinates_form_xyz(self, x, y, z, strategy=(0, 0, 0)):
        """
        根据给定的策略来确定用格子的哪个坐标(连续意义上的坐标)来代表格子本身，根据策略不同，坐标可分为27种，
        格子的坐标由3个值确定，每个值有3种取法，因此一共有3x3x3=27中取法。
        3种取法为：格子在某个方向的起始点，格子在某个方向的结束点，(格子在某个方向的起始点+格子在某个方向的结束点)/2

        :param x: 格子在x方向的离散点 min=1
        :param y: 格子在y方向的离散点
        :param z: 格子在z方向的离散点
        :param strategy: 0 中点,1 左边/下边/前边 ,2 右边/上边/后边
        :return: 用来代表格子的连续坐标
        """
        point = []
        x -= 1
        y -= 1
        # if z > 55:
        #     print(z)
        z -= 1  # z是倒着来的

        # """"""""""""计算用于代表格子位置的x坐标"""""""""""""""""#
        x_start = self.x_start_values[x]
        x_end = self.x_end_values[x]

        if strategy[0] == 0:
            point.append(0.5 * (x_start + x_end))

        elif strategy[0] == 1:
            point.append(x_start)

        elif strategy[0] == 2:
            point.append(x_end)
        else:
            raise Exception("当前策略不合法参考值:%s" % (self._strategy_set))

        # """"""""""""计算用于代表格子位置的y坐标"""""""""""""""""#
        y_start = self.y_start_values[y]
        y_end = self.y_end_values[y]
        if strategy[1] == 0:
            point.append(0.5 * (y_end + y_start))
        elif strategy[1] == 1:
            point.append(y_start)
        elif strategy[1] == 2:
            point.append(y_end)
        else:
            raise Exception("当前策略不合法参考值:%s" % (self._strategy_set))

        # """"""""""""计算用于代表格子位置的z坐标"""""""""""""""""#
        # try:
        #     z_start = self.z_start_values[z]
        # except:
        #     print(z)
        z_start = self.z_start_values[z]
        z_end = self.z_end_values[z]

        if strategy[2] == 0:
            point.append(0.5 * (z_end + z_start))

        elif strategy[2] == 1:
            point.append(z_start)

        elif strategy[2] == 2:
            point.append(z_end)
        else:
            raise Exception("当前策略不合法参考值:%s" % (self._strategy_set))

        return point

    def discretize_Physical_coordinates_x(self, x):
        """
        离散x值

        :param x: 连续x值
        :return: x所属区间的起始值、x所属区间的结束值、x本身的值、x离散化后的值
        """
        x_start = self.x_start_values
        x_end = self.x_end_values
        for i in range(len(x_start)):
            if x <= x_end[i] and x >= x_start[i]:
                return x_start[i], x_end[i], x, i + 1
        # 查询失败
        return -1, -1, x, -1

    def discretize_Physical_coordinates_y(self, y):
        """
        离散y值

        :param y: 连续y值
        :return: y所属区间的起始值、y所属区间的结束值、y本身的值、y离散化后的值
        """
        y_start = self.y_start_values
        y_end = self.y_end_values
        for i in range(len(y_start)):
            if y <= y_end[i] and y >= y_start[i]:
                return y_start[i], y_end[i], y, i + 1
        # 查询失败
        return -1, -1, y, -1

    def discretize_Physical_coordinates_z(self, z):
        """
        离散z值

        :param z: 连续z值
        :return: z所属区间的起始值、z所属区间的结束值、z本身的值、z离散化后的值
        """
        z_start = self.z_start_values
        z_end = self.z_end_values
        for i in range(len(z_start)):
            if z <= z_end[i] and z >= z_start[i]:
                return z_start[i], z_end[i], z, i + 1
        # 查询失败
        return -1, -1, z, -1

    def discretize_Physical_coordinates(self, x, y, z, detailed=False):
        """
        输入连续形式的坐标,返回离散坐标

        :param x: 连续x值
        :param y: 连续y值
        :param z: 连续z值
        :param detailed: 是否显示详情
        :return: 离散坐标(xyz方向第几个,可以直接计算j的坐标)
        """
        x = self.discretize_Physical_coordinates_x(x)
        y = self.discretize_Physical_coordinates_y(y)
        z = self.discretize_Physical_coordinates_z(z)
        if detailed:
            print("x的判断结果:", x)
            print("y的判断结果:", y)
            print("z的判断结果:", z)
        return x[-1], y[-1], z[-1]

    @staticmethod
    def generate_mesh_file(origin_of_coordinates,x_directions,y_directions,z_directions,mesh_file_path):
        """根据参数生成mesh文件并保存到mesh_file_path位置
        #!本接口不做类型和数据检查

        Args:
            origin_of_coordinates (list[1,3]): mesh文件空间的坐标原点
            x_directions (list[n,2]):  x方向格子的长度2维数组,每一个的第一个值是格子长度float,第二个值是这种格子的个数int
            y_directions (list[n,2]): y方向格子的长度2维数组,每一个的第一个值是格子长度float,第二个值是这种格子的个数int
            z_directions (list[n,2]): z方向格子的长度2维数组,每一个的第一个值是格子长度float,第二个值是这种格子的个数int
            mesh_file_path (string): mesh文件的保存位置
        return:
            MeshTools
        Simple code:
        
        from InvDataTools.MeshTools import MeshTools
        origin_of_coordinates=[0,1,-10]
        x_directions=[[1,10],[1.1,10],[1,10],[1.0,10]]
        y_directions=[[1,10],[1.1,10],[1,10],[1.0,10]]
        z_directions=[[1,10],[1.1,10],[1,10],[1.0,10]]
        mesh_file_path=r"M:\pycharm\Inversion\Temp\mesh"
        
        res=MeshTools.generate_mesh_file(origin_of_coordinates,x_directions,y_directions,z_directions,mesh_file_path)
        """
        file_obj=open(mesh_file_path,'w')
        x_count=0
        y_count=0
        z_count=0
        x_start=origin_of_coordinates[0]
        y_start=origin_of_coordinates[1]
        z_start=origin_of_coordinates[2]#!最终这个要的是end
        x_strig=""
        y_strig=""
        z_strig=""
        for data in x_directions:
            value,number=data
            x_count+=number
            for i in range(number):
                if x_strig =="":
                    x_strig+=str(value)
                else:
                    x_strig+=' '
                    x_strig+=str(value)
        for data in y_directions:
            value,number=data
            y_count+=number
            for i in range(number):
                if y_strig =="":
                    y_strig+=str(value)
                else:
                    y_strig+=' '
                    y_strig+=str(value)
        for data in z_directions:
            value,number=data
            z_count+=number
            for i in range(number):
                z_start+=value
                if z_strig =="":
                    z_strig+=str(value)
                else:
                    z_strig+=' '
                    z_strig+=str(value)
        file_obj.write(str(x_count))
        file_obj.write(" ")
        file_obj.write(str(y_count))
        file_obj.write(" ")
        file_obj.write(str(z_count))
        file_obj.write("\n")
        file_obj.write(str(x_start))
        file_obj.write(" ")
        file_obj.write(str(y_start))
        file_obj.write(" ")
        file_obj.write(str(z_start))
        file_obj.write("\n")
        file_obj.write(x_strig)
        file_obj.write("\n")
        file_obj.write(y_strig)
        file_obj.write("\n")
        file_obj.write(z_strig)
        file_obj.write("\n")
        file_obj.close()
        return MeshTools(mesh_file_path)

# %%
# if __name__ == '__main__':
#     mesh_tool = MeshTools(r"E:\vscode\Muon_Imaging_Algorithm\data\Input\real_data\45_58MaMian.msh")
#     print(mesh_tool.get_xs())
#     print(mesh_tool.get_ys())
#     print(mesh_tool.get_zs())
#     print(mesh_tool.get_shape())
#     print(mesh_tool.get_xyz_start())
# %%
