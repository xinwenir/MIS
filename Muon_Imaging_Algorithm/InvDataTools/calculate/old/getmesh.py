class Getmesh:

    def __init__(self, meshf):
        """
        得到格子边界的值

        :param meshf: 存放mesh文件的路径
        :return:
        """
        with open(meshf, 'r') as f:
            # mx,my,mz为模型在x,y,z方向被划分成的格子数量
            my, mx, mz = [int(val) for val in f.readline().strip().split()]

            # xmod,ymod,zmod的长度在划分格子的基础上加一,以xmod为例，第一个值存放x方向的最小值，后面的值存放沿x方向被划分的格子的长度，ymod和zmod同理
            xmod = [0.0] * (mx + 1)
            ymod = [0.0] * (my + 1)
            zmod = [0.0] * (mz + 1)
            # ymod[0], xmod[0], elev0 分别为模型的y方向最小值，x方向最小值，z方向最大值
            ymod[0], xmod[0], elev0 = [float(val) for val in f.readline().strip().split()]
            # 读取y,x,z方向的格子大小
            yy = f.readline().strip()
            xx = f.readline().strip()
            zz = f.readline().strip()
        # y方向均分
        if '*' in yy:
            num, size = yy.split('*')
            ymod[1:] = [float(size)] * int(num)
        # y方向不均分
        else:
            y_line=yy.split()
            ymod[1:] = [float(val) for val in y_line]

        # x方向均分
        if '*' in xx:
            num, size = xx.split('*')
            xmod[1:] = [float(size)] * int(num)
        # x方向不均分
        else:
            x_line=xx.split()
            xmod[1:] = [float(val) for val in x_line]

        # z方向均分
        if '*' in zz:
            num, size = zz.split('*')
            zmod[1:] = [float(size)] * int(num)
        # z方向不均分
        else:
            z_line=zz.split()
            zmod[1:] = [float(val) for val in z_line]
        self.my = my
        self.mx = mx
        self.mz = mz
        self.elev0 = elev0

        self.calc_cell_boundary_value(xmod,ymod,zmod)

    def calc_cell_boundary_value(self,xmod,ymod,zmod):

        # 除去第一个值外，后续的值更新(自身的值加上前一个的值)，也就是说，除去最后一个格子，ymod[i]的值是第i个格子的左边界值；
        # 假想最后面有一个哑格子，ymod[my]就是哑格子的左边界值;
        # 获取第i个格子的左右边值的方法：左边界->ymod[i],右边界->ymod[i+1].
        for i in range(1, self.my + 1):
            ymod[i] = ymod[i - 1] + ymod[i]

        for i in range(1, self.mx + 1):
            xmod[i] = xmod[i - 1] + xmod[i]

        # zmod[0]不需要声明为0了，因为定义的时候已经初始化为0了
        # 注意zmod是从最顶端开始计算的
        for i in range(1, self.mz + 1):
            zmod[i] = zmod[i - 1] + zmod[i]

        self.xmod = xmod
        self.ymod = ymod
        self.zmod = zmod

    def getmesh(self):
        return self.mx, self.my, self.mz, self.xmod, self.ymod, self.zmod, self.elev0