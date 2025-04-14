class Cellbnd:

    def __init__(self, mx, my, mz, xnode, ynode, znode):
        self.mx = mx
        self.my = my
        self.mz = mz
        self.xnode = xnode
        self.ynode = ynode
        self.znode = znode

    def getcellbnd(self, x0, y0, z0):
        """
        :param mx: 在x方向的格子数
        :param my: 在y方向的格子数
        :param mz: 在z方向的格子数
        :param x0: 探测器的x坐标
        :param y0: 探测器的y坐标
        :param z0: 探测器在z轴方向到最高点的距离
        :param xnode:
        :param ynode:
        :param znode:
        :return: xc为x0最近的且比x0小的点，yc,zc同理，也即探测器的离散坐标
        """

        # yc的起始值是-1，表示探测器位于模型最前边的外面，yc=0，表示探测器位于y方向第0个格子内，yc=my-1，表示探测器位于y方向最后一个格子内，yc=my，表示探测器位于模型最下边的外面
        yc = -1
        for i in range(self.my + 1):
            if y0 > self.ynode[i]:
                yc += 1
            else:
                break

        # xc的起始值是-1，表示探测器位于模型最左边的外面，xc=0，表示探测器位于x方向第0个格子内，xc=mx-1，表示探测器位于x方向最后一个格子内，xc=mx，表示探测器位于模型最右边的外面
        xc = -1
        for i in range(self.mx + 1):
            if x0 > self.xnode[i]:
                xc += 1
            else:
                break

        # zc的起始值是-1，表示探测器位于模型最上边的外面，zc=0，表示探测器位于z方向第0个格子内，zc=mz-1，表示探测器位于z方向最后一个格子内，zc=mz，表示探测器位于模型最下边的外面
        zc = -1
        for i in range(self.mz + 1):
            if z0 > self.znode[i]:
                zc += 1
            else:
                break

        if xc == -1 or yc == -1 or zc == -1 or xc == self.mx or yc == self.my or zc == self.mz:
            raise ValueError("探测器不在mesh边界之内")

        return xc, yc, zc