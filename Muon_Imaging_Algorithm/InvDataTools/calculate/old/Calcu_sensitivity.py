# 如果对数据要求精确的话，使用Decimal
# from decimal import Decimal
from math import pi
from math import sin, cos

from prettytable.prettytable import PrettyTable
from tqdm import tqdm

from InvDataTools.calculate.cellbnd import Cellbnd
from InvDataTools.calculate.getmesh import Getmesh


class Calcsensitivity:
    def __init__(self, obsf, meshf):
        """

        :param obsf: obs文件的路径
        :param meshf: mesh文件的路径
        """
        self.obs_abspath = obsf
        self.mesh_abspath = meshf
        self.cell_entity = None
        self.x_node = None
        self.y_node = None
        self.z_node = None
        self.mx = None
        self.my = None
        self.mz = None
        # 所有的异常射线的数据,包含三项[射线编号,等效长度,等效长度误差]
        self.h_list = []
        self.x_list = []
        self.y_list = []
        # 每个探测器的射线未被使用数
        self.h_unused_ray = 0
        self.x_unused_ray = 0
        self.y_unused_ray = 0
        # 总射线数
        self.n_dat = 0
        # 总的射线未被使用数
        # self.unused_ray = 0
        # 每个探测器总的射线数
        self.det_ray = []
        # z方向的最大值
        self.elev0 = 0
        # 所有被使用的射线穿过的格子
        self.used_ray_cross_cell = []
        # 每个探测器的被使用射线
        self.det_used_ray = [[]]
        # 读取obs文件
        self.file_data = None
        # 输出
        self.loger = None
        # 准备工作
        self.pre_work()

    def pre_work(self):
        """
        提前准备的工作

        """
        # 得到各数值
        self.mx, self.my, self.mz, self.x_node, self.y_node, self.z_node, self.elev0 = Getmesh(
            meshf=self.mesh_abspath).getmesh()
        # 实例化Cellbnd类，其中的getcellbnd方法是获取探测器的离散坐标
        self.cell_entity = Cellbnd(self.mx, self.my, self.mz, self.x_node, self.y_node, self.z_node)
        # 读obs数据,open在读取时，不指定newline,则默认开启Universal new line mode,所有的\n,\r,or \r\n被默认转换为\n,因此判断空行只需和\n比较,
        # 或者也可以用strip()和''比较
        with open(self.obs_abspath, 'r') as r:
            # 第一列是describe信息
            file_data = r.readlines()[1:]
        while file_data[-1] == '\n':
            file_data = file_data[:-1]
        self.file_data = file_data
        self.n_dat = len(self.file_data)

    def calc_single_ray(self, y0, x0, z0, theta0, phi0, need_consider_xy_ray, need_record_unused_ray=True, args=()):
        """
        !此方法会因为参数导致假执行 弃用
        计算单条射线

        为方便表述，下文中的xy方向均指x正方向、x负方向、y正方向、y负方向;x方向指x正方向、x负方向;y方向指y正方向、y负方向

        need_consider_xy_ray(是否考虑xy方向)和need_record_unused_ray(是否记录未使用的射线)两者之间的关系为：

        1 当need_consider_xy_ray=True时:
        1.1 need_record_unused_ray=True,记录近乎水平的射线
        1.2 need_record_unused_ray=False,不记录
        2 当need_consider_xy_ray=False时:
        2.1 need_record_unused_ray=True,记录近乎水平的射线、x方向的射线、y方向的射线
        2.2 need_record_unused_ray=False,不记录

        :param y0: 探测器的y坐标
        :param x0: 探测器的x坐标
        :param z0: 探测器的z坐标
        :param theta0: 射线的角度(和水平面xoy的夹角)
        :param phi0: 射线的角度(和y轴的夹角,因为相当于已经投影到水平面xoy了)
        :param need_record_unused_ray: 是否记录未使用的射线
        :param need_consider_xy_ray: 从xy方向穿出的射线是否需要考虑,为true时考虑,即结果包含此类射线,false不包含
        :param args:
        :return: [该条射线是否可用(False:直接进行下一条,True:该条可用)，射线穿过哪些格子，射线穿过格子的实际长度]
        """
        if need_record_unused_ray:
            if len(args) != 3:
                raise ValueError("args参数的长度不等于3,(射线编号,等效长度,等效长度的误差)")

        xnode = self.x_node
        ynode = self.y_node
        znode = self.z_node
        mx = self.mx
        my = self.my
        mz = self.mz

        ind, cross_cell_length = [], []  # ind存放射线穿过哪些格子，即j，sensrow存放射线穿过格子时，在格子内的路径长度

        # ux是x轴的，uy是y轴的，uz是z轴的，单位长度(可以看作是1)投影到各个轴上的长度，只需要乘以实际长度就可以得到实际的投影长度
        ux = abs(sin(phi0) * cos(theta0))
        uy = abs(cos(phi0) * cos(theta0))
        uz = sin(theta0)

        dirr = [0] * 2
        # dirr[0]是x轴的，dirr[1]是y轴的
        dirr[0] = -1 if sin(phi0) < 0 else 1
        dirr[1] = -1 if cos(phi0) < 0 else 1

        # 考虑射线近乎(贴近)水平时的情况，此类射线不参与计算
        if abs(uz) < 1e-5:
            if need_record_unused_ray:
                self.h_unused_ray += 1
                self.h_list.append(args)
            return False, [], []

        # 得到探测器的位置离散坐标
        xc, yc, zc = self.cell_entity.getcellbnd(x0, y0, z0)
        # 是否停止(结束)计算,done为True时继续计算,为False时停止(结束)计算
        done = True
        # 在不考虑xy方向的射线的情况下,当前射线是否应当被使用,is_next为True被使用,为False时不被使用
        is_next = True
        # 射线在每个格子停留时，必定是停留在格子的某个面上
        # 如果射线与z轴的夹角很小
        if abs(uz - 1) < 1e-5:
            # 如果探测器几乎紧挨着格子壁，x0和y0就减小一点
            if abs(x0 - xnode[xc + 1] < 1e-5): x0 -= (xnode[xc + 1] - xnode[xc]) / 10
            if abs(y0 - ynode[yc + 1] < 1e-5): y0 -= (ynode[yc + 1] - ynode[yc]) / 10

            while done:
                zbound = znode[zc]
                dr = z0 - zbound
                z1 = z0 - dr * uz
                xcel = xc * my * mz + yc * mz + zc + 1

                ind.append(xcel)
                cross_cell_length.append(dr)

                zc -= 1
                z0 = z1
                if zc == -1:
                    done = False
        # 探测器与z轴的夹角既不大也不小
        else:
            while done:
                # 找到射线此时所在的格子的z方向的上边界值
                zbound = znode[zc]
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
                xcel = xc * my * mz + yc * mz + zc + 1

                # 存放射线穿过的格子ind以及在格子内的长度sensrow,大小都是numz,即射线穿过的格子数量
                ind.append(xcel)
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
                    if zc == -1: done = False
                elif xc == -1 or xc == mx:
                    if need_record_unused_ray and not need_consider_xy_ray:
                        self.x_list.append(args)
                        self.x_unused_ray += 1
                    done = is_next = False
                else:
                    if need_record_unused_ray and not need_consider_xy_ray:
                        self.y_list.append(args)
                        self.y_unused_ray += 1
                    done = is_next = False

        if not need_consider_xy_ray and not is_next:
            return False, [], []

        return True, ind, cross_cell_length

    def calc_single_ray_from_number(self, number, need_consider_xy_ray=True):
        """
        obs文件的第几行，即第几条射线，从1开始

        :param number: 射线编号
        :param need_consider_xy_ray: 同calc_single_ray
        """
        data = self.file_data[number - 1].split()
        result = self.calc_single_ray(float(data[1]), float(data[2]), self.elev0 - float(data[3]),
                                      pi / 2 - float(data[4]), float(data[5]),
                                      need_consider_xy_ray)
        if result[0]:
            print(f"第{number}条射线的数据：")
            print("穿过的格子编号")
            for ind in result[1]:
                print(f"{ind} ", end='')
            print()
            print("对应的射线穿过格子的长度")
            for cross_cell_length in result[2]:
                print(f"{cross_cell_length} ", end='')
            print()
        else:
            print("该条射线被舍弃，未被参与计算")

    def calc_single_or_multiple_det(self, det_des_list):
        """
        单个或多个探测器

        :param det_des_list: 探测器列表,可单可多
        """
        pass

    def calcsensitivity(self, ijf, Ray_way_j_file, loger, need_consider_xy_ray=True, isprint=True):
        """
        计算射线穿过格子的长度ijg

        :param ijf: ij文件的路径，当ijf为None时，默认会存放在obsf所在的文件夹下的ij文件中，i和j均是从1开始的
        :param Ray_way_j_file: 正常射线穿过的格子编号
        :param loger: 输出的格式
        :param need_consider_xy_ray: 同calc_single_ray
        :param isprint: 是否打印
        """
        self.loger = loger
        # 第一个探测器的描述信息
        old_det_des = self.file_data[0].split()[0]
        # 初始索引值
        old_ray_ind = 1
        ray_ind = 1

        with open(ijf, 'w') as w:

            w.write(f'{self.n_dat} {self.mx * self.my * self.mz}\n')

            for data in tqdm(self.file_data):

                data_list = data.split()
                det_des, y0, x0, z0, theta0, phi0, d, d_err = data_list[0], float(data_list[1]), float(
                    data_list[2]), self.elev0 - float(data_list[3]), pi / 2 - float(data_list[4]), float(
                    data_list[5]), float(data_list[7]), float(data_list[8])

                # 另一种进度条形式，简单的
                # 最多输出方框的个数
                # num = n_dat // 50
                # print(f"\r进度:{ray_ind / n_dat * 100:3.0f}%: {'▓' * (ray_ind // num)}", end="")

                if isprint:
                    # 记录除最后一个探测器以外的其余探测器的异常射线的信息(最后一个探测器的异常射线信息在循环结束后记录)
                    if det_des != old_det_des:
                        old_det_des = det_des
                        self.det_ray.append(ray_ind - old_ray_ind)
                        self.recover()
                        old_ray_ind = ray_ind
                        # 到达下一个探测器,添加新的列表用来存放下一个探测器的数据
                        self.det_used_ray.append([])

                result = self.calc_single_ray(y0, x0, z0, theta0, phi0, need_consider_xy_ray,
                                              need_record_unused_ray=isprint, args=(ray_ind, d, d_err))

                if not result[0]:
                    ray_ind += 1
                    continue

                cell_ind, cross_cell_length = result[1], result[2]
                for jg, cell_length in zip(cell_ind, cross_cell_length):
                    w.write(f'{ray_ind} {jg} {cell_length}\n')
                # 被使用的射线穿过的体素
                self.used_ray_cross_cell.append([])
                for jg in cell_ind:
                    self.used_ray_cross_cell[-1].append(jg)

                if isprint:
                    self.det_used_ray[-1].append(ray_ind)

                ray_ind += 1

        self.storage_used_ray(Ray_way_j_file)
        if isprint:
            self.det_ray.append(ray_ind - old_ray_ind)
            self.recover()
            self.print_data(loger)

    def recover(self):
        self.h_unused_ray = 0
        self.x_unused_ray = 0
        self.y_unused_ray = 0

    def storage_used_ray(self, used_abspath):
        """
        正常射线穿过的格子编号

        :param used_abspath:
        :return:
        """
        with open(used_abspath, 'w') as w:
            for val in self.used_ray_cross_cell:
                for write_val in val:
                    w.write(f'{write_val} ')
                w.write('\n')

    def print_data(self, loger):
        try:
            self.loger_data(loger)
        except Exception as e:
            print(e)

    def loger_data(self, loger):
        # 存储在不同条件下排除射线的编号

        loger.write("未被使用的射线")
        loger.write('接近水平的射线(无论如何都不会被使用): ')
        h_str = ' '.join(str(h) for h in self.h_list)
        loger.write(h_str)
        loger.write('从x(正、负)方向射出的射线: ')
        x_str = ' '.join(str(x[0]) for x in self.x_list)
        loger.write(x_str)
        loger.write('从y(正、负)方向射出的射线: ')
        y_str = ' '.join(str(y) for y in self.y_list)
        loger.write(y_str)
        loger.write('未被使用的射线的具体信息,分别是"射线编号,等效长度,等效长度误差"')
        # 存储在不同条件下未被使用射线的具体信息(射线编号,等效长度,等效长度误差)
        for h in self.h_list:
            loger.write(f'{h[0]} {h[1]} {h[2]}')
        for x in self.x_list:
            loger.write(f'{x[0]} {x[1]} {x[2]}')
        for y in self.y_list:
            loger.write(f'{y[0]} {y[1]} {y[2]}')
        # 存储每个探测器的正常射线总数以及具体编号
        loger.write("每个探测器的被使用射线总数以及具体编号:")
        for used_ray in self.det_used_ray:
            loger.write(len(used_ray))
            ray_str = ' '.join(str(ray) for ray in used_ray)
            loger.write(ray_str)
        # 总计
        unused_ray_total = self.n_dat - len(self.used_ray_cross_cell)
        loger.write(f'共{self.n_dat}条射线------共{len(self.used_ray_cross_cell)}条射线被使用------共{unused_ray_total}条射线未被使用')
        loger.write(
            f'在未被使用的射线中，接近水平的有{len(self.h_list)}条，从x(正、负)方向射出的有{len(self.x_list)}条，从y(正、负)方向射出的有{len(self.y_list)}条')

        # 中文会引起对不齐，改成英文可以解决
        # tabulate也是一个不错的以表格形式展示数据的包
        table = PrettyTable(["探测器编号", "射线未被使用数", "射线总数", "未被使用百分比(未被使用数/总数)"])
        for no, ray in enumerate(self.det_ray):
            unused = ray - len(self.det_used_ray[no])
            table.add_row([no + 1, unused, ray, '{:.2%}'.format(unused / ray)])
        table.add_row(['总计', unused_ray_total, self.n_dat, '{:.2%}'.format(unused_ray_total / self.n_dat)])
        loger.write(table)
