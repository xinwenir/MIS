# author:高金磊
# datetime:2022/3/14 15:09


j_xyz_cache = {}


def getxyz_from_shape(shape, j):
    """
    代替jxyz文件的一种补充方案，之后可能代替jxyz文件的读写
    起始坐标是（1，1，1）与jxyz保持一致

    :param shape: 模型在x、y、z方向的格子数目
    :param j: range(1,shape[0]*shape[1]*shape[2]),如果是系统的j也就是从0开始的必须加1再传进来,比如newj_oldj
    :return: 格子编号对应的离散坐标
    """
    if j in j_xyz_cache.keys():
        return j_xyz_cache[j]
    j -= 1
    y = int(j / (shape[0] * shape[2]))
    j = j % (shape[0] * shape[2])
    x = int(j / shape[2])
    z = j % shape[2]
    xyz = [x + 1, y + 1, shape[2] - z]
    j_xyz_cache[j + 1] = xyz
    return xyz


def getj_from_xyz(shape, xyz):
    """
    根据格子的离散坐标x、y、z得到格子的编号

    :param shape: 模型在x、y、z方向上的格子数目
    :param xyz: 格子的离散坐标
    :return: 格子的编号 j range(1,shape[0]*shape[1]*shape[2])
    """
    x, y, z = xyz
    return shape[0] * shape[2] * (y - 1) + shape[2] * (x - 1) + (shape[2] - z + 1)


class Make_jxyz:
    """
    生成jxyz文件

    """
    def jxyz_form_shape(self, shape, jxyz_file):
        """
        根据shape自动生成jxyz文件,本文件的j包含空气，即模型分割后的所有格子

        :param shape: 模型在x、y、z方向上的格子数目 [xnum,ynum,znum]
        :param jxyz_file: 结果存放的文件路径
        """
        jxyz = open(jxyz_file, 'w')
        for x in range(shape[0]):
            for y in range(shape[1]):
                for z in range(shape[2]):
                    jxyz.write(
                        str(shape[0] * shape[2] * (y) + shape[2] * (x) + (shape[2] - z)) + " " + str(x + 1) + " " + str(
                            y + 1) + " " + str(z + 1) + "\n")
        # for y in range(shape[1]):
        #     for x in range(shape[0]):
        #         for z in range(shape[2]):
        #             jxyz.write(str(j) + " " + str(x + 1) + " " + str(y + 1) + " " + str(z + 1) + "\n")
        #             j += 1
        jxyz.flush()
        jxyz.close()

# if __name__ == '__main__':
# Make_jxyz().get_g()
# Make_jxyz().jxyz_form_shape([140,72,42])
# tool = show_ray_trace_tools([280, 144, 63])
# target=[i for i in range(26612,26621)]+[i  for i in range(1240,1249) ]
# target=[i for i in range(6714)]
# target=[100,300,555]
# target = [i for i in range(25770, 40000)]
# tool.mark_ray(target)
# InvDataTools.show_ray_trace_tools.show_ray_trace_tools([280, 144, 63]).mark_ray_all()

# tool.mark_ray([36929])
