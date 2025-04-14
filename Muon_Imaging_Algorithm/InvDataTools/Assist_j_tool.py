# author:高金磊
# datetime:2022/5/13 18:46

class Assist_j_tool:
    """
    将格子编号全部缩小某个倍数，缩小后的值小于1

    """
    @classmethod
    def make_assist_j_file(cls,shape,path):
        """
        将格子编号全部缩小某个倍数，缩小后的值小于1

        example: 格子数为n，则第i(格子编号从1开始)个格子经过处理后的结果为i/10**len(str(n))

        :param shape: 模型在x、y、z方向上的格子数
        :param path: 存放处理后的格子编号的文件路径
        """
        file=open(path,'w')
        n=shape[0]*shape[1]*shape[2]
        offset=10**(len(str(n)))
        for i in range(n):
            file.write(str((i+1)/offset))
            file.write('\n')
        file.close()

