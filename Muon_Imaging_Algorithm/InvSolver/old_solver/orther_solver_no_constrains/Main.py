# author:高金磊
# datetime:2021/11/21 9:49
from solver.old_solver.Setting import *
from scipy.sparse.linalg import lsqr

from DataTools.res_tools import res_tools

if __name__ == '__main__':
    print("开始计算.....%s" % (datetime.now()))
    res=lsqr(M,d,show=True,iter_lim=100)
    print(res)
    print("误差："+str(np.linalg.norm(M*res[0]-d)))
    print("计算完毕,正在读写文件...%s" % (datetime.now()))
    fo = open("../../DataTools/data/res", "w")
    for re in res[0]:
        fo.write(str(re))
        fo.write('\n')
    fo.close()
    print("正在重组结果...%s" % (datetime.now()))
    try:
        rt=res_tools()
        rt.remark_point(default_value=1000)
        rt.Conversion_2(shape=shape)
    except Exception:
        print("结果重组失败,请检查原因并尝试使用res_tools工具重组结果")
    print("程序运行完毕...%s" % (datetime.now()))