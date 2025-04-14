# author:高金磊
# datetime:2022/7/29 21:58
from pyinstrument import Profiler

def code_performance_show(test_fun, args=None):
    """
    测试代码性能
    :param test_fun:需要测试的方法名
    :param args: 方法对应的参数,不支持指定参数名和值,为元组
    :return: None,自动打印结果
    """
    profiler = Profiler()
    profiler.start()
    try:
        if args is None:
            test_fun()
        else:
            test_fun(*args)
    except Exception as e:
        print(e)

    # 要分析的代码
    profiler.stop()

    profiler.print()