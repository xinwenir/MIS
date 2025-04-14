# author:高金磊
# datetime:2022/8/1 8:56

import psutil
import os

def Chain_Process_tool(command:str,funs):
    """
    查询当前进程的情况,暂时不支持查询其他进程
    :param command:  命令"show process"
    :param fun: 未命中后执行的操作
    :return:
    """

    funs.pop(0)
    if command=="-h":
        return ["","show process\n"+str(funs[0](command,funs)[1])]
    if "show process"==command:
        return ["Chain_Process_tool",_get_Process_cpu_mem()]
    else:
        return funs[0](command,funs)

def _get_Process_cpu_mem(pid=None):
    res=""
    try:
        if pid is None:
            pid = os.getpid()
        p=psutil.Process(pid)
        res+="PID:%s 进程名称%s: 进程创建时间:%s cpu占用百分比:%s 内存占用百分比:%s"%(p.pid,p.name(),p.create_time(),p.cpu_percent(),p.memory_percent())
    except Exception as e:
        res+="获取失败将获取本进程的信息\n"+str(e)+"\n"+_get_Process_cpu_mem()

    return res