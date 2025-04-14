# author:高金磊
# datetime:2022/8/3 21:58
"""
远程性能监控,需要提前开启服务器

实用性较差-----暂停维护
"""
import datetime
import os
import time

import psutil
import threading

from MyTools.Socket_Tools.sockets import SocketServer, recver_message, SocketClient, send_message
from MyTools.myPrint import myPrint_Success


class _info_struct:
    def __init__(self,p):
        self.info={
            "PID":p.pid,
            "进程创建时间":p.create_time(),
            "cpu占用百分比":p.cpu_percent(),
            # "cpu占用核心数": p.cpu_num(),
            "内存占用百分比":p.memory_percent(),
            "本条记录时间": str(datetime.datetime.now())
        }
    def __str__(self):
        info=""
        for key in self.info:
            info+=str(key)+": "+str(self.info[key])+" "
        return info
    def getPID(self):
        return self.info["PID"]
    def getCPU(self):
        return self.info["cpu占用百分比"]
    def getMEM(self):
        return self.info["内存占用百分比"]

def _drow_cpu_mem(socket):
    import matplotlib.pyplot as plt
    plt.ion()
    cpu_infos=[0]
    mem_infos=[0]
    t_infos=[0]
    while 1:
        data= recver_message(socket)
        myPrint_Success("收到:", data)
        data=data.split(' ')
        cpu_infos.append(float(data[0]) * 100)
        mem_infos.append(float(data[1]) * 100)
        t_infos.append(t_infos[-1] + float(data[2]))
        plt.plot(t_infos, cpu_infos, c='r', ls='-', marker='o', mec='b', mfc='w')  ## 保存历史数据
        plt.pause(0.1)
def start_sever(is_multithreading=False):
    SocketServer().startup(fun=_drow_cpu_mem,is_multithreading=is_multithreading)
def process_monitoring(pid=os.getpid(),file="cache",interval=10):
    p = psutil.Process(pid)
    file_obj=open(file,'a')
    file_obj.write(str(datetime.datetime.now()))
    file_obj.write(' ')
    file_obj.write(str(pid))
    file_obj.write(' ')
    file_obj.write(str(p.name()))
    file_obj.write(' ')
    threading.Thread(target=timing_run,args=(interval,p,file_obj)).start()

def timing_run(interval,p,file_obj):
    infos=[]
    def fun(socket):
        while 1:
            time.sleep(interval)
            info=_info_struct(p)
            file_obj.write(str(info))
            file_obj.write("\n")
            infos.append(info)
            message=str(info.getCPU())+" "+str(info.getMEM())+" "+str(interval)
            print(message)
            send_message(socket=socket,data=message)
    SocketClient().connect(fun=fun)
        # print(info)




