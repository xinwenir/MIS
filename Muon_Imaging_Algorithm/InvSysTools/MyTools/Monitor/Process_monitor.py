# author:高金磊
# datetime:2022/8/3 14:59
import datetime
import os
import time

import psutil
import threading

class _info_struct:
    def __init__(self,p):
        self.info={
            "PID":p.pid,
            "进程创建时间":p.create_time(),
            "cpu占用百分比":p.cpu_percent()/len(p.cpu_affinity()),#####注意这里仅仅适用于python
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
    time.sleep(1)  # 给子线程一定时间以启动
def timing_run(interval,p,file_obj):
    infos=[]
    cpu_infos=[0.0 for i in range(60)]
    mem_infos=[0.0 for i in range(60)]
    t_infos=[-i for i in range(60)]
    t_infos.reverse()
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)
    ax1 = ax.twinx()
    plt.ion()
    all_mem=psutil.virtual_memory().total/1024/1024
    cpu_max=0
    mem_max=0
    while 1:
        # time.sleep(interval)
        info=_info_struct(p)
        file_obj.write(str(info))
        file_obj.write("\n")
        infos.append(info)
        # print(info)
        cpu_infos.append(float(info.getCPU()))
        mem_infos.append(float(info.getMEM())*all_mem/100)
        t_infos.append(t_infos[-1]+interval)
        cpu_max=max(cpu_infos[-1],cpu_max)
        mem_max=max(mem_max,mem_infos[-1])
        ax.cla()
        ax1.cla()
        ax.set_ylim(bottom=0,top=(cpu_max*1.1))
        ax1.set_ylim(bottom=0,top=mem_max*1.1)
        ax.plot(t_infos, cpu_infos, c='r', ls='-' ,label="CPU: %.2f %%"%(cpu_infos[-1]))  ## 保存历史数据
        ax.set_title('cpu_avg: %.2f %%  memory_avg: %.2f MB \n cpu_max: %.2f %% mem_max: %.2f MB' % (
        sum(cpu_infos) / (len(cpu_infos) - mem_infos.count(0.0)),
        sum(mem_infos) / (len(cpu_infos) - mem_infos.count(0.0)), cpu_max, mem_max), fontsize=12, color='b')
        ax.set_title('cpu_avg: %.2f %%  memory_avg: %.2f MB \n cpu_max: %.2f %% mem_max: %.2f MB'%(sum(cpu_infos)/(len(cpu_infos)-mem_infos.count(0.0)),sum(mem_infos)/(len(cpu_infos)-mem_infos.count(0.0)),cpu_max,mem_max), fontsize=12, color='b')
        ax.legend(loc="upper right")
        ax1.plot(t_infos, mem_infos, c='g', ls='-',label="Mem: %.2f MB"%(mem_infos[-1]))
        ax1.legend(loc="upper left")
        if len(t_infos)>60:
            t_infos.pop(0)
            cpu_infos.pop(0)
            mem_infos.pop(0)
        plt.pause(interval=interval)



def test():
    process_monitoring(interval=2)
    middle=[]
    for i in range(1000):
        middle.append([i*"12345678"])
        time.sleep(0.2)
# test()
