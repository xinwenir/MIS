# author:高金磊
# datetime:2022/8/3 22:04
import time

from InvSysTools.MyTools.Monitor.remote.Process_monitor_Remotes import process_monitoring


def test():
    process_monitoring(interval=5)
    middle=[]
    for i in range(100):
        middle.append([i*"12345678"])
        time.sleep(0.2)

test()