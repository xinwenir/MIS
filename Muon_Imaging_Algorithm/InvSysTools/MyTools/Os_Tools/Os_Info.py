# author:高金磊
# datetime:2022/7/31 20:29
# coding=utf-8

# 获取系统基本信息

import sys

import psutil

import time

import os
from Shell_Tools.online_remote_shell.Result import Result
from Shell_Tools.online_remote_shell.Tools import next_Cpu
import Data
#获取当前时间

time_str = time.strftime( "%Y-%m-%d", time.localtime( ) )

def get_os_info(command:str,funs):
    _all_fun = {
                   "os_mem":_mem, "os_cpu":_cpu, "os_user_info":_os_user, "os_disk":_disk
    }
    
    if command=="help":
        res=""
        for re in _all_fun.keys():
            res+=str(re)
            res+='\n'
        return Result.getResult(processor="get_os_info",output=res,returncode=Data.HELP)
    try:
        if command not in _all_fun.keys():
            return next_Cpu(command,funs)
        else:
            return Result.getResult(processor="get_os_info",output=_all_fun[command](),returncode=Data.SUCESS)
    except Exception as e:
        return next_Cpu(command,funs)

#获取系统内存使用情况
def _mem():

    memory_convent = 1024 * 1024

    mem = psutil.virtual_memory()

    print_str = " 内存状态如下:\n"

    print_str = print_str + " 系统的内存容量为: "+str( mem.total/( memory_convent ) ) + " MB\n"

    print_str = print_str + " 系统的内存以使用容量为: "+str( mem.used/( memory_convent ) ) + " MB\n"

    print_str = print_str + " 系统可用的内存容量为: "+str( mem.total/( memory_convent ) - mem.used/( 1024*1024 )) + "MB\n"

    # print_str = print_str + " 内存的buffer容量为: "+str( mem.buffers/( memory_convent ) ) + " MB\n"

    # print_str = print_str + " 内存的cache容量为:" +str( mem.cached/( memory_convent ) ) + " MB\n"
    return print_str
#获取cpu的相关信息

def _cpu():


    print_str = " CPU状态如下:\n"

    cpu_status = psutil.cpu_times()

    print_str = print_str + " user = " + str( cpu_status.user ) + "\n"

    # print_str = print_str + " nice = " + str( cpu_status.nice ) + "\n"

    print_str = print_str + " system = " + str( cpu_status.system ) + "\n"

    print_str = print_str + " idle = " + str ( cpu_status.idle ) + "\n"

    # print_str = print_str + " iowait = " + str ( cpu_status.iowait ) + "\n"

    # print_str = print_str + " irq = " + str( cpu_status.irq ) + "\n"

    # print_str = print_str + " softirq = " + str ( cpu_status.softirq ) + "\n"

    # print_str = print_str + " steal = " + str ( cpu_status.steal ) + "\n"

    # print_str = print_str + " guest = " + str ( cpu_status.guest ) + "\n"
    return print_str
#查看硬盘基本信息

def _disk():

    print_str = " 硬盘信息如下:\n"

    disk_status = psutil.disk_partitions()

    for item in disk_status :

        print_str = print_str + " "+ str( item ) + "\n"
    return print_str
#查看当前登录的用户信息

def _os_user():

    print_str = " 登录用户信息如下:\n "

    user_status = psutil.users()

    for item in user_status :

        print_str = print_str + " "+ str( item ) + "\n"

        print_str += "---------------------------------------------------------------\n"

        print ( print_str )
    return print_str

