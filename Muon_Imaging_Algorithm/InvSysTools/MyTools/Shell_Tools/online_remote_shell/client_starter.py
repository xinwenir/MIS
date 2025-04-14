# author:高金磊
# datetime:2022/6/22 10:08
import os
import sys


BASE_DIR = os.path.dirname(__file__)
for i in range(2):###距离项目根目录的距离____终端使用
    BASE_DIR=os.path.dirname(BASE_DIR)
    sys.path.append(BASE_DIR)
    print(BASE_DIR)

#接后续代码

import time

import cowsay
from Shell_Tools.online_remote_shell.sockets import SocketClient,recver_obj, send_message
from Shell_Tools.online_remote_shell.Result import Result
from myPrint import *
from Data import *
cwd=None
animation=0.02
animation_multiple=1
cowsay_fun=myPrint_Hint
def remote_command(command,sock):
    global cwd,animation,animation_multiple,cowsay_fun
    if command == "exit":
        send_message(sock,command)
        time.sleep(1)
        sock.close()
    else:
        send_message(sock,command)

def client_command(command):
    global cwd,animation,animation_multiple,cowsay_fun
    if command=="":
        return False
    elif command == "setting show":
        animation_multiple=input("请输入客户端动画效果（0--10）：")
        try:
            if int(animation_multiple)>=0 and int(animation_multiple)<=10:
                myPrint_Success("设置已经生效")
                animation_multiple = int(animation_multiple)
        except Exception:
            myPrint_Err("设置失败，请输入0到10的整数")
            animation_multiple=1
        return False
    elif command=="setting cowsay":
        list_fun=[[name,cowsay.char_funcs[name]] for name in cowsay.char_funcs.keys()]
        list_fun.append(["紫色打印",myPrint_Hint])
        
        list_fun_copy=copy(list_fun)
        while 1:
            char,char_func=list_fun[0]
            char_func(f'Hi! I am {char}')
            middle=input("使用这个输入y，查看下一个输入n,查看上一个输入u,退出用E,显示所有用A,或者输入你想使用的编号(0,%s):  "%(len(list_fun)-1))
            if middle=='y':
                cowsay_fun=char_func
                cowsay_fun("设置成功")
                break
            elif middle=='n':
                list_fun.append(list_fun.pop(0))
            elif middle=='u':
                list_fun.insert(0,list_fun.pop(-1))
            elif middle=='E':
                break
            elif middle=='A' or middle=='a':
                for j in range(len(list_fun_copy)):
                    char,char_func= list_fun_copy[j]
                    char_func(f'Hi! I am {char},id is {j}')
            else:
                try:
                    middle = int(middle)
                    if(middle<len(list_fun_copy) and middle>=0):
                        cowsay_fun=list_fun_copy[middle][1]
                        cowsay_fun("设置成功")
                        break
                    else:
                        myPrint_Err("无效的命令")
                except Exception as e:
                    myPrint_Err("设置失败")
        return False
    return True
def _run( sock):
    global cwd,animation,animation_multiple,cowsay_fun
    while 1:
        if cwd is None:
            command=Data("win").pwd_cd()
        else:
            command = input(f'({Host}){cwd}\>')
        try:
            if not client_command(command):
                pass
            else:
                remote_command(command,sock)
                res:Result = recver_obj(sock)
                res_cpu(res)
        except Exception as e:
            if str(e).endswith("远程主机强迫关闭了一个现有的连接。"):
                #重连机制
                SocketClient().connect(_run)
                return


def res_cpu(res:Result):
    global cwd,animation,animation_multiple,cowsay_fun
    cwd=res.cwd
    if res.err!="" or cowsay_fun==myPrint_Hint or res.command=='help' or res.command=='who':
        for i in res.output.split('\n'):
            if animation_multiple>0:
                time.sleep(animation_multiple*animation)
                myPrint_Hint(i)
        for i in res.err.split('\n'):
            if animation_multiple>0:
                time.sleep(animation_multiple*animation)
                myPrint_Err(i)
    else:
        txt=res.output
        try:
            cowsay_fun(txt)
        except Exception:
            pass
if __name__ == '__main__':

    SocketClient().connect(_run)
