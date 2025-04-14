# author:高金磊
# datetime:2022/6/22 10:08
import os
import sys
BASE_DIR = os.path.dirname(__file__)
for i in range(3):###距离项目根目录的距离____终端使用
    BASE_DIR=os.path.dirname(BASE_DIR)
sys.path.append(BASE_DIR)

#接后续代码

import time

import cowsay

from MyTools.Socket_Tools.sockets import SocketClient, EXIT, Sever_Message_start, send_message, Sever_Message_confirm, recver_message, \
    Sever_Message_end, recv_blog
from MyTools.myPrint import *


def _run( socket):
    animation=0.02
    animation_multiple=1
    cowsay_fun=cowsay.dragon
    while 1:
        command = input('gjl:\>：')
        if command == "exit":
            socket.send(EXIT.encode())
            socket.close()
            break
        elif command == "setting show":
            animation_multiple=input("请输入客户端动画效果（0--10）：")
            try:
                if int(animation_multiple)>=0 and int(animation_multiple)<=10:
                    myPrint_Success("设置已经生效")
                    animation_multiple = int(animation_multiple)
            except Exception:
                myPrint_Err("设置失败，请输入0到10的整数")
                animation_multiple=1
            continue
        elif command=="setting cowsay":
            list_fun=[[i,cowsay.chars[i]] for i in cowsay.chars.keys()]
            list_fun.append(["紫色打印",myPrint_Hint])
            while 1:
                char,char_func=list_fun[0]
                char_func(f'Hi! I am {char}')
                middle=input("使用这个输入y，查看下一个输入n,查看上一个输入u,退出用E:  ")
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

        try:
            socket.send(command.encode())
        except Exception as e:
            if str(e).endswith("远程主机强迫关闭了一个现有的连接。"):
                #重连机制
                SocketClient().connect(_run)
                return

        while 1:
            data = socket.recv(1024).decode('utf-8','ignore')
            if data == "":
                time.sleep(0.5)
                continue

            if data == Sever_Message_start:
                data = recv_blog(socket)
            if cowsay_fun==myPrint_Hint or command=='-h' or command=='who':
                for i in data.split('\n'):
                    if animation_multiple>0:
                        time.sleep(animation_multiple*animation)
                    myPrint_Hint(i)
            else:
                txt=data

                try:
                    cowsay_fun(txt)
                except Exception:
                    pass
            break
if __name__ == '__main__':

    SocketClient().connect(_run)
