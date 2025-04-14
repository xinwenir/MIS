# author:高金磊
# datetime:2022/6/22 9:48
import threading
import os
import sys
import traceback
BASE_DIR = os.path.dirname(__file__)
for i in range(3):###距离项目根目录的距离____终端使用
    BASE_DIR=os.path.dirname(BASE_DIR)
    sys.path.append(BASE_DIR)
import socket
import time

import tqdm
import pickle
from Data import *

from myPrint import *
class SocketServer(object):

    def __init__(self, port=PORT):
        self.port = port

    def startup(self,fun,is_multithreading=True):
        sock_server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock_server.bind((Host, self.port))
        sock_server.listen(0)
    
        myPrint_Hint("===========等待连接================")
        while True:
            try:
                sock, address = sock_server.accept()
                myPrint_Hint("===========接收客户端密匙================")
                data = recver_message(sock)
                
                if data!=cilent_key:
                    myPrint_Wran("客户端keyid不正确", data, "请检查版本是否一致")
                    send_message(sock,sever_key)
                    send_message(sock,"拒绝连接")
                    sock.close()
                else:
                    myPrint_Hint("===========客户端密匙验证成功================")
                    myPrint_Hint("===========发送服务端密匙================")
                    send_message(sock,sever_key)
                    if is_multithreading:
                        myPrint_Success("===========连接成功================")
                        thread = threading.Thread(target=fun, args=(sock,))
                        thread.setDaemon(True)
                        thread.start()
                    else:
                        fun(sock)

            except Exception as e:
                traceback.print_exc()
                myPrint_Err(e)



    def state(self):
        pass


def send_message(socket,data):
    socket.sendall(data.encode())
    
    socket.sendall(Sever_Message_end.encode())
    
def send_blog_message(socket, data_list):
    """
    发送批量消息
    Args:
        socket (_type_): _description_
        data_list (_type_): _description_
    """
        
    socket.sendall(Sever_Message_start.encode())
    recver_message_block(socket, Sever_Message_confirm)
    for i in data_list:
        socket.sendall(i.encode())
        recver_message_block(socket, Sever_Message_confirm)
        print(i)
    socket.sendall(Sever_Message_end.encode())
    recver_message_block(socket, Sever_Message_confirm)

def send_obj(socket:socket,obj):
    socket.sendall(pickle.dumps(obj=obj))
    time.sleep(0.5)
    


class SocketClient(object):

    def __init__(self, host=Host, port=PORT):
        self.host = host
        self.port = port

    def connect(self,fun,time_out=20):
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        # sock.settimeout(4)
        for i in range(10):
            try:
                myPrint_Hint("=============连接中========================")
                sock.connect((self.host, self.port))
                time.sleep(0.5)
                myPrint_Success("=============发送客户端密匙========================")
                send_message(sock,cilent_key)
                time.sleep(0.5)
                myPrint_Success("=============接收服务端密匙========================")
                data=recver_message(sock)
                break
            except Exception as e:
                # if str(e).endswith("在一个已经连接的套接字上做了一个连接请求。"):
                #     break
                myPrint_Wran("连接失败,正在进行尝试........")
                myPrint_Hint(e)
                if i==time_out-1:
                    myPrint_Err("=============连接超时========================")
                    return -1
                time.sleep(2)
       
        myPrint_Success("=============验证服务端密匙========================")
        if data!=sever_key:
            myPrint_Wran("服务器keyid不正确",data,"请检查版本是否一致")
            send_message(sock,"拒绝连接")
            myPrint_Success("=============连接断开========================")
            sock.close()
        else:
            myPrint_Success("=============服务器允许接入========================")
            myPrint_Success("=============连接成功========================")
            myPrint_Hint("服务器keyid:", data)
            fun(sock)
        return sock






def recver_message_block(socket,expect_message):
    while recver_message(socket) !=expect_message:
        pass
    return True
def recver_message(socket):
    data=""
    while 1:
        middle=socket.recv(1024).decode('utf-8','ignore')
        if middle.endswith(Sever_Message_end):
            middle=middle.replace(Sever_Message_end,"")
            data+=middle
            break
        else:
            data+=middle
    return data
def recv_blog(socket):
    blog_data = ""
    while 1:
        # 接收批量结果
        data =recver_message(socket)

        if data == Sever_Message_end:
            send_message(socket,Sever_Message_confirm)
            return blog_data
        else:
            blog_data +="\n"+data
            send_message(socket,Sever_Message_confirm)
def recver_obj(socket):
    # time.sleep(1)
    data=b""
    # while 1:
    #     middle=socket.recv(1024)
    #     if middle.endswith(Sever_Message_end):
    #         middle=middle.replace(Sever_Message_end,"")
    #         data+=middle
    #         break
    #     else:
    #         data+=middle
    while 1:
        data+=socket.recv(1024)
        try:
            obj=pickle.loads(data)
            return obj
        except Exception as e:
            if(str(e)=="pickle data was truncated"):
                continue
            else:
                import Result
                import Data
                return Result.Result.getResult(processor="net",output="网络波动，数据异常",returncode=Data.NETWORK_ERROR)