# author:高金磊
# datetime:2022/6/22 10:08
import os
import sys
BASE_DIR = os.path.dirname(__file__)
for i in range(3):###距离项目根目录的距离____终端使用
    BASE_DIR=os.path.dirname(BASE_DIR)
sys.path.append(BASE_DIR)
from copy import copy

from MyTools.Image_Tools.Image_char import get_txt_image
from MyTools.Os_Tools.Os_Info import get_os_info
from MyTools.Os_Tools.process_tool import Chain_Process_tool
from MyTools.myPrint import myPrint_Success,myPrint_Err

def general_manage(data,funs):
    if data =='-h':
        data='help'
    f = os.popen(data, 'r', )
    d = f.readlines()
    f.close()
    if len(funs) != 1:
        myPrint_Err("责任链顺序有误,请把本方法置于责任链最后")
    res=" "
    for i in d:
        res+="%s\n"%(i)
    return ["",res]
def other_command(command,funs):
    def hello():
       return "  _   _      _ _    __        __         _     _ \n" + \
        "| | | | ___| | | __\\ \\      / /__  _ __| | __| |\n" + \
        "| |_| |/ _ \\ | |/ _ \\ \\ /\\ / / _ \\| '__| |/ _` |\n" + \
        "|  _  |  __/ | | (_) \\ V  V / (_) | |  | | (_| |\n" + \
        "|_| |_|\\___|_|_|\\___/ \\_/\\_/ \\___/|_|  |_|\\__,_|"
    _all_fun = {
       "hello":hello,
        "who":get_txt_image
    }
    funs.pop(0)
    if command == "-h":
        res = get_txt_image()+"\n"
        for re in _all_fun.keys():
            res += str(re)
            res += '\n'
        return ["", res + funs[0](command, funs)[1]]
    try:
        if command not in _all_fun.keys():
            return funs[0](command,funs)
        else:
            return ["get_os_info",_all_fun[command]()]
    except Exception as e:
        return funs[0](command,funs)

def _run(socket):
    funs=[other_command,get_os_info,Chain_Process_tool,general_manage]
    while 1:
        data = recver_message(socket)
        myPrint_Success("收到:", data)
        if data == "exit":
            socket.close()
            return 1
        else:
            if data=="":
                continue

        response_message=funs[0](data,copy(funs))
        myPrint_Success(response_message)
        send_blog_message(socket, [response_message[1]])

if __name__ == '__main__':
    from  MyTools.Socket_Tools.sockets  import SocketServer, recver_message, send_blog_message

    SocketServer().startup(_run)