# author:高金磊
# datetime:2022/6/22 10:08
import os
import subprocess
import sys


BASE_DIR = os.path.dirname(__file__)
for i in range(3):###距离项目根目录的距离____终端使用
    BASE_DIR=os.path.dirname(BASE_DIR)
    sys.path.append(BASE_DIR)
from copy import copy
from Result import Result
from Image_Tools.Image_char import get_txt_image
from Os_Tools.Os_Info import get_os_info
from Os_Tools.process_tool import Chain_Process_tool
from myPrint import myPrint_Success,myPrint_Err,myPrint_Hint
import Data
from Shell_Tools.online_remote_shell.Tools import next_Cpu
import pickle
cwd=os.getcwd()
debug=True
data_str=Data.Data("win")
def general_manage(command,funs)->Result:
    global cwd
    if cwd is None:
        res = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    else: 
        res = subprocess.Popen(command,cwd=cwd, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    code=res.returncode
    if command == "help" or command.split(" ")[-1]=="-h":
        code=Data.HELP
    d = res.stdout.readlines()
    e = res.stderr.readlines()
    if len(funs) !=0:
        myPrint_Err("责任链顺序有误,请把本方法置于责任链最后")
    res=""
    if d==None:
        d=""
    if e==None:
        e=""
    for i in d:
        res+="%s"%(i.decode("gbk"))
    err_res=""
    for i in e:
        err_res+="%s"%(i.decode("gbk"))
    return Result.getResult(processor="remote_cmd",returncode=code,output=res,err=err_res) 

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
    if command == "help":
        res = get_txt_image(Data.QIAN_MING_IMAGE)+"\n"
        for re in _all_fun.keys():
            res += str(re)
            res += '\n'
        return Result.getResult(processor="",output=res,returncode=Data.HELP) 
    try:
        if command not in _all_fun.keys():
            return next_Cpu(command,funs)
        else:
            return Result.getResult(processor="os_info",returncode=1,output=_all_fun[command]()) 
    except Exception as e:
        return next_Cpu(command,funs)

def exit(command,funs):
        if command == "exit" or command == "logout":  
            return Result.getResult(output="good bye",returncode=Data.LOG_OUT)
        else:
           return next_Cpu(command,funs)
def help(command : str,funs):
    """
    !todo 实现所有方法的help

    Args:
        command (_type_): _description_
        funs (_type_): _description_

    Returns:
        _type_: _description_
    """
    if command == "help" or command.split(" ")[-1]=="help": 
        output="" 
        for fun in funs:
            res:Result=fun(command,[])
            if res.returncode==Data.HELP:
                output+=res.output+res.err+"\n"
        return Result.getResult(output=output)
                    
    else:
        return next_Cpu(command,funs)
def fliter(command,funs):
    if command == "help" or command.split(" ")[-1]=="help":
        return Result.getResult(processor="fliter",output="这是fliter，暂时没有帮助文档\n",returncode=Data.HELP)
    if command == "":  
        return Result.getResult(processor="fliter")
    elif command == "python":
        return Result.getResult(processor="fliter",output="不能使用交互式python终端",returncode=2)
    else:
        return next_Cpu(command,funs)
def _run(socket):
    funs=[help,fliter,cd_pwd,exit,other_command,get_os_info,Chain_Process_tool,general_manage]
    global cwd
    
    while 1:
        data = recver_message(socket)
        myPrint_Success("收到:", data)

        res:Result=next_Cpu(data,copy(funs))
        res.cwd=cwd
        res.command=data
        # myPrint_Success(res)
        send_obj(socket, res)
def cd_pwd(command:str,funs:list):
    global cwd
    if command==data_str.pwd_cd():
        return Result.getResult(processor="cd_pwd",output=cwd)
    elif command.split(" ")[0] == data_str.pwd_cd():
        path=command.split(" ")[1]
        try:
            middle=os.path.abspath(os.path.join(cwd,path))
            if os.path.exists(middle):
                cwd=middle
            else:
                raise Exception("目录不存在") 
        except Exception as e:
            return Result.getResult(processor="cd_pwd",err=path+"  不存在",returncode=2)
        return Result.getResult(processor="cd_pwd",err=cwd,returncode=1)
    return next_Cpu(command,funs)
if __name__ == '__main__':
    from  Shell_Tools.online_remote_shell.sockets  import SocketServer, recver_message, send_obj

    SocketServer().startup(_run)