
from Shell_Tools.online_remote_shell.Result import Result
from myPrint import *
def aop_preprocessing(command,fun):
    myPrint_Hint(f"指令：{command},下一个处理器：{fun.__name__}")

cache=None

def aop_postprocessing(res:Result):
    global cache
    if cache!=res:
        myPrint_Success(res)
        cache=res
    else:
        cache=res
def next_Cpu(command,funs):
    
    if len(funs)==0:
        res=Result.getResult(processor="next_cpu")
        aop_postprocessing(res)
        return res
    else:
        fun=funs.pop(0)
        aop_preprocessing(command,fun)
        res= fun(command,funs)
        aop_postprocessing(res)
        return res