class Result:
    def __init__(self) -> None:
        self.processor=""
        self.returncode=0
        self.output=""
        self.err=""
        self.cwd=""
        self.command=""
        
    def __str__(self) -> str:
        return  f"cwd:{self.cwd} \n processor: {self.processor}\n "+ f"returncode:{self.returncode}\n "+ f"output: {self.output}\n "+ f"err:{self.err}"
    @staticmethod
    def getResult(processor="",returncode=1,output="",err="",cwd="",command=""):
        if processor is None:
            processor=""
        if returncode is None:
            returncode=1
        if output is None:
            output=""
        if err is None:
            err=""
        if cwd is None:
            cwd=""
        if command is None:
            command=""
        
        res=Result()
        res.processor=processor
        res.returncode=returncode
        res.err=err
        res.output=output
        res.cwd=cwd
        res.command=command
        return res
        