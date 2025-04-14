class Data:
    def __init__(self,OS_TYPE) -> None:
        self.OS_TYPE=OS_TYPE
    def pwd_cd(self):
        return {"win":"cd","linux":"pwd"}[self.OS_TYPE]

Host="127.0.0.1"
# Host="56s815617a.goho.co"
PORT = 12522
EOF = '0x00'
QIAN_MING_IMAGE="E:\pycharm\MyTools\Image_Tools\签名.png"

sever_key="0x26916166291"#双方程序需要
cilent_key="0x3637126198"

Sever_Message_start="0x217732080"
Sever_Message_end="0x0q89888079"
Sever_Message_confirm="0x327092"

LOG_OUT="0x10003"
SUCESS="0x10200"
NOT_FOUND="0x10404"
NETWORK_ERROR="0x10400"
HELP=0x20001
NOT_HELP=0x20002

