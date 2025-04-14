# author:高金磊
# datetime:2021/8/4 11:08
import warnings
from copy import copy

from colorama import Style, Fore


def _build(data: tuple):
    """
    格式化数据，使之尽可能接近实际打印()---不推荐使用
    :param data:
    :return:
    """
    warnings.warn("不再使用，效果没有内置函数好", DeprecationWarning)
    res = ''
    if len(data) == 0:
        return res
    for datum in data:
        ty = type(datum)
        if ty != int and ty != str:
            res += " " + str(datum) + " "
        else:
            res += str(datum)
    return res


def myPrint_Err(*s):
    """
    打印错误提示信息,显示为红色,功能等同print
    :param s:  要打印的数据,字符串,数组,多个字段等所有 print可以打印的
    :return: None
    """
    print(Fore.LIGHTRED_EX, *s, Style.RESET_ALL)


def myPrint_Wran(*s):
    """
    打印警告提示信息,显示为黄色,功能等同print
    :param s:  要打印的数据,字符串,数组,多个字段等所有 print可以打印的
    :return: None
    """
    print(Fore.YELLOW, *s, Style.RESET_ALL)


def myPrint_Hint(*s):
    """
    打印较弱的提示信息,显示为紫色,功能等同print
    :param s:  要打印的数据,字符串,数组,多个字段等所有 print可以打印的
    :return: None
    """
    print(Fore.LIGHTMAGENTA_EX, *s, Style.RESET_ALL)


def myPrint_Success(*s):
    """
    打印成功提示信息,显示为绿色,功能等同print
    :param s:  要打印的数据,字符串,数组,多个字段等所有 print可以打印的
    :return: None
    """
    print(Fore.GREEN, *s, Style.RESET_ALL)


class Loger_Print:
    """
    集成日志和控制台输出,如有需要重写log的所有方法
    """

    def __init__(self, loger, printer):
        """
        构造Loger_Print对象
        :param loger: 任何支持 write,flush,close方法的对象,如loger,file,网络连接等IO对象,也可以是None
        :param printer: print对象
        """
        self.loger = loger
        self.printer = printer

    def write(self, txt, loger=None, printer=None):
        """
        将txt输出到控制台和loger中。为保持一致loger也自动换行。
        因为loger对象可能有缓存，所以会导致loger和printer在某时刻的内容不同。
        :param txt: 要输出的数据
        :param loger: 参考__init__，默认是None，将使用self中的默认值
        :param printer: 参考__init__，默认是None，将使用self中的默认值
        :return: bool 时候成功
        """
        if loger is None:
            loger = self.loger
        if printer is None:
            printer = self.printer
        try:  # 日志记录不应该影响到主程序运算
            printer(txt)
            if loger is not None:
                loger.write(str(txt))
                loger.write('\n')
                self.flush()
            return True
        except Exception as e:
            myPrint_Err(e)
            return False

    def err(self, txt):

        return self.write(txt, printer=myPrint_Err)

    def success(self, txt):

        return self.write(txt, printer=myPrint_Success)

    def waring(self, txt):
        return self.write(txt, printer=myPrint_Wran)

    def info(self, txt):
        return self.write(txt, printer=print)

    def important_info(self, txt):
        return self.write(txt, printer=myPrint_Hint)

    def flush(self):
        """
        清空loger的缓存区
        :return:
        """
        self.loger.flush()

    def close(self):
        """
        关闭loger
        :return:
        """
        self.loger.close()


if __name__ == '__main__':
    myPrint_Err("1", "2", [1, 2])
    print(1)
    myPrint_Wran("333")
    myPrint_Hint(444)
    print("尽可能接近实际打印")
    myPrint_Hint([1, 2], 1, [2, 3])
    print([1, 2], 1, [2, 3])
