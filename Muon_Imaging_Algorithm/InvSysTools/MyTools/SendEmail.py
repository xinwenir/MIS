# author:高金磊
# datetime:2021/11/11 21:26
import smtplib
from email.mime.text import MIMEText
from email.utils import formataddr
password=None
my_nick=None
sender=None
receivers=None
to_nick=None
# 自定义发送邮件的函数
def mail(subject,mail_msg,my_pass=password,my_nick=my_nick,to_nick=to_nick,my_sender=sender,to_user=receivers):
    """
    配置发邮件所需的基础信息
    :param my_sender: 配置发件人邮箱地址
    :param my_pass: 配置发件人邮箱密码
    :param to_user: 配置收件人邮箱地址
    :param my_nick: 配置发件人昵称
    :param to_nick: 配置收件人昵称
    :param mail_msg: 配置邮件内容
    :return: 状态
    """
    try:
        for rec in to_user:
            # 必须将邮件内容做一次MIME的转换 -- 这是发送含链接的邮件
            msg=MIMEText(mail_msg,'html','utf-8')
            # 配置发件人名称和邮箱地址
            msg['From']=formataddr([my_nick,my_sender])
            # 配置收件件人名称和邮箱地址
            msg['To']=formataddr([to_nick,rec])
            # 配置邮件主题（标题）
            msg['Subject']=subject
            # 配置Python与邮件的SMTP服务器的连接通道（如果不是QQ邮箱，SMTP服务器是需要修改的）
            server=smtplib.SMTP_SSL("smtp.qq.com", 465)
            # 模拟登陆
            server.login(my_sender, my_pass)
            # 邮件内容的发送
            server.sendmail(my_sender,[rec,],msg.as_string())
            # 关闭连接通道
            server.quit()
    except Exception:
        return False
    return True