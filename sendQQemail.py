#/usr/bin/env python3
import sys
import smtplib
from pathlib import Path
from email.mime.text import MIMEText
from email.utils import formataddr

def sendQQemail(content):
    ### 邮箱内容配置 ###
    # 邮箱文本


    msg = MIMEText(content, "html", "utf-8")
    # 邮件上显示发件人
    msg['From'] = formataddr(['吴彦祖的小管家', 'myth620137018@163.com'])
    # 邮件显示的主题
    msg['Subject'] = "任务结束请查看"


    ### 发送邮件
    server = smtplib.SMTP_SSL("smtp.163.com")
    server.login("myth620137018@163.com", "YFKMHHHUFYFACJCE")
    server.sendmail("myth620137018@163.com", "1157421359@qq.com", msg.as_string())
    server.quit()

if __name__ == "__main__":

    meachine_name = sys.argv[1]
    task_name     = sys.argv[2]
    task_path     = Path.cwd()

    content = f'''机器名称：{meachine_name}
    任务名称：{task_name}
    任务路径：{task_path}
    '''
    sendQQemail(content)