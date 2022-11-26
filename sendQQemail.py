#!/usr/bin/env python3
import sys
import os
import time
import re
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

def check_taskID(task_id):
    while True:
        res1 = os.popen(" qstat | awk '{print $1}' | cut -d . -f1").read()
        res2 = os.popen("  ps -ef | awk '{print $2}' ").read()
        all_taskid = re.findall(r"\d+", res1) + re.findall(r"\d+", res2)
        if task_id not in all_taskid:
            return True
        time.sleep(5)


if __name__ == "__main__":

    print("使用帮助: sendQQemail.py 机器名称 任务名称 任务ID号 &")
    print("例如: sendQQemail.py 超硬机器 计算200GPaNp2B2 200332 &")

    meachine_name = sys.argv[1]
    task_name     = sys.argv[2]
    task_path     = Path.cwd()
    task_id       = sys.argv[3]

    content = f'''机器名称：{meachine_name}
    任务名称：{task_name}
    任务路径：{task_path}
    任务ID号: {task_id}
    '''
    if check_taskID(task_id):
        sendQQemail(content)