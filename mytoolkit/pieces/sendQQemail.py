#!/usr/bin/env python3
import sys
import os
import time
import re
import smtplib
from pathlib import Path
from email.mime.text import MIMEText
from email.utils import formataddr

def sendQQemail(previous_content, content):
    ### 邮箱内容配置 ###
    # 邮箱文本
    total_content = "上一次的任务检查结果" + "<br/>" + previous_content + "<br/><br/><br/>" + "当前的任务检查结果<br/>" + content 

    msg = MIMEText(total_content, "html", "utf-8")
    # 邮件上显示发件人
    msg['From'] = formataddr(['吴彦祖的小管家', 'myth620137018@163.com'])
    # 邮件显示的主题
    msg['Subject'] = "任务结束请查看"


    ### 发送邮件
    server = smtplib.SMTP_SSL("smtp.163.com")
    server.login("myth620137018@163.com", "YFKMHHHUFYFACJCE")
    server.sendmail("myth620137018@163.com", "1157421359@qq.com", msg.as_string())
    server.quit()

def check_slurm(node_num=4):
    # 检查正在运行的任务数量
    rtasks = os.popen('squeue -u may -h --format "%A    %t    %R    %M    %Z     %s" | grep " R "').read().split('\n')
    running = [ t for t in rtasks if 'may' in t]

    # 检查正在排队的任务数量
    wtasks = os.popen('squeue -u may -h --format "%A    %t    %R    %M    %Z     %s" | grep " PD "').read().split('\n')
    waiting = [ t for t in wtasks if 'may' in t]

    # 准备发送的内容
    content = "Running" + "<br/>"
    for string in running:
        string = string + "<br/>"
        content+=string
    content = content + "Waiting" + "<br/>"
    for string in waiting:
        string = string + "<br/>"
        content+=string

    if len(running) < node_num:
        return True, content
    else:
        return False, content
    
def check_pbs(node_num=10):

    # 获得正在运行的任务的任务号 qstat -r
    tasknumber = os.popen("qstat -r | grep liuhy | awk '{print $1}' | cut -d . -f 1").read().split()
    # 遍历任务号，找到正在运行的任务
    running = []
    for num in tasknumber:
        tdata = os.popen('qstat -f ' + num + ' |grep Output_Path -A 5').read().split('\n')
        state = os.popen('qstat' + '| grep ' + num + " | awk '{print $5}'").read().strip('\n')
        for i in range(0, len(tdata)):
            if 'Priority' in tdata[i]:
                j = i
        ndata = tdata[:j]
        oneline = ''
        # 制作路径
        for iterm in ndata:
            iterm = iterm.strip()
            oneline += iterm
        predata = oneline.split('/')[1:-1]
        running.append("{}    {}    {}".format(num , state, '/' + '/'.join(predata)))

    # 获得正在排队的任务的任务号 qstat -i
    waiting = []
    tasknumber = os.popen("qstat -i | grep liuhy | awk '{print $1}' | cut -d . -f 1").read().split()
    for num in tasknumber:
        tdata = os.popen('qstat -f ' + num + ' |grep Output_Path -A 5').read().split('\n')
        state = os.popen('qstat' + '| grep ' + num + " | awk '{print $5}'").read().strip('\n')
        for i in range(0, len(tdata)):
            if 'Priority' in tdata[i]:
                j = i
        ndata = tdata[:j]
        oneline = ''
        for iterm in ndata:
            iterm = iterm.strip()
            oneline += iterm
        predata = oneline.split('/')[1:-1]
        waiting.append("{}    {}    {}".format(num , state, '/' + '/'.join(predata)))

    # 准备发送的内容
    content = "Running" + "<br/>"
    for string in running:
        string = string + "<br/>"
        content+=string
    content = content + "Waiting" + "<br/>"
    for string in waiting:
        string = string + "<br/>"
        content+=string

    if len(running) < node_num:
        return True, content
    else:
        return False, content


if __name__ == "__main__":


    meachine_name = "楼下集群"
    submit_system = "slurm"

    previous_content  = ""
    while True:
        if submit_system == "slurm":
            send_flag, content = check_slurm(6)
            if send_flag: # 发送邮件，就把content， previous_content 一起发送
                sendQQemail(previous_content, content)
                sys.exit(0)
            else: #  不发送邮件，就把content保存下来
                previous_content = content
        elif submit_system == "pbs":
            send_flag, content = check_pbs(10)
            if send_flag: # 发送邮件，就把content， previous_content 一起发送
                sendQQemail(previous_content, content)
                sys.exit(0)
            else: #  不发送邮件，就把content保存下来
                previous_content = content
        else:
            print("当前机器用的提交任务的脚本既不是slurm也不是pbs, 退出程序")
            sys.exit(0)
        time.sleep(1800)
