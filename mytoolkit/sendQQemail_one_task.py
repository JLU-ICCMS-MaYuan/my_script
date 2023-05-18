#!/usr/bin/env python3
import os
import sys
import time
import smtplib
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

def check_slurm(job_ids=None):
    # 检查正在运行的任务数量
    tasks = []

    for job_id in job_ids:
        task = os.popen(f'squeue -u may -h --format "%A    %t    %R    %M    %Z     %s" | grep "{job_id}"').read().strip("\n")
        tasks.append(task)
        
    # 准备发送的内容
    content = "one_tasks" + "<br/>"
    for string in tasks:
        string = string + "<br/>"
        content+=string

    if len(tasks) > 1:
        return False, content
    else:
        return True, content
    
def check_pbs(job_ids=None):

    # 获得正在运行的任务的任务号 qstat -r
    tasks = []
    for job_id in job_ids:
        tdata = os.popen('qstat -f ' + job_id + ' |grep Output_Path -A 5').read().split('\n')
        state = os.popen('qstat' + '| grep ' + job_id + " | awk '{print $5}'").read().strip('\n')
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
        tasks.append("{}    {}    {}".format(job_id , state, '/' + '/'.join(predata)))

    # 准备发送的内容
    content = "one_tasks" + "<br/>"
    for string in tasks:
        string = string + "<br/>"
        content+=string
    content = content + "Waiting" + "<br/>"

    if len(tasks) > 1:
        return False, content
    else:
        return True, content


if __name__ == "__main__":


    meachine_name = "楼下集群"
    submit_system = "slurm"

    print("Note: --------------------")
    print("    使用方法: sendQQemail_one_tasks.py <job_ids>")
    job_ids = sys.argv[1:]
    previous_content  = ""
    while True:
        if submit_system == "slurm":
            send_flag, content = check_slurm(job_ids)
            if send_flag: # 发送邮件，就把content， previous_content 一起发送
                sendQQemail(previous_content, content)
                sys.exit(0)
            else: #  不发送邮件，就把content保存下来
                previous_content = content
        elif submit_system == "pbs":
            send_flag, content = check_pbs(job_ids)
            if send_flag: # 发送邮件，就把content， previous_content 一起发送
                sendQQemail(previous_content, content)
                sys.exit(0)
            else: #  不发送邮件，就把content保存下来
                previous_content = content
        else:
            print("当前机器用的提交任务的脚本既不是slurm也不是pbs, 退出程序")
            sys.exit(0)
        time.sleep(300)
