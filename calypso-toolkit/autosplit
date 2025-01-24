#!/usr/bin/python3

__author__ = 'tqc'

# 运行逻辑
#   第一步：
#       a = Autocalypso() 初始化类参数
#   第二步：
#       a.autorun() 运行实例方法
#   第三步：
#       在autorun（）中，
#           先判断是否存在相应的提作业的脚本: self.checkfiles()
#           然后读取input.dat:              self.readinput()
#           最后开始运行vasp：              self.control_vasp()

# control_vasp()的基本逻辑：
#   1. 根据input.dat中 NumberOfParallel = 10 的开关，确定单次提交的任务数(即：self.NP），一般它要小于PopSize
#   2. 如果发现当前所有任务（正在运行的任务数+排队的任务数）< self.NP，那么就补交新的任务，使总任务数=self.NP
#   注意：
#       1. 可以同时提交多个calpso结构预测任务, 它可以通过取集合交集的方式获得某个calypso预测已提交的总任务数
#          如果某个calypso预测已提交的总任务数 < NumberOfParallel，那么就补交新的任务。
#       2. 
# 关于续算的问题：
#   1. 不用开启input.dat中的PickUp， 保持 PickUp = F
#   2. 但是要保证有step这个文件，你如果要在./calypso.x后迭代出第n代的结构，就要保证step中记录的是n-1
#   3. 如果要续算的话，就把current_stepfile中的POSCAR_*, CONTCAR_*, OUTCAR_*拷贝到当前目录，然后清理好results中的pso_ini_*, pso_opt_*, pso_sor_*
#      然后执行一下 write_vsave.py <Popsize>
#   这时你就是准备好了续算的所有前期工作，可以执行 nohup python autosplit.py slurm > 2>&1 &
 

import os
import time
import sys
import glob

class Autocalypso(object):
    def __init__(self,submit = 'sbatch vasp.slurm',stat = 'squeue',rstat = 'squeue | grep R',delete = 'scancel',\
                calypath = './calypso.x',machine = 'slurm'):
        self.CalyPsoPath = calypath
        self.submit = submit
        self.stat = stat
        self.rstat = rstat
        self.delete = delete
        self.machine = machine
        self.MaxTime = 5400
        self.PopSize = 360
        self.MaxStep = 300
        self.NumberOfParallel = 16
        self.PickUp = False
        self.NumberOfLocalOptim = 4
        self.f = open('split_calypso.log','w')

    def readinput(self):
        f = open('input.dat','r')
        while True:
            line = f.readline()
            if len(line) == 0:
                break
            if 'PopSize' in line:
                self.StrNum = int(line.split('=')[1])
            if 'MaxStep' in line:
                self.GenNum = int(line.split('=')[1])
            if 'MaxTime' in line:
                self.MaxTime = int(line.split('=')[1])
            if 'NumberOfLocalOptim' in line:
                self.NumberOfLocalOptim = int(line.split('=')[1])
            if 'NumberOfParallel' in line:
                self.NP = int(line.split('=')[1])
                if self.StrNum < (self.NP):
                    self.NP = self.StrNum
            if 'PickUp' in line:
                #print line.strip().split('=')[1]
                pickup = line.strip().split('=')[1]
                #if pickup == 'T':
                if 'T' in pickup:
                    self.PickUp = True
                else:
                    self.PickUp = False
                #print self.PickUp
                #sys.exit(0)

    def checkfiles(self):
        if self.machine == 'pbs':
            if not os.path.exists(r'./vasp.pbs'):
                print('There is no vasp.pbs file!!!')
                #print 'We will generate a file, and you need check the PATH of VASP and the number of node!\n'
                print('We will generate another file, maybe you need change it!!!')
                self.writevasppbs()
                sys.exit(0)
        elif self.machine == 'lsf':
            if not os.path.exists(r'./run.lsf'):
                print('No run.lsf file!!!')
                print('We will generate another file, maybe you need change it!!!')
                self.writerunlsf()
                sys.exit(0)
        elif self.machine == 'yh':
            if not os.path.exists(r'./vasp.sh'):
                print('No vasp.sh file!!!')
                print('We will generate another file, maybe you need change it!!!')
                self.writeyh()
                sys.exit(0)
        elif self.machine == 'slurm':
            if not os.path.exists(r'./vasp.slurm'):
                print('No vasp.slurm file!!!')
                print('We will generate another file, maybe you need change it!!!')
                self.writeslurm()
                sys.exit(0)
        if not os.path.exists(r'./input.dat'):
            print('No input.dat')
            sys.exit(0)
        elif not os.path.exists(r'./POTCAR'):
            print('No POTCAR')
            sys.exit(0)
        elif len(glob.glob(r'./INCAR_*')) == 0:
            print('No INCAR files')
            sys.exit(0)
        else:
            print('Check files completed!!!')

    def lpickup(self):
        if not self.PickUp:
            os.system('rm step')

    def submit_vasp(self,n):
        self.f.write("%s\n" % "caly_auto_split call CALYPSO generare structure" )
        self.f.flush()
        # 执行 ./calypso.x 
        os.system(self.CalyPsoPath) 
        if n != self.GenNum:
            self.f.write("%s\n" % "caly_auto_split submit structure relax job" )
            self.f.flush()
            self.control_vasp()

    def autorun(self):
        self.checkfiles()
        self.readinput()
        # self.lpickup()
        i = 0
        while i < self.GenNum:
            self.f.write("%s %s %s\n" % ("=================",str(i),"ITERATION ==================" ))
            self.f.flush() # 立即写入缓冲区, 在这里，缓冲区中的数据已被写入文件，可以安全地关闭文件对象了
            self.backupfiles()
            self.submit_vasp(i)
            i += 1

    def backupfiles(self):
        # 备份当前代计算好的POSCAR_*, CONTCAR_*, OUTCAR_*
        if not os.path.exists("current_stepfiles"):
            os.mkdir("current_stepfiles")
        for i in range(self.StrNum):
            if os.path.exists(r'POSCAR_%s' % (str(i+1))):
                os.system('cp POSCAR_%s  current_stepfiles/POSCAR_%s' % (str(i+1),str(i+1)))
                print("finish backup POSCAR_%s" % (str(i+1)))
            if os.path.exists(r'CONTCAR_%s' % str(i+1)):
                os.system('cp CONTCAR_%s  current_stepfiles/CONTCAR_%s' % (str(i+1),str(i+1)))
                print("finish backup CONTCAR_%s" % (str(i+1)))
            if os.path.exists(r'OUTCAR_%s' % str(i+1)):
                os.system('cp OUTCAR_%s  current_stepfiles/OUTCAR_%s' % (str(i+1),str(i+1)))
                print("finish backup OUTCAR_%s" % (str(i+1)))
        time.sleep(3); print("finish backup")

    def control_vasp(self):
        for i in range(self.StrNum):
            os.system('rm -rf %s' %str(i+1))
            os.system('mkdir %s' % str(i+1))
            if self.machine == 'pbs':
                os.system('cp vasp.pbs POTCAR INCAR_*  %s' % str(i+1))
            elif self.machine == 'lsf':
                os.system('cp run.lsf POTCAR INCAR_*  %s' % str(i+1))
            elif self.machine == 'yh':
                os.system('cp vasp.sh POTCAR INCAR_*  %s' % str(i+1))
            elif self.machine == 'slurm':
                os.system('cp vasp.slurm POTCAR INCAR_*  %s' % str(i+1))
            os.system('cp POSCAR_%s  %s/POSCAR' % (str(i+1),str(i+1)))
            os.system('cp POSCAR_%s  %s/POSCAR_%s' % (str(i+1),str(i+1),str(i+1)))
        self.f.write("%s\n" % "Set structure relax jobs finished")
        totaljobid = []
        runjobid = []
        splitjobid = []
        split_run_jobid = []
        jobtime = {}
        for i in range(self.NP): #  self.NP=50 单次并行提交的任务数，
            if self.machine == 'pbs':
                id = int(os.popen(' cd %s; %s;cd ..' % (str(i+1), self.submit)).read().split('.')[0])
            elif self.machine == 'lsf':
                id = int(os.popen(' cd %s; %s;cd ..' % (str(i+1),self.submit)).read().split(' ')[1].\
                                            strip('<').strip('>'))
            elif self.machine == 'yh':
                id = int(os.popen(' cd %s; %s;cd ..' % (str(i+1),self.submit)).read().split(' ')[3])
            elif self.machine == 'slurm':
                id = int(os.popen(' cd %s; %s;cd ..' % (str(i+1),self.submit)).read().split(' ')[3])
            splitjobid.append(id) # 第1次给splitjobid赋值，self.NP就是splitjobid的长度
            jobtime[id] = 0 
        self.f.write("%s\n" % "Submit structure relax jobs finished")
        #print(splitjobid
        num = self.NP
        finnum = 0   
        tenode = 0
        nover = 0
        while finnum < self.StrNum: # self.StrNum 一代提交的任务数就是PopSize, finnum已经完成的任务数
            totaljobid = self.run_jobid(self.stat)
            runjobid = self.run_jobid(self.rstat)
            splitjobid = list(set(splitjobid).intersection(set(totaljobid))) # 第2次splitjobid赋值，找到所有正在运行和正在排队的任务数
            self.f.write("%s\t" % "Split job ID")
            for job_id in splitjobid:
                self.f.write("%d\t" % job_id)
            self.f.write("\n")
            split_run_jobid = list(set(splitjobid).intersection(set(runjobid)))
            self.f.write("%s\t" % "Split RUN job ID")
            for job_id in split_run_jobid:
                self.f.write("%d\t" % job_id)
            self.f.write("\n" )
            self.f.write("%s\n" % "----------------------------------" )
            for ii in split_run_jobid:
                jobtime[ii] += 60
                if jobtime[ii] > self.MaxTime:
                    print('scancel %d' % ii)
                    os.system('%s %s' % (self.delete, str(ii)))  
            #print(jobtime  
            enode = len(splitjobid)  # splitjobid=50 是一次性提交计算的任务数， 等于self.NP
                                     # enode 表示当前所有没算完的任务（包括真正计算的和在排队的任务）
            if enode == 0:
                nover += 1 
            if nover == 5:
                break
            if enode < tenode:
                finnum += (tenode-enode)
            tenode = enode #  如果是第一次进入该while循环，tenode=self.NP, tenode就是 一次提交的self.NP个并行的任务数量
                           #  后续每次循环， tenode都是上次循环统计出来的所有没算完的任务（包括真正计算的和在排队的任务）
            fnode = self.NP-enode # 总的并行任务数 - 当前没算完的任务 = 空置的任务数
            self.f.write("%d  %d   %d   %d\n" %(fnode, enode, self.NP, finnum))
            if fnode > 0:  # 空置的任务数 > 0
                if num < self.StrNum: # self.StrNum 一代提交的任务数
                    for i in range(fnode): # 根据你空置的任务数，补交新的任务，新补交的任务数与当前所有任务数的总和等于self.NP
                        if (num+i+1) <= self.StrNum:
                            if self.machine == 'pbs':
                                id = int(os.popen(' cd %s; %s ;cd ..' % (str(num+i+1),self.submit)).\
                                                read().split('.')[0])
                            elif self.machine == 'lsf':
                                id = int(os.popen(' cd %s; %s ;cd ..' % (str(num+i+1),self.submit)).\
                                        read().split(' ')[1].strip('<').strip('>'))
                            elif self.machine == 'yh':
                                id = int(os.popen(' cd %s; %s ;cd ..' % (str(num+i+1),self.submit)).\
                                                read().split(' ')[3])
                            elif self.machine == 'slurm':
                                id = int(os.popen(' cd %s; %s ;cd ..' % (str(num+i+1),self.submit)).\
                                                read().split(' ')[3])
                            splitjobid.append(id)
                            jobtime[id] = 0 
                    num += fnode
            self.f.write("%s\t" % "Run number ")
            self.f.write("%d\n" % num)
            time.sleep(60)
            #print('submitted',num, 'finished',finnum

        # 这里表示完成了当前代所有的结构优化，开始读取结构
        self.f.write("%s\n" % "Structure relax jobs finished")
        for i in range(self.StrNum):
            if os.path.exists(r'./%s/CONTCAR' % str(i+1)):
                os.system('cp %s/CONTCAR  CONTCAR_%s' % (str(i+1),str(i+1)))
            else:
                print(i + 1,'th VASP JOB WRONG NO CONTCAR')
            if os.path.exists(r'./%s/OUTCAR' % str(i+1)):
                os.system('cp %s/OUTCAR  OUTCAR_%s' % (str(i+1),str(i+1)))
            else:
                print(i + 1,'th VASP JOB WRONG NO OUTCAR')
                os.system('echo 0 > OUTCAR_%s' % str(i+1))
        self.f.write("%s\n" % "cp OUTCAR and CONTCAR to .. finished")

    def run_jobid(self,cmd):
        runjobid = []
        aa = os.popen('%s' % cmd).readlines()
        for line in aa:
            try:
                if self.machine == 'pbs':
                    runjobid.append(int(line.split('.')[0]))
                elif self.machine == 'lsf':
                    runjobid.append(int(line.split(' ')[0]))
                elif self.machine == 'yh':
                    runjobid.append(int(line.strip().split(' ')[0]))
                elif self.machine == 'slurm':
                    runjobid.append(int(line.strip().split(' ')[0]))
            except:
                continue
        return runjobid

    def split_jobid(self):
        splitjobid = []
        for i in range(self.NP):
            splitjobid.append(int(os.popen(' cd %s;sbatch vasp.slurm;cd ..' % (str(i+1))).read().split('.')[0]))
        return splitjobid
                
    def writevasppbs(self,nodes = '1',ppn = '24',vasppath = '/work/software/vasp.6.1.0/vasp_std'):

        f = open('vasp.pbs','w')
        f.write('#!/bin/sh -f\n')
        f.write('#PBS -N split\n')
        f.write('#PBS -q hxl\n')
        f.write('#PBS -l nodes=%s:ppn=%s\n' % (nodes,ppn))
        f.write('#PBS -l walltime=4800:00:00\n')
        f.write('#PBS -V\n')
        f.write('#PBS -S /bin/bash\n')
        f.write('source /work/env/intel2020\n')
        f.write('cd $PBS_O_WORKDIR\n')
        #f.write('num_in=`ls -l |grep 'INCAR_' |wc -l`\n')
        f.write('for(( i=1; i<=%s; i++ ));\n' %  str(self.NumberOfLocalOptim))
        f.write('do\n')
        f.write('\tcp INCAR_$i INCAR\n')
        f.write('\tcp CONTCAR POSCAR\n')
        f.write('\tmpirun -np %s  %s > vasp.log_$i\n' % (ppn,vasppath))
        f.write('done\n')
        f.write('rm -rf CHG* WAVECAR\n')
        f.close()

    def writerunlsf(self,nodes = '1',ppn = '24',vasppath = '/apps/vasp/vasp.5.3.2'):
        f = open('run.lsf','w')
        f.write('#!/bin/sh\n')
        f.write('#BSUB -n 24\n')
        f.write("#BSUB -R 'span[ptile=24]'\n")
        f.write('for(( i=1; i<=%s; i++ ));\n' %  str(self.NumberOfLocalOptim))
        f.write('do\n')
        f.write('\tcp INCAR_$i INCAR\n')
        f.write('\tcp CONTCAR POSCAR\n')
        f.write('\tmpirun -np %s  %s > out_$i\n' % (ppn,vasppath))
        f.write('done\n')
        f.close()

    def writeyh(self):  
        f = open('vasp.sh','w')
        f.write('#!/bin/bash\n')
        f.write('for(( i=1; i<=%s; i++ ));\n' %  str(self.NumberOfLocalOptim))
        f.write('do\n')
        f.write('\tcp INCAR_$i INCAR\n')
        f.write('\tcp CONTCAR POSCAR\n')
        f.write('\tyhrun -N 1 -n 24 -p TH_NET /vol-th/home/maym/software/vasp.5.4.1/bin/vasp_std > vasp.log 2>&1\n')
        f.write('done\n')
        f.close()
            
    def writeslurm(self):
        f = open('vasp.slurm','w')
        f.write('#!/bin/sh\n')
#        f.write('source $MODULESHOME/init/sh\n')
#        f.write('module purge\n')
#        f.write('module add app_env/slurm\n')
#        f.write('module add app_env/vasp-5.4.1\n')
#        f.write('for(( i=1; i<=%s; i++ ));\n' %  str(self.NumberOfLocalOptim))
#        f.write('do\n')
#        f.write('\tcp INCAR_$i INCAR\n')
#        f.write('\tcp CONTCAR POSCAR\n')
#        f.write('\tsrun vasp_parallel\n')
#        f.write('done\n')
#        f.close()
        f.write('#SBATCH  --job-name=ll1118\n')
        f.write('#SBATCH  --output=log.out.%j\n')
        f.write('#SBATCH  --error=log.err.%j\n')
        f.write('#SBATCH  --partition=public\n')
        f.write('#SBATCH  --nodes=1\n')
        f.write('#SBATCH  --ntasks=24\n')
        f.write('#SBATCH  --ntasks-per-node=24\n')
        f.write('#SBATCH  --cpus-per-task=1\n')
        f.write('source /work/env/intel2018\n')
        f.write('srun hostname | sort | uniq >> /tmp/nodefile.$$\n')
        f.write('NP=`srun hostname | wc -l`\n')
        #f.write('cd $PBS_O_WORKDIR\n')
        #f.write('num_in=`ls -l |grep 'INCAR_' |wc -l`\n')
        f.write('for(( i=1; i<=%s; i++ ));\n' %  str(self.NumberOfLocalOptim))
        f.write('do\n')
        f.write('\tcp INCAR_$i INCAR\n')
        f.write('\tcp CONTCAR POSCAR\n')
        f.write('\t mpirun -np 24 /work/software/vasp.6.1.0/vasp_std > vasp.log_${i} 2>&1\n' % (ppn,vasppath))
        f.write('done\n')
        f.write('rm -rf CHG* WAVECAR\n')
        f.close()


if __name__ == '__main__':
    if len(sys.argv) == 1:
        a = Autocalypso()
    elif 'pbs' in sys.argv[1].lower():
        submit = 'qsub vasp.pbs'
        stat = 'qstat'
        rstat = 'qstat | grep R'
        delete = 'qdel'
        a = Autocalypso(submit,stat,rstat,delete,machine = 'pbs')
    elif 'lsf' in sys.argv[1].lower():
        submit = 'bsub -J lx.116.2fu.split < run.lsf'
        stat = 'bjobs'
        rstat = 'bjobs | grep RUN'
        delete = 'bkill'
        a = Autocalypso(submit,stat,rstat,delete,machine = 'lsf')
    elif 'yh' in sys.argv[1].lower():
        submit = 'yhbatch -p TH_NET -N 1 -n 24 vasp.sh'
        stat = 'yhqueue'
        rstat = 'yhqueue | grep   "\<R\>" '
        delete = 'yhcancel'
        a = Autocalypso(submit,stat,rstat,delete,machine = 'yh')
    elif 'slurm' in sys.argv[1].lower():
        submit = 'sbatch vasp.slurm'
        stat = 'squeue'
        rstat = 'squeue | grep   "\<R\>" '
        delete = 'scancel'
        a = Autocalypso(submit,stat,rstat,delete,machine = 'slurm')
    a.autorun()
