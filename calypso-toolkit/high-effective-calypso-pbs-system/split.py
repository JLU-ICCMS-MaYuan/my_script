#!/usr/bin/python
__author__ = 'tqc'

import os
import time
import sys
import glob

class Autocalypso(object):
    def __init__(self,submit = 'qsub vasp.pbs',stat = 'qstat',rstat = 'qstat | grep R',delete = 'qdel',\
                calypath = './calypso.x',machine = 'pbs'):
        self.CalyPsoPath = calypath
        self.submit = submit
        self.stat = stat
        self.rstat = rstat
        self.delete = delete
        self.machine = machine
        self.MaxTime = 5000
        self.PopSize = 200
        self.MaxStep = 20
        self.NumberOfParallel = 6
        self.PickUp = False
        self.NumberOfLocalOptim = 3
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
                print 'There is no vasp.pbs file!!!'
                #print 'We will generate a file, and you need check the PATH of VASP and the number of node!\n'
                print 'We will generate another file, maybe you need change it!!!'
                self.writevasppbs()
                sys.exit(0)
        elif self.machine == 'lsf':
            if not os.path.exists(r'./run.lsf'):
                print 'No run.lsf file!!!'
                print 'We will generate another file, maybe you need change it!!!'
                self.writerunlsf()
                sys.exit(0)
        elif self.machine == 'yh':
            if not os.path.exists(r'./vasp.sh'):
                print 'No vasp.sh file!!!'
                print 'We will generate another file, maybe you need change it!!!'
                self.writeyh()
                sys.exit(0)
        elif self.machine == 'slurm':
            if not os.path.exists(r'./submit.sh'):
                print 'No submit.sh file!!!'
                print 'We will generate another file, maybe you need change it!!!'
                self.writeslurm()
                sys.exit(0)
        if not os.path.exists(r'./input.dat'):
            print 'No input.dat'
            sys.exit(0)
        elif not os.path.exists(r'./POTCAR'):
            print 'No POTCAR'
            sys.exit(0)
        elif len(glob.glob(r'./INCAR_*')) == 0:
            print 'No INCAR files'
            sys.exit(0)
        else:
            print 'Check files completed!!!'

    def lpickup(self):
        if not self.PickUp:
            os.system('rm step')

    def submit_vasp(self,n):
        self.f.write("%s\n" % "caly_auto_split call CALYPSO generare structure" )
        self.f.flush()
        os.system(self.CalyPsoPath)
        if n != self.GenNum:
            self.f.write("%s\n" % "caly_auto_split submit structure relax job" )
            self.f.flush()
            self.control_vasp()

    def autorun(self):
        self.checkfiles()
        self.readinput()
        self.lpickup()
        i = 0
        while i < self.GenNum:
            self.f.write("%s %s %s\n" % ("=================",str(i),"ITERATION ==================" ))
            self.f.flush()
            self.submit_vasp(i)
            i += 1

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
                os.system('cp submit.sh POTCAR INCAR_*  %s' % str(i+1))
            os.system('cp POSCAR_%s  %s/CONTCAR' % (str(i+1),str(i+1)))
        self.f.write("%s\n" % "Set structure relax jobs finished")
        totaljobid = []
        runjobid = []
        splitjobid = []
        split_run_jobid = []
        jobtime = {}
        for i in range(self.NP):
            if self.machine == 'pbs':
                id = int(os.popen(' cd %s; %s;cd ..' % (str(i+1),self.submit)).read().split('.')[0])
            elif self.machine == 'lsf':
                id = int(os.popen(' cd %s; %s;cd ..' % (str(i+1),self.submit)).read().split(' ')[1].\
                                            strip('<').strip('>'))
            elif self.machine == 'yh':
                id = int(os.popen(' cd %s; %s;cd ..' % (str(i+1),self.submit)).read().split(' ')[3])
            elif self.machine == 'slurm':
                id = int(os.popen(' cd %s; %s;cd ..' % (str(i+1),self.submit)).read().split(' ')[3])
            splitjobid.append(id)
            jobtime[id] = 0 
        self.f.write("%s\n" % "Submit structure relax jobs finished")
        #print splitjobid
        num = self.NP
        finnum = 0   
        tenode = 0
        nover = 0
        while finnum < self.StrNum:
            totaljobid = self.run_jobid(self.stat)
            runjobid = self.run_jobid(self.rstat)
            splitjobid = list(set(splitjobid).intersection(set(totaljobid)))
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
            for  ii in split_run_jobid:
                jobtime[ii] += 30
                if jobtime[ii] > self.MaxTime:
                    print 'qdel %d' % ii
                    os.system('%s %s' % (self.delete,str(ii)))  
            #print jobtime  
            enode = len(splitjobid) 
            if enode == 0:
                nover += 1
            if nover == 5:
                break
            if enode < tenode:
                finnum += (tenode-enode)
            tenode = enode
            fnode = self.NP-enode
            self.f.write("%d  %d   %d   %d\n" %(fnode,enode,self.NP,finnum))
            if fnode > 0:
                if num < self.StrNum:
                    for i in range(fnode):
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
            time.sleep(30)
            #print 'submitted',num, 'finished',finnum
        self.f.write("%s\n" % "Structure relax jobs finished")
        for i in range(self.StrNum):
            if os.path.exists(r'./%s/CONTCAR' % str(i+1)):
                os.system('cp %s/CONTCAR  CONTCAR_%s' % (str(i+1),str(i+1)))
            else:
                print i + 1,'th VASP JOB WRONG NO CONTCAR'
            if os.path.exists(r'./%s/OUTCAR' % str(i+1)):
                os.system('cp %s/OUTCAR  OUTCAR_%s' % (str(i+1),str(i+1)))
            else:
                print i + 1,'th VASP JOB WRONG NO OUTCAR'
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
            splitjobid.append(int(os.popen(' cd %s;qsub vasp.pbs;cd ..' % (str(i+1))).read().split('.')[0]))
        return splitjobid
                
    def writevasppbs(self,nodes = '1',ppn = '28',vasppath = '/public/apps/vasp.5.4.1/bin/vasp_std'):
        f = open('vasp.pbs','w')
        f.write('#!/bin/bash\n')
        f.write('#PBS -N  Li4Au2N5.1\n')
        f.write('#PBS -q  liuhy\n')
        f.write('#PBS -l nodes=%s:ppn=%s\n' % (nodes,ppn))
        f.write('#PBS -j oe\n')
        f.write('#PBS -V\n')
        f.write('ulimit -s unlimited\n')
        f.write('cd $PBS_O_WORKDIR\n')
        #f.write('num_in=`ls -l |grep 'INCAR_' |wc -l`\n')
        f.write('for(( i=1; i<=%s; i++ ));\n' %  str(self.NumberOfLocalOptim))
        f.write('do\n')
        f.write('\tcp INCAR_$i INCAR\n')
        f.write('\tcp CONTCAR POSCAR\n')
        f.write('\tmpirun -np %s  %s > vasp.log_$i\n' % (ppn,vasppath))
        f.write('done\n')
        f.close()

    def writerunlsf(self,nodes = '1',ppn = '36',vasppath = '/data/software/vasp.5.4.1/bin/vasp_std'):
        f = open('run.lsf','w')
        f.write('#!/bin/sh\n')
        f.write('#BSUB -n 36\n')
        f.write("#BSUB -R 'span[ptile=36]'\n")
        f.write("#BSUB -q mym1\n")
        f.write('for(( i=1; i<=%s; i++ ));\n' %  str(self.NumberOfLocalOptim))
        f.write('do\n')
        f.write('\tcp INCAR_$i INCAR\n')
        f.write('\tcp CONTCAR POSCAR\n')
        f.write('\tmpirun -np %s  %s > vasp.log_$i\n' % (ppn,vasppath))
        f.write('done\n')
        f.write('rm -rf CHG* WAVECAR\n')
        f.close()


    def writeyh(self):  
        f = open('vasp.sh','w')
        f.write('#!/bin/bash\n')
        f.write('for(( i=1; i<=%s; i++ ));\n' %  str(self.NumberOfLocalOptim))
        f.write('do\n')
        f.write('\tcp INCAR_$i INCAR\n')
        f.write('\tcp CONTCAR POSCAR\n')
        f.write('\tyhrun -N 1 -n 12 -p TH_NET /vol-th/home/maym/software/vasp.5.4.1/bin/vasp_std > vasp.log 2>&1\n')
        f.write('done\n')
        f.close()
            
    def writeslurm(self):
        f = open('submit.sh','w')
        f.write('#!/bin/bash\n')
        f.write('source $MODULESHOME/init/sh\n')
        f.write('module purge\n')
        f.write('module add app_env/slurm\n')
        f.write('module add app_env/vasp-5.4.1\n')
        f.write('for(( i=1; i<=%s; i++ ));\n' %  str(self.NumberOfLocalOptim))
        f.write('do\n')
        f.write('\tcp INCAR_$i INCAR\n')
        f.write('\tcp CONTCAR POSCAR\n')
        f.write('\tsrun vasp_parallel\n')
        f.write('done\n')
        f.close()

if __name__ == '__main__':
    if len(sys.argv) == 1:
        a = Autocalypso()
    elif 'pbs' in sys.argv[1].lower():
        a = Autocalypso()
    elif 'lsf' in sys.argv[1].lower():
        submit = 'bsub -J lx.116.2fu.split < run.lsf'
        stat = 'bjobs'
        rstat = 'bjobs | grep RUN'
        delete = 'bkill'
        a = Autocalypso(submit,stat,rstat,delete,machine = 'lsf')
    elif 'yh' in sys.argv[1].lower():
        submit = 'yhbatch -p TH_NET -N 1 -n 12 vasp.sh'
        stat = 'yhqueue'
        rstat = 'yhqueue | grep   "\<R\>" '
        delete = 'yhcancel'
        a = Autocalypso(submit,stat,rstat,delete,machine = 'yh')
    elif 'slurm' in sys.argv[1].lower():
        submit = 'sbatch -n 12  -p normal submit.sh'
        stat = 'squeue'
        rstat = 'squeue | grep   "\<R\>" '
        delete = 'scancel'
        a = Autocalypso(submit,stat,rstat,delete,machine = 'slurm')
    a.autorun()
