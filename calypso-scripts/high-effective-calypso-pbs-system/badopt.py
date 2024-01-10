#!/data/home/mym/soft/anaconda/ana/bin/python
import os
import glob
import time
import math
while True:
    time.sleep(5)
    WrongOpt= os.popen(''' grep '\*\*\*\*\*\*\*' OUTCAR ''').read()
    if (WrongOpt):
      os.system('killall -9 vasp_std')
      print 'Optgoodkill'
#      break
   
