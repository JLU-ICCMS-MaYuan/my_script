import sys,os,re
import unicodedata
from bz2 import BZ2File
from subprocess import call


# http://docs.python.org/library/optparse.html
from optparse import OptionParser

#Run python from /software/testing/python/2.6.5-gcc-numerical/bin/
import numpy
import numpy.linalg

def inv(a):
    return numpy.array(
        [
        [(-(a[1,2]*a[2,1])+a[1,1]*a[2,2])/
         (-(a[0,2]*a[1,1]*a[2,0])+a[0,1]*a[1,2]*a[2,0]+a[0,2]*a[1,0]*a[2,1]-
           a[0,0]*a[1,2]*a[2,1]-a[0,1]*a[1,0]*a[2,2]+a[0,0]*a[1,1]*a[2,2]),
         (a[1,2]*a[2,0]-a[1,0]*a[2,2])/
         (-(a[0,2]*a[1,1]*a[2,0])+a[0,1]*a[1,2]*a[2,0]+a[0,2]*a[1,0]*a[2,1]-
           a[0,0]*a[1,2]*a[2,1]-a[0,1]*a[1,0]*a[2,2]+a[0,0]*a[1,1]*a[2,2]),
         (-(a[1,1]*a[2,0])+a[1,0]*a[2,1])/
         (-(a[0,2]*a[1,1]*a[2,0])+a[0,1]*a[1,2]*a[2,0]+a[0,2]*a[1,0]*a[2,1]-
           a[0,0]*a[1,2]*a[2,1]-a[0,1]*a[1,0]*a[2,2]+a[0,0]*a[1,1]*a[2,2])]
        ,
        [(a[0,2]*a[2,1]-a[0,1]*a[2,2])/
         (-(a[0,2]*a[1,1]*a[2,0])+a[0,1]*a[1,2]*a[2,0]+a[0,2]*a[1,0]*a[2,1]-
           a[0,0]*a[1,2]*a[2,1]-a[0,1]*a[1,0]*a[2,2]+a[0,0]*a[1,1]*a[2,2]),
         (-(a[0,2]*a[2,0])+a[0,0]*a[2,2])/
         (-(a[0,2]*a[1,1]*a[2,0])+a[0,1]*a[1,2]*a[2,0]+a[0,2]*a[1,0]*a[2,1]-
           a[0,0]*a[1,2]*a[2,1]-a[0,1]*a[1,0]*a[2,2]+a[0,0]*a[1,1]*a[2,2]),
         (a[0,1]*a[2,0]-
           a[0,0]*a[2,1])/
         (-(a[0,2]*a[1,1]*a[2,0])+a[0,1]*a[1,2]*a[2,0]+a[0,2]*a[1,0]*a[2,1]-
           a[0,0]*a[1,2]*a[2,1]-a[0,1]*a[1,0]*a[2,2]+a[0,0]*a[1,1]*a[2,2])]
        ,
        [(-(a[0,2]*a[1,1])+a[0,1]*a[1,2])/
         (-(a[0,2]*a[1,1]*a[2,0])+a[0,1]*a[1,2]*a[2,0]+a[0,2]*a[1,0]*a[2,1]-
           a[0,0]*a[1,2]*a[2,1]-a[0,1]*a[1,0]*a[2,2]+a[0,0]*a[1,1]*a[2,2]),
         (a[0,2]*a[1,0]-
           a[0,0]*a[1,2])/
         (-(a[0,2]*a[1,1]*a[2,0])+a[0,1]*a[1,2]*a[2,0]+a[0,2]*a[1,0]*a[2,1]-
           a[0,0]*a[1,2]*a[2,1]-a[0,1]*a[1,0]*a[2,2]+a[0,0]*a[1,1]*a[2,2]),
         (-(a[0,1]*a[1,0])+a[0,0]*a[1,1])/
         (-(a[0,2]*a[1,1]*a[2,0])+a[0,1]*a[1,2]*a[2,0]+a[0,2]*a[1,0]*a[2,1]-
           a[0,0]*a[1,2]*a[2,1]-a[0,1]*a[1,0]*a[2,2]+a[0,0]*a[1,1]*a[2,2])]
        ])
def mul(a,b):
    return numpy.array([a[0,0]*b[0]+a[0,1]*b[1]+a[0,2]*b[2],
                        a[1,0]*b[0]+a[1,1]*b[1]+a[1,2]*b[2],
                        a[2,0]*b[0]+a[2,1]*b[1]+a[2,2]*b[2]])

# Default value for properties
inf = float('inf')


parser = OptionParser()
parser.add_option("-n", help="Write a header line",
                  action="store_true", dest="header", default=False)
parser.add_option("-p", "--posfile", help="Specify posision output file, defalult infile.positions",
                  action="store", type="string", dest="posfile", default="infile.positions")
parser.add_option("-c", "--poscar", help="Specify a poscar file to use as base, defalult infile.ssposcar",
                  action="store", type="string", dest="poscar", default="infile.ssposcar")
parser.add_option("-f", "--forcefile", help="Specify force output file, defalult infile.forces",
                  action="store", type="string", dest="forcefile", default="infile.forces")
parser.add_option("-s", "--statfile", help="Specify stat output file defalult infile.stat",
                  action="store", type="string", dest="statfile", default="infile.stat")
parser.add_option("-m", "--metafile", help="Specify meta output file defalult infile.meta",
                  action="store", type="string", dest="metafile", default="infile.meta")
parser.add_option("--magneticfile", help="Specify magnetic output file defalult infile.magnetic",
                  action="store", type="string", dest="magneticfile", default="infile.magnetic")
parser.add_option("--magnetic", help="Write a magnetic output",
                  action="store_true", dest="magnetic", default=False)
parser.add_option("--skip", help="Skip the first N timesteps from each OUTCAR",
                  action="store", dest="skip", default=0)
parser.add_option("--use-bzip2", help="Assume that supplied OUTCARs and ssposcar are bzip2-encrypted",
                  action="store_true", dest="bzip2", default=False)
parser.add_option("--project_moments", help="Project magnetic moments ourselves instead of using VASP output",
                  action="store_true", dest="project_moments", default=False)
(options, args) = parser.parse_args()
assert len(args) > 0, "You forgot to supply at least one OUTCAR"

print( "posfile:   ", options.posfile)
print( "forcefile: ", options.forcefile)
print( "statfile:  ", options.statfile)
print("metafile:  ", options.metafile)
if(options.magnetic): print("magneticfile:  ", options.magneticfile)

rePOTIM    = re.compile(r"\s*POTIM\s*=\s*(\d+\.?\d+)\s+time-step for ionic-motion")
reNIONS    = re.compile(r".*NIONS =\s+(\d+)")
# Get the lattice layout
reLatt     = re.compile(r"\s+Lattice vectors\:\s+")
reVect     = re.compile(r"\s*A? = \(\s*(-?\d+\.\d+),\s*(-?\d+\.\d+),\s*(-?\d+\.\d+)\s*\)")
# Get the lattice layout, if isym=0 and the above wont work
reLatt2    = re.compile(r"\s+direct lattice vectors\s+reciprocal lattice vectors")
reVect2    = re.compile(r"\s+(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*.*")
#reTemp     = re.compile(r"\s*kinetic Energy EKIN\s*=\s*(\d+\.\d+)\s* \(temperature\s*(\d*.\d*) K\)")
reTemp     = re.compile(r"\s*nergy EKIN\s*=\s*(\d+\.\d+)\s*")
reTemp2    = re.compile(r"\s*lattice  EKIN_LAT=\s*(\d+\.\d+)\s* \(temperature\s*(\d*.\d*) K\)")
reEtot     = re.compile(r"\s*total energy   ETOTAL =\s*(-?\d*\.\d*)\s*eV")
reEpot     = re.compile(r"%?\s*ion-electron   TOTEN  =\s*(-?\d*\.\d*)\s*.*")
reEpot2    = re.compile(r"%?\s*free  energy   TOTEN  =\s*(-?\d*\.\d*)\s*.*")
reStress   = re.compile(r"\s*in kB (.*)")
rePressure = re.compile(r"\s*external pressure =\s*(-?\d*\.\d*)\s*kB  Pullay stress =.*")
rePosition = re.compile(r"\s*POSITION\s*TOTAL-FORCE.*")
reTEBEG    = re.compile(r"\s*TEBEG\s*=\s*(\d+\.?\d+);\s*TEEND\s*=\s*(\d+\.\d+)\s*temperature during run")
reSigma    = re.compile(r"\s*ISMEAR\s*=\s*(\+?-?\d*);\s*SIGMA\s*=\s*(\d*\.\d*)\s*broadening.*")
reIbrion   = re.compile(r"\s*IBRION\s*=\s*(\+?-?\d*)\s*")
reMag      = re.compile(r"\s*magnetization \(x\)\s*")

potim=inf
nions=inf
tebeg=inf
ismear=inf
sigma=inf

lattice = numpy.array([
    [inf, inf, inf],
    [inf, inf, inf],
    [inf, inf, inf],
    ])

posfile  = open(options.posfile,'w')
forcefile  = open(options.forcefile,'w')
statfile = open(options.statfile,'w')
statfileplot = open(options.statfile+".gnuplot",'w')
metafile = open(options.metafile,'w')
if(options.magnetic): magneticfile = open(options.magneticfile,'w')

## read base from a POSCAR if availible
try:
    if options.bzip2:
        if options.poscar.endswith('bz2'):
            poscar = BZ2File(options.poscar, 'r')
        else:
            poscar = BZ2File(options.poscar + ".bz2", 'r')
    else:
        poscar = open(options.poscar, 'r')

    poscar.readline() # comment
    a = float(poscar.readline()) # latticeconst
    if a < 0:
        print("Damit don't use volume... in poscar, fix me...")
        sys.exit(-1)
    lattice[0] = list(map(float,poscar.readline().strip().split()))
    lattice[1] = list(map(float,poscar.readline().strip().split()))
    lattice[2] = list(map(float,poscar.readline().strip().split()))
    lattice=lattice*a

    poscar.close()

    reciprocal_lattice = inv(lattice)
    print("lattice from "+options.poscar+":")
    print(lattice)
    print("reciprocal lattice:")
    print(reciprocal_lattice)
        
except IOError:
    print("Could not find "+options.poscar+". Will use the basis from"+
    "the OUTCAR, note that this might not be the same as the one in the POSCAR")
    
if(options.header):
    statfile.write("Step, Time [fs], Energy [eV] tot, pot, kin, Temp [K], Pressure [GPa], Stress [GPa] x,y,z,xy,zy,zx\n")

if options.bzip2:
    f = BZ2File(args[0], 'r')
else:
    f = open(args[0],'r')
while potim == inf or nions == inf or lattice[1,1] == inf or tebeg == inf or sigma == inf or ismear == inf:
    line = f.readline()
    if not line:
        print("ERROR: Can't parse all data out of OUTCAR")
        sys.exit(-1)

    if(rePOTIM.match(line)):
        m = rePOTIM.search(line)
        potim = float( m.group(1) )
        print("POTIM: ",potim)
    
    if(reNIONS.match(line)):
        m = reNIONS.search(line)
        nions = int( m.group(1) )
        print("NIONS: ", nions)

    if(reTEBEG.match(line)):
        m = reTEBEG.search(line)
        tebeg = float( m.group(1) )
        print("TEBEG: ", tebeg)

    if(reSigma.match(line)):
        m = reSigma.search(line)
        ismear = int( m.group(1) )
        sigma = float( m.group(2) )
        print("ISMEAR: ", ismear)
        print("SIGMA: ", sigma)

    if(reIbrion.match(line)):
        m = reIbrion.search(line)
        ibrion = m.group(1)
        print("IBRION: ", ibrion)

    if(lattice[1,1]==inf and reLatt.match(line)):
        f.readline()
        
        for i in range(0,3):
            line = f.readline()
            m = reVect.search(line)
            lattice[i] =[float(m.group(1)),float(m.group(2)),float(m.group(3))]

        reciprocal_lattice = inv(lattice)
        print("lattice:")
        print(lattice)
        print("reciprocal lattice:")
        print(reciprocal_lattice)

    if(lattice[1,1]==inf and reLatt2.match(line)):
        for i in range(0,3):
            line = f.readline()
            m = reVect2.search(line)
            lattice[i] =[float(m.group(1)),float(m.group(2)),float(m.group(3))]

        reciprocal_lattice = inv(lattice)
        print("lattice:")
        print(lattice)
        print("reciprocal lattice:")
        print(reciprocal_lattice)
        
        
f.close()

skip=int(options.skip)
step = 0
magtmp="" # skip last mag printout
filecount = 0

call(["rm", "-f","infile.projected_moments"]) #remove it even if option not set (so can't have mismatch)

for file in args:
    filecount += 1
    if options.project_moments:
        mag_dir=os.path.dirname(file)
        mag_out=os.path.join(mag_dir, "outfile.magnetic_moments")
        if not os.path.isfile(mag_out):
            print(mag_out+" does NOT exist so making it:")
            print("[tip: consider running project_magnetic_moments manually (e.g. right after the DFT calculation) because then you can run all configurations independently in parallel, as opposed to serially here in this script.  ...And you can customize your settings]")
            mag_ss=os.path.join(mag_dir,"POSCAR")
            mag_ch=os.path.join(mag_dir,"CHGCAR")
            if not (os.path.isfile(mag_ss) and os.path.isfile(mag_ch)):
                print("cannot find POSCAR and CHGCAR in "+mag_dir)
                print("put them named as such there (with the OUTCAR) or run project_magnetic_moments there yourself")
                print("exiting.  FAILED!!")
                sys.exit()

            call(["project_magnetic_moments", "--infile_ssposcar", mag_ss, "--infile_chgcar", mag_ch, "--outfile_mag", mag_out ])

        assert os.path.isfile(mag_out)
        os.system("cat " + mag_out+ " >> infile.projected_moments")

    if filecount > 1:
        if skip > 0:
            print("No skipped steps for the next files")
            skip = 0

    filestep = 0
    if options.bzip2:
        f = BZ2File(file, 'r')
    else:
        f = open(file,'r')
    sys.stdout.write(file+":")
    sys.stdout.flush()
    magtmp=""
    while 1:
        line=f.readline()
        if not line: break
        
        # matcha i stess blocket
        if(reStress.match(line)):
            m = reStress.search(line)
            stresses = list(map(float,m.group(1).split()))
            line=f.readline()
            m = rePressure.search(line)
            pressure = float(m.group(1))
            # plussa
            filestep+=1

        elif(reEpot2.match(line)):
            if(ibrion!="0"):
                m = reEpot2.search(line)
                epot = float( m.group(1) )
                etot = epot
                temp = 0
                ekin = 0
                step+=1
                statfile.write(("{0:<6d} {1:<8.2f}  {2: .6f}  {3: .6f}  {4: .6f}  {5: .2f}  {6: .3f} "+
                                "{7: .3f}  {8: .3f}  {9: .3f}  {10: .3f}  {11: .3f}  {12: .3f}\n").format(
                    step,(step-1)*potim,etot,epot,ekin,temp,0.1*pressure,
                    0.1*stresses[0],0.1*stresses[1],0.1*stresses[2],0.1*stresses[3],0.1*stresses[4],0.1*stresses[5])
                               )

        # matcha i energi blocket
        elif(reEpot.match(line)):
            m = reEpot.search(line)
            epot = float( m.group(1) )
            line=f.readline()  # Ekin
            m = reTemp.search(line)
            ekin = float( m.group(1) )
            # temp = float( m.group(2) )
            # temp = float( 1.0 )
            line=f.readline()  # kin. lattice ...
            m = reTemp2.search(line)
            temp = float( m.group(2) )
            line=f.readline()  # nose...
            line=f.readline()  # nose...
            line=f.readline()  # ----
            line=f.readline()  # etotal
            m = reEtot.search(line)
            etot = float( m.group(1) )
            # only write if we have skipped enough steps
            if ( filestep > skip ):
                step += 1
                if( (step % 100 ) == 0 ):
                    sys.stdout.write(".")
                    sys.stdout.flush()
                statfile.write(("{0:<6d} {1:<8.2f}  {2: .6f}  {3: .6f}  {4: .6f}  {5: .2f}  {6: .3f} "+
                                "{7: .3f}  {8: .3f}  {9: .3f}  {10: .3f}  {11: .3f}  {12: .3f}\n").format(
                    step,(step-1)*potim,etot,epot,ekin,temp,0.1*pressure,
                    0.1*stresses[0],0.1*stresses[1],0.1*stresses[2],0.1*stresses[3],0.1*stresses[4],0.1*stresses[5])
                               )
                if(options.magnetic):
                    magneticfile.write(magtmp)
                    magtmp=""
            
        elif(rePosition.match(line)):
            f.readline()
            for i in range(nions):
                line = list(map(float,f.readline().strip().split()))
                line[0:3] = mul(reciprocal_lattice, line[0:3])
                if ( filestep > skip ):
                    posfile.write("{0: .11f} {1: .11f} {2: .11f}\n".format(line[0], line[1], line[2]))
                    forcefile.write("{0: .8f} {1: .8f} {2: .8f}\n".format(line[3],line[4],line[5]))


        elif(options.magnetic and reMag.match(line)):
            f.readline()  #
            f.readline()  # # of ion     s       p       d       tot
            f.readline()  # ----------------------------------------
            for i in range(nions):
                magtmp+=f.readline()
            
    f.close()
    if(skip>0):
        print(" "+str(filestep) + " steps read, skipped the first "+str(skip))
    else:
        print(" "+str(filestep) + " steps read")

metafile.write(str(nions)+"\n")
metafile.write(str(step)+"\n")
metafile.write(str(potim)+"\n")
metafile.write(str(tebeg)+"\n")

statfileplot.write("dir=system(\"pwd\")\n")
statfileplot.write("set terminal x11 size 1200,800 font \"CMU Serif,10\"\n")
statfileplot.write("#Uncomment the following two lines to save a png file of the plot\n")
statfileplot.write("#set terminal pngcairo size 1200,800 font \",12\"\n")
statfileplot.write("#set out \"statplot.png\"\n")
statfileplot.write("set fit logfile system(\"mktemp\")\n")
statfileplot.write("set multiplot\n")
statfileplot.write(" set border lw 0.5 \n")
statfileplot.write(" set ytics scale 0.5 \n")
statfileplot.write(" set xtics scale 0.5 \n")
statfileplot.write(" set mxtics 5 \n")
statfileplot.write(" set mytics 5 \n")
statfileplot.write("  \n")
statfileplot.write(" set style line 1 lc rgb \"#65CAAE\" lw 1 \n") 
statfileplot.write(" set style line 2 lc rgb \"#CC8F12\" lw 1 \n")
statfileplot.write(" set style line 3 lc rgb \"#88172B\" lw 1 \n")
statfileplot.write(" set style line 4 lc rgb \"#90AB3F\" lw 1 \n")
statfileplot.write(" set style line 5 lc rgb \"#FDF488\" lw 1 \n")
statfileplot.write(" set style line 6 lc rgb \"#528E7C\" lw 1 \n")
statfileplot.write(" set style increment user \n")
statfileplot.write("set title font \"CMU Serif,14\"\n")
statfileplot.write("set xlabel \"Time (ps)\"\n")
statfileplot.write("set size 1.0,1.0\n")
statfileplot.write("f2p(n) = real(n)/1000.0\n")

statfileplot.write("set origin 0.0, 0.0\n")
statfileplot.write("set size 0.5, 0.5\n")
statfileplot.write("set title \"Total Energy\"\n")
statfileplot.write("set ylabel \"Energy (eV)\"\n")
statfileplot.write("fe(x) = ae *x + be\n")
statfileplot.write("fit fe(x) \"infile.stat\" u (f2p($2)):3 via ae, be\n")
statfileplot.write("plot \"infile.stat\" u (f2p($2)):3 t \"Energy\" w line, \\\n")
statfileplot.write("fe(x) t sprintf(\"fit: %f*t + %f\", ae,be)\n")

statfileplot.write("set origin 0.0, 0.5\n")
statfileplot.write("set size 0.5, 0.5\n")
statfileplot.write("set title \"Temperature\"\n")
statfileplot.write("set ylabel \"Temperature (K)\"\n")
statfileplot.write("ft(x) = at *x + bt\n")
statfileplot.write("fit ft(x) \"infile.stat\" u (f2p($2)):6 via at, bt\n")
statfileplot.write("plot \"infile.stat\" u (f2p($2)):6 t \"Temperature\" w line, \\\n")
statfileplot.write("ft(x) t sprintf(\"fit: %f*t + %f\", at,bt)\n")

statfileplot.write("set origin 0.5, 0.0\n") 
statfileplot.write("set size 0.5, 0.5\n")
statfileplot.write("set title \"Pressure\"\n")
statfileplot.write("set ylabel \"Pressure (GPa)\"\n")
statfileplot.write("fp(x) = ap *x + bp\n")
statfileplot.write("fit fp(x) \"infile.stat\" u (f2p($2)):7 via ap, bp\n")
statfileplot.write("plot \"infile.stat\" u (f2p($2)):7 t \"Pressure\" w line,\\\n")
statfileplot.write("fp(x) t sprintf(\"fit: %f*t + %f\", ap,bp)\n")

statfileplot.write("set origin 0.5, 0.5\n")
statfileplot.write("set size 0.5, 0.5\n")
statfileplot.write("set title \"Stress\"\n")
statfileplot.write("set ylabel \"Stress (GPa)\"\n")
statfileplot.write("plot \"infile.stat\" u (f2p($2)):8 t \"x\" w line,\\\n")
statfileplot.write("\"infile.stat\" u (f2p($2)):9  t \"y\"  w line,\\\n")
statfileplot.write("\"infile.stat\" u (f2p($2)):10 t \"z\"  w line,\\\n") 
statfileplot.write("\"infile.stat\" u (f2p($2)):11 t \"xy\" w line,\\\n")
statfileplot.write("\"infile.stat\" u (f2p($2)):12 t \"zy\" w line,\\\n")
statfileplot.write("\"infile.stat\" u (f2p($2)):13 t \"zx\" w line\n") 
statfileplot.write("unset multiplot\n") 

posfile.close()
forcefile.close()
metafile.close()
statfile.close()
statfileplot.close()
if(options.magnetic): magneticfile.close()

if options.project_moments:
    # make sure have same # of lines:
    num_lines_pos = sum(1 for line in open(options.posfile))
    num_lines_mag = sum(1 for line in open("infile.projected_moments"))
    if not num_lines_pos == num_lines_mag:
        print("ERROR: mismatch between magnetic file and infile.positions. DONT TRUST OUTPUT.")
        print("exiting.  FAILED!!")
        sys.exit()
