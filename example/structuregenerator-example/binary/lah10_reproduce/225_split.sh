[splitwps]
spacegroup_number = 225
nameofatoms = ["Ar", "H"]
wyckoffpositions = {'4a':  False,
                    '4b':  False,
                    '8c':  False,
                    '24d': False,
                    '24e': True,
                    '32f': True,
                    '48g': True,
                    '48h': True,
                    '48i': True,
                    '96j': True,
                    '96k': True,
                    '192l': True,
                    }
nonH_upper_limit = '8c'
H_lower_limit    = '8c'
sitesoccupiedrange=[[1,1],
                    [1,3],]
popsize=300
maxlimit=150
distancematrix=[[2.014, 1.590],
                [1.590, 1.116],]
clathrate_ratio=1.0
[pso]
numberOflbest = 4
simthreshold = 0.06
fingerprint = "bcm"
lbest = 4
critic = "enthalpy"
maxstep= 50
pso_ltype=["cubic"]
pso_ratio=0.5
