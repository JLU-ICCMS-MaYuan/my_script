# irrep 使用文档 

## irrep --help 帮助手册
```shell
irrep --help
Usage: irrep [OPTIONS]

              # ###   ###   #####  ###
              # #  #  #  #  #      #  #
              # ###   ###   ###    ###
              # #  #  #  #  #      #
              # #   # #   # #####  #

  version 2.3.0
  
  Calculates the expectation values of symmetry operations
  <Psi_nk | T(g) | Psi_nk >  as well as irreducible representations,
  Wannier charge centers (1D) Zak phases and many more.

  Examples:
  irrep -Ecut=50 -code=abinit -fWFK=Bi_WFK -refUC=0,-1,1,1,0,-1,-1,-1,-1 -kpoints=11 -IBend=5 -kpnames="GM"
  irrep -Ecut=50 -code=espresso -prefix=Bi -refUC=0,-1,1,1,0,-1,-1,-1,-1 -kpoints=11 -IBend=5 -kpnames="GM"

  If you have interest in the code and not sure how to use it, do not hesitate
  to contact the author:

  > Stepan S. Tsirkin
  > stepan.tsirkin@epfl.ch

Options:
  -Ecut FLOAT                     Energy cut-off in eV used in the
                                  calculation. A value of 50 eV is
                                  recommended. If not set, will default to the
                                  cut-off used in the DFT calculation.
  -correct_Ecut0 FLOAT            In case of VASP, if you get an error like '
                                  computed ncnt=*** != input nplane=*** ', try
                                  to set this parameter to a small positive or
                                  negative value (usually of order  +- 1e-7)
  -fWAV TEXT                      Filename for wavefunction in VASP WAVECAR
                                  format. Only used if code is "vasp".
  -fPOS TEXT                      Filename for wavefunction in VASP POSCAR
                                  format. Only used if code is "vasp".
  -fWFK TEXT                      Filename for wavefunction in ABINIT WFK
                                  format. Only used if code is "abinit".
  -gpaw_calc TEXT                 Filename for gpaw calculator. Only used if
                                  code is "gpaw".
  -prefix TEXT                    Prefix used for Quantum Espresso
                                  calculations (data should be in prefix.save)
                                  or seedname of Wannier90 files.
  -from_sym_file TEXT             if present, the symmetry operations will be
                                  read from this file, instead of those
                                  computed by spglib
  -IBstart INTEGER                The first band to be considered. If <= 0
                                  starting from the lowest band (count from
                                  one).
  -IBend INTEGER                  The last band to be considered. If <=0 up to
                                  the highest band (count from one).
  -code [vasp|abinit|espresso|wannier90|gpaw]
                                  Set which electronic structure code to
                                  interface with. If using ABINIT, always use
                                  "istwfk=1".
  -spinor                         Indicate the wavefunctions are spinor. Only
                                  used if code is "vasp".
  -kpoints TEXT                   Comma-separated list of k-point indices
                                  (starting from 1).
  -kpnames TEXT                   Comma-separated list of k-point names (as in
                                  the tables) with one entry per each value in
                                  the k-points list. Important! K-points is
                                  assumed to be an ordered list!
  -refUC TEXT                     The lattice vectors of the reference unit
                                  cell (as given in the crystallographic
                                  tables) expressed in terms of the unit cell
                                  vectors used in the calculation. Nine comma-
                                  separated numbers.
  -shiftUC TEXT                   The vector to shift the calculated unit cell
                                  origin (in units of the calculated lattice),
                                  to get the unit cell as defined in
                                  crystallographic tables. Three comma-
                                  separated numbers.
  -isymsep TEXT                   Index of the symmetry to separate the
                                  eigenstates. with new method works with any
                                  code/pseudopotentialPreviously worked well
                                  only for norm-conserving potentials.
  -onlysym                        Only calculate the symmetry operations
  -unk_formatted                  expect UNK files to be formatted (only
                                  relevant when -code=wannier90 )
  -writesym                       write the prefix.sym file needed for the
                                  Wannier90 sitesym calculations
  -alat FLOAT                     for writesym - the alat parameter. For QE,
                                  it is read from the prefix.save/data-file-
                                  schema.xmlFor other codes needs to be
                                  provided. (To be honest, the .sym file is
                                  useless for other codes for now, but still
                                  ..)
  -ZAK                            Calculate Zak phase
  -WCC                            Calculate Wannier charge centres
  -plotbands                      Write gnuplottable files with all symmetry
                                  eigenvalues
  -EF TEXT                        Fermi energy to shift energy-levels.
                                  Default: 0.0. If it is set to a number, this
                                  value will be used to shift the energy-
                                  levels. If it is set to 'auto', the code
                                  will try to parse it from DFT files and set
                                  it to 0.0 if it could not do so.
  -degenThresh FLOAT              Threshold to decide whether energy-levels
                                  are degenerate. Default: 1e-4
  -symprec FLOAT                  Symmetry precision. Default: 1e-5. Passed to
                                  spglib get_symmetry
  -groupKramers                   Group wave-functions in pairs of Kramers.
                                  Default: False.
  -symmetries TEXT                Indices of symmetries to be printed.
                                  Default: all detected symmetries.
  -suffix TEXT                    Suffix to name files containing data for
                                  band plotting. Default: tognuplot
  -config PATH                    Define irrep inputs from a configuration
                                  file in YAML or JSON format.
  -searchcell                     Find transformation to conventional unit
                                  cell. If it is not specified, the
                                  transformation to the conventional unit cell
                                  will not be calculated and, if refUC or
                                  shiftUC were given, it will not be checked
                                  if they correctly lead to the conventional
                                  cell.Default: False, unless -kpnames was
                                  specified in CLI
  -trans_thresh FLOAT             Threshold to compare translational parts of
                                  symmetries.Default: 1e-5
  -magmom TEXT                    Path to magnetic moments' file. One row per
                                  atom, three coordinates in Cartesian.
  -v                              Verbosity flag. -vv: very detailed info,
                                  recommended when you get an error. -v
                                  (default for CLI): info about some
                                  decissions taken internaly through the
                                  code's execution, recommended when the code
                                  runs without errors, but the result is not
                                  what you expected. If you don't set this
                                  tag, you will get the basic info.
  -json_file TEXT                 File to save the output in JSON format.
                                  (without extension, the '.json' will be
                                  added automatically)
  --time-reversal                 Consider TRS a symmetry and use the
                                  corresponding gray group.
  --print-hs-kpoints              Print high-symmetry k-points in the
                                  calculation and reference cell.
  --symmetry-indicators           Compute symmetry indicators if they are non-
                                  trivial. Irreps must be identified in the
                                  process.
  --ebr-decomposition             Compute the EBR decomposition and
                                  topological classification according to TQC.
                                  Irreps must be identified in the process.
  --version                       Show the version and exit.
  --help                          Show this message and exit.
```

## irrep使用帮助
```shell
# 以编号得形式指定被研究的k点，此时输出的结果关于irreps不包含不可约表示名称的信息
irrep -Ecut=900 -code=vasp -kpoints=1 -EF=11.8734

# 以名称得形式指定被研究的k点，此时输出的结果关于irreps包含不可约表示名称的信息
# 一般来说kpoints和kpnames要一起配合使用，比如GM点的编号为1就写为-kpoints=1  -kpnames=GM
irrep -Ecut=900 -code=vasp  -kpoints=1  -kpnames=GM -EF=11.8734

```

## 报错集锦

### 报错1 `RuntimeError: the kpoint [0.5 0.  0. ] does not correspond to the point GM ([0. 0. 0.] in refUC / [0. 0. 0.] in primUC) in the table`

这是由于指定的-kpoints和-kpnames不匹配导致的。因为-kpnames=GM对应的-kpoints是1而不是100

一般来说kpoints的设置严格对应与KPOINTS中各个高对称点的顺序。

**比如在如下的KPOINTS中，-kpnames=GM对应-kpoints=100**
```shell
K-Path Generated by VASPKIT.
   100
Line-Mode
Reciprocal
   0.5000000000   0.2500000000   0.7500000000     W              
   0.5000000000   0.5000000000   0.5000000000     L              
 
   0.5000000000   0.5000000000   0.5000000000     L              
   0.0000000000   0.0000000000   0.0000000000     GAMMA          
 
   0.0000000000   0.0000000000   0.0000000000     GAMMA          
   0.5000000000   0.0000000000   0.5000000000     X              
 
   0.5000000000   0.0000000000   0.5000000000     X              
   0.5000000000   0.2500000000   0.7500000000     W              
 
   0.5000000000   0.2500000000   0.7500000000     W              
   0.3750000000   0.3750000000   0.7500000000     K          
```
**比如在如下的KPOINTS中，-kpnames=GM对应-kpoints=1**
```shell
KPATH
100
Line-Mode
Reciprocal
0.00000000 0.00000000 0.00000000  $G$
0.50000000 0.00000000 0.00000000  $M$

0.50000000 0.00000000 0.00000000  $M$
0.33333333 0.33333333 0.00000000  $K$

0.33333333 0.33333333 0.00000000  $K$
0.00000000 0.00000000 0.00000000  $G$

0.00000000 0.00000000 0.00000000  $G$
0.00000000 0.00000000 0.50000000  $A$

0.00000000 0.00000000 0.50000000  $A$
0.50000000 0.00000000 0.50000000  $L$

0.50000000 0.00000000 0.50000000  $L$
0.33333333 0.33333333 0.50000000  $H$

0.33333333 0.33333333 0.50000000  $H$
0.00000000 0.00000000 0.50000000  $A$
```