
! Copyright (C) 2002-2009 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modmain

!----------------------------!
!     lattice parameters     !
!----------------------------!
! lattice vectors stored column-wise
real(8) avec(3,3)
! inverse of lattice vector matrix
real(8) ainv(3,3)
! reciprocal lattice vectors
real(8) bvec(3,3)
! inverse of reciprocal lattice vector matrix
real(8) binv(3,3)
! unit cell volume
real(8) omega
! any vector with length less than epslat is considered zero
real(8) epslat

!--------------------------!
!     atomic variables     !
!--------------------------!
! maximum allowed species
integer, parameter :: maxspecies=8
! maximum allowed atoms per species
integer, parameter :: maxatoms=200
! number of species
integer nspecies
! number of atoms for each species
integer natoms(maxspecies)
! maximum number of atoms over all the species
integer natmmax
! total number of atoms
integer natmtot
! index to atoms and species
integer idxas(maxatoms,maxspecies)
! inverse atoms and species indices
integer idxis(maxatoms*maxspecies)
integer idxia(maxatoms*maxspecies)
! molecule is .true. is the system is an isolated molecule
logical molecule
! primcell is .true. if primitive unit cell is to be found automatically
logical primcell
! atomic positions in lattice coordinates
real(8) atposl(3,maxatoms,maxspecies)
! atomic positions in Cartesian coordinates
real(8) atposc(3,maxatoms,maxspecies)

!----------------------------------!
!     atomic species variables     !
!----------------------------------!
! species files path
character(256) sppath
! species filenames
character(256) spfname(maxspecies)
! species name
character(256) spname(maxspecies)
! species symbol
character(256) spsymb(maxspecies)
! species nuclear charge
real(8) spzn(maxspecies)
! ptnucl is .true. if the nuclei are to be treated as point charges, if .false.
! the nuclei have a finite spherical distribution
logical ptnucl
! species electronic charge
real(8) spze(maxspecies)
! species mass
real(8) spmass(maxspecies)
! smallest radial point for each species
real(8) sprmin(maxspecies)
! effective infinity for species
real(8) sprmax(maxspecies)
! number of radial points to effective infinity for each species
integer spnr(maxspecies)
! maximum spnr over all the species
integer spnrmax
! maximum allowed states for each species
integer, parameter :: maxspst=40
! number of states for each species
integer spnst(maxspecies)
! maximum spnst over all the species
integer spnstmax
! core-valence cut-off energy for species file generation
real(8) ecvcut
! semi-core-valence cut-off energy for species file generation
real(8) esccut
! state principle quantum number for each species
integer spn(maxspst,maxspecies)
! state l value for each species
integer spl(maxspst,maxspecies)
! state k value for each species
integer spk(maxspst,maxspecies)
! spcore is .true. if species state is core
logical spcore(maxspst,maxspecies)
! total number of core states
integer nstcr
! state eigenvalue for each species
real(8) speval(maxspst,maxspecies)
! state occupancy for each species
real(8) spocc(maxspst,maxspecies)
! species radial mesh
real(8), allocatable :: spr(:,:)
! species charge density
real(8), allocatable :: sprho(:,:)
! species self-consistent potential
real(8), allocatable :: spvr(:,:)

!---------------------------------------------------------!
! number of energy intervals in the DOS/optics function   !
!---------------------------------------------------------!
integer, parameter :: nwdos=500

!--------------------------!
!     timing variables     !
!--------------------------!
! initialisation
real(8) timeinit
! Hamiltonian and overlap matrix set up
real(8) timemat
! first-variational calculation
real(8) timefv
! second-variational calculation
real(8) timesv
! charge density calculation
real(8) timerho
! potential calculation
real(8) timepot
! force calculation
real(8) timefor

!-----------------------------!
!     numerical constants     !
!-----------------------------!
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: twopi=6.2831853071795864769d0
real(8), parameter :: fourpi=12.566370614359172954d0
! spherical harmonic for l=m=0
real(8), parameter :: y00=0.28209479177387814347d0
! complex constants
complex(8), parameter :: zzero=(0.d0,0.d0)
complex(8), parameter :: zone=(1.d0,0.d0)
complex(8), parameter :: zi=(0.d0,1.d0)
! array of i^l values
complex(8), allocatable :: zil(:)
! Pauli spin matrices:
! sigma_x = ( 0  1 )   sigma_y = ( 0 -i )   sigma_z = ( 1  0 )
!           ( 1  0 )             ( i  0 )             ( 0 -1 )
complex(8) sigmat(2,2,3)
data sigmat / (0.d0,0.d0), (1.d0,0.d0), (1.d0,0.d0), (0.d0,0.d0), &
              (0.d0,0.d0), (0.d0,1.d0),(0.d0,-1.d0), (0.d0,0.d0), &
              (1.d0,0.d0), (0.d0,0.d0), (0.d0,0.d0),(-1.d0,0.d0) /
! Boltzmann constant in Hartree/kelvin (CODATA 2006)
real(8), parameter :: kboltz=3.166815343d-6
! speed of light in atomic units (=1/alpha) (CODATA 2006)
real(8), parameter :: sol=137.035999679d0
! scaled speed of light
real(8) solsc
! electron g-factor (CODATA 2006)
real(8), parameter :: gfacte=2.0023193043622d0
! hartree in SI units (CODATA 2006)
real(8), parameter :: ha_si=4.35974394d-18
! Bohr radius in SI units (CODATA 2006)
real(8), parameter :: au_si=0.52917720859d-10
! Planck constant in SI units (CODATA 2006)
real(8), parameter :: hbar_si=1.054571628d-34
! electron charge in SI units (CODATA 2006)
real(8), parameter :: e_si=1.602176487d-19
! atomic unit of magnetic flux density in SI
real(8), parameter :: b_si=hbar_si/(e_si*au_si**2)
! atomic unit of time in SI
real(8), parameter :: t_si=hbar_si/ha_si
! mass of 1/12 times carbon-12 in electron masses (CODATA 2006)
real(8), parameter :: amu=1822.88848426d0

!---------------------------------!
!     miscellaneous variables     !
!---------------------------------!
! code version
integer version(3)
data version / 1,3,31 /
! maximum number of tasks
integer, parameter :: maxtasks=40
! number of tasks
integer ntasks
! task array
integer tasks(maxtasks)
! current task
integer task
! tlast is .true. if the calculation is on the last self-consistent loop
logical tlast
! number of self-consistent loops after which STATE.OUT is written
integer nwrite
! filename extension for files generated by gndstate
character(256) filext
! default file extension
data filext / '.OUT' /
! scratch space path
character(256) scrpath
! maximum number of note lines
integer, parameter :: maxnlns=20
! number of note lines
integer notelns
! notes to include in INFO.OUT
character(80) notes(maxnlns)

end module

