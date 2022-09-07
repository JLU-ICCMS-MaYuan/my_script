MODULE Constants
    !
    USE kinds, ONLY: DP, i4b
    !
    ! ... The constants needed everywhere
    !
    IMPLICIT NONE
    public
    !
    SAVE
    !
    ! ... Mathematical constants
    !
    real(DP), parameter :: pi = 3.14159265358979323846_dp
    !***************CONSTANT PARAMETERS***********************
!
    integer(i4b), parameter   :: max_struct = 10000          ! the max number of initial structures
    integer(i4b), parameter   :: max_randomstr = 5000      !
    integer(i4b), parameter   :: max_cluster = 500         ! the max number of the cluster
    integer(i4b), parameter   :: max_natom = 10000
    integer(i4b), parameter   :: max_ntyp = 10
    integer(i4b), parameter   :: max_lpso = 20
    integer(i4b), parameter   :: max_fu = 30
    integer(i4b), parameter   :: max_mol = 60
    integer(i4b), parameter   :: max_num = 5000
    integer(i4b), parameter   :: max_save = 10
    real(dp), parameter :: opoint_tol = 1d-3
    real(dp), parameter :: trans_tol = 1d-1
    !***************CONSTANT PARAMETERS***********************
    !integer(i4b),parameter   :: ilammps_reax = 12
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
    INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
    REAL(DP), PARAMETER :: PIO2 = 1.57079632679489661923132169163975144209858_DP
    REAL(DP), PARAMETER :: TWOPI = 6.283185307179586476925286766559005768394_DP
    REAL(DP), PARAMETER :: SQRT2 = 1.41421356237309504880168872420969807856967_DP
    REAL(DP), PARAMETER :: EULER = 0.5772156649015328606065120900824024310422_DP
    REAL(DP), PARAMETER :: HARTREE = 27.21138386_DP
    REAL(DP), PARAMETER :: PRE = 0.001_DP
    INTEGER, PARAMETER :: NCAL = 5
    INTEGER, PARAMETER :: MAX_HARMONIC = 10
    REAL(DP), PARAMETER :: R_CUTOFF = 3.0_DP
    INTEGER, PARAMETER :: MAX_DIS_CHECK = 1000000
    INTEGER, PARAMETER :: MAX_BOND_CHECK = 8000000
    INTEGER, PARAMETER :: MAX_TYPE_CHECK = 800000
    INTEGER, PARAMETER :: MAX_ECR_CHECK = 800000
    INTEGER, PARAMETER :: MAX_SIM_CHECK = 10000
    INTEGER, PARAMETER :: MAX_PSO_MOVE = 2000000
    REAL(DP), PARAMETER :: BAR = 100000.0_DP    !1.0 bar = 1.0E5 Pa
    REAL(DP), PARAMETER :: KBAR = 100000000.0_DP
    REAL(DP), PARAMETER :: ATM = 101325.0_DP    !1.0 atm = 101325 Pa
    INTEGER, PARAMETER :: MAX_SORT_DIM = 8 !if elecnt>this value, then use GA
    !else, use permutation
    INTEGER, PARAMETER :: CHAR_LEN = 8 !length of the atomic/molecule name
    REAL(DP), PARAMETER :: MIN_DIS_SCALE = 0.7_DP !Minimum allowed distance scale
    REAL(DP), PARAMETER :: MAX_DIS_SCALE = 1.2_DP !Maximum allowed distance scale
    REAL(DP), PARAMETER :: MAX_BOND_SCALE = 1.5_DP !Distance out this range will NOT
    !treated as "bonded"
    CHARACTER(LEN=24), PARAMETER::PYTHON = 'python'
    CHARACTER(LEN=256)  :: workdir
    INTEGER, PARAMETER :: LOGID = 5211   !For CALYPSO.log
    CHARACTER(LEN=24), PARAMETER::PARSEVER = '1.1'  !For CALYPSO.log
    INTEGER, PARAMETER :: MAX_CHAR_LEN = 10240  !maximum length of string
    INTEGER, PARAMETER :: LINE_WIDTH = 90     !Line width for output text
    REAL(DP), PARAMETER :: NOTOK_VALUE = 610612509.0_DP

!Ref http://en.wikipedia.org/wiki/Laplacian_matrix
    REAL(DP), ALLOCATABLE, DIMENSION(:, :)::lap_matrix !L=D-A

END MODULE constants
