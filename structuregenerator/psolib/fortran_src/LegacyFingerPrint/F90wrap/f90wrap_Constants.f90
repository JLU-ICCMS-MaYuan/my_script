! Module constants defined in file ./FPP/Constants.fpp

subroutine f90wrap_constants__get__pi(f90wrap_pi)
    use constants, only: constants_pi => pi
    implicit none
    real(8), intent(out) :: f90wrap_pi
    
    f90wrap_pi = constants_pi
end subroutine f90wrap_constants__get__pi

subroutine f90wrap_constants__get__max_struct(f90wrap_max_struct)
    use constants, only: constants_max_struct => max_struct
    implicit none
    integer(4), intent(out) :: f90wrap_max_struct
    
    f90wrap_max_struct = constants_max_struct
end subroutine f90wrap_constants__get__max_struct

subroutine f90wrap_constants__get__max_randomstr(f90wrap_max_randomstr)
    use constants, only: constants_max_randomstr => max_randomstr
    implicit none
    integer(4), intent(out) :: f90wrap_max_randomstr
    
    f90wrap_max_randomstr = constants_max_randomstr
end subroutine f90wrap_constants__get__max_randomstr

subroutine f90wrap_constants__get__max_cluster(f90wrap_max_cluster)
    use constants, only: constants_max_cluster => max_cluster
    implicit none
    integer(4), intent(out) :: f90wrap_max_cluster
    
    f90wrap_max_cluster = constants_max_cluster
end subroutine f90wrap_constants__get__max_cluster

subroutine f90wrap_constants__get__max_natom(f90wrap_max_natom)
    use constants, only: constants_max_natom => max_natom
    implicit none
    integer(4), intent(out) :: f90wrap_max_natom
    
    f90wrap_max_natom = constants_max_natom
end subroutine f90wrap_constants__get__max_natom

subroutine f90wrap_constants__get__max_ntyp(f90wrap_max_ntyp)
    use constants, only: constants_max_ntyp => max_ntyp
    implicit none
    integer(4), intent(out) :: f90wrap_max_ntyp
    
    f90wrap_max_ntyp = constants_max_ntyp
end subroutine f90wrap_constants__get__max_ntyp

subroutine f90wrap_constants__get__max_lpso(f90wrap_max_lpso)
    use constants, only: constants_max_lpso => max_lpso
    implicit none
    integer(4), intent(out) :: f90wrap_max_lpso
    
    f90wrap_max_lpso = constants_max_lpso
end subroutine f90wrap_constants__get__max_lpso

subroutine f90wrap_constants__get__max_fu(f90wrap_max_fu)
    use constants, only: constants_max_fu => max_fu
    implicit none
    integer(4), intent(out) :: f90wrap_max_fu
    
    f90wrap_max_fu = constants_max_fu
end subroutine f90wrap_constants__get__max_fu

subroutine f90wrap_constants__get__max_mol(f90wrap_max_mol)
    use constants, only: constants_max_mol => max_mol
    implicit none
    integer(4), intent(out) :: f90wrap_max_mol
    
    f90wrap_max_mol = constants_max_mol
end subroutine f90wrap_constants__get__max_mol

subroutine f90wrap_constants__get__max_num(f90wrap_max_num)
    use constants, only: constants_max_num => max_num
    implicit none
    integer(4), intent(out) :: f90wrap_max_num
    
    f90wrap_max_num = constants_max_num
end subroutine f90wrap_constants__get__max_num

subroutine f90wrap_constants__get__max_save(f90wrap_max_save)
    use constants, only: constants_max_save => max_save
    implicit none
    integer(4), intent(out) :: f90wrap_max_save
    
    f90wrap_max_save = constants_max_save
end subroutine f90wrap_constants__get__max_save

subroutine f90wrap_constants__get__opoint_tol(f90wrap_opoint_tol)
    use constants, only: constants_opoint_tol => opoint_tol
    implicit none
    real(8), intent(out) :: f90wrap_opoint_tol
    
    f90wrap_opoint_tol = constants_opoint_tol
end subroutine f90wrap_constants__get__opoint_tol

subroutine f90wrap_constants__get__trans_tol(f90wrap_trans_tol)
    use constants, only: constants_trans_tol => trans_tol
    implicit none
    real(8), intent(out) :: f90wrap_trans_tol
    
    f90wrap_trans_tol = constants_trans_tol
end subroutine f90wrap_constants__get__trans_tol

subroutine f90wrap_constants__get__I2B(f90wrap_I2B)
    use constants, only: constants_I2B => I2B
    implicit none
    integer, intent(out) :: f90wrap_I2B
    
    f90wrap_I2B = constants_I2B
end subroutine f90wrap_constants__get__I2B

subroutine f90wrap_constants__get__I1B(f90wrap_I1B)
    use constants, only: constants_I1B => I1B
    implicit none
    integer, intent(out) :: f90wrap_I1B
    
    f90wrap_I1B = constants_I1B
end subroutine f90wrap_constants__get__I1B

subroutine f90wrap_constants__get__PIO2(f90wrap_PIO2)
    use constants, only: constants_PIO2 => PIO2
    implicit none
    real(8), intent(out) :: f90wrap_PIO2
    
    f90wrap_PIO2 = constants_PIO2
end subroutine f90wrap_constants__get__PIO2

subroutine f90wrap_constants__get__TWOPI(f90wrap_TWOPI)
    use constants, only: constants_TWOPI => TWOPI
    implicit none
    real(8), intent(out) :: f90wrap_TWOPI
    
    f90wrap_TWOPI = constants_TWOPI
end subroutine f90wrap_constants__get__TWOPI

subroutine f90wrap_constants__get__SQRT2(f90wrap_SQRT2)
    use constants, only: constants_SQRT2 => SQRT2
    implicit none
    real(8), intent(out) :: f90wrap_SQRT2
    
    f90wrap_SQRT2 = constants_SQRT2
end subroutine f90wrap_constants__get__SQRT2

subroutine f90wrap_constants__get__EULER(f90wrap_EULER)
    use constants, only: constants_EULER => EULER
    implicit none
    real(8), intent(out) :: f90wrap_EULER
    
    f90wrap_EULER = constants_EULER
end subroutine f90wrap_constants__get__EULER

subroutine f90wrap_constants__get__HARTREE(f90wrap_HARTREE)
    use constants, only: constants_HARTREE => HARTREE
    implicit none
    real(8), intent(out) :: f90wrap_HARTREE
    
    f90wrap_HARTREE = constants_HARTREE
end subroutine f90wrap_constants__get__HARTREE

subroutine f90wrap_constants__get__PRE(f90wrap_PRE)
    use constants, only: constants_PRE => PRE
    implicit none
    real(8), intent(out) :: f90wrap_PRE
    
    f90wrap_PRE = constants_PRE
end subroutine f90wrap_constants__get__PRE

subroutine f90wrap_constants__get__NCAL(f90wrap_NCAL)
    use constants, only: constants_NCAL => NCAL
    implicit none
    integer, intent(out) :: f90wrap_NCAL
    
    f90wrap_NCAL = constants_NCAL
end subroutine f90wrap_constants__get__NCAL

subroutine f90wrap_constants__get__MAX_HARMONIC(f90wrap_MAX_HARMONIC)
    use constants, only: constants_MAX_HARMONIC => MAX_HARMONIC
    implicit none
    integer, intent(out) :: f90wrap_MAX_HARMONIC
    
    f90wrap_MAX_HARMONIC = constants_MAX_HARMONIC
end subroutine f90wrap_constants__get__MAX_HARMONIC

subroutine f90wrap_constants__get__R_CUTOFF(f90wrap_R_CUTOFF)
    use constants, only: constants_R_CUTOFF => R_CUTOFF
    implicit none
    real(8), intent(out) :: f90wrap_R_CUTOFF
    
    f90wrap_R_CUTOFF = constants_R_CUTOFF
end subroutine f90wrap_constants__get__R_CUTOFF

subroutine f90wrap_constants__get__MAX_DIS_CHECK(f90wrap_MAX_DIS_CHECK)
    use constants, only: constants_MAX_DIS_CHECK => MAX_DIS_CHECK
    implicit none
    integer, intent(out) :: f90wrap_MAX_DIS_CHECK
    
    f90wrap_MAX_DIS_CHECK = constants_MAX_DIS_CHECK
end subroutine f90wrap_constants__get__MAX_DIS_CHECK

subroutine f90wrap_constants__get__MAX_BOND_CHECK(f90wrap_MAX_BOND_CHECK)
    use constants, only: constants_MAX_BOND_CHECK => MAX_BOND_CHECK
    implicit none
    integer, intent(out) :: f90wrap_MAX_BOND_CHECK
    
    f90wrap_MAX_BOND_CHECK = constants_MAX_BOND_CHECK
end subroutine f90wrap_constants__get__MAX_BOND_CHECK

subroutine f90wrap_constants__get__MAX_TYPE_CHECK(f90wrap_MAX_TYPE_CHECK)
    use constants, only: constants_MAX_TYPE_CHECK => MAX_TYPE_CHECK
    implicit none
    integer, intent(out) :: f90wrap_MAX_TYPE_CHECK
    
    f90wrap_MAX_TYPE_CHECK = constants_MAX_TYPE_CHECK
end subroutine f90wrap_constants__get__MAX_TYPE_CHECK

subroutine f90wrap_constants__get__MAX_ECR_CHECK(f90wrap_MAX_ECR_CHECK)
    use constants, only: constants_MAX_ECR_CHECK => MAX_ECR_CHECK
    implicit none
    integer, intent(out) :: f90wrap_MAX_ECR_CHECK
    
    f90wrap_MAX_ECR_CHECK = constants_MAX_ECR_CHECK
end subroutine f90wrap_constants__get__MAX_ECR_CHECK

subroutine f90wrap_constants__get__MAX_SIM_CHECK(f90wrap_MAX_SIM_CHECK)
    use constants, only: constants_MAX_SIM_CHECK => MAX_SIM_CHECK
    implicit none
    integer, intent(out) :: f90wrap_MAX_SIM_CHECK
    
    f90wrap_MAX_SIM_CHECK = constants_MAX_SIM_CHECK
end subroutine f90wrap_constants__get__MAX_SIM_CHECK

subroutine f90wrap_constants__get__MAX_PSO_MOVE(f90wrap_MAX_PSO_MOVE)
    use constants, only: constants_MAX_PSO_MOVE => MAX_PSO_MOVE
    implicit none
    integer, intent(out) :: f90wrap_MAX_PSO_MOVE
    
    f90wrap_MAX_PSO_MOVE = constants_MAX_PSO_MOVE
end subroutine f90wrap_constants__get__MAX_PSO_MOVE

subroutine f90wrap_constants__get__BAR(f90wrap_BAR)
    use constants, only: constants_BAR => BAR
    implicit none
    real(8), intent(out) :: f90wrap_BAR
    
    f90wrap_BAR = constants_BAR
end subroutine f90wrap_constants__get__BAR

subroutine f90wrap_constants__get__KBAR(f90wrap_KBAR)
    use constants, only: constants_KBAR => KBAR
    implicit none
    real(8), intent(out) :: f90wrap_KBAR
    
    f90wrap_KBAR = constants_KBAR
end subroutine f90wrap_constants__get__KBAR

subroutine f90wrap_constants__get__ATM(f90wrap_ATM)
    use constants, only: constants_ATM => ATM
    implicit none
    real(8), intent(out) :: f90wrap_ATM
    
    f90wrap_ATM = constants_ATM
end subroutine f90wrap_constants__get__ATM

subroutine f90wrap_constants__get__MAX_SORT_DIM(f90wrap_MAX_SORT_DIM)
    use constants, only: constants_MAX_SORT_DIM => MAX_SORT_DIM
    implicit none
    integer, intent(out) :: f90wrap_MAX_SORT_DIM
    
    f90wrap_MAX_SORT_DIM = constants_MAX_SORT_DIM
end subroutine f90wrap_constants__get__MAX_SORT_DIM

subroutine f90wrap_constants__get__CHAR_LEN(f90wrap_CHAR_LEN)
    use constants, only: constants_CHAR_LEN => CHAR_LEN
    implicit none
    integer, intent(out) :: f90wrap_CHAR_LEN
    
    f90wrap_CHAR_LEN = constants_CHAR_LEN
end subroutine f90wrap_constants__get__CHAR_LEN

subroutine f90wrap_constants__get__MIN_DIS_SCALE(f90wrap_MIN_DIS_SCALE)
    use constants, only: constants_MIN_DIS_SCALE => MIN_DIS_SCALE
    implicit none
    real(8), intent(out) :: f90wrap_MIN_DIS_SCALE
    
    f90wrap_MIN_DIS_SCALE = constants_MIN_DIS_SCALE
end subroutine f90wrap_constants__get__MIN_DIS_SCALE

subroutine f90wrap_constants__get__MAX_DIS_SCALE(f90wrap_MAX_DIS_SCALE)
    use constants, only: constants_MAX_DIS_SCALE => MAX_DIS_SCALE
    implicit none
    real(8), intent(out) :: f90wrap_MAX_DIS_SCALE
    
    f90wrap_MAX_DIS_SCALE = constants_MAX_DIS_SCALE
end subroutine f90wrap_constants__get__MAX_DIS_SCALE

subroutine f90wrap_constants__get__MAX_BOND_SCALE(f90wrap_MAX_BOND_SCALE)
    use constants, only: constants_MAX_BOND_SCALE => MAX_BOND_SCALE
    implicit none
    real(8), intent(out) :: f90wrap_MAX_BOND_SCALE
    
    f90wrap_MAX_BOND_SCALE = constants_MAX_BOND_SCALE
end subroutine f90wrap_constants__get__MAX_BOND_SCALE

subroutine f90wrap_constants__get__PYTHON(f90wrap_PYTHON)
    use constants, only: constants_PYTHON => PYTHON
    implicit none
    character(24), intent(out) :: f90wrap_PYTHON
    
    f90wrap_PYTHON = constants_PYTHON
end subroutine f90wrap_constants__get__PYTHON

subroutine f90wrap_constants__get__workdir(f90wrap_workdir)
    use constants, only: constants_workdir => workdir
    implicit none
    character(256), intent(out) :: f90wrap_workdir
    
    f90wrap_workdir = constants_workdir
end subroutine f90wrap_constants__get__workdir

subroutine f90wrap_constants__set__workdir(f90wrap_workdir)
    use constants, only: constants_workdir => workdir
    implicit none
    character(256), intent(in) :: f90wrap_workdir
    
    constants_workdir = f90wrap_workdir
end subroutine f90wrap_constants__set__workdir

subroutine f90wrap_constants__get__LOGID(f90wrap_LOGID)
    use constants, only: constants_LOGID => LOGID
    implicit none
    integer, intent(out) :: f90wrap_LOGID
    
    f90wrap_LOGID = constants_LOGID
end subroutine f90wrap_constants__get__LOGID

subroutine f90wrap_constants__get__PARSEVER(f90wrap_PARSEVER)
    use constants, only: constants_PARSEVER => PARSEVER
    implicit none
    character(24), intent(out) :: f90wrap_PARSEVER
    
    f90wrap_PARSEVER = constants_PARSEVER
end subroutine f90wrap_constants__get__PARSEVER

subroutine f90wrap_constants__get__MAX_CHAR_LEN(f90wrap_MAX_CHAR_LEN)
    use constants, only: constants_MAX_CHAR_LEN => MAX_CHAR_LEN
    implicit none
    integer, intent(out) :: f90wrap_MAX_CHAR_LEN
    
    f90wrap_MAX_CHAR_LEN = constants_MAX_CHAR_LEN
end subroutine f90wrap_constants__get__MAX_CHAR_LEN

subroutine f90wrap_constants__get__LINE_WIDTH(f90wrap_LINE_WIDTH)
    use constants, only: constants_LINE_WIDTH => LINE_WIDTH
    implicit none
    integer, intent(out) :: f90wrap_LINE_WIDTH
    
    f90wrap_LINE_WIDTH = constants_LINE_WIDTH
end subroutine f90wrap_constants__get__LINE_WIDTH

subroutine f90wrap_constants__get__NOTOK_VALUE(f90wrap_NOTOK_VALUE)
    use constants, only: constants_NOTOK_VALUE => NOTOK_VALUE
    implicit none
    real(8), intent(out) :: f90wrap_NOTOK_VALUE
    
    f90wrap_NOTOK_VALUE = constants_NOTOK_VALUE
end subroutine f90wrap_constants__get__NOTOK_VALUE

subroutine f90wrap_constants__array__lap_matrix(dummy_this, nd, dtype, dshape, dloc)
    use kinds
    use constants, only: constants_lap_matrix => lap_matrix
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(constants_lap_matrix)) then
        dshape(1:2) = shape(constants_lap_matrix)
        dloc = loc(constants_lap_matrix)
    else
        dloc = 0
    end if
end subroutine f90wrap_constants__array__lap_matrix

! End of module constants defined in file ./FPP/Constants.fpp

