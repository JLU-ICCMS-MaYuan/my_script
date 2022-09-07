! Module kinds defined in file ./FPP/Kinds.fpp

subroutine f90wrap_kinds__get__i4b(f90wrap_i4b)
    use kinds, only: kinds_i4b => i4b
    implicit none
    integer, intent(out) :: f90wrap_i4b
    
    f90wrap_i4b = kinds_i4b
end subroutine f90wrap_kinds__get__i4b

subroutine f90wrap_kinds__get__DP(f90wrap_DP)
    use kinds, only: kinds_DP => DP
    implicit none
    integer, intent(out) :: f90wrap_DP
    
    f90wrap_DP = kinds_DP
end subroutine f90wrap_kinds__get__DP

subroutine f90wrap_kinds__get__CHAR_LEN(f90wrap_CHAR_LEN)
    use kinds, only: kinds_CHAR_LEN => CHAR_LEN
    implicit none
    integer, intent(out) :: f90wrap_CHAR_LEN
    
    f90wrap_CHAR_LEN = kinds_CHAR_LEN
end subroutine f90wrap_kinds__get__CHAR_LEN

! End of module kinds defined in file ./FPP/Kinds.fpp

