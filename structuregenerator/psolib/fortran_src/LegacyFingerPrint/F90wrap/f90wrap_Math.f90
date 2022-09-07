! Module linear defined in file ./FPP/Math.fpp

subroutine f90wrap_upper(matrix1, matrix2, n0, n1, n2, n3)
    use linear, only: upper
    implicit none
    
    real(8), dimension(n0,n1) :: matrix1
    real(8), dimension(n2,n3) :: matrix2
    integer :: n0
    !f2py intent(hide), depend(matrix1) :: n0 = shape(matrix1,0)
    integer :: n1
    !f2py intent(hide), depend(matrix1) :: n1 = shape(matrix1,1)
    integer :: n2
    !f2py intent(hide), depend(matrix2) :: n2 = shape(matrix2,0)
    integer :: n3
    !f2py intent(hide), depend(matrix2) :: n3 = shape(matrix2,1)
    call upper(matrix1=matrix1, matrix2=matrix2)
end subroutine f90wrap_upper

subroutine f90wrap_lower(matrix1, matrix2, n0, n1, n2, n3)
    use linear, only: lower
    implicit none
    
    real(8), dimension(n0,n1) :: matrix1
    real(8), dimension(n2,n3) :: matrix2
    integer :: n0
    !f2py intent(hide), depend(matrix1) :: n0 = shape(matrix1,0)
    integer :: n1
    !f2py intent(hide), depend(matrix1) :: n1 = shape(matrix1,1)
    integer :: n2
    !f2py intent(hide), depend(matrix2) :: n2 = shape(matrix2,0)
    integer :: n3
    !f2py intent(hide), depend(matrix2) :: n3 = shape(matrix2,1)
    call lower(matrix1=matrix1, matrix2=matrix2)
end subroutine f90wrap_lower

! End of module linear defined in file ./FPP/Math.fpp

