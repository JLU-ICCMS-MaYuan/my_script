! Module bondcharmatrix defined in file ./FPP/BondCharMatrix.fpp

subroutine f90wrap_bondcharmatrix_calbcm(structyp, nele, natom, ele_num, lat_matrix, tau, sim_func, n0, n1, n2, n3, n4)
    use bondcharmatrix, only: bondcharmatrix_calbcm
    implicit none
    
    integer(4), intent(in) :: structyp
    integer(4), intent(in) :: nele
    integer(4), intent(in) :: natom
    integer(4), intent(in), dimension(n0) :: ele_num
    real(8), intent(in), dimension(3,n1) :: lat_matrix
    real(8), intent(in), dimension(3,n2) :: tau
    real(8), intent(inout), dimension(11,n3,n4) :: sim_func
    integer :: n0
    !f2py intent(hide), depend(ele_num) :: n0 = shape(ele_num,0)
    integer :: n1
    !f2py intent(hide), depend(lat_matrix) :: n1 = shape(lat_matrix,1)
    integer :: n2
    !f2py intent(hide), depend(tau) :: n2 = shape(tau,1)
    integer :: n3
    !f2py intent(hide), depend(sim_func) :: n3 = shape(sim_func,1)
    integer :: n4
    !f2py intent(hide), depend(sim_func) :: n4 = shape(sim_func,2)
    call bondcharmatrix_calbcm(structyp=structyp, nele=nele, natom=natom, ele_num=ele_num, lat_matrix=lat_matrix, tau=tau, &
        sim_func=sim_func)
end subroutine f90wrap_bondcharmatrix_calbcm

subroutine f90wrap_bondcharmatrix_car2sphe(x, y, z, r, theta, psi)
    use bondcharmatrix, only: bondcharmatrix_car2sphe
    implicit none
    
    real(8) :: x
    real(8) :: y
    real(8) :: z
    real(8) :: r
    real(8) :: theta
    real(8) :: psi
    call bondcharmatrix_car2sphe(x=x, y=y, z=z, r=r, theta=theta, psi=psi)
end subroutine f90wrap_bondcharmatrix_car2sphe

subroutine f90wrap_bondcharmatrix_sphericalharmonic(l, m, angle1, angle2, y)
    use bondcharmatrix, only: bondcharmatrix_sphericalharmonic
    implicit none
    
    integer(4) :: l
    integer(4) :: m
    real(8) :: angle1
    real(8) :: angle2
    complex :: y
    call bondcharmatrix_sphericalharmonic(l=l, m=m, angle1=angle1, angle2=angle2, Y=y)
end subroutine f90wrap_bondcharmatrix_sphericalharmonic

! End of module bondcharmatrix defined in file ./FPP/BondCharMatrix.fpp

