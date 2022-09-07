subroutine f90wrap_bondcharmatrix_caldis(nele, sim_func1, ret_bondcharmatrix_caldis, sim_func2, n0, n1, n2, n3)
    implicit none
    external bondcharmatrix_caldis
    real(8) bondcharmatrix_caldis
    
    integer(4), intent(in) :: nele
    real(8), intent(in), dimension(11,n0,n1) :: sim_func1
    real(8), intent(out) :: ret_bondcharmatrix_caldis
    real(8), intent(in), dimension(11,n2,n3) :: sim_func2
    integer :: n0
    !f2py intent(hide), depend(sim_func1) :: n0 = shape(sim_func1,1)
    integer :: n1
    !f2py intent(hide), depend(sim_func1) :: n1 = shape(sim_func1,2)
    integer :: n2
    !f2py intent(hide), depend(sim_func2) :: n2 = shape(sim_func2,1)
    integer :: n3
    !f2py intent(hide), depend(sim_func2) :: n3 = shape(sim_func2,2)
    ret_bondcharmatrix_caldis = bondcharmatrix_caldis(nele, sim_func1, sim_func2)
end subroutine f90wrap_bondcharmatrix_caldis

subroutine f90wrap_cal_volume(ret_cal_volume, matrix, n0)
    implicit none
    external cal_volume
    real(8) cal_volume
    
    real(8), intent(out) :: ret_cal_volume
    real(8), intent(in), dimension(3,n0) :: matrix
    integer :: n0
    !f2py intent(hide), depend(matrix) :: n0 = shape(matrix,1)
    ret_cal_volume = cal_volume(matrix)
end subroutine f90wrap_cal_volume

subroutine f90wrap_cal_area(ret_cal_area, matrix, n0)
    implicit none
    external cal_area
    real(8) cal_area
    
    real(8), intent(out) :: ret_cal_area
    real(8), intent(in), dimension(3,n0) :: matrix
    integer :: n0
    !f2py intent(hide), depend(matrix) :: n0 = shape(matrix,1)
    ret_cal_area = cal_area(matrix)
end subroutine f90wrap_cal_area

subroutine f90wrap_norm(ret_norm, a)
    implicit none
    external norm
    real(8) norm
    
    real(8), intent(out) :: ret_norm
    real(8), dimension(3) :: a
    ret_norm = norm(a)
end subroutine f90wrap_norm

subroutine f90wrap_detinmath(ret_detinmath, matrix, n0)
    implicit none
    external detinmath
    real(8) detinmath
    
    real(8), intent(out) :: ret_detinmath
    real(8), intent(in), dimension(3,n0) :: matrix
    integer :: n0
    !f2py intent(hide), depend(matrix) :: n0 = shape(matrix,1)
    ret_detinmath = detinmath(matrix)
end subroutine f90wrap_detinmath

subroutine f90wrap_crossprod(a, b, x)
    implicit none
    external crossprod
    
    real(8), dimension(3) :: a
    real(8), dimension(3) :: b
    real(8), dimension(3) :: x
    call crossprod(a, b, x)
end subroutine f90wrap_crossprod

subroutine f90wrap_lat_inv(matrix3, matrix2, n0, n1)
    implicit none
    external lat_inv
    
    real(8), dimension(3,n0) :: matrix3
    real(8), dimension(3,n1) :: matrix2
    integer :: n0
    !f2py intent(hide), depend(matrix3) :: n0 = shape(matrix3,1)
    integer :: n1
    !f2py intent(hide), depend(matrix2) :: n1 = shape(matrix2,1)
    call lat_inv(matrix3, matrix2)
end subroutine f90wrap_lat_inv

subroutine f90wrap_mul_vm(r, mat, p, n0)
    implicit none
    external mul_vm
    
    real(8), dimension(3) :: r
    real(8), dimension(3,n0) :: mat
    real(8), dimension(3) :: p
    integer :: n0
    !f2py intent(hide), depend(mat) :: n0 = shape(mat,1)
    call mul_vm(r, mat, p)
end subroutine f90wrap_mul_vm

subroutine f90wrap_quick_sort(a, n, s, e, n0)
    implicit none
    external quick_sort
    
    real(8), dimension(n0) :: a
    integer(4) :: n
    integer(4) :: s
    integer(4) :: e
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,0)
    call quick_sort(a, n, s, e)
end subroutine f90wrap_quick_sort

subroutine f90wrap_sortint_id(a, n, idx, n0, n1)
    implicit none
    external sortint_id
    
    integer(4), intent(in), dimension(n0) :: a
    integer(4), intent(in) :: n
    integer(4), intent(inout), dimension(n1) :: idx
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,0)
    integer :: n1
    !f2py intent(hide), depend(idx) :: n1 = shape(idx,0)
    call sortint_id(a, n, idx)
end subroutine f90wrap_sortint_id

subroutine f90wrap_int_to_char(ret_int_to_char, int_bn)
    implicit none
    external int_to_char
    character(6) int_to_char
    
    character(6), intent(out) :: ret_int_to_char
    integer, intent(in) :: int_bn
    ret_int_to_char = int_to_char(int_bn)
end subroutine f90wrap_int_to_char

subroutine f90wrap_eign(a, n, r, n0, n1, n2, n3)
    implicit none
    external eign
    
    real(8), intent(inout), dimension(n0,n1) :: a
    integer(4), intent(in) :: n
    real(8), intent(inout), dimension(n2,n3) :: r
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,0)
    integer :: n1
    !f2py intent(hide), depend(a) :: n1 = shape(a,1)
    integer :: n2
    !f2py intent(hide), depend(r) :: n2 = shape(r,0)
    integer :: n3
    !f2py intent(hide), depend(r) :: n3 = shape(r,1)
    call eign(a, n, r)
end subroutine f90wrap_eign

subroutine f90wrap_randomscale(lb, ret_randomscale, ub)
    implicit none
    external randomscale
    integer(4) randomscale
    
    integer(4) :: lb
    integer(4), intent(out) :: ret_randomscale
    integer(4) :: ub
    ret_randomscale = randomscale(lb, ub)
end subroutine f90wrap_randomscale

subroutine f90wrap_isprime(ret_isprime, num)
    implicit none
    external isprime
    logical isprime
    
    logical, intent(out) :: ret_isprime
    integer(4) :: num
    ret_isprime = isprime(num)
end subroutine f90wrap_isprime

subroutine f90wrap_matvectmul(mat, vect, mvec, n0)
    implicit none
    external matvectmul
    
    real(8), dimension(3,n0) :: mat
    real(8), dimension(3) :: vect
    real(8), dimension(3) :: mvec
    integer :: n0
    !f2py intent(hide), depend(mat) :: n0 = shape(mat,1)
    call matvectmul(mat, vect, mvec)
end subroutine f90wrap_matvectmul

subroutine f90wrap_rscale(lb, ret_rscale, ub)
    implicit none
    external rscale
    integer(4) rscale
    
    integer(4) :: lb
    integer(4), intent(out) :: ret_rscale
    integer(4) :: ub
    ret_rscale = rscale(lb, ub)
end subroutine f90wrap_rscale

subroutine f90wrap_dab(a, ret_dab, b)
    implicit none
    external dab
    real(8) dab
    
    real(8), dimension(3) :: a
    real(8), intent(out) :: ret_dab
    real(8), dimension(3) :: b
    ret_dab = dab(a, b)
end subroutine f90wrap_dab

subroutine f90wrap_det_i3(ret_det_i3, a, n0)
    implicit none
    external det_i3
    integer(4) det_i3
    
    integer(4), intent(out) :: ret_det_i3
    integer(4), dimension(3,n0) :: a
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,1)
    ret_det_i3 = det_i3(a)
end subroutine f90wrap_det_i3

subroutine f90wrap_dis_mat(a, b, dist, n0, n1)
    implicit none
    external dis_mat
    
    real(8), dimension(3,n0) :: a
    real(8), dimension(3,n1) :: b
    real(8) :: dist
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,1)
    integer :: n1
    !f2py intent(hide), depend(b) :: n1 = shape(b,1)
    call dis_mat(a, b, dist)
end subroutine f90wrap_dis_mat

subroutine f90wrap_inv_mat(a, inva, n0, n1)
    implicit none
    external inv_mat
    
    real(8), dimension(3,n0) :: a
    real(8), dimension(3,n1) :: inva
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,1)
    integer :: n1
    !f2py intent(hide), depend(inva) :: n1 = shape(inva,1)
    call inv_mat(a, inva)
end subroutine f90wrap_inv_mat

subroutine f90wrap_lat2matrix(lat, matrix, iflag, n0)
    implicit none
    external lat2matrix
    
    real(8), intent(inout), dimension(6) :: lat
    real(8), intent(inout), dimension(3,n0) :: matrix
    integer(4), intent(in) :: iflag
    integer :: n0
    !f2py intent(hide), depend(matrix) :: n0 = shape(matrix,1)
    call lat2matrix(lat, matrix, iflag)
end subroutine f90wrap_lat2matrix

subroutine f90wrap_check_lat(logic, lat1_matrix, n0)
    implicit none
    external check_lat
    
    logical :: logic
    real(8), dimension(3,n0) :: lat1_matrix
    integer :: n0
    !f2py intent(hide), depend(lat1_matrix) :: n0 = shape(lat1_matrix,1)
    call check_lat(logic, lat1_matrix)
end subroutine f90wrap_check_lat

subroutine f90wrap_cal_vect(a, lat, n0)
    implicit none
    external cal_vect
    
    real(8), dimension(3,n0) :: a
    real(8), dimension(6) :: lat
    integer :: n0
    !f2py intent(hide), depend(a) :: n0 = shape(a,1)
    call cal_vect(a, lat)
end subroutine f90wrap_cal_vect

