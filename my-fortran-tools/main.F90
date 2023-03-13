program myprogram
    ! implicit none
    use subF1 ! use语句必须放在implicit none的前面
    implicit none
    integer :: global_var
   
    global_var = 42
    call init_values(global_var)
    call get_spg_info()
    call gen_sg_num()
    write(*,*) "Global variable value is ", global_var
end program myprogram