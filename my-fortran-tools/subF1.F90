module Structure
    implicit none
    type CrystalParameters ! 将spgnum定义为模块变量，因此可以在模块中的所有子例程中访问它
        integer :: space_number
        character(len=:), allocatable :: elements(:)
        integer, allocatable :: numberofatoms(:)
        real, allocatable :: coords(3,:)
    end type
    
    type(CrystalParameters) :: CrystParas

    contains

        subroutine init_values(spgnum)
            implicit none
            integer,intent(in) :: spgnum
            spgdata%space_number = spgnum
        end subroutine init_values

        subroutine spg_set_list()
            write(*,*) "This is subroutine spg_set_list"
        end subroutine spg_set_list
        
        subroutine get_spg_info()
            write(*,*) "This is subroutine get_spg_info"
        end subroutine get_spg_info

        subroutine gen_sg_num()
            write(*,*) spgdata%space_number
            write(*,*) "This is subroutine gen_sg_num"
        end subroutine gen_sg_num

end module SpgModule


