MODULE FingerPrints

    USE KINDS, ONLY: DP, i4b
    USE InputModule, ONLY: T_input_paras

    implicit none

    real(DP), allocatable :: SF1(:, :, :, :)
    real(DP), allocatable :: SF4(:, :)
    integer(i4b) :: nsimilar = 0

CONTAINS

    subroutine AddStructF1(simf1, same_flag)
        USE KINDS, ONLY: DP, i4b
        USE InputModule, ONLY: T_input_paras

        real(DP) :: ss(11, T_input_paras%number_of_species, T_input_paras%number_of_species)
        real(DP) :: simf1(11, T_input_paras%number_of_species, T_input_paras%number_of_species)
        real(DP) :: BondCharMatrix_CalDis

        integer(i4b)        :: nele, i
        logical             :: same_flag
        real(DP)            :: dist, tol
        tol = 0.01d0
        nele = T_input_paras%number_of_species

        same_flag = .false.
        if (nsimilar == 0) then
            nsimilar = nsimilar + 1
            SF1(:, :, :, nsimilar) = SimF1
        else
            do i = 1, nsimilar
                ss = SF1(:, :, :, i)
                dist = BondCharMatrix_CalDis(nele, ss, simf1)
                if (dist < Tol) then
                    same_flag = .true.
                    exit
                end if
            end do
            if (.not. same_flag) then
                nsimilar = nsimilar + 1
                SF1(:, :, :, nsimilar) = SimF1
            end if
        end if

    end subroutine

    subroutine AddStructF4(simf4, same_flag)

        USE KINDS, ONLY: DP, i4b
        USE InputModule, ONLY: T_input_paras
        USE BondCRT, ONLY: sn, calstrdis

        real(DP) :: simf4(T_input_paras%number_of_species*(T_input_paras%number_of_species + 1)/2.0*sn)
        real(DP) :: tmpf4(T_input_paras%number_of_species*(T_input_paras%number_of_species + 1)/2.0*sn)
        integer(i4b)        :: nele, i
        logical             :: same_flag
        real(DP)            :: dist, tol
        tol = 0.075d0

        nele = T_input_paras%number_of_species

        same_flag = .false.
        if (nsimilar == 0) then
            nsimilar = nsimilar + 1
            write (*, *) "tag1"
            SF4(:, nsimilar) = SimF4
            write (*, *) SimF4
            write (*, *) "tag2"
            write (*, *) "SF4(:, i)"
            write (*, *) SF4(1, 1)
        else
            write (*, *) "tag3"
            write (*, *) SF4(1, 1)
            do i = 1, nsimilar
                write (*, *) "tag4"
                write (*, *) "SF4(:, i)", i
                write (*, *) SF4(:, i)
                tmpf4 = SF4(:, i)
                write (*, *) "tag5"
                dist = calstrdis(tmpf4, simf4, 1, 1, nele)
                write (*, *) "tag6"
                if (dist < Tol) then
                    same_flag = .true.
                    exit
                end if
            end do
            if (.not. same_flag) then
                nsimilar = nsimilar + 1
                SF4(:, nsimilar) = SimF4
            end if
        end if

    end subroutine
END MODULE
