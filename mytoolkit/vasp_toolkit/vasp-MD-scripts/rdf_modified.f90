!November 10, 2015  by Peihao Huang --<<Computer simulation of liquids>> page 183
program main
    implicit none
    integer :: k, i, j, Natoms, step, Nstep, n, n1, n2, BIN, MAXBIN, reject, Nreject
    real(kind=8) :: factor, a(3,3), r1, r2, r3, Rx, Ry, Rz, Rsq, R, DELR, RMAX, CONST, RHO, RLOWER, RUPPER, NIDEAL, RMIDDLE, Volume, Total, Coor, Rmin
    real(kind=8), allocatable :: x(:), y(:), z(:), GR(:)
    integer, allocatable :: HIST(:)
    real, parameter :: PI = 3.14159265
    character(len=256) :: input_file
    ! These parameters should be input by hand
    DELR = 0.1
    RMAX = 5.9
    Nreject = 6000
    Nstep = 10477

    ! Get input file name from command line arguments
    call get_command_argument(1, input_file)

    ! Initialize
    MAXBIN = INT(ANINT(RMAX / DELR))
    allocate(HIST(MAXBIN))
    allocate(GR(MAXBIN))
    do k = 1, MAXBIN
        HIST(k) = 0
    end do

    ! Read the structural information
    open(unit=10, file=trim(adjustl(input_file)))   ! Open the file specified by the command line argument
    read(10, *)
    read(10, *) factor
    read(10, *) a(1,1), a(1,2), a(1,3)
    read(10, *) a(2,1), a(2,2), a(2,3)
    read(10, *) a(3,1), a(3,2), a(3,3)
    do i = 1, 3
        do j = 1, 3
            a(i, j) = a(i, j) * factor
        end do
    end do
    read(10, *)
    read(10, *) Natoms
    allocate(x(Natoms))
    allocate(y(Natoms))
    allocate(z(Natoms))

    ! Calculate the average density of the system
    Volume = a(1,2) * a(2,3) * a(3,1) + a(1,3) * a(2,1) * a(3,2) + a(1,1) * a(2,2) * a(3,3) &
            - a(1,3) * a(2,2) * a(3,1) - a(1,1) * a(2,3) * a(3,2) - a(1,2) * a(2,1) * a(3,3)
    RHO = REAL(Natoms) / Volume
    write(*, *) "The total Volume is ", Volume

    ! Reject some steps before equilibrium
    do reject = 1, Nreject
        read(10, *)
        do n = 1, Natoms
            read(10, *)
        end do
    end do

    ! Read the equilibrium steps
    do step = 1, Nstep
        read(10, *)
        do n = 1, Natoms
            read(10, *) x(n), y(n), z(n)
        end do
        do n1 = 1, Natoms-1
            do n2 = n1 + 1, Natoms
                ! Calculate minimum image distances
                r1 = x(n2) - x(n1)
                r1 = r1 - ANINT(r1)
                r2 = y(n2) - y(n1)
                r2 = r2 - ANINT(r2)
                r3 = z(n2) - z(n1)
                r3 = r3 - ANINT(r3)
                ! Fractional to cartesian
                Rx = r1 * a(1,1) + r2 * a(2,1) + r3 * a(3,1)
                Ry = r1 * a(1,2) + r2 * a(2,2) + r3 * a(3,2)
                Rz = r1 * a(1,3) + r2 * a(2,3) + r3 * a(3,3)
                Rsq = Rx**2 + Ry**2 + Rz**2
                R = SQRT(Rsq)
                ! Sort the distances between the N atoms into histogram
                BIN = INT(R / DELR) + 1
                if (BIN <= MAXBIN) then
                    HIST(BIN) = HIST(BIN) + 2
                end if
            end do
        end do
    end do

    ! Normalize
    CONST = 4.0 * PI * RHO / 3
    do BIN = 1, MAXBIN
        RLOWER = REAL(BIN - 1) * DELR
        RUPPER = RLOWER + DELR
        NIDEAL = CONST * (RUPPER**3 - RLOWER**3)
        GR(BIN) = REAL(HIST(BIN)) / REAL(Nstep) / REAL(Natoms) / NIDEAL
    end do
    close(10)

    ! Write output to file
    open(unit=20, file='rdf.dat')
    do BIN = 1, MAXBIN
        RMIDDLE = REAL(BIN - 1) * DELR + DELR / 2
        write(20, *) RMIDDLE, GR(BIN)
    end do
    close(20)

    ! Check (RMAX should be set to sqrt(3)/2 * L for cubic box)
    Total = 0
    do BIN = 1, MAXBIN
        RLOWER = REAL(BIN - 1) * DELR
        RUPPER = RLOWER + DELR
        NIDEAL = CONST * (RUPPER**3 - RLOWER**3)
        Total = Total + GR(BIN) * NIDEAL
    end do
    write(*, *) "The total number of atoms is ", Total

    ! Calculate coordination number (method1: input by hand)
    Coor = 0
    do BIN = 1, 44
        RLOWER = REAL(BIN - 1) * DELR
        RUPPER = RLOWER + DELR
        NIDEAL = CONST * (RUPPER**3 - RLOWER**3)
        Coor = Coor + GR(BIN) * NIDEAL
        Rmin = REAL(BIN - 1) * DELR + DELR / 2
    end do
    write(*, *) "The coordination number is ", Coor, ", and rmin is ", Rmin

    ! (method 2)
    Coor = 0
    do BIN = 1, MAXBIN
        RLOWER = REAL(BIN - 1) * DELR
        RUPPER = RLOWER + DELR
        NIDEAL = CONST * (RUPPER**3 - RLOWER**3)
        Coor = Coor + GR(BIN) * NIDEAL
        Rmin = REAL(BIN - 1) * DELR + DELR / 2
        if ((BIN > 1) .and. (BIN < MAXBIN) .and. (GR(BIN) < GR(BIN - 1)) .and. (GR(BIN) < GR(BIN + 1))) then
            write(*, *) "The coordination number is ", Coor, ", and rmin is ", Rmin
            exit
        end if
    end do

    ! Deallocate arrays
    deallocate(x)
    deallocate(y)
    deallocate(z)
    deallocate(HIST)
    deallocate(GR)
end program main
