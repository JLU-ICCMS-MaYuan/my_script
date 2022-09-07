module WyckoffModule
    USE kinds, ONLY: DP, i4b
    USE SpgDataInfo
    USE SpgModule
    implicit none
    type(spgstruct), allocatable, dimension(:) :: Myspgs
    !Myspgs is used for substitute origin derived tyor array in setWyckoff
contains
    subroutine createMyspgs(length)
        integer, intent(in)::length
        if (allocated(Myspgs)) deallocate (Myspgs)
        if (.not. allocated(Myspgs)) allocate (Myspgs(1:length))
    end subroutine createMyspgs
    subroutine SetWyckoff(ntyp, typ, natom, latmatrix, dist_atom, u_natom, ocmatrix, pos, wycflag)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !> Filename:SetWyckoff.F90
        !>      Programmer            Email                   Web
        !>      ==================    ======================= ========================
        !>      Yanchao Wang          wyc@calypso.cn          www.calypso.cn
        !> Record of revisions:
        !>      Date          Programmer        Description of changes
        !>      ==========    ==============    ======================================
        !>      2017.05.10    Yanchao Wang      First build the module
        !> Discription:
        !>      This is designed to generate atomic positions
        !> Parameters:
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        USE kinds, ONLY: DP, i4b
        USE SpgDataInfo
        USE SpgModule
        !
        implicit none
        !
        !input and output
        !
        integer(i4b), intent(in)        :: ntyp
        integer(i4b), intent(in)        :: typ(ntyp)
        integer(i4b), intent(in)        :: natom
        real(DP), intent(in)        :: latmatrix(3, 3), dist_atom(ntyp, ntyp)
        integer(i4b), intent(out)       :: u_natom
        integer(i4b), intent(out)       :: ocmatrix(27, ntyp)
        real(DP), intent(inout)     :: pos(3, 4*natom)
        logical, intent(out)       :: wycflag

        integer(i4b)                     :: randomscale
        !local variables
        integer(i4b)  :: maxtrial
        integer(i4b)  :: i, j, k, jj, m, n, tmp, bp, ndiff
        integer(i4b)  :: iloop, ii
        integer(i4b), external    :: rscale
        integer(i4b)  :: iatom
        real(DP)      :: resf(3, 192)
        logical       :: same_ocflag, findsign, l_flag, flag
        integer(i4b)  :: remain
        integer(i4b)  :: ele_num(ntyp)

        ! type(spgstruct), dimension(4*natom) :: spgs

        call createMyspgs(4*natom)

        u_natom = natom*SpgData%R
        maxtrial = 2000*SpgData%R
        !!!!!!!!!!!!!!!!!!!!!!! choice is wyckoff sequence for each element!!!!!
        flag = .false.
        iloop = 0
        do while (.not. flag)
            iloop = iloop + 1
            call GenOcMatrix(ntyp, typ, Ocmatrix, flag)
            if (iloop > 500) then
                wycflag = .false.
                return
            end if
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        iatom = 0
        ele_num = typ*SpgData%R

        do i = 1, ntyp !!! for each elements
            do tmp = 1, spgdata%NumMulti !!!!! for all wyckoff positions
                do k = 1, ocmatrix(tmp, i) !!! how many times occupy for wyckoff
                    iloop = 0
                    do while (.true.)
                        ! write(*,*) "Tag1"
                        call GetCoords(tmp, ndiff, resf)
                        !write(*,*) "Tag2"
                        call dist_check_cry(latmatrix, resf, ndiff, resf, ndiff, dist_atom(i, i), findsign)
                        !write(*,*) "Tag3"
                        if ((findsign .eqv. .false.) .or. (ndiff /= spgdata%multi(tmp))) then
                            if (iloop > maxtrial) then
                                wycflag = .false.
                                return
                            end if
                            iloop = iloop + 1
                            cycle
                        else
                            if (iatom == 0) then
                                do ii = 1, ndiff
                                    iatom = iatom + 1
                                    ! spgs(iatom)%itype = i
                                    ! spgs(iatom)%pos(:) = resf(:, ii)
                                    !write(*,*) "Tag4"
                                    Myspgs(iatom)%itype = i
                                    !write(*,*) "Tag5"
                                    !write(*,*) "resf(:, ii)",resf(:, ii)
                                    !write(*,*) "Myspgs(iatom)%pos(:)",Myspgs(iatom)%pos(:)
                                    Myspgs(iatom)%pos(:) = resf(:, ii)
                                    !write(*,*) "Tag6"
                                    pos(:, iatom) = resf(:, ii)
                                    !write(*,*) "Tag7"
                                end do
                                exit
                            else
                                !write(*,*) "Tag8"
                                call dist_check_tau(ntyp, natom, dist_atom, latmatrix, iatom, &
                                                    i, resf, ndiff, l_flag)
                                if (l_flag) then
                                    do ii = 1, ndiff
                                        iatom = iatom + 1
                                        if (iatom > 4*natom) then
                                            wycflag = .false.
                                            !pause
                                            return
                                        end if
                                        Myspgs(iatom)%itype = i
                                        Myspgs(iatom)%pos(:) = resf(:, ii)
                                        ! spgs(iatom)%itype = i
                                        ! spgs(iatom)%pos(:) = resf(:, ii)
                                        pos(:, iatom) = resf(:, ii)
                                    end do
                                    exit
                                else
                                    iloop = iloop + 1
                                end if
                            end if
                        end if
                        if (iloop > maxtrial) then
                            wycflag = .false.
                            return
                        end if
                    end do
                end do
            end do
        end do
        if (iatom == u_natom) then
            !        u_natom=iatom
            !        print *, "dd=",u_natom
            wycflag = .true.
        end if
    end subroutine
    subroutine dist_check_tau(nele, natom, dist_atom, latmat, nstr, itype, tau, ndif, l_flag)
        USE kinds, ONLY: DP, i4b
        USE SpgModule, ONLY: spgstruct

        implicit none

        integer(i4b), intent(in)  :: nele, natom, itype, nstr, ndif
        real(DP), intent(in)  :: dist_atom(nele, nele), tau(3, 192)
        real(DP), intent(in)  :: latmat(3, 3)
        logical, intent(out) :: l_flag

        integer(i4b)    :: i, ii, iix, iiy, iiz, ix, iy, iz
        real(DP)        :: frac_diff(3), cart_diff(3), distance
        !pause
        l_flag = .true.
        do i = 1, nstr
            do ix = -1, 1
                do iy = -1, 1
                    do iz = -1, 1
                        do ii = 1, ndif
                            frac_diff(1) = Myspgs(i)%pos(1) - tau(1, ii) + real(ix)
                            frac_diff(2) = Myspgs(i)%pos(2) - tau(2, ii) + real(iy)
                            frac_diff(3) = Myspgs(i)%pos(3) - tau(3, ii) + real(iz)
                            call mul_vm(frac_diff, latmat, cart_diff)
                            distance = cart_diff(1)**2 + cart_diff(2)**2 + cart_diff(3)**2
                            distance = sqrt(distance)
                            if (distance < dist_atom(Myspgs(i)%itype, itype)) then
                                l_flag = .false.
                                return
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end subroutine
    subroutine dist_check_cry(latmatrix, pos1, natom1, pos2, natom2, length, dflag)
        !
        USE kinds, ONLY: DP, i4b
        !
        implicit none
        !
        !input and output
        !
        real(DP), intent(in) :: pos1(3, 192), pos2(3, 192)
        real(DP), intent(in) :: latmatrix(3, 3)
        real(DP), intent(in) :: length
        integer(i4b), intent(in) :: natom1, natom2
        logical, intent(out):: dflag
        !

        !local variables
        !
        integer(i4b) :: i, j, k, jat, l, i1, i2, i3, n
        integer(i4b) :: s, iix, iiy, iiz !(extend to supercell)
        real(DP)     :: dist, dif(3), xx(3), pp(3), p(3), dif_cal(3)
        dflag = .true.
        iix = 1
        iiy = 1
        iiz = 1

        do i = 1, natom1
            do j = 1, 3
                xx(j) = pos1(j, i)
            end do
            do i1 = -iix, iix
                do i2 = -iiy, iiy
                    do i3 = -iiz, iiz
                        p(1) = i1
                        p(2) = i2
                        p(3) = i3
                        do jat = 1, natom2 !(the i-th atom in cell )
                            dist = 0.0
                            if (i == jat) then
                                dist = 100
                            else
                                do l = 1, 3
                                    pp(l) = pos2(l, jat) + p(l)
                                    dif(l) = xx(l) - pp(l)
                                end do
                                call mul_vm(dif, latmatrix(:, :), dif_cal)
                                dist = dif_cal(1)*dif_cal(1) + dif_cal(2)*dif_cal(2) + &
                                       dif_cal(3)*dif_cal(3)
                                dist = sqrt(dist)
                            end if
                            if (dist < length) then
                                dflag = .false.
                                return
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end subroutine dist_check_cry
    subroutine reduce_symmetry_pos(tmpb, nbs, ndiff, resf)
        !
        ! author : Yanchao Wang
        ! email  : wyc@calypso.cn
        ! date   : 8 May 2017
        ! address: HongKong University
        !
        USE kinds, ONLY: DP, i4b

        implicit none

        real(DP), intent(in) :: tmpb(3, 192)
        integer, intent(in) :: nbs
        integer, intent(out):: ndiff
        real(DP), intent(out):: resf(3, 192)

        ! local parameters
        real(DP) :: ci(3), cm(3)
        real(DP) :: dis_c(3), dist
        integer(i4b) :: i, j, k, l, m, n, p, q, s
        real(DP) :: tmp(3, nbs), c(3, nbs), b(3, nbs)
        logical :: dstat

        q = 0
        tmp = 0.0
        c = 0.0
        b = 0.0
        k = 0
        do i = 1, nbs
            do j = 1, 3
                tmp(j, i) = tmpb(j, i)
            end do
        end do

        l = nbs
        do i = 1, l

            do p = 1, l
                do j = 1, 3
                    c(j, p) = tmp(j, p)
                end do
            end do

            do s = 1, i
                do j = 1, 3
                    b(j, s) = c(j, s)
                end do
            end do

            k = i

            do m = i + 1, l
                do n = 1, 3
                    ci(n) = c(n, i)
                    cm(n) = c(n, m)
                    ci(n) = ci(n) - floor(ci(n))
                    cm(n) = cm(n) - floor(cm(n))
                    if (ci(n) > 0.9999) ci(n) = 0.0
                    if (cm(n) > 0.9999) cm(n) = 0.0
                end do
                dist = 0.0
                do n = 1, 3
                    dist = dist + abs(ci(n) - cm(n))
                end do
                if (dist < 0.0005) then
                    dstat = .true.
                else
                    dstat = .false.
                end if
                if (dstat) then
                    q = q + 1
                else if (.not. dstat) then
                    k = k + 1
                    do j = 1, 3
                        b(j, k) = c(j, m)
                    end do
                end if
            end do

            l = k

            do p = 1, k
                do j = 1, 3
                    tmp(j, p) = b(j, p)
                end do
            end do
        end do

        resf = 0.0

        ndiff = k - q

        do i = 1, ndiff
            do j = 1, 3
                resf(j, i) = tmp(j, i)
            end do
        end do
    end subroutine
    subroutine GetCoords(iwck, ndiff, resf)
        USE kinds, ONLY: DP, i4b
        USE SpgDataInfo
        USE SpgModule
        !
        implicit none
        !
        !   input and output
        !

        integer(i4b), intent(in)  :: iwck
        integer(i4b), intent(out) :: ndiff

        real(DP)      :: tau(3), tol, restmp, trans
        real(DP)      :: resp(3, 192), resf(3, 192), res_tran_tmp(3, 192)
        logical       :: same_ocflag, findsign, l_flag
        integer(i4b)  :: tmp, m, n, j, k, bp
        tol = 0.01
        tmp = iwck

        call random_number(tau(1))
        call random_number(tau(2))
        call random_number(tau(3))
        do j = 1, 3
            if (abs(spgData%Coor(j, tmp)) < 1.) then
                tau(j) = spgdata%Coor(j, tmp)
            else
                if (abs(spgdata%Coor(j, tmp) - 10.d0) < tol) then
                    tau(j) = tau(1)
                else if (abs(spgdata%Coor(j, tmp) - 10.d0) < 1.0) then
                    trans = spgdata%Coor(j, tmp) - 10.d0
                    tau(j) = tau(1) + trans
                else if (abs(spgdata%Coor(j, tmp) + 10.d0) < tol) then
                    tau(j) = 1 - tau(1)
                else if (abs(spgdata%Coor(j, tmp) + 10.d0) < 1.0) then
                    trans = spgdata%Coor(j, tmp) + 10.d0
                    tau(j) = -tau(1) + trans
                else if (abs(spgdata%Coor(j, tmp) - 15.d0) < tol) then
                    tau(j) = 2*tau(1)
                end if
                if (j == 3) then
                    if (abs(spgdata%Coor(j, tmp) - 20.d0) < tol) then
                        tau(j) = tau(2)
                    else if (abs(spgdata%Coor(j, tmp) - 20.d0) < 1.0) then
                        trans = spgdata%Coor(j, tmp) - 20.d0
                        tau(j) = tau(2) + trans
                    else if (abs(spgdata%Coor(j, tmp) + 20.d0) < tol) then
                        tau(j) = 1 - tau(2)
                    else if (abs(spgdata%Coor(j, tmp) + 20.d0) < 1.0) then
                        trans = spgdata%Coor(j, tmp) + 20.d0
                        tau(j) = -tau(2) + trans
                    end if
                end if
            end if
            tau(j) = tau(j) - floor(tau(j))
        end do
        bp = spgdata%Multi(1)
        do m = 1, bp
            do n = 1, 3
                restmp = 0.
                do j = 1, 3
                    restmp = restmp + spgdata%pgop(j, n, m)*tau(j)
                end do
                restmp = restmp - floor(restmp)
                resp(n, m) = restmp
            end do
        end do
        do m = 1, bp
            res_tran_tmp(:, m) = resp(:, m) + spgdata%tgop(:, m)
        end do
        do m = 1, bp
            do n = 1, 3
                res_tran_tmp(n, m) = res_tran_tmp(n, m) - floor(res_tran_tmp(n, m))
            end do
        end do
        call reduce_symmetry_pos(res_tran_tmp, bp, ndiff, resf)
    end subroutine GetCoords
    subroutine GenOcMatrix(nele, typ, Omat, flag)

        USE kinds, ONLY: DP, i4b
        USE SpgModule

        implicit none

        integer(i4b), intent(in)  :: nele, typ(nele)
        integer(i4b), intent(out) :: Omat(27, nele)
        logical, intent(out) :: flag

        integer(i4b), external  :: rscale, randomscale

        integer(i4b) :: i1, i, j, k, iloop
        integer(i4b) :: ele_num(nele)
        logical      :: cflag, piece, wck_flag(27)
        logical      :: flag_tp
        integer(i4b) :: OCm(27), occheck(27)
        integer(i4b) :: num, mu, remain, wck_bound, n_var
        integer(i4b) :: tmp
        Omat = 0
        num = SpgData%NumMulti
        mu = SpgData%Multi(num)
        ele_num = typ*SpgData%R
        OCm = 0

        do i = 1, nele
            remain = ele_num(i)
            do while (remain > 0)
                if (remain < SpgData%Multi(SpgData%NumMulti)) then
                    flag = .false.
                    return
                end if
                flag_tp = .false.
                do while (.not. flag_tp)
                    tmp = randomscale(1, SpgData%NumMulti)
                    if (SpgData%multi(tmp) <= remain) then
                        flag_tp = .true.
                    end if

                end do
                remain = remain - SpgData%multi(tmp)
                Omat(tmp, i) = Omat(tmp, i) + 1
            end do
        end do

        occheck = 0
        do i = 1, nele
            do j = 1, 27
                occheck(j) = occheck(j) + Omat(j, i)
            end do
        end do
        cflag = .true.
        do i = 1, SpgData%NumMulti
            if (spgdata%ocflag(i) == 0) then
                if (occheck(i) > 1) then
                    cflag = .false.
                    exit
                end if
            end if
        end do
        if (cflag) then
            flag = .true.
        else
            flag = .false.
        end if
    end subroutine GenOcMatrix

end module WyckoffModule
