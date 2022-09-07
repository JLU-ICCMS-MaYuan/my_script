!   This module could calculate the distance between structures.
!   Author: Chuanxun Su
!   State Key Lab of Superhard Materials, Jilin University, Changchun.
!   Email: scx@calypso.cn
!   CREATED IN CHINA, Beijing
!   Date: 2014-7-10 21:21
!   Weather: Fine, 30 centigrade

!   Revision log:
!     10/07/2014  Program tested and worked. (Chuanxun Su)
!     09/04/2015  Add function storecfr and calstrdis. (Chuanxun Su)
!     15/05/2015  Optimize the algorithm. The efficiency is about 6 times
!                 higher than the old version, especially for big system with high symmetry.
!                 Improve the formula for character function around rcut(when ifr=2). (Chuanxun Su)
!     28/05/2015  A little improvement was made. The efficiency is about 15% higher than the
!                 last version, especially for elemental solids. Change the formula for
!                 character function a little. (Chuanxun Su)
!     02/06/2015  A small bug was fixed. (Chuanxun Su)

module BondCRT
    implicit none

    integer, parameter :: sn = 400
    real(8) :: rcut = 9.0d0, &
               dr = 0.006, dr2, &
               outr = 1.0d0, inr = 0.5d0, &
               pi = 3.1415926535897932384626433832795d0, delta, deltar, mtr(sn), tpfr(sn)
    real(8), allocatable :: minbd(:), wfscc(:, :), wfscc2(:, :), bondinf(:, :)
    integer, allocatable :: ibondinf(:, :)
    integer :: etypenub, atomnub, ifr = 2, ctnt
    logical :: lsetup = .false., l2d = .false., lct = .false.

contains

    subroutine storecfr(clycellv, atompsdr, atomid, cfr)

        ! Be very careful with clycellv, vector a is clycellv(1,:).
        real(8), intent(in) :: clycellv(3, 3), atompsdr(:, :)
        integer, intent(in) :: atomid(:)
        real(8), intent(out) :: cfr(:)
        real(8):: cellv(3, 3), lth
        integer :: i1, k, bk, i, binfn

        cellv = transpose(clycellv)
        if (.not. lsetup) then
            deltar = (rcut + outr)/(sn - 1)
            forall (i1=1:sn) mtr(i1) = deltar*(i1 - 1)
            lsetup = .true.
        end if

        call calcbond(cellv, atompsdr, atomid)

        binfn = size(bondinf, 1)
        bk = (etypenub + 1)*etypenub/2
        if (.not. allocated(minbd)) allocate (minbd(bk))
        minbd = 1000.0d0
        do k = 1, bk
            do i = 2, binfn
                if (bondinf(i, k) > 1.0d-10) then
                    minbd(k) = dr2*(i - 1)
                    exit
                end if
            end do
        end do

        cfr = 0.0d0
        do k = 1, bk
        do i = nint(minbd(k)/dr2 + 1), binfn
            if (bondinf(i, k) < 1.0d-10) cycle
            lth = dr2*(i - 1)
            if (ifr == 1) then
                call gasf(funcwr2(lth/minbd(k))*bondinf(i, k), lth, tpfr)
            else
                call gasf(funcwr3(lth)*bondinf(i, k), lth, tpfr)
            end if
            cfr((k - 1)*sn + 1:k*sn) = cfr((k - 1)*sn + 1:k*sn) + tpfr
        end do
        end do
        deallocate (bondinf, ibondinf, minbd)
        return
    end subroutine storecfr

    function calstrdis(cfr1, cfr2, spgit1, spgit2, ne)

        real(8) :: calstrdis, cfr1(:), cfr2(:)
        integer :: ne, i, j, loc, spgit1, spgit2

        if (spgit1 /= spgit2 .and. spgit1 /= 0 .and. spgit1 /= 999 .and. spgit2 /= 0 .and. spgit2 /= 999) then
            calstrdis = 1.0d0
            return
        end if
        ctnt = ctnt + 1
        if (.not. allocated(wfscc)) allocate (wfscc(ne, ne))
        wfscc = 0.0d0
        loc = 1
        do i = 1, ne
            do j = i, ne
                wfscc(j, i) = sum(cfr1(loc:loc + sn - 1))
                loc = loc + sn
            end do
        end do
        if (abs(sum(wfscc)) < 1.0d-20) then
            wfscc = 2.0d0/ne/(ne + 1)
        else
            wfscc = wfscc/sum(wfscc)
        end if

        if (.not. allocated(wfscc2)) allocate (wfscc2(ne, ne))
        wfscc2 = 0.0d0
        loc = 1
        do i = 1, ne
            do j = i, ne
                wfscc2(j, i) = sum(cfr2(loc:loc + sn - 1))
                loc = loc + sn
            end do
        end do
        if (abs(sum(wfscc2)) < 1.0d-20) then
            wfscc2 = 2.0d0/ne/(ne + 1)
        else
            wfscc2 = wfscc2/sum(wfscc2)
        end if
        wfscc = (wfscc + wfscc2)/2.0d0
        loc = 1
        calstrdis = 0.0d0
        do i = 1, ne
            do j = i, ne
                calstrdis = calstrdis + wfscc(j, i)*scc(cfr1(loc:loc + sn - 1), cfr2(loc:loc + sn - 1), sn)
                loc = loc + sn
            end do
        end do
        calstrdis = 1.0d0 - calstrdis
        deallocate (wfscc, wfscc2)
        return

    end function calstrdis

    function funcwr2(x)

        real(8) :: funcwr2, x, wfa = 2.20, wfb = 1.0

        funcwr2 = exp(-wfa*(x - wfb))
    end function funcwr2

    function funcwr3(x)

        real(8) :: funcwr3, x

        if (x < rcut - inr) then
            funcwr3 = 1.0d0
        else
            funcwr3 = exp(-3.0d0*(x - rcut + inr)/inr)
            return
        end if
        return
    end function

    subroutine gasf(wf, b, gr)

        real(8) :: wf, a = 60
        real(8) :: b
        real(8), dimension(:) :: gr
        integer :: i1

        forall (i1=1:sn) gr(i1) = wf*sqrt(a/pi)*exp(-a*(mtr(i1) - b)**2)
    end subroutine gasf

    subroutine calcbond(cellv, atompsdr, atomid)

        real(8) :: atompsdr(:, :), cellv(3, 3), recpv(3, 3), temp
        real(8), allocatable :: atomps(:, :)
        integer :: atomid(:), nabc(3), na, nb, nc, binfn, i, j, k, i1, i2, tpn
        integer, allocatable :: elenub(:)

        etypenub = size(atomid) - 1
        allocate (elenub(etypenub))
        elenub = atomid(2:etypenub + 1) - atomid(1:etypenub)
        atomnub = atomid(etypenub + 1) - 1
        binfn = rcut/dr + 2.0d0
        dr2 = rcut/(binfn - 1)
      if (.not. allocated(bondinf)) allocate (ibondinf(binfn, etypenub*(etypenub + 1)/2), bondinf(binfn, etypenub*(etypenub + 1)/2))
        ibondinf = 0
        recpv = recipvector(cellv)
        nabc(1) = ceiling(rcut*vectorlength(recpv(:, 1))/pi/2)
        nabc(2) = ceiling(rcut*vectorlength(recpv(:, 2))/pi/2)
        nabc(3) = ceiling(rcut*vectorlength(recpv(:, 3))/pi/2)
        !if(l2d) nabc(3)=0
        if (lct) nabc = 0
        allocate (atomps(3, atomnub))

        forall (i1=1:atomnub) atomps(:, i1) = cellv(:, 1)*atompsdr(1, i1) &
                                              + cellv(:, 2)*atompsdr(2, i1) + cellv(:, 3)*atompsdr(3, i1)
        k = 1
        do i = 1, etypenub
            do na = -nabc(1), nabc(1)
            do nb = -nabc(2), nabc(2)
            do nc = -nabc(3), nabc(3)
                if (na == 0 .and. nb == 0 .and. nc == 0) cycle
                temp = vectorlength(na*cellv(:, 1) + nb*cellv(:, 2) + nc*cellv(:, 3))
                if (temp > rcut) cycle
                tpn = nint(temp/dr2) + 1
                !print *, 'temp',temp,dr2,rcut,dr
                !print *, 'tpn1',tpn,k
                !if ( tpn<0 ) then
                !print *, 'cellv',cellv
                !print *, 'na',na,nb,nc
                !endif
                ibondinf(tpn, k) = ibondinf(tpn, k) + atomid(i + 1) - atomid(i)
            end do
            end do
            end do
            do na = -nabc(1), nabc(1)
            do nb = -nabc(2), nabc(2)
            do nc = -nabc(3), nabc(3)
                do i1 = atomid(i), atomid(i + 1) - 1
                    do i2 = i1 + 1, atomid(i + 1) - 1
                        temp = vectorlength(atomps(:, i1) - atomps(:, i2) - na*cellv(:, 1) - nb*cellv(:, 2) - nc*cellv(:, 3))
                        if (temp > rcut) cycle
                        tpn = nint(temp/dr2) + 1
                        !print *, 'tpn2',tpn,k
                        ibondinf(tpn, k) = ibondinf(tpn, k) + 2
                    end do
                end do
            end do
            end do
            end do
            k = k + 1
            do j = i + 1, etypenub
                do na = -nabc(1), nabc(1)
                do nb = -nabc(2), nabc(2)
                do nc = -nabc(3), nabc(3)
                    do i1 = atomid(i), atomid(i + 1) - 1
                    do i2 = atomid(j), atomid(j + 1) - 1
                        temp = vectorlength(atomps(:, i1) - atomps(:, i2) - na*cellv(:, 1) - nb*cellv(:, 2) - nc*cellv(:, 3))
                        if (temp > rcut) cycle
                        tpn = nint(temp/dr2) + 1
                        !print *, 'tpn3',tpn,k
                        ibondinf(tpn, k) = ibondinf(tpn, k) + 2
                    end do
                    end do
                end do
                end do
                end do
                k = k + 1
            end do
        end do
        bondinf = ibondinf/dfloat(2*atomnub)
        deallocate (elenub, atomps)

    end subroutine calcbond

    function vectorlength(vc)
        real(8) :: vc(3), vectorlength
        vectorlength = sqrt(vc(1)**2 + vc(2)**2 + vc(3)**2)
    end function

    function recipvector(lat)
        real(8), intent(in) :: lat(:, :)
        real(8) :: recipvector(3, 3)

        recipvector(:, 1) = crossp(lat(:, 2), lat(:, 3))
        recipvector(:, 2) = crossp(lat(:, 3), lat(:, 1))
        recipvector(:, 3) = crossp(lat(:, 1), lat(:, 2))
        recipvector = recipvector/volume(lat)*pi*2._8

    end function

    function volume(lat)
        real(8), intent(in) :: lat(:, :)
        real(8) :: volume

        volume = abs(sum(lat(:, 1)*crossp(lat(:, 2), lat(:, 3))))

    end function

    function crossp(va, vb)
        real(8), intent(in) :: va(3), vb(3)
        real(8) :: crossp(3)

        crossp(1) = va(2)*vb(3) - va(3)*vb(2)
        crossp(2) = va(3)*vb(1) - va(1)*vb(3)
        crossp(3) = va(1)*vb(2) - va(2)*vb(1)
    end function

    function scc(xi, yi, npt)

        !USE kinds,             ONLY : DP,i4b
        implicit none
        integer, parameter :: dp = 8, i4b = 4
        !!! pearson correlation coefficient

        real(DP)     :: xi(npt), yi(npt)
        integer(i4b) :: npt, i
        real(DP)     :: r, sum1, sum2, sum3, avx, avy
        real(DP)     :: scc
        avx = sum(xi)/npt
        avy = sum(yi)/npt
        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0
        do i = 1, npt
            sum1 = sum1 + (xi(i) - avx)*(yi(i) - avy)
            sum2 = sum2 + (xi(i) - avx)**2
            sum3 = sum3 + (yi(i) - avy)**2
        end do

        if (sum2 + sum3 == 0.0d0) then
            scc = 1.0d0
            return
        else if (sum2 == 0.0d0 .or. sum3 == 0.0d0) then
            scc = 0.0d0
            return
        end if
        scc = sum1/sqrt((sum2)*(sum3))

    end function scc

end module
