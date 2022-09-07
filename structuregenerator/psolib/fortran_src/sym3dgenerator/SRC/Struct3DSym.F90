subroutine Struct3DSym(ntype, type_num, natom, volume, dist_atom, spesg1, spesg2, sgid, ocmatrix, latmatrix, tau, l_flag)
    !v, l_flag)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> FiLEName: Struct3DSym.F90
    !>      Programmer            Email                   Web
    !>      ==================    ======================= ========================
    !>      Yanchao Wang          wyc@calypso.cn          www.calypso.cn
    !>      Date          Programmer        Description of changes
    !>      ==========    ==============    ======================================
    !>      2017.05.17    WYC               create
    !> Discription:
    !>      This module is designed to generate the crystal structure witm symmetry
    !> Parameters:
    !>
    !> References:
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    USE SpgDataInfo
    USE SpgModule
    USE SpgCryLat
    USE WyckoffModule

    implicit none
    !input paramters
    integer(i4b), intent(in)  :: ntype, type_num(ntype), natom
    real(DP), intent(in)  :: volume
    real(DP), intent(in)  :: dist_atom(ntype, ntype)
    INTEGER(i4b), intent(in) :: spesg1, spesg2 !!! edited By YuXin
    !output parameters
    integer(i4b), intent(out) :: sgid
    integer(i4b), intent(out) :: ocmatrix(27, ntype)
    real(DP), intent(out) :: latmatrix(3, 3)
    real(DP), intent(out) :: tau(3, natom)
    ! real(DP), intent(out) :: v(3, natom)
    logical, intent(out) :: l_flag

    integer(i4b)          :: u_type(ntype)

    real(DP) :: pos(3, 4*natom)
    integer :: i, j, k, iloop, unatom

    logical  :: wyc_flag, flag
    pos = 0.0
    call spg_set_list(ntype, type_num, spesg1, spesg2)
    wyc_flag = .false.
    iloop = 0
    l_flag = .true.
    do while (.not. wyc_flag)
        if (iloop > 1000) then
            l_flag = .false.
            return
        end if
        iloop = iloop + 1
        call gen_sg_num(sgid)
        call SpgGenLat(volume, sgid, latmatrix)
        call get_spg_info(sgid)
        call SetWyckoff(ntype, type_num, natom, latmatrix, dist_atom, unatom, ocmatrix, pos, wyc_flag)
        !write (*, *) "latmatrix", latmatrix
    end do
    ! write (*, *) "iloop", iloop
    !call OutputStruct(ntype,type_num*spgdata%R,4*natom,latmatrix,pos)
    call Unit2Prim(natom, spgdata%C, latmatrix, pos, tau)
    !print *, 'sym OK'
    ! v = 0.d0
    ! write (*, *) "l_flag", l_flag
end

subroutine Struct3DSymSgid(ntype, type_num, natom, volume, dist_atom, sgid, ocmatrix, latmatrix, tau, l_flag)
    !v, l_flag)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !given an spacegroup index ,generate an structure
    !!!!!!!!
    USE SpgDataInfo
    USE SpgModule
    USE SpgCryLat
    USE WyckoffModule

    implicit none
    !input paramters
    integer(i4b), intent(in)  :: ntype, type_num(ntype), natom
    real(DP), intent(in)  :: volume
    real(DP), intent(in)  :: dist_atom(ntype, ntype)
    INTEGER(i4b) :: spesg1, spesg2 !!! edited By YuXin
    !output parameters
    integer(i4b), intent(in) :: sgid
    integer(i4b), intent(out) :: ocmatrix(27, ntype)
    real(DP), intent(out) :: latmatrix(3, 3)
    real(DP), intent(out) :: tau(3, natom)
    ! real(DP), intent(out) :: v(3, natom)
    logical, intent(out) :: l_flag
    integer :: iloop
    real(DP) :: pos(3, 4*natom)
    integer :: unatom

    logical  :: wyc_flag
    spesg1 = 0
    spesg2 = 0
    pos = 0.0
    call spg_set_list(ntype, type_num, spesg1, spesg2)
    if (spglist(sgid) .eqv. .true.) Then
        wyc_flag = .false.
        iloop = 0
        l_flag = .true.
        do while (.not. wyc_flag)
            if (iloop > 1000) then
                l_flag = .false.
                return
            end if
            iloop = iloop + 1
            call SpgGenLat(volume, sgid, latmatrix)
            call get_spg_info(sgid)
            call SetWyckoff(ntype, type_num, natom, latmatrix, dist_atom, unatom, ocmatrix, pos, wyc_flag)
            call Unit2Prim(natom, spgdata%C, latmatrix, pos, tau)
            ! v = 0.d0
        end do
    else
        write (*, *) "Given Space Group index inapprociate"
        l_flag = .false.
    end if
end subroutine Struct3DSymSgid

subroutine GetLegalSpgIndex(ntype, type_num, LegalSpgIndex)
    USE SpgModule
    integer(i4b), intent(in)  :: ntype, type_num(ntype)
    integer(i4b), dimension(230), intent(out) :: LegalSpgIndex
    integer :: iloop, count
    LegalSpgIndex = 0
    call spg_set_list(ntype, type_num, 0, 0)
    count = 1
    do iloop = 1, 230
        if (spglist(iloop) .eqv. .true.) Then
            LegalSpgIndex(count) = iloop
            count = count + 1
        end if
    end do
end subroutine GetLegalSpgIndex
