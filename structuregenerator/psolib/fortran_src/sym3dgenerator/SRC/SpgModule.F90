MODULE SpgModule
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Filename:SpgModule.F90
    !>      Programmer            Email                   Web
    !>      ==================    ======================= ========================
    !>      Yanchao Wang          wyc@calypso.cn          www.calypso.cn
    !> Record of revisions:
    !>      Date          Programmer        Description of changes
    !>      ==========    ==============    ======================================
    !>      2017.05.10    Yanchao Wang      First build the module
    !> Discription:
    !>      This module is designed to generate structure by symmetry constraint.
    !> Parameters:
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    USE KINDS, ONLY: DP, i4b
    USE SpgDataInfo

    implicit none

    type spg_struct
        integer(i4b) :: NumMulti
        integer(i4b) :: R
        character    :: C
        integer(i4b), allocatable :: Multi(:)
        integer(i4b), allocatable :: Ocflag(:)
        real(DP), allocatable :: oc_probability(:)
        real(DP), allocatable :: pgop(:, :, :)
        real(DP), allocatable :: tgop(:, :)
        real(DP), allocatable :: Coor(:, :)
    end type

!type :: wyckoff
!   integer :: itype
!   integer :: wlist(27)
!   type(wyckoff),pointer :: next
!end type

    type spgstruct
        integer(i4b) :: itype
        real(DP)     :: pos(3)
    end type

!integer, allocatable :: nlist(:)
!integer, allocatable :: wyc_oc_matrix(:,:,:)
    logical              :: spglist(230)

    type(spg_struct)     :: spgdata

!END Module
contains

    subroutine spg_set_list(ntyp, typ, spesg1, spesg2)
        implicit none

        integer, intent(in) :: ntyp, typ(ntyp)
        INTEGER(i4b), intent(in) :: spesg1, spesg2

!local parameters
        integer :: minatom, i, j, m_mod
        integer :: elen(ntyp)

        spglist(:) = .true.
        if (spesg1 /= 0 .or. spesg2 /= 0) then
            spglist(:) = .false.
            do i = spesg1, spesg2
                spglist(i) = .true.
            end do
        end if

        do i = 1, 230
            elen(:) = typ(:)*spacdata(i)%R
            minatom = minval(elen)
            !print *, "i=",i,spacdata(i)%nummulti
            if (minatom < spacdata(i)%multi(spacdata(i)%nummulti)) then
                spglist(i) = .false.
            end if
            do j = 1, ntyp
                m_mod = mod(typ(j), spacdata(i)%multi(spacdata(i)%nummulti))
                if (m_mod /= 0) then
                    spglist(i) = .false.
                    exit
                end if
            end do
        end do

    end subroutine

    subroutine get_spg_info(sgindex)

        !USE SpgModule
        !USE SpgDataInfo

        implicit none

        integer, intent(in) :: sgindex

        ! local parameters
        integer :: numwycset, numop
        !print *, "sgindex=",sgindex
        numwycset = spacdata(sgindex)%nummulti
        numop = spacdata(sgindex)%multi(1)
        if (allocated(spgdata%Multi)) deallocate (spgdata%Multi)
        if (allocated(spgdata%ocflag)) deallocate (spgdata%ocflag)
        if (allocated(spgdata%pgop)) deallocate (spgdata%pgop)
        if (allocated(spgdata%tgop)) deallocate (spgdata%tgop)
        if (allocated(spgdata%coor)) deallocate (spgdata%coor)
        if (allocated(spgdata%oc_probability)) deallocate (spgdata%oc_probability)
        if (.not. allocated(spgdata%Multi)) allocate (spgdata%Multi(numwycset))
        if (.not. allocated(spgdata%ocflag)) allocate (spgdata%ocflag(numwycset))
        if (.not. allocated(spgdata%oc_probability)) allocate (spgdata%oc_probability(numwycset))
        if (.not. allocated(spgdata%pgop)) allocate (spgdata%pgop(3, 3, numop))
        if (.not. allocated(spgdata%tgop)) allocate (spgdata%tgop(3, numop))
        if (.not. allocated(spgdata%coor)) allocate (spgdata%coor(3, numwycset))

        spgdata%NumMulti = SpacData(sgindex)%NumMulti
        spgdata%Multi = SpacData(sgindex)%Multi(1:spgdata%NumMulti)
        spgdata%ocflag = SpacData(sgindex)%ocflag(1:spgdata%NumMulti)
        spgdata%pgop = SpacData(sgindex)%symematrix(1:3, 1:3, 1:numop)
        spgdata%tgop = SpacData(sgindex)%tranmatrix(1:3, 1:numop)
        spgdata%coor = SpacData(sgindex)%coor(1:3, 1:spgdata%NumMulti)
        spgdata%C = SpacData(sgindex)%C
        spgdata%R = SpacData(sgindex)%R
        spgdata%oc_probability = 1.0
    end subroutine

!subroutine set_wyck_list(ntype,type_num,spgdata,nlist,wyc_oc_matrix,l_flag)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!          modified file        !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!     author: Yanchao Wang       !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!     date  : 2015-4-21          !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!     email : wyc@calypso.cn  !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!USE kinds,           ONLY : DP,i4b
!
!implicit none
!
!input and output
!
!integer(i4b), intent (in)        :: ntype,type_num(ntype)
!logical, intent(out)             :: l_flag
!integer(i4b), intent (out)       :: nlist(:)
!integer(i4b), intent (out)       :: wyc_oc_matrix(:,:,:)
!

!type(spg_struct) :: spgdata

!integer(i4b), external  :: rscale

!local variables
!
!integer(i4b)     :: i,j,k,diff,nat,maxcount,icount
!
!integer(i4b)              :: dnatom,inatom,icol,maxnat
!integer(i4b)              :: row,col,numpossible
!integer(i4b)              :: ocmatrix(27,ntype)
!integer(i4b)              :: ocmax(27),ocvalue(27)
!integer(i4b),target       :: ocm(27)
!integer(i4b),pointer      :: ii
!logical  :: sameflag
!maxcount=1000
!ocmax=0
!ocmatrix=0
!maxnat=maxval(type_num)
!maxnat=maxnat*spgdata%R

!l_flag= .true.
!do i=1,spgdata%NumMulti
!        print *, maxnat,spgdata%Multi(i)
!        if (spgdata%ocflag(i)==0) then
!                ocmax(i)=1
!                ocm(i)=2
!        else
!                ocmax(i)=floor(real(maxnat)/spgdata%multi(i))
!                ocm(i)=ocmax(i)+1
!        endif
!        if(ocm(i)==0) ocm(i)=1
!        print *, ocm(i)
!        pause
!enddo
!call GenOcMatrix(spgdata,ntyp,type_num,Ocmatrix,flag)

!        numpossible=1
!        do i=1,spgdata%NumMulti
!                ii=>ocm(i)
!                numpossible=numpossible*ii
!        enddo
    !print *, (spgdata%multi(i),i=1,spgdata%NumMulti)
!        print *, "ocm=",(ocmax(i),i=1,spgdata%NumMulti)
!        pause
!        if(.not.allocated(ocmatrix)) allocate(ocmatrix(27,numpossible))
!        ocmatrix=0
!        col=0
!        nlist=0
!        icount=1
!        wyc_oc_matrix=0
!        do inatom=1,ntype
!                nat=type_num(inatom)*spgdata%R
    !        print *, "nat=",nat
    !        print *, (spgdata%multi(i),i=1,spgdata%NumMulti)
    !        pause
!                icol=0
!                do i=1,numpossible
!                        dnatom=0
!                        do j=1,spgdata%NumMulti
!                                !                        print *, "ocmatrix",ocmatrix(j,i),spgdata%multi(j)
!                                dnatom=dnatom+ocmatrix(j,i)*spgdata%multi(j)
    !                        print *, "dnatom=",dnatom
!                        enddo
!                        if(dnatom==nat) then
!                                icol=icol+1
!                                wyc_oc_matrix(:,icol,inatom)=ocmatrix(:,i)
    !print *, inatom,"type",dnatom,"itype=",wyc_oc_matrix(:,icol,inatom)
    !pause
!                        endif
!                        if(icol>1000) exit
!                enddo
!                nlist(inatom)=icol
    !print *, "inatom=",inatom,icol
    !pause
!        enddo
!        do i=1,ntype
!                if (nlist(i)==0) then
!                        l_flag=.false.
!                        return
!                endif
!        enddo
!        l_flag= .true.
!        if(allocated(ocmatrix)) deallocate(ocmatrix)
!        end subroutine

    !subroutine add_node(pos,itype,wlist)

    !        USE SpgModule

    !        type(wyckoff),pointer ::pos,tmp
    !        integer,intent(in) :: wlist(27),itype
    !        allocate(tmp)
    !        tmp%itype=itype
    !        tmp%wlist(:)=wlist(:)
    !        if(associated(pos%next)) then
    !                tmp%next=>pos%next
    !                pos%next=>tmp
    !        else
    !                nullify(tmp%next)
    !                pos%next=>tmp
    !        end if
    !end subroutine
    !
    !subroutine del_node(pos)
    !        !        USE SpgModule
    !        type(wyckoff),pointer ::pos,next
    !        next=>pos%next
    !        if(associated(next%next)) then
    !                pos%next=>next%next
    !                deallocate(next)
    !        else
    !                nullify(pos%next)
    !                deallocate(next)
    !        end if
    !end subroutine
    !
    !subroutine show_all(pos)
    !        !        USE SpgModule
    !        type(wyckoff),pointer ::pos,tmp
    !        integer :: cnt
    !        cnt=1
    !        tmp=>pos
    !        do while(associated(tmp))
    !                print*,cnt,"th "
    !                print*, tmp%itype
    !                print *, tmp%wlist(:)
    !                tmp=>tmp%next
    !                cnt=cnt+1
    !        end do
    !end subroutine
    !
    !subroutine list2matrix(pos,ntype)
    !        !        USE SpgModule
    !
    !        type(wyckoff),pointer ::pos,tmp
    !        !integer,allocatable :: wycmatrix(:,:,:)
    !        integer,allocatable :: cnt(:)
    !        integer :: i,j,ntype,cn
    !        if(.not.allocated(cnt)) allocate(cnt(ntype))
    !        if(.not.allocated(nlist)) allocate(nlist(ntype))
    !        cnt=1
    !        tmp=>pos
    !        do while(associated(tmp))
    !                cnt(tmp%itype)=cnt(tmp%itype)+1
    !                tmp=>tmp%next
    !        end do
    !        cn=maxval(cnt)
    !        if(.not.allocated(wyc_oc_matrix)) allocate(wyc_oc_matrix(ntype,cn,27))
    !        cnt=0
    !        wyc_oc_matrix=0
    !        tmp=>pos
    !        do while(associated(tmp))
    !                do i=1,ntype
    !                        if (i==tmp%itype) then
    !                                cnt(i)=cnt(i)+1
    !                                wyc_oc_matrix(i,cnt(i),:)=tmp%wlist(:)
    !                        endif
    !                enddo
    !                tmp=>tmp%next
    !        end do
    !        nlist=cnt
    !        do j=1,ntype
    !                print *, "itype=",j
    !                do i=1,cnt(j)
    !                        print *, wyc_oc_matrix(j,i,:)
    !                enddo
    !        enddo
    !end subroutine

    subroutine gen_sg_num(sgindex)
        !!!!! author: Yanchao Wang
        !!!!! date  : 2012-03-16
        USE kinds, ONLY: DP, i4b
        !        USE SpgModule,         ONLY : spglist
        implicit none

        real(DP)                 :: tmp, rtmp, t
        integer(i4b)             :: notoccupied(230), i, j, sgindex, sg_temp
        integer(i4b), external    :: rscale
        Character(len=100)       :: filename
        logical                  :: alive, test_spg
        test_spg = .false.
        do while (.true.)
            sg_temp = rscale(1, 230)
            !print *, "sg_temp=",sg_temp
            if (spglist(sg_temp) .eqv. .true.) then
                sgindex = sg_temp
                exit
            end if
        end do
    end subroutine

END MODULE SpgModule
