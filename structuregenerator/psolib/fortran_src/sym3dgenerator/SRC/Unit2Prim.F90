subroutine Unit2Prim(natom,C, lat, pos, tau)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!> FiLEName: u2p.F90                                                  
!!>      Programmer            Email                   Web  
!!>      ==================    ======================= ========================
!!>      Yanchao Wang          wyc@calypso.cn          www.calypso.cn
!!>      Date          Programmer        Description of changes                
!!>      ==========    ==============    ======================================
!!>      2017.05.15    WYC               create
!!> Discription:
!!>      This module is designed to transfer unit cell to primitive cell
!!> Parameters:
!!>       
!!> References:              
!!>          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

USE  kinds, ONLY : DP, i4b

implicit none

integer(i4b), intent(in)         :: natom
character, intent(in)     :: C

real(dp), dimension(3, 3)     :: lat,plat, invlat,invtau
real(dp), dimension(3, natom) :: tau
real(DP), dimension(3,4*natom):: pos,ptau
real(dp), dimension(3, 3) :: lp, ti, tf, ta, tmc, toc, li, lf, la, lmc, loc
integer(i4b) :: i

data lp /&
& 1.0, 0.0,  0.0,&
& 0.0, 1.0,  0.0,&
& 0.0, 0.0,  1.0/

data ti /&
& 0.0, 1.0,  1.0,&
& 1.0, 0.0,  1.0,&
& 1.0, 1.0,  0.0/

data tf /&
&-1.0, 1.0,  1.0,&
& 1.0,-1.0,  1.0,&
& 1.0, 1.0, -1.0/

data ta /&
& 0.0, 1.0,  1.0,&
& 0.0,-1.0,  1.0,&
& 1.0, 0.0,  0.0/

data tmc /&
& 1.0,-1.0,  0.0,&
& 1.0, 1.0,  0.0,&
& 0.0, 0.0,  1.0/

data toc /&
& 1.0, 1.0,  0.0,&
&-1.0, 1.0,  0.0,&
& 0.0, 0.0,  1.0/

data li /&
& -0.5,  0.5,  0.5, &
&  0.5, -0.5,  0.5, &
&  0.5,  0.5, -0.5/

data lf /&
& 0.0,   0.5,   0.5, &
& 0.5,   0.0,   0.5, &
& 0.5,   0.5,   0.0/

data la /&
& 0.0, 0.5, 0.5, &
& 0.0, -0.5,  0.5, &
& 1.0,  0.0,  0.0/

data lmc /&
&  0.5,   -0.5,   0.0, &
&  0.5,   0.5,   0.0, &
&  0.0,   0.0,   1.0/

data loc /&
&  0.5,   0.5,  0.0, &
&  -0.5,  0.5,  0.0, &
&  0.0,   0.0,   1.0/

select case (C)
case ('I')
	invlat=li
	invtau=ti
	!call unit_lat_2_prim_lat(li, ulat, plat)
	!call unit_tau_2_prim_tau(natom,plat, ti, utau, ptau)
case ('F')
	invlat=lf
	invtau=tf
case ('A')
	invlat=la
	invtau=ta
case ('C')
	invlat=loc
	invtau=toc
case ('M')
	invlat=lmc
	invtau=tmc
case ('P')
	invlat=lp
	invtau=lp
end select
call unit_lat_2_prim_lat(invlat, lat, plat)
call unit_tau_2_prim_tau(4*natom, invtau, pos, ptau)
lat=plat
do i=1,natom
	tau(:,i)=ptau(:,i)
enddo
end subroutine 

subroutine unit_lat_2_prim_lat(T, lat, plat)
use kinds, only : i4b, dp
implicit none
real(dp), dimension(3, 3) :: lat, plat, T
real(dp) :: tmp
integer(i4b) :: i, j, k

plat = 0.

do i = 1, 3
	do j = 1, 3
		do k = 1, 3
			plat(i, j) = plat(i, j) + T(k, i) * lat(k, j)
		end do
	end do
end do


end subroutine unit_lat_2_prim_lat

subroutine unit_tau_2_prim_tau(natom, T, utau, ptau)
use kinds, only : i4b, dp
implicit none

integer(i4b), intent(in)     :: natom

real(dp), dimension(3, 3)     ::  T
real(dp), dimension(3, natom) :: utau, ptau, p
real(DP) :: diff
logical :: flag
integer(i4b) :: i, ityp, iatom, j, k


do iatom = 1, natom 
	call mul_vm(utau(:, iatom), T, p(:, iatom))
	do i = 1, 3
		p(i, iatom) = p(i, iatom) - floor(p(i, iatom))
	end do
end do
iatom=0
do i=1,natom
	if (iatom==0) then
		iatom=iatom+1
		ptau(:,iatom)=p(:,i)
	else
		flag=.true.
		do j=1,iatom
			diff=0.d0
			do k=1,3
				diff=diff+abs(ptau(1,j)-p(1,i))+abs(ptau(2,j)-p(2,i))+abs(ptau(3,j)-p(3,i))
			enddo
			if(abs(diff)<0.000001d0) then
				flag=.false.
				exit
			endif
		enddo
		if(flag) then
			iatom=iatom+1
			ptau(:,iatom)=p(:,i)
		endif
	endif
enddo
end subroutine unit_tau_2_prim_tau
