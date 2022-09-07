Module SpgCryLat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Filename:SpgCryLat.F90                                                  
    !>      Programmer            Email                   Web  
    !>      ==================    ======================= ========================
    !>      Yanchao Wang          wyc@calypso.cn          www.calypso.cn
    !> Record of revisions:                                                       
    !>      Date          Programmer        Description of changes                
    !>      ==========    ==============    ======================================
    !>      2017.05.10    Yanchao Wang      First build the module
    !> Discription:
    !>      This module is designed to generate the lattice matrix of crystal structures 
    !> Parameters:
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

	subroutine SpgGenLat(vol,sgindex,lat_mat)
		
		USE kinds,            ONLY : DP, i4b
        USE SpgDataInfo

		implicit none

		! input and output parameters
		real(DP), intent(in)     :: vol
		integer(i4b), intent(in) :: sgindex 
		real(DP), intent(out)    :: lat_mat(3,3)

		! local parameters
		logical :: logic 
        real(DP) :: volume
		logic=.false.
		do while (.not.logic)
			volume=vol*spacdata(sgindex)%R
		!	print *, "vol=",volume
		!	print *, spacdata(sgindex)%R
			call GenSymLat(sgindex,volume,lat_mat)
			call CL_CheckLattice(lat_mat,logic)
		enddo
	end subroutine 

	subroutine GenSymLat(sgindex,Vol,lat_mat)

		USE kinds,             ONLY : DP,i4b 

		implicit none

		! input and output
		real(DP), intent(in)    :: Vol
		integer(i4b),intent(in) :: sgindex
		real(DP), intent(out)   :: lat_mat(3,3)

		! local parameter
		real(DP)     :: random_lat(6)
		!
		! determine the pg and transition operations
		!
		!  call random_seed
		lat_mat=0.0

		!  if (cluster) then
		!     lat_mat(1,1)=vol**0.33333
		!     lat_mat(2,2)=vol**0.33333
		!     lat_mat(3,3)=vol**0.33333
		!     return
		!  end if

		select case (sgindex)
		case ( 1:2 )
			do while(.true.)
				call random_number(random_lat(1:6))
				if (abs(random_lat(1)/random_lat(2))>0.1.or. abs(random_lat(1)/random_lat(2))<10.0& 
					.and. abs(random_lat(1)/random_lat(3))>0.1.or.abs(random_lat(1)/random_lat(3))<10.0 &
					.and. abs(random_lat(2)/random_lat(3))>0.1.or.abs(random_lat(2)/random_lat(3))<10.0 ) exit

			end do
			random_lat(4:6)=random_lat(4:6)*3.1415926/2.0
			call latconver(random_lat,lat_mat,vol)
		case ( 3:15  )
			do while(.true.)
				call random_number(random_lat(1:6))
				if (abs(random_lat(1)/random_lat(2))>0.1.or. abs(random_lat(1)/random_lat(2))<10.0& 
					.and. abs(random_lat(1)/random_lat(3))>0.1.or.abs(random_lat(1)/random_lat(3))<10.0 &
					.and. abs(random_lat(2)/random_lat(3))>0.1.or.abs(random_lat(2)/random_lat(3))<10.0 ) exit

			end do
			random_lat(4)=0.5*3.1415926
			random_lat(6)=0.5*3.1415926
			random_lat(5)=random_lat(5)*3.1415926/2.0
			call latconver(random_lat,lat_mat,vol)
		case ( 16:74 )
			do while(.true.)
				call random_number(random_lat(1:6))
				if (abs(random_lat(1)/random_lat(2))>0.1.or. abs(random_lat(1)/random_lat(2))<10.0& 
					.and. abs(random_lat(1)/random_lat(3))>0.1.or.abs(random_lat(1)/random_lat(3))<10.0 &
					.and. abs(random_lat(2)/random_lat(3))>0.1.or.abs(random_lat(2)/random_lat(3))<10.0 ) exit

			end do
			random_lat(4)=0.5
			random_lat(5)=0.5
			random_lat(6)=0.5
			random_lat(4:6)=random_lat(4:6)*3.1415926
			call latconver(random_lat,lat_mat,vol)
		case ( 75:142 )
			do while(.true.)
				call random_number(random_lat(1:6))
				if (abs(random_lat(1)/random_lat(2))>0.1.or. abs(random_lat(1)/random_lat(2))<10.0& 
					.and. abs(random_lat(1)/random_lat(3))>0.1.or.abs(random_lat(1)/random_lat(3))<10.0 &
					.and. abs(random_lat(2)/random_lat(3))>0.1.or.abs(random_lat(2)/random_lat(3))<10.0 ) exit

			end do
			random_lat(1)=random_lat(2)
			random_lat(4:6)=0.5
			random_lat(4:6)=random_lat(4:6)*3.1415926
			call latconver(random_lat,lat_mat,vol)
		case (143:167 )
			do while(.true.)
				call random_number(random_lat(1:6))
				if (abs(random_lat(1)/random_lat(2))>0.1.or. abs(random_lat(1)/random_lat(2))<10.0& 
					.and. abs(random_lat(1)/random_lat(3))>0.1.or.abs(random_lat(1)/random_lat(3))<10.0 &
					.and. abs(random_lat(2)/random_lat(3))>0.1.or.abs(random_lat(2)/random_lat(3))<10.0 ) exit

			end do
			random_lat(1)=random_lat(2)
			random_lat(5)=0.5
			random_lat(4)=random_lat(5)
			random_lat(6)=0.6666667
			random_lat(4:6)=random_lat(4:6)*3.1415926
			call latconver(random_lat,lat_mat,vol)
		case (168:194 )
			do while(.true.)
				call random_number(random_lat(1:6))
				if (abs(random_lat(1)/random_lat(2))>0.1.or. abs(random_lat(1)/random_lat(2))<10.0& 
					.and. abs(random_lat(1)/random_lat(3))>0.1.or.abs(random_lat(1)/random_lat(3))<10.0 &
					.and. abs(random_lat(2)/random_lat(3))>0.1.or.abs(random_lat(2)/random_lat(3))<10.0 ) exit

			end do
			random_lat(1)=random_lat(2)
			random_lat(5)=0.5
			random_lat(4)=random_lat(5)
			random_lat(6)=0.6666667
			random_lat(4:6)=random_lat(4:6)*3.1415926
			call latconver(random_lat,lat_mat,vol)
		case (195:230 )
			do while(.true.)
				call random_number(random_lat(1:6))
				if (abs(random_lat(1)/random_lat(2))>0.1.or. abs(random_lat(1)/random_lat(2))<10.0& 
					.and. abs(random_lat(1)/random_lat(3))>0.1.or.abs(random_lat(1)/random_lat(3))<10.0 &
					.and. abs(random_lat(2)/random_lat(3))>0.1.or.abs(random_lat(2)/random_lat(3))<10.0 ) exit

			end do
			random_lat(1)=random_lat(2)
			random_lat(3)=random_lat(2)
			random_lat(4:6)=0.5
			random_lat(4:6)=random_lat(4:6)*3.1415926
			call latconver(random_lat,lat_mat,vol)
		end select
		!if(D2C) then
		!    lat_mat(1,3)=0.0
		!    lat_mat(2,3)=0.0
		!    lat_mat(3,3)=Depth
		!    lat_mat(3,1)=0.0
		!    lat_mat(3,2)=0.0
		!endif
	end subroutine 

	subroutine LatConver(random_lat,lat_mat1,vol)
		USE kinds,            ONLY : DP,i4b

		implicit none

		real(DP) :: randV,lat1(6),random_lat(6),random_matrix(3,3),ratio,lat_mat1(3,3),rand_lat(6),vol

		rand_lat=random_lat
		call lat2matrix(rand_lat,random_matrix,1)
		randV=0.0
		randV=(random_matrix(1,2)*random_matrix(2,3)-random_matrix(1,3)*random_matrix(2,2))*random_matrix(3,1)
		randV=randV+(random_matrix(1,3)*random_matrix(2,1)-random_matrix(1,1)*random_matrix(2,3))*random_matrix(3,2)
		randV=randV+(random_matrix(1,1)*random_matrix(2,2)-random_matrix(1,2)*random_matrix(2,1))*random_matrix(3,3)

		ratio=vol/randV
		lat1(1)=random_lat(1)*(ratio**(0.3333333333))
		lat1(2)=random_lat(2)*(ratio**(0.3333333333))
		lat1(3)=random_lat(3)*(ratio**(0.3333333333))
		lat1(4)=random_lat(4)
		lat1(5)=random_lat(5)
		lat1(6)=random_lat(6)
		call lat2matrix(lat1,lat_mat1,1)
	end subroutine 
subroutine CL_CheckLattice(lat_matrix,logic)
  !
  USE kinds  ,     ONLY : i4b,DP
  !
  implicit none
  !
  integer(i4b)  :: i,j,k
  real(DP) :: lat_matrix(3,3),ra,rb,rc,alpha,beta,gama,cosinea,cosineb,cosinec
  logical       :: logic
  !
  !print *, lat_matrix
  !pause
  do i = 1, 3
     do j = 1, 3
        if (isnan(lat_matrix(i,j))) then
         logic=.false.
        return
        endif
     end do
  end do
     ra=sqrt(lat_matrix(1,1)**2+lat_matrix(1,2)**2+lat_matrix(1,3)**2)
     rb=sqrt(lat_matrix(2,1)**2+lat_matrix(2,2)**2+lat_matrix(2,3)**2)
     rc=sqrt(lat_matrix(3,1)**2+lat_matrix(3,2)**2+lat_matrix(3,3)**2)
     cosinea=(lat_matrix(2,1)*lat_matrix(3,1)+lat_matrix(2,2)*lat_matrix(3,2)+lat_matrix(2,3)*lat_matrix(3,3))/rb/rc
     cosineb=(lat_matrix(1,1)*lat_matrix(3,1)+lat_matrix(1,2)*lat_matrix(3,2)+lat_matrix(1,3)*lat_matrix(3,3))/ra/rc
     cosinec=(lat_matrix(1,1)*lat_matrix(2,1)+lat_matrix(1,2)*lat_matrix(2,2)+lat_matrix(1,3)*lat_matrix(2,3))/ra/rb
     alpha=(acos(cosinea)/3.1415926)*180.d0
     beta=(acos(cosineb)/3.1415926)*180.d0
     gama=(acos(cosinec)/3.1415926)*180.d0
!	 print *, alpha,beta,gama
     logic=.true.
     if (ra<1.0 .or. rb<1.0 .or. rc<1.0) then
        logic=.false.
        return
     end if
     if (alpha<20.0 .or. alpha>160.0) then
        logic=.false.
        return
     end if
     if (beta<20.0 .or. beta>160.0 ) then
        logic=.false.
        return
     end if
     if (gama< 20.0 .or.gama>160.0 ) then
        logic=.false.
        return
     end if
     if ( ra/rb >5.0 .or. ra/rb<0.2 ) then
        logic=.false.
        return
     end if
     if ( ra/rc >5.0 .or. ra/rc<0.2 ) then
        logic=.false.
        return
     end if
     if ( rb/rc >5.0 .or. rb/rc<0.2 ) then
        logic=.false.
        return
     end if
end subroutine CL_CheckLattice 

END MODULE
