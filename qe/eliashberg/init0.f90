
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: init0
! !INTERFACE:
subroutine init0
! !USES:
use modmain
!use modxcifc
!use modldapu
! !DESCRIPTION:
!   Performs basic consistency checks as well as allocating and initialising
!   global variables not dependent on the $k$-point set.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias
integer ist,l,m,lm,iv(3)
real(8) rsum,t1
real(8) ts0,ts1

!-------------------------------!
!     zero timing variables     !
!-------------------------------!
timeinit=0.d0
timemat=0.d0
timefv=0.d0
timesv=0.d0
timerho=0.d0
timepot=0.d0
timefor=0.d0
call timesec(ts0)

!------------------------------------!
!     index to atoms and species     !
!------------------------------------!
natmmax=0
ias=0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=ias+1
    idxas(ia,is)=ias
    idxis(ias)=is
    idxia(ias)=ia
  end do
! maximum number of atoms over all species
  natmmax=max(natmmax,natoms(is))
end do
! total number of atoms
natmtot=ias

return
end subroutine
!EOC

