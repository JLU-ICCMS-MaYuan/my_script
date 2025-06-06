
! Copyright (C) 2011 J. K. Dewhurst, A. Sanna, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mcmillan(w,a2f,lambda,wlog,wrms,tc,ndos)
use modmain
use modphonon
implicit none
! arguments
real(8), intent(in) :: w(nwdos)
real(8), intent(in) :: a2f(nwdos)
real(8), intent(out) :: lambda
real(8), intent(out) :: wlog
real(8), intent(out) :: wrms
real(8), intent(out) :: tc
! local variables
integer iw, ndos
real(8) l1,l2,f1,f2,t1
! allocatable arrays
real(8), allocatable :: f(:),g(:),cf(:,:)
allocate(f(ndos),g(ndos),cf(4,ndos))
! compute the total lambda
write(6,'("ndos =",i5)') ndos
do iw=1,ndos
  if (w(iw).gt.1.d-8) then
    f(iw)=a2f(iw)/w(iw)
  else
    f(iw)=0.d0
  end if
end do
call fderiv(-3,ndos,w,f,g,cf)
lambda=2.d0*g(ndos)
write(6,'("lambda =",f10.5)') lambda
! compute the logarithmic average frequency
do iw=1,ndos
  if (w(iw).gt.1.d-8) then
    f(iw)=a2f(iw)*log(w(iw))/w(iw)
  else
    f(iw)=0.d0
  end if
end do
call fderiv(-3,ndos,w,f,g,cf)
t1=(2.d0/lambda)*g(ndos)
wlog=exp(t1)
! compute < w^2 >^(1/2)
do iw=1,ndos
  if (w(iw).gt.1.d-8) then
    f(iw)=a2f(iw)*w(iw)
  else
    f(iw)=0.d0
  end if
end do
call fderiv(-3,ndos,w,f,g,cf)
t1=(2.d0/lambda)*g(ndos)
wrms=sqrt(abs(t1))
! compute McMillan-Allen-Dynes superconducting critical temperature
t1=(-1.04d0*(1.d0+lambda))/(lambda-mustar-0.62d0*lambda*mustar)
tc=(wlog/(1.2d0*kboltz))*exp(t1)
l1=2.46d0*(1.d0+3.8d0*mustar)
l2=1.82d0*(1.d0+6.3d0*mustar)*(wrms/wlog)
f1=(1.d0+(lambda/l1)**(3.d0/2.d0))**(1.d0/3.d0)
f2=1.d0+(wrms/wlog-1.d0)*(lambda**2)/(lambda**2+l2**2)
tc=tc*f1*f2
deallocate(f,g,cf)
return
end subroutine

