FUNCTION Cal_Volume(matrix)

USE kinds,             ONLY : DP,i4b

IMPLICIT NONE

real(DP), intent(in)  :: matrix(3,3)
real(dp), external    :: det
real(DP) :: Cal_Volume

Cal_Volume = abs(Det(matrix))

END FUNCTION Cal_Volume

Function cal_Area(matrix)

USE kinds,             ONLY : DP,i4b

IMPLICIT NONE

real(DP), intent(in)  :: matrix(3,3)

real(DP) :: Cal_Area
cal_area=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)

cal_area = abs(cal_area)

end function cal_area
Function norm(a)
USE kinds,             ONLY : DP,i4b

IMPLICIT NONE
real(DP) :: a(3)
real(DP) :: norm
norm=sqrt(a(1)**2+a(2)**2+a(3)**2)
end function norm

FUNCTION Det(matrix)

USE kinds,             ONLY : DP,i4b

IMPLICIT NONE

real(DP),  intent(in)  :: Matrix(3,3)
real(DP) :: Det


Det = Matrix(1,1)*(Matrix(2,2)*Matrix(3,3)-Matrix(2,3)*Matrix(3,2))&
-Matrix(1,2)*(Matrix(2,1)*Matrix(3,3)-Matrix(2,3)*Matrix(3,1))&
+Matrix(1,3)*(Matrix(2,1)*Matrix(3,2)-Matrix(3,1)*Matrix(2,2))

END FUNCTION Det

subroutine CROSSPROD(A,B,X)
USE kinds,             ONLY : DP,i4b

IMPLICIT NONE

REAL(DP) :: A(3),B(3)
REAL(DP) :: X(3)

X(1)=A(2)*B(3)-B(2)*A(3)
x(2)=A(3)*B(1)-B(3)*A(1)
x(3)=A(1)*B(2)-B(1)*A(2)
END subroutine CROSSPROD

module linear

implicit none

contains
subroutine upper(matrix1,matrix2)
USE kinds ,  ONLY : DP,i4b
implicit none
real(DP) :: matrix1(:,:)
real(DP) :: matrix2(:,:)
integer(i4b) :: i,j,m,n
real :: a 
m=size(matrix1,1)
n=size(matrix1,2)
do i=1,n
   do j=i+1,n
	  a=matrix1(j,i)/matrix1(i,i)
	  matrix1(j,:)=matrix1(j,:)-a*matrix1(i,:)
	  matrix2(j,:)=matrix2(j,:)-a*matrix2(i,:)
   end do
end do 
end subroutine

subroutine lower(matrix1,matrix2)

USE kinds ,  ONLY : DP,i4b

implicit none

real(DP) :: matrix1(:,:)
real(DP) :: matrix2(:,:)
integer(i4b) :: i,j,m,n
real :: a
m=size(matrix1,1)
n=size(matrix1,2)
do i=n,2,-1
   do j=i-1,1,-1
	  a=matrix1(j,i)/matrix1(i,i)
	  matrix1(j,:)=matrix1(j,:)-a*matrix1(i,:)
	  matrix2(j,:)=matrix2(j,:)-a*matrix2(i,:)
   end do
end do 
end subroutine

end module

subroutine lat_inv(matrix3,matrix2)
USE kinds ,  ONLY : DP,i4b
use linear
implicit none
real(DP) :: matrix1(3,3)
real(DP) :: matrix2(3,3),matrix3(3,3)
integer(i4b) :: i,j
matrix1=matrix3
do i=1,3
   do j=1,3
	  if(i==j)then
		 matrix2(i,j)=1
	  else
		 matrix2(i,j)=0
	  end if
   end do 
end do
call upper(matrix1,matrix2)
call lower(matrix1,matrix2)
do i=1,3
   matrix2(i,:)=matrix2(i,:)/matrix1(i,i)
end do
end subroutine

subroutine mul_vm(r,mat,p)

USE kinds,            ONLY : DP,i4b

implicit none
real(DP) r(3),mat(3,3),p(3),r1(3)
integer(i4b) i,j
!  p=matmul(r,mat)
do i=1,3
   r1(i)=0.0
   do j=1,3
	  r1(i)=r1(i)+mat(j,i)*r(j)
   enddo
enddo
p(1)=r1(1)
p(2)=r1(2)
p(3)=r1(3)
end subroutine mul_vm

recursive subroutine Quick_sort(A,n,S,E)
!!!!the sort the number
USE kinds,            ONLY : DP,i4b
implicit none

integer(i4b) :: N,S,E,L,R,i
real(DP)     :: A(N),K,tmp
logical      bol
L=S
R=E
if(R<=L) return
R=E+1
K=A(S)
do while (.true.)
   do while (.true.)
	  L=L+1
	  if((A(L)>K).or.(L>=E))exit
   enddo
   do while(.true.)
	  R=R-1
	  if((A(R)<K).or.(R<=S))exit
   enddo
   if(R<=L)exit
   tmp=A(L)
   A(L)=A(R)
   A(R)=tmp
enddo
tmp=A(S)
A(S)=A(R)
A(R)=tmp
call Quick_sort(A,N,S,R-1)
call Quick_sort(A,N,R+1,E)

return
end subroutine Quick_sort

subroutine sortint_id(a,n,idx)
USE kinds,             ONLY : DP,i4b
implicit none
integer(i4b), intent(in)          :: n
integer(i4b), intent(in)          :: a(n)
integer(i4b), intent(out)         :: idx(n)
integer(i4b)                      :: i,j,k,l,m
real(8), parameter :: eps=1.d-14
do i=1,n
   idx(i)=i
end do
if (n.eq.1) return
l=n/2+1
k=n
10 continue
if (l.gt.1) then
   l=l-1
   m=idx(l)
else
   m=idx(k)
   idx(k)=idx(1)
   k=k-1
   if (k.eq.1) then
	  idx(1)=m
	  return
   end if
end if
i=l
j=l+l
20 continue
if (j.le.k) then
   if (j.lt.k) then
	  if (a(idx(j)).lt.a(idx(j+1))+eps) j=j+1
   end if
   if (a(m).lt.a(idx(j))+eps) then
	  idx(i)=idx(j)
	  i=j
	  j=j+j
   else
	  j=k+1
   end if
   goto 20
end if
idx(i)=m
goto 10
end subroutine

FUNCTION int_to_char( int )
!-----------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: int
CHARACTER (LEN=6)   :: int_to_char
!
!
IF ( int <= -10 ) THEN
   WRITE( UNIT = int_to_char , FMT = "(I3)" ) int
ELSE IF(int == 0 ) THEN
   WRITE(UNIT = int_to_char,FMT= "(I1)") INT
ELSE IF(int< 0 ) THEN
   WRITE( UNIT = int_to_char , FMT = "(I2)" ) int
ELSE IF ( int < 10 ) THEN
   WRITE( UNIT = int_to_char , FMT = "(I1)" ) int
ELSE IF ( int < 100 ) THEN
   !
   WRITE( UNIT = int_to_char , FMT = "(I2)" ) int
   !
ELSE IF ( int < 1000 ) THEN
   !
   WRITE( UNIT = int_to_char , FMT = "(I3)" ) int
ELSE IF ( int < 10000 ) THEN
   !
   WRITE( UNIT = int_to_char , FMT = "(I4)" ) int
   !
ELSE
   !
   WRITE( UNIT = int_to_char , FMT = "(I5)" ) int
   !
END IF
!
RETURN
!
END FUNCTION int_to_char

subroutine eign(A,N,R)

USE kinds,             ONLY : DP,i4b

IMPLICIT NONE

!! input and output parameters
real(DP),intent(inout)     :: A(N,N)
real(DP),intent(out)       :: R(N,N)
integer(i4b),intent(in)    :: N

!!! local parameters
integer(i4b)  ::  I,J,P,Q
real(DP)      ::  AMAX,TEMP,ZEMP,COO,SII,CO,SI,APP,AQQ,APQ,API,AQI
real(DP)      ::  RIP,RIQ

DO I=1,N
   DO J=1,N
      R(I,J)=0
   ENDDO
   R(I,I)=1
ENDDO
do while (.true.)
   AMAX=ABS(A(2,1))
   P=2
   Q=1
   DO I=2,N
      DO J=1,I-1
	 IF(ABS(A(I,J)).GT.AMAX) THEN
	    AMAX=ABS(A(I,J))
	    P=I
	    Q=J
	 ENDIF
      ENDDO
   ENDDO
   if(AMAX.LE.1.0E-7) then
      exit
   else 
      TEMP=2*A(P,Q)/(A(P,P)-A(Q,Q)+1.0e-30)
      ZEMP=(A(P,P)-A(Q,Q))/(2*A(P,Q))

      IF(ABS(TEMP).LT.1.0) THEN
	 COO=(1+TEMP**2)**(-0.5)
	 SII=TEMP*(1+TEMP**2)**(-0.5)
      ELSE
	 COO=ABS(ZEMP)*(1+ZEMP**2)**(-0.5)
	 SII=SIGN(1.0_DP,ZEMP)*(1+ZEMP**2)**(-0.5)
      ENDIF

      CO=SQRT(0.5*(1+COO))
      SI=SII/(2*CO)
      DO I=1,N
	 RIP=R(I,P)*CO+R(I,Q)*SI
	 RIQ=-R(I,P)*SI+R(I,Q)*CO
	 R(I,P)=RIP
	 R(I,Q)=RIQ
      ENDDO
      APP=A(P,P)*CO**2+A(Q,Q)*SI**2+2*A(P,Q)*CO*SI
      AQQ=A(P,P)*SI**2+A(Q,Q)*CO**2-2*A(P,Q)*CO*SI
      APQ=0.5*(A(Q,Q)-A(P,P))*SII+A(P,Q)*COO
      A(P,P)=APP
      A(Q,Q)=AQQ
      A(P,Q)=APQ
      A(Q,P)=A(P,Q)

      DO I=1,N
	 IF(I.EQ.P.OR.I.EQ.Q) THEN
	 ELSE
	    API=A(P,I)*CO+A(Q,I)*SI
	    AQI=-A(P,I)*SI+A(Q,I)*CO
	    A(P,I)=API
	    A(Q,I)=AQI
	    A(I,P)=A(P,I)
	    A(I,Q)=A(Q,I)
	 ENDIF
      ENDDO
   endif
enddo
END subroutine eign

function randomscale(lb,ub)

USE kinds,            ONLY : DP,i4b
implicit none
integer(i4b) :: lb,ub
integer(i4b) :: len
integer(i4b) :: randomscale, rt
real(dp)     :: t
logical      :: isok
logical      :: isprime
if (lb == ub) then
   randomscale = lb
   return
else
   isok = .true.
   len = ub - lb
   call random_number(t)
   rt = nint(lb + len*t)
   randomscale=rt
end if
end function randomscale

logical function isprime(num)
!! 
USE kinds,            ONLY : DP,i4b 
implicit none
integer(i4b) :: num, i, j
real(DP) :: t

t = num
j = aint(sqrt(t))
if (num < 4) then
   isprime = .false.
   return
else
   do i = 2, j
      if ( mod(num, i) == 0 ) then
	 isprime = .false.
	 return
      end if
   end do
end if

isprime = .true.
return
end function isprime

subroutine MATVECTMUL(mat,vect,mvec)
USE kinds,             ONLY : DP,i4b

IMPLICIT NONE
INTEGER :: i
REAL(DP) :: mat(3,3)
REAL(DP) :: vect(3),mvec(3)

DO i=1,3
   mvec(i)=SUM(mat(i,:)*vect)
ENDDO
END subroutine MATVECTMUL

function rscale(lb,ub)

USE kinds,            ONLY : DP,i4b
implicit none
integer(i4b) :: lb,ub
integer(i4b) :: len
integer(i4b) :: rscale, rt
real(dp)     :: t
if (lb == ub) then
   rscale = lb
   return
else
   len=ub-lb
   call random_number(t)
   rt = lb + nint(len*t)
   rscale=rt
   return
end if
end function rscale

FUNCTION DAB(A,B)
!
USE kinds ,           ONLY : DP,i4b
!
implicit none

REAL(DP) :: A(3),B(3)
REAL(DP) :: DAB 
DAB=DOT_PRODUCT(A,B)
DAB=SQRT(dabs(DAB))
END
function DET_I3(a)
USE kinds,            ONLY : DP,i4b

implicit none

integer(i4b) :: a(3,3)
integer(i4b) :: det_I3
det_I3=0
det_I3=det_I3+ a(1,1) * (a(2,2) * a(3,3) - a(2,3) * a(3,2)) &
+ a(1,2) * (a(2,3) * a(3,1) - a(2,1) * a(3,3)) &
+ a(1,3) * (a(2,1) * a(3,2) - a(2,2) * a(3,1))
end

subroutine dis_mat(A,B,dist) 
USE KINDS,           ONLY : DP, i4b 
implicit none
integer(i4b)    :: i,j,k
real(DP)        :: A(3,3),B(3,3)
real(DP)        :: dist
dist=0.0
do i=1,3
   do j=1,3
      dist=dist+abs(A(i,j)-B(i,j))
   enddo
enddo
end subroutine

subroutine inv_mat(A,invA)

USE KINDS,           ONLY : DP, i4b

implicit none

integer(i4b)    :: i,j,k,l,is(3),js(3)
real(DP)        :: A(3,3),invA(3,3),A_old(3,3)
real(DP)        :: t,d
A_old=A
l=1
do k=1,3
   d=1D-3
   do i=k,3
      do j=k,3
	 if (abs(a(i,j)).gt.d) then
	    d=abs(a(i,j))
	    is(k)=i
	    js(k)=j
	 endif
      enddo
   enddo
   if ((d+1.0).eq.1.0) then
      l=0
   endif
   do j=1,3
      t=a(k,j)
      a(k,j)=a(is(k),j)
      a(is(k),j)=t
   enddo
   do i=1,3
      t=a(i,k)
      a(i,k)=a(i,js(k))
      a(i,js(k))=t
   enddo
   a(k,k)=1/a(k,k)
   do j=1,3
      if(j.ne.k) then
	 a(k,j)=a(k,j)*a(k,k)
      endif
   enddo
   do i=1,3
      if (i.ne.k) then
	 do j=1,3
	    if (j.ne.k) then
	       a(i,j)=a(i,j)-a(i,k)*a(k,j)
	    endif
	 enddo
      endif
   enddo
   do i=1,3
      if(i.ne.k) then
	 a(i,k)=-a(i,k)*a(k,k)
      endif
   enddo
enddo
do k=3,1,-1
   do j=1,3
      t=a(k,j)
      a(k,j)=a(js(k),j)
      a(js(k),j)=t
   enddo
   do i=1,3
      t=a(i,k)
      a(i,k)=a(i,is(k))
      a(i,is(k))=t
   enddo
enddo
invA=A
A=A_old
end subroutine

subroutine lat2matrix(lat,matrix,iflag)
!
use kinds    ,   only : i4b,DP
!
implicit none
!
integer(i4b),intent(in)    :: iflag ! if iflag==1, abc2matrix; iflag==2. matrix2abc
real(DP),    intent(inout) :: lat(6),matrix(3,3)

!local parameters
real(DP)     :: ra,rb,rc,cosinea,cosineb,cosinec,anglea,angleb,anglec

if (iflag==1) then
   matrix=0.0
   matrix(1,1) = lat(1)
   matrix(2,1) = lat(2)*cos(lat(6))
   matrix(2,2) = lat(2)*sin(lat(6))
   matrix(3,1) = lat(3)*cos(lat(5))
   matrix(3,2) = lat(3)*cos(lat(4))*sin(lat(6))-((lat(3)*cos(lat(5))&
   -lat(3)*cos(lat(4))*cos(lat(6)))/tan(lat(6)))
   matrix(3,3) = sqrt(lat(3)**2 -matrix(3,1)**2 - matrix(3,2)**2)
else
   lat=0.0
   ra=sqrt(matrix(1,1)**2+matrix(1,2)**2+matrix(1,3)**2)
   rb=sqrt(matrix(2,1)**2+matrix(2,2)**2+matrix(2,3)**2)
   rc=sqrt(matrix(3,1)**2+matrix(3,2)**2+matrix(3,3)**2)
   cosinea=(matrix(2,1)*matrix(3,1)+matrix(2,2)*matrix(3,2)+matrix(2,3)*matrix(3,3))/rb/rc
   cosineb=(matrix(1,1)*matrix(3,1)+matrix(1,2)*matrix(3,2)+matrix(1,3)*matrix(3,3))/ra/rc
   cosinec=(matrix(1,1)*matrix(2,1)+matrix(1,2)*matrix(2,2)+matrix(1,3)*matrix(2,3))/ra/rb
   anglea=acos(cosinea)
   angleb=acos(cosineb)
   anglec=acos(cosinec)
   lat(1)=ra
   lat(2)=rb
   lat(3)=rc
   lat(4)=anglea
   lat(5)=angleb
   lat(6)=anglec
endif
end subroutine lat2matrix
subroutine check_lat(logic,lat1_matrix)
!
use kinds  ,     only : i4b,DP
!
implicit none
!
integer(i4b)  :: i,j,k
real(kind=DP) :: lat1_matrix(3,3),lat_matrix_prim(3,3),lx,ly,lz,xy,xz,yz
real(kind=DP) :: lat_matrix(3,3),ra,rb,rc,alpha,beta,gama,cosinea,cosineb,cosinec
logical       :: logic
!

lat_matrix=lat1_matrix

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
!print *, alpha,beta,gama
!pause
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
!  if (Icode==9)  then
!    call abc_to_triclinic_box(ra,rb,rc,alpha*3.1415926/180.0,beta*3.1415926/180.0,gama*3.1415926/180.0,&
!                           & lx, ly, lz, xy, xz, yz)
!    if(xy>0.5*lx.or.xy<-0.5*lx) then
!       logic=.false.
!       return
!    endif
!    if(xz>0.5*lx.or.xz<-0.5*lx) then
!       logic=.false.
!       return
!    endif
!    if(yz>0.5*ly.or.yz<-0.5*ly) then
!       logic=.false.
!       return
!    endif
!   endif
end subroutine check_lat

subroutine cal_vect(a,lat)
!
USE kinds ,           ONLY : DP,i4b
!
implicit none
! This program is writed to calculate the relationship between the vectors in output_file of PWSCF
real(DP) :: a(3,3)!(input)
real(DP) :: lat(6)!(output)
!local parameter
real(DP) :: ra,rb,rc,cosinea,cosineb,cosinec,anglea,angleb,anglec

ra=sqrt(a(1,1)**2+a(1,2)**2+a(1,3)**2)
rb=sqrt(a(2,1)**2+a(2,2)**2+a(2,3)**2)
rc=sqrt(a(3,1)**2+a(3,2)**2+a(3,3)**2)
cosinea=(a(2,1)*a(3,1)+a(2,2)*a(3,2)+a(2,3)*a(3,3))/rb/rc
cosineb=(a(1,1)*a(3,1)+a(1,2)*a(3,2)+a(1,3)*a(3,3))/ra/rc
cosinec=(a(1,1)*a(2,1)+a(1,2)*a(2,2)+a(1,3)*a(2,3))/ra/rb
anglea=acos(cosinea)
angleb=acos(cosineb)
anglec=acos(cosinec)
lat(1)=ra
lat(2)=rb
lat(3)=rc
lat(4)=(anglea/3.1415926)*180.0
lat(5)=(angleb/3.1415926)*180.0
lat(6)=(anglec/3.1415926)*180.0
end subroutine cal_vect

