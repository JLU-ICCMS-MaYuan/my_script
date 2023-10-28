
! Copyright (C) 2011 A. Sanna and E. K. U. Gross
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readalpha2f(w,a2f,ndos)
use modmain
implicit none
! arguments
real(8), intent(out) :: w(nwdos)
real(8), intent(out) :: a2f(nwdos)
! local variables
integer iw,iostat,ndos
real(8) a2fw, ww
open(50,file='ALPHA2F.OUT',action='READ',form='FORMATTED',status='OLD', &
 iostat=iostat) !如果iostat的值为0，表示文件成功打开，否则，将有相应的错误代码表示文件打开失败。
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readalpha2f): error opening ALPHA2F.OUT")')
  write(*,*)
  stop
end if
iw = 0
100  read(50,*,end=200,iostat=iostat) ww  ,a2fw   
!*：这表示读取文件时采用默认的格式。在这种情况下，不指定特定的格式。
! iostat=iostat：这用于获取有关文件操作状态的信息。
! iostat是一个整数变量，用于记录读取操作是否成功。
! 如果iostat的值为0，表示读取操作成功，否则，它将包含有关读取错误的信息。  
  iw = iw + 1
  w(iw) = ww
  a2f(iw) = a2fw
  if (iostat.ne.0) then
    write(*,*)
    write(*,'("Error(readalpha2f): error reading from ALPHA2F.OUT")')
    write(*,'(" for frequency ",I6)') iw
    write(*,*)
    stop
  end if
go to 100
200 ndos = iw
    write(*,'("Number of points in Eliashberg spectral function =",I5)') ndos
close(50)
return
end subroutine

