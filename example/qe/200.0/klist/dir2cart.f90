program dir2cart

implicit none

real*8 cart(3,153),dir(3,153),vec(3,3)
integer i,j,k,natom

natom=8

open(101, file="input.dat")

do i=1,3
    read(101,*)(vec(j,i),j=1,3)
!    print *, "vec", (vec(j,i),j=1,3)
end do
!print *, "vec", vec
do i=1,natom
    read(101,*)(dir(j,i),j=1,3)
end do

do j=1,natom
    do i=1,3
        cart(i,j)=dir(1,j)*vec(i,1)+dir(2,j)*vec(i,2)+dir(3,j)*vec(i,3)
    end do
end do

print *, "Below cart:"
do i=1,natom
    write(*,"(3f12.7)")cart(1,i),cart(2,i),cart(3,i)
end do

close(101)
stop
end
