       program kpath_high_symmetry_dire
       implicit real*8 (a-h,o-z)
cc       dimension path(20,3,10)
       dimension path(20,20,3)
       dimension hs_k(21,3),delta_K(20,3)
c       character*8 skn(21)
       open(unit=1,file="kpath.ini",status="unknown")
       open(unit=11,file="kpath.dat",status="unknown")
cc     reading the number of high symmetry directions in the BZ.       
       write(*,*) "please input the number of high symmetry directions"
       read(*,*)n
cc     reading in the high symmetry (n+1) k points 
       do i=1,n+1
         read(1,*)hs_k(i,1),hs_k(i,2),hs_k(i,3)
c         write(*,*)hs_k(i,1),hs_k(i,2),hs_k(i,3)
         if(i.gt.1) then 
          delta_K(i-1,1)=hs_k(i,1)-hs_k(i-1,1)
          delta_K(i-1,2)=hs_k(i,2)-hs_k(i-1,2)
          delta_K(i-1,3)=hs_k(i,3)-hs_k(i-1,3)      
         endif
       enddo

       write(11,100)hs_k(1,1),hs_k(1,2),hs_k(1,3)
       do i=1,n
          
          do j=1,20
          
          path(i,j,1)=delta_K(i,1)/20.0*(j)+hs_k(i,1)  
          path(i,j,2)=delta_K(i,2)/20.0*(j)+hs_k(i,2)
          path(i,j,3)=delta_K(i,3)/20.0*(j)+hs_k(i,3)
c         if(j.eq.10) then 
c         write(2,*)path(i,j,1),path(i,j,2),path(i,j,2),skn(i+1)
c         else
c          write(2,100)path(i,j,1),path(i,j,2),path(i,j,3),8.0
          write(11,100)path(i,j,1),path(i,j,2),path(i,j,3)
c         endif
          enddo
       enddo 
  100  format(5x,3(100f16.8))
       stop

       end
