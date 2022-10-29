program read_band
  !      This program is to get the suitable data for ploting bands of energy
  real kx(300),ky(300),kz(300)
  real y(300,51),gam(300,51),lamb(300,51)
  character ch1,ch2,ch3,ch4
  real Rytocm1,RytoTHz,epsw,pi
  Rytocm1 = 109737.57990d0
  RytoTHz = 3.289828D3
  epsw = 1/Rytocm1
  pi = 3.14159265358979d0
  dosfit = 7.444133
  open(1,file='.freq')
  open(2,file='bands.dat')
  open(4,file='bb.dat')
  !      This file intend to generate the high symmetry point
  open(3,file='spoint.dat')
  read(1,*)ch1,ch2,nband,ch3,nks,ch4
  print *, nks
  print *, nband 
  !       read(1,*)emin,emax

  print *,  "读入 .freq"
  ! kx ky kz 代表着 nks 的三个方向
  ! y 代表一个k点的 nband 个能带
  do i=1,nks
     read(1,*)kx(i),ky(i),kz(i)
     read(1,*)(y(i,j),j=1,nband)
  enddo
  !       write(*,*)"nks=",nks,"nband=",nband

  print *, "读入 bb.dat,
其中nks代表着总的k点个数，一般为高对称的k点数+插入的k点个数，例如有8个高对称的k点，每两个k点之间插入20个k点，那么总的nks=1+20*7=141.
这里nks=", nks
  do i=1,nks
     read(4,*)tmp,(gam(i,j),j=1,nband)
     !print *, "tmp=", tmp
     !print *, "gam=", (gam(i,j),j=1,nband)
     !write(*,*)i
     !pause
  enddo
  print *, "结束了bb.dat" 

  print *, "将gam转化为THz"
  do i = 1,nks
     do j = 1,nband
        if (gam(i,j) > 0) then
           lamb(i,j) = gam(i,j)/pi/RytoTHz/dosfit/(y(i,j)/Rytocm1)**2
           !print *, "lamb(i, j)", lamb(i, j) 
           !print *, i,j
           !pause
           lamb(i,j) = sqrt(lamb(i,j)/pi)
        else
           lamb(i,j) = 0.0d0
!           print *, "lamb(i, j)", lamb(i, j) 
!           pause
           gam(i,j) = -1*gam(i,j)
        endif
     enddo
  enddo
  !      This part intend to chose the energy range
  !       do i=1,nks
  !         do j=1,nband
  !           if (y(i,j).lt.emin) then
  !             y(i,j)=emin
  !           else if (y(i,j).gt.emax) then
  !             y(i,j)=emax
  !           endif
  !         enddo
  !       enddo

  print *, "输出x"
  do j = 1, nband
     x=0.0
     x_tot=0.0
     do i=1,nks
        if(i.gt.1) then
           x=sqrt((kx(i)-kx(i-1))**2+(ky(i)-ky(i-1))**2+(kz(i)-kz(i-1))**2)
        endif
        x_tot=x_tot+x
        if(mod(i,21).eq.0) then
           write(3,'(100F10.4)')x_tot
        endif
        !         gam(i,j)=gam(i,j)/RytoTHz
        write(2,'(10(F12.6,2x))')x_tot,y(i,j),lamb(i,j),gam(i,j)
     enddo
     write(2,*)'                                         '
  enddo
  stop
end program read_band

