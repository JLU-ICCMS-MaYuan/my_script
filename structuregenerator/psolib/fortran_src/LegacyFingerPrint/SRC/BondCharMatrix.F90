MODULE BondCharMatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Filename: BondCharMatrix.F90
!>      Programmer            Email                   Web
!>      ==================    ======================= ========================
!>      Jian Lv               lvjian@calypso.cn       www.calypso.cn
!> Record of revisions:
!>      Date          Programmer        Description of changes
!>      ==========    ==============    ======================================
!>      2017.05.18    Jian Lv           Create
!> Discription:
!>      Module for structure characterizaiton using the BCM.
!> Parameters:
!>
!> References:
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    USE KINDS, ONLY: DP, I4B
    USE CONSTANTS, ONLY: PI

    IMPLICIT NONE

CONTAINS

    SUBROUTINE BondCharMatrix_CalBCM(structyp, nele, natom, ele_num, lat_matrix, tau, sim_func)
        !
        !input and output
        !
        INTEGER(I4B), INTENT(IN)      :: structyp   ! The type of input structures.
        ! 1 for crystal
        ! 2 for cluster
        INTEGER(I4B), INTENT(IN)      :: nele
        INTEGER(I4B), INTENT(IN)      :: natom
        INTEGER(I4B), INTENT(IN)      :: ele_num(nele)
        REAL(DP), INTENT(IN)      :: lat_matrix(3, 3)
        REAL(DP), INTENT(IN)      :: tau(3, natom)
        REAL(DP), INTENT(OUT)     :: sim_func(11, nele, nele)
        !
        !local variables
        !
        LOGICAL                  :: alive
        INTEGER(I4B)             :: nstru
        INTEGER(I4B)             :: i, j, i1, i2, i3, j1, j2, j3, k, n, l, m, ix, iy, iz
        INTEGER(I4B)             :: counter(nele, nele), super_ele_num(nele)
        REAL(DP)                 :: rms, x, y, z, r, theta, phi
        REAL(DP)                 :: mo(nele, nele)
        REAL(DP)                 :: pos(3, natom*500, nele), dis(nele, nele), tau_car(3, natom)
        REAL(DP)                 :: eps, eps1, tmp, temp(3)
        COMPLEX                  :: Yx(nele, nele), sumY(nele, nele)

        eps = 1e-1
        eps1 = 1e-4
        Yx = (0.0, 0.0)
        sumY = (0.0, 0.0)
        mo = 0.0
        pos = 0.0

        !WRITE(*,*) structyp, nele, natom, ele_num
        !DO i = 1, 3
        !  WRITE(*,*) (lat_matrix(i,j), j = 1, 3)
        !END DO
        !DO i = 1, natom
        !  WRITE(*,*) tau(:,i)
        !END DO

        IF (structyp == 2) THEN
            ix = 0
            iy = 0
            iz = 0
        ELSE IF (structyp == 1) THEN
            ix = 3
            iy = 3
            iz = 3
        END IF

        super_ele_num = 0
        k = 0
        DO i = 1, nele
            DO j = 1, ele_num(i)
                k = k + 1
                DO i1 = -ix, ix
                    DO i2 = -iy, iy
                        DO i3 = -iz, iz
                            super_ele_num(i) = super_ele_num(i) + 1

                            temp(1) = tau(1, k) + i1
                            temp(2) = tau(2, k) + i2
                            temp(3) = tau(3, k) + i3
                        pos(1, super_ele_num(i), i) = temp(1)*lat_matrix(1, 1) + temp(2)*lat_matrix(2, 1) + temp(3)*lat_matrix(3, 1)
                        pos(2, super_ele_num(i), i) = temp(1)*lat_matrix(1, 2) + temp(2)*lat_matrix(2, 2) + temp(3)*lat_matrix(3, 2)
                        pos(3, super_ele_num(i), i) = temp(1)*lat_matrix(1, 3) + temp(2)*lat_matrix(2, 3) + temp(3)*lat_matrix(3, 3)
                        END DO
                    END DO
                END DO
            END DO
        END DO
        k = 0
        DO i = 1, nele
            DO j = 1, ele_num(i)
                k = k + 1
                tau_car(1, k) = tau(1, k)*lat_matrix(1, 1) + tau(2, k)*lat_matrix(2, 1) + tau(3, k)*lat_matrix(3, 1)
                tau_car(2, k) = tau(1, k)*lat_matrix(1, 2) + tau(2, k)*lat_matrix(2, 2) + tau(3, k)*lat_matrix(3, 2)
                tau_car(3, k) = tau(1, k)*lat_matrix(1, 3) + tau(2, k)*lat_matrix(2, 3) + tau(3, k)*lat_matrix(3, 3)
            END DO
        END DO

        dis = 100.0
        DO l = 0, 10
            mo = 0.0
            DO m = -l, l
                sumY = (0.0, 0.0)
                counter = 0
                k = 0
                DO i1 = 1, nele
                    DO j1 = 1, ele_num(i1)
                        k = k + 1
                        DO i2 = 1, nele
                            DO j2 = 1, super_ele_num(i2)
                                x = tau_car(1, k) - pos(1, j2, i2)
                                y = tau_car(2, k) - pos(2, j2, i2)
                                z = tau_car(3, k) - pos(3, j2, i2)
                                call BondCharMatrix_Car2Sphe(x, y, z, r, theta, phi)
                                IF (r > eps .and. r < 3.5) THEN
                                    counter(i1, i2) = counter(i1, i2) + 1
                                    call BondCharMatrix_SphericalHarmonic(l, m, theta, phi, Yx(i1, i2))
                                    IF (structyp == 2) THEN
                                        sumY(i1, i2) = sumY(i1, i2) + Yx(i1, i2)
                                    ELSE
                                        !sumY(i1,i2) = sumY(i1,i2) + exp(-0.001*r)*Yx(i1,i2)
                                        sumY(i1, i2) = sumY(i1, i2) + Yx(i1, i2)
                                    END IF
                                END IF
                            END DO
                        END DO
                    END DO
                END DO
                DO i1 = 1, nele
                    DO i2 = 1, nele
                        IF (counter(i1, i2) == 0) THEN
                            sumY(i1, i2) = 0.0
                        ELSE
                            sumY(i1, i2) = sumY(i1, i2)/REAL(counter(i1, i2))
                        END IF
                        mo(i1, i2) = mo(i1, i2) + cabs(sumY(i1, i2))**2
                    END DO
                END DO
            END DO

            DO i1 = 1, nele
                DO i2 = 1, nele
                    sim_func(l + 1, i1, i2) = sqrt(mo(i1, i2)*4*pi/(2*l + 1))
                END DO
            END DO
        END DO

        DO l = 1, 10
            IF (mod(l, 2) == 0) sim_func(l, :, :) = 0
        END DO
        !
    END SUBROUTINE BondCharMatrix_CalBCM

    Subroutine BondCharMatrix_Car2Sphe(x, y, z, r, theta, psi)

        real(DP)      :: x, y, z, r, theta, psi, tol

        tol = 0.01

        r = sqrt(x*x + y*y + z*z)

        if (abs(x) < tol .and. abs(y) < tol) then
            if (abs(z) < tol) then
                theta = 0.0
            else if (z < 0.0) then
                theta = pi
            else
                theta = 0.0
            end if
        else
            theta = acos(z/r)
        end if
        if (abs(x) < tol .and. abs(y) < tol) then
            psi = 0.0
        else if (abs(x) < tol .and. y > 0.0) then
            psi = pi/2.0
        else if (abs(x) < tol .and. y < 0.0) then
            psi = 3*pi/2.0
        else if (abs(y) < tol .and. x > 0.0) then
            psi = 0.0
        else if (abs(y) < tol .and. x < 0.0) then
            psi = pi
        else if (x > 0 .and. y > 0) then
            psi = atan(y/x)
        else if (x < 0 .and. y > 0) then
            psi = atan(y/x) + pi
        else if (x < 0 .and. y < 0) then
            psi = atan(y/x) + pi
        else if (x > 0 .and. y < 0) then
            psi = atan(y/x) + 2*pi
        end if

    end subroutine BondCharMatrix_Car2Sphe

    subroutine BondCharMatrix_SphericalHarmonic(l, m, angle1, angle2, Y)
        ! angle1 is theta  0-pi
        ! angle2 is phi    0-2pi

        integer(i4b) l, m
        real(DP) angle1, angle2, sina, cosa, sinb, cosb, c, d, topi
        complex i, Y

        i = (0, 1)
        sina = sin(angle1)
        cosa = cos(angle1)
        topi = 2.0*pi
        if (l == 0) then
            if (m == 0) Y = 0.5*sqrt(1.0/pi)
        elseif (l == 1) then
            if (m == -1) then
                Y = 0.5*sqrt(3.0/(topi))*(exp(-i*angle2))*sin(angle1)
            elseif (m == 0) then
                Y = 0.5*sqrt(3.0/pi)*cos(angle1)
            elseif (m == 1) then
                Y = -0.5*sqrt(3.0/(topi))*exp(i*angle2)*sin(angle1)
            end if
        else if (l == 2) then
            if (m == -2) then
                Y = 0.25*sqrt(15.0/(topi))*exp(-i*2.0*angle2)*sin(angle1)**2
            elseif (m == -1) then
                Y = 0.5*sqrt(15.0/(topi))*exp(-i*angle2)*sin(angle1)*cos(angle1)
            elseif (m == 0) then
                Y = 0.25*sqrt(5.0/pi)*(3.0*(cos(angle1))**2 - 1.0)
            elseif (m == 1) then
                Y = -0.5*sqrt(15.0/(topi))*exp(i*angle2)*sin(angle1)*cos(angle1)
            elseif (m == 2) then
                Y = 0.25*sqrt(15.0/topi)*exp(i*2.0*angle2)*sina**2
            end if
        else if (l == 3) then
            if (m == -3) then
                Y = 1.0/8.0*sqrt(35.0/pi)*exp(-i*3.0*angle2)*(sin(angle1))**3
            elseif (m == -2) then
                Y = 1.0/4.0*sqrt(105.0/topi)*exp(-i*2.0*angle2)*(sin(angle1))**2*cos(angle1)
            elseif (m == -1) then
                Y = 1.0/8.0*sqrt(21.0/pi)*exp(-i*angle2)*sin(angle1)*(5*(cos(angle1))**2 - 1)
            elseif (m == 0) then
                Y = 1.0/4.0*sqrt(7.0/pi)*(5.0*(cos(angle1))**3 - 3.0*cos(angle1))
            elseif (m == 1) then
                Y = -1.0/8.0*sqrt(21.0/pi)*exp(i*angle2)*sin(angle1)*(5.0*(cos(angle1))**2 - 1)
            elseif (m == 2) then
                Y = 1.0/4.0*sqrt(105.0/(topi))*exp(2.0*i*angle2)*(sin(angle1))**2*cos(angle1)
            elseif (m == 3) then
                Y = -1.0/8.0*sqrt(35.0/pi)*exp(3*i*angle2)*(sin(angle1))**3
            end if
        else if (l == 4) then
            if (m == -4) then
                Y = 3.0/16.0*sqrt(35.0/topi)*exp(-4.0*i*angle2)*(sin(angle1))**4
            elseif (m == -3) then
                Y = 3.0/8.0*sqrt(35.0/pi)*exp(-3*i*angle2)*(sin(angle1))**3*cos(angle1)
            elseif (m == -2) then
                Y = 3.0/8.0*sqrt(5.0/topi)*exp(-2.0*i*angle2)*(sin(angle1))**2*(7.0*(cos(angle1))**2 - 1)
            elseif (m == -1) then
                Y = 3.0/8.0*sqrt(5.0/pi)*exp(-i*angle2)*sin(angle1)*(7.0*(cos(angle1))**3 - 3.0*cos(angle1))
            elseif (m == 0) then
                Y = 3.0/16.0*sqrt(1.0/pi)*(35.0*cosa**4 - 30.0*cosa**2 + 3.0)
            elseif (m == 1) then
                Y = -3.0/8.0*sqrt(5.0/pi)*exp(i*angle2)*sina*(7.0*cosa**3 - 3.0*cosa)
            elseif (m == 2) then
                Y = 3.0/8.0*sqrt(5.0/topi)*exp(2.0*i*angle2)*sina**2*(7.0*cosa**2 - 1.0)
            elseif (m == 3) then
                Y = -3.0/8.0*sqrt(35.0/pi)*exp(3.0*i*angle2)*sina**3*cosa
            elseif (m == 4) then
                Y = 3.0/16.0*sqrt(35.0/topi)*exp(4.0*i*angle2)*sina**4
            end if
        elseif (l == 5) then
            if (m == -5) then
                Y = 3.0/32.0*sqrt(77.0/pi)*exp(-5.0*i*angle2)*sina**5
            elseif (m == -4) then
                Y = 3.0/16.0*sqrt(385.0/topi)*exp(-4.0*i*angle2)*sina**4*cosa
            elseif (m == -3) then
                Y = 1.0/32.0*sqrt(385.0/pi)*exp(-3.0*i*angle2)*sina**3*(9.0*cosa**2 - 1.0)
            elseif (m == -2) then
                Y = 1.0/8.0*sqrt(1155.0/topi)*exp(-2.0*i*angle2)*sina**2*(3.0*cosa**3 - cosa)
            elseif (m == -1) then
                Y = 1.0/16.0*sqrt(165.0/topi)*exp(-i*angle2)*sina*(21.0*cosa**4 - 14.0*cosa**2 + 1.0)
            elseif (m == 0) then
                Y = 1.0/16.0*sqrt(11.0/pi)*(63.0*cosa**5 - 70.0*cosa**3 + 15.0*cosa)
            elseif (m == 1) then
                Y = -1.0/16.0*sqrt(165.0/topi)*exp(i*angle2)*sina*(21.0*cosa**4 - 14.0*cosa**2 + 1.0)
            elseif (m == 2) then
                Y = 1.0/8.0*sqrt(1155.0/topi)*exp(i*angle2)*sina**2*(3.0*cosa**3 - cosa)
            elseif (m == 3) then
                Y = -1.0/32.0*sqrt(385.0/pi)*exp(i*3.0*angle2)*sina**3*(9.0*cosa**2 - 1.0)
            elseif (m == 4) then
                Y = 3.0/16.0*sqrt(385.0/topi)*exp(i*4.0*angle2)*sina**4*cosa
            elseif (m == 5) then
                Y = -3.0/32.0*sqrt(77.0/pi)*exp(5.0*i*angle2)*sina**5
            end if
        elseif (l == 6) then
            if (m == -6) then
                Y = 1.0/64.0*sqrt(3003.0/pi)*exp(-6.0*i*angle2)*sina**6
            elseif (m == -5) then
                Y = 3.0/32.0*sqrt(1001.0/pi)*exp(-5.0*i*angle2)*sina**5*cosa
            elseif (m == -4) then
                Y = 3.0/32.0*sqrt(91.0/topi)*exp(-4.0*i*angle2)*sina**4*(11.0*cosa**2 - 1.0)
            elseif (m == -3) then
                Y = 1.0/32.0*sqrt(1365.0/pi)*exp(-3.0*i*angle2)*sina**3*(11.0*cosa**3 - 3.0*cosa)
            elseif (m == -2) then
                Y = 1.0/64.0*sqrt(1365.0/pi)*exp(-2.0*i*angle2)*sina**2*(33.0*cosa**4 - 18.0*cosa**2 + 1.0)
            elseif (m == -1) then
                Y = 1.0/16.0*sqrt(273.0/topi)*exp(-i*angle2)*sina*(33.0*cosa**5 - 30.0*cosa**3 + 5.0*cosa)
            elseif (m == 0) then
                Y = 1.0/32.0*sqrt(13.0/pi)*(231.0*cosa**6 - 315.0*cosa**4 + 105.0*cosa**2 - 5.0)
            elseif (m == 1) then
                Y = -1.0/16.0*sqrt(273.0/topi)*exp(i*angle2)*sina*(33.0*cosa**5 - 30.0*cosa**3 + 5.0*cosa)
            elseif (m == 2) then
                Y = 1.0/64.0*sqrt(1365.0/pi)*exp(i*2.0*angle2)*sina**2*(33.0*cosa**4 - 18.0*cosa**2 + 1.0)
            elseif (m == 3) then
                Y = -1.0/32.0*sqrt(1365.0/pi)*exp(3.0*i*angle2)*sina**3*(11.0*cosa**3 - 3.0*cosa)
            elseif (m == 4) then
                Y = 3.0/32.0*sqrt(91.0/topi)*exp(4.0*i*angle2)*sina**4*(11.0*cosa**2 - 1.0)
            elseif (m == 5) then
                Y = -3.0/32.0*sqrt(1001.0/pi)*exp(5.0*i*angle2)*sina**5*cosa
            elseif (m == 6) then
                Y = 1.0/64.0*sqrt(3003.0/pi)*exp(6.0*i*angle2)*sina**6
            end if
        elseif (l == 7) then
            if (m == -7) then
                Y = 3.0/64.0*sqrt(715.0/topi)*exp(-7.0*i*angle2)*sina**7
            elseif (m == -6) then
                Y = 3.0/64.0*sqrt(5005.0/pi)*exp(-6.0*i*angle2)*sina**6*cosa
            elseif (m == -5) then
                Y = 3.0/64.0*sqrt(385.0/topi)*exp(-5.0*i*angle2)*sina**5*(13.0*cosa**2 - 1.0)
            elseif (m == -4) then
                Y = 3.0/32.0*sqrt(385.0/topi)*exp(-4.0*i*angle2)*sina**4*(13.0*cosa**3 - 3.0*cosa)
            elseif (m == -3) then
                Y = 3.0/64.0*sqrt(35.0/topi)*exp(-3.0*i*angle2)*sina**3*(143.0*cosa**4 - 66.0*cosa**2 + 3.0)
            elseif (m == -2) then
                Y = 3.0/64.0*sqrt(35.0/pi)*exp(-2.0*i*angle2)*sina**2*(143.0*cosa**5 - 110.0*cosa**3 + 15.0*cosa)
            elseif (m == -1) then
                Y = 1.0/64.0*sqrt(105.0/topi)*exp(-i*angle2)*sina*(429.0*cosa**6 - 495.0*cosa**4 + 135.0*cosa**2 - 5.0)
            elseif (m == 0) then
                Y = 1.0/32.0*sqrt(15.0/pi)*(429.0*cosa**7 - 693.0*cosa**5 + 315.0*cosa**3 - 35.0*cosa)
            elseif (m == 1) then
                Y = -1.0/64.0*sqrt(105.0/topi)*exp(i*angle2)*sina*(429.0*cosa**6 - 495.0*cosa**4 + 135.0*cosa**2 - 5.0)
            elseif (m == 2) then
                Y = 3.0/64.0*sqrt(35.0/pi)*exp(2.0*i*angle2)*sina**2*(143.0*cosa**5 - 110.0*cosa**3 + 15.0*cosa)
            elseif (m == 3) then
                Y = -3.0/64.0*sqrt(35.0/topi)*exp(3.0*i*angle2)*sina**3*(143.0*cosa**4 - 66.0*cosa**2 + 3.0)
            elseif (m == 4) then
                Y = 3.0/32.0*sqrt(385.0/topi)*exp(4.0*i*angle2)*sina**4*(13.0*cosa**3 - 3.0*cosa)
            elseif (m == 5) then
                Y = -3.0/64.0*sqrt(385.0/topi)*exp(5.0*i*angle2)*sina**5*(13.0*cosa**2 - 1.0)
            elseif (m == 6) then
                Y = 3.0/64.0*sqrt(5005.0/pi)*exp(6.0*i*angle2)*sina**6*cosa
            elseif (m == 7) then
                Y = -3.0/64.0*sqrt(715.0/topi)*exp(7.0*i*angle2)*sina**7
            end if
        elseif (l == 8) then
            if (m == -8) then
                Y = 3.0/256.0*sqrt(12155.0/topi)*exp(-8.0*i*angle2)*sina**8
            elseif (m == -7) then
                Y = 3.0/64.0*sqrt(12155.0/topi)*exp(-7.0*i*angle2)*sina**7*cosa
            elseif (m == -6) then
                Y = 1.0/128.0*sqrt(7293.0/pi)*exp(-6.0*i*angle2)*sina**6*(15.0*cosa**2 - 1.0)
            elseif (m == -5) then
                Y = 3.0/64.0*sqrt(17017.0/topi)*exp(-5.0*i*angle2)*sina**5*(5.0*cosa**3 - cosa)
            elseif (m == -4) then
                Y = 3.0/128.0*sqrt(1309.0/topi)*exp(-4.0*i*angle2)*sina**4*(65.0*cosa**4 - 26.0*cosa**2 + 1.0)
            elseif (m == -3) then
                Y = 1.0/64.0*sqrt(19635.0/topi)*exp(-3.0*i*angle2)*sina**3*(39.0*cosa**5 - 26.0*cosa**3 + 3.0*cosa)
            elseif (m == -2) then
                Y = 3.0/128.0*sqrt(595.0/pi)*exp(-2.0*i*angle2)*sina**2*(143.0*cosa**6 - 143.0*cosa**4 + 33.0*cosa**2 - 1.0)
            elseif (m == -1) then
                Y = 3.0/64.0*sqrt(17.0/topi)*exp(-i*angle2)*sina*(715.0*cosa**7 - 1001.0*cosa**5 + 385.0*cosa**3 - 35.0*cosa)
            elseif (m == 0) then
                Y = 1.0/256.0*sqrt(17.0/pi)*(6435.0*cosa**8 - 12012.0*cosa**6 + 6930.0*cosa**4 - 1260.0*cosa**2 + 35.0)
            elseif (m == 1) then
                Y = -3.0/64.0*sqrt(17.0/topi)*exp(i*angle2)*sina*(715.0*cosa**7 - 1001.0*cosa**5 + 385.0*cosa**3 - 35.0*cosa)
            elseif (m == 2) then
                Y = 3.0/128.0*sqrt(595.0/pi)*exp(2.0*i*angle2)*sina**2*(143.0*cosa**6 - 143.0*cosa**4 + 33.0*cosa**2 - 1.0)
            elseif (m == 3) then
                Y = -1.0/64.0*sqrt(19635.0/topi)*exp(3.0*i*angle2)*sina**3*(39.0*cosa**5 - 26.0*cosa**3 + 3.0*cosa)
            elseif (m == 4) then
                Y = 3.0/128.0*sqrt(1309.0/topi)*exp(4.0*i*angle2)*sina**4*(65.0*cosa**4 - 26.0*cosa**2 + 1.0)
            elseif (m == 5) then
                Y = -3.0/64.0*sqrt(17017.0/topi)*exp(5.0*i*angle2)*sina**5*(5.0*cosa**3 - cosa)
            elseif (m == 6) then
                Y = 1.0/128.0*sqrt(7293.0/pi)*exp(6.0*i*angle2)*sina**6*(15.0*cosa**2 - 1.0)
            elseif (m == 7) then
                Y = -3.0/64.0*sqrt(12155.0/topi)*exp(7.0*i*angle2)*sina**7*cosa
            elseif (m == 8) then
                Y = 3.0/256.0*sqrt(12155.0/topi)*exp(8.0*i*angle2)*sina**8
            end if
        elseif (l == 9) then
            if (m == -9) then
                Y = 1.0/512.0*sqrt(230945.0/pi)*exp(-9.0*i*angle2)*sina**9
            elseif (m == -4) then
                Y = 3.0/128.0*sqrt(95095.0/topi)*exp(-4.0*i*angle2)*sina**4*(17.0*cosa**5 - 10.0*cosa**3 + cosa)
            elseif (m == -5) then
                Y = 3.0/256.0*sqrt(2717.0/pi)*exp(-5.0*i*angle2)*sina**5*(85.0*cosa**4 - 30.0*cosa**2 + 1.0)
            elseif (m == -6) then
                Y = 1.0/128.0*sqrt(40755.0/pi)*exp(-6.0*i*angle2)*sina**6*(17.0*cosa**3 - 3.0*cosa)
            elseif (m == -7) then
                Y = 3.0/512.0*sqrt(13585.0/pi)*exp(-7.0*i*angle2)*sina**7*(17.0*cosa**2 - 1.0)
            elseif (m == -8) then
                Y = 3.0/256.0*sqrt(230945.0/topi)*exp(-8.0*i*angle2)*sina**8*cosa
            elseif (m == -3) then
                Y = 1.0/256.0*sqrt(21945.0/pi)*exp(-3.0*i*angle2)*sina**3*(221.0*cosa**6 - 195.0*cosa**4 + 39.0*cosa**2 - 1.0)
            elseif (m == -2) then
                Y = 3.0/128.0*sqrt(1045.0/pi)*exp(-2.0*i*angle2)*sina**2*(221.0*cosa**7 - 273.0*cosa**5 + 91.0*cosa**3 - 7.0*cosa)
            elseif (m == -1) then
          Y = 3.0/256.0*sqrt(95.0/topi)*exp(-i*angle2)*sina*(2431.0*cosa**8 - 4004.0*cosa**6 + 2002.0*cosa**4 - 308.0*cosa**2 + 7.0)
            elseif (m == 0) then
                Y = 1.0/256.0*sqrt(19.0/pi)*(12155.0*cosa**9 - 25740.0*cosa**7 + 18018.0*cosa**5 - 4620.0*cosa**3 + 315.0*cosa)
            elseif (m == 1) then
          Y = -3.0/256.0*sqrt(95.0/topi)*exp(i*angle2)*sina*(2431.0*cosa**8 - 4004.0*cosa**6 + 2002.0*cosa**4 - 308.0*cosa**2 + 7.0)
            elseif (m == 2) then
                Y = 3.0/128.0*sqrt(1045.0/pi)*exp(2.0*i*angle2)*sina**2*(221.0*cosa**7 - 273.0*cosa**5 + 91.0*cosa**3 - 7.0*cosa)
            elseif (m == 3) then
                Y = -1.0/256.0*sqrt(21945.0/pi)*exp(3.0*i*angle2)*sina**3*(221.0*cosa**6 - 195.0*cosa**4 + 39.0*cosa**2 - 1.0)
            elseif (m == 4) then
                Y = 3.0/128.0*sqrt(95095.0/topi)*exp(4.0*i*angle2)*sina**4*(17.0*cosa**5 - 10*cosa**3 + cosa)
            elseif (m == 5) then
                Y = -3.0/256.0*sqrt(2717.0/pi)*exp(5.0*i*angle2)*sina**5*(85.0*cosa**4 - 30.0*cosa**2 + 1.0)
            elseif (m == 6) then
                Y = 1.0/128.0*sqrt(40755.0/pi)*exp(6.0*i*angle2)*sina**6*(17.0*cosa**3 - 3.0*cosa)
            elseif (m == 7) then
                Y = -3.0/512.0*sqrt(13585.0/pi)*exp(7.0*i*angle2)*sina**7*(17.0*cosa**2 - 1.0)
            elseif (m == 8) then
                Y = 3.0/256.0*sqrt(230945.0/topi)*exp(8.0*i*angle2)*sina**8.0*cosa
            elseif (m == 9) then
                Y = -1.0/512.0*sqrt(230945.0/pi)*exp(9.0*i*angle2)*sina**9
            end if
        elseif (l == 10) then
            if (m == -10) then
                Y = 1.0/1024.0*sqrt(969969.0/pi)*exp(-10.0*i*angle2)*sina**10
            elseif (m == -9) then
                Y = 1.0/512.0*sqrt(4849845.0/pi)*exp(-9.0*i*angle2)*sina**9*cosa
            elseif (m == -8) then
                Y = 1.0/512.0*sqrt(255255.0/topi)*exp(-8.0*i*angle2)*sina**8*(19.0*cosa**2 - 1.0)
            elseif (m == -7) then
                Y = 3.0/512.0*sqrt(85085.0/pi)*exp(-7.0*i*angle2)*sina**7*(19.0*cosa**3 - 3.0*cosa)
            elseif (m == -6) then
                Y = 3.0/1024.0*sqrt(5005.0/pi)*exp(-6.0*i*angle2)*sina**6*(323.0*cosa**4 - 102.0*cosa**2 + 3.0)
            elseif (m == -5) then
                Y = 3.0/256.0*sqrt(1001.0/pi)*exp(-5.0*i*angle2)*sina**5*(323.0*cosa**5 - 170.0*cosa**3 + 15*cosa)
            elseif (m == -4) then
                Y = 3.0/256.0*sqrt(5005.0/topi)*exp(-4.0*i*angle2)*sina**4*(323.0*cosa**6 - 255.0*cosa**4 + 45.0*cosa**2 - 1.0)
            elseif (m == -3) then
                Y = 3.0/256.0*sqrt(5005.0/pi)*exp(-3.0*i*angle2)*sina**3*(323.0*cosa**7 - 357.0*cosa**5 + 105.0*cosa**3 - 7.0*cosa)
            elseif (m == -2) then
  Y = 3.0/512.0*sqrt(385.0/topi)*exp(-2.0*i*angle2)*sina**2*(4199.0*cosa**8 - 6188.0*cosa**6 + 2730.0*cosa**4 - 364.0*cosa**2 + 7.0)
            elseif (m == -1) then
 Y = 1.0/256.0*sqrt(1155.0/topi)*exp(-i*angle2)*sina*(4199.0*cosa**9 - 7956.0*cosa**7 + 4914.0*cosa**5 - 1092.0*cosa**3 + 63.0*cosa)
            elseif (m == 0) then
       Y = 1.0/512.0*sqrt(21.0/pi)*(46189.0*cosa**10 - 109395.0*cosa**8 + 90090.0*cosa**6 - 30030.0*cosa**4 + 3465.0*cosa**2 - 63.0)
            elseif (m == 1) then
 Y = -1.0/256.0*sqrt(1155.0/topi)*exp(i*angle2)*sina*(4199.0*cosa**9 - 7956.0*cosa**7 + 4914.0*cosa**5 - 1092.0*cosa**3 + 63.0*cosa)
            elseif (m == 2) then
   Y = 3.0/512.0*sqrt(385.0/topi)*exp(2.0*i*angle2)*sina**2*(4199.0*cosa**8 - 6188.0*cosa**6 + 2730.0*cosa**4 - 364.0*cosa**2 + 7.0)
            elseif (m == 3) then
                Y = -3.0/256.0*sqrt(5005.0/pi)*exp(3.0*i*angle2)*sina**3*(323.0*cosa**7 - 357.0*cosa**5 + 105.0*cosa**3 - 7.0*cosa)
            elseif (m == 4) then
                Y = 3.0/256.0*sqrt(5005.0/topi)*exp(4.0*i*angle2)*sina**4*(323.0*cosa**6 - 255.0*cosa**4 + 45.0*cosa**2 - 1.0)
            elseif (m == 5) then
                Y = -3.0/256.0*sqrt(1001.0/pi)*exp(5.0*i*angle2)*sina**5*(323.0*cosa**5 - 170.0*cosa**3 + 15.0*cosa)
            elseif (m == 6) then
                Y = 3.0/1024.0*sqrt(5005.0/pi)*exp(6.0*i*angle2)*sina**6*(323.0*cosa**4 - 102.0*cosa**2 + 3.0)
            elseif (m == 7) then
                Y = -3.0/512.0*sqrt(85085.0/pi)*exp(7.0*i*angle2)*sina**7*(19.0*cosa**3 - 3.0*cosa)
            elseif (m == 8) then
                Y = 1.0/512.0*sqrt(255255.0/topi)*exp(8.0*i*angle2)*sina**8*(19.0*cosa**2 - 1.0)
            elseif (m == 9) then
                Y = -1.0/512.0*sqrt(4849845.0/pi)*exp(9.0*i*angle2)*sina**9*cosa
            elseif (m == 10) then
                Y = 1.0/1024.0*sqrt(969969.0/pi)*exp(10.0*i*angle2)*sina**10
            end if
        end if
    END SUBROUTINE BondCharMatrix_SphericalHarmonic

END MODULE BondCharMatrix

FUNCTION BondCharMatrix_CalDis(nele, sim_func1, sim_func2)

    USE KINDS, ONLY: DP, I4B

    INTEGER(I4B), INTENT(IN)      :: nele
    REAL(DP), INTENT(IN)      :: sim_func1(11, nele, nele)
    REAL(DP), INTENT(IN)      :: sim_func2(11, nele, nele)

    ! local variables

    INTEGER(I4B)                  :: i, j, k
    REAL(DP)  :: BondCharMatrix_CalDis

    dis = 0.0

    DO i = 1, nele
        DO j = 1, nele
            DO k = 0, 10, 2
                dis = dis + (sim_func1(k + 1, i, j) - sim_func2(k + 1, i, j))**2
            END DO
        END DO
    END DO

    BondCharMatrix_CalDis = sqrt(dis)

END FUNCTION
