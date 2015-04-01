!*************************************************************************!
!                                                                         !
!  Subroutine:   CALC_JACOBIAN.F90                                        !
!                                                                         !
!  Programmer:   Julian M. Toumey                                         !
!                Madison, WI                                              !
!                                                                         !
!  Date:         March 2015                                               !
!                                                                         !
!  Language:     FORTRAN90                                                !
!                                                                         !
!  Description:  Calculates the Jacobian of a system of finite difference !
!                equations that approximated heat conduction in a 1-D     !
!                bar. The Jacobian in this case is a tri-diagonal matrix, !
!                but this subroutine returns three column vectors         !
!                containing the sub diagonal, main diagonal, and sup      !
!                diagonal. 
!                                                                         !
!                  J_a   Sub diagonal                                     !
!                  J_b   Main diagonal                                    !
!                  J_c   Super diagonal                                   !
!                  d     Right side of Linear System                      !
!                  phi   Solution vector                                  !
!                                                                         !
!*************************************************************************!
SUBROUTINE CALC_JACOBIAN(n,dx,a,b,h,T_inf,T(n),J_a,J_b,J_c)
IMPLICIT NONE
!
integer n,ii ! n from main, ii used only in this fcn
real dx ! from main
real a,b,h ! from main
double precision T_inf,T(n)  ! from main
double precision J_a(n-1),J_b(n-1),J_c(n-1)  ! the subroutine returns these three
real c6 ! used only in this fcn
real T1,T2,T3,T4,T5 ! used only in this fcn
!
!...Construct Jacobian
!   Left boundary, Dirichlet b.c.
ii = 1
J_a(ii) = 0.
J_b(ii) = b*(T(ii+2)   - T(ii))**2/(2.*dx**2) - 2.*(b*T(ii+1)**2 &
   + a)/dx**2 + (2.*T(ii+1)*b*(T(ii) - 2.*T(ii+1) + T(ii+2)))/dx**2
J_c(ii) = (b* T(ii+1)**2 + a)/dx**2 - T(ii+1)*b*2.*(T(ii+2) &
   - T(ii))/(2.*dx**2)
!   
! Interior points
do ii = 2,n-2
J_a(ii) = (b* T(ii+1)**2 + a)/dx**2 - T(ii+1)*b*2.*(T(ii+2) &
      - T(ii))/(2.*dx**2)
J_b(ii) =  b*(T(ii+2)   - T(ii))**2/(2.*dx**2) &
      -  2.*(b*T(ii+1)**2 + a)/dx**2 + (2.*T(ii+1)*b*(T(ii) &
      -  2.*T(ii+1) + T(ii+2)))/dx**2
J_c(ii) = (b* T(ii+1)**2 + a)/dx**2 - T(ii+1)*b*2.*(T(ii+2) &
      - T(ii))/(2.*dx**2)
end do
!
!   Right boundary; Robin bc
ii = n-1
J_a(ii) = 2.*(b*T(ii+1)**2 + a)/(dx**2)
!
!   Coefficients of the derivative
C6 = b*T(ii+1)**2 + a
T1 = 2.*b*h**2*(T(ii+1) - T_inf)**2/C6**2
T2 = 2.*T(ii+1)*b*(2.*T(ii+1) - 2.*T(ii) + 2.*dx*h*(T(ii+1) &
- T_inf)/C6)/dx**2
T3 = C6*(2.*dx*h/C6 - (4.*T(ii+1)*b*dx*h*(T(ii+1) - T_inf))/C6**2 &
+ 2.)/dx**2
T4 = (2.*T(ii+1)*b*h**2*(2.*T(ii+1) - 2.*T_inf))/C6**2
T5 = 8.*T(ii+1)**2*b**2*h**2*(T(ii+1) - T_inf)**2/C6**3
J_b(ii) = T1 - T2 - T3 + T4 - T5
J_c(ii) = 0.



END SUBROUTINE CALC_JACOBIAN