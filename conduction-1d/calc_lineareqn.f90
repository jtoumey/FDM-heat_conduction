!*************************************************************************!
!                                                                         !
!  Subroutine:   CALC_LINEAREQN.F90                                       !
!                                                                         !
!  Programmer:   Julian M. Toumey                                         !
!                Madison, WI                                              !
!                                                                         !
!  Date:         April 2015                                               !
!                                                                         !
!  Language:     FORTRAN90                                                !
!                                                                         !
!  Description:  Calculates the linear system of discrete equations that  !
!                approximate the governing equation for heat conduction   !
!                in a 1-D bar for this problem.                           !
!                                                                         !
!                Inputs:                                                  !
!                  n     Number of grid pts                               !
!                  dx    grid spacing [m]                                 !
!                  a     physical parameter [W/m-K]                       !
!                  b     physical parameter [W/m-K^2]                     !
!                  h     physical parameter [W/K-m^2]                     !
!                  Q     physical parameter, source [W/m^3]               !
!                  T_inf ambient temp [K]                                 !
!                  T     size n; current temperature distribution [K]     !
!                                                                         !
!                Outputs:                                                 !
!                  f     size n-1; linear discrete equations              !
!                                                                         !
!*************************************************************************!
SUBROUTINE CALC_LINEAREQN(n,dx,a,b,h,Q,T_inf,T,f)
IMPLICIT NONE
!
integer n,ii
real dx
real a,b,Q,h
double precision T_inf,T(n)
double precision f(n-1)
real c1,c2,c3,c4,c5
real T1,T2,T3,T4,T5
real k
!
!...Construct solution vector f
!   interior points
do ii = 2,n-1
   !
   !...calculate the source term for MMS
   !
   call calc_source(ii,dx)
   !
   ! coefficients of f: c1, c2, c3, c4
   c1 = a + b*T(ii)**2
   c2 = (T(ii+1) - 2.*T(ii) + T(ii-1))/dx**2
   c3 = 2.*b*T(ii)
   c4 = ((T(ii+1) - T(ii-1))/(2.*dx))**2
   f(ii-1) = c1*c2 + c3*c4 + Q
end do
!   Right boundary, Robin b.c.
!
ii = n
k  = a + b*T(ii)**2
c5 = -2.*dx*h*(T(ii) - T_inf)/k + T(ii-1) 
f(ii-1) = (a + b*T(ii)**2)*( c5 - 2.*T(ii  ) +T(ii-1))/(dx**2) &
   + (2. * b*T(ii)   )*((c5 -   T(ii-1))/(2.*dx))**2 + Q
!
END SUBROUTINE CALC_LINEAREQN

SUBROUTINE CALC_SOURCE(ii,dx)
IMPLICIT NONE
!
integer ii
real dx
!
Q = -(     )



END SUBROUTINE CALC_SOURCE
