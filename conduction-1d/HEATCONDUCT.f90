!*************************************************************************!
!                                                                         !
!  File:         HEATCONDUCT.f                                            !
!                                                                         !
!  Programmer:   Julian M. Toumey                                         !
!                Madison, WI                                              !
!                                                                         !
!  Date:         March 2015                                               !
!                                                                         !
!  Language:     FORTRAN 90                                               !
!                                                                         !
!  Description:  This module solves a nonlinear boundary value problem    !
!                that describes heat conduction in a 1-D bar. Uses        !
!                the finite difference method with second-order accurate  !
!                discretization in space. The non-linear system is        !
!                solved by Newton-Raphson iteration, and the linear       !
!                system is solved by the Thomas Algorithm.                !
!                                                                         !
! Output Files:  TEMP  Two-column output: x-location vs. Temperature      !
!                                                                         !
!                The algorithms are taken from the GNU Scientific         !
!                Library (q.v.).                                          !
!                                                                         !
!*************************************************************************!
PROGRAM HEATCONDUCT
IMPLICIT NONE
!
integer n,ii
parameter (n=81) ! number of grid points
integer n_iter
real dx,L,x(n)
real a,b,Q,h
double precision T_0,T_L,T_inf,T(n),dT,T_old(n)
double precision f(n-1),J(n-1,n-1)
real tol
real c1,c2,c3,c4,c5
real k
!
!...INPUT SECTION
!
!WRITE *,'ENTER NUMBER OF GRID NODES: '
!READ(*,*) N
!
!...System parameters
!   Temperatures
T_0    = 350. ! [K] temp at the left wall (given)
T_L    = 500. ! [K] temp at the right wall (estimate)
T_inf  = 500. ! [K] ambient temp (given)
n_iter = 1
tol    = 1.
!   physical parameters
a = 0.01  ! [W/m-K]
b = 0.001 ! [W/m-K^2]
Q = 1.e6  ! [W/m^3]
h = 50.   ! [W/K-m^2]
!
!...Grid setup
!
L  = 0.2                  ! [m] Length of domain
dx = L/float(n-1)         ! [m]
dT = (T_L-T_0)/float(n-1) ! [K] temp increment for IC
!
!...grid vector
!
do ii = 1,n
   x(ii) = (ii-1)*dx
   T(ii) = T_0+(ii-1)*dT
end do
!
!...Begin Newton iteration
!
do while (tol .gt. 1.e-3)
   ! save IC to check convergence
   T_old = T; 
   ! write current cycle, error, temperature at right bndry
   write(6,401)n_iter,tol,T(n)
   !
   !...Construct solution vector f
   !   interior points
   do ii = 2,n-1
      ! interior pts
      c1 = a + b*T(ii)**2
      c2 = (T(ii+1) - 2*T(ii) + T(ii-1))/dx**2
      c3 = 2*b*T(ii)
      c4 = ((T(ii+1) - T(ii-1))/(2*dx))**2
      f(ii-1) = c1*c2 + c3*c4 + Q
   end do
   !   Right boundary, Robin b.c.
   ii = n
   k  = a + b*T(ii)**2
   c5 = -2*dx*h*(T(ii) - T_inf)/k + T(ii-1) 
   f(ii-1) = (a + b*T(ii)**2)*( c5 - 2*T(ii  ) +T(ii-1))/(dx**2) &
           + (2 * b*T(ii)   )*((c5 -   T(ii-1))/(2*dx))**2 + Q
   !...Construct Jacobian
   !   Left boundary, Dirichlet b.c.
   ii = 1
   J(ii,1) = 0
   !   Test convergence
   tol = maxval(T - T_old)
   n_iter = n_iter + 1
end do
!
401 format(3x,'***  Iteration : ',i8,3x,'Residual : ',f14.7,3x,'T_L = ',f14.7,'  ***')
!
END  
