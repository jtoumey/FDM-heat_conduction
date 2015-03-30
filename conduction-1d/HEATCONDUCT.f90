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
integer n,ii,jj,kk
parameter (n=81) ! number of grid points
integer n_iter
real dx,L,x(n)
real a,b,Q,h
double precision T_0,T_L,T_inf,T(n),dT_ic,T_old(n)
double precision f(n-1),J_a(n-1),J_b(n-1),J_c(n-1),phi(n-1)
real tol
real c1,c2,c3,c4,c5,c6
real T1,T2,T3,T4,T5
real k
!
!...INPUT SECTION
!
!WRITE *,'ENTER NUMBER OF GRID NODES: '
!READ(*,*)n
!
!...System parameters
!   Temperatures
T_0    = 350. ! [K] temp at the left wall (given)
T_L    = 500. ! [K] temp at the right wall (estimate)
T_inf  = 500. ! [K] ambient temp (given)
n_iter = 0
tol    = 1.
!   physical parameters
a = 0.01  ! [W/m-K]
b = 0.001 ! [W/m-K^2]
Q = 1.e6  ! [W/m^3]
h = 50.   ! [W/K-m^2]
!
!...Grid setup
!
L     = 0.2                  ! [m] Length of domain
dx    = L/float(n-1)         ! [m]
dT_ic = (T_L-T_0)/float(n-1) ! [K] temp increment for IC
!
!...grid vector
!
do ii = 1,n
   x(ii) = (ii-1)*dx
   T(ii) = T_0+(ii-1)*dT_ic ! temperature IC vector
end do
!
!...Begin Newton iteration
!
do while (tol .gt. 1.e-6)
   ! save IC to check convergence
   !
   T_old = T;
   !
   !...Construct solution vector f
   !   interior points
   do ii = 2,n-1
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
   !   
   !...Solve for dT: J(T_0)*dT = -F(T_0)  
   !
   call thomas(n-1,J_a,J_b,J_c,-f,phi)
   !
   !...Update T
   !
   T(1) = T(1) + 0. ! Because the left value does not change (dirichlet)
   do jj = 1,n-1
      T(jj+1) = T(jj+1) + phi(jj)
   end do
   !
   !   Test convergence
   !
   tol = maxval(T - T_old)
   n_iter = n_iter + 1
   !
   ! write current cycle, error, temperature at right bndry
   write(6,401)n_iter,tol,T(n)
end do
!
!...call subroutine to write file
!
call write_results(n,x,T)
!
!...format statements
!
401 format(3x,'***  Iteration : ',i8,3x,'Residual : ',f14.7,3x,'T_R = ',f14.7,'  ***')
!
501 format(3x,f20.7,3x,f20.7,3x,f20.7)
!
END
