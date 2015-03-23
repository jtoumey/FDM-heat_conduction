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
real dx,L,x(n)
real a,b,Q,h
double precision T_0,T_L,T_inf,T_IC
!
!...INPUT SECTION
!
!WRITE *,'ENTER NUMBER OF GRID NODES: '
!READ(*,*) N
!
!...System parameters
!   Temperatures
T_0   = 350. ! [K] temp at the left wall (given)
T_L   = 500. ! [K] temp at the right wall (estimate)
T_inf = 500. ! [K] ambient temp (given)
!   physical parameters
a = 0.01  ! [W/m-K]
b = 0.001 ! [W/m-K^2]
Q = 1.e6  ! [W/m^3]
h = 50.   ! [W/K-m^2]
!
!...Grid setup
!
L  = 0.2          ! [m] Length of domain
dx = L/float(n-1) ! [m] grid spacing

do ii = 1,n+1
   x(ii) = (ii-1)*dx
   write(*,*)x(ii)
end do


END  
