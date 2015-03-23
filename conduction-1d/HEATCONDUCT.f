!*************************************************************************!
!                                                                         !
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
INTEGER N
PARAMETER (N=81)
DOUBLE PRECISION T_0,T_IC
DOUBLE PRECISION a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11
!
!...INPUT SECTION
!
WRITE *,'ENTER NUMBER OF GRID NODES: '
READ(*,*) N


dx   = 1.0d0/dfloat(N-1)
!dt   = 0.90d0*Mach*dx/dsqrt(2.d0)



END  
