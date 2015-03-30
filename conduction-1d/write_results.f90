!*************************************************************************!
!                                                                         !
!  Subroutine:   WRITE_RESULTS.F90                                        !
!                                                                         !
!  Programmer:   Julian M. Toumey                                         !
!                Madison, WI                                              !
!                                                                         !
!  Date:         29 March 2015                                            !
!                                                                         !
!  Language:     FORTRAN90                                                !
!                                                                         !
!  Description:  Prints formatted results in two-column format: distance  !
!                along bar (in meters) and temperature (in Kelvin)        !
!                                                                         !
!  Inputs:                                                                !
!                  n     Number of grid nodes                             !
!                  x     Vector of length n containing the location of    !
!                        nodes                                            !
!                  T     Temperature                                      !
!                                                                         !
!*************************************************************************!
SUBROUTINE WRITE_RESULTS(n,x,T)
!
IMPLICIT NONE
integer n,ii
real x(n)
double precision T(n)

open(unit=7,file='temp_distr.dat')
!
!   write header
write(7,101)
! loop over and write results
do ii = 1,n
   write(7,201)x(ii),T(ii)
end do
!
!...Format statements
!
101 format(5x,'____x(n)____',5x,'____T(n)____')
201 format(3x,f12.5,5x,f12.5)
!
end subroutine write_results