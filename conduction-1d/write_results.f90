subroutine write_results(n,x,T)
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
201 format(3x,f12.5,3x,f12.5)
!
end subroutine write_results
