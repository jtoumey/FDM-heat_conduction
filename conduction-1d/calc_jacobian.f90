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
SUBROUTINE CALC_JACOBIAN()


END SUBROUTINE CALC_JACOBIAN