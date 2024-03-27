MODULE PoissonFFTSolverRoutines
USE Global
IMPLICIT NONE
!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE PoissonFFTSolver(u,f,h)
USE Global
USE, INTRINSIC:: iso_c_binding
!
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), INTENT(IN):: h
!
INCLUDE "fftw3.f03"
TYPE(C_PTR):: plan_forward, plan_backward
INTEGER:: i, j, mh
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: xlen, ylen, tmp
REAL(KIND=r8), DIMENSION(:), ALLOCATABLE:: eigx, eigy
REAL(KIND=r8), DIMENSION(:,:), ALLOCATABLE:: ftmp
COMPLEX(KIND=r8), DIMENSION(:,:), ALLOCATABLE:: aux
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
mh = mx(1)/2+1
!
ALLOCATE(aux( 1:mh   ,1:mx(2)))
ALLOCATE(ftmp(1:mx(1),1:mx(2)))
ALLOCATE(eigx(0:mx(1)),eigy(0:mx(2)))
!
xlen = xupper(1)-xlower(1)
ylen = xupper(2)-xlower(2)
!
! Eigenvalues
!
! x-direction
!
DO i = 0, mx(1)
  tmp = SIN(REAL(i,KIND=r8)*pi*h/xlen)
  eigx(i) = 4.0_r8*tmp*tmp/(h*h)
END DO
!
! y-direction
!
DO j = 0, mx(2)
  tmp = SIN(REAL(j,KIND=r8)*pi*h/ylen)
  eigy(j) = 4.0_r8*tmp*tmp/(h*h)
END DO
!
ftmp(1:mx(1),1:mx(2)) = f(1:mx(1),1:mx(2))
!
plan_forward = FFTW_Plan_DFT_R2C_2D(mx(2),mx(1),ftmp(1:mx(1),1:mx(2)),aux(1:mh,1:mx(2)),FFTW_ESTIMATE)
CALL FFTW_Execute_DFT_R2C(plan_forward,ftmp(1:mx(1),1:mx(2)),aux(1:mh,1:mx(2)))
CALL FFTW_Destroy_Plan(plan_forward)
!
DO j = 1, mx(2)
  DO i = 1, mh
    IF(i == 1 .AND. j == 1) THEN
      aux(i,j) = CMPLX(0.0_r8,KIND=r8)
    ELSE
      aux(i,j) = aux(i,j)/(eigx(i-1)+eigy(j-1))
    END IF
  END DO
END DO
!
plan_backward = FFTW_Plan_DFT_C2R_2D(mx(2),mx(1),aux(1:mh,1:mx(2)),u(1:mx(1),1:mx(2)),FFTW_ESTIMATE)
CALL FFTW_Execute_DFT_C2R(plan_backward,aux(1:mh,1:mx(2)),u(1:mx(1),1:mx(2)))
CALL FFTW_Destroy_Plan(plan_backward)
!
! Normalize:
u(1:mx(1),1:mx(2)) = u(1:mx(1),1:mx(2))/REAL(mx(1)*mx(2),KIND=r8)
! 
!Periodic boundary conditions: Order is important.
u(      0,1:mx(2)) = u(mx(1),1:mx(2))
u(mx(1)+1,1:mx(2)) = u(    1,1:mx(2))

u(0:mx(1)+1,      0) = u(0:mx(1)+1,mx(2))
u(0:mx(1)+1,mx(2)+1) = u(0:mx(1)+1,    1)
!
DEALLOCATE(eigy,eigx,aux,ftmp)
!
END SUBROUTINE PoissonFFTSolver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE PoissonFFTSolverRoutines
