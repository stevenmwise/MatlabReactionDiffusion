MODULE PoissonFFTSolverRoutines
USE Global
IMPLICIT NONE
!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION InvNegSPLap(f) RESULT(invresult)
USE Global
USE, INTRINSIC:: iso_c_binding
!
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), DIMENSION(1:SIZE(f,1),1:SIZE(f,2)):: invresult
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
mx(1) = SIZE(f,1); mx(2) = SIZE(f,2)
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
DO i = 0, mh-1
  eigx(i) = 4.0_r8*REAL(i,KIND=r8)*REAL(i,KIND=r8)*pi*pi/(xlen*xlen) 
END DO
!
! y-direction
!
DO j = 0, mh-1
  eigy(j) = 4.0_r8*REAL(j,KIND=r8)*REAL(j,KIND=r8)*pi*pi/(ylen*ylen) 
END DO
DO j = mh, mx(2)
  eigy(j) = 4.0_r8*REAL(mx(2)-j,KIND=r8)*REAL(mx(2)-j,KIND=r8)*pi*pi/(ylen*ylen) 
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
plan_backward = FFTW_Plan_DFT_C2R_2D(mx(2),mx(1),aux(1:mh,1:mx(2)),invresult(1:mx(1),1:mx(2)),FFTW_ESTIMATE)
CALL FFTW_Execute_DFT_C2R(plan_backward,aux(1:mh,1:mx(2)),invresult(1:mx(1),1:mx(2)))
CALL FFTW_Destroy_Plan(plan_backward)
!
! Normalize:
invresult(1:mx(1),1:mx(2)) = invresult(1:mx(1),1:mx(2))/REAL(mx(1)*mx(2),KIND=r8)
!
DEALLOCATE(eigy,eigx,aux,ftmp)
!
END Function InvNegSPLap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION SPLap(f) RESULT(lapresult)
USE Global
USE, INTRINSIC:: iso_c_binding
!
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), DIMENSION(1:SIZE(f,1),1:SIZE(f,2)):: lapresult
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
mx(1) = SIZE(f,1); mx(2) = SIZE(f,2)
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
DO i = 0, mh-1
  eigx(i) = -4.0_r8*REAL(i,KIND=r8)*REAL(i,KIND=r8)*pi*pi/(xlen*xlen) 
END DO
!
! y-direction
!
DO j = 0, mh-1
  eigy(j) = -4.0_r8*REAL(j,KIND=r8)*REAL(j,KIND=r8)*pi*pi/(ylen*ylen) 
END DO
DO j = mh, mx(2)
  eigy(j) = -4.0_r8*REAL(mx(2)-j,KIND=r8)*REAL(mx(2)-j,KIND=r8)*pi*pi/(ylen*ylen) 
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
    aux(i,j) = aux(i,j)*(eigx(i-1)+eigy(j-1))
  END DO
END DO
!
plan_backward = FFTW_Plan_DFT_C2R_2D(mx(2),mx(1),aux(1:mh,1:mx(2)),lapresult(1:mx(1),1:mx(2)),FFTW_ESTIMATE)
CALL FFTW_Execute_DFT_C2R(plan_backward,aux(1:mh,1:mx(2)),lapresult(1:mx(1),1:mx(2)))
CALL FFTW_Destroy_Plan(plan_backward)
!
! Normalize:
lapresult(1:mx(1),1:mx(2)) = lapresult(1:mx(1),1:mx(2))/REAL(mx(1)*mx(2),KIND=r8)
!
DEALLOCATE(eigy,eigx,aux,ftmp)
!
END Function SPLap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE PoissonFFTSolverRoutines
