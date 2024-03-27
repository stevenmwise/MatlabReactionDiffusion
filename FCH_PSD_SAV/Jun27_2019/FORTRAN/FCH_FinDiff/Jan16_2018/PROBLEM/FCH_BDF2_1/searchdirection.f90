MODULE SearchDirectionRoutines
USE Global
IMPLICIT NONE
!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE SearchDirection(res,h)
USE Global
USE ProblemDef
USE, INTRINSIC:: iso_c_binding
!
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(INOUT):: res
REAL(KIND=r8), INTENT(IN):: h
!
INCLUDE "fftw3.f03"
TYPE(C_PTR):: plan_forward, plan_backward
INTEGER:: i, j, mh
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: xlen, ylen, tmp0, tmp1, tmp2, tmp3, tmp4
REAL(KIND=r8), DIMENSION(:), ALLOCATABLE:: eigx, eigy
REAL(KIND=r8), DIMENSION(:,:), ALLOCATABLE:: faux
COMPLEX(KIND=r8), DIMENSION(:,:), ALLOCATABLE:: aux1, aux2
!
mx(1) = SIZE(res,1); mx(2) = SIZE(res,1)
mh = mx(1)/2+1
!
ALLOCATE(aux1(1:mh   ,1:mx(2)), &
         aux2(1:mh   ,1:mx(2)))
ALLOCATE(faux(1:mx(1),1:mx(2)))
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
  tmp0 = SIN(REAL(i,KIND=r8)*pi*h/xlen)
  eigx(i) = 4.0_r8*tmp0*tmp0/(h*h)
END DO
!
! y-direction
!
DO j = 0, mx(2)
  tmp0 = SIN(REAL(j,KIND=r8)*pi*h/ylen)
  eigy(j) = 4.0_r8*tmp0*tmp0/(h*h)
END DO
!
faux(1:mx(1),1:mx(2)) = res(1:mx(1),1:mx(2))
!
plan_forward = FFTW_Plan_DFT_R2C_2D(mx(2),mx(1),faux(1:mx(1),1:mx(2)),aux1(1:mh,1:mx(2)),FFTW_ESTIMATE)
CALL FFTW_Execute_DFT_R2C(plan_forward,faux(1:mx(1),1:mx(2)),aux1(1:mh,1:mx(2)))
CALL FFTW_Destroy_Plan(plan_forward)
!
tmp1 = 4.0_r8*gamma+eta+6.0_r8*gamma*eps*eps+4.0_r8*split*gamma*eps*eps
tmp2 = 0.0_r8 ! 4.0_r8*split*gamma*eps*eps+6.0_r8*gamma*eps*eps
tmp3 = gamma*eps*eps*eps*eps
!
DO j = 1, mx(2)
  DO i = 1, mh
    IF(i == 1 .AND. j == 1) THEN
      aux2(i,j) = CMPLX(0.0_r8,KIND=r8)
    ELSE
      tmp0 = -(eigx(i-1)+eigy(j-1))
      tmp4 = -1.0_r8/(2.0_r8*dt/3.0_r8*tmp0)+tmp1-tmp2*tmp0+tmp3*tmp0*tmp0
      aux2(i,j) = aux1(i,j)/tmp4
    END IF
  END DO
END DO
!
plan_backward = FFTW_Plan_DFT_C2R_2D(mx(2),mx(1),aux2(1:mh,1:mx(2)),res(1:mx(1),1:mx(2)),FFTW_ESTIMATE)
CALL FFTW_Execute_DFT_C2R(plan_backward,aux2(1:mh,1:mx(2)),res(1:mx(1),1:mx(2)))
CALL FFTW_Destroy_Plan(plan_backward)
!
! Normalize:
res(1:mx(1),1:mx(2)) = res(1:mx(1),1:mx(2))/REAL(mx(1)*mx(2),KIND=r8)
!
DEALLOCATE(eigy,eigx,aux1,aux2,faux)
!
END SUBROUTINE SearchDirection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
END MODULE SearchDirectionRoutines
