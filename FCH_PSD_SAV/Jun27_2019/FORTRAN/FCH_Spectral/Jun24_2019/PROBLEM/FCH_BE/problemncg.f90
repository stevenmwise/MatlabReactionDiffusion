MODULE ProblemNCG
USE Global
!
! FCH Equation
!
! Backward Euler formulation.
!
! Preconditioned Nonlinear Conjugate Gradient Solver.
!
CONTAINS
!
SUBROUTINE SetProblem
USE Global
USE ProblemDef
IMPLICIT NONE
!
INTEGER:: ierror
NAMELIST/problemdata/ eps, eta1, eta2, gamma
!
OPEN(UNIT=75,FILE='problemdata.dat',STATUS='OLD',ACTION='READ',IOSTAT=ierror)
IF(ierror/=0) THEN
  PRINT *,'Error opening input file problemdata.dat. Program stop.'
  STOP
END IF
READ(75,NML=problemdata)
CLOSE(75)
OPEN(UNIT=76,FILE='output.dat',STATUS='OLD',ACTION='WRITE',FORM='FORMATTED', &
     POSITION='APPEND')
WRITE(76,NML=problemdata)
CLOSE(76)
!
eps2 = eps*eps
!
END SUBROUTINE SetProblem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Initialize(u,h)
USE Global
USE Utilities, ONLY: Sqr
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(OUT):: u
REAL(KIND=r8), INTENT(IN):: h
!
INTEGER:: i, j
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: x, y, r
REAL(KIND=r8), DIMENSION(1:2):: c, p
!
mx(1) = SIZE(u,1)-2
mx(2) = SIZE(u,2)-2
!
p(1:2) = mx(1:2)*h
c(1:2) = (xlower(1:2)+(xlower(1:2)+p(1:2)))/2.0_r8
!
u(1:mx(1),1:mx(2)) = SmoothedAnnulus(p,mx)
!
CALL BoundaryConditions(u)
!
END SUBROUTINE Initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION TrueEnergy(u,f,h) RESULT(energyresult)
USE Global
USE ProblemDef
IMPLICIT NONE
!
! True energy.
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
!REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: q
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8):: energyresult
!
energyresult = 0.0_r8
!
END FUNCTION TrueEnergy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION OpEnergy(u,uo,uoo,s,h) RESULT(energyresult)
USE Global
USE ProblemDef
IMPLICIT NONE
!
! Operator energy.  Not necessary for computation
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uoo
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: s
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8):: energyresult
!
energyresult = 0.0_r8
!
END FUNCTION OpEnergy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Operator(u,uo,uoo,h) RESULT(operatorresult)
USE Global
USE PoissonFFTSolverRoutines, ONLY: InvNegSPLap, SPLap
USE ProblemDef
IMPLICIT NONE
!
! Operator:
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uoo
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,1:SIZE(u,2)-2):: operatorresult
!
INTEGER:: i, j
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: h2
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2  ,1:SIZE(u,2)-2  ):: kap, finv
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2  ,1:SIZE(u,2)-2  ):: f, res
!
mx(1) = SIZE(u,1)-2
mx(2) = SIZE(u,2)-2
!
h2 = h*h
!
CALL BoundaryConditions(u(0:mx(1)+1,0:mx(2)+1))
!
f(1:mx(1),1:mx(2)) = (u(1:mx(1),1:mx(2))-uo(1:mx(1),1:mx(2)))/dt
f(1:mx(1),1:mx(2)) = f(1:mx(1),1:mx(2))-SUM(f(1:mx(1),1:mx(2))) &
                   / REAL(mx(1)*mx(2),KIND=r8)
!
finv(1:mx(1),1:mx(2)) = InvNegSPLap(f)
!
kap(1:mx(1),1:mx(2)) = eps2*SPLap(u(1:mx(1),1:mx(2)))-DF(u(1:mx(1),1:mx(2)))
!
operatorresult(1:mx(1),1:mx(2)) &
  = finv(1:mx(1),1:mx(2)) &
  + eps2*SPLap(kap(1:mx(1),1:mx(2))) &
  - kap(1:mx(1),1:mx(2))*(D2F(u(1:mx(1),1:mx(2)))-eta1) &
  + (eta1-eta2)*DF(u(1:mx(1),1:mx(2)))
!
END FUNCTION Operator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE SetSource(u,uo,uoo,f,s,h)
USE Global
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: u
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uoo
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: f
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN OUT):: s
REAL(KIND=r8), INTENT(IN):: h
!
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2,0:SIZE(u,2)-2):: r
!
INTEGER:: i, j
INTEGER, DIMENSION(1:2):: mx
!
s = 0.0_r8
!
END SUBROUTINE SetSource
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE BoundaryConditions(u)
USE Global
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
!
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(u,1)-2
mx(2) = SIZE(u,2)-2
!
! Periodic boundary conditions: Order is important!
u(      0,1:mx(2)) = u(mx(1),1:mx(2))
u(mx(1)+1,1:mx(2)) = u(    1,1:mx(2))
!
u(0:mx(1)+1,      0) = u(0:mx(1)+1,mx(2))
u(0:mx(1)+1,mx(2)+1) = u(0:mx(1)+1,    1)
!
END SUBROUTINE BoundaryConditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION DF(p) RESULT(dfresult)
USE Global
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: p
REAL(KIND=r8):: dfresult
!
dfresult = (p-1.0_r8)*(p+1.0_r8)*(p+0.5*gamma)
!
END FUNCTION DF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION D2F(p) RESULT(d2fresult)
USE Global
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: p
REAL(KIND=r8):: d2fresult
!
d2fresult = (p-1.0_r8)*(p+1.0_r8)               &
          + (p-1.0_r8)*           (p+0.5*gamma) &
          +            (p+1.0_r8)*(p+0.5*gamma)
!
END FUNCTION D2F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION SmoothedAnnulus(p,mx) RESULT(uresult)
USE Global
USE Utilities, ONLY: Sqr
USE ProblemDef
USE, INTRINSIC:: iso_c_binding
!
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(1:2), INTENT(IN):: p
INTEGER, DIMENSION(1:2), INTENT(IN):: mx
REAL(KIND=r8), DIMENSION(1:mx(1),1:mx(2)):: uresult
!
INCLUDE "fftw3.f03"
TYPE(C_PTR):: plan_forward, plan_backward
INTEGER:: i, j, mmh
INTEGER, DIMENSION(1:2):: factor, mmx
REAL(KIND=r8):: alpha, r, x, y
REAL(KIND=r8), DIMENSION(1:2):: hm
REAL(KIND=r8), DIMENSION(:), ALLOCATABLE:: eigx, eigy
REAL(KIND=r8), DIMENSION(:,:), ALLOCATABLE:: um
COMPLEX(KIND=r8), DIMENSION(:,:), ALLOCATABLE:: aux1, aux2
!
! Ultrafine Grid:
mmx(1:2) = 1024
mmh = mmx(1)/2+1
!
IF(MODULO(mmx(1),mx(1)) == 0) THEN
  factor(1) = mmx(1)/mx(1)
ELSE
  PRINT *, 'There is an error.'
  STOP
END IF
!
IF(MODULO(mmx(2),mx(2)) == 0) THEN
  factor(2) = mmx(2)/mx(2)
ELSE
  PRINT *, 'There is an error.'
  STOP
END IF
!
ALLOCATE(um(  1:mmx(1),1:mmx(2)))
ALLOCATE(aux1(1:mmh   ,1:mmx(2)), &
         aux2(1:mmh   ,1:mmx(2)))
ALLOCATE(eigx(0:mmx(1)),eigy(0:mmx(2)))
!
um = 0.0_r8
hm(1:2) = p(1:2)/mmx(1:2)
!
DO i = 1, mmx(1)
  x = REAL(i,KIND=r8)*hm(1)
  DO j = 1, mmx(2)
    y = REAL(j,KIND=r8)*hm(2)
    r = SQRT(Sqr(x-p(1)/2.0_r8)+Sqr(y-p(2)/2.0_r8)/2.0_r8)
    IF(r>p(1)/4.0_r8+0.2_r8) THEN
      um(i,j) = -1.0_r8
    ELSE IF(r<p(1)/4.0_r8-0.2_r8) THEN
      um(i,j) = -1.0
    ELSE
      um(i,j) =  1.0
    END IF
  END DO
END DO
!
alpha = 50.0_r8*LOG(10.0_r8)
!
! Eigenvalues
!
! x-direction
DO i = 0, mmh-1
  eigx(i) = Sqr(2.0_r8*REAL(i,KIND=r8)/mmx(1))
END DO
!
! y-direction
DO j = 0, mmh-1
  eigy(j) = Sqr(2.0_r8*REAL(j,KIND=r8)/mmx(2))
END DO
DO j = mmh, mmx(2)
  eigy(j) = Sqr(2.0_r8*REAL(mmx(2)-j,KIND=r8)/mmx(2))
END DO
!
plan_forward = FFTW_Plan_DFT_R2C_2D(mmx(2),mmx(1),um(1:mmx(1),1:mmx(2)),aux1(1:mmh,1:mmx(2)),FFTW_ESTIMATE)
CALL FFTW_Execute_DFT_R2C(plan_forward,um(1:mmx(1),1:mmx(2)),aux1(1:mmh,1:mmx(2)))
CALL FFTW_Destroy_Plan(plan_forward)
!
DO j = 1, mmx(2)
  DO i = 1, mmh
    aux2(i,j) = EXP(-(eigx(i-1)+eigy(j-1))*alpha)*aux1(i,j)
  END DO
END DO
!
plan_backward = FFTW_Plan_DFT_C2R_2D(mmx(2),mmx(1),aux2(1:mmh,1:mmx(2)),um(1:mmx(1),1:mmx(2)),FFTW_ESTIMATE)
CALL FFTW_Execute_DFT_C2R(plan_backward,aux2(1:mmh,1:mmx(2)),um(1:mmx(1),1:mmx(2)))
CALL FFTW_Destroy_Plan(plan_backward)
!
! Normalize:
um(1:mmx(1),1:mmx(2)) = um(1:mmx(1),1:mmx(2))/REAL(mmx(1)*mmx(2),KIND=r8)
!
DO i = 1, mx(1)
  DO j = 1, mx(2)
    uresult(i,j) = um(i*factor(1),j*factor(2))
  END DO
END DO
!
DEALLOCATE(eigy,eigx,aux1,aux2,um)
!
END FUNCTION SmoothedAnnulus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE ProblemNCG
