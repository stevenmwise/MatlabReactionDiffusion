MODULE ProblemNCG
USE Global
!
! FCH Equation
!
! BDF2 formulation.
!
! Nonlinear Conjugate Gradient Solver.
!
CONTAINS
!
SUBROUTINE SetProblem
USE Global
USE ProblemDef
IMPLICIT NONE
!
INTEGER:: ierror
NAMELIST/problemdata/ eps, eta, gamma, split
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
PRINT *, 'eps:', eps, 'gamma:', gamma
PRINT *, 'eta:', eta, 'split:', split
!
END SUBROUTINE SetProblem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Initialize(u,h)
USE Global
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(OUT):: u
REAL(KIND=r8), INTENT(IN):: h
!
INTEGER:: i, j
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: x, y
REAL(KIND=r8), DIMENSION(1:2):: c, p
!
mx(1) = SIZE(u,1)-2
mx(2) = SIZE(u,2)-2
!
p(1:2) = mx(1:2)*h
c(1:2) = (xlower(1:2)+(xlower(1:2)+p(1:2)))/2.0_r8
!
DO j = 1, mx(2)
  y = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(2)
  DO i = 1, mx(1)
    x = (REAL(i,KIND=r8)-0.5_r8)*h+xlower(1)
!
    IF(x > SIN(4.0*pi*y/12.8)+6.4+0.34) THEN
      u(i,j) = -1.0
    ELSE IF(x < SIN(4.0*pi*y/12.8)+6.4-0.34) THEN
      u(i,j) = -1.0
    ELSE
      u(i,j) = 1.0
    END IF      
!     
  END DO
END DO
!
CALL BoundaryConditions(u)
!
END SUBROUTINE Initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION TrueEnergy(u,f,h) RESULT(energyresult)
USE Global
USE Utilities, ONLY: ULap2D
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
INTEGER:: i, j
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: f1, f2, h2, tmp1, tmp2, tmp3, tmp4
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2  ,1:SIZE(u,2)-2  ):: plapenergy, &
                                                             lapenergy, &
                                                             kapenergy
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2  ,0:SIZE(u,2)-2  ):: r
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2  ,1:SIZE(u,2)-2  ):: dx1u
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2  ,0:SIZE(u,2)-2  ):: dx2u
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2+1,0:SIZE(u,2)-2+1):: kap
!
energyresult = 0.0_r8
!
END FUNCTION TrueEnergy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION OpEnergy(u,uo,uoo,s,h) RESULT(energyresult)
USE Global
USE Utilities, ONLY: ULap2D
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
REAL(KIND=r8):: energyresult, energyminimum
!
INTEGER:: i, j, iter
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: f1, f2, h2, tmp1, tmp2, tmp3, tmp4
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,1:SIZE(u,2)-2):: plapenergy, kapenergy, f, res
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2,0:SIZE(u,2)-2):: r
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2+1,0:SIZE(u,2)-2+1):: lap
!
energyresult = 0.0_r8
!
END FUNCTION OpEnergy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Operator(u,uo,uoo,h) RESULT(operatorresult)
USE Global
USE Utilities, ONLY: ULap2D, ULap2DDiag, ULap2DDiagSpec, ULap2DSpec, Sqr, Cube, Quad, &
                     Quin, C2VAve, V2CAve
USE PoissonFFTSolverRoutines, ONLY: PoissonFFTSolver
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
INTEGER:: i, j, iter
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: h2, g1, g2, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8
!REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2  ,0:SIZE(u,2)-2  ):: r
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2  ,1:SIZE(u,2)-2  ):: m1
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2  ,0:SIZE(u,2)-2  ):: m2
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2+1,0:SIZE(u,2)-2+1):: kap, finv
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2  ,1:SIZE(u,2)-2  ):: f
!
mx(1) = SIZE(u,1)-2
mx(2) = SIZE(u,2)-2
!
h2 = h*h
tmp1 = 1.0_r8/(2.0_r8*h)
tmp2 = 3.0_r8*gamma
tmp3 = gamma+eta
tmp4 = 4.0_r8*split*gamma*eps2
tmp5 = gamma*Quad(eps)
tmp6 = 6.0_r8*gamma*eps2
tmp7 = 4.0_r8*gamma+eta+4.0_r8*split*gamma*eps2
tmp8 = (2.0_r8*gamma+eta)*eps2
!
CALL BoundaryConditions(u(0:mx(1)+1,0:mx(2)+1))
!
kap(1:mx(1),1:mx(2)) = ULap2D(u(0:mx(1)+1,0:mx(2)+1))/h2
!
CALL BoundaryConditions(kap(0:mx(1)+1,0:mx(2)+1))
!
f(1:mx(1),1:mx(2)) = (u(1:mx(1),1:mx(2))-(4.0_r8/3.0_r8)*uo (1:mx(1),1:mx(2))  &
                   +                     (1.0_r8/3.0_r8)*uoo(1:mx(1),1:mx(2))) &
                   / ((2.0_r8/3.0_r8)*dt)
f(1:mx(1),1:mx(2)) = f(1:mx(1),1:mx(2))-SUM(f(1:mx(1),1:mx(2))) &
                   / REAL(mx(1)*mx(2),KIND=r8)
!
CALL PoissonFFTSolver(finv,f,h)
!
!DO i = 0, mx(1)
!  DO j = 0, mx(2)
!
!    g1 = tmp1*(u(i+1,j+1)-u(i,j)-u(i,j+1)+u(i+1,j))
!    g2 = tmp1*(u(i,j+1)-u(i+1,j)+u(i+1,j+1)-u(i,j))
!    r(i,j) = g1*g1+g2*g2
!
!  END DO
!END DO
!
m1(0:mx(1),1:mx(2)) = (Sqr(u(0:mx(1),1:mx(2)))+Sqr(u(1:mx(1)+1,1:mx(2)  )))/2.0_r8
m2(1:mx(1),0:mx(2)) = (Sqr(u(1:mx(1),0:mx(2)))+Sqr(u(1:mx(1)  ,1:mx(2)+1)))/2.0_r8
!
operatorresult(1:mx(1),1:mx(2)) &
  = finv(1:mx(1),1:mx(2)) &
  + tmp2*Quin(u(1:mx(1),1:mx(2))) &
  + tmp3*u(1:mx(1),1:mx(2)) &
  + tmp4*Cube(u(1:mx(1),1:mx(2))) &
  + tmp5*ULap2D(kap(0:mx(1)+1,0:mx(2)+1))/h2 &
!  + tmp6*u(1:mx(1),1:mx(2))*V2CAve(GradSqr(u(0:mx(1)+1,0:mx(2)+1),h)) &
!  - tmp6*ULap2DDiagSpec(C2VAve(Sqr(u(0:mx(1)+1,0:mx(2)+1))),u(0:mx(1)+1,0:mx(2)+1))/h2 &
  + tmp6*u(1:mx(1),1:mx(2))*CGradSqr(u(0:mx(1)+1,0:mx(2)+1),h) &
  - tmp6*ULap2DSpec(m1(0:mx(1),1:mx(2)),m2(1:mx(1),0:mx(2)),u(0:mx(1)+1,0:mx(2)+1))/h2 &
!
! This is added for energy stability:
!  - tmp4*ULap2DDiagSpec(           r(0:mx(1)  ,0:mx(2)  )  ,u(0:mx(1)+1,0:mx(2)+1))/h2 &
!
! This part is usually explicit:
  - tmp7*Cube(u(1:mx(1),1:mx(2))) &
!  + tmp8*ULap2DDiag(u(0:mx(1)+1,0:mx(2)+1))/h2
  + tmp8*ULap2D(u(0:mx(1)+1,0:mx(2)+1))/h2
!
END FUNCTION Operator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE SetSource(u,uo,uoo,f,s,h)
USE Global
USE Utilities, ONLY: ULap2D, ULap2DDiag, ULap2DDiagSpec, Cube
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
REAL(KIND=r8):: g1, g2, h2, tmp1, tmp2, tmp3, tmp4
!
!mx(1) = SIZE(u,1)-2
!mx(2) = SIZE(u,2)-2
!
!h2 = h*h
!tmp1 = 1.0_r8/(2.0_r8*h)
!tmp2 = 4.0_r8*gamma+eta+4.0_r8*split*gamma*eps2
!tmp3 = (2.0_r8*gamma+eta)*eps2
!tmp4 = 4.0_r8*split*gamma*eps2
!
!DO i = 0, mx(1)
!  DO j = 0, mx(2)
!
!    g1 = tmp1*(uo(i+1,j+1)-uo(i,j)-uo(i,j+1)+uo(i+1,j))
!    g2 = tmp1*(uo(i,j+1)-uo(i+1,j)+uo(i+1,j+1)-uo(i,j))
!    r(i,j) = g1*g1+g2*g2
!
!  END DO
!END DO
!
!s(1:mx(1),1:mx(2)) = tmp2*Cube(uo(1:mx(1),1:mx(2))) &
!                   - tmp3*ULap2DDiag(uo(0:mx(1)+1,0:mx(2)+1))/h2 &
!                   - tmp4*ULap2DDiagSpec(r(0:mx(1),0:mx(2)),uo(0:mx(1)+1,0:mx(2)+1))/h2
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
FUNCTION GradSqr(u,h) RESULT(grad2)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: u
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2,0:SIZE(u,2)-2):: grad2
!
INTEGER:: i, j
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: g1, g2, tmp
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
tmp = 0.5_r8/h
!
DO i = 0, mx(1)
  DO j = 0, mx(2)
!
    g1 = tmp*(u(i+1,j+1)-u(i,j+1)+u(i+1,j)-u(i,j))
    g2 = tmp*(u(i+1,j+1)-u(i+1,j)+u(i,j+1)-u(i,j))
    grad2(i,j) = g1*g1+g2*g2
!
  END DO
END DO
!
END FUNCTION GradSqr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION CGradSqr(u,h) RESULT(grad2)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: u
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,1:SIZE(u,2)-2):: grad2
!
INTEGER:: i, j
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: g1, g2
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2,1:SIZE(u,2)-2):: grad_x_2
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,0:SIZE(u,2)-2):: grad_y_2
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
DO i = 0, mx(1)
  DO j = 1, mx(2)
    g1 = (u(i+1,j)-u(i,j))/h
    grad_x_2(i,j) = g1*g1
  END DO
END DO
!
DO i = 1, mx(1)
  DO j = 0, mx(2)
    g2 = (u(i,j+1)-u(i,j))/h
    grad_y_2(i,j) = g2*g2
  END DO
END DO
!
DO i = 1, mx(1)
  DO j = 1, mx(2)
    grad2(i,j) = (grad_x_2(i,j)+grad_x_2(i-1,j)+grad_y_2(i,j)+grad_y_2(i,j-1))/2.0_r8
  END DO
END DO
!
END FUNCTION CGradSqr
END MODULE ProblemNCG
