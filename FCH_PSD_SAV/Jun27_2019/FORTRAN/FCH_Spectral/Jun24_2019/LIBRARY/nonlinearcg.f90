MODULE NonLinearCGRoutines
USE Global
IMPLICIT NONE
!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE NonLinearCG(u,uo,uoo,f,r,h,iter,energyvalue)
USE Global
USE Utilities, ONLY: Norm2D
USE SearchDirectionRoutines
USE ProblemNCG, ONLY: TrueEnergy, SetSource
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uoo
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN OUT):: f
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(OUT):: r
REAL(KIND=r8), INTENT(IN):: h
INTEGER, INTENT(OUT):: iter
REAL(KIND=r8), INTENT(OUT):: energyvalue
!
INTEGER:: its
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: h2, residual, res_search
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2  ,1:SIZE(u,2)-2):: d, s
!
mx(1) = SIZE(u,1)-2
mx(2) = SIZE(u,2)-2
!
h2 = h*h
!
CALL SetSource(  u(0:mx(1)+1,0:mx(2)+1), &
                uo(0:mx(1)+1,0:mx(2)+1), &
               uoo(0:mx(1)+1,0:mx(2)+1), &
                 f(1:mx(1)  ,1:mx(2)  ), &
                 s(1:mx(1)  ,1:mx(2)  ),h)
!
r(1:mx(1),1:mx(2)) = OpEnergyGradient(  u(0:mx(1)+1,0:mx(2)+1), &
                                       uo(0:mx(1)+1,0:mx(2)+1), &
                                      uoo(0:mx(1)+1,0:mx(2)+1), &
                                        s(1:mx(1)  ,1:mx(2)  ),h)
!
r(1:mx(1),1:mx(2)) = r(1:mx(1),1:mx(2))-SUM(r(1:mx(1),1:mx(2)))/(mx(1)*mx(2))
!
d(1:mx(1),1:mx(2)) = SearchDirection(r(1:mx(1),1:mx(2)))
!
! Probably not necessary:
d(1:mx(1),1:mx(2)) = d(1:mx(1),1:mx(2))-SUM(d(1:mx(1),1:mx(2)))/(mx(1)*mx(2))

DO its = 1, ncgmaxits
!
  iter = its
!
  CALL GetLineMinimizer(  u(0:mx(1)+1,0:mx(2)+1), &
                         uo(0:mx(1)+1,0:mx(2)+1), &
                        uoo(0:mx(1)+1,0:mx(2)+1), &
                          d(1:mx(1)  ,1:mx(2)  ), &
                          s(1:mx(1)  ,1:mx(2)  ),h,energyvalue)
!
  r(1:mx(1),1:mx(2)) = OpEnergyGradient(  u(0:mx(1)+1,0:mx(2)+1), &
                                         uo(0:mx(1)+1,0:mx(2)+1), &
                                        uoo(0:mx(1)+1,0:mx(2)+1), &
                                          s(1:mx(1)  ,1:mx(2)  ),h)                                      
!
  r(1:mx(1),1:mx(2)) = r(1:mx(1),1:mx(2))-SUM(r(1:mx(1),1:mx(2)))/(mx(1)*mx(2))
!
  residual = MAXVAL(ABS(r)) ! Norm2D(r)
  res_search = MAXVAL(ABS(d)) ! Norm2D(d)
!
  WRITE(*,1000) its, 'PSD Residual =', residual, 'PSD Search Res =', res_search
  1000 FORMAT(I4,3X,A17,1X,ES13.5,A19,1X,ES13.5)
!
  IF(residual < ncgtol) RETURN
!
  IF(its == ncgmaxits) EXIT
!
  d(1:mx(1),1:mx(2)) = SearchDirection(r(1:mx(1),1:mx(2)))
!   
  d(1:mx(1),1:mx(2)) = d(1:mx(1),1:mx(2))-SUM(d(1:mx(1),1:mx(2)))/(mx(1)*mx(2))
!
END DO
!
PRINT *, 'NonLinearCG: maximum iterations exceeded.'
!
END SUBROUTINE NonLinearCG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE GetLineMinimizer(u,uo,uoo,r,s,h,energyminimum)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uoo
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN OUT):: r
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: s
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), INTENT(OUT):: energyminimum
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: x0, x1, xmin
!
mx(1) = SIZE(u,1)-2
mx(2) = SIZE(u,2)-2
!
! Starting points. Are there better ones?
x0 = 1.0_r8
x1 = 0.9_r8
!
! Secant Method:
!energyminimum = SecantMethod(  u(0:mx(1)+1,0:mx(2)+1), &
!                              uo(0:mx(1)+1,0:mx(2)+1), &
!                             uoo(0:mx(1)+1,0:mx(2)+1), &
!                               r(1:mx(1)  ,1:mx(2)  ), &
!                               s(1:mx(1)  ,1:mx(2)  ),h,x0,x1,xmin)
!
xmin = 0.7_r8
r(1:mx(1),1:mx(2)) = xmin*r(1:mx(1),1:mx(2))
u(1:mx(1),1:mx(2)) = u(1:mx(1),1:mx(2))+r(1:mx(1),1:mx(2))
!
END SUBROUTINE GetLineMinimizer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION SecantMethod(u,uo,uoo,r,s,h,x0,x1,xmin) RESULT(secantresult)
USE Global
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: u
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uoo
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN OUT):: r
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: s
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), INTENT(IN):: x0
REAL(KIND=r8), INTENT(IN):: x1
REAL(KIND=r8), INTENT(OUT):: xmin
REAL(KIND=r8):: secantresult
!
INTEGER:: iter
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: p0, p1, p2, q, q0, q1
!
mx(1) = SIZE(u,1)-2
mx(2) = SIZE(u,2)-2
!
p0 = x0
p1 = x1
!
q0 = LineOpEnergyGradient(p0,u(0:mx(1)+1,0:mx(2)+1), &
                            uo(0:mx(1)+1,0:mx(2)+1), &
                           uoo(0:mx(1)+1,0:mx(2)+1), &
                             r(1:mx(1)  ,1:mx(2)  ), &
                             s(1:mx(1)  ,1:mx(2)  ),h)
q1 = LineOpEnergyGradient(p1,u(0:mx(1)+1,0:mx(2)+1), &
                            uo(0:mx(1)+1,0:mx(2)+1), &
                           uoo(0:mx(1)+1,0:mx(2)+1), &
                             r(1:mx(1)  ,1:mx(2)  ), &
                             s(1:mx(1)  ,1:mx(2)  ),h)
!
DO iter = 1, secantmaxits
!
q = q1-q0
!
IF(q == 0.0_r8) THEN
  p2 = p1
  EXIT
ELSE
  p2 = p1-q1*(p1-p0)/(q1-q0)
END IF
!
IF(debug) PRINT *, 'Secant iteration:', iter, p2
!
! Do at least 2 iterations:
IF(iter >= 2 .AND. ABS(p2-p1) < secanttol) EXIT
!
p0 = p1
q0 = q1
p1 = p2
q1 = LineOpEnergyGradient(p2,u(0:mx(1)+1,0:mx(2)+1), &
                            uo(0:mx(1)+1,0:mx(2)+1), &
                           uoo(0:mx(1)+1,0:mx(2)+1), &
                             r(1:mx(1)  ,1:mx(2)  ), &
                             s(1:mx(1)  ,1:mx(2)  ),h)
!
END DO
!
xmin = p2
!secantresult = LineOpEnergy(p2,u(0:mx(1)+1,0:mx(2)+1), &
!                              uo(0:mx(1)+1,0:mx(2)+1), &
!                             uoo(0:mx(1)+1,0:mx(2)+1), &
!                               r(1:mx(1)  ,1:mx(2)  ), &
!                               s(1:mx(1)  ,1:mx(2)  ),h)
!
secantresult = 0.0_r8
!
END FUNCTION SecantMethod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION LineOpEnergy(alpha,u,uo,uoo,r,s,h) RESULT(lineenergyresult)
USE Global
USE ProblemNCG, ONLY: OpEnergy
IMPLICIT NONE
!
! OpEnergy along the line p+alpha*r, where alpha is a scalar.
!
REAL(KIND=r8), INTENT(IN):: alpha
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: u
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uoo
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN OUT):: r
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: s
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8):: lineenergyresult
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2+1,0:SIZE(u,2)-2+1):: utmp
!
mx(1) = SIZE(u,1)-2
mx(2) = SIZE(u,2)-2
!
! B.C.'s for utmp are imposed in OpEnergy:
utmp(1:mx(1),1:mx(2)) = u(1:mx(1),1:mx(2))+alpha*r(1:mx(1),1:mx(2))
!
lineenergyresult = OpEnergy(utmp(0:mx(1)+1,0:mx(2)+1), &
                              uo(0:mx(1)+1,0:mx(2)+1), &
                             uoo(0:mx(1)+1,0:mx(2)+1), &
                               s(1:mx(1)  ,1:mx(2)  ),h)
!
END FUNCTION LineOpEnergy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION LineOpEnergyGradient(alpha,u,uo,uoo,r,s,h) RESULT(linegradientresult)
USE Global
IMPLICIT NONE
!
! OpEnergy gradient along the line p+alpha*r, where alpha is a scalar.
!
REAL(KIND=r8), INTENT(IN):: alpha
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: u
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uoo
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN OUT):: r
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: s
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8):: linegradientresult
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2+1,0:SIZE(u,2)-2+1):: utmp
!
mx(1) = SIZE(u,1)-2
mx(2) = SIZE(u,2)-2
!
! B.C.'s for utmp are imposed in OpEnergyGradient:
utmp(1:mx(1),1:mx(2)) = u(1:mx(1),1:mx(2))+alpha*r(1:mx(1),1:mx(2))
!
linegradientresult = SUM(OpEnergyGradient(utmp(0:mx(1)+1,0:mx(2)+1), &
                                            uo(0:mx(1)+1,0:mx(2)+1), &
                                           uoo(0:mx(1)+1,0:mx(2)+1), &
                                             s(1:mx(1)  ,1:mx(2)  ),h) &
                   *                         r(1:mx(1)  ,1:mx(2)  ))
!
END FUNCTION LineOpEnergyGradient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION OpEnergyGradient(u,uo,uoo,s,h) RESULT(energygradientresult)
USE Global
USE ProblemNCG, ONLY: Operator 
IMPLICIT NONE
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN OUT):: u
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uo
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: uoo
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN):: s
REAL(KIND=r8), INTENT(IN):: h
REAL(KIND=r8), DIMENSION(1:SIZE(u,1)-2,1:SIZE(u,2)-2):: energygradientresult
!
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(u,1)-2
mx(2) = SIZE(u,2)-2
!
! B.C.'s for u are imposed in Operator:
energygradientresult(1:mx(1),1:mx(2)) =            s(1:mx(1)  ,1:mx(2)  ) &
                                      - Operator(  u(0:mx(1)+1,0:mx(2)+1), &
                                                  uo(0:mx(1)+1,0:mx(2)+1), &
                                                 uoo(0:mx(1)+1,0:mx(2)+1),h)
!
!
END FUNCTION OpEnergyGradient
!
END MODULE NonLinearCGRoutines
