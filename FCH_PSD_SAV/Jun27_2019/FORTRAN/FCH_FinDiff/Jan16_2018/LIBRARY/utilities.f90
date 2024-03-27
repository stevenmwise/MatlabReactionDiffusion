MODULE Utilities
USE Global
IMPLICIT NONE
!
CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Norm2D(g) RESULT(normresult)
USE Global
IMPLICIT NONE
REAL(KIND=r8), DIMENSION(1:,1:), INTENT(IN) :: g
REAL(KIND=r8) :: normresult
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: tmp
!
mx(1) = SIZE(g,1); mx(2) = SIZE(g,2)
tmp = REAL(mx(1)*mx(2),KIND=r8)
!
normresult = SQRT(SUM(g*g)/tmp)
!
END FUNCTION Norm2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION Norm2DVars(g) RESULT(normresult)
USE Global
IMPLICIT NONE
REAL(KIND=r8), DIMENSION(1:,1:,1:), INTENT(IN) :: g
REAL(KIND=r8) :: normresult
!
INTEGER:: numvars_loc
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: tmp
!
mx(1) = SIZE(g,1); mx(2) = SIZE(g,2); numvars_loc = SIZE(g,3) 
tmp = REAL(mx(1)*mx(2)*numvars_loc,KIND=r8)
!
normresult = SQRT(SUM(g*g)/tmp)
!
END FUNCTION Norm2DVars
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION C2VAve(u) RESULT(averesult)
USE Global
IMPLICIT NONE
!
! Center-to-vertex average. The variable v(0:mx(1)+1,0:mx(2)+1) is a cell-
! centered function. The result averesult(0:mx(1),0:mx(2)) is vertex-centered.
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: u
REAL(KIND=r8), DIMENSION(0:SIZE(u,1)-2,0:SIZE(u,2)-2):: averesult
!
INTEGER:: i, j
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(u,1)-2; mx(2) = SIZE(u,2)-2
!
DO i = 0, mx(1)
  DO j = 0, mx(2)
!
    averesult(i,j) = 0.25_r8*(u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))
!
  END DO
END DO
!
END FUNCTION C2VAve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION V2CAve(v) RESULT(averesult)
USE Global
IMPLICIT NONE
!
! Vertex-to-center average. The variable v(0:mx(1),0:mx(2)) is a vertex-
! centered function. The result averesult(1:mx(1),1:mx(2)) is cell-centered.
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: v
REAL(KIND=r8), DIMENSION(1:SIZE(v,1)-1,1:SIZE(v,2)-1):: averesult
!
INTEGER:: i, j
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(v,1)-1; mx(2) = SIZE(v,2)-1
!
DO i = 1, mx(1)
  DO j = 1, mx(2)
!
    averesult(i,j) = 0.25_r8*(v(i,j)+v(i-1,j)+v(i,j-1)+v(i-1,j-1))
!
  END DO
END DO
!
END FUNCTION V2CAve
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION UDiv2D(f1,f2) RESULT(udivresult)
USE Global
IMPLICIT NONE
!
! UNDIVIDED divergence of the 2D flux function 
!
!  f = (f1(0:mx(1),1:mx(2)),f2(1:mx(1),0:mx(2))).  
!
! The result is stored in udivresult(1:mx(1),1:mx(2)).
!
REAL(KIND=r8), DIMENSION(0:,1:), INTENT(IN):: f1
REAL(KIND=r8), DIMENSION(1:,0:), INTENT(IN):: f2
REAL(KIND=r8), DIMENSION(1:SIZE(f2,1), &
                         1:SIZE(f1,2)):: udivresult
!
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(f2,1)
mx(2) = SIZE(f1,2)
!
udivresult(1:mx(1),1:mx(2)) = f1(1:mx(1)  ,1:mx(2)  ) &
                            - f1(0:mx(1)-1,1:mx(2)  ) &
                            + f2(1:mx(1)  ,1:mx(2)  ) &
                            - f2(1:mx(1)  ,0:mx(2)-1)
!
END FUNCTION UDiv2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION ULap2D(a) RESULT(ulapresult)
USE Global
IMPLICIT NONE
!
! 2D UNDIVIDED laplacian of a(0:mx(1)+1,0:mx(2)+1). The result is stored in
! ulapresult(1:mx(1),1:mx(2)).
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2, &
                         1:SIZE(a,2)-2):: ulapresult
!
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(a,1)-2
mx(2) = SIZE(a,2)-2
!
ulapresult(1:mx(1),1:mx(2)) =        a(2:mx(1)+1,1:mx(2)  ) &
                            +        a(0:mx(1)-1,1:mx(2)  ) &
                            +        a(1:mx(1)  ,2:mx(2)+1) &
                            +        a(1:mx(1)  ,0:mx(2)-1) &
                            - 4.0_r8*a(1:mx(1)  ,1:mx(2)  )
!
END FUNCTION ULap2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION ULap2DSpec(m1,m2,a) RESULT(ulapresult)
USE Global
IMPLICIT NONE
!
! 2D UNDIVIDED laplacian of a(0:mx(1)+1,0:mx(2)+1) with non-constant mobility 
! m. The result is stored in ulapresult(1:mx(1),1:mx(2)).
!
! m is edge-centered, a is cell centered.
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(0:,1:), INTENT(IN):: m1
REAL(KIND=r8), DIMENSION(1:,0:), INTENT(IN):: m2
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2,1:SIZE(a,2)-2):: ulapresult
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8), DIMENSION(0:SIZE(a,1)-2,1:SIZE(a,2)-2):: d1
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2,0:SIZE(a,2)-2):: d2
!
mx(1) = SIZE(a,1)-2
mx(2) = SIZE(a,2)-2
!
d1(0:mx(1),1:mx(2)) = m1(0:mx(1),1:mx(2))*(a(1:mx(1)+1,1:mx(2)  ) &
                    -                      a(0:mx(1)  ,1:mx(2)  ))
d2(1:mx(1),0:mx(2)) = m2(1:mx(1),0:mx(2))*(a(1:mx(1)  ,1:mx(2)+1) &
                    -                      a(1:mx(1)  ,0:mx(2)  ))
!
ulapresult(1:mx(1),1:mx(2)) = d1(1:mx(1),1:mx(2))-d1(0:mx(1)-1,1:mx(2)  ) &
                            + d2(1:mx(1),1:mx(2))-d2(1:mx(1)  ,0:mx(2)-1)
!
END FUNCTION ULap2DSpec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION ULap2DDiag(a) RESULT(ulapresult)
USE Global
IMPLICIT NONE
!
! 2D UNDIVIDED diagonal laplacian of a(0:mx(1)+1,0:mx(2)+1). The result is 
! stored in ulapresult(1:mx(1),1:mx(2)).
!
! a is cell-centered
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2,1:SIZE(a,2)-2):: ulapresult
!
INTEGER, DIMENSION(1:2):: mx
!
mx(1) = SIZE(a,1)-2
mx(2) = SIZE(a,2)-2
!
ulapresult(1:mx(1),1:mx(2)) = 0.5_r8*(a(2:mx(1)+1,2:mx(2)+1)+a(0:mx(1)-1,0:mx(2)-1) &
                            +         a(0:mx(1)-1,2:mx(2)+1)+a(2:mx(1)+1,0:mx(2)-1) &
                            -         4.0_r8*                a(1:mx(1)  ,1:mx(2)  ))
!
END FUNCTION ULap2DDiag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
FUNCTION ULap2DDiagSpec(r,a) RESULT(ulapresult)
USE Global
IMPLICIT NONE
!
! 2D UNDIVIDED diagonal laplacian of a(0:mx(1)+1,0:mx(2)+1) with non-constant 
! mobility r. The result is stored in ulapresult(1:mx(1),1:mx(2)).
!
! r(0:mx(1),0:mx(2)) is vertex-centered, a(0:mx(1)+1,0:mx(2)+1) is cell centered.
!
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: a
REAL(KIND=r8), DIMENSION(0:,0:), INTENT(IN):: r
REAL(KIND=r8), DIMENSION(1:SIZE(a,1)-2,1:SIZE(a,2)-2):: ulapresult
!
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8), DIMENSION(0:SIZE(r,1)-1,0:SIZE(r,2)-1):: d1, d2
!
mx(1) = SIZE(a,1)-2
mx(2) = SIZE(a,2)-2
!
d1(0:mx(1),0:mx(2)) = r(0:mx(1),0:mx(2))*(a(1:mx(1)+1,1:mx(2)+1) &
                    -                     a(0:mx(1)  ,0:mx(2)  ))
d2(0:mx(1),0:mx(2)) = r(0:mx(1),0:mx(2))*(a(0:mx(1)  ,1:mx(2)+1) &
                    -                     a(1:mx(1)+1,0:mx(2)  ))
!
ulapresult(1:mx(1),1:mx(2)) = 0.5_r8*(d1(1:mx(1)  ,1:mx(2))-d1(0:mx(1)-1,0:mx(2)-1) &
                            +         d2(0:mx(1)-1,1:mx(2))-d2(1:mx(1)  ,0:mx(2)-1))
!
END FUNCTION ULap2DDiagSpec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Sqr(p) RESULT(sqrresult)
USE Global
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: p
REAL(KIND=r8):: sqrresult
!
sqrresult = p*p
!
END FUNCTION Sqr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Cube(p) RESULT(cuberesult)
USE Global
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: p
REAL(KIND=r8):: cuberesult
!
cuberesult = p*p*p
!
END FUNCTION Cube
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Quad(p) RESULT(quadresult)
USE Global
USE ProblemDef
IMPLICIT NONE
!
REAL(KIND=r8), INTENT(IN):: p
REAL(KIND=r8):: quadresult
!
quadresult = p*p*p*p
!
END FUNCTION Quad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
ELEMENTAL FUNCTION Quin(p) RESULT(quinresult)
USE Global
USE ProblemDef
IMPLICIT NONE

REAL(KIND=r8), INTENT(IN):: p
REAL(KIND=r8):: quinresult
!
quinresult = p*p*p*p*p
!
END FUNCTION Quin
END MODULE Utilities
