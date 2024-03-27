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
