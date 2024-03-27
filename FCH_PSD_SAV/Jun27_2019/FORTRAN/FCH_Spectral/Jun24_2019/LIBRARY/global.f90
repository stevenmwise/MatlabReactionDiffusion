MODULE Global
IMPLICIT NONE
!
LOGICAL:: debug, getenergy
INTEGER, PARAMETER:: r8 = SELECTED_REAL_KIND(15,307)
INTEGER:: ncgmaxits, secantmaxits
!
REAL(KIND=r8):: dt, ncgtol, secanttol, time
REAL(KIND=r8), PARAMETER:: pi = 3.1415926535897932_r8
REAL(KIND=r8), DIMENSION(1:2):: xlower, xupper
!
END MODULE Global
