PROGRAM NLCGDriver
USE Global
USE Utilities, ONLY: Ulap2D
USE NonLinearCGRoutines, ONLY: NonLinearCG
USE ProblemNCG, ONLY: SetProblem, Initialize, BoundaryConditions, TrueEnergy
IMPLICIT NONE
!
LOGICAL:: restart,convtest
INTEGER:: i, ierror, iter, j, k, maxtimeiterations, numprintouts, &
          printiterations, restartiteration
INTEGER, DIMENSION(1:2):: mx
REAL(KIND=r8):: energyvalue, dt_save, dx1, dx2, h, h2, start, finish, mass
REAL(KIND=r8), DIMENSION(1:2):: p
REAL(KIND=r8), DIMENSION(:,:), ALLOCATABLE:: f, kap, res, u, uo, uoo
!
CHARACTER(LEN=13):: file_name
CHARACTER(LEN=4):: string
!
NAMELIST/inputdata/ncgtol, mx, maxtimeiterations, dt, numprintouts, ncgmaxits, &
                   xlower, xupper, restart, restartiteration, &
                   debug, convtest, secanttol, secantmaxits
!
OPEN(UNIT=75,FILE='input.dat',STATUS='OLD',ACTION='READ',IOSTAT=ierror)
IF(ierror/=0) THEN
  PRINT *,'Error opening input file input.dat. Program stop.'
  STOP
END IF
READ(75,NML=inputdata)
CLOSE(75)
OPEN(UNIT=76,FILE='output.dat',STATUS='UNKNOWN',ACTION='WRITE', &
     FORM='FORMATTED', POSITION='APPEND')
WRITE(76,NML=inputdata)
CLOSE(76)
!
dt_save = dt
!
p(1:2) = xupper(1:2)-xlower(1:2)
dx1 = p(1)/REAL(mx(1),KIND=r8)
dx2 = p(2)/REAL(mx(2),KIND=r8)
!
IF(ABS(dx1-dx2)>1.0E-10_r8) THEN
  PRINT *, 'dx1 \= dx2.'
  STOP
ELSE
  h = dx1
END IF
!
h2 = h*h
IF (convtest) dt = 0.1_r8*h2
!
ALLOCATE(  f(1:mx(1)  ,1:mx(2)  ),res(1:mx(1)  ,1:mx(2)  ), &
           u(0:mx(1)+1,0:mx(2)+1),kap(0:mx(1)+1,0:mx(2)+1), &
          uo(0:mx(1)+1,0:mx(2)+1),uoo(0:mx(1)+1,0:mx(2)+1))
!
u = 0.0_r8
f = 0.0_r8
!
WRITE(string,'(i4)') mx(1)
DO i = 1, 4
  IF(string(i:i)==' ') string(i:i) = '0'
END DO
file_name = 'OUT/i'//string//'.dat'
!
IF(restart) THEN
  IF(restartiteration < numprintouts) THEN
    CALL ReadIn(restartiteration)
  ELSE
    PRINT *, 'restartiteration >= numprintouts'
    STOP
  END IF
ELSE
  restartiteration = 0
  CALL Initialize(u(0:mx(1)+1,0:mx(2)+1),h)
  CALL PrintOut(0)
END IF
CALL SetProblem
!
uo = u
uoo = u
!
OPEN(UNIT=199,FILE=file_name, STATUS='REPLACE', ACTION='WRITE')
printiterations = maxtimeiterations/numprintouts
!
! On the first step do Backward Euler. Need to adjust the timestep to do it:
dt = 1.5_r8*dt_save
!
DO i = restartiteration+1, numprintouts
  DO k = 1, printiterations
!
    time = REAL(k+(i-1)*printiterations,KIND=r8)*dt
!
    PRINT *, ' ' 
    PRINT *, 'Time = ', time, k+(i-1)*printiterations, maxtimeiterations
!
    uoo = uo
    uo = u
    u = 2.0_r8*uo-uoo   
!
    CALL CPU_TIME(start)
    CALL NonLinearCG(u,uo,uoo,f,res,h,iter,energyvalue)
    CALL CPU_TIME(finish)
!
    energyvalue = TrueEnergy(u(0:mx(1)+1,0:mx(2)+1), &
                             f(1:mx(1)  ,1:mx(2)  ),h)
!
    mass = h2*SUM(u(1:mx(1),1:mx(2)))
    WRITE(199,'(4(F25.12),I5)') time, energyvalue, finish-start, mass, iter
    
!    print '("CPU time of NCG = ",f6.3," seconds.")',finish-start
!
  dt = dt_save
!
  END DO
!
  kap(1:mx(1),1:mx(2)) = ULap2D(u(0:mx(1)+1,0:mx(2)+1))/h2
!  CALL PrintOut(mx(1))
  CALL PrintOut(i)
!
END DO
!
DEALLOCATE(f,kap,res,u,uo,uoo)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReadIn(k)
IMPLICIT NONE
!
INTEGER, INTENT(IN):: k
!
CHARACTER(LEN=13):: file_name
CHARACTER(LEN=4):: string
INTEGER:: i, j, idummy
INTEGER, DIMENSION(1:2):: mmx
REAL(KIND=r8):: dummy, hh, ttime, x, y
REAL(KIND=r8), DIMENSION(1:2):: xxlower, xxupper
!
WRITE(string,'(i4)') k
DO i = 1, 4
  IF(string(i:i)==' ') string(i:i) = '0'
END DO
file_name = 'OUT/m'//string//'.dat'
!
OPEN(UNIT=9, FILE=file_name, STATUS='UNKNOWN', ACTION='READ')
READ(9,'(F25.12)') dummy
READ(9,'(I8)') idummy
READ(9,'(F25.12)') hh
READ(9,'(F25.12)') dummy
READ(9,'(F25.12)') xxlower(1)
READ(9,'(F25.12)') xxlower(2)
READ(9,'(F25.12)') xxupper(1)
READ(9,'(F25.12)') xxupper(2)
READ(9,'(I8)') mmx(1)
READ(9,'(I8)') mmx(2)
!
IF(hh /= h) THEN 
  PRINT *, 'hh /= h'
  STOP
END IF
IF(ANY(xxlower /= xlower)) THEN 
  PRINT *, 'xxlower /= xlower'
  STOP
END IF
IF(ANY(xxupper /= xupper)) THEN 
  PRINT *, 'xxupper /= xupper'
  STOP
END IF
IF(ANY(mmx /= mx)) THEN 
  PRINT *, 'mmx /= mx'
  STOP
END IF
!
DO j = 1, mx(2)
  DO i = 1, mx(1)
    READ(9,'(4(f25.17,1x),f25.17)') dummy, dummy, u(i,j), dummy, dummy
  END DO
END DO
!
CLOSE(9)
!
CALL BoundaryConditions(u(0:mx(1)+1,0:mx(2)+1))
!
END SUBROUTINE ReadIn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE PrintOut(k)
IMPLICIT NONE
!
INTEGER, INTENT(IN):: k
!
CHARACTER(LEN=13):: file_name
CHARACTER(LEN=4):: string
INTEGER:: i, j
REAL(KIND=r8):: x, y
!
WRITE(string,'(i4)') k
DO i = 1, 4
  IF(string(i:i)==' ') string(i:i) = '0'
END DO
file_name = 'OUT/m'//string//'.dat'
!
OPEN(UNIT=9, FILE=file_name, STATUS='REPLACE', ACTION='WRITE')
WRITE(9,'(F25.12)') time
WRITE(9,'(I8)') 3
WRITE(9,'(F25.12)') h
WRITE(9,'(F25.12)') h
WRITE(9,'(F25.12)') xlower(1)
WRITE(9,'(F25.12)') xlower(2)
WRITE(9,'(F25.12)') xupper(1)
WRITE(9,'(F25.12)') xupper(2)
WRITE(9,'(I8)') mx(1)
WRITE(9,'(I8)') mx(2)
!
DO j = 1, mx(2)
  y = (REAL(j,KIND=r8)-0.5_r8)*h+xlower(2)
  DO i = 1, mx(1)
    x = (REAL(i,KIND=r8)-0.5_r8)*h+xlower(1)
    WRITE(9,'(4(f25.17,1x),f25.17)') x, y, u(i,j), kap(i,j), res(i,j)
  END DO
END DO
!
CLOSE(9)
!
END SUBROUTINE PrintOut
!
END PROGRAM NLCGDriver
