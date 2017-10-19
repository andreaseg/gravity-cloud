PROGRAM data_generator
!
!   Purpose:
!     To generate a set of data representing initial-conditions of a dust cloud in zero gravity.
!     The particles will be normal-distributed around origo with speed of similiar distribution
!
!
IMPLICIT NONE

! Data dictionary
INTEGER :: i,j                    ! Loop index
INTEGER :: num                    ! Number of particles
INTEGER :: errorcode              ! Error code catcher
INTEGER :: iunit                  ! IO unit
CHARACTER (len=32) :: filename    ! Name of file to write to
REAL, ALLOCATABLE, DIMENSION(:) :: randomNumbers       ! Array to hold random numbers
REAL :: mass                      ! Particle information
REAL, DIMENSION(3) :: pos, vel    ! Particle information
REAL :: mass_mean, mass_std, pos_mean, pos_std, vel_mean, vel_std ! Particle statistical information
REAL :: position_gauss_mag, velocity_gauss_mag ! Magnitude for Gauss-distribution
REAL :: start, finish             ! For CPU-timing
CHARACTER (len=8) :: name

iunit = 9 ! IO socket

! Constants for generating particles
mass_mean = 1
mass_std = 0.15
pos_mean = 0
pos_std = 10
vel_mean = 1
vel_std = 0.3

! Set initial values
WRITE (*,*) 'Please enter number of particles for simulation'
READ (*,*) num
WRITE (*,1000) num
1000 FORMAT (1X,'Number of particles set to ',I7)

WRITE (*,*) 'Please enter file to write to.'
READ (*,*) filename
WRITE (*,1020) filename
1020 FORMAT (1X,'Filename set to ',A32)

! Start program timing
CALL cpu_time(start)

! Open file for writing
OPEN (UNIT=iunit, FILE=filename, STATUS='NEW',ACTION='WRITE',IOSTAT=errorcode)
openif: IF (errorcode == 0) THEN

  ! Writing header row to file
  WRITE (iunit,1015) 'NAME','MASS','X-POS','Y-POS','Z-POS','X-VEL','Y-VEL','Z-VEL'
  1015 FORMAT (1X,A8,7A12)

  ! Allocates the array to hold random numbers
  WRITE (*,1025) num
  1025 FORMAT (1X,'Attempting to allocate array of size ',I7)
  ALLOCATE (randomNumbers(9*num), STAT=errorcode)
  arrayif: IF (errorcode == 0) THEN

    ! Initializes the random number generator
    WRITE (*,*) 'Array allocated, generating random numbers'
    CALL RANDOM_NUMBER(randomNumbers)

    ! Writes the numbers to file
    WRITE (*,*) 'Random numbers generated, writing to file'
    loop: DO i = 1, num
      mass = normal_quantile(randomNumbers(i),mass_mean,mass_std)
      inner: DO j =1,3
        pos(j) = randomNumbers(2*j*num+i)
        vel(j) = randomNumbers((2*j-1)*num+i)
      END DO inner
      position_gauss_mag = randomNumbers(7*num+i)
      velocity_gauss_mag = randomNumbers(8*num+i)
      CALL rotatedRandom(pos,position_gauss_mag,pos_mean,pos_std)
      CALL rotatedRandom(vel,velocity_gauss_mag,vel_mean,vel_std)
      WRITE(name,'(Z8)') i
      WRITE (iunit,1100) name,mass,pos(1),pos(2),pos(3),vel(1),vel(2),vel(3)
      1100 FORMAT (1X,A8,7F12.4)

    END DO loop
  DEALLOCATE (randomNumbers,STAT=errorcode)
  END IF arrayif

ELSE openif
  WRITE (*,1200) errorcode
  1200 FORMAT (1X,'Error opening file: IOSTAT = ',I6'. Perhaps file already exits.')
END IF openif

! Close file
CLOSE (UNIT=iunit)

! End program timing
CALL cpu_time(finish)
WRITE (*,1300) (finish-start)
1300 FORMAT (1x,'(Time = ',F6.3,' seconds.)')

CONTAINS

REAL FUNCTION normal_quantile (num, mean, std)
IMPLICIT NONE
REAL, INTENT(IN) :: mean
REAL, INTENT(IN) :: std
REAL, INTENT(IN) :: num
normal_quantile = mean + std * sqrt(2.) * inverf(2*num-1)
END FUNCTION normal_quantile

REAL FUNCTION inverf (x)
IMPLICIT NONE
REAL, INTENT(IN) :: x
REAL :: pi
pi = 3.14159265359
inverf = sqrt(pi)*(.5*x + 1./24.*pi*x**3 + 7./960.*pi**2*x**5 + 127./80640.*pi**3*x**7) ! se http://mathworld.wolfram.com/InverseErf.html
END FUNCTION inverf

SUBROUTINE rotatedRandom (vec,mag,mean,std)
IMPLICIT NONE
REAL, DIMENSION(3), INTENT(INOUT) :: vec
REAL, INTENT(IN) :: mag
REAL, INTENT(IN) :: mean
REAL, INTENT(IN) :: std
REAL :: magnitude,length
REAL :: pi
pi = 3.14159265359
magnitude = normal_quantile(mag,mean,std)
length = sqrt( vec(1)**2 + vec(2)**2 + vec(3)**2 )

vec(1) = magnitude*vec(1)/length
vec(2) = magnitude*vec(2)/length
vec(3) = magnitude*vec(3)/length
END SUBROUTINE rotatedRandom

END PROGRAM data_generator
