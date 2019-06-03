PROGRAM heatflux_format
!This program rewrites the log.lammps file by truncating the end columns so that it can be used with einst.x program
!
!---Sanghamitra Neogi---
!---Mainz, June, 2014---
!-----------------------

  IMPLICIT NONE
  
  CHARACTER(LEN=100) :: buffer,title
  CHARACTER(LEN=100) :: fname
  INTEGER :: time, numcolumn, time0, initime
  INTEGER :: n, p
  DOUBLE PRECISION :: heat1, heat2, heat3, timestep
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: extracolumn
  CHARACTER(LEN=8),DIMENSION (:), ALLOCATABLE    :: name
  
  p=iargc()

  IF(p<3)THEN
     WRITE(*,*)"USAGE"
     WRITE(*,*)"heatflux_format [inputfile numextracolumns timestep initime]"
     STOP
  END IF

  CALL getarg(1,buffer)
  READ(buffer,*) fname
  CALL getarg(2,buffer)
  READ(buffer,*) numcolumn
  CALL getarg(3,buffer)
  READ(buffer,*) timestep
  if (p>3)then
    CALL getarg(4,buffer)
    READ(buffer,*) initime
  endif

  ALLOCATE(extracolumn(numcolumn))

  OPEN(7,file=fname,status='old')

  OPEN(8,file='HEATFLUX_formatted')

  n=0
  DO WHILE (.TRUE.)
     READ(7,*,END=100)time, heat1, heat2, heat3, extracolumn(1:numcolumn)
     if (n.eq.0) then
       if (p>3)then
          time0 = initime
       else
          time0 = time
       endif
       print*, time0
     endif
     if (mod(n,100000).eq.0) then
       print*, (time-time0)*timestep, heat1, heat2, heat3, extracolumn(1:numcolumn)
     endif
     WRITE(8,'(4f20.12)')(time-time0)*timestep, heat1, heat2, heat3
     n=n+1
  END DO
 
100 CONTINUE

print*,'total number of lines in file',n

END PROGRAM
