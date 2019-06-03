 PROGRAM therm

 !This program computes thermal conductivity using Einstein's relation
 !
 !------Davide Donadio------
 !
 !-----------------------
 
 IMPLICIT NONE 
 INTEGER, PARAMETER                :: ndim = 3
 INTEGER                           :: i, j, k, l, ndata, ndata2, norm, ncorr, ncorrmax
 INTEGER                           :: iargc, nspace, npt, dimens =3
 REAL*8,ALLOCATABLE,DIMENSION(:,:) :: corr, errcorr
 REAL*8,ALLOCATABLE,DIMENSION(:,:) :: cprim, cdot, cdot2
 REAL*8,ALLOCATABLE,DIMENSION(:)   ::corrtmp, datareal, dataset, errtmp, datasetx, datasety, datasetz
 COMPLEX*16,ALLOCATABLE,DIMENSION(:) :: datatmp
 REAL*8                            :: tmax, volume, temp, dummy, dt, kT2vol,dt0,dt1
 REAL*8, DIMENSION(ndim)           :: average, conduct, sd, sd2, omega2, integ
 REAL*8                            :: cell(3,3), kTsqrt,thick,twopi, radius
 REAL*8, PARAMETER                 :: Kboltz=5.39054d-8 ! conversion from eV, A, ps 
 REAL*8, PARAMETER                 :: kjmol =  0.01036410d0 ! kjmol->eV
 CHARACTER*60                      :: buffer, fname, dum
 LOGICAL                           :: ltube = .false., l2D=.false., isotropic =.true.
 
 dt = .001
 twopi = 2.d0*ACOS(-1.)
 thick = 3.4d0

 if (iargc().lt.1) then 
   write(*,*)  "USAGE: conductance.x <ncorr>       (length of the fit in data points)"
   write(*,*)  "                     <temperature> (in Kelvin)            "
   write(*,*)  "                     <filename>    (deafult is HEATFLUX)     "
   write(*,*)  "                     <stride>      (sampling stride)          " 
   write(*,*)  "                     <dim>         (dimensionality, default 3D-iso)"
   write(*,*)  "                     <thickness>   (diameter/thickness, only for 1D or 2D systems)"
   write(*,*)  "NOTES:"
   write(*,*)  "The heatflux file is assumed in the format <time(ps)  jx  jy  jz (eV*A/ps)>."
   Write(*,*)  "The volume is computed by reading a CONFIG file (dlpoly format), and it assumes a triangular cell matrix."
   Write(*,*)  "1D systems are assumed to be periodic in the z direction (jx, jy are ignored)."
   Write(*,*)  "If 3D is specified, then it is assumed ANISOTROPIC."
   stop
 endif
 CALL GETARG(1,buffer)
 READ(buffer,*) ncorrmax
 CALL GETARG(2,buffer)
 READ(buffer,*) temp
 if (iargc().lt.3) then 
   fname = "HEATFLUX"
 else
   CALL GETARG(3,fname)
 endif
 if (iargc().lt.4) then
   nspace = 1000
 else 
   CALL GETARG(4,buffer)
   READ(buffer,*) nspace
 endif
 if (iargc().ge.5) then
   CALL GETARG(5,buffer)
   READ(buffer,*) dimens
 endif
 if (iargc().ge.6) then
   CALL GETARG(6,buffer)
   READ(buffer,*) thick
 endif

 if (dimens==1) ltube=.true.
 if (dimens==2) l2D=.true.
 if (dimens==3) isotropic = .false.
 
 write(*,*) "Stride: ", nspace
 write(*,*) "Temperature: ", temp

! read the volume of the system

 OPEN (unit=20,file='CONFIG',status='old')
 read(20,*)
 read(20,*)
 DO i = 1,3
   read(20,*) cell(:,i)
 ENDDO
 close(20)
 volume = cell(1,1)*cell(2,2)*cell(3,3)
 
!OPEN (unit=20,file='CONTROL',status='old')
!do while (.true.)
!   read(20,*,END=10) dum
!   if (index(dum,'tube').NE.0) ltube = .true.
!enddo
!close(20)
10 continue

 if (ltube) then
    radius = thick/2.d0
    write(*,*) "It's a 1D system! radius = ",radius
    volume = 0.5*twopi*radius**2*cell(3,3)
 endif
 if (l2D) then
   write(*,*) "It's a 2D system! thickness = ",thick
   volume = thick*cell(1,1)*cell(2,2)
 endif
 
 OPEN (unit=21,file=fname)

 READ(21,*) dt0
 READ(21,*) dt1
 dt = dt1-dt0
 REWIND(21)

 print*, "Data timestep: ", dt
 
 i = 0 
 DO WHILE(.true.)
   READ(21,*,END=100)  
   i = i+1
 ENDDO
 
 100 CONTINUE

 ndata = i
 print *,'ndata', ndata
 ALLOCATE(cprim (3,ndata))
 ALLOCATE(cdot  (3,ndata))
 ALLOCATE(cdot2 (3,ndata))
 ALLOCATE(corr(3,0:ncorrmax))
 ALLOCATE(errcorr(3,0:ncorrmax))
 ALLOCATE(corrtmp(0:ncorrmax))
 ALLOCATE(errtmp(0:ncorrmax))
 npt = ncorrmax/nspace+1
 ALLOCATE(dataset(npt))
 ALLOCATE(datasetx(npt))
 ALLOCATE(datasety(npt))
 ALLOCATE(datasetz(npt))

 tmax = dt*dble(ncorrmax)

 REWIND(21)
 
 DO i = 1, ndata
   READ(21,*) dummy, cdot(:,i)
   cdot(:,i) = cdot(:,i) !*kjmol
 ENDDO
 
 kT2vol =  Kboltz * Temp**2 * volume 
 kTsqrt = sqrt(kT2vol)!/2.d0
 cdot = cdot/kTsqrt

 do i = 1, 3
   if ((ltube .and. i.lt.3 ).or.(l2d .and. i==3)) then
     average(i) = 0.d0
     sd(i) = 0.d0
     cycle
   endif
   average(i) = SUM(cdot(i,:))/ndata
   sd2(i)     = SUM(cdot(i,:)*cdot(i,:))/ndata - average(i)**2
   cdot(i,:) = cdot(i,:) - average(i)
 enddo
! Compute the primitive and the derivative of the flux

 print *,'Computing the primitive'
 cprim(:,1) = 0.d0
 l=3
 DO i = 1, ndata - 1
   k=1 
   if(ltube) k = 3
   if(l2D) l = 2
   cprim(k:l,i+1) = cprim(k:l,i)  + (cdot(k:l,i) + cdot(k:l,i+1))*dt/2.d0
 ENDDO
 print *,'DONE'
  
 write(*,'(a16,3f15.3)')'Average         ', average
 write(*,'(a16,3f15.3)')'Av. Sq. fluct.  ', sd2
 write(*,'(a16, f15.3)')'Volume    (A^3) ', volume
 write(*,'(a16, f15.3)')'Temperature (K) ', temp
 write(*,'(a16, f15.3,i8)')'Max corr time   ', tmax, ncorrmax

! --- compute the correlation time for the integral

! CALL corrtime(cprim,ndim,ndata,average,sd,dt)

! --- compute Einstein diffusion coefficient

  ndata  = ndata - 1
  ALLOCATE(datareal(ndata))
  DO i = 1, ndim
    print *, '  JJ(' , i,')'
    if (ltube .and. i.lt.3 ) then 
       corr(i,:) = 0.d0
       cycle
    endif
    DO j = 1, ndata 
       datareal(j) = cprim(i,j)
    ENDDO
    call einstein(datareal,ndata,corrtmp,errtmp,ncorrmax,dt,nspace)
    corr(i,:) = corrtmp
    errcorr(i,:) = errtmp
  ENDDO

! open(unit =22, file='conduct_Einstein.dat')
  open(unit =12, file='diffusion.dat')
  j = 0
  DO i = 0, ncorrmax-1,nspace
    if (i==0) then 
      conduct =0.
    else  
      !conduct = corr(:,i)/2.d0/(dble(i)*dt)
      conduct = (corr(:,i+nspace) - corr(:,i))/2.d0/dt/nspace
    endif
!   write(22,'(f10.4,3g16.4)') dble(i)*dt, conduct(:)
!   write(23,'(f10.4,3g16.4)') dble(i)*dt, corr(:,i)
    j=j+1
    if (ltube) then 
        dataset(j) = corr(3,i)
        errtmp(j) = errcorr(3,i)
    else if (l2D) then
        datasetx(j) = corr(1,i)
        datasety(j) = corr(2,i)   
        errtmp(j) = (errcorr(1,i)+errcorr(2,i))/2.d0
    else if (.not.isotropic) then
        datasetx(j) = corr(1,i)
        datasety(j) = corr(2,i)
        datasetz(j) = corr(3,i)
        errtmp(j) = ( errcorr(1,i)+errcorr(2,i)+errcorr(3,i) )/3.d0
    else 
        dataset(j) = (corr(1,i)+corr(2,i)+corr(3,i)) /3.d0
        errtmp(j) = sum(errcorr(:,i))/3.d0
    endif
  ENDDO

  
  npt = npt -1
  dt = dt*DBLE(nspace)
  if (l2D) then ! we have to separate x and y
    call MCfit(npt,datasetx,dt,errtmp)
    call MCfit(npt,datasety,dt,errtmp)
  else if (.not.isotropic) then
    call MCfit(npt,datasetx,dt,errtmp)
    call MCfit(npt,datasety,dt,errtmp)
    call MCfit(npt,datasetz,dt,errtmp)
  else
    call MCfit(npt,dataset,dt,errtmp)
  endif
  
 END PROGRAM therm
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE einstein(cdata,npt,corr,errc,nc,dt,nspace)
  IMPLICIT NONE

! in/out arguments

  INTEGER    :: nspace, nc !max corr points
  INTEGER    :: npt        !size of the dataset
  REAL*8     :: dt
  REAL*8     :: corr(0:nc), errc(0:nc)
  REAL*8     :: cdata(npt)

! local parameters

  INTEGER    :: i, ncorr, norm

  DO ncorr = 0, nc, nspace
    norm = (npt-ncorr)/nspace
    DO i = 1, npt - ncorr, nspace
      corr(ncorr) = corr(ncorr) + (cdata(i)-cdata(i+ncorr))**2
      errc(ncorr) = errc(ncorr) + (cdata(i)-cdata(i+ncorr))**4
    ENDDO
    corr(ncorr) = corr(ncorr)/ norm
    errc(ncorr) = sqrt( (errc(ncorr)/norm - corr(ncorr)**2)/norm )
    IF(MOD(ncorr,nc/10).EQ.0)print'(i3," % done")',100*(ncorr)/(nc)
  ENDDO

  RETURN
  END SUBROUTINE 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE corrtime(cdata,ndim,ndata,average,sd2,dt)
  IMPLICIT NONE

! in/out arguments

  INTEGER    :: ndim, ndata  !max corr points
  REAL*8     :: cdata(ndim,ndata),sd2(ndim),dt,average(ndim)

! local variables

  INTEGER, PARAMETER :: nptbmin = 1000
  INTEGER    :: nblocks, nptblock, nbloops
  INTEGER    :: i, j, idm, n_ini, n_fin

  REAL*8     :: taublock, sd2blk
  REAL*8, ALLOCATABLE, DIMENSION(:)   :: aveblk
  REAL*8, ALLOCATABLE, DIMENSION(:,:) :: taus

  OPEN (unit=24, file='tau.dat')
  nbloops = ndata/nptbmin
  ALLOCATE(taus(ndim,nbloops))

  i=0
  WRITE(*,*) 'Computing the statistical inefficiency'

  DO
    i = i+1
    nptblock = i*nptbmin
    taublock = DBLE(nptblock)*dt
    nblocks = ndata/nptblock-1
    IF (nblocks<10) exit
    ALLOCATE(aveblk(nblocks))
    DO idm = 1, ndim

       DO j = 1, nblocks
         n_ini= (j-1)*nptblock+1
         n_fin= j*nptblock
         aveblk(j) = SUM(cdata(idm,n_ini:n_fin))/nptblock
       ENDDO
       sd2blk = SUM( (aveblk - average(idm))**2 )/nblocks
       taus(idm,i) = taublock * sd2blk/sd2(idm)

    ENDDO !dimension loop
    DEALLOCATE(aveblk)
    WRITE(24,'(f10.3,4g16.6,i6)') taublock, taus(:,i),sum(taus(:,i))/3,nblocks 

  ENDDO ! do while loop
  close(24)

  return
  END SUBROUTINE corrtime
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE MCfit(npt, dataset, dtp, errc)

  INTEGER npt
  REAL*8 dataset(npt), dtp, errc(npt)

  INTEGER              :: i, j, k
  INTEGER              :: nstep = 4, npar = 2 , nseek =50000
  REAL*8, DIMENSION(2) :: par, parmax, parmin, deltapar, parnew
  REAL*8               :: chisq, chinew, rand, sigmapar, deltax,dt,sumdt

  chisq = 1.0e30
  chinew= 1.0e30
  parmax(1) = 1.d5
  parmin(1) = 0.d0
  parmax(2) = 1.d4
  parmin(2) = 0.d0

  par(2)=dataset(npt)/2./npt/dtp
  print*,'guess k=', par(2) 
  par(1)=1.d2 ! tau

  deltapar(1) = 1.d2
  deltapar(2) = 1.d2

  call compchi(npar,par,chisq,npt,dataset,dtp,errc)

  do j =1, nstep
    do i = 1,nseek

!     generate new parameters

      do k = 1,npar
        do
          call random_number(rand)
          parnew(k) = par(k) + (0.5 - rand)*deltapar(k)
          if(parnew(k)<parmax(k) .and. parnew(k)>parmin(k)) exit
        enddo
      enddo
      call compchi(npar,parnew,chinew,npt,dataset,dtp,errc)
      if(chinew.lt.chisq) then
         do n = 1,npar
            par(n) = parnew(n)
         enddo
         chisq = chinew
!        print*,'chi',chinew
      endif
    enddo
    do n = 1,npar
       deltapar(n) = deltapar(n)/5.d0
    enddo
  enddo
! compute error on kappa
 
  deltax = 0.d0
  sumdt = 0.d0
  do i =1, npt
    dt = dble(i-1)*dtp  
    deltax = deltax+ dble(npt)*dt**2 
    sumdt = sumdt+dt
  enddo
  deltax = deltax - sumdt**2
  sigmapar = sqrt(dble(npt)*chisq/deltax)
  write(*,'("kappa= ",f12.2," tau= ",f12.2," chisq=",2g12.4)') par(2), par(1), chisq, sigmapar
  do i = 2, npt
    write(12,'(g14.4,4g16.8)') (i-1)*dtp, (dataset(i)-dataset(i-1))/2.d0/dtp, dataset(i), & 
                               errc(i), brown((i-1.)*dtp, par(1), par(2))
                              
  enddo
  RETURN

  END SUBROUTINE MCfit
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  SUBROUTINE compchi(npar,par,chisq,npt,dataset,dtp,errc)

  INTEGER npt
  REAL*8 dataset(npt), errc(npt), dtp

  INTEGER                 :: npar
  REAL*8, DIMENSION(npar) :: par
  REAL*8                  :: chisq

  INTEGER :: i, j
  REAL*8                  :: tau, kappa, time

  tau   = par(1)
  kappa = par(2)

  chisq = 0.d0
  do i = 1, npt

    time = DBLE(i-1)*dtp

    if(abs(dataset(i))>1.d-8)   & 
      chisq = chisq + ((brown(time, tau, kappa) - dataset(i))**2+errc(i)**2)!/dataset(i)**2
  enddo
  chisq  = chisq/ (npt-2)
  RETURN
  END SUBROUTINE compchi
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  DOUBLE PRECISION FUNCTION brown(t, tau, kappa)

  REAL*8 t, tau, kappa

  brown = 2.d0*kappa*( t + tau*(exp(-t/tau)-1.d0) )

  END

