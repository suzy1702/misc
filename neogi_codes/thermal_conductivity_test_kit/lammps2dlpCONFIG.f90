program lammps2dlpCONFIG
!This program converts lammps file to CONFIG file, it writes out the CONFIG for every time step till the end timestep. 
!
!---Sanghamitra Neogi---
!---Mainz, April, 2014---
!-----------------------

  implicit none
  
  character(len=100) :: buffer,title,text
  character(len=100) :: fname,fn
  integer :: j, k, outunit=44,tempcount,tempnum,particletype,numtype,keytrj,perboundarykey
  integer :: n, natms, nstep, nstep0,numbins, direction, tempj, scaletag, imagetag
  double precision, dimension (3) :: lo, hi, cell
  double precision, allocatable, dimension (:,:) :: pos, vel, force
  double precision                               :: timestep
  integer, allocatable, dimension (:) :: num, atom_type
  integer, allocatable, dimension (:,:) :: image 
!  double precision, allocatable, dimension (:) :: masstype
  character(len=2), allocatable, dimension (:) :: nametype
  
  n=iargc()

  if(n<5)THEN
     write(*,*)"USAGE"
     write(*,*)"lammps2dlpCONFIG [inputfile timestep keytrj(0/1/2) perboundarykey scaletag(1/0,y/n) imagetag(1/0,y/n) particletype nametype]"
     STOP
  end if

  call getarg(1,buffer)
  read(buffer,*) fname
  call getarg(2,buffer)
  read(buffer,*)timestep
  call getarg(3,buffer)
  read(buffer,*) keytrj
  call getarg(4,buffer)
  read(buffer,*)perboundarykey
  call getarg(5,buffer)
  read(buffer,*)scaletag
  call getarg(6,buffer)
  read(buffer,*)imagetag
  if (n>6) then
     call getarg(7,buffer)
     read(buffer,*)particletype
  else
     particletype = 1
  endif
  allocate(nametype(particletype))
  do j = 1,particletype
     call getarg(7+j,buffer)
     read(buffer,*)nametype(j)
  end do

  open(7,file=fname,status='old')

  open(8,file='CONFIG')

  read(7,'(a)')title
  write(8,'(a)')'CONFIG'
  write(*,'(a)')'File header: ',title
  read(7,*)nstep
  nstep0 = nstep
  write(*,*)nstep-nstep0

  read(7, *)text
  read(7,*)natms
  print*,'total number of atoms in file',natms

  write(8,'(2i10, 2i10)')0, perboundarykey

  read(7,*)text
!  write(*,*)text
  read(7,*)lo(1), hi(1)
  read(7,*)lo(2), hi(2)
  read(7,*)lo(3), hi(3)
  
  cell = hi - lo
  print*, cell

!  write(8,*)'timestep',nstep-nstep0,natms,keytrj,perboundarykey,timestep*(nstep-nstep0)
!  write(8,*)lo(1), hi(1)
!  write(8,*)lo(2), hi(2)
!  write(8,*)lo(3), hi(3)

  read(7,*)text

!  write(*,*)text
  allocate(pos(3,natms),num(natms),atom_type(natms))
  allocate(vel(3,natms))
  if (imagetag.eq.1)then
     allocate(image(3,natms))
  endif

  do while (.true.)
!     tempcount = 0
     do j=1,natms
        if (keytrj.eq.0)then
           if (imagetag.eq.1) then
              read(7,*,end=100)tempnum, atom_type(tempnum), pos(:,tempnum), image(:,tempnum)
           else
              read(7,*,end=100)tempnum, atom_type(tempnum), pos(:,tempnum)
           endif
        elseif (keytrj.eq.1)then
           if (imagetag.eq.1) then
              read(7,*,end=100)tempnum, atom_type(tempnum), pos(:,tempnum), image(:,tempnum), vel(:,tempnum)
           else
              read(7,*,end=100)tempnum, atom_type(tempnum), pos(:,tempnum), vel(:,tempnum)
           endif
        endif
!        print*, 'before check',tempnum, atom_type(tempnum), pos(:,tempnum), vel(:,tempnum)
!        if (atom_type(tempnum).eq.particletype) then
!           print*, 'tempcount', tempcount
!           tempcount = tempcount + 1
!           num(tempcount) = tempnum
!           pos(:,num(tempcount)) = pos(:,tempnum)
!           vel(:,num(tempcount)) = vel(:,tempnum)
!        else
!           continue
!        endif
!        numtype = tempcount
!        write(*,*)j, num(j), atom_type(num(j)), pos(:,num(j)), vel(:,num(j))
     enddo

!!     print*, 'numtype',numtype
!
!!     do k=1,numtype
!       do j=1,natms
!!         if(num(j).eq.k) then
!!           tempj = num(j)
!!	   Output format to be used with c code and tfreq package
!!	   To match with the output of Quantum-ESPRESSO package 
!!          write(8,*)vel(:,tempj),tag
!!           write(8,*)tempj,num(j),num(k)
!           if(nstep.gt.0) then
!              write(8,'(a2,i10,2f12.6)')nametype(atom_type(j)),j,masstype(atom_type(j)),0.0
!!              write(*,*)nametype(atom_type(j)),j,masstype(atom_type(j)),0.0
!              write(8,'(1p,3e12.4)')pos(1,j)-lo(1),pos(2,j)-lo(2),pos(3,j)-lo(3)
!!              write(*,*)pos(:,j)
!              if(keytrj.eq.1)then
!                write(8,'(1p,3e12.4)')vel(:,j)
!              endif
!           endif 
!       end do
!!     end do
     read(7,'(a)',end=100)title
     read(7,*)nstep
     write(*,*)nstep-nstep0

     read(7, *)text
     read(7,*)natms

     read(7,*)text
     read(7,*)lo(1), hi(1)
     read(7,*)lo(2), hi(2)
     read(7,*)lo(3), hi(3)

!     if(nstep.gt.0) then
!        write(8,'(a8,4i10,f12.6)')'timestep',nstep-nstep0,natms,keytrj,perboundarykey,timestep*(nstep-nstep0)
!        write(8,'(3g12.4)')hi(1)-lo(1),0.0,0.0
!        write(8,'(3g12.4)')0.0,hi(2)-lo(2),0.0
!        write(8,'(3g12.4)')0.0,0.0,hi(3)-lo(3)
!     endif
     read(7,*)text
  end do

100 continue  

  write(8,'(3f20.12)')hi(1)-lo(1),0.0,0.0
  write(8,'(3f20.12)')0.0,hi(2)-lo(2),0.0
  write(8,'(3f20.12)')0.0,0.0,hi(3)-lo(3)

  do j=1,natms
!    if(num(j).eq.k) then
!    tempj = num(j)
!    Output format to be used with c code and tfreq package
!    To match with the output of Quantum-ESPRESSO package
!    write(8,*)vel(:,tempj),tag
!    write(8,*)tempj,num(j),num(k)
    if(nstep.gt.0) then
!      write(8,'(a2,i10,2f12.6)')nametype(atom_type(j)),j
     write(8,*)nametype(atom_type(j)),j
      if (scaletag.eq.0) then
!         write(8,'(3f20.12)')pos(1,j)-lo(1),pos(2,j)-lo(2), pos(3,j)-lo(3)
         write(8,'(3f20.12)')pos(1,j),pos(2,j), pos(3,j)
      else
         write(8,'(3f20.12)')pos(1,j)*cell(1)-lo(1),pos(2,j)*cell(2)-lo(2), pos(3,j)*cell(3)-lo(3)
      endif
!     write(*,*)pos(:,j)
!      if(keytrj.eq.1)then
!        write(8,'(1p,3e12.4)')vel(:,j)
!      endif
    endif
  end do


close(8)

end program
