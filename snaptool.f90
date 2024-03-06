! Snaptool - an utility to visualize cosmological simulation snapshots
! 
! Copyright (C) 2024 Sergey Pilipenko
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.



program snaptool
! modules and libraries
 use gif_util
 use octtree

! interpreter global vars
character(len=800) :: arg
integer:: infile

! command-specific global vars

! set map
logical:: was_set_map = .false.
integer:: isize = 500, jsize = 500, iaxis = 3, ia=1, ja=2, ka=3
real, allocatable:: map(:,:)

! read pos
logical:: was_read_pos = .false., plusread = .false.
real, allocatable:: pos(:,:)
integer:: Ntot, numfiles, Npos, Nall_arr(6), gformat, Npart(6)
byte, allocatable:: ptype(:)
real:: Box, tsnap
integer*8:: Nall
real*8:: massarr(6)
character(len=4):: secname

! set read every
logical:: was_set_read_every = .false.
integer*8:: every = 1

! map

! set center
logical:: was_set_center = .false.
real:: center(3), oldcenter(3), shift(3)

! set size
logical:: was_set_size = .false.
real:: regionsize

! set cut
logical:: set_cut = .false.

! set N pos max
integer:: Nposmax = 2000000000

! set format
character(len=3):: outformat

! set outfile
character(len=800):: basename

! add suffix z
logical:: suffix_z = .false.

! add suffix num
logical:: suffix_num = .false.
integer:: numout = 0

! add suffix axes
logical:: suffix_axes = .false.

! set minimum brightness
integer:: mincol = 30

! select sphere
logical,allocatable:: is_selected(:)

! set stats file
logical:: set_stats = .false.
character(len=800):: statsfile
integer:: fstat = 19

! dens
real, allocatable:: dens(:)
type(tree) octree

! vel
real, allocatable:: vel(:,:)
logical:: readvel = .false.
real*8::velrms

! set drange
real:: drange = 5

! set ntree
integer::ntree = 30

! set maxalloc
logical:: maxalloc = .false.

! rot
real:: rot(3,3)
logical:: was_rot = .false.

! AHF
character(len=800):: ahf_file
integer,parameter:: N_halo_max = 100000
integer:: nh = 0
real:: halo_pos(3,N_halo_max), halo_r(N_halo_max)

! palette
integer:: palette = 1

! grid
real:: grid_step = 1.0
logical:: do_grid = .false.

! read ascii
integer:: Nmax_ascii = 1000000

! evolve center
logical::was_set_center_evolution = .false.
real::startcenter(3), endcenter(3)
real:: tstart, tend

! main

rot = 0.0
rot(1,1) = 1.0
rot(2,2) = 1.0
rot(3,3) = 1.0

shift = 0.0

if(iargc().lt.1) then
  print*,'Usage: snaptool in-file'
  print*,'or: snapplot --help'
  return
endif

call getarg(1,arg)

if(arg.eq.'--help') then
  call help()
  return
endif

basename = 'snapplot'
infile = 20
call interpret(arg)


contains

recursive subroutine interpret(fin)
  character(len=800):: fin, finclude
  character(len=32):: command
  open(infile,file=trim(fin),status='old')
5 read(infile,'(A)',end=10,err=10) command
  command = trim(command)
  
  if(command.eq.'set map') then
    call setmap()
    goto 5
  endif
  if(command.eq.'read pos') then
    plusread = .false.
    call readpos()
    goto 5
  endif
  if(command.eq.'+read pos') then
    plusread = .true.
    call readpos()
    goto 5
  endif
  if(command.eq.'end') then
    return
  endif
  if(command.eq.'set read every') then
    call set_read_every()
    goto 5
  endif
  if(command.eq.'map') then
    call fillmap()
    goto 5
  endif
  if(command.eq.'set center') then
    read(infile,*) center
    was_set_center = .true.
    goto 5
  endif
  if(command.eq.'move center') then
    call setcenter()
    goto 5
  endif
  if(command.eq.'set size') then
    call setsize()
    goto 5
  endif
  if(command.eq.'set cut') then
    set_cut = .true.
    goto 5
  endif
  if(command.eq.'unset cut') then
    set_cut = .false.
    goto 5
  endif
  if(command.eq.'set N pos max') then
    read(infile,*) Nposmax
    goto 5
  endif
  if(command.eq.'add suffix z') then
    suffix_z = .true.
    goto 5
  endif
  if(command.eq.'add suffix num') then
    suffix_num = .true.
    goto 5
  endif
  if(command.eq.'add suffix axes') then
    suffix_axes = .true.
    goto 5
  endif
  if(command.eq.'set outfile') then
    read(infile,'(A)') basename
    goto 5
  endif
  if(command.eq.'plot') then
    call plot()
    goto 5
  endif
  if(command.eq.'include') then
    read(infile,'(A)') finclude
    infile = infile+1
    call interpret(finclude)
    infile = infile-1
    goto 5
  endif
  if(command.eq.'select sphere') then
    call select_sphere()
    goto 5
  endif
  if(command.eq.'select type') then
    call select_type()
    goto 5
  endif
  if(command.eq.'unselect') then
    call unselect()
    goto 5
  endif
  if(command.eq.'rotate') then
    call rotate()
    goto 5
  endif
  if(command.eq.'set minimum brightness') then
    read(infile,*) mincol
    goto 5
  endif
  if(command.eq.'set stats file') then
    read(infile,'(A)') statsfile
    open(fstat,file=trim(statsfile))
    set_stats = .true.
    goto 5
  endif
  if(command.eq.'make density') then
    call makedens()
    goto 5
  endif
  if(command.eq.'make velocity') then
    call makevel()
    goto 5
  endif
  if(command.eq.'make energy density') then
    call makeenergydens()
    goto 5
  endif
  if(command.eq.'make phase density') then
    call makeentropy()
    goto 5
  endif
  if(command.eq.'dens') then
    call densmap()
    goto 5
  endif
  if(command.eq.'set vel') then
    readvel = .true.
    goto 5
  endif
  if(command.eq.'set drange') then
    read(infile,*) drange
    goto 5
  endif
  if(command.eq.'set ntree') then
    read(infile,*) ntree
    goto 5
  endif
  if(command.eq.'set maxalloc') then
    maxalloc = .true.
    goto 5
  endif
  if(command.eq.'AHF') then
    read(infile,'(A)') ahf_file
    call ahf()
    goto 5
  endif
  if(command.eq.'set palette') then
    read(infile,*) palette
    goto 5
  endif
  if(command.eq.'set grid') then
    read(infile,*) grid_step
    do_grid = .true.
    goto 5
  endif
   if(command.eq.'read ascii') then
    plusread = .false.
    call readascii()
    goto 5
  endif
  if(command.eq.'set N ascii max') then
    read(infile,*) Nmax_ascii
    goto 5
  endif
  if(command.eq.'load ramses map') then
    call load_ramses()
    goto 5
  endif
  if(command.eq.'print dist') then
    call print_dist()
    goto 5
  endif
  if(command.eq.'print mass') then
    call print_mass()
    goto 5
  endif
  if(command.eq.'evolve center') then
    call setcenterevolution()
    goto 5
  endif
  if(command.eq.''.or.command(1:1).eq.'#') then
    goto 5
  endif
  print*,'Unknown command: ',command
10 close(infile)
end subroutine interpret

! commands

subroutine setmap()
  character(len=2):: caxis
  read(infile,*) isize, jsize, caxis
  if(caxis.eq.'xy') iaxis = 3
  if(caxis.eq.'xz') iaxis = 2
  if(caxis.eq.'yz') iaxis = 1
  if(allocated(map)) deallocate(map)
  allocate(map(isize,jsize))
  map = 0
  if(iaxis.eq.1) then
   ia=2
   ja=3
   ka=1
  endif
  if(iaxis.eq.2) then
   ia=1
   ja=3
   ka=2
  endif
  if(iaxis.eq.3) then
   ia=1
   ja=2
   ka=3
  endif
  was_set_map = .true.
end subroutine setmap

subroutine readpos()
  character(len=800):: fsnap
  real, allocatable:: buf(:,:), vbuf(:,:)
  print*, 'Reading...'
  numfiles = 1
  nfilesread = 0
  if(plusread) then
	ipos = Npos
  else
    ipos = 1
  endif
  if(allocated(pos).and..not.plusread) then
    deallocate(pos)
    deallocate(is_selected)
    deallocate(ptype)
    if(allocated(vel)) deallocate(vel)
  endif
  do while(numfiles.gt.nfilesread)
    read(infile,'(A)') fsnap
    call read_header(fsnap)
    Npos = min(Nall/every+every+numfiles, Nposmax)
!    print*, Npos
    if(maxalloc) Npos = Nposmax
    if(nfilesread.eq.0.and..not.plusread) then
      allocate(pos(3,Npos))
      allocate(is_selected(Npos))
      allocate(ptype(Npos))
      is_selected = .true.
      if(readvel) allocate(vel(3,Npos))
    endif
    allocate(buf(3,Ntot))
    if(gformat.eq.2) read(2) secname
    read(2) buf
    if(readvel) then
		allocate(vbuf(3,Ntot))
		if(gformat.eq.2) read(2) secname
		read(2) vbuf
	endif
    close(2)
    call set_types(ipos)
    
    if(.not.was_set_size) regionsize = Box
    if(.not.was_set_center) center(1:3)=Box/2
    if(.not.plusread) oldcenter(1:3)=Box/2
    do k=1,3
      do i=1,Ntot
        buf(k,i)=buf(k,i)-center(k)
        if(buf(k,i).gt.Box/2) buf(k,i)=buf(k,i)-Box
        if(buf(k,i).lt.-Box/2) buf(k,i)=buf(k,i)+Box
      enddo
    enddo
    
    do i=1,Ntot,every
      if(.not.set_cut) then
        pos(1:3,ipos) = buf(1:3,i)
        if(readvel) vel(1:3,ipos) = vbuf(1:3,i)
        ipos = ipos+1
      else
        if(abs(buf(1,i)).lt.regionsize/2.and. &
           abs(buf(2,i)).lt.regionsize/2.and. &
           abs(buf(3,i)).lt.regionsize/2) then
             pos(1:3,ipos) = buf(1:3,i)
             if(readvel) vel(1:3,ipos) = vbuf(1:3,i)
			 ipos = ipos+1
        endif
      endif
      if(ipos.gt.Npos) then
        print*,'Warning: increase N pos max!!!'
        exit
      endif
    enddo
    deallocate(buf)
    if(readvel) deallocate(vbuf)
    if(allocated(map)) map = 0
    nfilesread = nfilesread+1
  enddo
  print*,'read snapshot at z = ', 1/tsnap - 1
  Npos = ipos-1
  print*,'total ', Npos, ' particles'	
  if(readvel) velrms = sqrt(sum(vel(1,:)*vel(1,:))+sum(vel(2,:)*vel(2,:))+sum(vel(3,:)*vel(3,:)))
  if(was_set_center_evolution) call evolvecenter()
  was_read_pos = .true.
end subroutine readpos

subroutine set_types(ipos)
  integer:: npart_sum(6)
  npart_sum(1) = Npart(1)
  do k=2,6
    npart_sum(k) = Npart_sum(k-1) + Npart(k)
  enddo
  do i=1,Ntot, every
    if(i.le.Npart_sum(1)) ptype(i+ipos-1) = 0
    do k=1,5
      if(i.gt.Npart_sum(k).and.i.le.Npart_sum(k+1)) ptype(i+ipos-1) = k
    enddo
  enddo
end subroutine set_types


subroutine readascii
  character(len=800):: fsnap
  read(infile,'(A)') fsnap
  if(allocated(pos)) then
	deallocate(pos)
	deallocate(is_selected)
  endif
  allocate(pos(3,Nmax_ascii))
  allocate(is_selected(Nmax_ascii))
  is_selected = .true.
  open(2,file=fsnap,status='old')
  Npos = 1
  do while(.true.)
    read(2,*,end=12,err=12) (pos(k,Npos),k=1,3)
    Npos = Npos + 1
  enddo
 12 close(2)
  Npos = Npos - 1
  was_read_pos = .true.
  if(.not.was_set_center) center(1:3)=0
  if(.not.plusread) oldcenter(1:3)=0
end subroutine readascii

subroutine set_read_every()
  was_set_read_every = .true.
  read(infile,*) every
end subroutine set_read_every

subroutine fillmap()
logical:: zfit
mijhalf = (isize-jsize)/2
!$OMP PARALLEL PRIVATE(i,mi,mj,zfit,bx,by,tx,ty,mip,mjp), REDUCTION(+ : map)
!$OMP DO
do i=1,Npos
  if(is_selected(i)) then
  mi=int(normpos(ia,pos(:,i))*jsize+1+mijhalf)
  mj=int(normpos(ja,pos(:,i))*jsize+1)
  zfit=(normpos(ka,pos(:,i)).ge.0).and.(normpos(ka,pos(:,i)).le.1)
  if(mi.ge.1.and.mj.ge.1.and.mi.le.isize.and.mj.le.jsize.and.zfit) then
     bx = (normpos(ia,pos(:,i))*jsize+1+mijhalf) - mi
     by = (normpos(ja,pos(:,i))*jsize+1) - mj
     tx = 1 - bx
     ty = 1 - by
     mip = min(mi+1,isize)
     mjp = min(mj+1,jsize)
     map(mi,mj)=map(mi,mj)+tx*ty
     map(mip,mj)=map(mip,mj)+bx*ty
     map(mi,mjp)=map(mi,mjp)+tx*by
     map(mip,mjp)=map(mip,mjp)+bx*by
  endif
  endif
enddo
!$OMP END PARALLEL
end subroutine fillmap

subroutine load_ramses()
  character(len=800):: fmap
  read(infile,*) fmap
  open(9,file=trim(fmap),status='old')
  do i = 1,isize
    do j = 1,jsize
      read(9,*) tx,ty,map(i,j)
    enddo
  enddo
  close(9)
end subroutine load_ramses


subroutine makevel()
if(allocated(dens)) deallocate(dens)
allocate(dens(Npos))
do i=1,Npos
	dens(i) = sum(vel(1:3,i)*vel(1:3,i))
enddo
end subroutine makevel

subroutine makeenergydens()
call makedens()
do i=1,Npos
	dens(i) = dens(i)*sum(vel(1:3,i)*vel(1:3,i))
enddo
end subroutine makeenergydens

subroutine makedens()
real:: boxtree(3), centtree(3)
! add something to clear tree
if(allocated(dens)) deallocate(dens)
if(octree%np.ne.0) call tree_delete(octree)
boxtree = regionsize
centtree = 0
octree = tree_new(boxtree,centtree,Npos,ntree)
print*, 'Making tree from ',Npos,' points...'
call tree_construct(octree, pos)
print*, 'done (',octree%nlevels,' levels)'
allocate(dens(Npos))
print*, 'Calculating density...'
!$OMP PARALLEL DO
do i=1,Npos
	dens(i)=tree_dens(pos(1,i),octree)
enddo
!$OMP END PARALLEL DO
print*, 'done'
end subroutine makedens


subroutine makeentropy()
real:: boxtree(3), centtree(3), v(3), vv(3)
type(node), pointer:: p
if(allocated(dens)) deallocate(dens)
if(octree%np.ne.0) call tree_delete(octree)
boxtree = regionsize
centtree = 0
octree = tree_new(boxtree,centtree,Npos,ntree)
print*, 'Making tree from ',Npos,' points...'
call tree_construct(octree, pos)
call tree_makelist(octree, pos)
print*, 'done (',octree%nlevels,' levels)'
allocate(dens(Npos))
print*, 'Calculating phase density...'
dens = -1
do i=1,Npos
	if(dens(i).lt.0) then
	p => tree_find(pos(1,i),octree,lev)
	sv=0
	v=0
	vv=0
	if(p%np.lt.3) then
		sv = velrms
		dens(i)=tree_dens(pos(1,i),octree)/sv/sv/sv
	else
		ind = p%start
		do while(ind.ne.-1)
			v=v+vel(1:3,ind)
			vv = vv+vel(1:3,ind)*vel(1:3,ind)
			ind = octree%list(ind)	
		enddo
		vv = vv/p%np
		v = (v/p%np)*(v/p%np)
		sv = sqrt(sum(vv-v))
		ind = p%start
		d=tree_dens(pos(1,i),octree)/sv/sv/sv
		do while(ind.ne.-1)
			dens(ind)=d
			ind = octree%list(ind)
		enddo
	endif
	endif
enddo
print*, 'done'
end subroutine makeentropy



subroutine densmap()
logical:: zfit

do i=1,Npos
  if(is_selected(i)) then
  mi=int(normpos(ia,pos(:,i))*isize+1)
  mj=int(normpos(ja,pos(:,i))*jsize+1)
  zfit=(normpos(ka,pos(:,i)).ge.0).and.(normpos(ka,pos(:,i)).le.1)
  if(mi.ge.1.and.mj.ge.1.and.mi.le.isize.and.mj.le.jsize.and.zfit) then
     bx = (normpos(ia,pos(:,i))*isize+1) - mi
     by = (normpos(ja,pos(:,i))*jsize+1) - mj
     tx = 1 - bx
     ty = 1 - by
     mip = min(mi+1,isize)
     mjp = min(mj+1,jsize)
     map(mi,mj)= max(map(mi,mj),tx*ty*dens(i))
     map(mip,mj)= max(map(mip,mj),bx*ty*dens(i))
     map(mi,mjp)= max(map(mi,mjp),tx*by*dens(i))
     map(mip,mjp)= max(map(mip,mjp),bx*by*dens(i))
  endif
  endif
enddo
end subroutine densmap


real function normpos(iax,oldpos)
  real:: newpos(3), oldpos(3)
  if(was_rot) then
    do k=1,3
      newpos(k) = sum(oldpos(1:3)*rot(k,1:3))
    enddo
  else
    newpos = oldpos(1:3)
  endif
  normpos = (newpos(iax))/regionsize+0.5
end function normpos

subroutine setcenter()
  real:: newcenter(3)
  read(infile,*) newcenter
  call setcenter_core(newcenter)
end subroutine setcenter

subroutine setcenterevolution()
  read(infile,*) startcenter, zstart
  read(infile,*) endcenter, zend
  tstart = 1/(1+zstart)
  tend = 1/(1+zend)
  was_set_center_evolution = .true.
end subroutine setcenterevolution

subroutine evolvecenter()
  real:: newcenter(3)
  newcenter = startcenter + (endcenter - startcenter)/(tend - tstart) * (tsnap - tstart)
  call setcenter_core(newcenter)
  print*,'center: ',newcenter
end subroutine evolvecenter

subroutine setcenter_core(newcenter)
  real:: newcenter(3)
  oldcenter = center
  center = newcenter
  was_set_center = .true.
  if(allocated(pos).and.Npos.gt.0) then
!$OMP PARALLEL PRIVATE(k,i)
!$OMP DO
  do i=1,Npos
    do k=1,3
      pos(k,i)=pos(k,i)-center(k)+oldcenter(k)
      if(pos(k,i).lt.-Box/2) pos(k,i)=pos(k,i)+Box
      if(pos(k,i).gt.Box/2) pos(k,i)=pos(k,i)-Box
    enddo
  enddo
!$OMP END PARALLEL
  endif
  shift = shift-center+oldcenter
end subroutine setcenter_core

subroutine setsize()
  read(infile,*) regionsize
  was_set_size = .true.
end subroutine setsize

subroutine plot()
integer, parameter :: ncol=256
integer 	   :: pix(isize,jsize), colmap(3,0:ncol-1)
character(len=820)::fout
character(len=6):: suf
logical:: zfit

! setting output filename
fout = trim(basename)
ilast = scan(fout,' ')-1
if(suffix_z) then
  z = 1/tsnap - 1
  if(z.lt.10) write(suf,'(F6.3)') z
  if(z.gt.10.and.z.lt.100) write(suf,'(F6.2)') z
  if(z.gt.100) write(suf,'(F6.1)') z
  fout=fout(1:ilast)//'_z'//suf(2:6)
  ilast = ilast+7
endif
if(suffix_axes) then
  if(iaxis.eq.1) fout=fout(1:ilast)//'_yz'
  if(iaxis.eq.2) fout=fout(1:ilast)//'_xz'
  if(iaxis.eq.3) fout=fout(1:ilast)//'_xy'
  ilast = ilast+3
endif
if(suffix_num) then
  write(suf,'(I6)') 100000+ numout
  fout=fout(1:ilast)//'_'//suf(2:6)
  ilast = ilast+6
endif

amax = 0
amin = 1e20
do mi = 1, isize
do mj = 1,jsize
  amax = max(amax,map(mi,mj))
  if(map(mi,mj).gt.0) amin = min(amin,map(mi,mj))
enddo
enddo

sm = 0
sm2 = 0
nmap = 0
do mi = 1, isize
do mj = 1,jsize
  if(map(mi,mj).gt.0) then
	  sm = sm + log(map(mi,mj))
	  sm2 = sm2 + log(map(mi,mj))*log(map(mi,mj))
	  nmap = nmap + 1
  endif
enddo
enddo
sm = sm/nmap
sm2 = sqrt(sm2/nmap - sm*sm)
amin = max(log(amin) - sm2, sm - sm2*drange/2)
amax = sm + sm2*drange
!print*, sm, sm2, amin, amax

map = max(amin, min(amax,log(map)))

map = (map-amin)/(amax-amin)
map = 3*map**2 - 2*map**3

pix = int(map*255)


!!!!! draw ahf halos !!!!!
pscale = isize/regionsize
do i=1, nh 
  mi0=int(normpos(ia,halo_pos(:,i)-center)*isize+1)
  mj0=int(normpos(ja,halo_pos(:,i)-center)*jsize+1)
  zfit=(normpos(ka,halo_pos(:,i)-center).ge.0).and.(normpos(ka,halo_pos(:,i)-center).le.1)
  mj = mj0
  do mi = mi0-int(halo_r(i)*pscale), mi0+int(halo_r(i)*pscale)
	if(mi.ge.1.and.mj.ge.1.and.mi.le.isize.and.mj.le.jsize.and.zfit) pix(mi,mj)=255
  enddo
  mi = mi0
  do mj = mj0-int(halo_r(i)*pscale), mj0+int(halo_r(i)*pscale)
	if(mi.ge.1.and.mj.ge.1.and.mi.le.isize.and.mj.le.jsize.and.zfit) pix(mi,mj)=255
  enddo
enddo

!!!!! draw grid !!!!!
if(do_grid) then
  ngrid = Box/grid_step
  do igrid = 1, ngrid
	mi = int(((igrid*grid_step-center(ia))/regionsize+0.5)*isize+1)
	if(mi.ge.1.and.mi.le.isize) pix(mi,:) = 128	
  enddo
  do jgrid = 1, ngrid
	mj = int(((jgrid*grid_step-center(ja))/regionsize+0.5)*jsize+1)
	if(mj.ge.1.and.mj.le.jsize) pix(:,mj) = 128	
  enddo
endif


do mi=1,isize
  do mj=1,jsize
    if(pix(mi,mj).eq.0) pix(mi,mj)=pix(mi,mj)+1
  enddo
enddo

if(palette.eq.1) then
	colmap=0
	do i=2,128
	  colmap(1,i)=i*((255-mincol)/128.0)+mincol
	enddo
	do i=129,255
	  colmap(2,i)=(i-129)*2
	enddo
	 colmap(1,129:255)=255
else
	colmap = 255
	do i=2,255
		colmap(:,i)=255-(i*((255-mincol)/256.0)+mincol)
	enddo
endif

 call writegif(fout(1:ilast)//'.gif', pix, colmap)

print*,'Written file ', fout(1:ilast)//'.gif'

map = 0
numout = numout+1
end subroutine plot

subroutine select_sphere()
  read(infile,*) radius
  if(.not.allocated(pos)) then
    print*,'Read something before selecting.'
    return
  endif
  Nsph = Npos
  do i=1,Npos
    r = sqrt(pos(1,i)*pos(1,i)+pos(2,i)*pos(2,i)+pos(3,i)*pos(3,i))
    if(r.gt.radius) then
      is_selected(i) = .false.
      Nsph = Nsph-1
    endif
  enddo
  vol = 4*3.141593/3*radius**3
  overdens = (Nsph/vol)/(Nall/Box**3)
  write(fstat,*) 'center=',center,' radius=',radius,' overdens=',overdens
end subroutine select_sphere

subroutine select_type()
  read(infile,*) itype
  if(.not.allocated(pos)) then
    print*,'Read something before selecting.'
    return
  endif
  is_selected = .false.
  do i=1,Npos
    is_selected(i) = itype.eq.ptype(i)
  enddo
!  !!! Warning: this does not work for multiple GADGET files
!  istart = 1+sum(Nall_arr(1:itype+1))-Nall_arr(itype+1)
!  iend = 1+sum(Nall_arr(1:itype+1))
!  is_selected = .false.
!  is_selected(istart:iend) = .true.
end subroutine select_type

subroutine unselect()
  is_selected = .true.
end subroutine unselect

subroutine print_mass()
  real*8:: fmass
  integer*8:: nsel
  fmass = 0d0
  do i = 1, Npos
      if(is_selected(i)) then
       fmass = fmass + massarr(ptype(i)+1)
      endif
  enddo
  print*, 'Total mass of selection: ',fmass
end subroutine print_mass

subroutine print_dist()
  rmin = 1e18
  do i=1,Npos
    if(is_selected(i)) then
      r = sqrt(pos(1,i)*pos(1,i)+pos(2,i)*pos(2,i)+pos(3,i)*pos(3,i))
      rmin = min(r,rmin)
    endif
  enddo
  print*,'Minimal distance to selection: ',rmin
end subroutine print_dist

subroutine matprod(mat)
  real:: mat(3,3), new(3,3)
  do i=1,3
    do j=1,3
      new(i,j) = sum(rot(i,1:3)*mat(1:3,j))
    enddo
  enddo
  rot(:,:) = new(:,:)
end subroutine matprod

subroutine rotate()
  real:: mat(3,3), norm
  read(infile,*) angle, x, y, z
  if(.not.allocated(pos)) then
    print*,'Read something before rotating.'
    return
  endif
  w = cos(angle/2/180*3.14159265)
  norm = sqrt((x*x+y*y+z*z)/(1-w*w))
  x=x/norm
  y=y/norm
  z=z/norm
  mat(1,1) = 1 - 2*y*y - 2*z*z
  mat(1,2) = 2*x*y - 2*z*w
  mat(1,3) = 2*x*z + 2*y*w
  mat(2,1) = 2*x*y + 2*z*w
  mat(2,2) = 1 - 2*x*x - 2*z*z
  mat(2,3) = 2*y*z - 2*x*w
  mat(3,1) = 2*x*z - 2*y*w
  mat(3,2) = 2*y*z + 2*x*w
  mat(3,3) = 1 - 2*x*x - 2*y*y
  call matprod(mat)
  was_rot = .true.
end subroutine rotate

subroutine ahf
	integer*8 :: hostid
	open(9,file=trim(ahf_file),status='old')
	read(9,*)
	nh=1
	do while(.true.)
	  read(9,*,end=90,err=90) a1,hostid,a3,a4,a5,(halo_pos(k,nh),k=1,3), a9,a10,a11, halo_r(nh)
	  nh=nh+1
	enddo
	90 close(2)
	nh=nh-1
	halo_pos = halo_pos/1e3
	halo_r = halo_r/1e3
	print*,'Halos read: ',nh
end subroutine ahf

! other subroutines

subroutine read_header(fname)
integer:: bytesleft((256-6*4-8*8-48)/4)
real*8:: time, redshift, boxsize
character(len=800):: fname

open(2,file=trim(fname),form='unformatted',status='old')

gformat=1
read(2,err=120) secname
if(secname(1:1).eq.'H')gformat=2
120 continue
if(gformat.eq.1) rewind(2)

read(2) npart, massarr, time, redshift, tmp, tmp, Nall_arr, tmp, numfiles, boxsize, bytesleft

Ntot=sum(npart)
Box=boxsize
tsnap = time
Nall = sum(Nall_arr)
end subroutine read_header

! help

subroutine help()
  print*
  print*,'     ---== SNAPTOOL HELP ==---'
  print*
  print*,'In-file contains a sequence of commands.'
  print*,'The line with the command name should not contain anything else:'
  print*,'the arguments, if required, are specified on subsequent line(s).'
  print*,'The in-file should finish with the end command.'
  print*
  print*,'List of commands:'
  print*
  print*,'# the rest of line is a comment'
  print*
  print*,'set map'
  print*,'isize jsize plane'
  print*,'  sets the map size in pixels and the projection plane,'
  print*,'  which is one of: xy, xz, yz.'
  print*,'  Default: 500 500 xy'
  print*
  print*,'read pos'
  print*,'path/to/snapshot_xxx.0'
  print*,'...'
  print*,'path/to/shapshot_xxx.last'
  print*,'  reads positions of particles from either a single-file or'
  print*,'  multiple-file GADGET snapshot. The number of files to read'
  print*,'  is determined automatically from the header of the first'
  print*,'  file in the list. Clears the map.'
  print*
  print*,'map'
  print*,'  fills the density map using particles currently read into memory.'
  print*,'  The previous contents are not cleared.'
  print*
  print*,'plot'
  print*,'  saves the map as an image, clears the map.'
  print*
  print*,'set read every'
  print*,'n'
  print*,'  read only every n-th particle from snapshots. Default: 1.'
  print*
  print*,'set center'
  print*,'X Y Z'
  print*,'  sets center of the region to map and plot. Default: Box/2 Box/2 Box/2'
  print*
  print*,'move center'
  print*,'X Y Z'
  print*,'  the same as "set center", but affects only the positions of particles'
  print*,'  that are already read, and does not affects next "read pos".'
  print*
  print*,'evolve center'
  print*,'X1 Y1 Z1 zstart'
  print*,'X2 Y2 Z2 zend'
  print*,'  defines the movement of the map center between two positions for each'
  print*,'  snapshot.'
  print*
  print*,'set size'
  print*,'S'
  print*,'  sets size of the working region in all three dimensions. Default: Box.'
  print*
  print*,'set cut'
  print*,'  or'
  print*,'unset cut'
  print*,'  chooses whether to cut the region of the size S and center X, Y, Z'
  print*,'  during reading the snapshot. This is useful to decrease memory usage.'
  print*
  print*,'set N pos max'
  print*,'N'
  print*,'  sets the maximum number of particles which can be allocated.'
  print*,'  Default: 2 000 000 000.'
  print*
  print*,'set outfile'
  print*,'filename'
  print*,'  sets the output filename base. Default: snapplot'
  print*
  print*,'add suffix z'
  print*,'  adds _zX.XXX to the output filenames.'
  print*,'add suffix axes'
  print*,'  adds _xy, _xz or _yz to the output filenames.'
  print*,'add suffix num'
  print*,'  adds _NNNNN to the output filenames. The number is increased'
  print*,'  each time writing occures.'
  print*
  print*,'include'
  print*,'filename'
  print*,'  run the sequence of commands from the given file,'
  print*,'  than return to the current file.'
  print*
  print*,'select sphere'
  print*,'R'
  print*,'  marks all particles with distance >R from center'
  print*,'  as non-plottable.'
  print*
  print*,'select type'
  print*,'i'
  print*,'  selects only particles of type i=0..5'
  print*
  print*,'print dist'
  print*,'  prints distance from previously set center to the'
  print*,'  closest particle from current selection.'
  print*
  print*,'rotate'
  print*,'angle x y z'
  print*,'  rotate everything on the given angle around the axis given'
  print*,'  by the three coordinates.'
  print*
  print*,'AHF'
  print*,'file'
  print*,'  display halos from AHF _halos catalog'
  print*
  print*,'set palette'
  print*,'N'
  print*,'  N=1 (default) - black-red-yellow. N>1 - white-black.'
  print*
  print*,'set grid'
  print*,'step'
  print*,'  plots a grid with given step size and origin at box origin.'
  print*
  print*,'load ramses map'
  print*,'filename'
  print*,'  load ascii file created by RAMSES part2map utility.'
end subroutine help

end program snaptool
