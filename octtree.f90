module octtree
	
	type node
		type(node), pointer :: next(:,:,:), parent
		integer:: start=-1
		integer:: np=0
	end type node
	
	type tree
		type(node), pointer :: top
		real :: box(3), center(3)
		logical :: periodic
		integer :: nmax, np=0, nlevels
		integer, pointer :: list(:)
	end type tree
contains
type(tree) function tree_new(box,center,np,nm)
	real:: box(3), center(3)
	logical :: periodic
	tree_new%box = box
	tree_new%center = center
	tree_new%periodic = .true.
	tree_new%nmax = nm
	tree_new%np = np
	allocate(tree_new%top)
end function tree_new

subroutine tree_delete(tr)
	type(tree):: tr
	if(associated(tr%list)) nullify(tr%list)
	if(associated(tr%top)) call tree_rmnodes(tr%top,tr%nmax)
end subroutine tree_delete

recursive subroutine tree_rmnodes(nd,nmax)
	type(node):: nd
	if(nd%np.gt.nmax) then
		do i=1,2
		do j=1,2
		do k=1,2
			call tree_rmnodes(nd%next(i,j,k),nmax)
		enddo
		enddo
		enddo
		nullify(nd%next)
	endif
end subroutine tree_rmnodes

subroutine tree_construct(tr,xyz)
	type(tree):: tr
	real:: xyz(1:3,tr%np)
	logical:: donext
	tr%top%np=0
	tr%top%start=-1
	level = 1
	donext = .true.
	do while (donext)
		donext=.false.
		do i=1,tr%np
			call tree_mark(xyz(1,i),tr,level,donext)
		enddo
		level=level+1
	enddo
	tr%nlevels=level
end subroutine tree_construct

subroutine tree_makelist(tr,xyz)
	type(tree):: tr
	real:: xyz(1:3,tr%np)
	!if(.not.associated(tr%list)) 
	!if(associated(tr%list)) nullify(tr%list)
	allocate(tr%list(tr%np))
	do i=1,tr%np
		call tree_addtolist(xyz(1,i),tr,i)
	enddo
end subroutine tree_makelist

function tree_find(xyz,tr,level) result(p)
	type(tree):: tr
	real:: xyz(3), center(3), dc(3)
	integer:: addr(3)
	integer, intent(out)::level
	type(node), pointer:: p
	center = tr%center
	dc = tr%box/4
	p => tr%top
	level = 1
	do while(p%np.ge.tr%nmax)
		do k=1,3
			if(xyz(k).gt.center(k)) then
				addr(k)=2
				center(k)=center(k)+dc(k)
			else
				addr(k)=1
				center(k)=center(k)-dc(k)
			endif
		enddo
		dc=dc/2
		p=>p%next(addr(1),addr(2),addr(3))
		level = level+1
	enddo
end function tree_find

subroutine tree_addtolist(xyz,tr,ind)
	type(tree):: tr
	real:: xyz(3)
	type(node), pointer:: p
	p => tree_find(xyz,tr,lev)
	tr%list(ind) = p%start
	p%start = ind
end subroutine tree_addtolist

function tree_dens(xyz,tr)
	type(tree):: tr
	real:: xyz(3)
	type(node), pointer:: p
	p => tree_find(xyz,tr,lev)
	if(p%np.lt.1) then
		p => p%parent
		lev = lev-1
	endif
	tree_dens = p%np/(tr%box(1)*tr%box(2)*tr%box(3)/2**(3*(lev-1)))
end function tree_dens

	
subroutine tree_mark(xyz,tr,level,donext)
	type(tree):: tr
	real:: xyz(3), center(3), dc(3)
	integer:: addr(3)
	type(node), pointer:: p
	logical:: donext
	center = tr%center
	do k=1,3
		if(xyz(k).lt.center(k)-tr%box(k)/2 .or. xyz(k).gt.center(k)+tr%box(k)/2) return
	enddo
	dc = tr%box/4
	p => tr%top
	do l = 1,level-1
		if(p%np.lt.tr%nmax) return
		do k=1,3
			if(xyz(k).gt.center(k)) then
				addr(k)=2
				center(k)=center(k)+dc(k)
			else
				addr(k)=1
				center(k)=center(k)-dc(k)
			endif
		enddo
		dc=dc/2
		p=>p%next(addr(1),addr(2),addr(3))
	enddo
	p%np = p%np+1
	if(p%np.eq.tr%nmax) then
		allocate(p%next(2,2,2))
		!allocate(p%next(1:2,1:2,1:2)%parent)
		p%next(1:2,1:2,1:2)%np = 0
		p%next(1:2,1:2,1:2)%start = -1
		do i=1,2
		do j=1,2
		do k=1,2
			p%next(i,j,k)%parent => p
		enddo
		enddo
		enddo
		donext=.true.
	endif
end subroutine tree_mark
end module octtree
