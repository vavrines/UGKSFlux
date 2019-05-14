
! grid.f90
!-------------------------------------------------------------------

module gridmod

use datamod

implicit none

real(srk),pointer :: coordUS(:,:)  !unstructured grid
real(srk),pointer :: coordS(:,:,:) !structured grid

type gridtype
	integer :: numNp,numElem,numEdge,nodePerElem

	! mapping between Elem and Np
	integer,pointer :: mapEN(:,:) !每个单元对应的节点
	integer,pointer :: edgeMapNp(:,:) !每条边对应两个节点
	integer,pointer :: edgeNeighborElem(:,:) !每条边对应的单元
	integer,pointer :: elemNeighbor(:,:) !单元相邻的单元
	integer,pointer :: elemMapEdge(:,:) !!每个单元对应的边

	! index of structured node points/elements
	integer :: q1index, q2index

	! map from 3D list to 1D list
	integer,pointer :: mapStoU(:,:) !结构网格点对应非结构网格点编号
end type gridtype

type(gridtype) :: grid	

integer, parameter :: selfgnt = 1
integer, parameter :: loadmsh = 2

contains

!-------------------------------------------------------------------

subroutine initGrid

logical :: found

select case(initype)

	case(selfgnt)

		write(*,*) "Generating grid"

		call gntCartesian
		
	case(loadmsh)
	
		write(*,*) "Loading grid"

		inquire(file=gridfile,exist=found)

		if (.not. found) then
			write(*,*) 'Grid file "',trim(gridfile),'" does not exist'
			stop
		else
			write(*,*) 'Grid file is ',trim(gridfile)
		end if

end select

end subroutine initGrid

!-------------------------------------------------------------------

subroutine gntCartesian

real(srk) :: dq1,dq2
real(srk) :: s1,s2
integer :: i,j

! grid type
grid%q1index = 46
grid%q2index = 46
grid%numNp = grid%q1index*grid%q2index
grid%numElem = (grid%q1index-1)*(grid%q2index-1)
grid%nodePerElem = 4

! element
q1imin = 1
q1imax = grid%q1index-1
q2imin = 1
q2imax = grid%q2index-1

write(*,*) "x cell number: ",q1imax-q1imin+1
write(*,*) "y cell number: ",q2imax-q2imin+1

dq1 = xlength/(q1imax-q1imin+1)
dq2 = ylength/(q2imax-q2imin+1)

if(associated(coordS)) then
	deallocate(coordS)
endif
allocate(coordS(2,grid%q1index,grid%q2index))

allocate(ctr(q1imin-2:q1imax+2,q2imin-2:q2imax+2)) 
allocate(a1face(grid%q1index,grid%q2index-1))
allocate(a2face(grid%q1index-1,grid%q2index))

do i=1,grid%q1index
do j=1,grid%q2index
	coordS(1,i,j) = (i-1)*dq1
	coordS(2,i,j) = (j-1)*dq2
end do
end do

do i=q1imin-2,q1imax+2
do j=q2imin-2,q2imax+2
	ctr(i,j)%x = (i-0.5d0)*dq1
	ctr(i,j)%y = (j-0.5d0)*dq2
	ctr(i,j)%area = dq1*dq2
	ctr(i,j)%length(1) = dq1
	ctr(i,j)%length(2) = dq2
end do
end do

do i=1,grid%q1index
do j=1,grid%q2index-1
	a1face(i,j)%length = dq2
	a1face(i,j)%cosa = 1.d0
	a1face(i,j)%sina = 0.d0
end do
end do

do i=1,grid%q1index-1
do j=1,grid%q2index
	a2face(i,j)%length = dq1
	a2face(i,j)%cosa = 0.d0
	a2face(i,j)%sina = 1.d0
end do
end do

end subroutine gntCartesian

!-------------------------------------------------------------------

end module gridmod