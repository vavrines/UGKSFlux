
! postprocess.f90
!-------------------------------------------------------------------

module postmod

use datamod

use gridmod

use kineticmod

implicit none

contains

!-------------------------------------------------------------------

subroutine postjob

write(*,*) 'simulation time: ', simtime
write(*,*) 'iteration: ', iter

call writeMacro(iter)
call writeMicro(iter,q1imax/2,q2imax/2)
call writeNeq(iter)

end subroutine postjob

!-------------------------------------------------------------------

subroutine writeMacro(loop)

integer,intent(in) :: loop
character(len=7) :: ctemp
real(srk),dimension(:,:,:),allocatable :: csolution
real(srk) :: prim(dim+2)
integer :: i,j

allocate(csolution(7,q1imin:q1imax,q2imin:q2imax))

do i=1,q1imax
do j=1,q2imax
	prim = get_primitive(ctr(i,j)%w)

	csolution(1:3,i,j) = prim(1:3)
	csolution(4,i,j) = 1.d0/prim(4)
	csolution(5,i,j) = 0.5*prim(1)/prim(4)
	csolution(6:7,i,j) = get_heat_flux(ctr(i,j)%h,ctr(i,j)%b,uspace,vspace,prim)
end do
end do

write(ctemp,'(i7)') loop
open(unit=fileoutid,file=rstfile//trim(adjustl(ctemp))//'.dat',status="replace",action="write")
write(fileoutid,*) "VARIABLES = X, Y, Density, U, V, Temperature, Pressure, Qx, Qy"

select case(output_method)
	case(pointsvalue)
	    write(fileoutid,*) "ZONE  I = ",grid%q1index,", J = ",grid%q2index,"DATAPACKING=BLOCK, VARLOCATION=([3-9]=CELLCENTERED)"

	    !write geometry (node value)
	    write(fileoutid,"(6(ES23.16,2X))") coordS(1,:,:)
		write(fileoutid,"(6(ES23.16,2X))") coordS(2,:,:)
	case(centervalue)
	    write(fileoutid,*) "ZONE  I = ",grid%q1index-1,", J = ",grid%q2index-1,"DATAPACKING=BLOCK"

	    !write geometry (cell centered value)
	    write(fileoutid,"(6(ES23.16,2X))") ctr(q1imin:q1imax,q2imin:q2imax)%x
	    write(fileoutid,"(6(ES23.16,2X))") ctr(q1imin:q1imax,q2imin:q2imax)%y
end select

do i=1,7
	write(fileoutid,"(6(ES23.16,2X))") csolution(i,:,:)
end do

close(fileoutid)

end subroutine writeMacro

!-------------------------------------------------------------------

subroutine writeMicro(loop,q1idx,q2idx)

integer,intent(in) :: loop
integer,intent(in) :: q1idx,q2idx
character(len=7) :: ctemp
integer :: i,j

write(ctemp,'(i7)') loop
open(unit=fileoutid,file=MSTFILE//trim(adjustl(cTemp))//'.dat',status="replace",action="write")
write(fileoutid,*) "VARIABLES = U, V, H, B"
write(fileoutid,'(A,I5,A,I5,A)') 'ZONE I =',unum, ', J =',vnum, ', DATAPACKING=BLOCK'

write(fileoutid,"(6(ES23.16,2X))") uspace(:,:)
write(fileoutid,"(6(ES23.16,2X))") vspace(:,:)
write(fileoutid,"(6(ES23.16,2X))") ctr(q1idx,q2idx)%h
write(fileoutid,"(6(ES23.16,2X))") ctr(q1idx,q2idx)%b

close(fileoutid)

end subroutine writeMicro

!-------------------------------------------------------------------

subroutine writeNeq(loop)

integer,intent(in) :: loop
character(len=7) :: ctemp

real(srk),dimension(:,:),allocatable :: h,b
real(srk),dimension(:,:,:),allocatable :: neq
real(srk) :: prim(dim+2)
integer :: i,j
integer :: k,l

allocate(neq(2,q1imin:q1imax,q2imin:q2imax))
allocate(h(unum,vnum))
allocate(b(unum,vnum))

neq = 0.d0

do i=1,q1imax
do j=1,q2imax
	prim = get_primitive(ctr(i,j)%w)
    call reduced_maxwell(h,b,uspace,vspace,prim)
	
    do k=1,unum
        do l=1,vnum
            neq(1,i,j) = neq(1,i,j)+abs(ctr(i,j)%h(k,l)-h(k,l))
            neq(2,i,j) = neq(2,i,j)+abs(ctr(i,j)%b(k,l)-b(k,l))
        end do
    end do
end do
end do

write(ctemp,'(i7)') loop
open(unit=fileoutid,file=neqfile//trim(adjustl(ctemp))//'.dat',status="replace",action="write")
write(fileoutid,*) "VARIABLES = X, Y, Dh, Db"

select case(output_method)
	case(pointsvalue)
	    write(fileoutid,*) "ZONE  I = ",grid%q1index,", J = ",grid%q2index,"DATAPACKING=BLOCK, VARLOCATION=([3-9]=CELLCENTERED)"

	    !write geometry (node value)
	    write(fileoutid,"(6(ES23.16,2X))") coordS(1,:,:)
		write(fileoutid,"(6(ES23.16,2X))") coordS(2,:,:)
	case(centervalue)
	    write(fileoutid,*) "ZONE  I = ",grid%q1index-1,", J = ",grid%q2index-1,"DATAPACKING=BLOCK"

	    !write geometry (cell centered value)
	    write(fileoutid,"(6(ES23.16,2X))") ctr(q1imin:q1imax,q2imin:q2imax)%x
	    write(fileoutid,"(6(ES23.16,2X))") ctr(q1imin:q1imax,q2imin:q2imax)%y
end select

do i=1,2
	write(fileoutid,"(6(ES23.16,2X))") neq(i,:,:)
end do

close(fileoutid)

end subroutine writeNeq

!-------------------------------------------------------------------

end module postmod