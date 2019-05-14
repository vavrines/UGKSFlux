
! postprocess.f90
!-------------------------------------------------------------------

module postmod

use datamod

use kineticmod

implicit none

contains

!-------------------------------------------------------------------

subroutine postjob

write(*,*) 'simulation time: ', simtime
write(*,*) 'iteration: ', iter

call writeMacro(iter)

end subroutine postjob

!-------------------------------------------------------------------

subroutine writeMacro(loop)

integer,intent(in):: loop
character(len=8) :: ctemp
real(srk),dimension(:,:),allocatable :: csolution
real(srk) :: prim(DIM+2)
integer :: i
real(srk) :: rmid,xmid

allocate(csolution(10,q1imax))

do i=q1imin,q1imax
	prim = get_primitive(ctr(i)%w)

	csolution(1:2,i) = prim(1:2)
	csolution(3,i) = 1.d0/prim(3)
	csolution(4,i) = 0.5d0*csolution(3,i)*csolution(1,i)
	csolution(5,i) = get_heat_flux(ctr(i)%h,ctr(i)%b,prim)

	csolution(6,i) = ctr(i)%rgtype
end do

!find middle location - the location of average density
rmid = 0.5*(ctr(q1imin-1)%prim(1)+ctr(q1imax+1)%prim(1))

do i=q1imin,q1imax
	if ((ctr(i)%prim(1)-rmid)*(ctr(i+1)%prim(1)-rmid)<=0) then
	    xmid = ctr(i)%x+(ctr(i+1)%x-ctr(i)%x)/(ctr(i+1)%prim(1)-ctr(i)%prim(1))*(rmid-ctr(i)%prim(1))
	end if
end do

!normalize
!csolution(1,:) = (csolution(1,:)-ctr(q1imin-1)%prim(1))/(ctr(q1imax+1)%prim(1)-ctr(q1imin-1)%prim(1))
!csolution(3,:) = (csolution(3,:)-1.d0/ctr(q1imin-1)%prim(3))/(1.d0/ctr(q1imax+1)%prim(3)-1.d0/ctr(q1imin-1)%prim(3))

write(ctemp,'(i8)') loop
open(unit=fileOutId,file=RSTFILE//trim(adjustl(cTemp))//'.dat',status="replace",action="write")
write(fileOutId,*) "VARIABLES = X, Density, Velocity, Temperature, Pressure, Qx, Scheme"

write(fileOutId,'(A,I5,A,I5,A)') 'ZONE I =',q1imax-q1imin+1,', DATAPACKING=BLOCK'
write(fileOutId,"(6(ES23.16,2X))") ctr(q1imin:q1imax)%x!-xmid
do i=1,6
	write(fileOutId,"(6(ES23.16,2X))") csolution(i,:)
end do

close(fileOutId)

end subroutine writeMacro

!-------------------------------------------------------------------

end module postmod