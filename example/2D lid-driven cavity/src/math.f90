
! math.f90
!-------------------------------------------------------------------

module mathmod

use datamod

implicit none

contains

!-------------------------------------------------------------------

function global_frame(w,cosa,sina)

real(srk),intent(in) :: w(dim+2)
real(srk),intent(in) :: cosa,sina
real(srk) :: global_frame(dim+2)

global_frame(1) = w(1)
global_frame(2) = w(2)*cosa-w(3)*sina
global_frame(3) = w(2)*sina+w(3)*cosa
global_frame(4) = w(4)

end function global_frame


function local_frame(w,cosa,sina)
real(srk),intent(in) :: w(dim+2)
real(srk),intent(in) :: cosa,sina
real(srk) :: local_frame(dim+2)

local_frame(1) = w(1)
local_frame(2) = w(2)*cosa+w(3)*sina
local_frame(3) = w(3)*cosa-w(2)*sina
local_frame(4) = w(4)

end function local_frame

!-------------------------------------------------------------------

subroutine interp_boundary(cell_L,cell_N,cell_R,dirc)

type(cell_center),intent(inout) :: cell_N
type(cell_center),intent(inout) :: cell_L,cell_R
integer,intent(in) :: dirc

cell_N%sh(:,:,dirc) = (cell_R%h-cell_L%h)/(0.5*cell_R%length(dirc)+0.5*cell_L%length(dirc))
cell_N%sb(:,:,dirc) = (cell_R%b-cell_L%b)/(0.5*cell_R%length(dirc)+0.5*cell_L%length(dirc))

cell_N%sw(:,dirc) = (cell_R%w-cell_L%w)/(0.5*cell_R%length(dirc)+0.5*cell_L%length(dirc))

end subroutine interp_boundary


subroutine interp_inner(cell_L,cell_N,cell_R,dirc)

type(cell_center),intent(in) :: cell_L,cell_R
type(cell_center),intent(inout) :: cell_N
integer,intent(in) :: dirc
real(srk),allocatable,dimension(:,:) :: sL,sR
real(srk) :: swL(dim+2),swR(dim+2)

allocate(sL(unum,vnum))
allocate(sR(unum,vnum))

sL = (cell_N%h-cell_L%h)/(0.5*cell_N%length(dirc)+0.5*cell_L%length(dirc))
sR = (cell_R%h-cell_N%h)/(0.5*cell_R%length(dirc)+0.5*cell_N%length(dirc))
cell_N%sh(:,:,dirc) = (sign(up,sR)+sign(up,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)

sL = (cell_N%b-cell_L%b)/(0.5*cell_N%length(dirc)+0.5*cell_L%length(dirc))
sR = (cell_R%b-cell_N%b)/(0.5*cell_R%length(dirc)+0.5*cell_N%length(dirc))
cell_N%sb(:,:,dirc) = (sign(up,sR)+sign(up,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)

swL = (cell_N%w-cell_L%w)/(0.5*cell_N%length(dirc)+0.5*cell_L%length(dirc))
swR = (cell_R%w-cell_N%w)/(0.5*cell_R%length(dirc)+0.5*cell_N%length(dirc))
cell_N%sw(:,dirc) = (sign(up,swR)+sign(up,swL))*abs(swR)*abs(swL)/(abs(swR)+abs(swL)+smv)

end subroutine interp_inner

!-------------------------------------------------------------------

end module mathmod