
! math.f90
!-------------------------------------------------------------------

module mathmod

use datamod

implicit none

contains

!-------------------------------------------------------------------

subroutine interp_boundary(cell_L,cell_N,cell_R)

type(cell_center),intent(inout):: cell_N
type(cell_center),intent(inout):: cell_L,cell_R

cell_N%sh = (cell_R%h-cell_L%h)/(0.5*cell_R%length+0.5*cell_L%length)
cell_N%sb = (cell_R%b-cell_L%b)/(0.5*cell_R%length+0.5*cell_L%length)
cell_N%sw = (cell_R%w-cell_L%w)/(0.5*cell_R%length+0.5*cell_L%length)

end subroutine interp_boundary

!-------------------------------------------------------------------

subroutine interp_micro(cell_L,cell_N,cell_R)

type(cell_center),intent(in) :: cell_L,cell_R
type(cell_center),intent(inout) :: cell_N
real(srk),allocatable,dimension(:) :: sL,sR

allocate(sL(unum))
allocate(sR(unum))

sL = (cell_N%h-cell_L%h)/(0.5*cell_N%length+0.5*cell_L%length)
sR = (cell_R%h-cell_N%h)/(0.5*cell_R%length+0.5*cell_N%length)
cell_N%sh = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)

sL = (cell_N%b-cell_L%b)/(0.5*cell_N%length+0.5*cell_L%length)
sR = (cell_R%b-cell_N%b)/(0.5*cell_R%length+0.5*cell_N%length)
cell_N%sb = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)

end subroutine interp_micro

!-------------------------------------------------------------------

subroutine interp_macro(cell_L,cell_N,cell_R)

type(cell_center),intent(in) :: cell_L,cell_R
type(cell_center),intent(inout) :: cell_N
real(srk) :: sL(dim+2),sR(dim+2)

sL = (cell_N%w-cell_L%w)/(0.5*cell_N%length+0.5*cell_L%length)
sR = (cell_R%w-cell_N%w)/(0.5*cell_R%length+0.5*cell_N%length)
cell_N%sw = (sign(UP,sR)+sign(UP,sL))*abs(sR)*abs(sL)/(abs(sR)+abs(sL)+SMV)

end subroutine interp_macro

!-------------------------------------------------------------------

end module mathmod