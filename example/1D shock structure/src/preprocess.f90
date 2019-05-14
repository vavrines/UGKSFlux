
! preprocess.f90
!-------------------------------------------------------------------

module premod

use datamod

implicit none

contains

!-------------------------------------------------------------------

subroutine prejob

call preHead

call preSet

end subroutine prejob

!-------------------------------------------------------------------

subroutine preHead

character(len=20) :: separator20

separator20 = '--------------------'

write(*, '(4a20)') separator20, separator20, separator20, separator20

write(*,*) 'Welcome to Cloud CFD 1D (version 0.1)'
write(*,*) 'A Software Package in Scientific and Engineering Simulation'

write(*, '(4a20)') separator20, separator20, separator20, separator20

end subroutine preHead

!-------------------------------------------------------------------

subroutine preSet

isNewrun = .true.

gridfile = 'CloudCFD.neu'
initype = 1 ! 1:generate 2:read

rstfile = "macro"
mstfile = "micro"
neqfile = "noneq"

hstfile = "CloudCFD.hst"

end subroutine preSet

!-------------------------------------------------------------------

end module premod