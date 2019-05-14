
!-------------------------------------------------------------------!
!                                                                   !
!                    Hybrid Gas-Kinetic Scheme                      !
!                                                                   !
!     CFD Software Package in Scientific and Engineering Research   !
!                                                                   !
! 					Developed by Tianbai Xiao (2019)                !  
!                                                                   !
! 			main.f90 - Main program for CFD simulation frame        !  
!                                                                   !
!-------------------------------------------------------------------!

program main

use premod

use initmod

use loopmod

use postmod

implicit none

!-------------------------------------------------------------------

!- set header information, could add software usage protection

call prejob

!- set initial status

call initjob

!- start main loop

call loopjob

!- postprocessing

call postjob

!-------------------------------------------------------------------

end program main