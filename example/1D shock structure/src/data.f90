
! data.f90
!-------------------------------------------------------------------

module datamod

implicit none

! Data kind
integer, parameter :: srk = selected_real_kind (8)

! Math
integer, parameter :: dim = 1
integer, parameter :: inside = 1 
integer, parameter :: outside = 2 

real(srk), parameter :: pi = 3.1415926535d0
real(srk), parameter :: smv = 1.d-8
real(srk), parameter :: eps = 1.d-7
real(srk), parameter :: up = 1.d0
integer, parameter :: rn = 1 !no frame rotation
integer, parameter :: ry = -1 !with frame rotation

integer, parameter :: firstorder = 1
integer, parameter :: secondorder = 2 
integer, parameter :: thirdorder = 3 
integer, parameter :: fourthorder = 4 
integer, parameter :: centervalue = 1
integer, parameter :: pointsvalue = 2
integer, parameter :: continuum = 1
integer, parameter :: transition = 2
integer, parameter :: rarefied = 3
integer,parameter :: unigrid =1
integer,parameter :: refinegrid = 2

! Init
integer :: initype
logical :: isNewrun

! Computation
integer :: interp_method
integer :: grid_method

integer :: iter
real(srk) :: cfl
real(srk) :: simtime,maxtime
real(srk) :: dt
real(srk) :: res(dim+2)

! File
integer, parameter :: fileinid = 17
integer, parameter :: fileoutid = 18
integer, parameter :: fileinitid = 19
integer, parameter :: filehstid = 20
character(len=15) :: inifile
character(len=15) :: gridfile
character(len=15) :: hstfile, rstfile, mstfile
character(len=15) :: neqfile

! Geometry
real(srk) :: xlength
integer :: q1imin,q1imax

! Physics
real(srk), parameter :: R0 = 8.314d0
real(srk), parameter :: bz = 1.3806488d-23

! Gas property
real(srk) :: inK
real(srk) :: namass !mole mass
real(srk) :: molediameter !molecule diameter
real(srk) :: molemass !molecule mass
real(srk) :: gasR
real(srk) :: gamma

! Non-dimensional parameters
real(srk) :: prandtl
real(srk) :: knudsen
real(srk) :: mach

! Reference state
real(srk) :: rhoref,velref,Tref,lambdaref
real(srk) :: momref,Eref
real(srk) :: sosref

real(srk) :: mfp
real(srk) :: omega !HS model for tau
real(srk) :: alpha_ref !coefficient in VHS model
real(srk) :: omega_ref !coefficient in VHS model
real(srk) :: muref

! Force
real(srk) :: phi
real(srk),allocatable,dimension(:) :: pot

! Boundary
real(srk) :: bcL(dim+2),bcR(dim+2)

! Shock tube
real(srk) :: rhoLeft,rhoRight
real(srk) :: TLeft,TRight
real(srk) :: lambdaLeft,lambdaRight
real(srk) :: ELeft,ERight

! Velocity space
integer :: unum
real(srk) :: umax
real(srk) :: umin
real(srk) :: lengthU
real(srk) :: centerU
real(srk) :: du
real(srk), allocatable, dimension(:) :: uspace
real(srk), allocatable, dimension(:) :: weight

! Class
type :: cell_center
	real(srk) :: x
	real(srk) :: length
	real(srk) :: w(dim+2)
	real(srk) :: sw(dim+2)
	real(srk) :: prim(dim+2)
	real(srk),allocatable,dimension(:) :: h,b
	real(srk),allocatable,dimension(:) :: sh,sb
	!flow regime
	integer :: rgtype
	!weno reconstruction
	real(srk),allocatable,dimension(:) :: hl,bl
	real(srk),allocatable,dimension(:) :: hr,br
end type cell_center

type :: cell_interface
	integer :: bctype 
	real(srk) :: flux(dim+2)
	real(srk),allocatable,dimension(:) :: flux_h,flux_b
	real(srk) :: qf
end type cell_interface

type(cell_center),allocatable,dimension(:) :: ctr
type(cell_interface),allocatable,dimension(:) :: face

end module datamod