
! data.f90
!-------------------------------------------------------------------

module datamod

implicit none

! Data kind
integer, parameter :: srk = selected_real_kind (8)

! Math
integer, parameter :: dim = 2
integer, parameter :: q1dirc = 1 
integer, parameter :: q2dirc = 2 
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

integer,parameter :: recgeometry =1
integer,parameter :: polargeometry = 2
integer,parameter :: unigrid =1
integer,parameter :: refinegrid = 2

! Init
integer :: initype
logical :: isNewrun

! Computation
integer :: interp_method
integer :: output_method
integer :: geometry_method
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
real(srk) :: ylength

real(srk) :: rlength
real(srk) :: thetalength
real(srk) :: rin

integer :: npe
integer :: q1imin,q1imax
integer :: q2imin,q2imax

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
real(srk) :: rhoref,velref(dim),Tref,lambdaref
real(srk) :: momref(dim),Eref
real(srk) :: sosref

real(srk) :: mfp
real(srk) :: omega !HS model for tau
real(srk) :: alpha_ref !coefficient in VHS model
real(srk) :: omega_ref !coefficient in VHS model
real(srk) :: muref

real(srk) :: tau_ref

! Force
real(srk) :: phi

! Boundary
real(srk) :: rhowall,velwall,Twall,lambdawall
real(srk) :: bc_move(dim+2),bc_static(dim+2)
real(srk) :: bc_hot(dim+2),bc_cold(dim+2)
real(srk) :: bc_uniform(dim+2)
real(srk) :: bc_inverse(dim+2)
real(srk), allocatable, dimension(:,:) :: bc_flex

! Velocity space
integer :: unum,vnum
real(srk) :: umax,vmax
real(srk) :: umin,vmin
real(srk) :: lengthU,lengthV
real(srk) :: centerU,centerV
real(srk) :: du,dv
real(srk), allocatable, dimension(:,:) :: uspace,vspace
real(srk), allocatable, dimension(:,:) :: weight

! Class
type :: cell_center
	real(srk) :: x,y
	real(srk) :: r,theta
	real(srk) :: area
	real(srk) :: length(dim)
	real(srk) :: w(dim+2)
	real(srk) :: sw(dim+2,2)
	real(srk) :: prim(dim+2)
	real(srk),allocatable,dimension(:,:) :: h,b
	real(srk),allocatable,dimension(:,:,:) :: sh,sb
end type cell_center

type :: cell_interface
	integer :: bctype 
	real(srk) :: length
	real(srk) :: cosa,sina
	real(srk) :: flux(dim+2)
	real(srk),allocatable,dimension(:,:) :: flux_h,flux_b
	real(srk) :: qf(dim)
end type cell_interface

type(cell_center),allocatable,dimension(:,:) :: ctr
type(cell_interface),allocatable,dimension(:,:) :: a1face,a2face

end module datamod