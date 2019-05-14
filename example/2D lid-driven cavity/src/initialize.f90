
! initialize.f90
!-------------------------------------------------------------------

module initmod

use datamod

use gridmod

use kineticmod

use postmod

implicit none

contains

!-------------------------------------------------------------------

subroutine initjob

if(isNewrun) then
	call initNew
else
	!call initRestart
end if

end subroutine initjob

!-------------------------------------------------------------------

subroutine initNew

call initPara

call initGeometry

call initGrid

call initGas

call initVelspace(velref,lambdaref,8.d0)

call initBoundary

call initExforce

call initFlowfield

call initIo

end subroutine initNew

!-------------------------------------------------------------------

subroutine initPara

interp_method = secondorder
output_method = pointsvalue

cfl = 0.8d0
simtime = 0.d0
maxtime = 1.0d-5
iter = 1

write(*,*) "CFL number: ",cfl

end subroutine initPara

!-------------------------------------------------------------------

subroutine initGeometry

! Cartesian
xlength = 1.d0
ylength = 1.d0

write(*,*) "xlength: ",xlength
write(*,*) "ylength: ",ylength

end subroutine initGeometry

!-------------------------------------------------------------------

subroutine initGas

! Argon gas
prandtl = 2.d0/3.d0
knudsen = 0.075d0

namass = 39.95d0
molediameter = 4.17d-10
molemass = 6.63d-26
gasR = R0/namass*1.d3 

! Hard sphere model
inK = 1.d0
gamma = get_gamma(inK)

omega = 0.72d0
alpha_ref = 1.d0
omega_ref = 0.5d0
muref = get_mu(knudsen,alpha_ref,omega_ref)

write(*,*) "Prandtl number: ",prandtl
write(*,*) "Knudsen number: ",knudsen
write(*,*) "gamma: ",gamma
write(*,*) "dynamic viscosity: ",muref

! Reference state
rhoref = 1.d0
velref(1) = 0.d0
velref(2) = 0.d0
Tref = 1.d0
lambdaref = 1.d0/Tref
sosref = 1.d0

momref = rhoref*velref
Eref = 0.5d0*sum(momref*momref)/rhoref+&
       0.5d0*rhoref/lambdaref/(gamma-1.0)

end subroutine initGas

!-------------------------------------------------------------------

subroutine initVelspace(vel,lambda,alpha)

real(srk),intent(in) :: vel(dim),lambda,alpha
real(srk) :: referVel

lengthU = 2.d0*alpha*dsqrt(0.5d0*gamma/lambda)
centerU = vel(1)

lengthV = 2.d0*alpha*dsqrt(0.5d0*gamma/lambda)
centerV = vel(2)

umax = centerU+0.5d0*lengthU
umin = centerU-0.5d0*lengthU
vmax = centerV+0.5d0*lengthV
vmin = centerV-0.5d0*lengthV

!call initNewton
call initGauss(velref(1),velref(2))

end subroutine initVelspace

!-------------------------------------------------------------------
	
subroutine initNewton

integer :: i,j

unum = 41
vnum = 41

write(*,*) "Newton-Cotes velocity space"
write(*,*) "velocity points: ",unum,vnum
write(*,*) "U range: ",umin,umax
write(*,*) "V range: ",vmin,vmax

allocate(uspace(unum,vnum))
allocate(vspace(unum,vnum))
allocate(weight(unum,vnum))

du = (umax-umin)/(unum-1)
dv = (vmax-vmin)/(vnum-1)

forall(i=1:unum,j=1:vnum)
	uspace(i,j) = umin+(i-1)*du
	vspace(i,j) = vmin+(j-1)*dv

	weight(i,j) = (newton_coeff(i,unum)*du)*(newton_coeff(j,vnum)*dv)
end forall

contains

  pure function newton_coeff(idx,num)
      integer,intent(in) :: idx,num
      real(srk) :: newton_coeff

      if (idx==1 .or. idx==num) then
          newton_coeff = 14.0/45.0
      else if (mod(idx-5,4)==0) then
          newton_coeff = 28.0/45.0
      else if (mod(idx-3,4)==0) then
          newton_coeff = 24.0/45.0
      else
          newton_coeff = 64.0/45.0
      end if
  end function newton_coeff

end subroutine initNewton

!-------------------------------------------------------------------

subroutine initGauss(umid,vmid)

real(srk),intent(in) :: umid,vmid
real(srk) :: vcoords(28), weights(28)
integer :: i,j

write(*,*) "Gauss velocity space"

vcoords=[ -0.5392407922630E+01, -0.4628038787602E+01, -0.3997895360339E+01, -0.3438309154336E+01,&
          -0.2926155234545E+01, -0.2450765117455E+01, -0.2007226518418E+01, -0.1594180474269E+01,&
          -0.1213086106429E+01, -0.8681075880846E+00, -0.5662379126244E+00, -0.3172834649517E+00,&
          -0.1331473976273E+00, -0.2574593750171E-01, +0.2574593750171E-01, +0.1331473976273E+00,&
          +0.3172834649517E+00, +0.5662379126244E+00, +0.8681075880846E+00, +0.1213086106429E+01,&
          +0.1594180474269E+01, +0.2007226518418E+01, +0.2450765117455E+01, +0.2926155234545E+01,&
          +0.3438309154336E+01, +0.3997895360339E+01, +0.4628038787602E+01, +0.5392407922630E+01 ]

weights=[ +0.2070921821819E-12, +0.3391774320172E-09, +0.6744233894962E-07, +0.3916031412192E-05,&
          +0.9416408715712E-04, +0.1130613659204E-02, +0.7620883072174E-02, +0.3130804321888E-01,&
          +0.8355201801999E-01, +0.1528864568113E+00, +0.2012086859914E+00, +0.1976903952423E+00,&
          +0.1450007948865E+00, +0.6573088665062E-01, +0.6573088665062E-01, +0.1450007948865E+00,&
          +0.1976903952423E+00, +0.2012086859914E+00, +0.1528864568113E+00, +0.8355201801999E-01,&
          +0.3130804321888E-01, +0.7620883072174E-02, +0.1130613659204E-02, +0.9416408715712E-04,&
          +0.3916031412192E-05, +0.6744233894962E-07, +0.3391774320172E-09, +0.2070921821819E-12 ]

unum = 28
vnum = 28

write(*,*) "velocity points: ",unum,vnum

allocate(uspace(unum,vnum)) 
allocate(vspace(unum,vnum)) 
allocate(weight(unum,vnum))

forall(i=1:unum,j=1:vnum)
    uspace(i,j) = umid+vcoords(i)
    vspace(i,j) = vmid+vcoords(j)
    weight(i,j) = (weights(i)*exp(vcoords(i)**2))*(weights(j)*exp(vcoords(j)**2))
end forall

!find the maximum micro velocity for the determination of dt
umax = maxval(abs(uspace(:,1)))
vmax = maxval(abs(vspace(1,:)))
umin = -umax
vmin = -vmax
du = (umax-umin)/unum
dv = (vmax-vmin)/vnum

write(*,*) "U range: ",umin,umax
write(*,*) "V range: ",vmin,vmax

end subroutine initGauss

!-------------------------------------------------------------------

subroutine initBoundary

bc_uniform = [ rhoref,0.15d0,0.d0,lambdaref ]
bc_move = [ rhoref,0.15d0,0.d0,lambdaref ]
bc_static = [ rhoref,0.d0,0.d0,lambdaref ]

end subroutine initBoundary

!-------------------------------------------------------------------

subroutine initExforce

phi = 0.d0

write(*,*) "external acceleration: ",phi

end subroutine initExforce

!-------------------------------------------------------------------

subroutine initFlowfield

real(srk),allocatable,dimension(:,:) :: h0,b0
integer :: i,j

do i=q1imin-2,q1imax+2
do j=q2imin-2,q2imax+2
  ctr(i,j)%prim(1) = rhoref

  ctr(i,j)%prim(2) = 0.d0
  ctr(i,j)%prim(3) = 0.d0
  ctr(i,j)%prim(4) = lambdaref

  ctr(i,j)%w = get_conserved(ctr(i,j)%prim)
end do
end do

allocate(h0(unum,vnum))
allocate(b0(unum,vnum))

do i=q1imin-2,q1imax+2
do j=q2imin-2,q2imax+2
  allocate(ctr(i,j)%h(unum,vnum))
  allocate(ctr(i,j)%b(unum,vnum))
  allocate(ctr(i,j)%sh(unum,vnum,dim))
  allocate(ctr(i,j)%sb(unum,vnum,dim))
end do
end do

do i=q1imin,q1imax+1
do j=q2imin,q2imax
  allocate(a1face(i,j)%flux_h(unum,vnum))
  allocate(a1face(i,j)%flux_b(unum,vnum))
end do
end do

do i=q1imin,q1imax
do j=q2imin,q2imax+1
  allocate(a2face(i,j)%flux_h(unum,vnum))
  allocate(a2face(i,j)%flux_b(unum,vnum))
end do
end do

do i=q1imin-2,q1imax+2
do j=q2imin-2,q2imax+2
  call reduced_maxwell(h0,b0,uspace,vspace,ctr(i,j)%prim)

  ctr(i,j)%h = h0
  ctr(i,j)%b = b0
  ctr(i,j)%sh = 0.d0
  ctr(i,j)%sb = 0.d0
end do
end do

end subroutine initFlowfield

!-------------------------------------------------------------------

subroutine initIo

call writeMacro(iter)

open(unit=filehstid,file=hstfile,status="replace",action="write")

end subroutine initIo

!-------------------------------------------------------------------

end module initmod