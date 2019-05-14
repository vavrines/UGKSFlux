
! initialize.f90
!-------------------------------------------------------------------

module initmod

use datamod

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

call initGas

call initVelspace(velref,lambdaref,8.d0)

call initBoundary

call initExforce

call initFlowfield

call initIo

end subroutine initNew

!-------------------------------------------------------------------

subroutine initPara

grid_method = unigrid
interp_method = secondorder

cfl = 0.1d0
simtime = 0.d0
maxtime = 0.2d0
iter = 1

write(*,*) "CFL number: ",cfl

end subroutine initPara

!-------------------------------------------------------------------

subroutine initGeometry

integer :: i
real(srk) :: dq1

xlength = 1.d0
q1imin = 1
q1imax = 100
dq1 = xlength/(q1imax-q1imin+1)

write(*,*) "xlength: ",xlength
write(*,*) "dx: ",dq1

! allocate
allocate(ctr(q1imin-2:q1imax+2))
allocate(face(q1imin-2:q1imax+3))
allocate(pot(q1imin-2:q1imax+3)) !external potential on nodes

! geometry
do i=q1imin-2,q1imax+2
  ctr(i)%x = (i-0.5d0)*dq1
  ctr(i)%length = dq1
  ctr(i)%rgtype = continuum
end do
!ctr(q1imax+1:q1imax+2)%rgtype = rarefied
!ctr(q1imin-2:q1imin-1)%rgtype = rarefied

do i=q1imin-2,q1imax+3
  face(i)%bcType = 1
end do

end subroutine initGeometry

!-------------------------------------------------------------------

subroutine initGas

! Argon gas
prandtl = 2.d0/3.d0
knudsen = 0.0001d0
mach = 8.d0

namass = 39.95d0
molediameter = 4.17d-10
molemass = 6.63d-26
gasR = R0/namass*1.d3

! Hard sphere model
inK = 2.d0
gamma = get_gamma(inK)

omega = 0.81d0
alpha_ref = 1.d0
omega_ref = 0.5d0
muref = get_mu(knudsen,alpha_ref,omega_ref)

write(*,*) "Prandtl number: ",prandtl
write(*,*) "Knudsen number: ",knudsen
write(*,*) "gamma: ",gamma
write(*,*) "dynamic viscosity: ",muref

! Reference state
rhoref = 1.d0
velref = 0.d0
Tref = 1.d0
lambdaref = 1.d0/Tref
sosref = 1.d0

momref = rhoref*velref
Eref = 0.5d0*momref*momref/rhoref+&
	     0.5d0*rhoref/lambdaref/(gamma-1.0)

end subroutine initGas

!-------------------------------------------------------------------

subroutine initVelspace(vel,lambda,alpha)

real(srk),intent(in) :: vel,lambda,alpha
real(srk) :: referVel

lengthU = 2.d0*alpha*dsqrt(0.5d0*gamma/lambda)
centerU = vel

umax = centerU+0.5d0*lengthU
umin = centerU-0.5d0*lengthU

call initNewton

end subroutine initVelspace

!-------------------------------------------------------------------
	
subroutine initNewton

integer :: i,j

unum = 101

write(*,*) "Newton-Cotes velocity space"
write(*,*) "velocity points: ",unum
write(*,*) "U range: ",umin,umax

allocate(uspace(unum))
allocate(weight(unum))

du = (umax-umin)/(unum-1)

forall(i=1:unum)
  uspace(i) = umin+(i-1)*du
  weight(i) = (newton_coeff(i,unum)*du)
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

subroutine initBoundary

bcL = [1.d0,0.d0,0.5d0]
bcR = [1.d0,0.d0,0.5d0]

end subroutine initBoundary

!-------------------------------------------------------------------

subroutine initExforce

integer :: i
real(srk) :: dq1

phi = -0.d0
dq1 = xlength/(q1imax-q1imin+1)

do i=q1imin-2,q1imax+3
  pot(i) = -phi*dq1*(i-q1imin)
end do

write(*,*) "external acceleration: ",phi

end subroutine initExforce

!-------------------------------------------------------------------

subroutine initFlowfield

real(srk),allocatable,dimension(:) :: h0,b0
integer :: i
real(srk) :: ratio_T
real(srk) :: ma_R
real(srk) :: prim_L(dim+2),prim_R(dim+2)

ratio_T = (1+(gamma-1)/2*Mach**2)*(2*gamma/(gamma-1)*Mach**2-1)/(Mach**2*(2*gamma/(gamma-1)+(gamma-1)/2))
Ma_R = sqrt((Mach**2*(gamma-1)+2)/(2*gamma*Mach**2-(gamma-1)))

prim_L(1) = 1.0 !density
prim_L(2) = Mach*sqrt(gamma/2.0) !velocity
prim_L(3) = 1.0 !lambda=1/temperature

prim_R(1) = prim_L(1)*(gamma+1)*Mach**2/((gamma-1)*Mach**2+2)
prim_R(2) = Ma_R*sqrt(gamma/2.0)*sqrt(ratio_T)
prim_R(3) = prim_L(3)/ratio_T

!do i=q1imin-2,q1imax+2
!  if(i<=q1imax/2) then
!    ctr(i)%prim = prim_L
!  else
!    ctr(i)%prim = prim_R
!  end if

!  ctr(i)%w = get_conserved(ctr(i)%prim)
!end do


do i=q1imin-2,q1imax+2
  if(i<=q1imax/2) then
    ctr(i)%prim(1) = 1.d0
    ctr(i)%prim(2) = 0.d0
    ctr(i)%prim(3) = 0.5d0
  else
    ctr(i)%prim(1) = 0.125d0
    ctr(i)%prim(2) = 0.d0
    ctr(i)%prim(3) = 1.d0/1.6d0
  end if

  ctr(i)%w = get_conserved(ctr(i)%prim)
end do


allocate(h0(unum))
allocate(b0(unum))

do i=q1imin-2,q1imax+2
  allocate(ctr(i)%h(unum))
  allocate(ctr(i)%b(unum))
  allocate(ctr(i)%sh(unum))
  allocate(ctr(i)%sb(unum))
end do

do i=q1imin-2,q1imax+3
    allocate(face(i)%flux_h(unum))
    allocate(face(i)%flux_b(unum))
  end do

do i=q1imin-2,q1imax+2
  call reduced_maxwell(h0,b0,ctr(i)%prim)
  ctr(i)%h = h0
  ctr(i)%b = b0
  ctr(i)%sh = 0.d0
  ctr(i)%sb = 0.d0
end do

end subroutine initFlowfield

!-------------------------------------------------------------------

subroutine initIo

call writeMacro(iter)

open(unit=filehstid,file=hstfile,status="replace",action="write")

end subroutine initIo

!-------------------------------------------------------------------

end module initmod
