
! kinetic.f90
!-------------------------------------------------------------------

module kineticmod

use datamod

implicit none

real(srk),pointer :: duh(:)
real(srk),pointer :: dub(:)

contains

!-------------------------------------------------------------------

subroutine reduced_maxwell(H,B,prim)

real(srk),dimension(:),intent(out) :: H,B
real(srk),intent(in) :: prim(DIM+2)

H = prim(1)*(prim(3)/PI)**(1.0/2.0)*exp(-prim(3)*(uspace-prim(2))**2)
B = h*ink/(2.0*prim(3))

end subroutine reduced_maxwell

!-------------------------------------------------------------------

subroutine shakhov_part(H,B,qf,prim,H_plus,B_plus)

real(srk),dimension(:),intent(in) :: H,B
real(srk),intent(in) :: qf
real(srk),intent(in) :: prim(DIM+2)
real(srk),dimension(:),intent(out) :: H_plus,B_plus

H_plus = 0.8*(1-prandtl)*prim(3)**2/prim(1)*&
         (uspace-prim(2))*qf*(2*prim(3)*(uspace-prim(2))**2+ink-5)*H
B_plus = 0.8*(1-prandtl)*prim(3)**2/prim(1)*&
         (uspace-prim(2))*qf*(2*prim(3)*(uspace-prim(2))**2+ink-3)*B

end subroutine shakhov_part

!-------------------------------------------------------------------

function get_moment_conserved(h,b) result(w)

real(srk), dimension(:), intent(in) :: h,b
real(srk) :: w(dim+2)
integer :: i

w(1) = sum(weight*h)
w(2) = sum(weight*uspace*h)
w(3) = 0.5*(sum(weight*uspace**2*h)+sum(weight*b))

end function get_moment_conserved

!-------------------------------------------------------------------

function get_conserved(prim)

real(srk),intent(in) :: prim(DIM+2)
real(srk) :: get_conserved(DIM+2)

get_conserved(1) = prim(1)
get_conserved(2) = prim(1)*prim(2)
get_conserved(3) = 0.5d0*prim(1)/prim(3)/(gamma-1.0)+0.5*prim(1)*prim(2)**2

end function get_conserved


function get_primitive(w)

real(srk),intent(in) :: w(DIM+2)
real(srk) :: get_primitive(DIM+2)

get_primitive(1) = w(1)
get_primitive(2) = w(2)/w(1)
get_primitive(3) = 0.5*w(1)/(gamma-1.0)/(w(3)-0.5*w(2)**2/w(1))

end function get_primitive

!-------------------------------------------------------------------

function get_gamma(K)

real(srk),intent(in) :: K
real(srk) :: get_gamma

get_gamma = (K+3.0d0)/(K+1.0d0)

end function get_gamma

!-------------------------------------------------------------------

function get_sos(lam)

real(srk),intent(in) :: lam
real(srk) :: get_sos 

get_sos = sqrt(0.5d0*gamma/lam)

end function get_sos

!-------------------------------------------------------------------

! calculate collision time: VHS model
function get_tau(prim)
real(srk),intent(in) :: prim(DIM+2)
real(srk) :: get_tau

!get_tau = muref/prim(1)
get_tau = muref*2*prim(3)**(1-omega)/prim(1)
!get_tau = 0.01*dt

end function get_tau

!-------------------------------------------------------------------

! calculate nondimensionalized viscosity coefficient
! alpha,omega :index related to HS/VHS/VSS model
function get_mu(kn,alpha,omega)
real(srk),intent(in) :: kn,alpha,omega
real(srk) :: get_mu

!get_mu = 0.5d0*sqrt(PI)*kn

! Hard Sphere model(HS model)
get_mu = 5.0d0*(alpha+1.0d0)*(alpha+2.0d0)*sqrt(PI)/&
        (4*alpha*(5-2*omega)*(7-2*omega))*kn

end function get_mu

! ------------------------------------------------------------

function get_heat_flux(h,b,prim)
real(srk),dimension(:),intent(in) :: h,b
real(srk),intent(in) :: prim(DIM+2)
real(srk) :: get_heat_flux

get_heat_flux = 0.5d0*(sum(weight*(uspace-prim(2))*(uspace-prim(2))**2*h)+sum(weight*(uspace-prim(2))*b))

end function get_heat_flux

!-------------------------------------------------------------------

subroutine cal_fdv(cell)

type(cell_center),intent(in) :: cell
real(srk) :: du
integer :: i

if(associated(duh)) deallocate(duh)
allocate(duh(unum))
if(associated(dub)) deallocate(dub)
allocate(dub(unum))

du = lengthU/(unum-1)

duh(1) = (cell%h(2)-cell%h(1))/du
dub(1) = (cell%b(2)-cell%b(1))/du
duh(unum) = (cell%h(unum)-cell%h(unum-1))/du
dub(unum) = (cell%b(unum)-cell%b(unum-1))/du

do i=2,unum-1
    duh(i)=(cell%h(i+1)-cell%h(i))/du
    dub(i)=(cell%b(i+1)-cell%b(i))/du
end do

end subroutine cal_fdv

!-------------------------------------------------------------------

end module kineticmod