
! kinetic.f90
!-------------------------------------------------------------------

module kineticmod

use datamod

implicit none

real(srk),pointer :: duh(:,:)
real(srk),pointer :: dub(:,:)
real(srk),pointer :: dvh(:,:)
real(srk),pointer :: dvb(:,:)

contains

!-------------------------------------------------------------------

subroutine reduced_maxwell(H,B,vn,vt,prim)

real(srk),dimension(:,:),intent(out) :: H,B
real(srk),dimension(:,:),intent(in) :: vn,vt
real(srk),dimension(:,:),intent(in) :: prim(dim+2)

H = prim(1)*(prim(4)/PI)*exp(-prim(4)*((vn-prim(2))**2+(vt-prim(3))**2))
B = H*inK/(2.0*prim(4))

end subroutine reduced_maxwell

!-------------------------------------------------------------------

subroutine shakhov_part(H,B,vn,vt,qf,prim,H_plus,B_plus)

real(srk),dimension(:,:),intent(in) :: H,B
real(srk),dimension(:,:),intent(in) :: vn,vt
real(srk),intent(in) :: qf(dim)
real(srk),intent(in) :: prim(dim+2)
real(srk),dimension(:,:),intent(out) :: H_plus,B_plus

H_plus = 0.8*(1-prandtl)*prim(4)**2/prim(1)*&
         ((vn-prim(2))*qf(1)+(vt-prim(3))*qf(2))*(2*prim(4)*((vn-prim(2))**2+(vt-prim(3))**2)+inK-5)*H
B_plus = 0.8*(1-prandtl)*prim(4)**2/prim(1)*&
         ((vn-prim(2))*qf(1)+(vt-prim(3))*qf(2))*(2*prim(4)*((vn-prim(2))**2+(vt-prim(3))**2)+inK-3)*B

end subroutine shakhov_part

!-------------------------------------------------------------------

function get_moment_conserved(h,b,vn,vt) result(w)

real(srk), intent(in), dimension(:,:) :: vn, vt
real(srk), dimension(:,:), intent(in) :: h,b
real(srk) :: w(dim+2)

w(1) = sum(weight*h)
w(2) = sum(weight*vn*h)
w(3) = sum(weight*vt*h)
w(4) = 0.5*(sum(weight*(vn**2+vt**2)*h)+sum(weight*b))

end function get_moment_conserved

!-------------------------------------------------------------------

function get_conserved(prim)

real(srk),intent(in) :: prim(dim+2)
real(srk) :: get_conserved(dim+2)

get_conserved(1) = prim(1)
get_conserved(2) = prim(1)*prim(2)
get_conserved(3) = prim(1)*prim(3)
get_conserved(4) = 0.5*prim(1)/prim(4)/(gamma-1.0)+0.5*prim(1)*(prim(2)**2+prim(3)**2)

end function get_conserved

 
function get_primitive(w)

real(srk),intent(in) :: w(dim+2)
real(srk) :: get_primitive(dim+2)

get_primitive(1) = w(1)
get_primitive(2) = w(2)/w(1)
get_primitive(3) = w(3)/w(1)
get_primitive(4) = 0.5*w(1)/(gamma-1.0)/(w(4)-0.5*(w(2)**2+w(3)**2)/w(1))

end function get_primitive

!-------------------------------------------------------------------

function get_gamma(K)

real(srk),intent(in) :: K
real(srk) :: get_gamma

get_gamma = (K+4.0d0)/(K+2.0d0)

end function get_gamma

!-------------------------------------------------------------------

function get_sos(lam)
real(srk),intent(in) :: lam
real(srk) :: get_sos 

get_sos = sqrt(0.5d0*gamma/lam)

end function get_sos

!-------------------------------------------------------------------

function get_tau(prim)

real(srk),intent(in) :: prim(dim+2)
real(srk) :: get_tau

get_tau = muref*2*prim(4)**(1-omega)/prim(1)

end function get_tau

!-------------------------------------------------------------------

! VSS model

function get_mu(kn,alpha,omega)

real(srk),intent(in) :: kn,alpha,omega
real(srk) :: get_mu
    
get_mu = 5.0d0*(alpha+1.0d0)*(alpha+2.0d0)*sqrt(PI)/&
         (4.d0*alpha*(5.d0-2.d0*omega)*(7.d0-2.d0*omega))*kn

end function get_mu

!-------------------------------------------------------------------

function get_heat_flux(h,b,vn,vt,prim)

real(srk),dimension(:,:),intent(in) :: h,b
real(srk),dimension(:,:),intent(in) :: vn,vt
real(srk),intent(in) :: prim(dim+2)
real(srk) :: get_heat_flux(dim)

get_heat_flux(1) = 0.5*(sum(weight*(vn-prim(2))*((vn-prim(2))**2+(vt-prim(3))**2)*h)+sum(weight*(vn-prim(2))*b)) 
get_heat_flux(2) = 0.5*(sum(weight*(vt-prim(3))*((vn-prim(2))**2+(vt-prim(3))**2)*h)+sum(weight*(vt-prim(3))*b)) 

end function get_heat_flux

!-------------------------------------------------------------------

function get_temperature(h,b,vn,vt,prim)

real(srk),dimension(:,:),intent(in) :: h,b
real(srk),dimension(:,:),intent(in) :: vn,vt
real(srk),intent(in) :: prim(dim+2)
real(srk) :: get_temperature

get_temperature = 2.0*(sum(weight*((vn-prim(2))**2+(vt-prim(3))**2)*h)+sum(weight*b))/(inK+2)/prim(1)

end function get_temperature

!-------------------------------------------------------------------

function get_pressure(h,b,vn,vt,prim)

real(srk),dimension(:,:),intent(in) :: h,b
real(srk),dimension(:,:),intent(in) :: vn,vt
real(srk),intent(in) :: prim(dim+2)
real(srk) :: get_pressure

get_pressure = (sum(weight*((vn-prim(2))**2+(vt-prim(3))**2)*h)+sum(weight*b))/(inK+2)

end function get_pressure

!------------------------------------------------------------------- 

function get_stress(h,b,vn,vt,prim)

real(srk),dimension(:,:),intent(in) :: h,b
real(srk),dimension(:,:),intent(in) :: vn,vt
real(srk),intent(in) :: prim(dim+2)
real(srk) :: get_stress(3)

get_stress(1) = sum(weight*(vn-prim(2))*(vn-prim(2))*h)
get_stress(2) = sum(weight*(vn-prim(2))*(vt-prim(3))*h)
get_stress(3) = sum(weight*(vt-prim(3))*(vt-prim(3))*h)

end function get_stress

!------------------------------------------------------------------- 

function micro_slope(prim,sw)

real(srk),intent(in) :: prim(dim+2),sw(dim+2)
real(srk) :: micro_slope(dim+2)

micro_slope(4) = 4.0*prim(4)**2/(inK+2)/prim(1)*(2.0*sw(4)-2.0*prim(2)*sw(2)-&
                 2.0*prim(3)*sw(3)+sw(1)*(prim(2)**2+prim(3)**2-0.5*(inK+2)/prim(4)))

micro_slope(3) = 2.0*prim(4)/prim(1)*(sw(3)-prim(3)*sw(1))-prim(3)*micro_slope(4)
micro_slope(2) = 2.0*prim(4)/prim(1)*(sw(2)-prim(2)*sw(1))-prim(2)*micro_slope(4)
micro_slope(1) = sw(1)/prim(1)-prim(2)*micro_slope(2)-prim(3)*micro_slope(3)-&
                 0.5*(prim(2)**2+prim(3)**2+0.5*(inK+2)/prim(4))*micro_slope(4)

end function micro_slope

!-------------------------------------------------------------------
    
subroutine calc_moment(prim,Mu,Mv,Mxi,Mu_L,Mu_R)

real(srk),intent(in) :: prim(dim+2)
real(srk),intent(out) :: Mu(0:6),Mu_L(0:6),Mu_R(0:6)
real(srk),intent(out) :: Mv(0:6)
real(srk),intent(out) :: Mxi(0:2)
integer :: i

!moments of normal velocity
Mu_L(0) = 0.5*erfc(-sqrt(prim(4))*prim(2))
Mu_L(1) = prim(2)*Mu_L(0)+0.5*exp(-prim(4)*prim(2)**2)/sqrt(PI*prim(4))
Mu_R(0) = 0.5*erfc(sqrt(prim(4))*prim(2))
Mu_R(1) = prim(2)*Mu_R(0)-0.5*exp(-prim(4)*prim(2)**2)/sqrt(PI*prim(4))

do i=2,6
    Mu_L(i) = prim(2)*Mu_L(i-1)+0.5*(i-1)*Mu_L(i-2)/prim(4)
    Mu_R(i) = prim(2)*Mu_R(i-1)+0.5*(i-1)*Mu_R(i-2)/prim(4)
end do

Mu = Mu_L+Mu_R

!moments of tangential velocity
Mv(0) = 1.0
Mv(1) = prim(3)

do i=2,6
    Mv(i) = prim(3)*Mv(i-1)+0.5*(i-1)*Mv(i-2)/prim(4)
end do

!moments of \xi
Mxi(0) = 1.0 !<\xi^0>
Mxi(1) = 0.5*inK/prim(4) !<\xi^2>
Mxi(2) = (inK**2+2.0*inK)/(4.0*prim(4)**2) !<\xi^4>

end subroutine calc_moment

!-------------------------------------------------------------------

function moment_uv(Mu,Mv,Mxi,alpha,beta,theta)

real(srk),intent(in) :: Mu(0:6),Mv(0:6),Mxi(0:2)
integer,intent(in) :: alpha,beta,theta
real(srk) :: moment_uv(dim+2)

moment_uv(1) = Mu(alpha)*Mv(beta)*Mxi(theta/2)
moment_uv(2) = Mu(alpha+1)*Mv(beta)*Mxi(theta/2)
moment_uv(3) = Mu(alpha)*Mv(beta+1)*Mxi(theta/2)
moment_uv(4) = 0.5*(Mu(alpha+2)*Mv(beta)*Mxi(theta/2)+&
               Mu(alpha)*Mv(beta+2)*Mxi(theta/2)+Mu(alpha)*Mv(beta)*Mxi((theta+2)/2))

end function moment_uv


function moment_au(a,Mu,Mv,Mxi,alpha,beta)
real(srk),intent(in) :: a(dim+2)
real(srk),intent(in) :: Mu(0:6),Mv(0:6),Mxi(0:2)
integer,intent(in) :: alpha,beta
real(srk) :: moment_au(dim+2)

moment_au = a(1)*moment_uv(Mu,Mv,Mxi,alpha+0,beta+0,0)+&
            a(2)*moment_uv(Mu,Mv,Mxi,alpha+1,beta+0,0)+&
            a(3)*moment_uv(Mu,Mv,Mxi,alpha+0,beta+1,0)+&
            0.5*a(4)*moment_uv(Mu,Mv,Mxi,alpha+2,beta+0,0)+&
            0.5*a(4)*moment_uv(Mu,Mv,Mxi,alpha+0,beta+2,0)+&
            0.5*a(4)*moment_uv(Mu,Mv,Mxi,alpha+0,beta+0,2)

end function moment_au

!-------------------------------------------------------------------

subroutine cal_fdv(cell)

type(cell_center),intent(in) :: cell
integer :: i,j
real(srk) :: sL,sR  

if(associated(duh)) deallocate(duh)
allocate(duh(unum,vnum))
if(associated(dub)) deallocate(dub)
allocate(dub(unum,vnum))
if(associated(dvh)) deallocate(dvh)
allocate(dvh(unum,vnum))
if(associated(dvb)) deallocate(dvb)
allocate(dvb(unum,vnum))

!boundary
duh(1,:) = (cell%h(2,:)-cell%h(1,:))/(uspace(2,1)-uspace(1,1))
dub(1,:) = (cell%b(2,:)-cell%b(1,:))/(uspace(2,1)-uspace(1,1))
dvh(:,1) = (cell%h(:,2)-cell%h(:,1))/(vspace(1,2)-vspace(1,1))
dvb(:,1) = (cell%b(:,2)-cell%b(:,1))/(vspace(1,2)-vspace(1,1))
duh(unum,:) = (cell%h(unum,:)-cell%h(unum-1,:))/(uspace(unum,1)-uspace(unum-1,1))
dub(unum,:) = (cell%b(unum,:)-cell%b(unum-1,:))/(uspace(unum,1)-uspace(unum-1,1))
dvh(:,vnum) = (cell%h(:,vnum)-cell%h(:,vnum-1))/(vspace(1,vnum)-vspace(1,vnum-1))
dvb(:,vnum) = (cell%b(:,vnum)-cell%b(:,vnum-1))/(vspace(1,vnum)-vspace(1,vnum-1))

!inner
do i=2,unum-1
do j=1,vnum
    sl=(cell%h(i,j)-cell%h(i-1,j))/(uspace(i,1)-uspace(i-1,1))
    sr=(cell%h(i+1,j)-cell%h(i,j))/(uspace(i+1,1)-uspace(i,1))
    !duh(i,j)=sr
    duh(i,j)=0.5*(sl+sr)

    sl=(cell%b(i,j)-cell%b(i-1,j))/(uspace(i,1)-uspace(i-1,1))
    sr=(cell%b(i+1,j)-cell%b(i,j))/(uspace(i+1,1)-uspace(i,1))
    !dub(i,j)=sr
    dub(i,j)=0.5*(sl+sr)
end do
end do

do i=1,unum
do j=2,vnum-1
    sl=(cell%h(i,j)-cell%h(i,j-1))/(vspace(1,j)-vspace(1,j-1))
    sr=(cell%h(i,j+1)-cell%h(i,j))/(vspace(1,j+1)-vspace(1,j))
    !dvh(i,j)=sr
    dvh(i,j)=0.5*(sl+sr)
    
    sl=(cell%b(i,j)-cell%b(i,j-1))/(vspace(1,j)-vspace(1,j-1))
    sr=(cell%b(i,j+1)-cell%b(i,j))/(vspace(1,j+1)-vspace(1,j))
    !dvb(i,j)=sr
    dvb(i,j)=0.5*(sl+sr)
end do
end do

end subroutine cal_fdv

!-------------------------------------------------------------------

end module kineticmod