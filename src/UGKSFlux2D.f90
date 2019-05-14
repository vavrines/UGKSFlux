
! flux.f90
!-------------------------------------------------------------------

module fluxmod

implicit none

real(kind=8), parameter :: PI = 3.141592654d0

integer,parameter :: MNUM = 6 !number of normal velocity moments
integer,parameter :: MTUM = 6 !number of tangential velocity moments

contains

!-------------------------------------------------------------------

subroutine flux_ugks2d(wL, hL, bL, shL, sbL, lenL, fluxw, fluxh, fluxb, lenFace, wR, hR, bR, shR, sbR, lenR, &
                       unum, vnum, uspace, vspace, weight, &
                       ink, gamma, muref, omega, prandtl, &
                       dt, dirc, cosa, sina)
    
real(kind=8), dimension(:), intent(in) :: wL, wR !// conservative variables
real(kind=8), dimension(:,:), intent(in) :: hL, bL, hR, bR !// distribution function
real(kind=8), dimension(:,:,:), intent(in) :: shL, sbL, shR, sbR !// slopes of distribution function
real(kind=8), dimension(:), intent(inout) :: fluxw !// interface fluxes of conservative variables
real(kind=8), dimension(:,:), intent(inout) :: fluxh, fluxb !// interface fluxes of distribution function
real(kind=8), intent(in) :: lenFace, lenL, lenR !// lengths
integer, intent(in) :: unum, vnum !// number of velocity grids
real(kind=8), dimension(:,:), intent(in) :: uspace, vspace, weight !// velocity quadrature points and weights
real(kind=8), intent(in) :: ink, gamma !// internal degrees of freedom of gas and Poisson ratio
real(kind=8), intent(in) :: muref, omega, prandtl !// reference viscosity, VHS model index, and Prandtl number
real(kind=8), intent(in) :: dt !// time step
integer,intent(in) :: dirc !// direction
real(kind=8), intent(in) :: cosa, sina !// directional cosine and sine

!Heaviside step function
integer, allocatable, dimension(:,:) :: delta

!Local micro velocity space
real(kind=8), allocatable, dimension(:,:) :: vn, vt

!interface variable
real(kind=8), allocatable, dimension(:,:) :: h, b
real(kind=8), allocatable, dimension(:,:) :: H0, B0
real(kind=8), allocatable, dimension(:,:) :: H_plus, B_plus
real(kind=8), allocatable, dimension(:,:) :: sh, sb 
real(kind=8) :: w(4), prim(4)
real(kind=8) :: qf(2)
real(kind=8) :: sw(4)
real(kind=8) :: aL(4), aR(4), aT(4)

!moments variable
real(kind=8) :: Mu(0:MNUM), Mu_L(0:MNUM), Mu_R(0:MNUM), Mv(0:MTUM), Mxi(0:2)
real(kind=8) :: Muv(4)
real(kind=8) :: Mau_L(4), Mau_R(4)
real(kind=8) :: Mau_T(4)
real(kind=8) :: tau
real(kind=8) :: Mt(5)

!--------------------------------------------------
! initialize
!--------------------------------------------------
allocate(vn(unum,vnum))
allocate(vt(unum,vnum))
allocate(delta(unum,vnum))
allocate(h(unum,vnum))
allocate(b(unum,vnum))
allocate(sh(unum,vnum))
allocate(sb(unum,vnum))
allocate(H0(unum,vnum))
allocate(B0(unum,vnum))
allocate(H_plus(unum,vnum))
allocate(B_plus(unum,vnum))

vn = uspace * cosa + vspace * sina
vt = vspace * cosa - uspace * sina

delta = (sign(1.d0,vn) + 1.d0) / 2.d0

!--------------------------------------------------
! reconstruct initial distribution
!--------------------------------------------------
h = (hL + 0.5 * lenL * shL(:,:,DIRC)) * delta + &
    (hR - 0.5 * lenR * shR(:,:,DIRC)) * (1.d0 - delta)
b = (bL + 0.5 * lenL * sbL(:,:,DIRC)) * delta + &
    (bR - 0.5 * lenR * sbR(:,:,DIRC)) * (1.d0 - delta)

sh = shL(:,:,DIRC) * delta + shR(:,:,DIRC) * (1.d0 - delta)
sb = sbL(:,:,DIRC) * delta + sbR(:,:,DIRC) * (1.d0 - delta)

!--------------------------------------------------
! obtain macroscopic variables at interface
!--------------------------------------------------
!conservative variables W_0 
w(1) = sum(weight * h)
w(2) = sum(weight * vn * h)
w(3) = sum(weight * vt * h)
w(4) = 0.5d0 * (sum(weight * (vn**2 + vt**2) * h) + sum(weight * b))
            
!convert to primary variables
prim = prim_variable(w, gamma)

!heat flux
qf = heat_flux(h, b, prim, vn, vt, weight) 

!--------------------------------------------------
! calculate a^L,a^R
!--------------------------------------------------
sw = (w - local_frame(wL, cosa, sina)) / (0.5 * lenL)
aL = micro_slope(prim, sw, ink)

sw = (local_frame(wR, cosa, sina) - w) / (0.5 * lenR)
aR = micro_slope(prim, sw, ink)

!--------------------------------------------------
! calculate time slope of W and A
!--------------------------------------------------
!<u^n>,<v^m>,<\xi^l>,<u^n>_{>0},<u^n>_{<0}
call calc_moment(prim, Mu, Mv, Mxi, Mu_L, Mu_R, ink) 

Mau_L = moment_au(aL, Mu_L, Mv, Mxi, 1, 0)
Mau_R = moment_au(aR, Mu_R, Mv, Mxi, 1, 0)

sw = -prim(1) * (Mau_L + Mau_R)
aT = micro_slope(prim, sw, ink)

!--------------------------------------------------
! calculate collision time and some time integration terms
!--------------------------------------------------
tau = collision_time(prim, muref, omega)

Mt(4) = tau*(1.d0 - exp(-dt / tau))
Mt(5) = -tau * dt * exp(-dt / tau) + tau * Mt(4)
Mt(1) = dt - Mt(4)
Mt(2) = -tau * Mt(1) + Mt(5)
Mt(3) = dt**2 / 2.d0 - tau * Mt(1)

!--------------------------------------------------
! calculate the flux of conservative variables related to g0
!--------------------------------------------------
Muv = moment_uv(Mu, Mv, Mxi, 1, 0, 0)
Mau_L = moment_au(aL, Mu_L, Mv, Mxi, 2, 0)
Mau_R = moment_au(aR, Mu_R, Mv, Mxi, 2, 0)
Mau_T = moment_au(aT, Mu, Mv, Mxi, 1, 0)

fluxw = Mt(1) * prim(1) * Muv + Mt(2) * prim(1) * (Mau_L + Mau_R) + Mt(3) * prim(1) * Mau_T

!--------------------------------------------------
! calculate the flux of conservative variables related to g+ and f0
!--------------------------------------------------
!Maxwellian distribution H0 and B0
call maxwell(H0, B0, prim, vn, vt, ink)

!Shakhov part H+ and B+
call shakhov(H0, B0, qf, prim, H_plus, B_plus, vn, vt, ink, prandtl)

!macro flux related to g+ and f0
fluxw(1) = fluxw(1) + Mt(1) * sum(weight * vn * H_plus) + Mt(4) * sum(weight * vn * h) - Mt(5) * sum(weight * vn**2 * sh)
fluxw(2) = fluxw(2) + Mt(1) * sum(weight * vn * vn * H_plus) + Mt(4) * sum(weight * vn * vn * h) - Mt(5) * sum(weight * vn * vn**2 * sh)
fluxw(3) = fluxw(3) + Mt(1) * sum(weight * vt * vn * H_plus) + Mt(4) * sum(weight * vt * vn * h) - Mt(5) * sum(weight * vt * vn**2 * sh)
fluxw(4) = fluxw(4) + &
           Mt(1) * 0.5d0 * (sum(weight * vn *(vn**2 + vt**2) * H_plus) + sum(weight * vn * B_plus)) + &
           Mt(4) * 0.5d0 * (sum(weight * vn *(vn**2 + vt**2) * h) + sum(weight * vn * b)) - &
           Mt(5) * 0.5d0 * (sum(weight * vn**2 * (vn**2 + vt**2) * sh) + sum(weight * vn**2 * sb))

!--------------------------------------------------
! calculate flux of distribution function
!--------------------------------------------------
fluxh = Mt(1) * vn * (H0 + H_plus) + &
        Mt(2) * vn**2 * (aL(1) * H0 + aL(2) * vn * H0 + aL(3) * vt * H0 + 0.5d0 * aL(4) * ((vn**2 + vt**2) * H0 + B0)) * delta + & 
        Mt(2) * vn**2 * (aR(1) * H0 + aR(2) * vn * H0 + aR(3) * vt * H0 + 0.5d0 * aR(4) * ((vn**2 + vt**2) * H0 + B0)) * (1.d0 - delta) + &
        Mt(3) * vn * (aT(1) * H0 + aT(2) * vn * H0 + aT(3) * vt * H0 + 0.5d0 * aT(4) * ((vn**2 + vt**2) * H0 + B0)) + &
        Mt(4) * vn * h - Mt(5) * vn**2 * sh

fluxb = Mt(1) * vn * (B0 + B_plus) + &
        Mt(2) * vn**2 * (aL(1) * B0 + aL(2) * vn * B0 + aL(3) * vt * B0 + 0.5d0 * aL(4) * ((vn**2 + vt**2) * B0 + Mxi(2) * H0)) * delta + &
        Mt(2) * vn**2 * (aR(1) * B0 + aR(2) * vn * B0 + aR(3) * vt * B0 + 0.5d0 * aR(4) * ((vn**2 + vt**2) * B0 + Mxi(2) * H0)) * (1.d0 - delta) + &
        Mt(3) * vn * (aT(1) * B0 + aT(2) * vn * B0 + aT(3) * vt * B0 + 0.5d0 * aT(4) * ((vn**2 + vt**2) * B0 + Mxi(2) * H0)) + &
        Mt(4) * vn * b - Mt(5) * vn**2 * sb

!--------------------------------------------------
! final flux
!--------------------------------------------------
! convert to global frame
fluxw = global_frame(fluxw, cosa, sina)

! total flux
fluxw = lenFace * fluxw 
fluxh = lenFace * fluxh
fluxb = lenFace * fluxb

end subroutine flux_ugks2d

!-------------------------------------------------------------------

subroutine flux_isothermal(bc, wCell, hCell, bCell, shCell, sbCell, lenCell, fluxw, fluxh, fluxb, lenFace, &
                           unum, vnum, uspace, vspace, weight, &
                           ink, gamma, muref, omega, &
                           dt, dirc, rot, cosa, sina)

real(kind=8), intent(in) :: bc(4) !// boundary condition
real(kind=8), dimension(:), intent(in) :: wCell !// conservative variables
real(kind=8), dimension(:,:), intent(in) :: hCell, bCell !// distribution function
real(kind=8), dimension(:,:,:), intent(in) :: shCell, sbCell !// slopes of distribution function
real(kind=8), dimension(:), intent(inout) :: fluxw !// interface fluxes of conservative variables
real(kind=8), dimension(:,:), intent(inout) :: fluxh, fluxb !// interface fluxes of distribution function
real(kind=8), intent(in) :: lenFace, lenCell !// lengths
integer, intent(in) :: unum, vnum !// number of velocity grids
real(kind=8), dimension(:,:), intent(in) :: uspace, vspace, weight !// velocity quadrature points and weights
real(kind=8), intent(in) :: ink, gamma !// internal degrees of freedom of gas and Poisson ratio
real(kind=8), intent(in) :: muref, omega !// reference viscosity, VHS model index
real(kind=8), intent(in) :: dt !// time step
integer,intent(in) :: dirc !// direction
integer,intent(in) :: rot !// frame rotation
real(kind=8), intent(in) :: cosa, sina !// directional cosine and sine

!Heaviside step function
integer, allocatable, dimension(:,:) :: delta 

!local micro velocity space
real(kind=8), allocatable, dimension(:,:) :: vn,vt 

!interface variable
real(kind=8), allocatable, dimension(:,:) :: h, b
real(kind=8), allocatable, dimension(:,:) :: H0, B0
real(kind=8) :: prim(4) 
real(kind=8) :: R1, R2

!--------------------------------------------------
! initialize
!--------------------------------------------------
!allocate array
allocate(vn(unum,vnum))
allocate(vt(unum,vnum))
allocate(delta(unum,vnum))
allocate(h(unum,vnum))
allocate(b(unum,vnum))
allocate(H0(unum,vnum))
allocate(B0(unum,vnum))

!convert the micro velocity to local frame
vn = uspace * cosa + vspace * sina
vt = vspace * cosa - uspace * sina

!Heaviside step function. The rotation accounts for the right wall
delta=(sign(1.d0, vn) * rot + 1.d0) / 2.d0

!boundary condition in local frame
prim = local_frame(bc, cosa, sina)

!--------------------------------------------------
!obtain h^{in} and b^{in}, rotation accounts for the right wall
!--------------------------------------------------
h = hCell - rot * 0.5d0 * lenCell * shCell(:,:,DIRC)
b = bCell - rot * 0.5d0 * lenCell * sbCell(:,:,DIRC)

!--------------------------------------------------
!calculate wall density and Maxwellian distribution
!--------------------------------------------------
R1 = sum(weight * vn * h * (1.d0 - delta))
R2 = (prim(4) / PI)*sum(weight * vn * exp(-prim(4) * ((vn - prim(2))**2 + (vt - prim(3))**2)) * delta)

prim(1) = -R1 / R2

call maxwell(H0, B0, prim, vn, vt, ink)

!--------------------------------------------------
!distribution function at the boundary interface
!--------------------------------------------------
h = H0 * delta + h * (1.d0 - delta)
b = B0 * delta + b * (1.d0 - delta)

!--------------------------------------------------
!calculate flux
!--------------------------------------------------
fluxw(1) = sum(weight * vn * h)
fluxw(2) = sum(weight * vn * vn * h)
fluxw(3) = sum(weight * vn * vt * h)
fluxw(4) = 0.5d0 * sum(weight * vn * ((vn**2 + vt**2) * h + b))

fluxh = vn * h
fluxb = vn * b

!--------------------------------------------------
!final flux
!--------------------------------------------------
!convert to global frame
fluxw = global_frame(fluxw, cosa, sina) 

!total flux
fluxw = dt * lenFace * fluxw
fluxh = dt * lenFace * fluxh
fluxb = dt * lenFace * fluxb

end subroutine flux_isothermal

!-------------------------------------------------------------------

function global_frame(w, cosx, cosy)

    real(kind=8), intent(in) :: w(4)
    real(kind=8), intent(in) :: cosx, cosy
    real(kind=8) :: global_frame(4)

    global_frame(1) = w(1)
    global_frame(2) = w(2) * cosx - w(3) * cosy
    global_frame(3) = w(2) * cosy + w(3) * cosx
    global_frame(4) = w(4)
    
end function global_frame


function local_frame(w, cosx, cosy)

    real(kind=8), intent(in) :: w(4)
    real(kind=8), intent(in) :: cosx, cosy
    real(kind=8) :: local_frame(4)

    local_frame(1) = w(1)
    local_frame(2) = w(2) * cosx + w(3) * cosy
    local_frame(3) = w(3) * cosx - w(2) * cosy
    local_frame(4) = w(4)
    
end function local_frame

!-------------------------------------------------------------------

function prim_variable(w, gamma)

    real(kind=8), intent(in) :: w(4), gamma
    real(kind=8) :: prim_variable(4)

    prim_variable(1) = w(1)
    prim_variable(2) = w(2) / w(1)
    prim_variable(3) = w(3) / w(1)
    prim_variable(4) = 0.5d0 * w(1) / (gamma - 1.d0) / (w(4) - 0.5d0 * (w(2)**2 + w(3)**2) / w(1))

end function prim_variable

!-------------------------------------------------------------------

subroutine maxwell(h, b, prim, vn, vt, ink)

    real(kind=8), dimension(:,:), intent(out) :: h, b
    real(kind=8), dimension(:,:), intent(in) :: vn, vt
    real(kind=8), intent(in) :: prim(4), ink

    h = prim(1) * (prim(4) / PI) * exp(-prim(4) * ((vn - prim(2))**2 + (vt - prim(3))**2))
    b = h * ink / (2.d0 * prim(4))
    
end subroutine maxwell

!-------------------------------------------------------------------

subroutine shakhov(H, B, qf, prim, H_plus, B_plus, vn, vt, ink, pr)
            
    real(kind=8), dimension(:,:), intent(in) :: H, B
    real(kind=8), intent(in) :: qf(2)
    real(kind=8), intent(in) :: prim(4)
    real(kind=8), dimension(:,:), intent(out) :: H_plus, B_plus
    real(kind=8), dimension(:,:), intent(in) :: vn, vt
    real(kind=8), intent(in) :: ink, pr

    H_plus = 0.8d0 * (1.d0 - pr) * prim(4)**2 / prim(1) * &
             ((vn - prim(2)) * qf(1) + (vt - prim(3)) * qf(2)) * (2.d0 * prim(4) * ((vn - prim(2))**2 + (vt - prim(3))**2) + ink - 5.d0) * H
    B_plus = 0.8d0 * (1.d0 - pr) * prim(4)**2 / prim(1) * &
             ((vn - prim(2)) * qf(1) + (vt - prim(3)) * qf(2)) * (2.d0 * prim(4) * ((vn - prim(2))**2 + (vt - prim(3))**2) + ink - 3.d0) * B
            
end subroutine shakhov

!-------------------------------------------------------------------

function collision_time(prim, muref, omega)

    real(kind=8), intent(in) :: prim(4), muref, omega
    real(kind=8) :: collision_time

    collision_time = muref * 2.d0 * prim(4)**(1.d0 - omega) / prim(1)

end function collision_time

!-------------------------------------------------------------------

function heat_flux(h, b, prim, vn, vt, weight)

    real(kind=8), dimension(:,:), intent(in) :: h, b
    real(kind=8), dimension(:,:), intent(in) :: vn, vt, weight
    real(kind=8), intent(in) :: prim(4)
    real(kind=8) :: heat_flux(2)

    heat_flux(1) = 0.5d0 * (sum(weight * (vn - prim(2)) * ((vn - prim(2))**2 + (vt - prim(3))**2) * h) + sum(weight * (vn - prim(2)) * b))
    heat_flux(2) = 0.5d0 * (sum(weight * (vt - prim(3)) * ((vn - prim(2))**2 + (vt - prim(3))**2) * h) + sum(weight * (vt - prim(3)) * b))
    
end function heat_flux

!-------------------------------------------------------------------

function micro_slope(prim, sw, ink)
            
    real(kind=8), intent(in) :: prim(4), sw(4), ink
    real(kind=8) :: micro_slope(4)

    micro_slope(4) = 4.d0 * prim(4)**2 / (ink + 2.d0) / prim(1) * (2.d0 * sw(4) - 2.d0 * prim(2) * sw(2) - 2.d0 * prim(3) * sw(3) + sw(1) * (prim(2)**2 + prim(3)**2 - 0.5d0 * (ink + 2.d0) / prim(4)))
    micro_slope(3) = 2.d0 * prim(4) / prim(1) * (sw(3) - prim(3) * sw(1)) - prim(3) * micro_slope(4)
    micro_slope(2) = 2.d0 * prim(4) / prim(1) * (sw(2) - prim(2) * sw(1)) - prim(2) * micro_slope(4)
    micro_slope(1) = sw(1) / prim(1) - prim(2) * micro_slope(2) - prim(3) * micro_slope(3) - 0.5d0 * (prim(2)**2 + prim(3)**2 + 0.5d0 * (ink + 2.d0) / prim(4)) * micro_slope(4)

end function micro_slope

!-------------------------------------------------------------------

subroutine calc_moment(prim, Mu, Mv, Mxi, Mu_L, Mu_R, ink)

    real(kind=8), intent(in) :: prim(4)
    real(kind=8), intent(out) :: Mu(0:MNUM), Mu_L(0:MNUM), Mu_R(0:MNUM)
    real(kind=8), intent(out) :: Mv(0:MTUM)
    real(kind=8), intent(out) :: Mxi(0:2)
    real(kind=8), intent(in) :: ink
    integer :: i

    !moments of normal velocity
    Mu_L(0) = 0.5d0 * erfc(-sqrt(prim(4)) * prim(2))
    Mu_L(1) = prim(2) * Mu_L(0) + 0.5d0 * exp(-prim(4) * prim(2)**2) / sqrt(PI * prim(4))
    Mu_R(0) = 0.5d0 * erfc(sqrt(prim(4)) * prim(2))
    Mu_R(1) = prim(2) * Mu_R(0) - 0.5d0 * exp(-prim(4) * prim(2)**2) / sqrt(PI * prim(4))

    do i=2,MNUM
        Mu_L(i) = prim(2) * Mu_L(i-1) + 0.5d0 * (i-1) * Mu_L(i-2) / prim(4)
        Mu_R(i) = prim(2) * Mu_R(i-1) + 0.5d0 * (i-1) * Mu_R(i-2) / prim(4)
    end do

    Mu = Mu_L + Mu_R

    !moments of tangential velocity
    Mv(0) = 1.0
    Mv(1) = prim(3)

    do i=2,MTUM
        Mv(i) = prim(3) * Mv(i-1) + 0.5d0 * (i-1) * Mv(i-2) / prim(4)
    end do

    !moments of \xi
    Mxi(0) = 1.d0 !<\xi^0>
    Mxi(1) = 0.5d0 * ink / prim(4) !<\xi^2>
    Mxi(2) = (ink**2 + 2.d0 * ink) / (4.d0 * prim(4)**2) !<\xi^4>
    
end subroutine calc_moment

!-------------------------------------------------------------------

function moment_uv(Mu,Mv,Mxi,alpha,beta,delta)

    real(kind=8), intent(in) :: Mu(0:MNUM), Mv(0:MTUM), Mxi(0:2)
    integer, intent(in) :: alpha, beta, delta
    real(kind=8) :: moment_uv(4)

    moment_uv(1) = Mu(alpha) * Mv(beta) * Mxi(delta/2)
    moment_uv(2) = Mu(alpha+1) * Mv(beta) * Mxi(delta/2)
    moment_uv(3) = Mu(alpha) * Mv(beta+1) * Mxi(delta/2)
    moment_uv(4) = 0.5d0 * (Mu(alpha+2) * Mv(beta) * Mxi(delta/2) + Mu(alpha) * Mv(beta+2) * Mxi(delta/2) + Mu(alpha) * Mv(beta) * Mxi((delta+2)/2))

end function moment_uv

!-------------------------------------------------------------------

function moment_au(a,Mu,Mv,Mxi,alpha,beta)
    
    real(kind=8), intent(in) :: a(4)
    real(kind=8), intent(in) :: Mu(0:MNUM), Mv(0:MTUM), Mxi(0:2)
    integer, intent(in) :: alpha, beta
    real(kind=8) :: moment_au(4)

    moment_au = a(1) * moment_uv(Mu,Mv,Mxi,alpha+0,beta+0,0) + &
                a(2) * moment_uv(Mu,Mv,Mxi,alpha+1,beta+0,0) + &
                a(3) * moment_uv(Mu,Mv,Mxi,alpha+0,beta+1,0) + &
                0.5d0 * a(4) * moment_uv(Mu,Mv,Mxi,alpha+2,beta+0,0) + &
                0.5d0 * a(4) * moment_uv(Mu,Mv,Mxi,alpha+0,beta+2,0) + &
                0.5d0 * a(4) * moment_uv(Mu,Mv,Mxi,alpha+0,beta+0,2)
    
end function moment_au

!-------------------------------------------------------------------

end module fluxmod