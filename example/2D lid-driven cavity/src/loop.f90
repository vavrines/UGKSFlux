
! loop.f90
!-------------------------------------------------------------------

module loopmod

use datamod

use gridmod

use mathmod

use kineticmod

use fluxmod

use postmod

use bcmod

implicit none

contains

!-------------------------------------------------------------------

subroutine loopjob

write(*,*) "main loop"

do while(.true.)

    call timestep

    call interpolation

    call evolution

    call update

    !if(simtime >= maxtime) exit

    if(all(res<eps)) exit

    if(mod(iter,10)==0) then
    	write(*,"(A24,I7,6E15.7)") "iteration,simtime,dt,res:",iter,simtime,dt,res
        write(filehstid,"(I15,6E15.7)") iter,simtime,dt,res
    end if

    if(mod(iter,500)==0) then
        call writeMacro(iter)
        !call writeNeq(iter)
    end if

    iter = iter+1
    simtime = simtime+dt

end do

write(*,*) 'end loop'

end subroutine loopjob

!-------------------------------------------------------------------

subroutine timestep

real(srk) :: tmin !minimum dt allowed
real(srk) :: temp !max 1/dt allowed
real(srk) :: sos !speed of sound
real(srk) :: prim(dim+2)
integer :: i,j

!set initial value
tmin = 0.0
temp = 0.0

!$omp parallel 

!$omp do private(i,j,sos,prim) reduction(max:temp)
do j=q2imin,q2imax
do i=q1imin,q1imax
    prim = get_primitive(ctr(i,j)%w)
    sos = get_sos(prim(4))

    !maximum velocity(plus speed of sound)
    prim(2) = max(umax,abs(prim(2)))+sos
    prim(3) = max(vmax,abs(prim(3)))+sos

    !maximum 1/dt allowed
    temp = max(temp,(ctr(i,j)%length(2)*prim(2)+ctr(i,j)%length(1)*prim(3))/ctr(i,j)%area,abs(phi)/(vspace(1,15)-vspace(1,14)))
end do
end do 
!$omp end do

!$omp end parallel

tmin = 1.d0/temp
dt = cfl*tmin

end subroutine timestep

!-------------------------------------------------------------------

subroutine interpolation

integer :: i,j

!first order NO NEED interpolation of slope
if(interp_method==firstorder) return 

!$omp parallel

!$omp do
do j=q2imin,q2imax
    call interp_boundary(ctr(q1imin,j),ctr(q1imin,j),ctr(q1imin+1,j),q1dirc)
    call interp_boundary(ctr(q1imax-1,j),ctr(q1imax,j),ctr(q1imax,j),q1dirc)
end do
!$omp end do nowait

!$omp do
do j=q2imin,q2imax
do i=q1imin+1,q1imax-1
    call interp_inner(ctr(i-1,j),ctr(i,j),ctr(i+1,j),q1dirc)
end do
end do
!$omp end do nowait

!$omp do
do i=q1imin,q1imax
    call interp_boundary(ctr(i,q2imin),ctr(i,q2imin),ctr(i,q2imin+1),q2dirc)
    call interp_boundary(ctr(i,q2imax-1),ctr(i,q2imax),ctr(i,q2imax),q2dirc)
end do
!$omp end do nowait

!$omp do
do j=q2imin+1,q2imax-1
do i=q1imin,q1imax
    call interp_inner(ctr(i,j-1),ctr(i,j),ctr(i,j+1),q2dirc)
end do
end do
!$omp end do nowait

!$omp end parallel

end subroutine interpolation

!-------------------------------------------------------------------

subroutine evolution

integer:: i,j

!$omp parallel

!$omp do
do j=q2imin,q2imax
    call flux_isothermal(bc_static, ctr(q1imin,j)%w, ctr(q1imin,j)%h, ctr(q1imin,j)%b, ctr(q1imin,j)%sh, ctr(q1imin,j)%sb, ctr(q1imin,j)%length(q1dirc), &
                         a1face(q1imin,j)%flux, a1face(q1imin,j)%flux_h, a1face(q1imin,j)%flux_b, a1face(q1imin,j)%length, &
                         unum, vnum, uspace, vspace, weight, &
                         ink, gamma, muref, omega, &
                         dt, q1dirc, RN, a1face(q1imin,j)%cosa, a1face(q1imin,j)%sina)
    call flux_isothermal(bc_static, ctr(q1imax,j)%w, ctr(q1imax,j)%h, ctr(q1imax,j)%b, ctr(q1imax,j)%sh, ctr(q1imax,j)%sb, ctr(q1imax,j)%length(q1dirc), &
                         a1face(q1imax+1,j)%flux, a1face(q1imax+1,j)%flux_h, a1face(q1imax+1,j)%flux_b, a1face(q1imax+1,j)%length, &
                         unum, vnum, uspace, vspace, weight, &
                         ink, gamma, muref, omega, &
                         dt, q1dirc, RY, a1face(q1imax+1,j)%cosa, a1face(q1imax+1,j)%sina)
end do
!$omp end do nowait

!$omp do
do j=q2imin,q2imax
do i=q1imin+1,q1imax
    call flux_ugks2d(ctr(i-1,j)%w, ctr(i-1,j)%h, ctr(i-1,j)%b, ctr(i-1,j)%sh, ctr(i-1,j)%sb, ctr(i-1,j)%length(q1dirc), &
                     a1face(i,j)%flux, a1face(i,j)%flux_h, a1face(i,j)%flux_b, a1face(i,j)%length, &
                     ctr(i,j)%w, ctr(i,j)%h, ctr(i,j)%b, ctr(i,j)%sh, ctr(i,j)%sb, ctr(i,j)%length(q1dirc), &
                     unum, vnum, uspace, vspace, weight, &
                     ink, gamma, muref, omega, prandtl, &
                     dt, q1dirc, a1face(i,j)%cosa, a1face(i,j)%sina)
end do
end do
!$omp end do nowait

!$omp do
do i=q1imin,q1imax
    call flux_isothermal(bc_static, ctr(i,q2imin)%w, ctr(i,q2imin)%h, ctr(i,q2imin)%b, ctr(i,q2imin)%sh, ctr(i,q2imin)%sb, ctr(i,q2imin)%length(q2dirc), &
                         a2face(i,q2imin)%flux, a2face(i,q2imin)%flux_h, a2face(i,q2imin)%flux_b, a2face(i,q2imin)%length, &
                         unum, vnum, uspace, vspace, weight, &
                         ink, gamma, muref, omega, &
                         dt, q2dirc, RN, a2face(i,q2imin)%cosa, a2face(i,q2imin)%sina)
    call flux_isothermal(bc_move, ctr(i,q2imax)%w, ctr(i,q2imax)%h, ctr(i,q2imax)%b, ctr(i,q2imax)%sh, ctr(i,q2imax)%sb, ctr(i,q2imax)%length(q2dirc), &
                         a2face(i,q2imax+1)%flux, a2face(i,q2imax+1)%flux_h, a2face(i,q2imax+1)%flux_b, a2face(i,q2imax+1)%length, &
                         unum, vnum, uspace, vspace, weight, &
                         ink, gamma, muref, omega, &
                         dt, q2dirc, RY, a2face(i,q2imax+1)%cosa, a2face(i,q2imax+1)%sina)
end do
!$omp end do nowait

!$omp do
do j=q2imin+1,q2imax
do i=q1imin,q1imax
    call flux_ugks2d(ctr(i,j-1)%w, ctr(i,j-1)%h, ctr(i,j-1)%b, ctr(i,j-1)%sh, ctr(i,j-1)%sb, ctr(i,j-1)%length(q2dirc), &
                     a2face(i,j)%flux, a2face(i,j)%flux_h, a2face(i,j)%flux_b, a2face(i,j)%length, &
                     ctr(i,j)%w, ctr(i,j)%h, ctr(i,j)%b, ctr(i,j)%sh, ctr(i,j)%sb, ctr(i,j)%length(q2dirc), &
                     unum, vnum, uspace, vspace, weight, &
                     ink, gamma, muref, omega, prandtl, &
                     dt, q2dirc, a2face(i,j)%cosa, a2face(i,j)%sina)
end do
end do
!$omp end do nowait

!$omp end parallel 

end subroutine evolution

!-------------------------------------------------------------------

subroutine update

integer :: i,j

!t=t^n
real(srk),allocatable,dimension(:,:) :: H_old,B_old !equilibrium distribution at t=t^n
real(srk) :: w_old(dim+2)
real(srk) :: prim_old(dim+2)
real(srk) :: tau_old

!t=t^n+1
real(srk),allocatable,dimension(:,:) :: H,B !equilibrium distribution at t=t^{n+1}
real(srk),allocatable,dimension(:,:) :: H_plus,B_plus !Shakhov part

real(srk) :: prim(dim+2) !primary variables at t^n and t^{n+1}
real(srk) :: qf(dim)
real(srk) :: tau !collision time and t^n and t^{n+1}
real(srk) :: sum_res(dim+2),sum_avg(dim+2)

!external force
real(srk) :: fx,fy
integer :: p,q
real(srk),allocatable,dimension(:,:) :: hsource
real(srk),allocatable,dimension(:,:) :: bsource

allocate(H_old(unum,vnum))
allocate(B_old(unum,vnum))
allocate(H(unum,vnum))
allocate(B(unum,vnum))
allocate(H_plus(unum,vnum))
allocate(B_plus(unum,vnum))
allocate(hsource(unum,vnum))
allocate(bsource(unum,vnum))

!set initial value
res = 0.0
sum_res = 0.0
sum_avg = 0.0

fx = 0.d0
fy = phi

do j=q2imin,q2imax
do i=q1imin,q1imax 
    !--------------------------------------------------
    ! store W^n and calculate H^n,B^n,\tau^n
    !--------------------------------------------------
    w_old = ctr(i,j)%w
    prim_old = get_primitive(w_old)

    call reduced_maxwell(H_old,B_old,uspace,vspace,prim_old)

    tau_old=get_tau(prim_old)

    !--------------------------------------------------
    ! update W^{n+1} and calculate H^{n+1},B^{n+1},\tau^{n+1}
    !--------------------------------------------------
    ctr(i,j)%w = ctr(i,j)%w+(a1face(i,j)%flux-a1face(i+1,j)%flux+a2face(i,j)%flux-a2face(i,j+1)%flux)/ctr(i,j)%area
    prim = get_primitive(ctr(i,j)%w)

    call reduced_maxwell(H,B,uspace,vspace,prim)
    tau = get_tau(prim)

    !--------------------------------------------------
    ! record residual
    !--------------------------------------------------
    sum_res = sum_res+(w_old-ctr(i,j)%w)**2
    sum_avg = sum_avg+abs(ctr(i,j)%w)

    !--------------------------------------------------
    ! Shakhov part
    !--------------------------------------------------
    !heat flux at t=t^n
    qf = get_heat_flux(ctr(i,j)%h,ctr(i,j)%b,uspace,vspace,prim_old) 

    !h^+ = H+H^+ at t=t^n
    call shakhov_part(H_old,B_old,uspace,vspace,qf,prim_old,H_plus,B_plus)
    H_old = H_old+H_plus !h^+
    B_old = B_old+B_plus !b^+

    !h^+ = H+H^+ at t=t^{n+1}
    call shakhov_part(H,B,uspace,vspace,qf,prim,H_plus,B_plus)
    H = H+H_plus
    B = B+B_plus

    !--------------------------------------------------
    ! update distribution function
    !--------------------------------------------------
    ctr(i,j)%h = (ctr(i,j)%h+(a1face(i,j)%flux_h-a1face(i+1,j)%flux_h+&
                 a2face(i,j)%flux_h-a2face(i,j+1)%flux_h)/ctr(i,j)%area+&
                 0.5*dt*(H/tau+(H_old-ctr(i,j)%h)/tau_old))/(1.0+0.5*dt/tau)
    ctr(i,j)%b = (ctr(i,j)%b+(a1face(i,j)%flux_b-a1face(i+1,j)%flux_b+&
                 a2face(i,j)%flux_b-a2face(i,j+1)%flux_b)/ctr(i,j)%area+&
                 0.5*dt*(B/tau+(B_old-ctr(i,j)%b)/tau_old))/(1.0+0.5*dt/tau)

end do
end do

!final residual
res = sqrt(grid%numElem*sum_res)/(sum_avg+SMV)

end subroutine update

!-------------------------------------------------------------------

end module loopmod