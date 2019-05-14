
! loop.f90
!-------------------------------------------------------------------

module loopmod

use datamod

use kineticmod

use mathmod

use UGKSFlux1D

use postmod

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

    if(simtime >= maxtime) exit
    !if(all(res<eps)) exit

    if(mod(iter,10)==0) then
    	write(*,"(A24,I8,6E15.7)") "iteration,simtime,dt,res:",iter,simtime,dt,res
        write(filehstid,"(I15,6E15.7)") iter,simtime,dt,res
    end if

    if(mod(iter,1000)==0) then
        call writeMacro(iter)
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
!$omp do private(i,sos,prim) reduction(max:temp)
do i=q1imin,q1imax
    prim = get_primitive(ctr(i)%w)
    sos = get_sos(prim(3))

    !maximum velocity(plus speed of sound)
    prim(2) = max(umax,abs(prim(2)))+sos

    !maximum 1/dt allowed
    temp = max(temp,prim(2)/ctr(i)%length)
end do
!$omp end do
!$omp end parallel

tmin = 1.d0/temp
dt = cfl*tmin

end subroutine timestep

!-------------------------------------------------------------------

subroutine interpolation

integer :: i,j

call interp_boundary(ctr(q1imin),ctr(q1imin),ctr(q1imin+1))
call interp_boundary(ctr(q1imax-1),ctr(q1imax),ctr(q1imax))

!$omp parallel

!$omp do
do i=q1imin+1,q1imax-1
    call interp_macro(ctr(i-1),ctr(i),ctr(i+1))
    call interp_micro(ctr(i-1),ctr(i),ctr(i+1))
end do
!$omp end do nowait

!$omp end parallel

end subroutine interpolation

!-------------------------------------------------------------------

subroutine evolution

integer :: i

!$omp parallel

!$omp do
do i=q1imin,q1imax+1
    !call flux_ugks1d(ctr(i-1),face(i),ctr(i))
    call flux_ugks1d(ctr(i-1)%w, ctr(i-1)%h, ctr(i-1)%b, ctr(i-1)%sh, ctr(i-1)%sb, ctr(i-1)%length, face(i)%flux, face(i)%flux_h, face(i)%flux_b, ctr(i)%w, ctr(i)%h, ctr(i)%b, ctr(i)%sh, ctr(i)%sb, ctr(i)%length, unum, uspace, weight, ink, gamma, muref, omega, prandtl, dt)
end do
!$omp end do nowait

!$omp end parallel

end subroutine evolution

!-------------------------------------------------------------------

subroutine update

integer :: i

!t=t^n
real(srk),allocatable,dimension(:) :: H_old,B_old !equilibrium distribution at t=t^n
real(srk) :: w_old(dim+2)
real(srk) :: prim_old(dim+2)
real(srk) :: tau_old

!t=t^n+1
real(srk),allocatable,dimension(:) :: H,B !equilibrium distribution at t=t^{n+1}
real(srk),allocatable,dimension(:) :: H_plus,B_plus !Shakhov part

real(srk) :: prim(dim+2) !primary variables at t^n and t^{n+1}
real(srk) :: qf
real(srk) :: tau !collision time and t^n and t^{n+1}
real(srk) :: sum_res(dim+2),sum_avg(dim+2)

allocate(H_old(unum))
allocate(B_old(unum))
allocate(H(unum))
allocate(B(unum))
allocate(H_plus(unum))
allocate(B_plus(unum))

!set initial value
res = 0.0
sum_res = 0.0
sum_avg = 0.0

do i=q1imin,q1imax
    !--------------------------------------------------
    ! store W^n and calculate H^n,B^n,\tau^n
    !--------------------------------------------------
    w_old = ctr(i)%w

    prim_old = get_primitive(w_old)
    call reduced_maxwell(H_old,B_old,prim_old)
    tau_old = get_tau(prim_old)

    !--------------------------------------------------
    ! update W^{n+1} and calculate H^{n+1},B^{n+1},\tau^{n+1}
    !--------------------------------------------------
    ctr(i)%w = ctr(i)%w+(face(i)%flux-face(i+1)%flux)/ctr(i)%length
    ctr(i)%prim = get_primitive(ctr(i)%w)

    prim = get_primitive(ctr(i)%w)
    call reduced_maxwell(H,B,prim)
    tau = get_tau(prim)

    !--------------------------------------------------
    ! record residual
    !--------------------------------------------------
    sum_res = sum_res+(w_old-ctr(i)%w)**2
    sum_avg = sum_avg+abs(ctr(i)%w)

    !--------------------------------------------------
    ! Shakhov part
    !--------------------------------------------------
    !heat flux at t=t^n
    qf = get_heat_flux(ctr(i)%h,ctr(i)%b,prim_old) 

    !h^+ = H+H^+ at t=t^n
    call shakhov_part(H_old,B_old,qf,prim_old,H_plus,B_plus) !H^+ and B^+
    H_old = H_old+H_plus !h^+
    B_old = B_old+B_plus !b^+

    !h^+ = H+H^+ at t=t^{n+1}
    call shakhov_part(H,B,qf,prim,H_plus,B_plus)
    H = H+H_plus
    B = B+B_plus

    !--------------------------------------------------
    ! update distribution function
    !--------------------------------------------------
    ctr(i)%h = (ctr(i)%h+(face(i)%flux_h-face(i+1)%flux_h)/ctr(i)%length+&
                0.5*dt*(H/tau+(H_old-ctr(i)%h)/tau_old))/(1.0+0.5*dt/tau)
    ctr(i)%b = (ctr(i)%b+(face(i)%flux_b-face(i+1)%flux_b)/ctr(i)%length+&
                0.5*dt*(B/tau+(B_old-ctr(i)%b)/tau_old))/(1.0+0.5*dt/tau)

end do

!final residual
res = sqrt(q1imax*sum_res)/(sum_avg+SMV)

end subroutine update

!-------------------------------------------------------------------

end module loopmod
