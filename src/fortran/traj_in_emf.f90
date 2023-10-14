module traj_in_emf_m
 use ode_solver_m, only : ode_solution, ode_solver
 use vector_ops_mod, only : vector_ops
 use phys_const, only : e_charge, erg_eV, c_light, me_ev, pi, mprot
 use util_test, only : cwd_parent_dir

 implicit none
 private
 save

 public :: main_calc


 real(8) :: mc2, emc2, eloss, emc1
 real(8) :: bm0_ch, rc0_ch, chi0_ch(3)
 real(8) :: bm0, rho, chi(3), rot_num, accuracy
 real(8) :: start_pos(3), start_vel(3), start_lorf
 real(8) :: u1_const, u2_const, u3_const
 character(2) :: field_type
 character(10) :: particle
 character(500) :: fname(4)
 logical :: loss_on
 
 
 contains


subroutine initiate_constants
!=================================================
! Constants for ODE system
!=================================================

  ! c
  u1_const = c_light
  ! e/mc
  u2_const = e_charge*c_light/(me_ev*erg_eV)
  ! (2/3)(e^2/mc)(1/c^2)
  u3_const = u2_const*e_charge*2d0/(3d0*c_light**2)

end subroutine initiate_constants

function efield(xc)
  real(8), intent(in) :: xc(3)
  real(8):: efield(3)
 
  efield = 1d0
 
end function efield

subroutine em_field(xc, e_field, b_field)
  real(8), intent(in) :: xc(3)
  real(8), intent(out) :: e_field(3), b_field(3)
  type(vector_ops) :: vop

  e_field = [0d0, 0d0, 0d0]
  b_field = [0d0, 0d0, 1d3]

end subroutine em_field  


subroutine motion_ode(x, y, dydx)
  real(8), intent(in) :: x, y(:)
  real(8), intent(out) :: dydx(:)
  real(8) :: xc(3), vc(3), gamma, em(3), bm(3), vde, vxb(3), dbdt2 
  type(vector_ops) :: vop
  
  xc = y(1:3)
  vc = y(4:6)
  gamma = y(7)

  ! 1st eq: dr/dt = c*beta
  dydx(1:3)=u1_const*vc

  call em_field(xc, em, bm)
  vde=dot_product(vc, em)
  vxb = vop%cross_prod(vc, bm)
  ! 2nd eq: dbeta/dt = e/(mc*gamma)*(E - (beta*E)*beta + [beta x B])
  dydx(4:6)=(u2_const/gamma)*(em - vde*vc + vxb)


  dbdt2 = dot_product(dydx(4:6), dydx(4:6))
  ! 3rd eq: dgamma/dt = (e/mc)*(beta*E) 
  ! - 2/3(e^2/mc)(gamma^4/c2)*(dbeta/dt)^2
  dydx(7) = u2_const*vde - u3_const*dbdt2*gamma**4

 end subroutine motion_ode

subroutine cond1(this)
!================================================
! Condition for fast calculation
!================================================
  type(ode_solver) :: this
  type(vector_ops) :: vop
  ! real(8) :: rad, vel
  ! real(8) :: dydx(6), min_vel=1d-2, min_dist=1d0*psec
  real(8) :: min_vel=1d-2, maxvel=1d0 + 1d-5, vnorm
  integer :: i
  integer, save :: nstep = 0
  real(8), save :: gamma_init

  if (nstep == 0) then
    gamma_init = this%yc(7)
  end if
  
  ! Stop condition:
  ! If initial gamma decreased to 1%
  if (this%yc(7)/gamma_init < 1e-2) then
    this%stop_flag = .true.
  end if  
  nstep=nstep+1

  ! normalize velocity
  vnorm = vop%vec_norm(this%yc(4:6))
  if (vnorm > maxvel) then
    this%yc(4:6) = this%yc(4:6)/vnorm
  end if  
  
  ! this%yeps_flag=.false.
  if (mod(nstep,100000)==0) then
  write(*,*) nstep, this%hc, this%yc(7)
  end if

end subroutine cond1


 subroutine calc_syst
  type(ode_solution) :: sol
  type(ode_solver) :: ode
  integer, parameter :: ns=7
  real(8) :: y(ns)
  real(8) :: t1, t2, eps
 
  call initiate_constants 
  y(1:3)=[0d0, 0d0, 0d0]
  y(4:6)=[1d0, 0d0, 0d0]
  y(7) = 1d6

  t1 = 0d0
  t2 = 1d1
  eps = 1d-6
 
  call ode%ode_sys(ns, motion_ode)
  call ode%init_cond(t1,y)
  call ode%solve_to(t2)
  call ode%sol_out_to(sol)
  call ode%toler(eps)
  call ode%init_step(1d-8)
  call ode%add_cond(1,cond1)
  call ode%solve_ode
 
  write(*,*) "Solved with np = ", sol%np
  call write_results(sol)
 
 
 end subroutine calc_syst
 
 
 subroutine write_results(sol)
  integer :: i
  type(ode_solution) :: sol
  character(500) :: fname(4)

  fname(1)=trim(cwd_parent_dir(2))//'/data/test_res2.dat'
  write(*,*) "fname = ", fname(1)
  open(1,file=fname(1))
  
  do i=1,sol%np
  !  write(*,'(A,Es14.6,A,Es14.6)') "time = ", sol%y(0, i), " gamma = ", sol%y(7, i)
   write(1,*) sol%y(0:7, i)
  end do
  
  close(1)
  write(*, *) "Total points = ", sol%np
 
 
 end subroutine write_results


subroutine main_calc

 call calc_syst

end subroutine main_calc


end module traj_in_emf_m



program main
 use traj_in_emf_m, only : main_calc
 
 call main_calc

end program main




