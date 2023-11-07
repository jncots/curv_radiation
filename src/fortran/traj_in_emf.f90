module traj_in_emf_m
 use ode_solver_m, only : ode_solution, ode_solver
 use vector_ops_mod, only : vector_ops
 use phys_const, only : e_charge, erg_eV, c_light, me_ev, pi, mprot
 use phys_const, only : hbar
 use util_test, only : cwd_parent_dir
 use mag_sph_m, only : mfields, norm_field, solid_rot, mf_rlc


 implicit none
 private
 save

 public :: main_calc


 real(8) :: rlight_cylinder
 real(8) :: mc2, emc2, eloss, emc1
 real(8) :: bm0_ch, rc0_ch, chi0_ch(3)
 real(8) :: bm0, rho, chi(3), rot_num, accuracy
 real(8) :: start_pos(3), start_vel(3), start_lorf
 real(8) :: u1_const, u2_const, u3_const
 real(8) :: dnde_const, ec_const, loss_const, rgc_const, mc2e_const
 character(2) :: field_type
 character(10) :: particle
 character(500) :: fname(4)
 logical :: loss_on
 
 
 contains


subroutine initiate_constants
  real(8) :: alpha
!=================================================
! Constants for ODE system
!=================================================

  ! c
  u1_const = c_light
  ! e/mc
  u2_const = e_charge*c_light/(me_ev*erg_eV)
  ! (2/3)(e^2/mc)(1/c^2)
  u3_const = u2_const*e_charge*2d0/(3d0*c_light**2)

  ! For radiation:
  ! (1/sqrt(3)*pi)*(alpha/hbar)
  dnde_const = (e_charge/hbar)**2/(sqrt(3d0)*pi*c_light)
  dnde_const = dnde_const*erg_eV
  ! sqrt(3/2*c*mc^2/e^2)

  ec_const = sqrt(3*me_ev*erg_eV*c_light/(2*e_charge**2))
  loss_const = 1/ec_const
  ! Characteristic energy in eV:
  ec_const = ec_const*3*hbar/(2*erg_eV)

  rgc_const = (me_ev*erg_eV)/(e_charge*c_light)
  mc2e_const = (me_ev*erg_eV)/(e_charge)

end subroutine initiate_constants

subroutine initiate_fields

  call norm_field
  call solid_rot(90d0)
  ! rlight_cylinder = 1.5e7

  ! write(*,'(A,Es14.6)') "rlight_cylinder = ", rlight_cylinder

end subroutine initiate_fields  




subroutine em_field(xc, e_field, b_field)
  real(8), intent(in) :: xc(3)
  real(8), intent(out) :: e_field(3), b_field(3)
  real(8) :: xc1(3)
  type(vector_ops) :: vop


  xc1 = xc/mf_rlc

  ! write(*,'(A, 3Es20.10)') 'xc1 = ', xc1
  call mfields(xc1,b_field,e_field, 'fff')
  ! write(*,'(A, 3Es20.10)') 'Bfield = ', b_field

  e_field = [0d0, 0d0, 0d0]
  ! b_field = [0d0, 0d0, 1d3]
  ! read(*,*)

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
  ! dydx(7) = 0d0

 end subroutine motion_ode


subroutine rad_val(x, y, acc, loss, ec, loss_pr, acc_norm)
  real(8), intent(in) :: x, y(:)
  real(8), intent(out) :: acc, loss, ec, loss_pr, acc_norm
  real(8) :: dydx(7), bmmod
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

  acc =  u2_const*vde
  loss = u3_const*dbdt2*gamma**4
  ! dydx(7) = u2_const*vde - u3_const*dbdt2*gamma**4
  ! dydx(7) = 0d0

  ! Characteristic energy in eV
  ec = ec_const*y(7)*sqrt(loss)

  ! gamma change per rotation
  loss_pr = loss_const*sqrt(loss)*y(7)**2


  bmmod = sqrt(dot_product(bm, bm))
  acc_norm = rgc_const*y(7)*sqrt(dbdt2)/bmmod

end subroutine rad_val

subroutine cond1(this)
!================================================
! Condition for fast calculation
!================================================
  type(ode_solver) :: this
  type(vector_ops) :: vop
  ! real(8) :: rad, vel
  ! real(8) :: dydx(6), min_vel=1d-2, min_dist=1d0*psec
  real(8) :: min_vel=1d-2, maxvel=1d0 + 1d-5, vnorm 
  real(8) :: rg, em(3), bm(3), bm_val
  integer :: i
  integer, save :: nstep = 0
  real(8), save :: gamma_init

  if (nstep == 0) then
    gamma_init = this%yc(7)
  end if
  
  ! Stop condition:
  ! If initial gamma decreased to 0.1%
  if (this%yc(7)/gamma_init < 1e-4) then
    this%stop_flag = .true.
  end if
  
  call em_field(this%yc(1:3), em, bm)
  bm_val = vop%vec_norm(bm)
  if (bm_val < 1e0) then
    this%stop_flag = .true.
    write(*,*) "Magnetic field small = ", bm_val
  end if  
  
  nstep=nstep+1

  ! normalize velocity
  vnorm = vop%vec_norm(this%yc(4:6))
  if (vnorm > maxvel) then
    this%yc(4:6) = this%yc(4:6)/vnorm
  end if  
  
  ! this%yeps_flag=.false.
  if (mod(nstep,10000)==0) then
    vnorm = vop%vec_norm(this%yc(1:3))
    call em_field(this%yc(1:3), em, bm)
    bm_val = vop%vec_norm(bm)
    rg = mc2e_const*this%yc(7)/bm_val
    write(*,'(I10,10Es14.6)') nstep, this%xc, vnorm, this%yc(7), bm_val, &
                              rg, rg/c_light

    ! read(*,*)                          
  end if

end subroutine cond1


 subroutine calc_syst
  type(ode_solution) :: sol, sol_red
  type(ode_solver) :: ode
  type(vector_ops) :: vop
  integer, parameter :: ns=7
  real(8) :: y(ns)
  real(8) :: t1, t2, eps
  real(8) :: bm(3), em(3)
 
  call initiate_constants
  call initiate_fields
  
  write(*,'(A,Es14.6)') "mf_rlc =", mf_rlc
  ! read(*,*)


  y(1:3)=[4d-1, 4d-1, 4d-1]*mf_rlc

  call em_field(y(1:3), em, bm)
  y(4:6)=vop%vec_dir([1d0, 1d0, 1d0])
  y(7) = 1d7

  t1 = 0d0
  t2 = 1d0
  eps = 1d-6
 
  call ode%ode_sys(ns, motion_ode)
  call ode%init_cond(t1,y)
  call ode%solve_to(t2)
  call ode%sol_out_to(sol)
  call ode%toler(eps)
  call ode%init_step(1d-20)
  call ode%add_cond(1,cond1)
  call ode%solve_ode
 
  call sol%reduce_to_npoints(10000,sol_red) 
  write(*,*) "Solved with np = ", sol%np
  call write_results(sol_red)
 
 
 end subroutine calc_syst
 
 
 subroutine write_results(sol)
  integer :: i
  type(ode_solution) :: sol
  real(8) :: bm(3), em(3)
  real(8) :: acc, loss, ec, loss_pr, acc_norm
  character(500) :: fname(4)

  fname(1)=trim(cwd_parent_dir(3))//'/results/data/traj_in_fff.dat'
  write(*,*) "fname = ", fname(1)
  open(1, file=fname(1))
  
  do i=1,sol%np
  !  write(*,'(A,Es14.6,A,Es14.6)') "time = ", sol%y(0, i), " gamma = ", sol%y(7, i)
    call em_field(sol%y(1:3, i), em, bm)
    call rad_val(sol%y(0, i), sol%y(1:, i), acc, &
    loss, ec, loss_pr, acc_norm)
    write(1,'(20Es40.30)') sol%y(0:7, i), ec, acc, &
    loss, loss_pr, loss_pr/sol%y(7, i), sol%y(7, i)/loss, acc_norm 
  end do
  
  close(1)
  write(*, *) "Total points = ", sol%np

  call write_int_spectr(sol)
 
 
 end subroutine write_results


 subroutine write_int_spectr(sol)
  use utools_mod, only : utools
  type(utools) :: ut
  type(ode_solution) :: sol
  real(8), allocatable :: egam(:), spectr(:)
  real(8), allocatable :: ec(:), gm(:), tm(:)
  real(8) :: acc1, loss1, ec1, loss_pr1, acc_norm
  integer :: i, ngam, nc
  character(500) :: fname(4)

  ngam = 100
  call ut%grid(egam, 1d6, 1d18, ngam,'log')
  call ut%grid(spectr, 1d8, 1d13, ngam,'log')

  nc = sol%np
  allocate(ec(1:nc))
  allocate(gm(1:nc))
  allocate(tm(1:nc))

  do i=1,sol%np
      call rad_val(sol%y(0, i), sol%y(1:, i), &
      acc1, loss1, ec1, loss_pr1, acc_norm)
      ec(i) = ec1
      gm(i) = sol%y(7, i)
      tm(i) = sol%y(0, i)
  end do 

  call int_spectr(ec, gm, tm, nc, egam, ngam, spectr)

  fname(1)=trim(cwd_parent_dir(3))//'/results/data/spectr.dat'
  write(*,*) "fname = ", fname(1)
  open(1,file=fname(1))

  do i = 1, ngam
    write(1,*) egam(i), dnde_const*spectr(i)
  end do
  
  close(1)

 
 end subroutine write_int_spectr  


 subroutine show_field
 ! Function for plotting field lines  
  use mag_sph_m, only : mfields, norm_field, solid_rot
  use utools_mod, only : utools
  use util_test, only : cwd_parent_dir
  type(utools) :: ut
  real(8) :: bmf(3), emf(3), xc(3), cube_size
  real(8), allocatable ::xgrid(:), ygrid(:)
  integer :: ix, iy, nx, ny
  character(500) :: fname

  call norm_field
 call solid_rot(0d0)

 nx = 201
 ny = 201
 cube_size = 2d0
  call ut%grid(xgrid,-cube_size,cube_size,nx,'lin')
  call ut%grid(ygrid,-cube_size,cube_size,ny,'lin')

  fname = trim(cwd_parent_dir(3))//'/results/data/mag_field/bm_201.dat'
  open(1, file=fname)

  do ix=1,nx
    do iy=1,ny
      xc = [xgrid(ix), ygrid(iy), 0d0]
      call mfields(xc,bmf,emf,'fff')
      write(1, *) xc, bmf
    end do
  end do    
  close(1)
!  write(*, '(A, 3Es14.6)') "bmf = ", bmf
!  write(*, '(A, 3Es14.6)') "emf = ", emf

 end subroutine show_field

 subroutine show_field3d
  ! Function for plotting field lines  
   use mag_sph_m, only : mfields, norm_field, solid_rot
   use utools_mod, only : utools
   use util_test, only : cwd_parent_dir
   type(utools) :: ut
   real(8) :: bmf(3), emf(3), xc(3), cube_size
   real(8), allocatable ::xgrid(:), ygrid(:), zgrid(:)
   integer :: ix, iy, iz,  nx, ny, nz
   character(500) :: fname
 
  call norm_field
  call solid_rot(0d0)
 
  nx = 201
  ny = 201
  nz = 201
  cube_size = 3d0
   call ut%grid(xgrid,-cube_size,cube_size,nx,'lin')
   call ut%grid(ygrid,-cube_size,cube_size,ny,'lin')
   call ut%grid(zgrid,-cube_size,cube_size,nz,'lin')
 
   fname = trim(cwd_parent_dir(3))//'/results/data/mag_field/bm_3d_dip.dat'
   open(1, file=fname)
 
   do ix=1,nx
     do iy=1,ny
      do iz=1,nz
       xc = [xgrid(ix), ygrid(iy), zgrid(iz)]
       call mfields(xc,bmf,emf)
       write(1, *) xc, bmf
      end do 
     end do
   end do    
   close(1)
 !  write(*, '(A, 3Es14.6)') "bmf = ", bmf
 !  write(*, '(A, 3Es14.6)') "emf = ", emf
 
  end subroutine show_field3d


 subroutine int_spectr(ec, gm, tm, nc, egam, ngam, spectr)
  ! e_critical, gamma, time
  real(8), intent(in) :: ec(:), gm(:), tm(:), egam(:)
  integer, intent(in) :: nc, ngam
  real(8), intent(out) :: spectr(:)
  integer iegam, iec
  real(8) :: res, fsyn, xc, res_step

  do iegam = 1, ngam
    ! spectr(iegam) = 0d0

    write(*,'(A,Es14.6)') "Calc spec for egam = ", egam(iegam)
    res = 0d0
    do iec = 1, nc
      xc = egam(iegam)/ec(iec)
      call f_synch(xc, fsyn)
      fsyn = fsyn/(xc*gm(iec)**2)
      if (iec.ne.nc) then
        res_step =  fsyn*(tm(iec+1) - tm(iec))
      else
        res_step = fsyn*(tm(iec) - tm(iec-1)) 
      end if
      if (res_step == res_step) res = res + res_step
    end do  
    spectr(iegam) = res
  end do
 
 end subroutine int_spectr

 subroutine f_synch(x,res)
  !
  ! Emmisivity function in homogeneous magnetic field
  !
   real(8), intent(in) :: x
   real(8), intent(out) :: res
   real(8) :: x13, x23, x43, p1, p2, p3
   
   x13=x**(1d0/3d0)
   x23=x13**2
   x43=x23**2
   
   p1=1d0+0.884d0*x23+0.471d0*x43
   p2=1d0+1.64d0*x23+0.974d0*x43
   p3=2.15d0*x13*(1d0+3.06*x)**(1d0/6d0)
   
   res=p3*(p1/p2)*exp(-x)
  end subroutine f_synch


subroutine main_calc

!  call calc_syst
!  call show_field
  call show_field3d

end subroutine main_calc


end module traj_in_emf_m



program main
 use traj_in_emf_m, only : main_calc
 use mag_sph_m, only : mfields, norm_field, solid_rot
 use utools_mod, only : utools
 use util_test, only : cwd_parent_dir
 type(utools) :: ut
 real(8) :: bmf(3), emf(3), xc(3)
 real(8), allocatable ::xgrid(:), ygrid(:)
 integer :: ix, iy, nx, ny
 character(500) :: fname
 
 
 call main_calc

end program main




