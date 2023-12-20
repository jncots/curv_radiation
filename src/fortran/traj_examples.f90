!*************************************************************************
!      Calculation trajectories in Y-point field
!*************************************************************************
! Author     : Anton Prosekin
! Date       : 2023-12-20
! Description:
!    I need to output trajectories with directions of velocity
!    * plot trajectories
!    * plot direction of radiation
!    * plot max point of ec
!    * plot ec evolution
! Usage:
!
! Notes:
!
!
! Log:
!   Date       | Modification Description
!   -----------|-------------------------------------
!   2023-12-20 | Initial version created
!
!*************************************************************************

module traj_in_emf_m
   use ode_solver_m, only : ode_solution, ode_solver
   use vector_ops_mod, only : vector_ops
   use utools_mod, only : utools
   use phys_const, only : e_charge, erg_eV, hbar, c_light, me_ev
   use phys_const, only : pi
   use util_test, only : cwd_parent_dir
   use hyper_field_m, only : hyper_field


   implicit none
   private
   save

   public :: main_calc

   real(8) :: u1_const, u2_const, u3_const
   real(8) :: dnde_const, ec_const, loss_const, rgc_const, mc2e_const
   real(8) :: bfield_val
   real(8) :: ec_min0, ec_max0
   character(500) :: traj_file_name, spectrum_file_name
   type(hyper_field) :: hfield
   type(vector_ops) :: vop
   type(utools) :: ut
   integer :: nstep_cond = 0


contains


   function root_out_dir()
      character(500) :: root_out_dir
      root_out_dir = trim(cwd_parent_dir(3))//'/results/data/05_traj_examples_ypoint/'
   end function root_out_dir

   function r_gyro(gamma, bm_field, sinp)
!================================================
!   returns gyroradius
!   gamma - Lorentz factor
!   bm_field - module of magnetic field
!   sinp - pitch angle
!================================================
      real(8), intent(in) :: gamma, bm_field, sinp
      real(8) :: r_gyro
      r_gyro = mc2e_const*gamma/(bm_field*sinp)
   end function r_gyro

   subroutine initiate_constants
!=================================================
! Constants for ODE system
!=================================================
      ! c
      u1_const = c_light
      ! e/mc
      u2_const = -e_charge*c_light/(me_ev*erg_eV)
      ! (2/3)(e^2/mc)(1/c^2)
      u3_const = -u2_const*e_charge*2d0/(3d0*c_light**2)

      ! write(*,*) "u1_const = ", u1_const
      ! write(*,*) "u2_const = ", u2_const
      ! write(*,*) "u3_const = ", u3_const
      ! read(*,*)

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
      ! mc**2/e
      mc2e_const = (me_ev*erg_eV)/(e_charge)

   end subroutine initiate_constants


   subroutine em_field(xc, e_field, b_field)
      real(8), intent(in) :: xc(3)
      real(8), intent(out) :: e_field(3), b_field(3)
      real(8) :: vx, vy, tpar, bpar

      ! write(*,'(A, 3Es20.10)') 'xc1 = ', xc1
      call hfield%field_ypoint(xc(1), xc(2), vx, vy, tpar, bpar)
      ! write(*,'(A, 10Es20.10)') 'Bfield = ', xc(1), xc(2), vx, vy, tpar, bpar

      ! call hfield%field(-1.8d0, 1d0, vx, vy, tpar, bpar)
      ! write(*,'(A, 10Es20.10)') 'Bfield = ', xc(1), xc(2), vx, vy, tpar, bpar

      ! read(*,*)

      b_field = bfield_val*[vx, vy, 0d0]
      e_field = [0d0, 0d0, 0d0]

   end subroutine em_field




   !  subroutine em_field(xc, e_field, b_field)
   !     real(8), intent(in) :: xc(3)
   !     real(8), intent(out) :: e_field(3), b_field(3)
   !     real(8) :: xc1(3)


   !     xc1 = xc/mf_rlc

   !     ! write(*,'(A, 3Es20.10)') 'xc1 = ', xc1
   !     call mfields(xc1,b_field,e_field, 'fff')
   !     ! write(*,'(A, 3Es20.10)') 'Bfield = ', b_field

   !     e_field = [0d0, 0d0, 0d0]
   !     ! b_field = [0d0, 0d0, 1d3]
   !     ! read(*,*)

   !  end subroutine em_field


   subroutine motion_ode(x, y, dydx)
      real(8), intent(in) :: x, y(:)
      real(8), intent(out) :: dydx(:)
      real(8) :: xc(3), vc(3), gamma, em(3), bm(3), vde, vxb(3), dbdt2

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


   subroutine rad_val(x, y, acc, loss, ec, loss_pr, acc_norm)
      real(8), intent(in) :: x, y(:)
      real(8), intent(out) :: acc, loss, ec, loss_pr, acc_norm
      real(8) :: dydx(7), bmmod
      real(8) :: xc(3), vc(3), gamma, em(3), bm(3), vde, vxb(3), dbdt2

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
      real(8) :: maxvel=1d0 + 1d-5, vnorm, rnorm
      real(8) :: rg, em(3), bm(3), bm_val
      real(8) :: acc, loss, ec, loss_pr, acc_norm
      integer, save :: nstep = 0
      real(8), save :: gamma_init
      real(8), save :: ec_max = 0, ec_min = 0

      if (nstep_cond == 0) then
         nstep = 0
         ec_max = 0
         ec_min = 0
         gamma_init = this%yc(7)
         nstep_cond = 1
      end if

      ! Stop condition:
      ! If initial gamma decreased to 0.1%
      if (this%yc(7)/gamma_init < 1e-5) then
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

      call rad_val(this%yc(0), this%yc(1:), acc, &
         loss, ec, loss_pr, acc_norm)

      if (ec > ec_max) then
         ec_max = ec
         ec_min = ec_max
      else
         if (ec < ec_min) then
            ec_min = ec
         end if
      end if


      if (ec < ec_max*1d-5) then
         this%stop_flag = .true.
         ec_min0 = ec_min
         ec_max0 = ec_max
         write(*,'(A,10Es14.6)') "ec = ", ec
         write(*,'(A,10Es14.6)') "ec_max = ", ec_max
         write(*,'(A,10Es14.6)') "ec_min = ", ec_min
      end if

      ! this%yeps_flag=.false.
      if (mod(nstep,1000000)==0) then
         rnorm = vop%vec_norm(this%yc(1:3))
         call em_field(this%yc(1:3), em, bm)
         bm_val = vop%vec_norm(bm)
         rg = mc2e_const*this%yc(7)/bm_val
         write(*,'(I10,10Es14.6)') nstep, this%xc, rnorm, this%yc(7), ec, &
            rg, rg/c_light

         ! read(*,*)
      end if

   end subroutine cond1


   subroutine calculate_trajectory
      integer, parameter :: ns=7
      real(8) :: y0(ns)
      integer :: nquad, branch
      real(8) :: gamma, cross_angle
      real(8) :: rmin, xstart, r1, r2
      integer :: i, nrmin, lunit
      real(8), allocatable :: rmin_grid(:)
      character(500) :: fname

      call initiate_constants

      ! Initial Lorentz factor
      gamma = 1d7

      ! Set field parameters
      !quadrant
      nquad = 2
      ! Crossing angle in degrees
      cross_angle = 10d0
      ! Magnetic field strength
      bfield_val = 1d6

      ! Set starting field line
      branch = 1

      ! nrmin = 100
      ! r1 = 1d0
      ! r2 = 1d8
      ! call ut%grid(rmin_grid, r1, r2, nrmin, "log")
      ! fname=trim(cwd_parent_dir(3))//'/results/data/04_ec_plots_ypoint_field/ec_max_g1e8_a10.dat'
      ! lunit = 77
      ! open(lunit, file = fname)


      ! do i=1,nrmin
      !    ! Semi-major axis
      !    rmin = rmin_grid(i)
      !    ! Starting position
      !    xstart = -rmin*1d2

      !    y0 = initial_cond(gamma, cross_angle, nquad, branch, rmin, xstart)
      !    call calc_syst(y0)
      !    ! write(lunit,"(10(A,Es14.6))") "rmin=", rmin, " ecmin=", ec_min0, " ecmax=", ec_max0
      !    ! write(lunit,*) rmin, ec_min0, ec_max0
      ! end do
      ! close(lunit)


      rmin = 1d7
      xstart = -rmin*1d2
      traj_file_name = "01_run/tgam_1e7_a10.dat"
      spectrum_file_name = "01_run/spectrum.dat"
      y0 = initial_cond(gamma, cross_angle, nquad, branch, rmin, xstart)
      call calc_syst(y0)


   end subroutine calculate_trajectory


   function initial_cond(gamma, cross_angle, nquad, branch, rmin, xstart)
      real(8), intent(in) :: gamma, cross_angle, rmin, xstart
      integer, intent(in) :: nquad, branch
      integer, parameter :: ns=7
      real(8) :: initial_cond(ns)
      real(8) :: y(ns), ystart
      real(8) :: bm(3), em(3)

      call hfield%set_angle(cross_angle)
      call hfield%ycoord(nquad, branch, rmin, xstart, ystart)
      y(1:3)=[xstart, ystart, 0d0]

      call em_field(y(1:3), em, bm)
      bm(2) = bm(2)*(1-1e-6)
      y(4:6) = vop%vec_dir(bm)

      ! Change direction of velocity, so that
      ! vx > 0, i.e. directed towards
      ! center of coordinates
      if (y(4) < 0) then
         y(4:6) = vop%vec_dir(-bm)
      end if
      y(7) = gamma
      initial_cond = y

      ! Print setup info
      call setup_info(nquad, rmin, xstart, ystart, gamma, bm, y)

   end function initial_cond


   subroutine setup_info(nquad, rmin, xstart, ystart, gamma, bm, y)
      integer, intent(in) :: nquad
      real(8), intent(in) :: rmin, xstart, ystart, gamma, bm(3), y(7)
      real(8) :: xstart_test, tpar, rg0
      real(8) :: nbx, nby, tpar_test, rmin_test

      call hfield%xcoord(nquad, ystart, rmin, xstart_test, tpar)
      rg0 = r_gyro(gamma, bfield_val, 1d0)
      call hfield%field(xstart, ystart, nbx, nby, tpar_test, rmin_test)

      ! write(*,'(A,10Es14.6)') "xmin=", -bpar*hfield%xmin_q2
      ! write(*,'(A,10Es14.6)') "x0=", x0
      write(*,'(A,10Es14.6)') "gamma =", gamma
      write(*,'(A,10Es14.6)') "xstart =", xstart
      write(*,'(A,10Es14.6)') "ystart =", ystart
      write(*,'(A,10Es14.6)') "rmin=", rmin
      ! write(*,'(A,10Es14.6)') "bpar1=", bpar1
      ! write(*,'(A,10Es14.6)') "tpar=", tpar
      ! write(*,'(A,10Es14.6)') "tpar1=", tpar1
      write(*,'(A,10Es14.6)') "rg0 =", rg0
      write(*,'(A,10Es14.6)') "rg0/rmin =", rg0/rmin
      write(*,*) "\n-----"
      ! write(*,'(A,10Es14.6)') "nbx, nby =", nbx, nby
      ! write(*,'(A,10Es14.6)') "bm=", bm
      ! write(*,'(A,10Es14.6)') "v =", y(4:6)
      ! write(*,'(A,10Es14.6)') "sin(pitch_angle) =", vop%ab_sin(y(4:6), bm)
      ! write(*,'(A,10Es14.6)') "y =", y
      ! read(*,*)

   end subroutine setup_info

   subroutine calc_syst(y0)
      type(ode_solution) :: sol, sol_red
      type(ode_solver) :: ode
      integer, parameter :: ns=7
      real(8), intent(in) :: y0(ns)
      real(8) :: t1, t2, eps

      nstep_cond = 0
      write(*,*) "nstep_cond = ", nstep_cond
      t1 = 0d0
      t2 = 4e2
      eps = 1d-6

      call ode%ode_sys(ns, motion_ode)
      call ode%init_cond(t1, y0)
      call ode%solve_to(t2)
      call ode%sol_out_to(sol)
      call ode%toler(eps)
      call ode%init_step(1d-20)
      call ode%add_cond(1, cond1)
      call ode%solve_ode

      write(*,*) "Solved with np = ", sol%np
      ! call sol%reduce_to_npoints(100000, sol_red)
      call write_results(sol)

   end subroutine calc_syst


   subroutine write_results(sol)
      integer :: i
      type(ode_solution), intent(in) :: sol
      type(ode_solution) :: sol_red
      real(8) :: bm(3), em(3)
      real(8) :: acc, loss, ec, loss_pr, acc_norm
      character(500) :: fname
      logical :: start_write

      fname=trim(root_out_dir())//trim(traj_file_name)
      write(*,*) "Trajectory saved in = ", trim(fname)
      open(1, file=fname)

      write(1,'(A1, A39, 20A40)') "#", "Time", "X", "Y", "Z", "VX", "VY", "VZ", "Gamma",&
         "Ec", "Acc", "Loss", "Loss per rot", "Relative loss per rot", &
         "Char loss time", "Acceleration=nu*sqrt(...)"
      start_write = .False.

      call sol_red%set(sol%nsys)
      do i=1,sol%np
         call em_field(sol%y(1:3, i), em, bm)
         call rad_val(sol%y(0, i), sol%y(1:, i), acc, &
            loss, ec, loss_pr, acc_norm)
         if ((.not.start_write)&
            .and.(abs(sol%y(7, i)-sol%y(7, 1))/sol%y(7,1) > 1e-10)) then
            start_write = .True.
         end if

         if (start_write) then
            write(1,'(20Es40.30)') sol%y(0:7, i), ec, acc, &
               loss, loss_pr, loss_pr/sol%y(7, i), sol%y(7, i)/loss, acc_norm
            call sol_red%add_res(sol%y(0,i),sol%y(1:,i))
         end if
      end do

      close(1)
      write(*, *) "Saved points = ", sol_red%np
      call write_int_spectr(sol_red)


   end subroutine write_results


   subroutine write_int_spectr(sol)
      type(ode_solution) :: sol
      real(8), allocatable :: egam(:), spectr(:)
      real(8), allocatable :: ec(:), gm(:), tm(:)
      real(8) :: acc1, loss1, ec1, loss_pr1, acc_norm
      integer :: i, ngam, nc
      character(500) :: fname

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

      fname=trim(root_out_dir())//trim(spectrum_file_name)
      write(*,*) "Spectrum saved in = ", trim(fname)
      open(1,file=fname)

      do i = 1, ngam
         write(1,*) egam(i), dnde_const*spectr(i)
      end do

      close(1)


   end subroutine write_int_spectr


   subroutine show_field
      ! Function for plotting field lines
      use mag_sph_m, only : mfields, norm_field, solid_rot
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

      fname = trim(root_out_dir())//'mag_field/bm_201.dat'
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

      fname = trim(root_out_dir())//'mag_field/bm_3d_dip.dat'
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

         ! write(*,'(A,Es14.6)') "Calc spec for egam = ", egam(iegam)
         res = 0d0
         do iec = 1, nc
            xc = egam(iegam)/ec(iec)
            call f_synch(xc, fsyn)
            fsyn = fsyn/(xc*gm(iec)**2)
            if (fsyn /= fsyn) continue
            if (iec.ne.nc) then
               res_step =  fsyn*(tm(iec+1) - tm(iec))
            else
               res_step = fsyn*(tm(iec) - tm(iec-1))
            end if
            res = res + res_step
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

      call calculate_trajectory
      !  call show_field
      ! call show_field3d

   end subroutine main_calc


end module traj_in_emf_m



program main
   use traj_in_emf_m, only : main_calc
   call main_calc

end program main




