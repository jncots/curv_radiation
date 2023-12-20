module hyper_field_test
   use hyper_field_m, only : hyper_field
   use utools_mod, only : utools
   use util_test, only : cwd_parent_dir

   implicit none
   private
   save

   public :: main_calc

contains


   subroutine test_field_lines
      type(utools) :: ut
      type(hyper_field) :: hfield
      real(8), allocatable :: tgrid(:)
      real(8) :: b, theta, tmax
      real(8) :: x1, y1, x2, y2, x3, y3, x4, y4
      integer :: i, nt
      character(500) :: fname


      theta = 100d0
      call hfield%set_angle(theta)
      nt = 100
      tmax = 1d2
      call ut%grid(tgrid, 1/tmax, tmax, nt, 'log')

      fname = trim(cwd_parent_dir(3))//'/results/data/hyper_field.dat'
      open(1,file=fname)

      b = 0.4
      do i=1,nt
         call hfield%field_line(1, tgrid(i), b, x1, y1)
         call hfield%field_line(2, tgrid(i), b, x2, y2)
         call hfield%field_line(3, tgrid(i), b, x3, y3)
         call hfield%field_line(4, tgrid(i), b, x4, y4)
         write(1,*) tgrid(i), b, x1, y1, x2, y2, x3, y3, x4, y4
      end do

   end subroutine test_field_lines


   subroutine test_xcoord
      type(utools) :: ut
      type(hyper_field) :: hfield
      real(8), allocatable :: ygrid(:), bgrid(:), theta_grid(:)
      real(8) :: b, theta, tmax
      real(8) :: x1, y1, x, t, tol
      integer :: i, j, k, k1, ny, nb, nquad, ntheta

      ntheta = 100
      call ut%grid(theta_grid, 0d0, 180d0, ntheta, 'lin')

      do k=1,ntheta
         theta = theta_grid(k)
         do k1 = 1,4
            nquad = k1
            tol = 1d-14
            call hfield%set_angle(theta)
            ny = 10
            nb = 10

            call ut%grid(theta_grid, 0d0, 180d0, 100, 'lin')
            call ut%grid(ygrid, 1d-5, 20d0, ny, 'lin')
            call ut%grid(bgrid, 1d-5, 20d0, nb, 'lin')


            do i=1,nb
               do j=1,ny
                  call hfield%xcoord(nquad, ygrid(j), bgrid(i), x, t)
                  call hfield%field_line(nquad, t, bgrid(i), x1, y1)

                  if (abs(y1 - ygrid(j)) > tol) then
                     write(*, *) "ygrid=", ygrid(j), " y=", y1
                  end if

                  if (abs(x1 - x) > tol) then
                     write(*, *) "x=", x, " x1=", x1
                  end if

               end do
            end do
         end do
      end do

   end subroutine test_xcoord


   subroutine test_field
      type(utools) :: ut
      type(hyper_field) :: hfield
      real(8), allocatable :: xgrid(:), ygrid(:)
      real(8) :: theta, vx, vy, t, b
      real(8) :: xmin, xmax, ymin, ymax
      integer :: i, j, nx, ny
      character(500) :: fname


      theta = 120d0
      call hfield%set_angle(theta)

      nx = 100
      ny = 100
      xmin = -10
      xmax = 10
      ymin = -10
      ymax = 10

      call ut%grid(xgrid, xmin, xmax, nx, 'lin')
      call ut%grid(ygrid, ymin, ymax, ny, 'lin')

      fname = trim(cwd_parent_dir(3))//'/results/03_test_ypoint_field/ypoint_vector_field.dat'
      open(1,file=fname)

      do i = 1, nx
         do j = 1, ny
            call hfield%field_ypoint(xgrid(i), ygrid(j), vx, vy, t, b)
            write(1, *) xgrid(i), ygrid(j), vx, vy, t, b
         end do
      end do

   end subroutine test_field


   subroutine test_yfield
      ! Output vector field of y-point field
      type(utools) :: ut
      type(hyper_field) :: hfield
      real(8), allocatable :: xgrid(:), ygrid(:)
      real(8) :: theta, vx, vy, t, b
      real(8) :: xmin, xmax, ymin, ymax
      integer :: i, j, nx, ny
      character(500) :: fname


      theta = 30d0
      call hfield%set_angle(theta)

      nx = 100
      ny = 100
      xmin = -10
      xmax = 10
      ymin = -10
      ymax = 10

      call ut%grid(xgrid, xmin, xmax, nx, 'lin')
      call ut%grid(ygrid, ymin, ymax, ny, 'lin')

      fname = trim(cwd_parent_dir(3))//'/results/data/03_test_ypoint_field/ypoint_vector_field.dat'
      open(1,file=fname)

      do i = 1, nx
         do j = 1, ny
            call hfield%field_ypoint(xgrid(i), ygrid(j), vx, vy, t, b)
            write(1, *) xgrid(i), ygrid(j), vx, vy, t, b
         end do
      end do

   end subroutine test_yfield



   subroutine main_calc

      !   call test_field_lines
      !   call test_xcoord
      call test_yfield

   end subroutine main_calc


end module hyper_field_test



program main
   use hyper_field_test, only : main_calc

   call main_calc


end program main
