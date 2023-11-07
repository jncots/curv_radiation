module hyper_field_m
   use phys_const, only : pi
   implicit none
   private
   save

   public :: hyper_field

   type hyper_field
      real(8) :: theta, theta_rad, cost, sint, tant
      real(8) :: sint2, cost2, isint2, icost2
   contains
      procedure :: set_angle, xcoord, field_line, field
   end type hyper_field

contains


   subroutine set_angle(this, theta)
      class(hyper_field) :: this
      real(8), intent(in) :: theta

      this%theta = theta
      this%theta_rad = theta*pi/180

      this%cost = cos(this%theta_rad)
      this%sint = sin(this%theta_rad)
      this%tant = this%sint/this%cost

      this%sint2 = sin(this%theta_rad/2)
      this%cost2 = cos(this%theta_rad/2)
      this%isint2 = 1/this%sint2
      this%icost2 = 1/this%cost2


   end subroutine set_angle


   subroutine xcoord(this, nquad, y, b, x, t)
      class(hyper_field) :: this
      integer :: nquad
      real(8), intent(in) :: y, b
      real(8), intent(out) :: x, t
      real(8) :: aa

      if (nquad==2) then

         aa = b*this%sint2
         t = y/aa
         x = -aa*(1/t + t*this%cost)/this%sint

      else if (nquad==4) then

         aa = b*this%sint2
         t = -y/aa
         x = aa*(1/t + t*this%cost)/this%sint

      else if (nquad==1) then

         aa = b*this%cost2
         t = y/aa
         x = aa*(1/t - t*this%cost)/this%sint

      else if (nquad==3) then

         aa = b*this%cost2
         t = -y/aa
         x = -aa*(1/t - t*this%cost)/this%sint
      else
         t = 0d0
         x = 0d0
      end if

   end subroutine xcoord

   subroutine field_line(this, nquad, t, b, x, y)
      class(hyper_field) :: this
      integer :: nquad
      real(8), intent(in) :: t, b
      real(8), intent(out) :: x, y
      real(8) :: aa

      ! nquad is quadrant number
      if (nquad==2) then

         aa = b*this%sint2
         x = -aa*(1/t + t*this%cost)/this%sint
         y = aa*t

      else if (nquad==4) then

         aa = b*this%sint2
         x = aa*(1/t + t*this%cost)/this%sint
         y = -aa*t

      else if (nquad==1) then

         aa = b*this%cost2
         x = aa*(1/t - t*this%cost)/this%sint
         y = aa*t

      else if (nquad==3) then

         aa = b*this%cost2
         x = -aa*(1/t - t*this%cost)/this%sint
         y = -aa*t

      else
         x=0d0
         y=0d0
      end if

   end subroutine field_line



   subroutine field(this, x, y, nx, ny, t, b)
      class(hyper_field) :: this
      real(8), intent(in) :: x, y
      real(8), intent(out) :: nx, ny, t, b
      real(8) :: max_y, l1, l2, it, t2, it2
      logical :: cond1, cond2, cond

      if (abs(y) < 1d-16) then
         ny = 0
         t = 0
         b = 0
         if (x < 0) then
            nx = -1
         else
            nx = 1
         end if
         return
      end if

      if ((abs(y*(x*this%sint + y*this%cost)) < 1d-16)&
         .or.(abs((x*this%tant + y)) < 1d-16)) then
         l1 = 1d0/sqrt(x**2 + y**2)
         nx = x*l1
         ny = y*l1
         b = 0
         t = 0

         ! if ((nx.ne.nx).or.(ny.ne.ny).or.(nx**2+ny**2 > 2)&
         !    .or.(nx.ne.nx).or.(ny.ne.ny)) then
         !    write(*,"(A,10Es14.6)") "x=", x
         !    write(*,"(A,10Es14.6)") "y=", y
         !    write(*,"(A,10Es14.6)") "nx=", nx
         !    write(*,"(A,10Es14.6)") "ny=", ny
         !    write(*,"(A,10Es14.6)") "l1 = ", l1
         !    write(*,"(A,10Es14.6)") "t = ", t
         !    read(*,*)
         ! end if

         return
      end if

      max_y = -x*this%tant
      cond1 = (y > max_y).and.(y < 0)
      cond2 = (y < max_y).and.(y > 0)
      cond = cond1.or.cond2

      if (this%tant < 0) then
         cond = .not.cond
      end if

      ! Field in 2nd and 4th quadrant
      if (cond) then
         l1 = sqrt(-y*(x*this%sint + y*this%cost))
         b = l1*this%isint2

         t = y/l1
         it = 1d0/t
         t2 = t*t
         it2 = it*it

         l2 = 1d0/sqrt(t2 + it2 - 2*this%cost)
         nx = (it - t*this%cost)*l2
         ny = t*this%sint*l2

         ! if ((nx.ne.nx).or.(ny.ne.ny).or.(nx**2+ny**2 > 2)) then
         !    write(*,"(A,10Es14.6)") "x=", x
         !    write(*,"(A,10Es14.6)") "y=", y
         !    write(*,"(A,10Es14.6)") "l1 = ", l1
         !    write(*,"(A,10Es14.6)") "t = ", t
         ! end if
      else
         ! Field in 1st and 3rd quadrant
         l1 = sqrt(y*(x*this%sint + y*this%cost))
         b = l1*this%icost2

         t = y/l1
         it = 1d0/t
         t2 = t*t
         it2 = it*it

         l2 = 1d0/sqrt(t2 + it2 + 2*this%cost)
         nx = -(it + t*this%cost)*l2
         ny = t*this%sint*l2
      end if

      ! if ((nx.ne.nx).or.(ny.ne.ny).or.(nx**2+ny**2 > 2)&
      !    .or.(nx.ne.nx).or.(ny.ne.ny)) then
      !    write(*,*) "cond1", cond1
      !    write(*,*) "cond2", cond2
      !    write(*,*) "cond", cond
      !    write(*,*) "max_y=", max_y
      !    write(*,*) "y=", y
      !    write(*,"(A,10Es14.6)") "l1**2=", y*(x*this%sint + y*this%cost)
      !    write(*,"(A,10Es14.6)") "x=", x
      !    write(*,"(A,10Es14.6)") "y=", y
      !    write(*,"(A,10Es14.6)") "nx=", nx
      !    write(*,"(A,10Es14.6)") "ny=", ny
      !    write(*,"(A,10Es14.6)") "l1 = ", l1
      !    write(*,"(A,10Es14.6)") "t = ", t
      !    read(*,*)
      ! end if

   end subroutine field


end module hyper_field_m





