module tricubic_interpol_m
     use utools_mod, only : utools
     use util_test, only : cwd_parent_dir
     use intpol_mod, only : arr_ind_short   


     implicit none
     private
     save
    
     public :: main_calc  
          
    
     type(utools) :: ut
     real(8) :: t3c_coef(4, 4)
     
     
     contains
    
    
     subroutine init_t3c_coef

     t3c_coef = 0.5d0*transpose(reshape((/0, 2, 0, 0,&
                                        -1, 0, 1, 0,&
                                        2, -5, 4, -1,&
                                        -1, 3, -3, 1/),&
                                        shape(t3c_coef)))

     end subroutine init_t3c_coef
     
     
     function t3cint(u, p)
          implicit none
          real(8), intent(in) :: u, p(4)
          real(8) :: t3cint
          real(8) :: u2, u3, uvec(4)
       
          u2 = u*u
          u3 = u2*u
          uvec = [1d0, u, u2, u3]

          write(*, *) "uvec = ", uvec
       
          uvec = matmul(uvec, t3c_coef)

          write(*, *) "uvec = ", uvec
          t3cint = dot_product(uvec, p)
          write(*, *) "t3cint  = ", t3cint 
       end function t3cint

     function tcint_3d(upoint, fun_vals)
          real(8), intent(in) :: upoint(3), fun_vals(4, 4, 4)
          real(8) :: tcint_3d
          real(8) :: tz(4, 4), ty(4)
          integer :: i, j

          write(*,*) "fun_vals = ", fun_vals

          do i=1,4
               do j=1,4
               tz(i, j) = t3cint(upoint(3),fun_vals(i, j, :))
               end do
          end do
          
          write(*,*) "tz = ", tz

          do i=1,4
               ty(i) = t3cint(upoint(2), tz(i, :))
          end do

          write(*,*) "ty = ", ty
          tcint_3d = t3cint(upoint(1), ty)
     end function tcint_3d
     
     function test_3dfun(x, y, z)
     real(8), intent(in) :: x, y, z
     real(8) :: test_3dfun
     
     test_3dfun = (x**(2d0/3d0)*(y-5d0) - y**2)*z**(2d0)

     end function test_3dfun     

     subroutine test_interpol  
        real(8), allocatable :: xgrid(:), ygrid(:), zgrid(:)
        real(8), allocatable :: fgrid(:, :, :)
        real(8) :: dx, dy, dz
        real(8) :: x, y, z, upoint(3), fun_vals(4, 4, 4)
        real(8) :: int_res, real_res 
        integer :: nn, nx, ny, nz
        integer :: ix, iy, iz, i, j, k
        integer :: ix1, ix2, iy1, iy2, iz1, iz2


        write(*,*) "OK1"
        nn = 10
        nx = nn
        ny = nn
        nz = nn
        write(*,*) "OK2"
     call ut%grid(xgrid, 0d0, 10d0, nx, 'lin')
     call ut%grid(ygrid, 0d0, 10d0, ny, 'lin')
     call ut%grid(zgrid, 0d0, 10d0, nz, 'lin') 
     allocate(fgrid(nx, ny, nz))
     write(*,*) "OK3"
     do ix=1,nx
          do iy=1,ny
               do iz=1,nz

                  fgrid(ix, iy, iz) = &
                   test_3dfun(xgrid(ix), ygrid(iy), zgrid(iz))   
                    
               end do     
          end do          
     end do
     write(*,*) "OK4"
     
     x = 8.33546d0
     y = 2.496585d0
     z = 5.2904d0


     call arr_ind_short(1,nx,xgrid, x,ix1, ix2)
     call arr_ind_short(1,ny,ygrid, y,iy1, iy2)
     call arr_ind_short(1,nz,zgrid, z,iz1, iz2)
     write(*,*) "OK5"
     dx = xgrid(2)-xgrid(1)
     dy = ygrid(2)-ygrid(1)
     dz = zgrid(2)-zgrid(1)
     write(*,*) "OK6"

     write(*,*) "OK7"
     upoint = [(x-xgrid(ix1))/dx, (y-ygrid(iy1))/dy, (z-zgrid(iz1))/dz]
     write(*,*) ix1, iy1, iz1

     do i = 1, 4
          do j = 1, 4
               do k = 1, 4
                  fun_vals(i, j, k) = &
                  fgrid(ix1-2+i, iy1-2+j, iz1-2+k)
               end do
          end do
     end do              
     write(*,*) "OK8"
     write(*,*) fun_vals
     write(*,*) "upoint = ", upoint
     int_res = tcint_3d(upoint, fun_vals)
     real_res = test_3dfun(x, y, z)

     write(*,*) "OK9"

     write(*,*) "int_res  = ", int_res
     write(*,*) "real_res = ", real_res
     write(*,*) "rel_error = ", abs(real_res - int_res)/real_res

     write(*,*) "OK10"

     end subroutine test_interpol     
    
    subroutine main_calc
    
     call init_t3c_coef
     call test_interpol 
    
    end subroutine main_calc     

end module tricubic_interpol_m    

program main
     use tricubic_interpol_m, only : main_calc     
     call main_calc  
end program main

! program test
!      use utools_mod, only : utools
!      use util_test, only : cwd_parent_dir
!      implicit none
!      type(utools) :: ut
!      real(8) :: m1(2,2), m2(2, 2), v(2)
!      real(8) :: t3c_coef(4, 4)
!      real(8) :: t3cint, val
!      real(8), allocatable :: ux(:)
!      integer :: i, nux
!      character(500) :: fname
!      real(8) :: u0(3), u1(3), u2(3), u3(3)

!      m1 = transpose(reshape((/1, 1, 3, 1/), shape(m1)))

!      t3c_coef = 0.5d0*transpose(reshape((/0, 2, 0, 0,&
!                                   -1, 0, 1, 0,&
!                                    2, -5, 4, -1,&
!                                   -1, 3, -3, 1/),&
!                                    shape(t3c_coef)))

!      ! do i=1, 4
!      !  write(*, *) t3_coef(i,:)

!      ! end do 
       
!      u0 = [1, 1, 1]                              
!      u1 = [1, 2, 3]
!      u2 = u1*u1
!      u3 = u2*u1
!      write(*,*) [u0, u1, u2, u3]   
!      ! write(*,*) uu*uu*uu               

!      ! nux = 100
!      ! call ut%grid(ux, 0d0, 1d0, nux,'lin')
     
!      ! fname=trim(cwd_parent_dir(2))//'/data/t3_interp.dat'
!      ! open(1,file=fname)
     
!      ! do i=1,100  
          
!      !    val = t3cint(t3c_coef, ux(i), [4d0, 5d0, 3d0, -3d0])  
!      !    write(1, *) ux(i)*100, val
        
!      ! end do   

! end

! function t3cint(t3c_coef, u, p)
!    implicit none
!    real(8), intent(in) :: t3c_coef(4, 4), u, p(4)
!    real(8) :: t3cint
!    real(8) :: u2, u3, uvec(4)

!    u2 = u*u
!    u3 = u2*u
!    uvec = [1d0, u, u2, u3]

!    uvec = matmul(uvec, t3c_coef)
!    t3cint = dot_product(uvec, p)

! end function t3cint


! function tcint_3d(upoint, fun_vals)

!      tz = 0d0
!      do i=1,4
!        do j=1,4
!          tz(i, j) = t3cint(upoint(3),fun_vals(i, j, :))
!        end do
!      end do
     
!      ty = 0d0
!      do i=1,4
!            tz(i) = t3cint(upoint(2), tz(i, :))
!      end do

!      ty = 0d0

!      val = t3cint(upoint(1), ty)








! end function tcint_3d    
   
   



