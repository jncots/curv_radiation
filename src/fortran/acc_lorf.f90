module acc_lorf_m
 use phys_const, only : me_ev, e_charge, erg_eV, c_light, pi, psec
 use phys_const, only : mprot
 use utools_mod, only : utools
 
 implicit none
 private
 save


 public :: main_calc


 real(8), parameter :: re=e_charge**2/(me_ev*erg_eV)
 real(8), parameter :: racc0=3*e_charge**(3d0/2d0)/(2*re**2)
 real(8), parameter :: racc0_mp=racc0*(mprot/me_ev)**2/psec
 real(8), parameter :: racc0_me=racc0/psec
 real(8), parameter :: bmax0_mp=(3d0)/(2d0*sqrt(e_charge))*racc0*(mprot/me_ev)**2



 real(8) :: mp_cur, zch_cur



 contains





subroutine acc_cond(bm,eta,res)
 real(8), intent(in) :: bm, eta
 real(8), intent(out) :: res

 res=c_light*zch_cur*e_charge*bm
 res=res/(me_ev*erg_eV)
 res=res*(me_ev/mp_cur)
 res=res*eta

end subroutine acc_cond


subroutine curv_loss(rcurv,gam,res)






end subroutine curv_loss










subroutine racc_fun(rs,bm,res)
 real(8), intent(in) :: rs, bm
 real(8), intent(out) :: res
 real(8) :: rmax

 rmax=racc0_mp*bm**(-3d0/2d0)
 res=rmax*tanh(rs/rmax)

end subroutine racc_fun


subroutine bracc_fun(rs,bm)
 real(8), intent(in) :: rs
 real(8), intent(out) :: bm
 real(8) :: rmax

 bm=(rs/racc0_mp)**(-2d0/3d0)

end subroutine bracc_fun


subroutine bgmax(gam,res)
 real(8), intent(in) :: gam
 real(8), intent(out) :: res

 res=bmax0_mp/gam**2

end subroutine bgmax



subroutine xi_acc(ap,zp,eta,res)
 real(8), intent(in) :: ap, zp, eta
 real(8), intent(out) :: res

 res=ap**(4d0/3d0)
 res=res*zp**(-5d0/3d0)
 res=res*eta**(-1d0/3d0)

end subroutine xi_acc


subroutine plot_racc

type(utools) :: ut
 real(8), allocatable :: bm(:), rs(:), res(:)
 real(8) :: rbm, rbmax
 integer :: i, j, nb

 nb=10
 call ut%grid(bm,1d-6,1d6,nb)
 call ut%grid(res,1d-6,1d6,nb)
 call ut%grid(rs,1d-10,1d12,200)


 open(1,file='racc.dat')

 do i=1,size(rs)
 
  call bracc_fun(rs(i),rbm)
  call bgmax(rs(i),rbmax)
  do j=1,size(bm)

   call racc_fun(rs(i),bm(j),res(j))

  end do
  write(1,*) rs(i), rbm, rbmax, (res(j),j=1,nb)

 end do
 
 close(1)
  

end subroutine plot_racc
 
  

 
subroutine main_calc


 call plot_racc


end subroutine main_calc




end module acc_lorf_m



program main
 use acc_lorf_m, only : main_calc

 call main_calc

end program main
