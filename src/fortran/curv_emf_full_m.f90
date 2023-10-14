module curv_emf_full_m
 use ode_solver_m, only : ode_solution, ode_solver
 use vector_ops_mod, only : vector_ops
 use phys_const, only : e_charge, erg_eV, c_light, me_ev, pi, mprot

 implicit none
 private
 save

 public :: main_calc


 real(8) :: mc2, emc2, eloss, emc1
 real(8) :: bm0_ch, rc0_ch, chi0_ch(3)
 real(8) :: bm0, rho, chi(3), rot_num, accuracy
 real(8) :: start_pos(3), start_vel(3), start_lorf
 character(2) :: field_type
 character(10) :: particle
 character(500) :: fname(4)
 logical :: loss_on
 
 
 contains



subroutine init_param
 character(500) :: fdir
 real(8) :: rg, nu, eqgam

!=================================================
! Model 1
!=================================================

! fdir='/home/anton/work/project/accl_rev/num_check_full/res1/'

! fname(1)=trim(fdir)//'bvel.dat'
! fname(2)=trim(fdir)//'acc_loss.dat'
! fname(3)=trim(fdir)//'coord.dat'

! field_type='wf' ! wire field ('cc' constant with local curvature)
! loss_on=.true.
! particle='electron' ! electron
! call part_const(particle)

! bm0=1d4 ! G
! rho=1d12 ! cm
! chi=[0.1d0,0.1d0,0.1d0]
! start_lorf=1.3364d9


! call larm_rad(start_lorf,bm0,rg)
! nu=rg/rho
!  
!! call par_drift(chi,nu,start_vel)
!! call perp_drift(chi,nu,start_vel)

! call tot_drift_vel(chi,nu,start_vel)
!! write(*,*) 'nu=' , nu, rg, rho
!! write(*,*) 'b_parallel=' , start_vel(1)
!! write(*,*) 'b_normal=  ' , start_vel(2)
!! write(*,*) 'b_binormal=' , start_vel(3)
!! write(*,*) 'beta=', sqrt(dot_product(start_vel,start_vel))
!! read(*,*)


! start_pos=[1.0d0,0.0d0,0d0]
!! start_vel=[1d0,1d0,1d0] ! velocity - beta, par, norm, binorm
!! start_lorf=1.3293d9
! rot_num=1d0
! accuracy=1d-10

!=================================================
! Model 2
!=================================================

 fdir='/mnt/e/work/project/accl_rev/calc/'
 fname(1)=trim(fdir)//'bvel.dat'
 fname(2)=trim(fdir)//'acc.dat'
 fname(3)=trim(fdir)//'coord.dat'
 fname(4)=trim(fdir)//'acc_comp_loss.dat'

 field_type='cc' ! wire field ('cc' constant with local curvature)
 loss_on=.true.
 particle='proton' ! electron
 call part_const(particle)

 bm0=1d4 ! G
 rho=3d12 ! cm
! chi=[0.1d0,0.1d0,-0.1d0]
 chi=[1d0,0.01d0,0.01d0]
! start_lorf=1.3370d9
 call eqvilib_gam(bm0,rho,chi,eqgam)
 start_lorf=eqgam


 call larm_rad(start_lorf,bm0,rg)
 nu=rg/rho

!tot_drift_vel(chi,nu,bp)
 call tot_drift_vel(chi,nu,start_vel)
! call perp_drift(chi,nu,start_vel)
! start_vel(2)=-0.089d0
! start_vel(3)=-0.1088d0
! start_vel(1)=sqrt(1-start_vel(2)**2-start_vel(3)**2)

! call eqvilib_gam(bm0,rho,chi,eqgam)

! start_vel=[0.1d0,0.1d0,0.1d0]

 start_vel=start_vel+[1d0,3d0,2d0]*1d-2
 write(*,*) 'nu=' , nu, rg, rho
 write(*,'(A,10Es14.6)') ' eqgam=', eqgam
 write(*,*) 'b_parallel=' , start_vel(1)
 write(*,*) 'b_normal=  ' , start_vel(2)
 write(*,*) 'b_binormal=' , start_vel(3)
 write(*,*) 'beta=', sqrt(dot_product(start_vel,start_vel))


 start_pos=[1.0d0,0.0d0,0d0]
 rot_num=0.1d0
 accuracy=1d-10

 rc0_ch=rho
 bm0_ch=bm0
 chi0_ch=chi


end subroutine init_param


subroutine eqvilib_gam(bm,rho,chi,res,nu)
!========================================================================
! gamma_0 is the equlibrium between acceleration and curvature radiation
!========================================================================
 real(8), intent(in) :: bm, rho, chi(3)
 real(8), intent(out) :: res
 real(8), intent(in), optional :: nu
 real(8) :: xi2, xi, xin, rho2, chi1, s1, sigma, nu1
 real(8) :: hc, xp, sig1

 if (present(nu)) then
  nu1=nu
 else
  nu1=0d0
 end if

 xi2=1-chi(2)**2-chi(3)**2
 xi=sqrt(xi2)
 xin=1-chi(2)**2
 rho2=rho**2/(xi2*xin)


  
 xp=chi(1)
 hc=dot_product(chi,chi) ! chi^2
 hc=(1-hc)/2

 sig1=xp/sqrt(hc+sqrt(hc**2+xp**2))

 chi1=chi(1)/xi 
 s1=(1-chi1**2/xi2)*chi(3)
 sigma=sig1+nu1*s1
! sigma=chi1+nu1*s1

 res=3*bm/(2*e_charge)
 res=res*sigma*rho2
 res=res**0.25d0
end subroutine eqvilib_gam


subroutine part_const(pt)
 character(10) :: pt

 if (trim(pt)=='proton') then
  mc2=mprot*erg_eV
 else if (trim(pt)=='electron') then
  mc2=me_ev*erg_eV
 else
  write(*,*) 'Particle type is unknow'
  mc2=me_ev*erg_eV
 end if


 emc2=e_charge/mc2
 emc1=emc2*c_light
 eloss=(2d0/3d0)*mc2*emc2**4

end subroutine part_const




subroutine calc_syst
 type(ode_solution) :: sol
 type(ode_solver) :: ode
 integer, parameter :: ns=7
 real(8) :: y(ns)
 real(8) :: t1, t2, eps, vel_init(3)
 type(vector_ops) :: vop

 call init_param
! rc0_ch=rho
! bm0_ch=bm0
! chi0_ch=chi

 y(1:3)=start_pos*rc0_ch
 vel_init=vop%vec_dir(start_vel)
 call from_cyl_sys(y(1:3),start_vel,y(4:6))
 y(7)=start_lorf
 
 t1=0d0
 t2=2*pi*rho*rot_num
 eps=accuracy

 call ode%ode_sys(ns,motion_ode)
 call ode%init_cond(t1,y)
 call ode%solve_to(t2)
 call ode%sol_out_to(sol)
 call ode%toler(eps)
! call ode%init_step(1d3)
! call ode%add_cond(1,cond1)
 call ode%solve_ode

 call write_results(sol)


end subroutine calc_syst



subroutine write_results(sol)
 type(ode_solution) :: sol
 integer :: i, j
 real(8) :: vloc(3), vloc_an(3)
 real(8) :: acceb, force(3), acc, undim_acc
 real(8) :: em(3), bm(3), nu, gaml, bmval, rg, rho0
 real(8) :: xy(3), dgam_dl, dgam_rho, f2
 real(8) :: chi2_perp, nu_perp, perp_fact
 real(8) :: dgam_dl_acc, dgam_dl_loss, dgam_dl_sigma, em_loc(3), eqgam
 real(8) :: phi0, phic, anl_loss, eps0, ct0, ct1, phic0, phicc

 open(1,file=fname(1))
 open(2,file=fname(2))
 open(3,file=fname(3))
 open(4,file=fname(4))

 phic0=0d0
 ct0=0d0

 do i=1,sol%np
  call cyl_sys(sol%y(1:3,i),sol%y(4:6,i),vloc)
  call emfield(sol%y(1:3,i),em,bm)
  call lor_force(em,bm,sol%y(4:6,i),acceb,force)
  gaml=sol%y(7,i)
  bmval=sqrt(dot_product(bm,bm))
  call larm_rad(gaml,bmval,rg)
  
  f2=dot_product(force,force)
  undim_acc=sqrt(f2)/bmval
  xy=[sol%y(1,i),sol%y(2,i),0d0]
  rho0=sqrt(dot_product(xy,xy))
  nu=rg/rho0

  chi2_perp=chi0_ch(2)**2+chi0_ch(3)**2
  perp_fact=sqrt((1-chi0_ch(2)**2)*(1-chi2_perp))
  nu_perp=nu*perp_fact


  dgam_dl_acc=emc2*acceb
  dgam_dl_loss=-eloss*f2*gaml**2
  dgam_dl=dgam_dl_acc+dgam_dl_loss
  dgam_rho=dgam_dl*rho0/gaml

  call acc_analyt(bmval,chi0_ch,nu,dgam_dl_sigma)

!  call cyl_sys(sol%y(1:3,i),em,em_loc)
!  write(*,*) 'acceb=', dot_product(sol%y(4:6,i),em)/bmval,&
!  dot_product(vloc,em_loc)/bmval
!  write(*,*) 'em=', em/bmval
!  write(*,*) 'em_loc=', em_loc/bmval
!  write(*,*) 'vv_loc=', vloc/(sqrt(dot_product(vloc,vloc)))
!  sqrt(dot_product(sol%y(4:6,i),sol%y(4:6,i))),(sqrt(dot_product(em,em)))
!  read(*,*)


  call eqvilib_gam(bmval,rho0,chi0_ch,eqgam)
  
  eps0=nu*190d0*0.9937687539
  phi0=4.7d0
  phicc=phi_gyro(sol%y(0,i),bmval,eqgam) ! equlibrium
!  phic=phi_gyro(sol%y(0,i),bmval,gaml)   ! current gamma

  ct1=sol%y(0,i)
  phic=phi_gyro_step(ct0,ct1,phic0,bmval,gaml)
  call acc_term_anal(chi0_ch,nu,eps0,phi0,phic,anl_loss) ! calculation of a^2 term
  phic0=phic
  ct0=ct1
!  write(*,'(A,10Es14.6)') 'phic=', phic, phicc
 
  call tot_drift_vel(chi0_ch,nu,vloc_an)
  write(1,*) sol%y(0,i)/(2*pi*rc0_ch), vloc, vloc_an
  write(2,*) sol%y(0,i)/(2*pi*rc0_ch), gaml, undim_acc, nu, dgam_dl, dgam_rho, &
  undim_acc/nu_perp, undim_acc/nu, perp_fact, dgam_dl_acc, dgam_dl_loss, dgam_dl_sigma,&
  eqgam, abs(gaml-eqgam)/gaml, sqrt(anl_loss)/nu_perp 
  write(3,*) sol%y(0,i)/(2*pi*rc0_ch), sol%y(1:3,i)/rc0_ch

  write(4,*) sol%y(0,i)/(2*pi*rc0_ch), sqrt(anl_loss)/nu_perp, undim_acc/nu_perp
 end do
 
 close(4)
 close(3)
 close(2)
 close(1)


end subroutine write_results

subroutine acc_analyt(bm,chi,nu,res)
 real(8), intent(in) :: bm, chi(3), nu
 real(8), intent(out) :: res
 real(8) :: xi

 xi=sqrt(1-chi(2)**2-chi(3)**2)

! res=nu*(chi(3)*(1-chi(1)**2)+2*xi*chi(2)*chi(1))
! res=res+chi(1)*xi
! res=chi(1)+nu*chi(3)
! write(*,'(A,10Es14.6)') 'chi(1), nu*chi(3)=', chi(1)*xi, nu*(chi(3)*(1-chi(1)**2)+2*xi*chi(2)*chi(1))
! write(*,'(A,10Es14.6)') 'xi, 2*xi*chi(2)*chi(1)', xi, 2*xi*chi(2)*chi(1)/chi(3)*(1-chi(1)**2)

! write(*,'(A,10Es14.6)') 'chi(3)*(1-chi(1)**2), 2*xi*chi(2)*chi(1)', chi(3)*(1-chi(1)**2), 2*xi*chi(2)*chi(1)
! read(*,*)
! res=chi(1)/xi*(1-((1-xi**2)/(2*xi**2))*(chi(1)/xi)**2)+nu*chi(3)

 res=dot_product(chi,chi)
 res=(1-res)/2
 res=sqrt(res**2+chi(1)**2)+res
 res=chi(1)/sqrt(res)
 res=res+nu*(chi(3)+(chi(1)/xi)*(chi(3)*(chi(1)/xi)+(xi**2-1)*chi(2)))
 res=res*emc2*bm
end subroutine acc_analyt



subroutine test_g_min_max
 use utools_mod, only : utools
 real(8) :: bm, gam, ct
 real(8) :: eps, chi(3), rho, phi, epss
 real(8) :: gmin, gmax
 type(utools) :: ut
 real(8), allocatable :: xct(:)
 integer :: i, ndev
 character(500) :: fname

 call init_param

! ct=rc0_ch/1d3
 bm=bm0_ch
 rho=rc0_ch
 chi=chi0_ch

 call eqvilib_gam(bm,rho,chi,gam)


! gam=1d7

! bm0_ch, rc0_ch, chi0_ch(3)

! write(*,'(A,Es14.6,A)') 'ct=', ct, ' cm'
! write(*,'(A,Es14.6,A)') 'bm=', bm, ' Gauss'
! write(*,'(A,Es14.6,A)') 'gam=', gam, ' '

! write(*,'(A,Es14.6,A)') 'phi=', phi_gyro(ct,bm,gam), ' times'
 
 
 eps=1d-1
! phi=phi_gyro(ct,bm,gam)
! write(*,*) 'chi=', chi
! write(*,'(A,Es14.6,A)') 'epsgam=', epsgam(eps,chi,bm,rho,phi), ' times'
 call ut%grid(xct,rc0_ch/1d5, rc0_ch/1d1,1000)
 

 fname='/mnt/e/work/project/accl_rev/calc/grange.dat' 
!  fname='/home/anton/work/project/accl_rev/num_check_full/gam_range/grange.dat'
 open(newunit=ndev,file=fname)

 do i=1,size(xct)
  phi=phi_gyro(xct(i),bm,gam)
  epss=epsgam(eps,chi,bm,rho,phi)
  call gam_min_max(epss,gam,gmin,gmax) 
!  write(*,'(A,10Es14.6)') 'xct, gmin, gmax=', xct(i), gmin, gmax

  write(ndev,*)  xct(i), gmin, gmax
 end do



end subroutine test_g_min_max




function phi_gyro(ct,bm,gam)
!===========================================================
! Gyrophase
!===========================================================
 real(8), intent(in) :: bm, gam, ct
 real(8) :: phi_gyro

 phi_gyro=emc2*bm*ct/gam

end function phi_gyro


function phi_gyro_step(ct0,ct1,phi0,bm,gam)
!===========================================================
! Gyrophase small step
!===========================================================
 real(8), intent(in) :: bm, gam, ct0, ct1, phi0
 real(8) :: phi_gyro_step

 phi_gyro_step=phi0+emc2*bm*(ct1-ct0)/gam

end function phi_gyro_step


function epsgam(eps,chi,bm,rho,phi)
!===========================================================
! Calculations of the 
! epsilon_{gamma}=1/2*epsilon_{hat}*rho_perp/r_cycl
! for the equation 
! |gamma^2 +- epsilon_{gamma}*gamma|-{gamma_0}^2=0 
!===========================================================
 real(8), intent(in) :: eps, chi(3), bm, rho, phi
 real(8) :: epsgam
 real(8) :: xin2, xi2, xin, xi, rho_perp, irc


 xin2=1-chi(2)**2
 xi2=xin2-chi(3)**2

 xin=sqrt(xin2)
 xi=sqrt(xi2)
 rho_perp=rho/(xin*xi)

 irc=emc2*bm   ! 1/r_cycl, where r_cycl=mc^2/eB - cyclotron radius

 epsgam=sqrt(1+chi(1)**2)*xi2*exp(-chi(1)*phi)
 epsgam=epsgam*eps*rho_perp*irc/2

end function epsgam


subroutine gam_min_max(epss,g0,gmin,gmax)
 real(8), intent(in) :: epss, g0
 real(8), intent(out) :: gmin, gmax
 real(8) :: alpha, al2, sal

 alpha=epss/g0 ! epss is synchrotron epsilon
 al2=alpha**2
 
 if (alpha>1d0) then
  gmin=g0/(sqrt(al2+1)+alpha)
  gmax=g0/(sqrt(al2-1)+alpha)
 else
  sal=sqrt(al2+1)
  gmin=g0/(sal+alpha)
  gmax=g0/(sal-alpha)
 end if

end subroutine gam_min_max



subroutine par_drift(chi,nu,bp)
 real(8), intent(in) :: chi(3), nu
 real(8), intent(out) :: bp(3)
 real(8) :: hh, gg, bp2, sigma

 hh=(1-chi(1)**2)/2
 gg=nu**2+chi(1)**2
 bp2=1/(hh+sqrt(hh**2+gg))
 bp(1)=sqrt(bp2)
 sigma=chi(1)*bp(1)
 bp(3)=nu*bp2/(1+sigma**2)
 bp(2)=-sigma*bp(3)


! write(*,'(A,10Es14.6)') 'bp=', bp

end subroutine par_drift




subroutine perp_drift(chi,nu,bp)
 real(8), intent(in) :: chi(3), nu
 real(8), intent(out) :: bp(3)
 real(8) :: dd, nn

 nn=nu-chi(2)
 dd=2/(sqrt(1+4*nu*nn)+1)

 bp(2)=chi(3)*dd
 bp(3)=nn*dd
 bp(1)=sqrt(1-bp(2)**2-bp(3)**2)

! write(*,'(A,10Es14.6)') 'bp=', bp

end subroutine perp_drift


subroutine tot_drift_vel(chi,nu,bp)
!
! Drift velocity for general case of parallel and perpendicular
! electric field in case of nu<<1
!
 real(8), intent(in) :: chi(3), nu
 real(8), intent(out) :: bp(3) ! velosity
 real(8) :: xp, xn, xb, xi, xpi, xi2
 real(8) :: b1, b1n, b2, b2n, b3, b3n
 real(8) :: hc, xp2, bb1, b12

 xp=chi(1)
 xn=chi(2)
 xb=chi(3)

 xi2=1-chi(2)**2-chi(3)**2
 if (xi2>0d0) then
  xi=sqrt(xi2)
 else
  xi=0d0
 end if

! bp(1) is b_parallel
! bp(2) is b_normal
! bp(3) is b_binormal

 xpi=xp/xi

 hc=dot_product(chi,chi) ! chi^2
 hc=(1-hc)/2
 xp2=xp**2

 b12=hc+sqrt(hc**2+xp2)
 b1=sqrt(b12)
 bb1=b1/(xp2+b12)
 
 b2=(xn*xp+xb*b1)*bb1
 b3=(xb*xp-xn*b1)*bb1
  


! b1=xi+((1-xi2)/(2*xi))*xpi**2 ! expantion for small xp
 b1n=xi*xn+((1-xi2)/xi)*xb*xpi+((1-3*xi2)/(2*xi))*xn*xpi**2

! b2=xb+(xn-xb*xpi)*xpi ! expantion for small xp
 b2n=xn*xb-(2*xb**2+xi2)*xpi-((1+3*xi2)/xi2)*xn*xb*xpi**2

! b3=-xn+(xb+xn*xpi)*xpi ! expantion for small xp
 b3n=1-xn**2+2*xn*xb*xpi-(xn**2-(2+xb**2/xi2))*xpi**2


 bp(1)=b1+nu*b1n 
 bp(2)=b2+nu*b2n
 bp(3)=b3+nu*b3n

! bp(1)=xi*(1+xn*nu)+nu*xpi*xb/xi+xpi**2/(2*xi)
! bp(2)=xb*(1+xn*nu)+xpi*(xn-nu-xpi*xb)
! bp(3)=-xn+(1-xn**2)*nu+xpi*(xb+xpi*xn)

end subroutine tot_drift_vel


subroutine acc_term_anal(chi,nu,eps,phi0,phi,res)
 real(8), intent(in) :: chi(3), nu, eps, phi0, phi
 real(8), intent(out) :: res
 real(8) :: chi_perp2, chi_perp, xi2, xi, phic, sxi
 real(8) :: ef, eta, nu_perp

 chi_perp2=chi(2)**2+chi(3)**2
 chi_perp=sqrt(chi_perp2)
 xi2=1-chi_perp2
 xi=sqrt(xi2)

 nu_perp=nu*sqrt(xi2*(1-chi(2)**2))
 eta=eps/nu_perp
 eta=eta*xi2*sqrt(1+chi(1)**2)
 phic=phi-phi0
 sxi=sin(xi*phic)
 eta=eta/(1-eps*chi_perp*sxi)
! ef=(chi(1)/xi+chi(3)*nu/2)*phic
! ef=(chi(1)/xi)*phi
 ef=(chi(1)/xi)*(phi)**(0.9)
 eta=eta*exp(-ef)

! write(*,*) chi(1)
! read(*,*)

 res=(1+2*eta*sxi+eta**2)
 res=res*nu_perp**2

end subroutine acc_term_anal


subroutine cyl_sys(x,v,vc)
 real(8), intent(in) :: x(3), v(3)
 real(8), intent(out) :: vc(3)
 real(8) :: r0, cc, ss
 

 r0=1d0/sqrt(x(1)**2+x(2)**2)
 cc=x(1)*r0
 ss=x(2)*r0

 vc(2)=-(v(1)*cc+v(2)*ss)     !-e_r, e_norm
 vc(1)=-v(1)*ss+v(2)*cc       ! e_phi, e_par
 vc(3)=v(3)                   ! e_z,  e_binorm

end subroutine cyl_sys


subroutine from_cyl_sys(x,v,vc)
 real(8), intent(in) :: x(3), v(3)
 real(8), intent(out) :: vc(3)
 real(8) :: r0, cc, ss
 

 r0=1d0/sqrt(x(1)**2+x(2)**2)
 cc=x(1)*r0
 ss=x(2)*r0

 vc(1)=-v(2)*cc-v(1)*ss    ! v(2) is along normal
 vc(2)=-v(2)*ss+v(1)*cc    ! v(1) is along tangent
 vc(3)=v(3)                ! v(3) is along binormal

end subroutine from_cyl_sys



subroutine lor_force(em,bm,vel,acc,force)
 real(8), intent(in) :: em(3), bm(3), vel(3)
 real(8), intent(out) :: acc, force(3)
 type(vector_ops) :: vop

 acc=dot_product(vel,em)
 force=em-acc*vel
 force=force+vop%cross_prod(vel,bm)

end subroutine lor_force

subroutine larm_rad(gam,bm,res)
 real(8) :: gam, bm
 real(8) :: res
 res=mc2*gam/(e_charge*bm)
end subroutine larm_rad


subroutine emfield(x,em,bm)
 real(8), intent(in) :: x(3)
 real(8), intent(out) :: em(3), bm(3)
 
 if (field_type=='wf') then
  call current_field(x,em,bm)
 else if (field_type=='cc') then
  call const_curvf(x,em,bm) 
 else
  write(*,*) 'Field is not defined'
  em=0d0
  bm=0d0
 end if
 
end subroutine emfield


subroutine const_curvf(x,em,bm)
!=========================================
! Magnetic field of constant strength
!=========================================
 real(8), intent(in) :: x(3)
 real(8), intent(out) :: em(3), bm(3)
 real(8) :: bm_val
 real(8) :: el_dir(3)
 real(8) :: bm_cyl_dir(3), bm_dir(3)
 
 bm_cyl_dir=[1d0,0d0,0d0]
 call from_cyl_sys(x,bm_cyl_dir,bm_dir)

 bm_val=bm0_ch
 bm=bm_val*bm_dir             ! B=B0*R0/r

 call from_cyl_sys(x,chi0_ch,el_dir)
 em=bm_val*el_dir

end subroutine const_curvf


subroutine current_field(x,em,bm)
!=========================================
! Magnetic field created by straight wire
!=========================================
 real(8), intent(in) :: x(3)
 real(8), intent(out) :: em(3), bm(3)
 real(8) :: bm_norm, r0, bm_val
 real(8) :: el_dir(3)
 real(8) :: bm_cyl_dir(3), bm_dir(3)
 
 bm_norm=bm0_ch*rc0_ch
 bm_cyl_dir=[1d0,0d0,0d0]
 call from_cyl_sys(x,bm_cyl_dir,bm_dir)

 r0=sqrt(x(1)**2+x(2)**2)
 bm_val=bm_norm/r0
 bm=bm_val*bm_dir             ! B=B0*R0/r

 call from_cyl_sys(x,chi0_ch,el_dir)
 em=bm_val*el_dir

end subroutine current_field


subroutine motion_ode(x,y,dydx)
 real(8), intent(in) :: x, y(:)
 real(8), intent(out) :: dydx(:)
 real(8) :: em(3), bm(3), force(3), acc, f2
  

 call emfield(y(1:3),em,bm)
 call lor_force(em,bm,y(4:6),acc,force)
 f2=dot_product(force,force)

 dydx(1:3)=y(4:6)
 dydx(4:6)=(emc2/y(7))*force

 if (loss_on) then
  dydx(7)=emc2*acc-eloss*f2*y(7)**2
 else
  dydx(7)=0d0
 end if

end subroutine motion_ode


subroutine main_calc

 call calc_syst

 call test_g_min_max

end subroutine main_calc


end module curv_emf_full_m



program main
 use curv_emf_full_m, only : main_calc

 call main_calc


end program main




