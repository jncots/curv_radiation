module mag_sph_m
 use field_3dlm, only : field_3dl
 use vector_rot3d_m, only : vector_rot3d
 use mag_dipole_m, only : mag_dipole
 use ret_dipole_mod, only : ret_dipole
 use phys_const, only : c_light, pi


 implicit none
 private
 save

 public :: norm_field, mfields, rot_incl, solid_rot
 public :: mf_rlc, mf_omega

 real(8) :: bnorm=1d0, al_incl, mf_rlc, mf_omega
 real(8), parameter :: nom0(3)=[0d0,0d0,1d0] ! along z diretion
 real(8) :: nom(3) ! axis of rotation

 type(field_3dl) :: fd
 type(mag_dipole) :: mdip
 type(ret_dipole) :: rdip
 type(vector_rot3d) :: rf_dip, rf_dat, rf_solid, rf_incl


 contains



subroutine set_field
 integer :: i
! Reading data
 fd%fname='/home/antonpr/work/pulsar_rad/data_files/CrabLowRes1.35_ri.dat'
 call fd%read

 ! do i=1,size(fd%x)
 !   write(*,'(I5,10Es14.6)') i, fd%x(i), fd%y(i), fd%z(i)
 ! end do
 ! read(*,*)

! Rotate data
 call rf_dat%set(1)
 call rf_dat%set_rot(1,3,126d0)
!  call rf_dat%set_rot(1,3,0d0)
! call rf_dat%set_rot(1,2,90d0)
! call rf_dat%set_rot(1,3,-140d0)
! call rf_dat%set_rot(1,3,0d0)

! Retarded dipole 
  rdip%bm_star=2.1d12
!  rdip%bm_star=1d11
  rdip%r_star=1.5d6
  call rdip%set
  call rdip%inclin(90d0,'deg')
  call rdip%phase(0d0,'deg')
  mf_omega=rdip%omega
  mf_rlc=rdip%r_unit

! Rotate dipole
 call rf_dip%set(2)
 call rf_dip%set_rot(1,2,0d0) ! mu along z
 call rf_dip%set_rot(2,3,0d0) ! mu along z
! call rf_dip%calc_rm


! Dipole setting and rotating
! mdip%r_star=1d6 ! 10 km
 mdip%r_star=1.5d6 ! 15 km
 mdip%bm_star=2.1d12
! mdip%bm_star=1d10

 call mdip%set('lc')
 mf_omega=mdip%omega
 mf_rlc=mdip%r_unit

 write(*,'(A, Es14.6)') "r_surface = ", mdip%rsu

! Rotate dipole
 call rf_dip%set(1)
 call rf_dip%set_rot(1,2,90d0)
 
! call rf_dip%calc_rm
! write(*,'(3Es14.6)') rf_dip%trm
! read(*,*)


end subroutine set_field


subroutine norm_field
 real(8) :: x(3), bdat(3), bdip(3), bdat0, bdip0

 call set_field


 x=[0d0,1d0,0d0]
 x=x/sqrt(dot_product(x,x))
 x=0.35d0*x

 bdat=bm_dat(x)
 bdip=bm_dip(x)

 bdat0=sqrt(dot_product(bdat,bdat))
 bdip0=sqrt(dot_product(bdip,bdip))
 bnorm=bdip0/bdat0
! write(*,*) 'bnorm=', bnorm

! bdat=bdat/sqrt(dot_product(bdat,bdat))
! bdip=bdip/sqrt(dot_product(bdip,bdip))
!
! write(*,*) bdat
! write(*,*) bdip
! write(*,*) 'ang=',(acos(dot_product(bdip,bdat)))*180/pi


end subroutine norm_field



function bm_dat(x)
 real(8), intent(in) :: x(3)
 real(8) :: bm_dat(3)
 real(8) :: x1(3), x0(3), b0(3)

! Shift of data
 x1=x
 x1(1)=x1(1)+0.035d0
 x1(2)=x1(2)-0.0245d0
 x1(3)=x1(3)+0.04d0

 call rf_dat%rotate(x1,x0)
 call fd%cfield(x0(1),x0(2),x0(3),b0)
!  call fd%cfield_cubic(x0(1),x0(2),x0(3),b0)
 call rf_dat%rotate_back(b0,bm_dat)
 bm_dat=bnorm*bm_dat

end function bm_dat



function bm_ret(x)
 real(8), intent(in) :: x(3)
 real(8) :: bm_ret(3)
 real(8) :: x1(3), x0(3), b0(3)


 call rdip%calc(x,bm_ret)
! x1=x
! call rf_dip%rotate_back(x1,x0)
! call rdip%calc(x0,b0)

! call rf_dip%rotate(b0,bm_ret)

end function bm_ret



function bm_dip(x)
 real(8), intent(in) :: x(3)
 real(8) :: bm_dip(3)
 real(8) :: x1(3), x0(3), b0(3)

 x1=x
 call rf_dip%rotate_back(x1,x0)
 call mdip%calc(x0,b0)

 call rf_dip%rotate(b0,bm_dip)

end function bm_dip




subroutine mfields(x,bmf,emf,typ)
 real(8), intent(in) :: x(3)
 real(8), intent(out) :: bmf(3)
 real(8), optional, intent(out) :: emf(3)
 real(8), parameter :: rjoint=0.35d0
 real(8) :: xa, x0(3), bm0(3), rb, nb
 character(3), intent(in), optional :: typ
 character(3) :: typ_field

 if (present(typ)) then
  typ_field=typ
 else
  typ_field='dip'
 end if

 call rf_solid%rotate_back(x,x0)
 xa=sqrt(dot_product(x0,x0))

 if (xa<=rjoint) then
  bm0=bm_dip(x0)
!  bm0=bm_ret(x0)

 else
  if (typ_field=='fff') then
   bm0=bm_dat(x0)
  else
   bm0=bm_dip(x0)
  end if
!  bm0=bm_ret(x0)
 end if
 call rf_solid%rotate(bm0,bmf)
! If x is in light cylinder units, then
! knowing 'nom' - direction of rotation
! one can find E=-(vxB)/c as:

! write(*,*) nom
! read(*,*)

 if (present(emf)) then
  rb=dot_product(x,bmf)
  nb=dot_product(nom,bmf)
  emf=rb*nom-nb*x
 end if

end subroutine mfields


subroutine rot_incl(alpha)
!==========================================
! Setting angle 'al_incl' of inclination
! of rotational axis relative z axis
! in zx plane and calculation of the
! direction of angular velosity of rotation
! 'nom'
!==========================================
 real(8), intent(in) :: alpha

 al_incl=alpha
 call rf_incl%set(1)
 call rf_incl%set_rot(1,2,al_incl)
 call rf_incl%rotate(nom0,nom)

end subroutine rot_incl


subroutine solid_rot(phase)
!==========================================
! Calculation of transformations due to
! solid rotation of pulsar
!==========================================
 real(8), intent(in) :: phase

 write(*,*) "al_incl = ", al_incl
 call rf_solid%set(2)
 call rf_solid%set_rot(1,2,0d0)
 call rf_solid%set_rot(2,3,phase)

end subroutine solid_rot



end module mag_sph_m
