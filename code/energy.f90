
module ENERGY
interface

subroutine calc_div_B(h,mu_hat,div_B)
use parameters
	implicit none
	complex(prc),pointer :: h(:,:,:,:)
	complex(prc),pointer :: mu_hat(:,:,:,:,:)
	complex(prc),pointer :: div_B(:,:,:)
end subroutine calc_div_B

subroutine calc_div_D(e,div_D)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)	
	complex(prc),pointer :: div_D(:,:,:)
end subroutine calc_div_D

subroutine calc_energy_density(E_cur,E_prev,H_prev,rho,eps_inv,zi,zf,bbox)
use parameters
	implicit none
	complex(prc),pointer :: E_cur(:,:,:,:)
	complex(prc),pointer :: E_prev(:,:,:,:)
	complex(prc),pointer :: H_prev(:,:,:,:)	
	complex(prc),pointer                    :: eps_inv(:,:,:,:,:)	
	real(prc),pointer :: rho(:,:,:)
	integer,intent(in)::zi
    integer,intent(in)::zf
	real(prc),pointer ::bbox(:,:,:)
end subroutine calc_energy_density

subroutine calc_current(J,e_cur,h_cur,ix,iy,iz)
use parameters
	implicit none
	complex(prc),intent(out) :: J(3)
	complex(prc),pointer :: e_cur(:,:,:,:),h_cur(:,:,:,:)
	integer,intent(in) :: ix,iy,iz
end subroutine calc_current

subroutine calc_current_OLD(J,e_cur,h_cur,h_prev,ix,iy,iz)
use parameters
	implicit none
	complex(prc),intent(out) :: J(3)
	complex(prc),pointer :: e_cur(:,:,:,:),h_cur(:,:,:,:),h_prev(:,:,:,:)
	integer,intent(in) :: ix,iy,iz
end subroutine calc_current_OLD
end interface

end module ENERGY


! -----------------------------------------------------------------------------
! Subroutine calc_div_B
! Calculates the divergence of B integrated over each mesh cell
! -----------------------------------------------------------------------------
subroutine calc_div_B(h,mu_hat,div_B)
use parameters
use files
implicit none

complex(prc),pointer :: h(:,:,:,:)
complex(prc),pointer :: mu_hat(:,:,:,:,:)
complex(prc),pointer :: div_B(:,:,:)

integer :: ix,iy,iz,ixp1,iyp1,izp1,i

! We actually calculate (div B)/Q0 here

div_B=0.0

do iz=pmlzi+1,pmlzf-1
   do iy=1,iymax
      do ix=1,ixmax

! Boundary conditions for mu. Assumes periodic but this will also
! work for metal BC's as the field will be zero anyway.

         ixp1=ix+1
         iyp1=iy+1
         izp1=iz+1
         if (ix==ixmax) ixp1=1
         if (iy==iymax) iyp1=1
         if (iz==izmax) izp1=1

! Calculate the divergence

         div_B(ix,iy,iz)=(0.0,0.0)
         do i=1,3
            div_B(ix,iy,iz)=div_B(ix,iy,iz)                        &
&                          +mu_hat(1,i,ixp1,iy,iz)*h(i,ix+1,iy,iz) &
&                          -mu_hat(1,i,ix,iy,iz)*h(i,ix,iy,iz)     &
&                          +mu_hat(2,i,ix,iyp1,iz)*h(i,ix,iy+1,iz) &
&                          -mu_hat(2,i,ix,iy,iz)*h(i,ix,iy,iz)     &
&                          +mu_hat(3,i,ix,iy,izp1)*h(i,ix,iy,iz+1) &
&                          -mu_hat(3,i,ix,iy,iz)*h(i,ix,iy,iz)
         enddo
      enddo
   enddo
enddo

return
end subroutine calc_div_B


! -----------------------------------------------------------------------------
! Subroutine calc_div_D
! Calculates the divergence of D integrated over each mesh cell
! -----------------------------------------------------------------------------
subroutine calc_div_D(e,div_D)
use parameters
use files
implicit none

complex(prc),pointer :: e(:,:,:,:)
complex(prc),pointer :: div_D(:,:,:)

integer :: ix,iy,iz,ixm1,iym1,izm1,i

! We actually calculate (div D)/Q0 here

div_D=0.0

!do iz=pmlzi+1,pmlzf-1
!   do iy=pmlyi+1,pmlyf-1
!      do ix=pmlxi+1,pmlxf-1

do iz=pmlzi+1,pmlzf-1
   do iy=1,iymax
      do ix=1,ixmax

! Boundary conditions for eps. Assumes periodic but this will also
! work for metal BC's as the field will be zero anyway.

         ixm1=ix-1
         iym1=iy-1
         izm1=iz-1
         if (ix==1) ixm1=ixmax
         if (iy==1) iym1=iymax
         if (iz==1) izm1=izmax

! Calculate the divergence

         div_D(ix,iy,iz)=(0.0,0.0)         
		    !>>> Fórmula que no contiene eps solo válida para el caso de vacío
            div_D(ix,iy,iz)=div_D(ix,iy,iz)                         &
&                          +real(e(1,ix,iy,iz))     &
&                          -real(e(1,ix-1,iy,iz)) &
&                          +real(e(2,ix,iy,iz))     &
&                          -real(e(2,ix,iy-1,iz)) &
&                          +real(e(3,ix,iy,iz))     &
&                          -real(e(3,ix,iy,iz-1))
             !    Fórmula que no contiene eps solo válida para el caso de vacío <<<         
      enddo
   enddo
enddo

return
end subroutine calc_div_D


! -----------------------------------------------------------------------------
! Subroutine calc_energy_density
! Calculates the stored energy in the system excluding the absorbing layers
! -----------------------------------------------------------------------------
subroutine calc_energy_density(e_cur,e_prev,h_prev,rho,eps_inv,zi,zf,bbox)
use physconsts
use parameters
use files
implicit none

complex(prc),pointer :: e_cur(:,:,:,:)
complex(prc),pointer :: e_prev(:,:,:,:)
complex(prc),pointer :: h_prev(:,:,:,:)
complex(prc),pointer                 :: eps_inv(:,:,:,:,:)

real(prc),pointer :: rho(:,:,:)
integer,intent(in)::zi
integer,intent(in)::zf
real(prc),pointer ::bbox(:,:,:)

integer :: ix,iy,iz
integer :: i
real(precision) :: dtcq2

! Calculate the energy stored (rho) in the mesh cell at ix,iy,iz
! Actually calc rho/((Q0*eps0)/(2*Q1*Q2*Q3*|u1.u2xu3|))
! but the constant doesn't matter

dtcq2=(dt*c0/Q)**2

rho(:,:,:)=(0.0,0.0)

do iz=zi,hi-1
 do iy=1,iymax
  do ix=1,ixmax
   do i=1,3   
    rho(ix,iy,iz)=rho(ix,iy,iz)+dble(1/eps_inv(i,i,ix,iy,iz))*(conjg(E_cur(i,ix,iy,iz))*E_prev(i,ix,iy,iz))
    rho(ix,iy,iz)=rho(ix,iy,iz)+dble(1/mu_inv)*(conjg(H_prev(i,ix,iy,iz))*H_prev(i,ix,iy,iz))/dtcq2   
   enddo
  enddo
 enddo
enddo

do iz=hi,hf
 do iy=1,iymax
  do ix=1,ixmax
   do i=1,3   
    rho(ix,iy,iz)=rho(ix,iy,iz)+(1.0 - bbox(ix,iy,iz))*dble(1/eps_inv(i,i,ix,iy,iz))*(conjg(E_cur(i,ix,iy,iz))*E_prev(i,ix,iy,iz))
    rho(ix,iy,iz)=rho(ix,iy,iz)+(1.0 - bbox(ix,iy,iz))*dble(1/mu_inv)*(conjg(H_prev(i,ix,iy,iz))*H_prev(i,ix,iy,iz))/dtcq2   
   enddo
  enddo
 enddo
enddo

do iz=hf+1,zf
 do iy=1,iymax
  do ix=1,ixmax
   do i=1,3   
    rho(ix,iy,iz)=rho(ix,iy,iz)+dble(1/eps_inv(i,i,ix,iy,iz))*(conjg(E_cur(i,ix,iy,iz))*E_prev(i,ix,iy,iz))
    rho(ix,iy,iz)=rho(ix,iy,iz)+dble(1/mu_inv)*(conjg(H_prev(i,ix,iy,iz))*H_prev(i,ix,iy,iz))/dtcq2   
   enddo
  enddo
 enddo
enddo

return
end subroutine calc_energy_density

! -----------------------------------------------------------------------------
! Subroutine calc_current
! Calculates the three components of the current flowing away from the point
! at ix,iy,iz
! -----------------------------------------------------------------------------
subroutine calc_current(J,e_cur,h_cur,ix,iy,iz)
use parameters
implicit none

complex(prc),intent(out) :: J(3)
complex(prc),pointer :: e_cur(:,:,:,:),h_cur(:,:,:,:)
complex(prc),pointer :: h_prev(:,:,:,:)
integer,intent(in) :: ix,iy,iz

!For complex-valued Plane Waves

J(1)=e_cur(2,ix,iy,iz)*conjg(h_cur(3,ix,iy,iz)) - e_cur(3,ix,iy,iz)*conjg(h_cur(2,ix,iy,iz))

J(2)=e_cur(3,ix,iy,iz)*conjg(h_cur(1,ix,iy,iz)) - e_cur(1,ix,iy,iz)*conjg(h_cur(3,ix,iy,iz))

J(3)=e_cur(1,ix,iy,iz)*conjg(h_cur(2,ix,iy,iz)) - e_cur(2,ix,iy,iz)*conjg(h_cur(1,ix,iy,iz))

return
end subroutine calc_current

!******************************************************************!

subroutine calc_current_OLD(J,e_cur,h_cur,h_prev,ix,iy,iz)
use parameters
implicit none

complex(prc),intent(out) :: J(3)
complex(prc),pointer :: e_cur(:,:,:,:),h_cur(:,:,:,:)
complex(prc),pointer :: h_prev(:,:,:,:)
integer,intent(in) :: ix,iy,iz

J(1)=-(e_cur(3,ix+1,iy,iz)*conjg(h_cur(2,ix,iy,iz)) &
&   -e_cur(2,ix+1,iy,iz)*conjg(h_cur(3,ix,iy,iz))   &
&   +conjg(e_cur(3,ix+1,iy,iz))*h_prev(2,ix,iy,iz)  &
&   -conjg(e_cur(2,ix+1,iy,iz))*h_prev(3,ix,iy,iz))

J(2)=-(e_cur(1,ix,iy+1,iz)*conjg(h_cur(3,ix,iy,iz)) &
&   -e_cur(3,ix,iy+1,iz)*conjg(h_cur(1,ix,iy,iz))   &
&   +conjg(e_cur(1,ix,iy+1,iz))*h_prev(3,ix,iy,iz)  &
&   -conjg(e_cur(3,ix,iy+1,iz))*h_prev(1,ix,iy,iz))

J(3)=-(e_cur(2,ix,iy,iz+1)*conjg(h_cur(1,ix,iy,iz)) &
&   -e_cur(1,ix,iy,iz+1)*conjg(h_cur(2,ix,iy,iz))   &
&   +conjg(e_cur(2,ix,iy,iz+1))*h_prev(1,ix,iy,iz)  &
&   -conjg(e_cur(1,ix,iy,iz+1))*h_prev(2,ix,iy,iz))

return
end subroutine calc_current_OLD