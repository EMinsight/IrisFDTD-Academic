!____________________________________________________________________________________________________________________________________________________________________________________________ 
!
! IrisFDTD Academic is an implementation of the Finite-Difference Time-Domain (FDTD) method, one of the most widely used Maxwell’s equation solvers in the field of computational electrodynamics for 
! photonics and nanotechnology.
!
! IrisFDTD Academic is licensed under the AGPL and it is free to use. However, if you are using, or plan to use it, specially if it is for research or academic purposes, please send an email with your 
! name, institution and a brief description of your interest for this program. If you use IrisFDTD in a work that leads to a scientific or academic publication, we would appreciate it if you 
! would kindly cite IrisFDTD in your manuscript as:
!
! > Sergio G. Rodrigo, Optical Properties of Nanostructured Metallic Systems: Studied with the Finite-Difference Time-Domain Method, Springer-Verlag, Berlin, (2012).
!
! Commercial applications may also acquire a commercial license. Please contact Sergio G Rodrigo sergut@unizar.es
!
! Copyright (C) 2005-2020 Sergio G Rodrigo sergut@unizar.es
!
! IrisFDTD Academic is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of
! the License,or (at your option) any later version.
! 
! IrisFDTD Academic is distributed in the hope that it will be useful for research or/and academic purpouses, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
! FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details. You should have received a copy of the GNU Affero General Public License along with IrisFDTD. If not,see 
! http://www.gnu.org/licenses/.
!______________________________________________________________________________________________________________________________________________________________________________________________ 
! 

module TFSF
interface

subroutine tfsf_update_2DArray(it,e_cur,h_cur,tfsf_data)
	use physconsts
	use parameters
	implicit none	
	complex(prc),pointer :: e_cur(:,:,:,:)
	complex(prc),pointer :: h_cur(:,:,:,:)
	complex(prc),pointer :: tfsf_data(:,:,:,:)
	integer :: it
end subroutine

subroutine tfsf_update_3DArray(it,e_cur,h_cur,tfsf_data)
	use physconsts
	use parameters
	implicit none	
	complex(prc),pointer :: e_cur(:,:,:,:)
	complex(prc),pointer :: h_cur(:,:,:,:)
	complex(prc),pointer :: tfsf_data(:,:,:,:)
	integer :: it
end subroutine

subroutine tfsf_gauss(it,izo,Eamp,wlength_o,spectral_width,tfsf_data)
	use parameters
	use physconsts
	use files
	implicit none
	integer,intent(in)::it
	real(precision) :: Eamp,izo
	real(precision) :: wlength_o
	real(precision) :: spectral_width
	complex(prc),pointer :: tfsf_data(:,:,:,:)
end subroutine

subroutine tfsf_plane_wave(it,Eamp,wlength_o,it_o,it_f,tfsf_data)
	use parameters
	use physconsts
	use files
	implicit none
	real(precision)       :: Eamp,wlength_o,it_o,it_f
	integer,intent(in)    :: it
	complex(prc),pointer       :: tfsf_data(:,:,:,:)
end subroutine

end interface
end module


! ****************************************************************** !

subroutine tfsf_update_3DArray(it,e_cur,h_cur,tfsf_data)
use physconsts
use parameters
implicit none

complex(prc),pointer :: e_cur(:,:,:,:)
complex(prc),pointer :: h_cur(:,:,:,:)
complex(prc),pointer :: tfsf_data(:,:,:,:)
integer :: it

integer :: ix,iy,iz
real(precision) :: dtcq2
logical,parameter::rigth=.true.

dtcq2=(dt*c0/Q)**2

!>>> Normal Incidence
 !>>> E->Z=Zmin ; H->Z=Zmin-1
  !>>> E-fields
do iy=iy_min,iy_max
do ix=ix_min,ix_max-1
    e_cur(1,ix,iy,iz_min)=e_cur(1,ix,iy,iz_min) + eps_in(1)*tfsf_data(5,ix,iy,iz_min-1)
end do 
end do

do iy=iy_min,iy_max-1
 do ix=ix_min,ix_max
    e_cur(2,ix,iy,iz_min)=e_cur(2,ix,iy,iz_min) - eps_in(1)*tfsf_data(4,ix,iy,iz_min-1)     
 end do  
end do
  !    E-fields <<<
  !>>> H-fields
do iy=iy_min,iy_max-1
 do ix=ix_min,ix_max
    h_cur(1,ix,iy,iz_min-1)=h_cur(1,ix,iy,iz_min-1) - dtcq2*mu_inv*tfsf_data(2,ix,iy,iz_min)
 end do  
end do

do iy=iy_min,iy_max
do ix=ix_min,ix_max-1
    h_cur(2,ix,iy,iz_min-1)=h_cur(2,ix,iy,iz_min-1) + dtcq2*mu_inv*tfsf_data(1,ix,iy,iz_min) 
end do 
end do
  !    H-fields <<<
 !    E->Z=Zmin ; H->Z=Zmin-1 <<<

if(.false.)then
!>>> E->Z=Zmax ; H->Z=Zmax
  !>>> E-fields
do iy=iy_min,iy_max
do ix=ix_min,ix_max-1
    e_cur(1,ix,iy,iz_max)=e_cur(1,ix,iy,iz_max) - eps_in(1)*tfsf_data(5,ix,iy,iz_max)
end do 
end do

do iy=iy_min,iy_max-1
 do ix=ix_min,ix_max
    e_cur(2,ix,iy,iz_max)=e_cur(2,ix,iy,iz_max) + eps_in(1)*tfsf_data(4,ix,iy,iz_max)     
 end do  
end do
  !    E-fields <<<
  !>>> H-fields
do iy=iy_min,iy_max-1
 do ix=ix_min,ix_max
    h_cur(1,ix,iy,iz_max)=h_cur(1,ix,iy,iz_max) + dtcq2*mu_inv*tfsf_data(2,ix,iy,iz_max)
 end do  
end do

do iy=iy_min,iy_max
do ix=ix_min,ix_max-1
    h_cur(2,ix,iy,iz_max)=h_cur(2,ix,iy,iz_max) - dtcq2*mu_inv*tfsf_data(1,ix,iy,iz_max) 
end do 
end do
  !    H-fields <<< 
 !    E->Z=Zmax ; H->Z=Zmax <<< 
end if
return
end subroutine 


!*******************************************************************************************!




!********************************************************************************************!

subroutine tfsf_update_2DArray(it,e_cur,h_cur,tfsf_data)
use physconsts
use parameters
implicit none

complex(prc),pointer :: e_cur(:,:,:,:)
complex(prc),pointer :: h_cur(:,:,:,:)
complex(prc),pointer :: tfsf_data(:,:,:,:)
integer :: it

integer :: ix,iy,iz
real(precision) :: dtcq2,eps_medium_source
logical,parameter::rigth=.true.

dtcq2=(dt*c0/Q)**2
if(iz_min>hi)then
eps_medium_source=eps(3)
else
eps_medium_source=eps(1)
end if
do iy=iy_min,iy_max
 do ix=ix_min,ix_max

if(.true.)then     
!>>> Left TF/SF boundary
! E-fields
    e_cur(1,ix,iy,iz_min)=e_cur(1,ix,iy,iz_min) + tfsf_data(5,ix,iy,iz_min)/eps_medium_source
    e_cur(2,ix,iy,iz_min)=e_cur(2,ix,iy,iz_min) - tfsf_data(4,ix,iy,iz_min-1)/eps_medium_source
    e_cur(3,ix,iy,iz_min)=e_cur(3,ix,iy,iz_min) 

! H-fields
! s-polarization
    h_cur(1,ix,iy,iz_min-1)=h_cur(1,ix,iy,iz_min-1)-dtcq2*mu_inv*tfsf_data(2,ix,iy,iz_min)
    h_cur(2,ix,iy,iz_min-1)=h_cur(2,ix,iy,iz_min-1)
    h_cur(3,ix,iy,iz_min-1)=h_cur(3,ix,iy,iz_min-1)
! p-polarization
    h_cur(1,ix,iy,iz_min)=h_cur(1,ix,iy,iz_min)
    h_cur(2,ix,iy,iz_min)=h_cur(2,ix,iy,iz_min)+dtcq2*mu_inv*tfsf_data(1,ix,iy,iz_min)
    h_cur(3,ix,iy,iz_min)=h_cur(3,ix,iy,iz_min)
!    Left TF/SF boundary <<<
end if
if(.false.)then
!>>> Right TF/SF boundary
! E-fields
! s-polarization
    e_cur(1,ix,iy,iz_max)=e_cur(1,ix,iy,iz_max)
    e_cur(2,ix,iy,iz_max)=e_cur(2,ix,iy,iz_max)+tfsf_data(4,ix,iy,iz_max)/eps_medium_source
    e_cur(3,ix,iy,iz_max)=e_cur(3,ix,iy,iz_max)
! p-polarization
    e_cur(1,ix,iy,iz_max+1)=e_cur(1,ix,iy,iz_max+1)-tfsf_data(5,ix,iy,iz_max)/eps_medium_source
    e_cur(2,ix,iy,iz_max+1)=e_cur(2,ix,iy,iz_max+1)
    e_cur(3,ix,iy,iz_max+1)=e_cur(3,ix,iy,iz_max+1)

! H-fields
    h_cur(1,ix,iy,iz_max)=h_cur(1,ix,iy,iz_max)+dtcq2*mu_inv*tfsf_data(2,ix,iy,iz_max)
    h_cur(2,ix,iy,iz_max)=h_cur(2,ix,iy,iz_max)-dtcq2*mu_inv*tfsf_data(1,ix,iy,iz_max+1)
    h_cur(3,ix,iy,iz_max)=h_cur(3,ix,iy,iz_max)
!    Right TF/SF boundary <<<
end if
 enddo 
enddo

!    Top TF/SF boundary <<<
return
end subroutine 

! ****************************************************************** !

subroutine tfsf_gauss(it,izo,Eamp,wlength_o,spectral_width,tfsf_data)
use parameters
use physconsts
use files


implicit none
integer,intent(in)::it
real(precision) :: Eamp,izo
real(precision) :: wlength_o
real(precision) :: spectral_width
complex(prc),pointer :: tfsf_data(:,:,:,:)
integer :: ix,iy,iz

real(precision)::ko	 
real(precision)::Dk  
real(precision) :: kzz(1:2)   
complex(precision) :: ex,ey,ez,hx,hy,hz
real(precision) :: e,h,win
real(precision) :: eps_medium_source,sign_kz
integer:: ito,ixo
type K_DIRECTION
     real(precision) x
	 real(precision) y
	 real(precision) z
end type K_DIRECTION   

real(precision) :: DELTA 

type(K_DIRECTION) u

real(precision)::sum

ixo = pmlzi-50
!izo = pmlzi-(3*GAP)
ito =  0

Dk =  spectral_width*Q !GAP*Q                        
ko = 2*Pi/(wlength_o*f_abohr)  !2*Pi/(lf*f_abohr)  
if(iz_min>hi)then
eps_medium_source=eps(3)
sign_kz=-1.0d0
else
eps_medium_source=eps(1)
sign_kz=1.0d0
end if

!tfsf_data = 0.0

DELTA = RATIO_ANGLE * Pi/2.0  !>> RATIO_ANGLE = 1.0 - Incidencia Normal <<<!  	   

u%x = dcos(DELTA)
u%y = 0.0
u%z = dsin(DELTA) !Normal incidence

FTwindow = 0
select case(FTwindow)
 case(0)
   win = 1.0
 case(1)
   win = 0.42 - 0.5*cos(2*pi*it/(texit-1))+0.08*cos(4*pi*it/(texit-1))   !>>> Blackman window <<<!
 case(2)
   win = exp(-(1.0/2.0)*((it-(texit-1)/2.0)/(0.5*(texit-1)/2.0))**2)     !>>> Gauss window <<<!
 case(3)
   win = 0.53836 - 0.46164*cos(2*pi*it/(texit-1))                        !>>> Hamming window  <<<!
 case(4)
   win = 0.5*(1-cos(2*pi*it/(texit-1)))                                  !>>> Hann window  <<<!
end select


select case(i_pol)
 case(1) !>>> E//y <<<!  
   do iy=iy_min-1,iy_max+1
	do ix=ix_min-1,ix_max+1
	 do iz=iz_min-1,iz_max+1
	  kzz(1) = ko*(sqrt(eps_medium_source)*sign_kz*(Q*(ix-ixo)*u%x+Q*(iz-izo)*u%z) - c0*dt*(it-ito))  	  

	  h = EXP(-(((Q*(ix-ixo)*u%x+Q*(iz-izo)*u%z) - c0*dt*(it-ito)/sqrt(eps_medium_source))/Dk)**2) ! (z-vt)^2 v=c/n
	  hx =-sign_kz*(1.0,0.0)*EXP(ci*kzz(1))*h*u%z 
	  hz =(1.0,0.0)*EXP(ci*kzz(1))*h*u%x
	  ey =(1.0,0.0)*(u%x*hz-u%z*sign_kz*hx) !sqrt(eps_medium_source)*

	  if(cosine==1)then
	  ey = (ey+conjg(ey))/2.0
	  hz = (hz+conjg(hz))/2.0
	  hx = (hx+conjg(hx))/2.0
	  end if
      
	  tfsf_data(4,ix,iy,iz)= tfsf_data(4,ix,iy,iz)+Eamp*Q*(dt*c0/Q)*hx
	  tfsf_data(6,ix,iy,iz)=tfsf_data(6,ix,iy,iz)+Eamp*Q*(dt*c0/Q)*hz
	  tfsf_data(2,ix,iy,iz)=tfsf_data(2,ix,iy,iz)+Eamp*Q*ey	  
	  
	 end do
	end do
	end do

 case(2) !>>> E//x <<<!
   do iy=iy_min-1,iy_max+1
	do ix=ix_min-1,ix_max+1
	 do iz=iz_min-1,iz_max+1
	  kzz(1) = ko*(sqrt(eps_medium_source)*sign_kz*(Q*(ix-ixo)*u%x+Q*(iz-izo)*u%z) - c0*dt*(it-ito))  	  

	  e = EXP(-(((Q*(ix-ixo)*u%x+Q*(iz-izo)*u%z) - c0*dt*(it-ito)/sqrt(eps_medium_source))/Dk)**2) 
	  ex = sign_kz*(1.0,0.0)*EXP(ci*kzz(1))*e*u%z 
	  ez =-(1.0,0.0)*EXP(ci*kzz(1))*e*u%x
	  hy =(1.0,0.0)*(u%z*sign_kz*ex-u%x*ez) !sqrt(eps_medium_source)*
	  
	  if(cosine==1)then
	  ex =(ex+conjg(ex))/2.0
	  ez =(ez+conjg(ez))/2.0
	  hy =(hy+conjg(hy))/2.0
	  end if
      
      
	  tfsf_data(1,ix,iy,iz)=tfsf_data(1,ix,iy,iz)+Eamp*Q*ex
	  tfsf_data(3,ix,iy,iz)=tfsf_data(3,ix,iy,iz)+Eamp*Q*ez
	  tfsf_data(5,ix,iy,iz)=tfsf_data(5,ix,iy,iz)+Eamp*Q*(dt*c0/Q)*hy	  
	 end do
	end do
	end do
end select

tfsf_data = win*tfsf_data
return
end subroutine

!*********************************************************************************!

subroutine tfsf_plane_wave(it,Eamp,wlength_o,it_o,it_f,tfsf_data)
use parameters
use physconsts
use files


implicit none
real(precision)          :: Eamp,wlength_o,it_o,it_f
integer,intent(in)::it
complex(prc),pointer :: tfsf_data(:,:,:,:)
integer :: ix,iy,iz

real(precision)::ko	 
real(precision)::Dk  
real(precision) :: kzz
complex(precision) :: ex,ey,ez,hx,hy,hz
real(precision) :: e,h,win
integer:: ito,izo,ixo
type K_DIRECTION
     real(precision) x
	 real(precision) y
	 real(precision) z
end type K_DIRECTION   

real(precision) :: DELTA 

type(K_DIRECTION) u

real(precision)::sum


ixo = pmlzi-50
if(ratio_angle>=0)then
izo = pmlzi-(3*GAP)
else
izo = izmax+(3*GAP)
end if
ito =  0

Dk =  GAP*Q                        
ko = 2*Pi/(wlength_o*f_abohr)  

tfsf_data = 0.0

DELTA = RATIO_ANGLE * Pi/2.0  !>> RATIO_ANGLE = 1.0 - Incidencia Normal <<<!  	   

u%x = dcos(DELTA)
u%y = 0.0
u%z = dsin(DELTA) !Incidencia Normal

FTwindow = 0
select case(FTwindow)
 case(0)
   win = 1.0
 case(1)
   win = 0.42 - 0.5*cos(2*pi*it/(texit-1))+0.08*cos(4*pi*it/(texit-1))   !>>> Blackman window <<<!
 case(2)
   win = exp(-(1.0/2.0)*((it-(texit-1)/2.0)/(0.5*(texit-1)/2.0))**2)     !>>> Gauss window <<<!
 case(3)
   win = 0.53836 - 0.46164*cos(2*pi*it/(texit-1))                        !>>> Hamming window  <<<!
 case(4)
   win = 0.5*(1-cos(2*pi*it/(texit-1)))                                  !>>> Hann window  <<<!
end select


select case(i_pol)
 case(1) !>>> E//y <<<!  
   do iy=iy_min-1,iy_max+1
	do ix=ix_min-1,ix_max+1
	 do iz=iz_min-1,iz_max+1
	  
	  kzz = ko*((Q*(ix-ixo)*u%x+Q*(iz-izo)*u%z) - c0*dt*(it-ito)/sqrt(eps(1)))  	 
    
	 !>>> switch on/off
	  if(.false.)then
		  if(c0*dt*(it-ito)/sqrt(eps(1))<=Q*abs(iz-izo))then
	  		h = EXP(-(Q*abs(iz-izo) - c0*dt*(it-ito)/sqrt(eps(1)))/Dk)
		  else		   
			h = 1.0		   
		  end if	
      else
		  if(c0*dt*(it-it_o)<=Q*abs(iz-izo))then		
			h = EXP(-(Q*abs(iz-izo)-c0*dt*(it-it_o))/Dk)
		  else
		   if(c0*dt*(it-it_f)>=Q*abs(iz-izo))then
			 h = EXP((Q*abs(iz-izo)-c0*dt*(it-it_f))/Dk)
 		   else
			 h = 1.0
		   end if
		  end if
	  end if
	  !<<< switch on/off

	  hx =-(1.0,0.0)*exp(ci*kzz)*h*u%z 
	  hz =(1.0,0.0)*exp(ci*kzz)*h*u%x
	  ey =sqrt(eps(1))*(1.0,0.0)*(u%x*hz-u%z*hx)

	  if(cosine==1)then
	  ey = (ey+conjg(ey))/2.0
	  hz = (hz+conjg(hz))/2.0
	  hx = (hx+conjg(hx))/2.0
	  end if
      
	  tfsf_data(4,ix,iy,iz)=Eoxy*Q*(dt*c0/Q)*hx
	  tfsf_data(6,ix,iy,iz)=Eoxy*Q*(dt*c0/Q)*hz
	  tfsf_data(2,ix,iy,iz)=Eoxy*Q*ey	  
	 end do
	end do
	end do

 case(2) !>>> E//x <<<!
  do iy=iy_min-1,iy_max+1
	do ix=ix_min-1,ix_max+1
	 do iz=iz_min-1,iz_max+1
	  
	  kzz = ko*((Q*(ix-ixo)*u%x+Q*(iz-izo)*u%z) - c0*dt*(it-ito)/sqrt(eps(1)))
	  
	  !>>> switch on/off
	  if(.false.)then
		  if(c0*dt*(it-ito)/sqrt(eps(1))<=Q*abs(iz-izo))then
	  		e = EXP(-(Q*abs(iz-izo) - c0*dt*(it-ito)/sqrt(eps(1)))/Dk)
		  else		   
			e = 1.0		   
		  end if	
      else
		  if(c0*dt*(it-it_o)<=Q*abs(iz-izo))then		
			e = EXP(-(Q*abs(iz-izo)-c0*dt*(it-it_o))/Dk)
		  else
		   if(c0*dt*(it-it_f)>=Q*abs(iz-izo))then
			 e = EXP((Q*abs(iz-izo)-c0*dt*(it-it_f))/Dk)
 		   else
			 e = 1.0
		   end if
		  end if
	  end if
	  !<<< switch on/off    
	   
  	  ex = (1.0,0.0)*EXP(ci*kzz)*e*u%z 
	  ez =-(1.0,0.0)*EXP(ci*kzz)*e*u%x
	  hy =sqrt(eps(1))*(1.0,0.0)*(u%z*ex-u%x*ez)

	  if(cosine==1)then
	  ex =(ex+conjg(ex))/2.0
	  ez =(ez+conjg(ez))/2.0
	  hy =(hy+conjg(hy))/2.0
	  end if
      
	  tfsf_data(1,ix,iy,iz)=Eoxy*Q*ex
	  tfsf_data(3,ix,iy,iz)=Eoxy*Q*ez
	  tfsf_data(5,ix,iy,iz)=Eoxy*Q*(dt*c0/Q)*hy	  	   	  
	  
	 end do
	end do
	end do
end select

tfsf_data = win*tfsf_data
return
end subroutine

!*********************************************************************************!
