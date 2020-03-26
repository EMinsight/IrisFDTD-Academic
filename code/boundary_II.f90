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

module	UPML_2DArray
interface

subroutine UPML_Array_e(e_cur,e_prev,h_prev,d_cur,d_prev,fpzmn,fmzmn,fpzmx,fmzmx)
use parameters
use physconsts
	implicit none
	complex(prc),pointer :: e_cur(:,:,:,:),e_prev(:,:,:,:),h_prev(:,:,:,:)	
	type(COMPLEX_POINTER_ARRAY_II),pointer :: d_cur(:),d_prev(:)        
	real(precision),pointer :: fpzmn(:,:),fmzmn(:,:)
	real(precision),pointer :: fpzmx(:,:),fmzmx(:,:)
end subroutine 

subroutine UPML_Array_h(e_cur,h_cur,h_prev,b_cur,b_prev,fpzmn,fmzmn,fpzmx,fmzmx)
use parameters
use physconsts
	implicit none
	complex(prc),pointer :: e_cur(:,:,:,:),h_cur(:,:,:,:),h_prev(:,:,:,:)	
    type(COMPLEX_POINTER_ARRAY_II),pointer :: b_cur(:),b_prev(:)		
	real(precision),pointer :: fpzmn(:,:),fmzmn(:,:)
	real(precision),pointer :: fpzmx(:,:),fmzmx(:,:)
end subroutine 

subroutine set_sz(fpzmn,fmzmn,fpzmx,fmzmx)
use parameters
use physconsts
	implicit none
	real(precision),pointer :: fpzmn(:,:),fmzmn(:,:)
	real(precision),pointer :: fpzmx(:,:),fmzmx(:,:)
end subroutine set_sz
end interface
end module 

! ************************************************************************** !
                                                          
subroutine UPML_Array_e(e_cur,e_prev,h_prev,d_cur,d_prev,fpzmn,fmzmn,fpzmx,fmzmx)
use parameters
use physconsts
	implicit none
	complex(prc),pointer :: e_cur(:,:,:,:),e_prev(:,:,:,:),h_prev(:,:,:,:)	
	type(COMPLEX_POINTER_ARRAY_II),pointer :: d_cur(:),d_prev(:)	
	real(precision),pointer :: fpzmn(:,:),fmzmn(:,:)
	real(precision),pointer :: fpzmx(:,:),fmzmx(:,:)

integer :: ix,iy,iz
complex(precision) :: curl(3)

! Update the face at iz=1 
do iz=1,pmlzi
  do iy=1,iymax
    do ix=1,ixmax
 ! Define Curl H
		 curl(1)=h_prev(3,ix,iy,iz)-h_prev(3,ix,iy-1,iz) &
&               -h_prev(2,ix,iy,iz)+h_prev(2,ix,iy,iz-1)
		 curl(2)=h_prev(1,ix,iy,iz)-h_prev(1,ix,iy,iz-1) &
&               -h_prev(3,ix,iy,iz)+h_prev(3,ix-1,iy,iz)
		 curl(3)=h_prev(2,ix,iy,iz)-h_prev(2,ix-1,iy,iz) &
&               -h_prev(1,ix,iy,iz)+h_prev(1,ix,iy-1,iz)
		
curl = (eps0**2)*curl

 d_cur(1)%pml(1,ix,iy,iz)=(fmzmn(2,iz)/fpzmn(2,iz))*d_prev(1)%pml(1,ix,iy,iz)+ &
                      &   ((2*eps(1))/fpzmn(2,iz))*curl(1)  
 d_cur(1)%pml(2,ix,iy,iz)=(fmzmn(3,iz)/fpzmn(3,iz))*d_prev(1)%pml(2,ix,iy,iz)+ &
                      &   ((2*eps(1))/fpzmn(3,iz))*curl(2)  
 d_cur(1)%pml(3,ix,iy,iz)=(fmzmn(1,iz)/fpzmn(1,iz))*d_prev(1)%pml(3,ix,iy,iz)+ &
                      &   ((2*eps(1))/fpzmn(1,iz))*curl(3)                           
 

 e_cur(1,ix,iy,iz)=(fmzmn(3,iz)/fpzmn(3,iz))*e_prev(1,ix,iy,iz) + &
                  & (1/(eps0*eps(1)*fpzmn(3,iz)))*(fpzmn(1,iz)*d_cur(1)%pml(1,ix,iy,iz) - fmzmn(1,iz)*d_prev(1)%pml(1,ix,iy,iz))
 e_cur(2,ix,iy,iz)=(fmzmn(1,iz)/fpzmn(1,iz))*e_prev(2,ix,iy,iz) + &
                  & (1/(eps0*eps(1)*fpzmn(1,iz)))*(fpzmn(2,iz)*d_cur(1)%pml(2,ix,iy,iz) - fmzmn(2,iz)*d_prev(1)%pml(2,ix,iy,iz))
 e_cur(3,ix,iy,iz)=(fmzmn(2,iz)/fpzmn(2,iz))*e_prev(3,ix,iy,iz) + &
                  & (1/(eps0*eps(1)*fpzmn(2,iz)))*(fpzmn(3,iz)*d_cur(1)%pml(3,ix,iy,iz) - fmzmn(3,iz)*d_prev(1)%pml(3,ix,iy,iz))


    enddo
  enddo
enddo

! Update the face at izmax
do iz=pmlzf,izmax    
   do iy=1,iymax
    do ix=1,ixmax       
 	     curl(1)=h_prev(3,ix,iy,iz)-h_prev(3,ix,iy-1,iz) &
&               -h_prev(2,ix,iy,iz)+h_prev(2,ix,iy,iz-1)
	     curl(2)=h_prev(1,ix,iy,iz)-h_prev(1,ix,iy,iz-1) &
&               -h_prev(3,ix,iy,iz)+h_prev(3,ix-1,iy,iz)
	     curl(3)=h_prev(2,ix,iy,iz)-h_prev(2,ix-1,iy,iz) &
&               -h_prev(1,ix,iy,iz)+h_prev(1,ix,iy-1,iz)

curl = (eps0**2)*curl

 d_cur(2)%pml(1,ix,iy,iz)=(fmzmx(2,iz)/fpzmx(2,iz))*d_prev(2)%pml(1,ix,iy,iz)+ &
                      &   ((2*eps(3))/fpzmx(2,iz))*curl(1)  
 d_cur(2)%pml(2,ix,iy,iz)=(fmzmx(3,iz)/fpzmx(3,iz))*d_prev(2)%pml(2,ix,iy,iz)+ &
                      &   ((2*eps(3))/fpzmx(3,iz))*curl(2)  
 d_cur(2)%pml(3,ix,iy,iz)=(fmzmx(1,iz)/fpzmx(1,iz))*d_prev(2)%pml(3,ix,iy,iz)+ &
                      &   ((2*eps(3))/fpzmx(1,iz))*curl(3)                           
 

 e_cur(1,ix,iy,iz)=(fmzmx(3,iz)/fpzmx(3,iz))*e_prev(1,ix,iy,iz) + &
                  & (1/(eps0*eps(3)*fpzmx(3,iz)))*(fpzmx(1,iz)*d_cur(2)%pml(1,ix,iy,iz) - fmzmx(1,iz)*d_prev(2)%pml(1,ix,iy,iz))
 e_cur(2,ix,iy,iz)=(fmzmx(1,iz)/fpzmx(1,iz))*e_prev(2,ix,iy,iz) + &
                  & (1/(eps0*eps(3)*fpzmx(1,iz)))*(fpzmx(2,iz)*d_cur(2)%pml(2,ix,iy,iz) - fmzmx(2,iz)*d_prev(2)%pml(2,ix,iy,iz))
 e_cur(3,ix,iy,iz)=(fmzmx(2,iz)/fpzmx(2,iz))*e_prev(3,ix,iy,iz) + &
                  & (1/(eps0*eps(3)*fpzmx(2,iz)))*(fpzmx(3,iz)*d_cur(2)%pml(3,ix,iy,iz) - fmzmx(3,iz)*d_prev(2)%pml(3,ix,iy,iz))
      enddo
    enddo
  enddo

return
end subroutine 

! ************************************************************************** !

subroutine UPML_Array_h(e_cur,h_cur,h_prev,b_cur,b_prev,fpzmn,fmzmn,fpzmx,fmzmx)
use parameters
use physconsts
	implicit none
	complex(prc),pointer :: e_cur(:,:,:,:),h_cur(:,:,:,:),h_prev(:,:,:,:)
	type(COMPLEX_POINTER_ARRAY_II),pointer :: b_cur(:),b_prev(:)	
	real(precision),pointer :: fpzmn(:,:),fmzmn(:,:)
	real(precision),pointer :: fpzmx(:,:),fmzmx(:,:)

integer :: ix,iy,iz
real(precision) :: wz_e,wz_h
real(precision) :: dtq2
complex(precision) :: curl(3)

dtq2=(dt/Q)**2

! Update the face at iz=1
do iz=1,pmlzi
  do iy=1,iymax
    do ix=1,ixmax
! Define Curl E
		 curl(1)=e_cur(3,ix,iy+1,iz)-e_cur(3,ix,iy,iz) &
&               -e_cur(2,ix,iy,iz+1)+e_cur(2,ix,iy,iz)
		 curl(2)=e_cur(1,ix,iy,iz+1)-e_cur(1,ix,iy,iz) &
&               -e_cur(3,ix+1,iy,iz)+e_cur(3,ix,iy,iz)
		 curl(3)=e_cur(2,ix+1,iy,iz)-e_cur(2,ix,iy,iz) &
&               -e_cur(1,ix,iy+1,iz)+e_cur(1,ix,iy,iz)
 
 b_cur(1)%pml(1,ix,iy,iz)=(fmzmn(2,iz)/fpzmn(2,iz))*b_prev(1)%pml(1,ix,iy,iz) - &
                      &   ((2*eps(1))/fpzmn(2,iz))*dtq2*curl(1)  
 b_cur(1)%pml(2,ix,iy,iz)=(fmzmn(3,iz)/fpzmn(3,iz))*b_prev(1)%pml(2,ix,iy,iz) - &
                      &   ((2*eps(1))/fpzmn(3,iz))*dtq2*curl(2)  
 b_cur(1)%pml(3,ix,iy,iz)=(fmzmn(1,iz)/fpzmn(1,iz))*b_prev(1)%pml(3,ix,iy,iz) - &
                      &   ((2*eps(1))/fpzmn(1,iz))*dtq2*curl(3)                           
 

 h_cur(1,ix,iy,iz)=(fmzmn(3,iz)/fpzmn(3,iz))*h_prev(1,ix,iy,iz) + &
                  & (1/(mu0*fpzmn(3,iz)))*(fpzmn(1,iz)*b_cur(1)%pml(1,ix,iy,iz) - fmzmn(1,iz)*b_prev(1)%pml(1,ix,iy,iz))
 h_cur(2,ix,iy,iz)=(fmzmn(1,iz)/fpzmn(1,iz))*h_prev(2,ix,iy,iz) + &
                  & (1/(mu0*fpzmn(1,iz)))*(fpzmn(2,iz)*b_cur(1)%pml(2,ix,iy,iz) - fmzmn(2,iz)*b_prev(1)%pml(2,ix,iy,iz))
 h_cur(3,ix,iy,iz)=(fmzmn(2,iz)/fpzmn(2,iz))*h_prev(3,ix,iy,iz) + &
                  & (1/(mu0*fpzmn(2,iz)))*(fpzmn(3,iz)*b_cur(1)%pml(3,ix,iy,iz) - fmzmn(3,iz)*b_prev(1)%pml(3,ix,iy,iz))

      enddo
    enddo
  enddo


! Update the face at izmax
  do iz=pmlzf,izmax    
   do iy=1,iymax
    do ix=1,ixmax
! Define Curl E
		 curl(1)=e_cur(3,ix,iy+1,iz)-e_cur(3,ix,iy,iz) &
&               -e_cur(2,ix,iy,iz+1)+e_cur(2,ix,iy,iz)
		 curl(2)=e_cur(1,ix,iy,iz+1)-e_cur(1,ix,iy,iz) &
&               -e_cur(3,ix+1,iy,iz)+e_cur(3,ix,iy,iz)
		 curl(3)=e_cur(2,ix+1,iy,iz)-e_cur(2,ix,iy,iz) &
&               -e_cur(1,ix,iy+1,iz)+e_cur(1,ix,iy,iz)

 b_cur(2)%pml(1,ix,iy,iz)=(fmzmx(2,iz)/fpzmx(2,iz))*b_prev(2)%pml(1,ix,iy,iz) - &
                      &   ((2*eps(3))/fpzmx(2,iz))*dtq2*curl(1)  
 b_cur(2)%pml(2,ix,iy,iz)=(fmzmx(3,iz)/fpzmx(3,iz))*b_prev(2)%pml(2,ix,iy,iz) - &
                      &   ((2*eps(3))/fpzmx(3,iz))*dtq2*curl(2)  
 b_cur(2)%pml(3,ix,iy,iz)=(fmzmx(1,iz)/fpzmx(1,iz))*b_prev(2)%pml(3,ix,iy,iz) - &
                      &   ((2*eps(3))/fpzmx(1,iz))*dtq2*curl(3)                           
 

 h_cur(1,ix,iy,iz)=(fmzmx(3,iz)/fpzmx(3,iz))*h_prev(1,ix,iy,iz) + &
                  & (1/(mu0*fpzmx(3,iz)))*(fpzmx(1,iz)*b_cur(2)%pml(1,ix,iy,iz) - fmzmx(1,iz)*b_prev(2)%pml(1,ix,iy,iz))
 h_cur(2,ix,iy,iz)=(fmzmx(1,iz)/fpzmx(1,iz))*h_prev(2,ix,iy,iz) + &
                  & (1/(mu0*fpzmx(1,iz)))*(fpzmx(2,iz)*b_cur(2)%pml(2,ix,iy,iz) - fmzmx(2,iz)*b_prev(2)%pml(2,ix,iy,iz))
 h_cur(3,ix,iy,iz)=(fmzmx(2,iz)/fpzmx(2,iz))*h_prev(3,ix,iy,iz) + &
                  & (1/(mu0*fpzmx(2,iz)))*(fpzmx(3,iz)*b_cur(2)%pml(3,ix,iy,iz) - fmzmx(3,iz)*b_prev(2)%pml(3,ix,iy,iz))
      enddo
    enddo
  enddo
return
end subroutine 

!******************************************************************************************!
                  
subroutine set_sz(fpzmn,fmzmn,fpzmx,fmzmx)
use parameters
use physconsts
use files
implicit none
 real(precision),pointer :: fpzmn(:,:),fmzmn(:,:)
 real(precision),pointer :: fpzmx(:,:),fmzmx(:,:)

 real(precision),pointer :: sg(:,:),k(:,:) 
 real(precision),parameter::pmlalfa = 1.8
 real(precision)::u
 integer :: i

write(logfile,*) '...UPML 2D Array running'

allocate(sg(1:3,pmlzi),k(1:3,pmlzi))

! >>> Z = 0 
sg=0.0
do i=1,pmlzi
  sg(3,i)=sqrt(eps(1))*wzPML*((pmlzi-i)**nPML/(dble(pmlzi-1))**nPML) 
enddo

k=1.0
do i=1,pmlzi
  k(3,i)=1 + sqrt(eps(1))*pmlalfa*((pmlzi-i)**nPML/(dble(pmlzi-1))**nPML)  
enddo

fpzmn = 2*eps0*eps(1)*k + dt*sg
fmzmn = 2*eps0*eps(1)*k - dt*sg
!     Z = 0 <<<

deallocate(sg,k)
allocate(sg(1:3,pmlzf:izmax),k(1:3,pmlzf:izmax))

! >>> Z = Zmax 
sg=0.0
do i=pmlzf,izmax
  sg(3,i)=sqrt(eps(3))*wzPML*((i-pmlzf))**nPML/(dble(izmax - pmlzf)**nPML)  
enddo

k=1.0
do i=pmlzf,izmax
  k(3,i)=1 + sqrt(eps(3))*pmlalfa*((i-pmlzf))**nPML/(dble(izmax - pmlzf)**nPML) 
enddo

fpzmx = 2*eps0*eps(3)*k + dt*sg
fmzmx = 2*eps0*eps(3)*k - dt*sg
!     Z = Zmax <<<

deallocate(sg,k)
return
end subroutine set_sz

