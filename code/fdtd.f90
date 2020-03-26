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

module FDTD_CODE
interface
subroutine e_fdtd(it,e_cur,e_prev,h_prev,eps_inv,bbox)
use parameters
use physconsts
implicit none
integer              :: it
	complex(prc),pointer :: e_cur(:,:,:,:)
	complex(prc),pointer :: e_prev(:,:,:,:)	
	complex(prc),pointer :: h_prev(:,:,:,:)
	
	complex(prc),pointer                    :: eps_inv(:,:,:,:,:)	
	real(prc),pointer ::bbox(:,:,:)			
end subroutine e_fdtd

subroutine h_fdtd(h_cur,h_prev,e_cur,dtcq2)
	use parameters
	use physconsts
	implicit none
    
	complex(prc),pointer :: h_cur (:,:,:,:)
	complex(prc),pointer :: h_prev(:,:,:,:)
	complex(prc),pointer :: e_cur (:,:,:,:)	
     
	real(precision),intent(in) :: dtcq2	
end subroutine h_fdtd

end interface
end module FDTD_CODE
          
!***************************************************************************************************!

subroutine e_fdtd(it,e_cur,e_prev,h_prev,eps_inv,bbox)
use parameters
use physconsts
implicit none
integer              :: it
complex(prc),pointer :: e_cur(:,:,:,:)
complex(prc),pointer :: e_prev(:,:,:,:)
complex(prc),pointer :: h_prev(:,:,:,:)

complex(prc),pointer :: eps_inv(:,:,:,:,:)
real(prc),pointer    :: bbox(:,:,:)

integer :: ix,iy,iz
complex(precision) :: curl(3)

! First the E-fields
do iz=1+pmlzi,hi-1 
  do iy=pmlyi_aux+1,pmlyf_aux-1
    do ix=pmlxi_aux+1,pmlxf_aux-1   	
! Define Curl H

         curl(1)=h_prev(3,ix,iy,iz)-h_prev(3,ix,iy-1,iz) &
&               -h_prev(2,ix,iy,iz)+h_prev(2,ix,iy,iz-1)
         curl(2)=h_prev(1,ix,iy,iz)-h_prev(1,ix,iy,iz-1) &
&               -h_prev(3,ix,iy,iz)+h_prev(3,ix-1,iy,iz)
         curl(3)=h_prev(2,ix,iy,iz)-h_prev(2,ix-1,iy,iz) &
&               -h_prev(1,ix,iy,iz)+h_prev(1,ix,iy-1,iz)

! Integrate fields in time
         e_cur(1,ix,iy,iz)=e_prev(1,ix,iy,iz)+ epsa_inv(1,1,1)*curl(1)+&
		                                     & epsa_inv(1,1,2)*curl(2)+&
											 & epsa_inv(1,1,3)*curl(3)  

         e_cur(2,ix,iy,iz)=e_prev(2,ix,iy,iz)+ epsa_inv(1,2,1)*curl(1)+&
		                                     & epsa_inv(1,2,2)*curl(2)+&
											 & epsa_inv(1,2,3)*curl(3) 

         e_cur(3,ix,iy,iz)=e_prev(3,ix,iy,iz)+ epsa_inv(1,3,1)*curl(1)+&
		                                     & epsa_inv(1,3,2)*curl(2)+&
											 & epsa_inv(1,3,3)*curl(3)

      enddo
   enddo
enddo
 
! The structure between [hi,hf]
  do iz=hi,hf
   do iy=pmlyi_aux+1,pmlyf_aux-1
     do ix=pmlxi_aux+1,pmlxf_aux-1       

	     	 curl(1)=h_prev(3,ix,iy,iz)-h_prev(3,ix,iy-1,iz) &
	&               -h_prev(2,ix,iy,iz)+h_prev(2,ix,iy,iz-1)
			 curl(2)=h_prev(1,ix,iy,iz)-h_prev(1,ix,iy,iz-1) &
	&               -h_prev(3,ix,iy,iz)+h_prev(3,ix-1,iy,iz)
			 curl(3)=h_prev(2,ix,iy,iz)-h_prev(2,ix-1,iy,iz) &
	&               -h_prev(1,ix,iy,iz)+h_prev(1,ix,iy-1,iz)
               
             			 
	 e_cur(1,ix,iy,iz)=e_prev(1,ix,iy,iz) + eps_inv(1,1,ix,iy,iz)*curl(1)
	 e_cur(2,ix,iy,iz)=e_prev(2,ix,iy,iz) + eps_inv(1,1,ix,iy,iz)*curl(2)
	 e_cur(3,ix,iy,iz)=e_prev(3,ix,iy,iz) + eps_inv(1,1,ix,iy,iz)*curl(3)  

     enddo
   enddo
  enddo 	 	
! The structure between [hi,hf]

do iz=hf+1,pmlzf-1 
  do iy=pmlyi_aux+1,pmlyf_aux-1
    do ix=pmlxi_aux+1,pmlxf_aux-1    
! Define Curl H

         curl(1)=h_prev(3,ix,iy,iz)-h_prev(3,ix,iy-1,iz) &
&               -h_prev(2,ix,iy,iz)+h_prev(2,ix,iy,iz-1)
         curl(2)=h_prev(1,ix,iy,iz)-h_prev(1,ix,iy,iz-1) &
&               -h_prev(3,ix,iy,iz)+h_prev(3,ix-1,iy,iz)
         curl(3)=h_prev(2,ix,iy,iz)-h_prev(2,ix-1,iy,iz) &
&               -h_prev(1,ix,iy,iz)+h_prev(1,ix,iy-1,iz)

! Integrate fields in time
         e_cur(1,ix,iy,iz)=e_prev(1,ix,iy,iz)+ epsa_inv(2,1,1)*curl(1)+&
		                                     & epsa_inv(2,1,2)*curl(2)+&
											 & epsa_inv(2,1,3)*curl(3)  

         e_cur(2,ix,iy,iz)=e_prev(2,ix,iy,iz)+ epsa_inv(2,2,1)*curl(1)+&
		                                     & epsa_inv(2,2,2)*curl(2)+&
											 & epsa_inv(2,2,3)*curl(3) 

         e_cur(3,ix,iy,iz)=e_prev(3,ix,iy,iz)+ epsa_inv(2,3,1)*curl(1)+&
		                                     & epsa_inv(2,3,2)*curl(2)+&
											 & epsa_inv(2,3,3)*curl(3)

      enddo
   enddo
enddo

end subroutine e_fdtd

!***************************************************************************************************!

subroutine h_fdtd(h_cur,h_prev,e_cur,dtcq2)
use parameters
use physconsts
implicit none

complex(prc),pointer :: h_cur(:,:,:,:)
complex(prc),pointer :: h_prev(:,:,:,:)
complex(prc),pointer :: e_cur(:,:,:,:)

real(precision),intent(in) :: dtcq2
integer :: zmx_eff
integer :: ix,iy,iz
complex(precision) :: curl(3)


if(pmlzf<hf) then
 zmx_eff = hf
else
 zmx_eff = pmlzf-1
end if

do iz=pmlzi,hi-1
	do iy=pmlyi_aux+1,pmlyf_aux-1
      do ix=pmlxi_aux+1,pmlxf_aux-1 
			! Define Curl E
					 curl(1)=e_cur(3,ix,iy+1,iz)-e_cur(3,ix,iy,iz) &
			&               -e_cur(2,ix,iy,iz+1)+e_cur(2,ix,iy,iz)
					 curl(2)=e_cur(1,ix,iy,iz+1)-e_cur(1,ix,iy,iz) &
			&               -e_cur(3,ix+1,iy,iz)+e_cur(3,ix,iy,iz)
					 curl(3)=e_cur(2,ix+1,iy,iz)-e_cur(2,ix,iy,iz) &
			&               -e_cur(1,ix,iy+1,iz)+e_cur(1,ix,iy,iz)

			! Integrate the fields in time

			 h_cur(1,ix,iy,iz)=h_prev(1,ix,iy,iz)-dtcq2*curl(1)                  

			 h_cur(2,ix,iy,iz)=h_prev(2,ix,iy,iz)-dtcq2*curl(2)     

			 h_cur(3,ix,iy,iz)=h_prev(3,ix,iy,iz)-dtcq2*curl(3)

  	   enddo
     enddo
   enddo 


  do iz=hi,hf
	do iy=pmlyi_aux+1,pmlyf_aux-1
      do ix=pmlxi_aux+1,pmlxf_aux-1 
			! Define Curl E
					 curl(1)=e_cur(3,ix,iy+1,iz)-e_cur(3,ix,iy,iz) &
			&               -e_cur(2,ix,iy,iz+1)+e_cur(2,ix,iy,iz)
					 curl(2)=e_cur(1,ix,iy,iz+1)-e_cur(1,ix,iy,iz) &
			&               -e_cur(3,ix+1,iy,iz)+e_cur(3,ix,iy,iz)
					 curl(3)=e_cur(2,ix+1,iy,iz)-e_cur(2,ix,iy,iz) &
			&               -e_cur(1,ix,iy+1,iz)+e_cur(1,ix,iy,iz)

			! Integrate the fields in time

			 h_cur(1,ix,iy,iz)=h_prev(1,ix,iy,iz)-dtcq2*mu_inv_bbox*curl(1)                  

			 h_cur(2,ix,iy,iz)=h_prev(2,ix,iy,iz)-dtcq2*mu_inv_bbox*curl(2)     

			 h_cur(3,ix,iy,iz)=h_prev(3,ix,iy,iz)-dtcq2*mu_inv_bbox*curl(3)

  	   enddo
     enddo
   enddo 

do iz=hf+1,zmx_eff
	do iy=pmlyi_aux+1,pmlyf_aux-1
      do ix=pmlxi_aux+1,pmlxf_aux-1 
			! Define Curl E
					 curl(1)=e_cur(3,ix,iy+1,iz)-e_cur(3,ix,iy,iz) &
			&               -e_cur(2,ix,iy,iz+1)+e_cur(2,ix,iy,iz)
					 curl(2)=e_cur(1,ix,iy,iz+1)-e_cur(1,ix,iy,iz) &
			&               -e_cur(3,ix+1,iy,iz)+e_cur(3,ix,iy,iz)
					 curl(3)=e_cur(2,ix+1,iy,iz)-e_cur(2,ix,iy,iz) &
			&               -e_cur(1,ix,iy+1,iz)+e_cur(1,ix,iy,iz)

			! Integrate the fields in time

			 h_cur(1,ix,iy,iz)=h_prev(1,ix,iy,iz)-dtcq2*curl(1)                  

			 h_cur(2,ix,iy,iz)=h_prev(2,ix,iy,iz)-dtcq2*curl(2)     

			 h_cur(3,ix,iy,iz)=h_prev(3,ix,iy,iz)-dtcq2*curl(3)

  	   enddo
     enddo
   enddo 
return
end subroutine h_fdtd
















