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

module MESH
interface
!#2
subroutine slit(bbox,eps_inv)
use parameters
	implicit none
	real(prc),pointer :: bbox(:,:,:)
	complex(prc),pointer                    :: eps_inv(:,:,:,:,:)
end subroutine slit

! **************************************** !

subroutine set_emconstants(eps_inv)
use parameters
	implicit none   
	complex(prc),pointer                    :: eps_inv(:,:,:,:,:)
end subroutine set_emconstants

subroutine mesh_graph(bbox,eps_inv)
	use parameters
	use files
	use SDL_MESSAGES
    use SDL_MEMBERS 
    implicit none
    real(prc),pointer :: bbox(:,:,:)	          	  
    complex(prc),pointer                    :: eps_inv(:,:,:,:,:)	
end subroutine mesh_graph

end interface
end module MESH

! **************************************** !
subroutine slit(bbox,eps_inv)
use parameters
use physconsts
implicit none
real(prc),pointer :: bbox(:,:,:)
complex(prc),pointer                    :: eps_inv(:,:,:,:,:)

bbox =0.0 ! bbox could be used to define metals (metal=1; otherwise=0)
eps_inv(:,:,aix:afx,:,hi:hf)=1.0d0 ! A slit of air in the dielectric
!Modify here eps_inv to build your own lossless structure between z=[hi,hf]

return
end subroutine slit

! **************************************** !
    
! **************************************** !

subroutine set_emconstants(eps_inv)
use parameters
use physconsts 
use files
	implicit none   
    complex(prc),pointer                    :: eps_inv(:,:,:,:,:)
	integer :: i,ix,iy,iz
	
     n(1) = sqrt(eps(1))
     n(2) = sqrt(eps(3))

	 !>>> eps anisotropic regions I and III
	 epsa(1,1,1)         = eps(1);epsa(1,1,2)         = 0.0   ;epsa(1,1,3)         = 0.0
	 epsa(1,2,1)         = 0.0   ;epsa(1,2,2)         = eps(1);epsa(1,2,3)         = 0.0
	 epsa(1,3,1)         = 0.0	 ;epsa(1,3,2)         = 0.0   ;epsa(1,3,3)         = eps(1)
	 epsa(2,1,1)         = eps(3);epsa(2,1,2)         = 0.0   ;epsa(2,1,3)         = 0.0
	 epsa(2,2,1)         = 0.0   ;epsa(2,2,2)         = eps(3);epsa(2,2,3)         = 0.0
	 epsa(2,3,1)         = 0.0	 ;epsa(2,3,2)         = 0.0   ;epsa(2,3,3)         = eps(3)
	 !>>> eps_inv anisotropic regions I and III
	 epsa_inv(1,1,1)         = 1.0d0/eps(1);epsa_inv(1,1,2)         = 0.0         ;epsa_inv(1,1,3)         = 0.0
	 epsa_inv(1,2,1)         = 0.0         ;epsa_inv(1,2,2)         = 1.0d0/eps(1);epsa_inv(1,2,3)         = 0.0
	 epsa_inv(1,3,1)         = 0.0	       ;epsa_inv(1,3,2)         = 0.0         ;epsa_inv(1,3,3)         = 1.0d0/eps(1)
	 epsa_inv(2,1,1)         = 1.0d0/eps(3);epsa_inv(2,1,2)         = 0.0         ;epsa_inv(2,1,3)         = 0.0
	 epsa_inv(2,2,1)         = 0.0         ;epsa_inv(2,2,2)         = 1.0d0/eps(3);epsa_inv(2,2,3)         = 0.0
	 epsa_inv(2,3,1)         = 0.0	       ;epsa_inv(2,3,2)         = 0.0         ;epsa_inv(2,3,3)         = 1.0d0/eps(3)
	
	
     allocate (eps_inv(1,1,0:ixmax+1,0:iymax+1,hi-1:hf+1)) !>>> hi-1;hf+1 porque en subroutinas fdtd.f90 hacen referencia a estos valores
                                                             !>>> lo mismo para 0:ixmax+1 0:iymax+1 (añadido 29nov16 para hyperbolic)
	 eps_inv               = 1.0/eps_bbox
	 eps_inv(:,:,:,:,hi-1) = 1.0/eps(1)
	 eps_inv(:,:,:,:,hf+1) = 1.0/eps(3)            	
	
	return
end subroutine set_emconstants

! **************************************** !

subroutine mesh_graph(bbox,eps_inv)
use parameters
use files
use SDL_MESSAGES
use SDL_MEMBERS 
implicit none
  real(prc),pointer :: bbox(:,:,:)	  
  complex(prc),pointer                    :: eps_inv(:,:,:,:,:)  
  complex(prc),pointer::sdl(:,:,:,:)
  real(prc),pointer::maxem(:,:)
  integer::ilayer
  integer i,j,k,r
  integer :: graphic,nolayer,il    
  logical :: flag
  integer       :: mxi2,mxf2,myi2,myf2,mzi2,mzf2 !To solve BUG: SERS input override by mesh_graph
  
  msg='' 
  rewind(unit=inFDTD)  
  do while(msg/=msg2)
     read(inFDTD,*) msg
  end do
  
read(inFDTD,*) graphic  
if(graphic==1)then
  open(unit=meshgraph,file='mesh_graph.dat',status='replace',action='write')
  read(inFDTD,*) LayerMESH
  mxi2=mxi;mxf2=mxf;myi2=myi;myf2=myf;mzi2=mzi;mzf2=mzi!To solve BUG: SERS input override by mesh_graph
  read(inFDTD,*) mxi,mxf,myi,myf,mzi,mzf
  read(inFDTD,*) nolayer    

if(.true.)then
  if(nolayer/=0)then   
   il = 1
   do while(il<=nolayer)
     read(inFDTD,*) k
     write(meshgraph,*) '#',k			  		  
	   select case (LayerMESH)
	     case (1) ! >>> Plano X=cte <<< 
		   do i=myi,myf
			 write(meshgraph,*) '#',i
             do j=mzi,mzf			   
		       write(meshgraph,*) (1.0-bbox(k,i,j))/real(eps_inv(1,1,k,i,j))			   
			 end do 
           end do
	      case (2) ! >>> Plano Y=cte <<<	
		   !do i=1,ixmax		    
		     do i=mxi,mxf
			 write(meshgraph,*) '#',i
             do j=mzi,mzf
		       write(meshgraph,*) (1.0-bbox(i,k,j))/real(eps_inv(1,1,i,k,j))			   			   
			 end do 
           end do       
	      case (3) ! >>> Plano Z=cte <<<
		   do i=mxi,mxf
			 write(meshgraph,*) '#',i
             do j=myi,myf
		       write(meshgraph,*) (1.0-bbox(i,j,k))/real(eps_inv(1,1,i,j,k))			   			   
			 end do 
           end do     
	   end select	                  	   			
     il = il+1
	end do  
  end if
write(logfile,*) "...mesh-graph generated"
print*,          "...mesh-graph generated"

else

do i=100,100
 do j=1,iymax
   flag = .true.
   do k=hi,hf
   if(flag)then
     if(bbox(i,j,k)==1.0)then
		write(meshgraph,*) i," ",j," ",hf-k
		flag = .false.	
     end if
   end if
   end do
 end do
end do

end if
	  
close(unit=meshgraph)
mxi=mxi2;mxf=mxf2;myi=myi2;myf=myf2;mzi=mzi2;mzf=mzi2!To solve BUG: SERS input override by mesh_graph
stop
end if
return  
         	  
end subroutine mesh_graph


