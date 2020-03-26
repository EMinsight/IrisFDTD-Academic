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
    
module MAIN_PROGRAM
use parameters
use FDTD_DRIVER
use PROCESO_DATOS
use MESH
use physconsts
use files
use COMMON_SUBROUTINES
interface
 subroutine main(i,j,k)
  implicit none
  integer,optional::i,j,k
 end subroutine
end interface
end module MAIN_PROGRAM

! ******************************************************************** !

subroutine main(i,j,k)
use parameters
use FDTD_DRIVER
use PROCESO_DATOS
use MESH
use physconsts
use files
use COMMON_SUBROUTINES
use SDL_MESSAGES 
implicit none


integer,optional       :: i,j,k

integer                :: ios

integer,pointer        :: store_pts(:,:)

complex(prc),pointer   :: fdata(:,:,:,:,:)

complex(prc),pointer   :: e(:,:,:,:)
complex(prc),pointer   :: h(:,:,:,:)

complex(prc),pointer   :: eps_inv(:,:,:,:,:)
complex(prc),pointer   :: bbox_aux(:,:,:,:,:)
real(prc),pointer      :: bbox(:,:,:)

complex(prc),pointer   :: fft_data(:,:)
complex(prc),pointer   :: fft_data_o(:,:)
integer                :: iw,ix,iy,iz,indx,indx_o,indx_f

!>>> Variables que se utilizan en proceso de datos en cada región
type(INTEGER_POINTER_ARRAY),pointer   :: k_wmax(:)
type(REAL_POINTER_ARRAY),pointer      :: k_w(:) 
integer,pointer                       :: no_propag(:)
complex(precision),pointer            :: kz_wmax(:)
!    Variables que se utilizan en proceso de datos en cada región <<< 

character               :: no_salir
character(100)           :: descripcion
CHARACTER(10)              hora
CHARACTER(8)               dia
character(50)              grating
integer                 :: dimk
integer                 :: near_field
integer                 :: doyouwanna

character(100)           :: RELEASE_VERSION
character(100)           :: LICENCE

real :: x

real(prc),pointer    :: maxem(:,:)
 

                                       !*******************************!
                                       !*******************************!
        RELEASE_VERSION = "Copyright (C) 2005-2020 Sergio G Rodrigo sergut@unizar.es"
        LICENCE =         "         (GNU Affero General Public License)             "
									   !*******************************!
                                       !*******************************!  


!******************************************** BEGIN PROGRAM ***********************************************************!

if (present(i))then
	descripcion = ''
else
   read(inFDTD,'(A100)') descripcion	   
   !descripcion = ''
end if

open(unit=logfile,file='log.dat',status='replace',action='write')
print*
write(logfile,*) trim(RELEASE_VERSION)
print*,          trim(RELEASE_VERSION)
write(logfile,*) trim(LICENCE)
print*,          trim(LICENCE)
print*
write(logfile,*) "Description     : ",trim(descripcion)
print*,          "Description     : ",trim(descripcion)
print*

! >>> Init parameters
print*,"-> Init parameters"
call setparam()
! >>> Initialise position dependent eps 
 call set_emconstants(eps_inv)
 call spectra_lossless_dielectric()
!     Initialise position dependent eps <<<

!>>> Running initializing  >>> bbox not initialized yet!!
    if(runmode==norunmode_2)then  
    select case(rwmode)
		case(0)
		   select case(emmode)
		   	 Case (0)
		   	   allocate(fdata(6,sxi:sxf,syi:syf,szi:szf,nof)) 			   
		   	 case (1)
		   	   allocate(fdata(1:3,sxi:sxf,syi:syf,szi:szf,nof)) 	
		   	 case (2)
		   	   allocate(fdata(4:6,sxi:sxf,syi:syf,szi:szf,nof)) 	
		   	end select
		   	fdata = 0.0		 	 		 	     	    
	end select		     
    end if
!   Running initializing <<<
     
!>>> Build structure
allocate(bbox(0:ixmax,0:iymax,hi:hf))  
bbox = 0.0 
print*,"...black box done"

select case(m_structure)
   case(nostr_2)		  
   call slit(bbox,eps_inv) 	
end select    	
print*,"...structure done" 		
!    Build structure <<< 

! >>> Allocate arrays for fields  
select case(ABCmodel)
case(0)  !>>> UPML <<<!
  allocate(e(3,0:ixmax+1,0:iymax+1,0:izmax+1))
  allocate(h(3,0:ixmax+1,0:iymax+1,0:izmax+1))
end select
e=(0.0,0.0)
h=(0.0,0.0)
!     Allocate arrays for fields  <<<

!>>> Transfer control to central driver
 call DATE_AND_TIME(DATE=dia,TIME = hora) 
 print*
 write(logfile,*) '...Starting main calculation ',dia(7:8),"/",dia(5:6),"/",dia(1:4),"-",hora(1:2),":",hora(3:4),":",hora(4:5)
 print*, '...Starting main calculation ',dia(7:8),"/",dia(5:6),"/",dia(1:4),"-",hora(1:2),":",hora(3:4),":",hora(4:5)
       
	call driver(e,h,eps_inv,bbox,fft_data,fft_data_o,store_pts,fdata)
		
 call DATE_AND_TIME(DATE=dia,TIME = hora)
 write(logfile,*) '...Main calculation end ',dia(7:8),"/",dia(5:6),"/",dia(1:4),"-",hora(1:2),":",hora(3:4),":",hora(4:5)  
 print*, '...Main calculation end ',dia(7:8),"/",dia(5:6),"/",dia(1:4),"-",hora(1:2),":",hora(3:4),":",hora(4:5)  
 write(logfile,*) '...done! (or not?)'
!   Transfer control to central driver <<<		

!>>> Tidy up
 deallocate(e,h)	 
 print*,"e,h deallocated...status=",stat_alloc

 deallocate(eps_inv)
 print*,"eps_inv deallocated...status=",stat_alloc
 
select case(runmode)
 case(2)    
  if((rwmode==0).or.(rwmode==1)) deallocate(fdata,STAT=stat_alloc)    
  print*,"fdata deallocated...status=",stat_alloc  
  deallocate(fw,                 STAT=stat_alloc)     	    	    	    
  print*,"fw deallocated...status=",stat_alloc
 case(3)  
  !deallocate(fft_data_o,         STAT=stat_alloc)
  !print*,"fft_data_o deallocated...status=",stat_alloc
  deallocate(fft_data,           STAT=stat_alloc)
  print*,"fft_data deallocated...status=",stat_alloc
  deallocate(store_pts,STAT=stat_alloc)    
  print*,"store_pts deallocated...status=",stat_alloc

 case default		  
end select 
! Tidy up <<<
return
!******************************************** END PROGRAM ***********************************************************!
end subroutine main




