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

module FDTD_DRIVER
 interface
subroutine driver(e,h,eps_inv,bbox,fft_data,fft_data_o,store_pts,fdata)
use boundary
use UPML_2DArray
use parameters
use files
use physconsts
use SDL_MESSAGES
use MESH
use FDTD_CODE
use PROCESO_DATOS
use TFSF
use SDL_MEMBERS
	
implicit none
 complex(prc),pointer                 :: e(:,:,:,:)
 complex(prc),pointer                 :: h(:,:,:,:)
 complex(prc),pointer                 :: eps_inv(:,:,:,:,:)
 real(prc),pointer                    :: bbox(:,:,:)
	
 complex(prc),pointer                 :: fft_data(:,:)	
 complex(prc),pointer                 :: fft_data_o(:,:)
 integer,pointer                      :: store_pts(:,:)
   	 
 type(INTEGER_POINTER_ARRAY),pointer  :: k_wmax(:)
 type(REAL_POINTER_ARRAY),pointer     :: k_w(:) 
 integer,pointer                      :: no_propag(:)    
	
 complex(prc),pointer                 :: fdata(:,:,:,:,:)	 
end subroutine driver

 end interface
end module FDTD_DRIVER


! -----------------------------------------------------------------------------
! Subroutine driver_data
! Heart of the calculation. This subroutine actually does the time integration
! loop. Should be made to be as flexible as possible
! -----------------------------------------------------------------------------
subroutine driver(e_cur,h_cur,eps_inv,bbox,fft_data,fft_data_o,store_pts,fdata)
use BOUNDARY
use UPML_2DArray
use parameters
use files
use physconsts
use SDL_MESSAGES
use MESH
use FDTD_CODE
use PROCESO_DATOS
use TFSF
use SDL_MEMBERS 

implicit none

complex(prc),pointer     :: e_cur(:,:,:,:)
complex(prc),pointer     :: h_cur(:,:,:,:)

complex(prc),pointer     :: temp(:,:,:,:),temp2(:,:,:,:,:)
complex(precision),pointer        :: tempr(:,:,:,:)
complex(prc),pointer     :: pmltemp(:,:,:,:)
complex(prc),pointer     :: ccomtemp(:,:,:,:)

complex(prc),pointer     :: eps_inv(:,:,:,:,:)
real(prc),pointer        :: bbox(:,:,:)

complex(prc),pointer     :: fft_data(:,:)
complex(prc),pointer     :: fft_data_o(:,:)
integer,pointer          :: store_pts(:,:)

complex(prc),pointer     :: fdata(:,:,:,:,:)
complex(prc),pointer     :: fdata_scattering(:,:,:,:,:)
complex(prc),pointer     :: fdata_source(:,:,:,:,:)

complex(prc),pointer    :: e_prev(:,:,:,:)
complex(prc),pointer    :: h_prev(:,:,:,:)
complex(prc),pointer :: daux_prev(:,:,:,:)
complex(prc),pointer :: daux_cur(:,:,:,:) 

complex(prc),pointer :: Pn1(:,:,:,:,:)
complex(prc),pointer :: Pn(:,:,:,:,:) 
complex(prc),pointer :: Pn_1(:,:,:,:,:) 
complex(precision),pointer    :: Npn(:,:,:,:)
complex(precision),pointer    :: Npn1(:,:,:,:)

complex(precision),pointer     :: dipoleJ(:)
complex(precision),pointer     :: FTdipoleJ(:,:)

!>>> POLARIZER
complex(prc)         :: fx,fy
!<<< POLARIZER

! <<<Función auxiliar para metales que debe actualizarse cada instante 
	complex(prc),pointer :: ud(:,:,:,:,:)	
	complex(prc),pointer :: ulr(:,:,:,:,:),uli(:,:,:,:,:)
! Función auxiliar para metales que debe actualizarse cada instante>>>

! <<<SDL
	integer              :: nArchivo
	integer              :: nos
	real(prc),pointer    :: maxem(:,:)
! SDL>>>


!>>> Scattering/Absorbption
  type(COMPLEX_POINTER_ARRAY_VI),pointer  :: emg(:)
  type(COMPLEX_POINTER_ARRAY_VII),pointer :: spts(:)
  type(COMPLEX_POINTER_ARRAY_VI),pointer  :: emg2(:)
  type(COMPLEX_POINTER_ARRAY_VII),pointer :: spts2(:)
!    Scattering/Absorbption <<<
  
! >>> TF/SF
complex(prc),pointer          :: tfsf_data(:,:,:,:), egauss(:,:,:)
! <<< TF/SF


! >>> UMPL
type(COMPLEX_POINTER_ARRAY_II),pointer :: d_cur(:)
type(COMPLEX_POINTER_ARRAY_II),pointer :: b_cur(:)
type(COMPLEX_POINTER_ARRAY_II),pointer :: d_prev(:)
type(COMPLEX_POINTER_ARRAY_II),pointer :: b_prev(:)
real(precision),pointer                :: fpzmn(:,:),fmzmn(:,:),fpzmx(:,:),fmzmx(:,:)
real(precision),pointer                ::fp(:,:,:,:),fm(:,:,:,:)
!     UMPL <<<


integer                  :: ix,iy,iz,it,i_pts
complex(precision)       :: curl(3)
real(precision)          :: dtcq2
real(precision)          :: w
real(prc),pointer        :: rho(:,:,:)
integer                  :: nt
integer                  :: i,ij,ip,id,iw
integer,dimension(8)     :: hora,hora2
CHARACTER(8)                dia
CHARACTER(10)               hora_char
real(precision)          :: ddt

!complex(prc)             :: J(3)
real(precision)          :: Jz
complex(prc),pointer     :: J(:,:)
real(precision)          :: Jy

!>>> TEST ABC´s 
real(prc)::e_ref(itmax/no_downloads)

logical                :: single_box

nt = texit	

! >>> UMPL  
	 allocate(d_cur(1:2))
	 allocate(d_cur(1)%pml(3,ixmax,iymax,1:pmlzi))
	 allocate(d_cur(2)%pml(3,ixmax,iymax,pmlzf:izmax))
	 allocate(b_cur(1:2))
	 allocate(b_cur(1)%pml(3,ixmax,iymax,1:pmlzi))
	 allocate(b_cur(2)%pml(3,ixmax,iymax,pmlzf:izmax))
	 allocate(d_prev(1:2))
	 allocate(d_prev(1)%pml(3,ixmax,iymax,1:pmlzi))
	 allocate(d_prev(2)%pml(3,ixmax,iymax,pmlzf:izmax))
	 allocate(b_prev(1:2))
	 allocate(b_prev(1)%pml(3,ixmax,iymax,1:pmlzi))
	 allocate(b_prev(2)%pml(3,ixmax,iymax,pmlzf:izmax))
	 allocate(fpzmn(1:3,pmlzi),fmzmn(1:3,pmlzi))
	 allocate(fpzmx(1:3,pmlzf:izmax),fmzmx(1:3,pmlzf:izmax))
	 d_cur(1)%pml  = (0.0,0.0)
	 b_cur(1)%pml  = (0.0,0.0)
	 d_prev(1)%pml = (0.0,0.0)
	 b_prev(1)%pml = (0.0,0.0)
	 d_cur(2)%pml  = (0.0,0.0)
	 b_cur(2)%pml  = (0.0,0.0)
	 d_prev(2)%pml = (0.0,0.0)
	 b_prev(2)%pml = (0.0,0.0)

	 call set_sz(fpzmn,fmzmn,fpzmx,fmzmx)
print*,"...UPML done"
	!     UPML <<<
    

if(runmode==norunmode_2)then   
   if(rwmode==6)then
   !>>> BOX one
     allocate(spts(1:3)) !>> i=1,3 : i->Cada par de superficies paralelas, pj. ixmin e ixmax
	 if(ixmax/=1)then
	  allocate(spts(1)%pto(2*6*(syf+1-syi)*(szf+1-szi),4))
	 else
      allocate(spts(1)%pto(1,1))
	 end if  	 
	 if(iymax/=1)then
	  allocate(spts(2)%pto(2*6*(sxf+1-sxi)*(szf+1-szi),4))
	 else
      allocate(spts(2)%pto(1,1))
	 end if  
	 allocate(spts(3)%pto(2*6*(sxf+1-sxi)*(syf+1-syi),4))
	 do i=1,3
	    spts(i)%pto = 1
     end do 
	 
	 call init_store_pts_fourier(spts,1)

     allocate(emg(1:3)) !>> i=1,3 : i->Cada par de superficies paralelas, pj. ixmin e ixmax	 
	 
	 if(ixmax/=1) then
	  allocate(emg(1)%fdata(2*6*(syf+1-syi)*(szf+1-szi),nof))
	 else
	  allocate(emg(1)%fdata(1,1))
	 end if

	 if(iymax/=1) then
	  allocate(emg(2)%fdata(2*6*(sxf+1-sxi)*(szf+1-szi),nof))
	 else
	  allocate(emg(2)%fdata(1,1))
	 end if
	 allocate(emg(3)%fdata(2*6*(sxf+1-sxi)*(syf+1-syi),nof))
	 do i=1,3
	    emg(i)%fdata = 0.0
     end do    
   !<<< BOX one

   !>>> BOX two
    single_box=(sxi==mxi).AND.(syi==myi).AND.(szi==mzi).AND.(sxf==mxf).AND.(syf==myf).AND.(szf==mzf) 
    if(not(single_box))then
     allocate(spts2(1:3)) !>> i=1,3 : i->Cada par de superficies paralelas, pj. ixmin e ixmax
	 if(ixmax/=1)then
	  allocate(spts2(1)%pto(2*6*(myf+1-myi)*(mzf+1-mzi),4))
	 else
      allocate(spts2(1)%pto(1,1))
	 end if  	 
	 if(iymax/=1)then
	  allocate(spts2(2)%pto(2*6*(mxf+1-mxi)*(mzf+1-mzi),4))
	 else
      allocate(spts2(2)%pto(1,1))
	 end if  
	 allocate(spts2(3)%pto(2*6*(mxf+1-mxi)*(myf+1-myi),4))
	 do i=1,3
	    spts2(i)%pto = 1
     end do 
	 
	 call init_store_pts_fourier(spts2,2)

     allocate(emg2(1:3)) !>> i=1,3 : i->Cada par de superficies paralelas, pj. ixmin e ixmax	 
	 
	 if(ixmax/=1) then
	  allocate(emg2(1)%fdata(2*6*(myf+1-myi)*(mzf+1-mzi),nof))
	 else
	  allocate(emg2(1)%fdata(1,1))
	 end if

	 if(iymax/=1) then
	  allocate(emg2(2)%fdata(2*6*(mxf+1-mxi)*(mzf+1-mzi),nof))
	 else
	  allocate(emg2(2)%fdata(1,1))
	 end if
	 allocate(emg2(3)%fdata(2*6*(mxf+1-mxi)*(myf+1-myi),nof))
	 do i=1,3
	    emg2(i)%fdata = 0.0
     end do    
   end if
   !<<< BOX two
   end if     
end if


if(runmode/=norunmode_1)then
if(SDL_ON) then !>>> SDL on Fourier Maps     
   call setparam_SDL(inFDTD,logfile,runmode_1)      
   allocate(maxem(1:2,texit/no_downloads))   
end if
end if

print*,"...run mode options done"

!>>> Initialize 2nd set of fields
select case(ABCmodel)
 case (0)
   allocate(e_prev(3,0:ixmax+1,0:iymax+1,0:izmax+1))
   allocate(h_prev(3,0:ixmax+1,0:iymax+1,0:izmax+1))
   e_prev=(0.0,0.0)
   h_prev=(0.0,0.0) 
end select
print*,"...e_prev,h_prev done"
!    Initialize 2nd set of fields <<<  

    ! >>> TF/SF initialization
	if(em_source(0)/=0)then
    
	     !>>> Para que funcione en Incidencia normal como si fuera pulso gaussiano el rango debe
		 !>>> ser entre 0 e Max+1	
	     iy_min=0
		 iy_max=iymax+1
	     ix_min=0
		 ix_max=ixmax+1	
		 if(tfsf_zi/=0)then
	     iz_min= tfsf_zi
		 iz_max= tfsf_zf
		 end if
		 pmlxi_aux = 0
		 pmlxf_aux = ixmax+1
		 pmlyi_aux = 0
		 pmlyf_aux = iymax+1
	     allocate(tfsf_data(1:6,ix_min-1:ix_max+1,iy_min-1:iy_max+1,iz_min-1:iz_max+1))		
         if(plane_wave==4) then
            allocate(egauss(ix_min-1:ix_max+1,iy_min-1:iy_max+1,iz_min-1:iz_min+1))        
            egauss=0.0
         end if    
	end if	
print*,"...TF/SF done"
! >>> TF/SF initialization
    
dtcq2=(dt*c0/Q)**2
nos = 1
! >>> FDTD core <<< !
time: do it=1,nt 
!>>> Information on the screen
call DATE_AND_TIME(VALUES = hora) 	  
call DATE_AND_TIME(DATE=dia,TIME = hora_char) 
if(MOD(it,OUTPUT_TIME_STEP).EQ.0) then
ddt=hora(5)*60*60*1000+hora(6)*60*1000+hora(7)*1000+hora(8)
ddt=ddt - (hora2(5)*60*60*1000+hora2(6)*60*1000+hora2(7)*1000+hora2(8))
ddt=ddt/1000 !Tiempo en segundos
        
print*, "  Date = ",dia(7:8),"/",dia(5:6),"/",dia(1:4),"-",hora_char(1:2),":",hora_char(3:4),":",hora_char(4:5)
print*, "  Iteration ",it,' of ',nt
print*, "  Iteration/second = ",OUTPUT_TIME_STEP/ddt 
print*, "  Finish in ",((nt-it)*(ddt/OUTPUT_TIME_STEP))/3600,' hours'
print*, "  Completed ",100*dble(it)/dble(nt),' %'
call DATE_AND_TIME(VALUES = hora2)             
end if	!call DATE_AND_TIME
!<<< Information on the screen 

!>>> UPML update
    !>>> Current fields become previous fields
	temp=>e_cur
	e_cur=>e_prev
	e_prev=>temp

	temp=>h_cur
	h_cur=>h_prev
	h_prev=>temp 	
	
	temp=>daux_cur
	daux_cur=>daux_prev
	daux_prev=>temp 	 
    
	pmltemp=>d_cur(1)%pml    
	d_cur(1)%pml=>d_prev(1)%pml
	d_prev(1)%pml=>pmltemp

	pmltemp=>d_cur(2)%pml    
	d_cur(2)%pml=>d_prev(2)%pml
	d_prev(2)%pml=>pmltemp

	pmltemp=>b_cur(1)%pml
	b_cur(1)%pml=>b_prev(1)%pml
	b_prev(1)%pml=>pmltemp

	pmltemp=>b_cur(2)%pml
	b_cur(2)%pml=>b_prev(2)%pml
	b_prev(2)%pml=>pmltemp 
	!    Current fields become previous fields <<<

	!>>> "Black box" update - E field
    call e_fdtd(it,e_cur,e_prev,h_prev,eps_inv,bbox)
    !    "Black box" update - E field <<<

	! >>> BOUNDARY CONDITIONS for the E-field
    !Update Berenger type PML equations for the E-field (Absorbing boundary layer)	
	 call UPML_Array_e(e_cur,e_prev,h_prev,d_cur,d_prev,fpzmn,fmzmn,fpzmx,fmzmx)  	  				 
	! En este punto se actualiza el campo en el contorno de la derecha ya que
	! se necesitará en el cálculo de H
	! Update the boundary conditions needed for the E-fields
    
	 call bc_xmax_bloch(e_cur,h_cur)
	 call bc_ymax_bloch(e_cur,h_cur)		 		
	 call bc_zmax_bloch(e_cur,h_cur)	
	! <<< BOUNDARY CONDITIONS for the E-field
     
	
	!>>> "Black box" update - H field
	call h_fdtd(h_cur,h_prev,e_cur,dtcq2) 
    !    "Black box" update - H field <<<
        
	! >>> BOUNDARY CONDITIONS for the H-field
	! Update Berenger type PML equations for the H-field (Absorbing boundary layer)	  
	 call UPML_Array_h(e_cur,h_cur,h_prev,b_cur,b_prev,fpzmn,fmzmn,fpzmx,fmzmx)		
	! Update the boundary conditions needed for the H-fields    
	 
     call bc_xmin_bloch(e_cur,h_cur)
	 call bc_ymin_bloch(e_cur,h_cur)	
	 call bc_zmin_bloch(e_cur,h_cur) !>>> Boundary conditions through 3D	
    ! <<< BOUNDARY CONDITIONS for the H-field    		 				


!>>> Total Field Scattering Field (TF/SF) source conditions
if(em_source(0)/=0)then
    tfsf_data=0.0
    select case (plane_wave) 
	case(0)	  	  	 
	  call tfsf_plane_wave(it,Epump,wlpump,tincr_PW,tdecay_PW,tfsf_data) 		  
    case(1)      
	  call tfsf_gauss(it,zopump,Epump,wlpump,Gpump,tfsf_data)       
	end select
end if

!>>>  ¿Por qué tanta subroutina? Verás:
!>>>  i) Las .._2D permiten cálculos en incidencia obliqua
!>>> ii) Las .._2D no funcionan en sistemas 3D
!>>>iii) Las .._3D no permiten cálculos en incidencia obliqua
!>>>iv) Las .._3D sirven en onda s y p en 2D (aunque no funcionan para onda s en 3D
!>>>     , de todas formas, en 3D puedes "rotar la estructura", de forma que puedes hacer
!>>>     cálculos para distintas polarizaciones
 		if((iymax/=1).and.(ixmax/=1)) then
			call tfsf_update_3DArray(it,e_cur,h_cur,tfsf_data) !>>> .._3D doesn't work with s-polarization in 3D<<<!
		else		
			 if(i_pol==2)   call tfsf_update_2DArray(it,e_cur,h_cur,tfsf_data) !>>> .._2D doesn't work with s-polarization <<<!
			 if(i_pol==1)   call tfsf_update_3DArray(it,e_cur,h_cur,tfsf_data)				
        end if  
!<<< Total Field Scattering Field (TF/SF) source conditions
 
! >>> Fourier Transformation on-the-fly 
if((it>=to_Fourier).and.(MOD(it,sample_Fourier)==0)) call scattering_fourier(emg,spts,e_cur,h_cur,it+1-to_Fourier,1)		   
single_box=(sxi==mxi).AND.(syi==myi).AND.(szi==mzi).AND.(sxf==mxf).AND.(syf==myf).AND.(szf==mzf) 
if(not(single_box))then
    if((it>=to_Fourier).and.(MOD(it,sample_Fourier)==0)) call scattering_fourier(emg2,spts2,e_cur,h_cur,it+1-to_Fourier,2)	
end if
! <<< Fourier Transformation on-the-fly 

! >>> On-the-fly data proccesing 
! >>> Results are saved in files no_downloads times
if(it>=to_Fourier)then
    if ((it==(to_fourier+(no_downloads*nos))).or.(it==nt)) then	  	             
          if(runmode==norunmode_2)then
		    select case(rwmode)		
			    case(6)
                call scattering_fourier_out(emg,spts,it,1)
                single_box=(sxi==mxi).AND.(syi==myi).AND.(szi==mzi).AND.(sxf==mxf).AND.(syf==myf).AND.(szf==mzf) 
		        if(not(single_box))then
			      call scattering_fourier_out(emg2,spts2,it,2)
                end if         						  
            end select			 
          end if	             
	       
          if(SDL_ON)then                !
            if(SDL_plane_on(1)==1) then
            do i=1,SDL_Xextralayers
               SDL_xi1=xi(1)
               xi(1)=xi(1)+SDL_deltaX(i)
		       call setSDL_FILE(e_cur,h_prev,maxem,nos+1,1)
               xi(1)=SDL_xi1
            end do
            end if               
		   
           if(SDL_plane_on(2)==1) then
            do i=1,SDL_Yextralayers
               SDL_yi2=yi(2)
               yi(2)=yi(2)+SDL_deltaY(i)
		       call setSDL_FILE(e_cur,h_prev,maxem,nos+1,2)
               yi(2)=SDL_yi2
            end do
           end if		   
                   
           if(SDL_plane_on(3)==1) then
            do i=1,SDL_Zextralayers
               SDL_zi3=zi(3)
               zi(3)=zi(3)+SDL_deltaZ(i)
		       call setSDL_FILE(e_cur,h_prev,maxem,nos+1,3)
               zi(3)=SDL_zi3
            end do
           end if
        end if           
    nos=nos+1	
    end if		              	  
end if  	  
! <<< On-the-fly data proccesing 
end do time

! >>> Free computer memory
deallocate(e_prev,h_prev,STAT=stat_alloc)
print*,"e_prev,h_prev deallocated...status=",stat_alloc

if((metalmodel==1).or.(metalmodel==4).or.(metalmodel==5).or.(metalmodel==6).or.(metalmodel==7))then
 deallocate(bbox,ud,ulr,uli,STAT=stat_alloc)
 print*,"bbox,ud,ulr,uli deallocated...status=",stat_alloc
end if

select case (runmode) 
 case(norunmode_2)
   select case(rwmode)
    case(6)
	   do i=1,3	   
	    deallocate(emg(i)%fdata,STAT=stat_alloc)
	    deallocate(spts(i)%pto,STAT=stat_alloc)
        single_box=(sxi==mxi).AND.(syi==myi).AND.(szi==mzi).AND.(sxf==mxf).AND.(syf==myf).AND.(szf==mzf) 
        if(not(single_box))then
		deallocate(emg2(i)%fdata,STAT=stat_alloc)
        deallocate(spts2(i)%pto,STAT=stat_alloc)
        end if
	   end do
	   deallocate(emg,spts,STAT=stat_alloc)
	   print*,"emg,spts deallocated...status=",stat_alloc  
   end select 
end select

if(em_source(0)/=0)then
  deallocate(tfsf_data,STAT=stat_alloc)
  print*,"  deallocated...status=",stat_alloc
end if

select case(ABCmodel)
 case(0)
 ! >>> UMPL
     do i = 1,2		
	  deallocate(d_cur(i)%pml,STAT=stat_alloc)
	  deallocate(b_cur(i)%pml,STAT=stat_alloc)
	  deallocate(d_prev(i)%pml,STAT=stat_alloc)
	  deallocate(b_prev(i)%pml,STAT=stat_alloc)	
	 end do

	 deallocate(d_cur,b_cur,d_prev,b_prev,STAT=stat_alloc)
 	 print*,"d_cur,b_cur,d_prev,b_prev deallocated...status=",stat_alloc

	 deallocate(fpzmn,fmzmn,fpzmx,fmzmx,STAT=stat_alloc)
	 print*,"fpzmn,fmzmn,fpzmx,fmzmx deallocated...status=",stat_alloc
 !     UMPL <<< 
 end select 
 ! <<< Free computer memory
close(unit=logfile)

return
end subroutine driver





