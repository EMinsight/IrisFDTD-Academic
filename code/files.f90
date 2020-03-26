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

module files
! Set the unit numbers for the files we want to access
integer,parameter :: logfile=10,wpfile=11
integer,parameter :: inFDTD=16
integer,parameter :: meshgraph=17
integer,parameter :: out1=18,out2=19,out3=20,out4=21,out5=22,out6=23,out7=24,out8=25,out9=41
integer,parameter :: aux1=26,aux2=27,testPML=28,spectra=29
logical::exists
end module files

! -----------------------------------------------------------------------------
! Subroutine setparam
! Set the key parameters from a config file
! -----------------------------------------------------------------------------
subroutine setparam()
use parameters
use physconsts
use files
implicit none
integer::i,j
character(100) grating,runname


rewind(unit=inFDTD)
print*, msg0
do while(msg/=msg0)
     read(inFDTD,*) msg
end do

!>> Lattice parameters
read(inFDTD,*) i_pol ! ipol=1 -> E//y ; ipol=2 -> E//x

if ((i_pol/=1).and.(i_pol/=2).and.(i_pol/=3).and.(i_pol/=4)) call err(2)
select case(i_pol)
  case(1) !>>> E//y <<<!
     write(logfile,'(A1,A20)') '#','Polarization = E//y'   
  case(2) !>>> E//x <<<!
     write(logfile,'(A1,A20)') '#','Polarization = E//x' 
end select 

read(inFDTD,*) Q
if (Q<=0.0) call err(2)
write(logfile,'(A1,A12,F15.5,A12,F15.5)') '#','q(nm)=',Q,'q(u.a)=',Q*f_abohr

read(inFDTD,*) ixmax
if (ixmax<1) call err(2)
write(logfile,'(A1,A12,I15,A12,F15.5)') '#','ixmax=',ixmax,'Lx(nm)=',ixmax*Q

read(inFDTD,*) iymax
if (iymax<1) call err(2)
write(logfile,'(A1,A12,I15,A12,F15.5)') '#','iymax=',iymax,'Ly(nm)=',iymax*Q

read(inFDTD,*) izmax
if (izmax<1) call err(2)
write(logfile,'(A1,A12,I15,A12,F15.5)') '#','izmax=',izmax,'Lz(nm)=',izmax*Q

Q = f_abohr*Q
dL    = ixmax*Q !>>> Lo más típico <<<

read(inFDTD,*) CFL
if (CFL<=0.0) call err(2)
write(logfile,'(A1,A12,F15.5)') '#','CFL =',CFL

dt = CFL/(c0*dsqrt(3.0/Q**2))
if (dt<=0.0) call err(2)

read(inFDTD,*) wzPML
if (wzPML<0.0) call err(2)
wzPML = f_abohr*(wzPML/Q)
write(logfile,'(A1,A15,F15.5)') '#','PML strength =',wzPML
write(logfile,'(A65)') "...BE CAREFUL WITH wzPML. Now, is aplied wzPML=wzPML(q=1nm)*Q"

read(inFDTD,*) pml_taflove
!read(inFDTD,*) nPML
nPML=2 !Note 21/nov/18: This value has not been changed for long so the reading line for it is reused for pml_taflove
if (nPML<0.0) call err(2)
write(logfile,'(A1,A15,I10)') '#','PML exponent =',nPML
    pmlxi = 0
	pmlxf = ixmax+1
	pmlyi = 0
	pmlyf = iymax+1
	pmlxi_aux = pmlxi
	pmlxf_aux = pmlxf
	pmlyi_aux = pmlyi
	pmlyf_aux = pmlyf

	read(inFDTD,*) pmlzi
	write(logfile,'(A1,A12,I10)') '#','pmlzi=',pmlzi

	read(inFDTD,*) pmlzf
	write(logfile,'(A1,A12,I10)') '#','pmlzf=',pmlzf

!>>> Reading optional parameters
rewind(unit=inFDTD)  
do while(msg/=msg1)
     read(inFDTD,*) msg
end do

read(inFDTD,*) runmode,SDL_ON
select case(runmode)  
  case(2)
     runname = runmode_2  
end select
write(logfile,'(A16,A20)') runname//' mode running...'    

read(inFDTD,*) m_structure 
   select case(m_structure)	  
	   case(nostr_2)
	    grating = str_2
	end select    	

read(inFDTD,*) no_downloads,OUTPUT_TIME_STEP

read(inFDTD,*) texit_fs
texit=int(texit_fs/(dt*uat))
itmax=texit
write(logfile,'(A1,A12,F15.8)') '#','dt(fs)=',dt*uat
write(logfile,'(A1,A26,I10,A26,F15.5,A7)') '#','Iteration max. number =',texit,' Simulation time =',uat*dt*texit,' (fs)'
no_downloads=int(texit/no_downloads)

read(inFDTD,*) hi
    write(logfile,'(A1,A12,I10)') '#','hi=',hi

if(SDL_ON==1) read(inFDTD,*) SDL_plane_on(1),SDL_plane_on(2),SDL_plane_on(3)
rewind(unit=inFDTD)  
 do while(msg/=msg5)
   read(inFDTD,*) msg
 end do
	 	
read(inFDTD,*) eps(1)
write(logfile,'(A4,A20,F10.5)') '#  -','Eps región I  =',eps(1)
eps_in(1) = 1.0/eps(1)
read(inFDTD,*) eps_bbox
write(logfile,'(A4,A20,F10.5,F10.5)') '   -','Eps región II =',eps_bbox
eps_in(2) = 1.0/eps_bbox
eps(2)=eps_bbox !This asignation is important in a few places (see drude.f90 drude_lorentz_hyperbolic()
read(inFDTD,*) eps(3)
write(logfile,'(A4,A20,F10.5)') '   -','Eps región III=',eps(3)
eps_in(3) = 1.0/eps(3)


rewind(unit=inFDTD)
do while(msg/=runmode_2)
 if(not(eof(inFDTD))) then 
 read(inFDTD,*) msg
 else
 exit
 end if
end do

if(runmode==norunmode_2)then
	read(inFDTD,*) rwmode	
	write(logfile,'(A1,A22,I1)') '#','F(r,wo) or F(ro,w)..=',rwmode
        read(inFDTD,*) sxi,mxi
        read(inFDTD,*) sxf,mxf
        read(inFDTD,*) syi,myi
        read(inFDTD,*) syf,myf
        read(inFDTD,*) szi,mzi
        read(inFDTD,*) szf,mzf  
    
		rewind(unit=inFDTD)
        do while(msg/='w=cte')
          read(inFDTD,*) msg
        end do		
        
		read(inFDTD,*) emmode
		write(logfile,'(A1,A22,I1)') '#','Electromagnectic mode=',emmode		

		read(inFDTD,*) file_direct
		write(logfile,'(A1,A22,I1)') '#','Direct wavelength reading=',file_direct 

		read(inFDTD,*) nof	
	
		select case (file_direct)
		case(0)
		 nowav=1
		 allocate(fw(nof))
		 i=1
		 read(inFDTD,*) dlamda
		 read(inFDTD,*) lno
         lno_aux=lno
         do while(i<=nof)		   
		   fw(i) =2.0*pi*c0*dt/(f_abohr*lno_aux)
		   write(wpfile,*),i,' = ',lno_aux
		   i= i+1		
		   lno_aux = lno_aux + dlamda              
         end do	 
           lnf= lno_aux
           nwo=1; nwf=nof; nopto=nof !for spectra_lossless_dielectric() in common.f90
		case(1)	
		 nowav=1	
		 allocate(fw(nof))
		 i=1
		 do while(i<=nof)
		   read(inFDTD,*) lno		  
		   fw(i) =2.0*pi*c0*dt/(f_abohr*lno)
		   print*,i,") ",lno,'-',fw(i) 		   		  	                
		   i= i+1		   
		 end do	 
		 case(2)		  
		  read(inFDTD,*) nowav !>>> number of points allocated at each wavelength (must be even)
		  read(inFDTD,*) dlamda		  
		  allocate(fw(nof*nowav))
		  i=1
		  do while(i<=nof)
		   read(inFDTD,*) lno		  
		   do j=1,nowav
		    fw((i-1)*nowav+j) =2.0*pi*c0*dt/(f_abohr*(lno-(nowav/2-j)*dlamda))
			write(78,*) lno-(nowav/2-j)*dlamda
		   end do
		   print*,i,") ",lno,'-',fw(i) 		   		  	                
		   i= i+1		   
		  end do	 		 
		  nof = nowav*nof !>>> This is important when running the subroutines remaining.
		 end select


   read(inFDTD,*) to_Fourier,sample_Fourier
   write(logfile,'(A1,A12,(I2))') '#','to_Fourier=',to_Fourier  	
	
end if


rewind(unit=inFDTD)  
do while(msg/=msg12)
if(not(eof(inFDTD))) then
     read(inFDTD,*) msg
else
	exit
end if
end do 

if(msg==msg12)then
em_source(0)=1 !TF/SF by default in this version 
if(em_source(0)/=0)then        
    RATIO_ANGLE= 1.0 !Normal incidence by default for this version
    read(inFDTD,*) plane_wave
    read(inFDTD,*) tincr_PW,tdecay_PW
    write(logfile,'(A1,A22,I1)') '#','Plane-Wave=',plane_wave
    read(inFDTD,*) Epump,wlpump,pump_DT,def_pulse       !>>> kV/cm,nm,fs,decay factor def_pulse=E^(-def_pulse)	
    read(inFDTD,*) tfsf_zi !It is redefined in most of the cases
    tfsf_zf=tfsf_zi+5 !An infinite half-space starting at z=tfsf_o by default in this version
    curr_zo=pmlzi+1;curr_zf=pmlzf-1 
	read(inFDTD,*) cosine    
    write(logfile,'(A1,A40,F15.5)') '#','TF/SF angle of incidence(x_direcction)=',ratio_angle*(pi/2.0)
           
    tincr_PW= tincr_PW/uat
    tdecay_PW=tdecay_PW/uat
	Epump=Epump/(f_e*1e-5) !>>> from kV/cm to u.a. electric field
	Eoxy=Epump
	pump_DT=pump_DT/uat
	Gpump=c0*pump_DT/(2.d0*q*sqrt(eps(1)*def_pulse)) !>>> in mesh units
	zopump=-sqrt(def_pulse)*Gpump-tfsf_zi

    write(logfile,*) Gpump,zopump,(c0*pump_DT/(2.d0*q*gap*sqrt(eps(1))))**2
end if
end if

!>>> Extra parameters for SDL: several layers in the three planes
SDL_Xextralayers=1; SDL_Yextralayers=1;SDL_Zextralayers=1 !By default
allocate(SDL_deltaX(1:SDL_Xextralayers),SDL_deltaY(1:SDL_Yextralayers),SDL_deltaZ(1:SDL_Zextralayers))
SDL_deltaX=0;SDL_deltaY=0;SDL_deltaZ=0; !By default

	rewind(unit=inFDTD)
	msgaux=trim(msg7)//'-SDLextended'

	 do while(msg/=msgaux)	 
	   if(not(eof(inFDTD))) then
		 read(inFDTD,*) msg   	
	   else   
		 exit
	   end if
	 end do
   if(msg==msgaux)then      
   read(inFDTD,*) SDL_Xextralayers,SDL_Yextralayers,SDL_Zextralayers !They may be 1 at least   
   deallocate(SDL_deltaX,SDL_deltaY,SDL_deltaZ)
   allocate(SDL_deltaX(1:SDL_Xextralayers),SDL_deltaY(1:SDL_Yextralayers),SDL_deltaZ(1:SDL_Zextralayers))
   do i=1,SDL_Xextralayers
   read(inFDTD,*) SDL_deltaX(i) !Regarding xi(1) for X=cte layer
   end do
   do i=1,SDL_Yextralayers
   read(inFDTD,*) SDL_deltaY(i)!Regarding yi(2) for Y=cte layer
   end do
   do i=1,SDL_Zextralayers
   read(inFDTD,*) SDL_deltaZ(i) !Regarding zi(3) for Z=cte layer
   end do
   end if

! >>> Leyendo parámetros de estructura
  rewind(unit=inFDTD)
  do while(msg/=grating)
    read(inFDTD,*) msg
  end do

  write(logfile,'(A1,A60)') '#',msg  
  read(inFDTD,*) Dh !>>> Structure thickness
  hf = (hi+Dh) - 1  !>>> hf final coordinate of the structure
  write(logfile,'(A1,A12,I10)') '#','System starts at z=hi=',hi
  write(logfile,'(A1,A12,I10)') '#','System ends   at z=hf=',hf
  write(logfile,'(A1,A12,F15.8)') '#','System thickness, h=',(Q/f_abohr)*(hf+1-hi)


  select case(m_structure)	        
	   case(nostr_2)   
	     read(inFDTD,*) dx	
         call struct_limit(ixmax,dx,aix,afx)
   	     write(logfile,'(A1,A12,I10)') '#','aix=',aix		
		 write(logfile,'(A1,A12,I10)') '#','afx=',afx
		 write(logfile,'(A1,A12,F15.8)') '#','    dx=',(Q/f_abohr)*(afx+1-aix)    
   end select  
return
end subroutine setparam

! ****************************************************** !

! -----------------------------------------------------------------------------
! Subroutine err
! Error trap. Write error message to log file and stop
! -----------------------------------------------------------------------------
subroutine err(ierr)
use files
use parameters
implicit none

integer,intent(in) :: ierr

error: select case(ierr)
   case(1)
      write(logfile,*) Err_1
   case(2)
      write(logfile,*) Err_2  
   case(4)
      write(logfile,*) Err_3
   case(5)
      write(logfile,*) Err_4
   case(6)
      write(logfile,*) Err_5
   case(9)
      write(logfile,*) Err_6
   case(10)
      write(logfile,*) Err_7   
   case(13)
      write(logfile,*) Err_8
   case default
      write(logfile,*) Err_0
end select error

stop
end subroutine err

!**************************************************!

subroutine struct_limit(N,a,ai,af)
implicit none
   
   integer,intent(in)::N,a
   integer,intent(out)::ai,af
 
   if(Mod(a,2)==0) then !>>> a par <<<!  
		if(Mod(N,2)==0)then !>>> N par <<<!
            ai = N/2 - (a/2-1)
			af = N/2 + (a/2)			
        else
		    ai = (N+1)/2 - a/2 
			af = (N+1)/2 + a/2 - 1			
		end if
      else
	    if(Mod(N,2)==0)then !>>> N par <<<!
            ai = N/2 - ((a+1)/2-1) 
			af = N/2 + ((a-1)/2)
        else
		    ai = (N+1)/2 - (a-1)/2
			af = (N+1)/2 + (a-1)/2			
		end if
      end if
return
end subroutine struct_limit
