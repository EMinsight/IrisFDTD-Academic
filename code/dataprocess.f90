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

module PROCESO_DATOS
interface 

 subroutine scattering_fourier(emg,spts,e,h,it,nobox)
	use parameters	
	use physconsts
	use COMMON_SUBROUTINES
	use files
	 implicit none	 
	 type(COMPLEX_POINTER_ARRAY_VI),pointer :: emg(:) 
	 type(COMPLEX_POINTER_ARRAY_VII),pointer :: spts(:)
	 complex(prc),pointer        :: e (:,:,:,:)
	 complex(prc),pointer        :: h(:,:,:,:)
	 integer,intent(in)          :: it  
	 integer,intent(in)          :: nobox
 end subroutine

  subroutine scattering_fourier_out(emg,spts,it,nobox)
	use parameters	
	use physconsts
	use COMMON_SUBROUTINES
	use files
	 type(COMPLEX_POINTER_ARRAY_VI),pointer :: emg(:) 
	 type(COMPLEX_POINTER_ARRAY_VII),pointer :: spts(:)
	 integer,intent(in)          :: it  
	 integer,intent(in)          :: nobox
end subroutine
  
subroutine init_store_pts_fourier(spts,nobox)
	use files
	use parameters
	implicit none

	type(COMPLEX_POINTER_ARRAY_VII),pointer :: spts(:)
	integer                                 :: nobox
end subroutine
  
end interface

end module PROCESO_DATOS
    

module NEAR_FIELD
interface
  subroutine load_near_field(fdata,path,no_input,nolambda)
	use parameters
	use physconsts	
	use files
	use COMMON_SUBROUTINES
	use SDL_MESSAGES
	 implicit none	 	
	 complex(prc),pointer         :: fdata(:,:,:,:,:)
	 character(50)                :: path
	  integer                     :: no_input
	 integer,pointer              :: nolambda(:)
  end subroutine

  subroutine unload_near_field(fdata,path,no_input,nolambda)
	use parameters
	use physconsts	
	use files
	use COMMON_SUBROUTINES
	use SDL_MESSAGES
	implicit none	 	
	complex(prc),pointer         :: fdata(:,:,:,:,:)
	character(50)                :: path
	integer                     :: no_input
	integer,pointer              :: nolambda(:)
  end subroutine

  subroutine load_near_field_epsilon(eps_inv)
	use parameters
	use physconsts	
	use files
	use COMMON_SUBROUTINES
	use SDL_MESSAGES
	implicit none	 	
	complex(prc),pointer         :: eps_inv(:,:,:,:,:)
  end subroutine  
  
end interface   
end module NEAR_FIELD

!*******************************************************************************************************!
subroutine scattering_fourier(emg,spts,e,h,it,nobox)
use parameters	
use physconsts
use COMMON_SUBROUTINES
use files
 implicit none	 
 type(COMPLEX_POINTER_ARRAY_VI),pointer  :: emg(:) 
 type(COMPLEX_POINTER_ARRAY_VII),pointer :: spts(:)
 complex(prc),pointer        :: e (:,:,:,:)
 complex(prc),pointer        :: h(:,:,:,:)
 integer,intent(in)          :: it  
 integer,intent(in)          :: nobox
 integer ::ix,iy,iz,j,i_pts,bxi,bxf,byi,byf,bzi,bzf
 integer::iw,icm,cmo
 real(precision) :: w
 integer  :: npts(3)

if(nobox==1)then
bxi=sxi;bxf=sxf;byi=syi;byf=syf;bzi=szi;bzf=szf
else
bxi=mxi;bxf=mxf;byi=myi;byf=myf;bzi=mzi;bzf=mzf
end if

npts(1) = 2*6*(byf+1-byi)*(bzf+1-bzi)
npts(2) = 2*6*(bxf+1-bxi)*(bzf+1-bzi)
npts(3) = 2*6*(bxf+1-bxi)*(byf+1-byi)

do iw=1,nof 
w=fw(iw)  
if((iymax/=1).and.(ixmax/=1))then
if(system_model==0)then
 do j = 1,3
  do i_pts=1,npts(j)   
   if (spts(j)%pto(i_pts,1)<4) then
    emg(j)%fdata(i_pts,iw)= emg(j)%fdata(i_pts,iw) + &
                              & + dt*e(spts(j)%pto(i_pts,1),spts(j)%pto(i_pts,2),spts(j)%pto(i_pts,3),spts(j)%pto(i_pts,4))*cdexp(ci*w*it)   
   else
    emg(j)%fdata(i_pts,iw)= emg(j)%fdata(i_pts,iw) + &
                  & + dt*h(spts(j)%pto(i_pts,1)-3,spts(j)%pto(i_pts,2),spts(j)%pto(i_pts,3),spts(j)%pto(i_pts,4))*cdexp(ci*w*it)   
   end if
  end do
 end do
else
  j = 3
  do i_pts=1,npts(j)   
   if (spts(j)%pto(i_pts,1)<3) then
    emg(j)%fdata(i_pts,iw)= emg(j)%fdata(i_pts,iw) + &
                              & + dt*e(spts(j)%pto(i_pts,1),spts(j)%pto(i_pts,2),spts(j)%pto(i_pts,3),spts(j)%pto(i_pts,4))*cdexp(ci*w*it)   
   end if
   if((spts(j)%pto(i_pts,1)>3).and.(spts(j)%pto(i_pts,1)<6)) then
    emg(j)%fdata(i_pts,iw)= emg(j)%fdata(i_pts,iw) + &
                  & + dt*h(spts(j)%pto(i_pts,1)-3,spts(j)%pto(i_pts,2),spts(j)%pto(i_pts,3),spts(j)%pto(i_pts,4))*cdexp(ci*w*it)   
   end if
  end do
 
end if
else
 if((iymax==1).and.(ixmax/=1))then
 do j=1,3,2 !>>> Si 2D-System entonces la reserva de puntos en ymin e ymax resulta totalmente inecesaria
   do i_pts=1,npts(j)   
   if (spts(j)%pto(i_pts,1)<4) then
    emg(j)%fdata(i_pts,iw)= emg(j)%fdata(i_pts,iw) + &
                              & + dt*e(spts(j)%pto(i_pts,1),spts(j)%pto(i_pts,2),spts(j)%pto(i_pts,3),spts(j)%pto(i_pts,4))*cdexp(ci*w*it)   
   else
    emg(j)%fdata(i_pts,iw)= emg(j)%fdata(i_pts,iw) + &
                  & + dt*h(spts(j)%pto(i_pts,1)-3,spts(j)%pto(i_pts,2),spts(j)%pto(i_pts,3),spts(j)%pto(i_pts,4))*cdexp(ci*w*it)   
   end if
  end do
  end do
 end if

 if((ixmax==1).and.(iymax/=1))then
 do j=2,3 !>>> Si 2D-System entonces la reserva de puntos en ymin e ymax resulta totalmente inecesaria
   do i_pts=1,npts(j)   
   if (spts(j)%pto(i_pts,1)<4) then
    emg(j)%fdata(i_pts,iw)= emg(j)%fdata(i_pts,iw) + &
                              & + dt*e(spts(j)%pto(i_pts,1),spts(j)%pto(i_pts,2),spts(j)%pto(i_pts,3),spts(j)%pto(i_pts,4))*cdexp(ci*w*it)   
   else
    emg(j)%fdata(i_pts,iw)= emg(j)%fdata(i_pts,iw) + &
                  & + dt*h(spts(j)%pto(i_pts,1)-3,spts(j)%pto(i_pts,2),spts(j)%pto(i_pts,3),spts(j)%pto(i_pts,4))*cdexp(ci*w*it)   
   end if
  end do
  end do
 end if

 if((iymax==1).and.(ixmax==1))then
 j=3 
   do i_pts=1,npts(j)   
   if (spts(j)%pto(i_pts,1)<4) then
    emg(j)%fdata(i_pts,iw)= emg(j)%fdata(i_pts,iw) + &
                              & + dt*e(spts(j)%pto(i_pts,1),spts(j)%pto(i_pts,2),spts(j)%pto(i_pts,3),spts(j)%pto(i_pts,4))*cdexp(ci*w*it)   
   else
    emg(j)%fdata(i_pts,iw)= emg(j)%fdata(i_pts,iw) + &
                  & + dt*h(spts(j)%pto(i_pts,1)-3,spts(j)%pto(i_pts,2),spts(j)%pto(i_pts,3),spts(j)%pto(i_pts,4))*cdexp(ci*w*it)   
   end if
  end do 
 end if
end if
end do
return
end subroutine 

!*****************************************************************

subroutine scattering_fourier_out(emg,spts,it,nobox)
use parameters	
use physconsts
use COMMON_SUBROUTINES
use files
 type(COMPLEX_POINTER_ARRAY_VI),pointer :: emg(:) 
 type(COMPLEX_POINTER_ARRAY_VII),pointer :: spts(:)
 integer,intent(in)          :: it  
 integer,intent(in)          :: nobox
 integer ::ix,iy,iz,i,j,k
 character(5) :: jAscii
 character(5) :: noboxAscii
 integer::iw,icm,cmo
 integer  :: npts(3)
 real(precision),pointer :: Jabs(:,:)
 real(precision)         :: Jw
 integer                 :: bxi,bxf,byi,byf,bzi,bzf
 
if(nobox==1)then
bxi=sxi;bxf=sxf;byi=syi;byf=syf;bzi=szi;bzf=szf
else
bxi=mxi;bxf=mxf;byi=myi;byf=myf;bzi=mzi;bzf=mzf
end if

call inttoASCII(nobox,noboxAscii)

npts(1) = 2*6*(byf+1-byi)*(bzf+1-bzi)
npts(2) = 2*6*(bxf+1-bxi)*(bzf+1-bzi)
npts(3) = 2*6*(bxf+1-bxi)*(byf+1-byi)

 allocate(Jabs(1:6,nof)) !>>> 1:6 -> Para cada una de las 6 "superficies de medida" <<<!
 Jabs = 0.0
 
 do iw=1,nof 
   do j = 1,3   
	 do i_pts=1,npts(j),6	   
	   select case(j)
	   case(1)
	   if(ixmax/=1)then
	    if(i_pts<=npts(j)/2) then
	     Jabs(1,iw)= Jabs(1,iw) - real(emg(j)%fdata(i_pts+1,iw)*conjg(emg(j)%fdata(i_pts+5,iw))-emg(j)%fdata(i_pts+2,iw)*conjg(emg(j)%fdata(i_pts+4,iw)))      
		else
		 Jabs(2,iw)= Jabs(2,iw) + real(emg(j)%fdata(i_pts+1,iw)*conjg(emg(j)%fdata(i_pts+5,iw))-emg(j)%fdata(i_pts+2,iw)*conjg(emg(j)%fdata(i_pts+4,iw)))      
		end if
	   end if
	   case(2)
	   if(iymax/=1)then
	    if(i_pts<=npts(j)/2) then
	     Jabs(3,iw)= Jabs(3,iw) - real(emg(j)%fdata(i_pts+2,iw)*conjg(emg(j)%fdata(i_pts+3,iw))-emg(j)%fdata(i_pts,iw)*conjg(emg(j)%fdata(i_pts+5,iw)))       
		else
		 Jabs(4,iw)= Jabs(4,iw) + real(emg(j)%fdata(i_pts+2,iw)*conjg(emg(j)%fdata(i_pts+3,iw))-emg(j)%fdata(i_pts,iw)*conjg(emg(j)%fdata(i_pts+5,iw)))       
		end if        
	   end if
	   case(3)
	    if(i_pts<=npts(j)/2) then
	     Jabs(5,iw)= Jabs(5,iw) - real(emg(j)%fdata(i_pts,iw)*conjg(emg(j)%fdata(i_pts+4,iw))-emg(j)%fdata(i_pts+1,iw)*conjg(emg(j)%fdata(i_pts+3,iw)))      
		else
		 Jabs(6,iw)= Jabs(6,iw) + real(emg(j)%fdata(i_pts,iw)*conjg(emg(j)%fdata(i_pts+4,iw))-emg(j)%fdata(i_pts+1,iw)*conjg(emg(j)%fdata(i_pts+3,iw)))      
		end if
	   end select
	 end do
   end do
end do

Jabs = (eps0*q/(2.0*dt))*Jabs

do j=1,6
  call inttoASCII(j,jAscii)
  open(unit=out6,file=trim(noboxAscii)//'_scatt' // trim(jAscii) // '.dat',status='unknown',action='write',position='append')    
  write(out6,*) '#',it

  do iw=1,nof
    write(out6,*) Jabs(j,iw)
  end do

  close(out6)
end do

!>>> Total current
open(unit=out6,file=trim(noboxAscii)//'_scatt.dat',status='unknown',action='write',position='append')    
 
  write(out6,*) '#',it
  do iw=1,nof
    write(out6,*) Jabs(1,iw)+Jabs(2,iw)+Jabs(3,iw)+Jabs(4,iw)+Jabs(5,iw)+Jabs(6,iw)
  end do

close(out6)
!    Total current <<<

deallocate(Jabs)
return
    end subroutine 
   

!**************************************************!

subroutine near_field_map(bbox)
use parameters
use physconsts	
use files
use COMMON_SUBROUTINES
use SDL_MESSAGES
use NEAR_FIELD
 implicit none	 
 real(prc),pointer            :: bbox(:,:,:)	
 complex(prc),pointer         :: fdata(:,:,:,:,:),fdata_aux(:,:,:,:,:)  
 type(EMR)                    :: emrw
 complex(precision)           :: ex,ey,ez,hx,hy,hz
 real(precision)              :: e,h,norm
 integer                      :: i,iaux,j,jaux,k,kaux,iw,icm,id ,xnorm,znorm,ynorm
 integer                      :: li
 integer                      :: store_r
 character(5)                 :: nfr
 character(5)                 :: nz 
 character(7)                 :: dir
 integer                      :: icase,ilayer
 integer                      :: nolayer 
 character(50)                :: name
 integer                      :: abcisa0,abcisaF,ordenada0,ordenadaF
 complex(precision)           :: fH,fE
 integer,pointer              :: nolambda(:)
 integer                      :: no_input
 integer                      :: xonorm,yonorm,zonorm 
 logical                      :: flag
 logical,parameter            :: default=.true.
 logical                      :: flagIP=.false. 
 integer                      :: indx

if(not(default))then
no_input=nof
xnorm=0
ynorm=0
znorm=0
end if 

print*,"No of lamdas( inputs)  ="
if(default) read*,no_input
if(flagIP)  read*,no_input
!no_input = 1
if((default).or.(flagIP))then
print*,"x for being normalized ="
 read*,xnorm
print*,"y for being normalized ="
read*,ynorm
print*,"z for being normalized ="
read*,znorm
end if


allocate(nolambda(no_input))

if(not(default)) then
do iw=1, no_input
 nolambda(iw)=iw
end do
end if


do iw=1,no_input
    print*,"Wavelenght label ="
	if(default) read*,nolambda(iw)	
	if(flagIP)  nolambda(iw)=iw
	wo=1/(2.0*pi*c0*dt/(fw(nolambda(iw))))        
	print*,"Wavelength (nm) =",1.0/(f_abohr*wo) 
end do


! >>> FORMAT LIST
1 format(I5,A1,I5,1F30.12)
2 format(I5,A1,I5,2F30.12)
3 format(I5,A1,I5,3F30.12)
4 format(I5,A1,I5,4F30.12)
5 format(I5,A1,I5,5F30.12)
6 format(I5,A1,I5,6F30.12)
7 format(I5,A1,I5,7F30.12)
8 format(I5,A1,I5,8F30.12)
!     FORMAT LIST

name=""

call load_near_field(fdata,name,no_input,nolambda)

! >>> Process DATA
 do icase=1,3
 select case(icase)
  case(1)
	 rewind(unit=inFDTD)  
	 do while(msg/=SERS_X)
		read(inFDTD,*) msg
	 end do		
	 read(inFDTD,*) nolayer 
  case(2)
     rewind(unit=inFDTD)  
	 do while(msg/=SERS_Y)
		read(inFDTD,*) msg
	 end do		
	 read(inFDTD,*) nolayer 
  case(3)
     rewind(unit=inFDTD)  
	 do while(msg/=SERS_Z)
		read(inFDTD,*) msg
	 end do		
	 read(inFDTD,*) nolayer 
 end select

ilayer=1
do while(ilayer<=nolayer)
  read(inFDTD,*) store_r  

 select case (icase)
   case (1) ! >>> Plano X=cte <<<      
      abcisa0   = myi
	  abcisaF   = myf
	  ordenada0 = mzi
	  ordenadaF = mzf
   case (2) ! >>> Plano Y=cte <<<	     
      abcisa0   = mxi
	  abcisaF   = mxf
	  ordenada0 = mzi
	  ordenadaF = mzf	   	      
   case (3) ! >>> Plano Z=cte <<<	      
      abcisa0   = mxi
	  abcisaF   = mxf
	  ordenada0 = myi
	  ordenadaF = myf	  
 end select	



if(.true.)then
! >>> E(r,w) Field en unidades atómicas de campo eléctrico <<< !
if((emmode==0).or.(emmode==1))then
!fE = dt/Q
fE = 1.0
 do iw=1,no_input 
  li = nolambda(iw)
  call inttoASCII(li,nfr)      
  
  k= store_r    
  call inttoASCII(k,nz)
     
    select case (icase)
	   case (1) ! >>> Plano X=cte <<< 
		dir = 'X(cte)='	    			  				   
	   case (2) ! >>> Plano Y=cte <<<		   	
		dir = 'Y(cte)='	      			  
	   case (3) ! >>> Plano Z=cte <<<				
		dir = 'Z(cte)='	     			  
	end select	

   	open(unit=out5,file=trim(idx_name)//'Ex'//nfr//dir//nz//'.dat',status='replace',action='write')    
	open(unit=out6,file=trim(idx_name)//'Ey'//nfr//dir//nz//'.dat',status='replace',action='write')
	open(unit=out7,file=trim(idx_name)//'Ez'//nfr//dir//nz//'.dat',status='replace',action='write')
	open(unit=out8,file=trim(idx_name)//'E'//nfr//dir//nz//'.dat',status='replace',action='write')
		 		   
	write(out5,*) 
	write(out5,*) '#',iw    
	write(out6,*) 
	write(out6,*) '#',iw    
	write(out7,*) 
	write(out7,*) '#',iw    
	!write(out8,*) 
	!write(out8,*) '#',iw    
	    if(xnorm/=0)then	    
		   ! >>> Plano Y=cte <<<		   	
		   ex =fE*fdata(1,xnorm,ynorm,znorm,iw)
		   ey =fE*fdata(2,xnorm,ynorm,znorm,iw)
		   ez =fE*fdata(3,xnorm,ynorm,znorm,iw)	
	       norm = dsqrt(cdabs(conjg(ex)*ex)+cdabs(conjg(ey)*ey)+cdabs(conjg(ez)*ez))
		else
		   norm = 1.0           
           print*,"E-field normalized by... ="; read*,norm
		end if
	    do i=abcisa0,abcisaF
		    do j=ordenada0,ordenadaF
			 select case (icase)
			   case (1) ! >>> Plano X=cte <<< 
			      ex =fE*fdata(1,k,i,j,iw)
				  ey =fE*fdata(2,k,i,j,iw)
				  ez =fE*fdata(3,k,i,j,iw)					  				  				   					  				   					  				   
			   case (2) ! >>> Plano Y=cte <<<		   	
			      ex =fE*fdata(1,i,k,j,iw)
				  ey =fE*fdata(2,i,k,j,iw)
				  ez =fE*fdata(3,i,k,j,iw)						  
			   case (3) ! >>> Plano Z=cte <<<				
			      ex =fE*fdata(1,i,j,k,iw)
				  ey =fE*fdata(2,i,j,k,iw)
				  ez =fE*fdata(3,i,j,k,iw)				 				  
			 end select	 
			  	   
			 e = cdabs(conjg(ex)*ex)+cdabs(conjg(ey)*ey)+cdabs(conjg(ez)*ez) ! So expresed, cdabs takes a real valued argument so cdabs = || so here e=e^2

			 e = e/norm**2
             ex=ex/norm
             ey=ey/norm
             ez=ez/norm

             write(out5,*)  dble(ex),imag(ex)
             write(out6,*)  dble(ey),imag(ey)
             write(out7,*)  dble(ez),imag(ez)
             
			 write(out8,*)  dsqrt(e)
	      end do
	   end do
     close(out5)	       
	 close(out6)
	 close(out7)
	 close(out8)
  end do
end if
	
! >>> H(r,w) Field <<< !
!>>> Fi(r,w) = dt*(c*mu0)*(eps0*Q/dt)*(1/Qi)*`Hi´ -> [F] = [E(r,w)]
if((emmode==0).or.(emmode==2))then
!fH = 1.0/c0
fH=1.0
 
 do iw=1,no_input  
  li = nolambda(iw)
  call inttoASCII(li,nfr)         
  
  k= store_r
  call inttoASCII(k,nz)

   	select case (icase)
	   case (1) ! >>> Plano X=cte <<< 
		dir = 'X(cte)='	    			  				   
	   case (2) ! >>> Plano Y=cte <<<		   	
		dir = 'Y(cte)='	      			  
	   case (3) ! >>> Plano Z=cte <<<				
		dir = 'Z(cte)='	     			  
	end select	

   	open(unit=out5,file=trim(idx_name)//'Hx'//nfr//dir//nz//'.dat',status='replace',action='write')    
	open(unit=out6,file=trim(idx_name)//'Hy'//nfr//dir//nz//'.dat',status='replace',action='write')
	open(unit=out7,file=trim(idx_name)//'Hz'//nfr//dir//nz//'.dat',status='replace',action='write')
	open(unit=out8,file=trim(idx_name)//'H'//nfr//dir//nz//'.dat',status='replace',action='write')
	 		   
	write(out5,*) 
	write(out5,*) '#',iw    
	write(out6,*) 
	write(out6,*) '#',iw    
	write(out7,*) 
	write(out7,*) '#',iw    
	write(out8,*) 
	write(out8,*) '#',iw   	
	   if(xnorm/=0)then	    
		   ! >>> Plano Y=cte <<<		   	
		   hx =fH*fdata(4,xnorm,ynorm,znorm,iw)
		   hy =fH*fdata(5,xnorm,ynorm,znorm,iw)
		   hz =fH*fdata(6,xnorm,ynorm,znorm,iw)	
	       norm = dsqrt(cdabs(conjg(hx)*hx)+cdabs(conjg(hy)*hy)+cdabs(conjg(hz)*hz))
		else
		   norm = 1.0
           print*,"H-field normalized by... ="; read*,norm
		end if
    
	    do i=abcisa0,abcisaF
		    do j=ordenada0,ordenadaF
			 select case (icase)
			   case (1) ! >>> Plano X=cte <<< 			     
				  hx =fH*fdata(4,k,i,j,iw)
				  hy =fH*fdata(5,k,i,j,iw)
				  hz =fH*fdata(6,k,i,j,iw)			     				 					   
			   case (2) ! >>> Plano Y=cte <<<				      
				  hx =fH*fdata(4,i,k,j,iw)
				  hy =fH*fdata(5,i,k,j,iw)
				  hz =fH*fdata(6,i,k,j,iw)				  				  
			   case (3) ! >>> Plano Z=cte <<<			      
				  hx =fH*fdata(4,i,j,k,iw)
				  hy =fH*fdata(5,i,j,k,iw)
				  hz =fH*fdata(6,i,j,k,iw)					  
			 end select			     	
			 h = cdabs(conjg(hx)*hx)+cdabs(conjg(hy)*hy)+cdabs(conjg(hz)*hz)				  
			 
			 !h= h/norm**2
			 hx=hx/norm
			 hy=hy/norm
			 hz=hz/norm

			 write(out5,2) i," ",j,dble(hx),imag(hx)
			 write(out6,2) i," ",j,dble(hy),imag(hy)
			 write(out7,2) i," ",j,dble(hz),imag(hz)
		 	   
			 write(out8,*)  dsqrt(h)/norm             
             
	      end do
	   end do
     close(out5)	       
	 close(out6)
	 close(out7)
	 close(out8)	      
	end do	
  end if

else
fE = dt/Q
do iw=1,no_input 
do i=1,ixmax
 do j=1,iymax
   flag = .true.
   do k=hi,hf
   if(flag)then
     if(bbox(i,j,k)==1.0)then
	    ex = fE*fdata(1,i,j,k,iw)
		ey = fE*fdata(2,i,j,k,iw)
		ez = fE*fdata(3,i,j,k,iw)
	    e  = cdabs(conjg(ex)*ex)+cdabs(conjg(ey)*ey)+cdabs(conjg(ez)*ez)
		write(meshgraph,*) i," ",j," "," ",hf-k,dsqrt(e)
		flag = .false.	
     end if
   end if
   end do
 end do
end do
end do

end if

  ilayer = ilayer+1
  end do
end do

!    Process DATA <<<
deallocate(fdata)
deallocate(nolambda)
return
end subroutine 

!****************************************************************+

subroutine load_near_field(fdata,path,no_input,nolambda)
use parameters
use physconsts	
use files
use COMMON_SUBROUTINES
use SDL_MESSAGES
 implicit none	 	
 complex(prc),pointer         :: fdata(:,:,:,:,:)
 character(50)                :: path
 integer                      :: no_input
 integer,pointer              :: nolambda(:)
 type(EMR)                    :: emrw 
 integer                      :: i,j,k,iw
 integer                      :: li 
 character(5)                 :: nfr  
 

select case(emmode)
 case (0)
   allocate(fdata(6,mxi:mxf,myi:myf,mzi:mzf,no_input)) 	
 case (1)
   allocate(fdata(1:3,mxi:mxf,myi:myf,mzi:mzf,no_input)) 	
 case (2)
   allocate(fdata(4:6,mxi:mxf,myi:myf,mzi:mzf,no_input)) 	
end select 

! >>> E(r,w) Field <<< !
if((emmode==0).or.(emmode==1))then
 print*, "...loading E-Field"
 do iw=1,no_input  
  li = nolambda(iw)
  call inttoASCII(li,nfr)
     
   open (unit=out5,file=trim(idx_name)//trim(path)//'E_l(nm)='//nfr//'.dat',action='read',FORM="UNFORMATTED",POSITION="REWIND")
	
	 do i=sxi,sxf
      do j=syi,syf
        do k=szi,szf
		    read(out5) emrw		   
			if((i>=mxi).and.(i<=mxf).and.(j>=myi).and.(j<=myf).and.(k>=mzi).and.(k<=mzf))then 
		    fdata(1,i,j,k,iw)=emrw%ex
			fdata(2,i,j,k,iw)=emrw%ey
			fdata(3,i,j,k,iw)=emrw%ez		   
			end if
	      end do
	   end do     	       
     end do
	 close(out5)
 end do
end if

! >>> H(r,w) Field <<< !
if((emmode==0).or.(emmode==2))then
 print*, "...loading H-Field"
 do iw=1,no_input  
  li = nolambda(iw)
  call inttoASCII(li,nfr)     	
	open (unit=out5,file=trim(idx_name)//trim(path)//'H_l(nm)='//nfr//'.dat',action='read',FORM="UNFORMATTED",POSITION="REWIND")
	    
	do i=sxi,sxf
      do j=syi,syf
        do k=szi,szf	
		 read(out5) emrw		
		 if((i>=mxi).and.(i<=mxf).and.(j>=myi).and.(j<=myf).and.(k>=mzi).and.(k<=mzf))then    
		 fdata(4,i,j,k,iw)=emrw%ex
		 fdata(5,i,j,k,iw)=emrw%ey
		 fdata(6,i,j,k,iw)=emrw%ez			 
		 end if
	     end do
	   end do
	end do
    close(out5)
 end do
end if
return
end subroutine 

!***********************************************************************************************!

subroutine unload_near_field(fdata,path,no_input,nolambda)
use parameters
use physconsts	
use files
use COMMON_SUBROUTINES
use SDL_MESSAGES
 implicit none	 	
 complex(prc),pointer         :: fdata(:,:,:,:,:)
 character(50)                :: path
 integer                      :: no_input
 integer,pointer              :: nolambda(:)
 type(EMR)                    :: emrw 
 integer                      :: i,j,k,iw
 integer                      :: li 
 character(5)                 :: nfr  
 
! >>> E(r,w) Field <<< !
nof=no_input
if((emmode==0).or.(emmode==1))then
 do iw=1,nof   
  li = iw
  call inttoASCII(li,nfr)       
	open (unit=out5,file=trim(idx_name)//trim(path)//'UnE_l(nm)='//nfr//'.dat',status='replace',action='write',FORM="UNFORMATTED",POSITION="REWIND")	  
	
	do i=mxi,mxf
      do j=myi,myf
        do k=mzi,mzf	    		   
			emrw%i=i	   		   
			emrw%j=j
			emrw%k=k
			emrw%ex=fdata(1,i,j,k,iw)
			emrw%ey=fdata(2,i,j,k,iw)
			emrw%ez=fdata(3,i,j,k,iw)	
		    write(out5) emrw	
	      end do
	   end do     	       
     end do
	 close(out5)
 end do
end if

! >>> H(r,w) Field <<< !
if((emmode==0).or.(emmode==2))then
 do iw=1,nof  
   li = iw
  call inttoASCII(li,nfr)  
      
  open (unit=out5,file=trim(idx_name)//trim(path)//'UnH_l(nm)='//nfr//'.dat',status='replace',action='write',FORM="UNFORMATTED",POSITION="REWIND")	 		   
	    
	do i=mxi,mxf
      do j=myi,myf
        do k=mzi,mzf	    		    		   		
		    emrw%i=i	   		   
			emrw%j=j
			emrw%k=k
			emrw%ex=fdata(4,i,j,k,iw)
			emrw%ey=fdata(5,i,j,k,iw)
			emrw%ez=fdata(6,i,j,k,iw)
		    write(out5) emrw					     		   
	     end do
	   end do
	end do
    close(out5)
 end do
end if
return
    end subroutine 

!********************************************************************************************!

subroutine init_store_pts_fourier(spts,nobox)
use files
use parameters
implicit none

type(COMPLEX_POINTER_ARRAY_VII),pointer :: spts(:)
integer                                 :: nobox
integer :: i,j,ix,iy,iz,ip,bxi,bxf,byi,byf,bzi,bzf

integer  :: npts(3),i_pts

write(logfile,*) '...Running - Subroutine: init_store_pts_fourier'

if(nobox==1)then
bxi=sxi;bxf=sxf;byi=syi;byf=syf;bzi=szi;bzf=szf
else
bxi=mxi;bxf=mxf;byi=myi;byf=myf;bzi=mzi;bzf=mzf
end if

!>>> X = bxi y bxf
if(ixmax/=1) then
j = 1
	do iy=1,(byf+1-byi)
		do iz=1,(bzf+1-bzi)
			do ip=1,6             
			  i=ip+6*(iy-1+(byf+1-byi)*(iz-1))
	          spts(j)%pto(i,1) = ip		
			  spts(j)%pto(i,2) = bxi
			  spts(j)%pto(i,3) = iy+(byi-1)		  
			  spts(j)%pto(i,4) = iz+(bzi-1)
			enddo
		enddo
	enddo 

	  ! Los siguientes 4NxNy puntos corresponden a la segunda superficie de cálculo
	  ! antes de atravesar el centro dispersor
     do iy=1,(byf+1-byi)
		do iz=1,(bzf+1-bzi)
			do ip=1,6
			  i=6*(byf+1-byi)*(bzf+1-bzi) + ip +6*(iy-1+(byf+1-byi)*(iz-1))
			  spts(j)%pto(i,1) = ip	
			  spts(j)%pto(i,2) = bxf
			  spts(j)%pto(i,3) = iy+(byi-1)		  
			  spts(j)%pto(i,4) = iz+(bzi-1)              
			enddo
		enddo
	enddo 
end if
!    X = bxi y bxf <<<

!>>> Y = byi y byf
if(iymax/=1) then
j = 2
	do ix=1,(bxf+1-bxi)
		do iz=1,(bzf+1-bzi)
			do ip=1,6
			  i=ip+6*(ix-1+(bxf+1-bxi)*(iz-1))
			  spts(j)%pto(i,1) = ip	
			  spts(j)%pto(i,2) = ix+(bxi-1)
			  spts(j)%pto(i,3) = byi			  
			  spts(j)%pto(i,4) = iz+(bzi-1)
			enddo
		enddo
	enddo 

	  ! Los siguientes 4NxNy puntos corresponden a la segunda superficie de cálculo
	  ! antes de atravesar el centro dispersor
  do ix=1,(bxf+1-bxi)
		do iz=1,(bzf+1-bzi)	
			do ip=1,6
			  i=6*(bxf+1-bxi)*(bzf+1-bzi) + ip+6*(ix-1+(bxf+1-bxi)*(iz-1))
			  spts(j)%pto(i,1) = ip	
			  spts(j)%pto(i,2) = ix+(bxi-1)
			  spts(j)%pto(i,3) = byf			                
			  spts(j)%pto(i,4) = iz+(bzi-1)
			enddo
		enddo
	enddo 
!    Y = byi y byf <<<   
end if

!>>> Z = bzi y bzf
j = 3
	do ix=1,(bxf+1-bxi)
		do iy=1,(byf+1-byi)
			do ip=1,6
			  i=ip+6*(ix-1+(bxf+1-bxi)*(iy-1))
			  spts(j)%pto(i,1) = ip	
			  spts(j)%pto(i,2)=ix+(bxi-1)
			  spts(j)%pto(i,3)=iy+(byi-1)				  
			  spts(j)%pto(i,4) = bzi 
			enddo
		enddo
	enddo 

	  ! Los siguientes 4NxNy puntos corresponden a la segunda superficie de cálculo
	  ! antes de atravesar el centro dispersor
	do ix=1,(bxf+1-bxi)
		do iy=1,(byf+1-byi)
			do ip=1,6
			  i=6*(bxf+1-bxi)*(byf+1-byi) + ip+6*(ix-1+(bxf+1-bxi)*(iy-1))
			  spts(j)%pto(i,1) = ip	
			  spts(j)%pto(i,2)=ix+(bxi-1)
			  spts(j)%pto(i,3)=iy+(byi-1)			                
			  spts(j)%pto(i,4) = bzf 
			enddo
		enddo
	enddo 
!    Z = bzi y bzf <<<
return
end subroutine init_store_pts_fourier