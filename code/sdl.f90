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

module SDL_MEMBERS 
use parameters
  integer         :: NO_FILES   
  integer         :: xi(3),xf(3),yi(3),yf(3),zi(3),zf(3)
  integer         :: SDLField(6)
  integer         :: SDLModuleE  
  integer         :: SDL_copy_flag
  integer         :: SDL_copy_no  
  integer         :: SDL_rotate_90,SDL_rotate_180,SDL_scale,scale_factor

  ! Rutas de acceso  
  integer, parameter::LONGITUD_ARCHIVO=100
  character(LONGITUD_ARCHIVO), save::RUTA_SDL  
end module SDL_MEMBERS

!>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<!

module SDL_PRIVATE
  interface
   subroutine SDL_output(ino,NoComponente,em,pInfoSDL,maxem,LayerSDL)
	  use parameters
	  use SDL_MEMBERS 
	  implicit none
	  integer, intent(IN)::ino,NoComponente,em
   	  complex(prc), pointer::pInfoSDL(:,:,:,:)	
	  real(prc),pointer::maxem(:,:)
	  integer          :: LayerSDL
   end subroutine
  
  subroutine SDL_input(ino,NoComponente,em,pInfoSDL,maxem,LayerSDL)
	  use parameters
	  use SDL_MEMBERS 
	  implicit none
	  integer, intent(IN)::ino,NoComponente,em
   	  complex(prc), pointer::pInfoSDL(:,:,:,:)	
	  real(prc),pointer::maxem(:,:)
	  integer          :: LayerSDL
   end subroutine
  end interface    
end module SDL_PRIVATE

! **************************************************** !

subroutine SDL_output(ino,NoComponente,em,pInfoSDL,maxem,LayerSDL)
use parameters
use physconsts
use SDL_MEMBERS 
use COMMON_SUBROUTINES
implicit none
  integer, intent(IN)::ino,NoComponente,em
  character(5) ::charno,charno2
  complex(prc), pointer::pInfoSDL(:,:,:,:)	
  real(prc),pointer::maxem(:,:)    
  integer ::LayerSDL
  integer ix,iy,iz,i,j  
  character(LONGITUD_ARCHIVO) Archivo  
  character(4) EXT,EXT2
  real(prc), pointer::Output(:,:)
  real(prc), pointer::Output_aux(:,:)
  integer           ::ixo,ixf,iyo,iyf
  real(prc) :: glog
  integer::Escala
  real(prc) :: SDLnorm

  select case(LayerSDL)
    case(1) !>>> Corte X=cte (LayerSDL=1)<<<!
		ixo =yi(1)
		ixf =yf(1)
		iyo =zi(1)
		iyf =zf(1)
		allocate(Output(ixo:ixf,iyo:iyf))
	case(2) !>>> Corte Y=cte (LayerSDL=2)<<<!	
		ixo =xi(2)
		ixf =xf(2)
		iyo =zi(2)
		iyf =zf(2)
		allocate(Output(ixo:ixf,iyo:iyf))
	case(3) !>>> Corte Z=cte (LayerSDL=3)<<<!		
		ixo =xi(3)
		ixf =xf(3)
		iyo =yi(3)
		iyf =yf(3)
		allocate(Output(ixo:ixf,iyo:iyf))
  end select     
    
  !em=1----> Campo Eléctrico
  !em=2----> Campo Magnético
  if (em==1) then
	   !e1 = Ex
	  if(NoComponente.EQ.1) then
	   EXT ='.e1'
	  end if
	   !e2 = Ey
	  if(NoComponente.EQ.2) then
	   EXT ='.e2'
	  end if
	   !e3 = Ez
	  if(NoComponente.EQ.3) then
	   EXT = '.e3'
      end if   
  else
  	   !h1 = Hx
	  if(NoComponente.EQ.1) then
	   EXT ='.h1'	   
	  end if
	   !h2 = Hy
	  if(NoComponente.EQ.2) then
	   EXT ='.h2'
	  end if
	   !h3 = Hz
	  if(NoComponente.EQ.3) then
	   EXT = '.h3'
      end if      
  end if
  
  EXT2=''
  if((em==1).or.(em==2))then
  select case(SDLField(NoComponente))
   case(1)
   EXT2='_R'!R=Real
   case(2)
   EXT2='_I'!I=Imaginary
   case(3)
   EXT2='_MC' !M=Modulus complex sqrt(real^2+imag^2)
   case(4)
   EXT2='_MR' !M=Modulus real sqrt(real^2)
   case(5)
   EXT2='_MI' !M=Modulus imag sqrt(imag^2)
   end select      
  end if

 if(em==0)then !>>> Modulo campo eléctrico <<<!
   EXT ='.e'
   if(SDLModuleE==1)then
   EXT2='_E'
   else 
   EXT2='_H'
   end if 
 end if
  
 if(em==-1)then
   EXT='.h3' !Warning (26nov18): Hz component is not used too much so I used this variable for charge density, but this must be fixed.
 end if

 call inttoASCII(ino,charno)

 if(RUTA_SDL/='cluster')then
 select case (LayerSDL)
 case (1)
   call inttoASCII(xi(1),charno2)
   Archivo='X_'//trim(charno2)//'layer'//trim(RUTA_SDL)//trim(EXT2)//trim(charno)//EXT 
 case (2)
   call inttoASCII(yi(2),charno2)
   Archivo='Y_'//trim(charno2)//'layer'//trim(RUTA_SDL)//trim(EXT2)//trim(charno)//EXT 
 case (3)
   call inttoASCII(zi(3),charno2)
   Archivo='Z_'//trim(charno2)//'layer'//trim(RUTA_SDL)//trim(EXT2)//trim(charno)//EXT 
 end select			
 else !>>> CLUSTER
   select case (LayerSDL)
   case (1)
   call inttoASCII(xi(1),charno2)
   Archivo='./X_SDL'//trim(charno2)//'layer'//trim(EXT2)//trim(charno)//EXT 
   case (2)
   call inttoASCII(yi(2),charno2)
   Archivo='./Y_SDL'//trim(charno2)//'layer'//trim(EXT2)//trim(charno)//EXT 
   case (3)
   call inttoASCII(zi(3),charno2)
   Archivo='./Z_SDL'//trim(charno2)//'layer'//trim(EXT2)//trim(charno)//EXT 
 end select	
 end if	 	

If((SDL_Xextralayers==1).and.(SDL_Yextralayers==1).and.(SDL_Zextralayers==1)) then !OLD mode 
 if(RUTA_SDL/='cluster')then
 select case (LayerSDL)
 case (1)   
   Archivo='X_'//trim(RUTA_SDL)//trim(EXT2)//trim(charno)//EXT 
 case (2)   
   Archivo='Y_'//trim(RUTA_SDL)//trim(EXT2)//trim(charno)//EXT 
 case (3)   
   Archivo='Z_'//trim(RUTA_SDL)//trim(EXT2)//trim(charno)//EXT 
 end select			
 else !>>> CLUSTER
   select case (LayerSDL)
   case (1)   
   Archivo='./X_SDL'//trim(EXT2)//trim(charno)//EXT 
   case (2)
   Archivo='./Y_SDL'//trim(EXT2)//trim(charno)//EXT 
   case (3)   
   Archivo='./Z_SDL'//trim(EXT2)//trim(charno)//EXT 
 end select	
 end if	 
 end if
 
!if(em/=0)then
Select case (em)
case(1,2)
 select case(LayerSDL)
    case(1) !>>> Corte X=cte <<<!
	  ix = xi(1)
      do iz=zi(1),zf(1)
        do iy=yi(1),yf(1)       		
		!Output(iy,iz) =  cdsqrt(dcmplx(conjg(pInfoSDL(NoComponente,ix,iy,iz))*pInfoSDL(NoComponente,ix,iy,iz)))
		!Output(iy,iz) =   sqrt(conjg(pInfoSDL(NoComponente,ix,iy,iz))*pInfoSDL(NoComponente,ix,iy,iz))
		!Output(iy,iz) =   dble(pInfoSDL(NoComponente,ix,iy,iz))!For real part of E-field & H-field
        select case(SDLField(NoComponente))
        case(1)
            Output(iy,iz) =  dble(pInfoSDL(NoComponente,ix,iy,iz)) !Real
        case(2)
            Output(iy,iz) = imag(pInfoSDL(NoComponente,ix,iy,iz))!Imag
        case(3)
            Output(iy,iz) =  sqrt(conjg(pInfoSDL(NoComponente,ix,iy,iz))*pInfoSDL(NoComponente,ix,iy,iz)) !Mod
        case(4) 
            Output(iy,iz) =  sqrt(dble(pInfoSDL(NoComponente,ix,iy,iz))*dble(pInfoSDL(NoComponente,ix,iy,iz))) !Abs(Real)
        case(5) 
            Output(iy,iz) =  sqrt(imag(pInfoSDL(NoComponente,ix,iy,iz))*imag(pInfoSDL(NoComponente,ix,iy,iz))) !Abs(Imag)
        end select
        end do
      end do	
	case(2) !>>> Corte Y=cte <<<!	 
	  iy = yi(2)
      do iz=zi(2),zf(2)
        do ix=xi(2),xf(2)		  
        !Output(ix,iz) = cdsqrt(dcmplx(conjg(pInfoSDL(NoComponente,ix,iy,iz))*pInfoSDL(NoComponente,ix,iy,iz))) 		
		!Output(ix,iz) =  dble(pInfoSDL(NoComponente,ix,iy,iz))!For real part of E-field & H-field
        select case(SDLField(NoComponente))
        case(1)
            Output(ix,iz) =  dble(pInfoSDL(NoComponente,ix,iy,iz)) !Real
        case(2)
            Output(ix,iz) = imag(pInfoSDL(NoComponente,ix,iy,iz))!Imag
        case(3)
            Output(ix,iz) =  sqrt(conjg(pInfoSDL(NoComponente,ix,iy,iz))*pInfoSDL(NoComponente,ix,iy,iz)) !Mod
        case(4) 
            Output(ix,iz) =  sqrt(dble(pInfoSDL(NoComponente,ix,iy,iz))*dble(pInfoSDL(NoComponente,ix,iy,iz))) !Abs(Real)
        case(5) 
            Output(ix,iz) =  sqrt(imag(pInfoSDL(NoComponente,ix,iy,iz))*imag(pInfoSDL(NoComponente,ix,iy,iz))) !Abs(Imag)
        end select
        end do
      end do
	case(3) !>>> Corte Z=cte <<<!
	  iz = zi(3)
      do iy=yi(3),yf(3)               
        do ix=xi(3),xf(3)				
		!Output(ix,iy) = cdsqrt(dcmplx(conjg(pInfoSDL(NoComponente,ix,iy,iz))*pInfoSDL(NoComponente,ix,iy,iz)))
		!Output(ix,iy) = dble(pInfoSDL(NoComponente,ix,iy,iz))!For real part of E-field & H-field        
        select case(SDLField(NoComponente))
        case(1)
            Output(ix,iy) =  dble(pInfoSDL(NoComponente,ix,iy,iz)) !Real
        case(2)
            Output(ix,iy) = imag(pInfoSDL(NoComponente,ix,iy,iz))!Imag
        case(3)
            Output(ix,iy) =  sqrt(conjg(pInfoSDL(NoComponente,ix,iy,iz))*pInfoSDL(NoComponente,ix,iy,iz)) !Mod
        case(4) 
            Output(ix,iy) =  sqrt(dble(pInfoSDL(NoComponente,ix,iy,iz))*dble(pInfoSDL(NoComponente,ix,iy,iz))) !Abs(Real)
        case(5) 
            Output(ix,iy) =  sqrt(imag(pInfoSDL(NoComponente,ix,iy,iz))*imag(pInfoSDL(NoComponente,ix,iy,iz))) !Abs(Imag)
        end select
        end do
      end do   	
    end select  
!else
case(0)
 select case(LayerSDL)
    case(1) !>>> Corte X=cte <<<!
	  ix = xi(1)
      do iz=zi(1),zf(1)
        do iy=yi(1),yf(1)
		glog = conjg(pInfoSDL(1,ix,iy,iz))*pInfoSDL(1,ix,iy,iz)
		glog = glog+conjg(pInfoSDL(2,ix,iy,iz))*pInfoSDL(2,ix,iy,iz)
		glog = glog+conjg(pInfoSDL(3,ix,iy,iz))*pInfoSDL(3,ix,iy,iz)
		glog= sqrt(glog)        
		Output(iy,iz) = glog		
        end do
      end do	
	case(2) !>>> Corte Y=cte <<<!	 
	  iy = yi(2)
      do iz=zi(2),zf(2)
        do ix=xi(2),xf(2)		
		glog = conjg(pInfoSDL(1,ix,iy,iz))*pInfoSDL(1,ix,iy,iz)
		glog = glog+conjg(pInfoSDL(2,ix,iy,iz))*pInfoSDL(2,ix,iy,iz)
		glog = glog+conjg(pInfoSDL(3,ix,iy,iz))*pInfoSDL(3,ix,iy,iz)
		glog= sqrt(abs(glog))
        Output(ix,iz) = glog 		
        end do
      end do
	case(3) !>>> Corte Z=cte <<<!
	  iz = zi(3)
      do iy=yi(3),yf(3)               
        do ix=xi(3),xf(3)
		glog = conjg(pInfoSDL(1,ix,iy,iz))*pInfoSDL(1,ix,iy,iz)
		glog = glog+conjg(pInfoSDL(2,ix,iy,iz))*pInfoSDL(2,ix,iy,iz)
		glog = glog+conjg(pInfoSDL(3,ix,iy,iz))*pInfoSDL(3,ix,iy,iz)
		glog= sqrt(glog)
		Output(ix,iy) = glog         
        end do
      end do   	
    end select 
case(-1) !charge density
 select case(LayerSDL)
    case(1) !>>> Corte X=cte <<<!
	  ix = xi(1)
      do iz=zi(1),zf(1)
        do iy=yi(1),yf(1)
		glog = pInfoSDL(1,ix,iy,iz)-pInfoSDL(1,ix-1,iy,iz)
		glog = glog+pInfoSDL(1,ix,iy,iz)-pInfoSDL(1,ix,iy-1,iz)
		glog = glog+pInfoSDL(1,ix,iy,iz)-pInfoSDL(1,ix,iy,iz-1)		
		Output(iy,iz) = glog		
        end do
      end do	
	case(2) !>>> Corte Y=cte <<<!	 
	  iy = yi(2)
      do iz=zi(2),zf(2)
        do ix=xi(2),xf(2)		
		glog = pInfoSDL(1,ix,iy,iz)-pInfoSDL(1,ix-1,iy,iz)
		glog = glog+pInfoSDL(1,ix,iy,iz)-pInfoSDL(1,ix,iy-1,iz)
		glog = glog+pInfoSDL(1,ix,iy,iz)-pInfoSDL(1,ix,iy,iz-1)		
        Output(ix,iz) = glog 		
        end do
      end do
	case(3) !>>> Corte Z=cte <<<!
	  iz = zi(3)
      do iy=yi(3),yf(3)               
        do ix=xi(3),xf(3)
		glog = pInfoSDL(1,ix,iy,iz)-pInfoSDL(1,ix-1,iy,iz)
		glog = glog+pInfoSDL(1,ix,iy,iz)-pInfoSDL(1,ix,iy-1,iz)
		glog = glog+pInfoSDL(1,ix,iy,iz)-pInfoSDL(1,ix,iy,iz-1)		
		Output(ix,iy) = glog         
        end do
      end do   	
    end select 
!end if
end select    

SDLnorm=1.0d0
if((em_source(0)==1))then !TFST or gaussian wave packet at t=0
SDLnorm=q*Epump!*(f_e*1e-5/q) field in E(kV/cm) or /(q*Epump) normalized to |Einc|
end if
    
If(SDLnorm==0.0) SDLnorm=1.0 !If no source is defined or there is somothing wrong

output=real(output)/SDLnorm 

!>>> Se gira la estructura 
if(SDL_rotate_180)then !>> 180º
   allocate(Output_aux(ixo:ixf,iyo:iyf))
   do ix=ixo,ixf
	 do iy=iyo,iyf
	    output_aux(ix,iy)    = output(ix,iy)	  
	 end do
	end do
		
	do ix=ixo,ixf
	 do iy=0,(iyf-iyo)	     
	     output(ix,iy+iyo) = output_aux(ix,iyf-iy)		 
	  end do 
	 end do
	 deallocate(output_aux)
end if   

if(SDL_rotate_90)then !>> 90º
   allocate(output_aux(iyo:iyf,ixo:ixf))
   do ix=ixo,ixf
	 do iy=iyo,iyf
	    output_aux(iy,ix)    = output(ix,iy)	  
	 end do
	end do
	!deallocate(output)
	!allocate(output(iyo:iyf,ixo:ixf))	
	!output=>output_aux	
	!deallocate(output_aux)
end if    
!<<< Se gira la estructura 

!>>> Varias copias
if(SDL_copy_flag)then 
if (LayerSDL==3)then !Array XY
    allocate(output_aux(ixo:SDL_copy_no*(ixf+1-ixo),iyo:SDL_copy_no*(iyf+1-iyo)))
    do i=0,(SDL_copy_no-1)
    do j=0,(SDL_copy_no-1)
     do ix=ixo,ixf
	   do iy=iyo,iyf	   
	    output_aux(ix+i*(ixf+1-ixo),iy+j*(iyf+1-iyo))    = output(ix,iy)	  
	   end do
	 end do
    end do
    end do
	!deallocate(output)
	!allocate(output(ixo:SDL_copy_no*(ixf+1-ixo),iyo:iyf))		
else !Array XZ,YZ 
   allocate(output_aux(ixo:SDL_copy_no*(ixf+1-ixo),iyo:iyf))
    do i=0,(SDL_copy_no-1)
     do ix=ixo,ixf
	   do iy=iyo,iyf	   
	    output_aux(ix+i*(ixf+1-ixo),iy)    = output(ix,iy)	  
	   end do
	 end do
    end do
	!deallocate(output)
	!allocate(output(ixo:SDL_copy_no*(ixf+1-ixo),iyo:iyf))		
end if
end if    
!<<< Varias copias


  if(SDL_scale)then
 !>>> Points to allocate
 ix=0;iy=0	
 do i=ixo,ixf,scale_factor              	   
   ix=ix+1 	       
 end do
 do j=iyo,iyf,scale_factor		   			  	 				
   iy=iy+1
 end do
 !<<< Points to allocate

 allocate(output_aux(1:ix,1:iy))
 ix=1
 iy=1	
 do i=ixo,ixf,scale_factor      
   do j=iyo,iyf,scale_factor		   			
  	 output_aux(ix,iy)=output(i,j)				
	 iy=iy+1
   end do
   iy=1		  	   
   ix=ix+1 	       
 end do
  end if

   !>>> Escala de colores
 Escala=0
if(Escala==1)then
  iz =1
  ix=1
  glog=0.0
  do iy=0,53
     !Output(:,iz:iz+10) = glog
    output(ix:ix+ixmax/54-1,:)=glog
	 glog=glog+1.d0/54.0	 
	 ix=ix+ixmax/54
  end do 
end if
!    Escala de colores

  
!maxem(em,ino) = maxval(output)
!output=output*(f_e*1e-5/q) !E(kV/cm)
!if(SDL_copy_flag) output_aux=output_aux*(f_e*1e-5/q) !E(kV/cm)  
  

 open (UNIT=4,FILE=Archivo,STATUS="REPLACE",FORM="UNFORMATTED",POSITION="REWIND",Action="write")		  
 
 
 if((SDL_copy_flag).or.(SDL_rotate_90).or.(SDL_scale).or.(SDL_rotate_180))then
 write (UNIT=4) real(output_aux)     
 else
 write (UNIT=4) real(output)
 end if


 !if(not(SDL_copy_flag))then
 !write (UNIT=4) real(output)
 !else
 !write (UNIT=4) real(output_aux)
 !end if
 !
 !if(not(SDL_rotate_90))then
 !write (UNIT=4) real(output)
 !else
 !write (UNIT=4) real(output_aux)
 !end if
 !
 !if(not(SDL_scale))then
 !write (UNIT=4) real(output)
 !else
 !write (UNIT=4) real(output_aux)
 !end if

close (UNIT=4)  

print*,trim(Archivo)  

!>>> Sutil: a Windows 7 no le gusta que un puntero "asociado" como => sea deallocatado 
!>>> Por eso la línea output=>output_aux está comentada y he modificado
!>>> el resto de líneas en consecuencia. La próxima vez que leas esto no sabrás de que
!>>> coño estás hablando (21abr11)
deallocate(output)   
!deallocate(output_aux)
if(SDL_copy_flag)  deallocate(output_aux)
if(SDL_rotate_90)  deallocate(output_aux)
if(SDL_scale)      deallocate(output_aux)

return
    end subroutine SDL_output	
    
!*****************************************************************************************+

    subroutine SDL_input(ino,NoComponente,em,pInfoSDL,maxem,LayerSDL)
use parameters
use physconsts
use SDL_MEMBERS 
use COMMON_SUBROUTINES
implicit none
  integer, intent(IN)::ino,NoComponente,em
  character(5) ::charno
  complex(prc), pointer::pInfoSDL(:,:,:,:)	
  real(prc),pointer::maxem(:,:)    
  integer ::LayerSDL
  integer ix,iy,iz,i,j  
  character(LONGITUD_ARCHIVO) Archivo  
  character(4) EXT
  real(4), pointer::Output(:,:)
  real(4), pointer::Output_aux(:,:)
  integer           ::ixo,ixf,iyo,iyf
  real(prc) :: glog
  integer:: output_file
  integer :: ixcc,izcc

  select case(LayerSDL)
    case(1) !>>> Corte X=cte (LayerSDL=1)<<<!
		ixo =yi(1)
		ixf =yf(1)
		iyo =zi(1)
		iyf =zf(1)
		allocate(Output(ixo:ixf,iyo:iyf))
	case(2) !>>> Corte Y=cte (LayerSDL=2)<<<!	
		ixo =xi(2)
		ixf =xf(2)
		iyo =zi(2)
		iyf =zf(2)
		allocate(Output(ixo:ixf,iyo:iyf))
	case(3) !>>> Corte Z=cte (LayerSDL=3)<<<!		
		ixo =xi(3)
		ixf =xf(3)
		iyo =yi(3)
		iyf =yf(3)
		allocate(Output(ixo:ixf,iyo:iyf))
  end select 
       
  !em=1----> Campo Eléctrico
  !em=2----> Campo Magnético
  if (em==1) then
	   !e1 = Ex
	  if(NoComponente.EQ.1) then
	   EXT ='.e1'
	  end if
	   !e2 = Ey
	  if(NoComponente.EQ.2) then
	   EXT ='.e2'
	  end if
	   !e3 = Ez
	  if(NoComponente.EQ.3) then
	   EXT = '.e3'
	  end if     
  else
  	   !h1 = Hx
	  if(NoComponente.EQ.1) then
	   EXT ='.h1'	   
	  end if
	   !h2 = Hy
	  if(NoComponente.EQ.2) then
	   EXT ='.h2'
	  end if
	   !h3 = Hz
	  if(NoComponente.EQ.3) then
	   EXT = '.h3'
	  end if  
  end if

 if(em==0)then !>>> Modulo campo eléctrico <<<!
   EXT ='.e'
 end if

 call inttoASCII(ino,charno)

 if(RUTA_SDL/='cluster')then
 select case (LayerSDL)
  case (1)
   Archivo='X_'//trim(RUTA_SDL)//trim(charno)//EXT 
  case (2)
   Archivo='Y_'//trim(RUTA_SDL)//trim(charno)//EXT 
  case (3)
   Archivo='Z_'//trim(RUTA_SDL)//trim(charno)//EXT 
 end select			
 else !>>> CLUSTER
   select case (LayerSDL)
  case (1)
   Archivo='./X_SDL'//trim(charno)//EXT 
  case (2)
   Archivo='./Y_SDL'//trim(charno)//EXT 
  case (3)
   Archivo='./Z_SDL'//trim(charno)//EXT 
 end select	
 end if	 	
 
 !>>> Se gira la estructura 
if(SDL_rotate_180)then !>> 180º
   allocate(Output_aux(ixo:ixf,iyo:iyf))
   do ix=ixo,ixf
	 do iy=iyo,iyf
	    output_aux(ix,iy)    = output(ix,iy)	  
	 end do
	end do
		
	do ix=ixo,ixf
	 do iy=0,(iyf-iyo)	     
	     output(ix,iy+iyo) = output_aux(ix,iyf-iy)		 
	  end do 
	 end do
	 deallocate(output_aux)
end if   

if(SDL_rotate_90)then !>> 90º
   allocate(output_aux(iyo:iyf,ixo:ixf))
   do ix=ixo,ixf
	 do iy=iyo,iyf
	    output_aux(iy,ix)    = output(ix,iy)	  
	 end do
	end do
	!deallocate(output)
	!allocate(output(iyo:iyf,ixo:ixf))	
	!output=>output_aux	
	!deallocate(output_aux)
end if    
!<<< Se gira la estructura 

!>>> Varias copias
if(SDL_copy_flag)then 
   allocate(output_aux(ixo:SDL_copy_no*(ixf+1-ixo),iyo:iyf))
    do i=0,(SDL_copy_no-1)
     do ix=ixo,ixf
	   do iy=iyo,iyf	   
	    output_aux(ix+i*(ixf+1-ixo),iy)    = output(ix,iy)	  
	   end do
	 end do
	end do
	!deallocate(output)
	!allocate(output(ixo:SDL_copy_no*(ixf+1-ixo),iyo:iyf))		
end if    
!<<< Varias copias

if(SDL_scale)then
 !>>> Points to allocate
 ix=0;iy=0	
 do i=ixo,ixf,scale_factor              	   
   ix=ix+1 	       
 end do
 do j=iyo,iyf,scale_factor		   			  	 				
   iy=iy+1
 end do
 !<<< Points to allocate

 allocate(output_aux(1:ix,1:iy))
 ix=1
 iy=1	
 do i=ixo,ixf,scale_factor      
   do j=iyo,iyf,scale_factor		   			
  	 output_aux(ix,iy)=output(i,j)				
	 iy=iy+1
   end do
   iy=1		  	   
   ix=ix+1 	       
 end do
end if

!maxem(em,ino) = maxval(output)
!output=output*(f_e*1e-5/q) !E(kV/cm)
!if(SDL_copy_flag) output_aux=output_aux*(f_e*1e-5/q) !E(kV/cm)

 open (UNIT=4,FILE=Archivo,FORM="UNFORMATTED",POSITION="REWIND",Action="read")		  
 
 !WARNING: output and output_aux has forced to be a REAL 4bytes to read. In SDL_output data (output/output_aux) is saved as with 4bytes.
 !if(not(SDL_copy_flag))then
 if((SDL_rotate_90==0).and.(SDL_rotate_90==0).and.(SDL_scale==0))then
 read (UNIT=4) output
 else
 read (UNIT=4) output_aux
 end if
 !
 !if(not(SDL_rotate_90))then
 !read (UNIT=4) output
 !else
 !read (UNIT=4) output_aux
 !end if
 !
 !if(not(SDL_scale))then
 !read (UNIT=4) output
 !else
 !read (UNIT=4) output_aux
 !end if
 close (UNIT=4) 
print*,trim(Archivo)  
 
!>>> OUTPUT


if(.false.)then
if(RUTA_SDL/='cluster')then
 select case (LayerSDL)
  case (1)
   Archivo='MX_'//trim(RUTA_SDL)//trim(charno)//EXT 
  case (2)
   Archivo='MY_'//trim(RUTA_SDL)//trim(charno)//EXT 
  case (3)
   Archivo='MZ_'//trim(RUTA_SDL)//trim(charno)//EXT 
 end select			
 else !>>> CLUSTER
   select case (LayerSDL)
  case (1)
   Archivo='./MX_SDL'//trim(charno)//EXT 
  case (2)
   Archivo='./MY_SDL'//trim(charno)//EXT 
  case (3)
   Archivo='./MZ_SDL'//trim(charno)//EXT 
 end select	
 end if	 	
 
  open (UNIT=4,FILE=Archivo,STATUS="REPLACE",FORM="UNFORMATTED",POSITION="REWIND",Action="write")	  
 
 
 if((SDL_copy_flag).or.(SDL_rotate_90).or.(SDL_scale).or.(SDL_rotate_180))then
 output_aux=-output_aux 
 write (UNIT=4) real(output_aux)     
 else
 output=-output

 !do ix=ixo,ixf
 !do iy=iyo,iyf	   
	!    output(ix,iy)    = (output(ix,iy)+output(ixmax-(ix-1),iy))/2.0	  !SIMETRIC PART XZ plane
 !end do
 !end do
 
 write (UNIT=4) real(output)
 end if

 close (UNIT=4) 
 else
  !em=1----> Campo Eléctrico
  !em=2----> Campo Magnético

    if(em==0)then !>>> Modulo campo eléctrico <<<!
       output_file = 100
    end if

  if (em==1) then
	   !e1 = Ex
	  if(NoComponente.EQ.1) then
	   output_file=102
	  end if
	   !e2 = Ey
	  if(NoComponente.EQ.2) then
	   output_file=104
	  end if
	   !e3 = Ez
	  if(NoComponente.EQ.3) then
	   output_file=105
	  end if     
  else
  	   !h1 = Hx
	  if(NoComponente.EQ.1) then
	   output_file=106 
	  end if
	   !h2 = Hy
	  if(NoComponente.EQ.2) then
	   output_file=108
	  end if
	   !h3 = Hz
	  if(NoComponente.EQ.3) then
	   output_file=110
	  end if  
  end if
 
 select case(LayerSDL)
    case(1) !>>> Corte X=cte <<<!
	  ix = xi(1)
      do iz=zi(1),zf(1)
        do iy=yi(1),yf(1)       				
		!do something
        end do
      end do	
	case(2) !>>> Corte Y=cte <<<!	 
	    iy = yi(2)        
        print*,"Fixed iy=",iy,"..readig ixo,izo for cross cuts [if (iz<0 --> 2D map else X Y Z)]"; read*,ixcc,izcc
        if(izcc>=0)then
            write(output_file,*) '#',ino
            do ix=xi(2),xf(2)		  
            write(output_file,*) output(ix,izcc)
            end do      
            write(output_file+1,*) '#',ino
            do iz=zi(2),zf(2)		  
            write(output_file+1,*) output(ixcc,iz)
            end do      
        else
            !Recovers the whole SDL map in text format
            do iz=zi(2),zf(2)
            write(121,*) '#',iz
            do ix=xi(2),xf(2)		  
            write(121,*) output(ix,iz)
            end do
            end do
      end if
	case(3) !>>> Corte Z=cte <<<!
	  iz = zi(3)
      do iy=yi(3),yf(3)               
        do ix=xi(3),xf(3)						
		!do something
        end do
      end do   	
    end select   
 
!>>> Sutil: a Windows 7 no le gusta que un puntero "asociado" como => sea deallocatado 
!>>> Por eso la línea output=>output_aux está comentada y he modificado
!>>> el resto de líneas en consecuencia. La próxima vez que leas esto no sabrás de que
!>>> coño estás hablando (21abr11)
 end if
 
deallocate(output)   
if(SDL_copy_flag)  deallocate(output_aux)
if(SDL_rotate_90)  deallocate(output_aux)
if(SDL_scale)      deallocate(output_aux)

return
    end subroutine SDL_input

    !>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<>>><<<!

module SDL_MESSAGES
 interface 
    subroutine setparam_SDL(inputSDL_FILE,log_FILE,SDL_name)
	  use SDL_MEMBERS	  
      integer,intent(in)::inputSDL_FILE
	  integer,intent(in)::log_FILE
      character(LONGITUD_ARCHIVO),intent(in)::SDL_name      
    end subroutine setparam_SDL

	subroutine setSDL_FILE(e,h,maxem,nt,LayerSDL)
     use SDL_MEMBERS
	 use SDL_PRIVATE
	 implicit none
	 complex(prc), pointer::e(:,:,:,:),h(:,:,:,:)
	 real(prc),pointer::maxem(:,:)
	 integer,intent(in)::nt
	 integer           :: LayerSDL
    end subroutine setSDL_File

    subroutine fromSDL_FILE(e,h,maxem,nt,LayerSDL)
     use SDL_MEMBERS
	 use SDL_PRIVATE
	 implicit none
	 complex(prc), pointer::e(:,:,:,:),h(:,:,:,:)
	 real(prc),pointer::maxem(:,:)
	 integer,intent(in)::nt
	 integer           :: LayerSDL
    end subroutine fromSDL_File
    
 end interface
end module SDL_MESSAGES   

! **************************************************** !

subroutine setparam_SDL(inputSDL_FILE,log_FILE,SDL_name) 
 use SDL_MEMBERS
  implicit none
 integer,intent(in)::inputSDL_FILE
 integer,intent(in)::log_FILE
 character(LONGITUD_ARCHIVO),intent(in)::SDL_name

  character(LONGITUD_ARCHIVO)::message
  message="" !>>> Tal vez lo hayas utilizado antes y por el motivo que sea queda en memoria <<<!
  rewind(unit=inputSDL_FILE)  
  do while(message/=SDL_name)
    read(inputSDL_FILE,*) message
  end do  
    
  read(inputSDL_FILE,*) RUTA_SDL
  if (RUTA_SDL=='') call err(13)
  write(log_FILE,'(A1,A12,A80)') '#','SDL path=',RUTA_SDL 
  write(log_FILE,'(A1,A28,I10,A6)') '#','Archivos ha generar..=',NO_FILES,' pasos'
  read(inputSDL_FILE,*) SDL_copy_flag,SDL_copy_no,SDL_rotate_90,SDL_rotate_180 !,SDL_scale,scale_factor
  read(inputSDL_FILE,*) xi(1),xi(2),xi(3)
  read(inputSDL_FILE,*) xf(1),xf(2),xf(3)
  read(inputSDL_FILE,*) yi(1),yi(2),yi(3)
  read(inputSDL_FILE,*) yf(1),yf(2),yf(3)
  read(inputSDL_FILE,*) zi(1),zi(2),zi(3)
  read(inputSDL_FILE,*) zf(1),zf(2),zf(3)
  read(inputSDL_FILE,*) SDLModuleE
  read(inputSDL_FILE,*) SDLField(1)
  read(inputSDL_FILE,*) SDLField(2)
  read(inputSDL_FILE,*) SDLField(3)
  read(inputSDL_FILE,*) SDLField(4)
  read(inputSDL_FILE,*) SDLField(5)
  read(inputSDL_FILE,*) SDLField(6)     
end subroutine setparam_SDL

!******************************************************************!

subroutine setSDL_FILE(e,h,maxem,nt,LayerSDL)
use SDL_MEMBERS
use SDL_PRIVATE
implicit none

complex(prc), pointer :: e(:,:,:,:)
complex(prc), pointer :: h(:,:,:,:)
real(prc),pointer     :: maxem(:,:)
integer,intent(in)    :: nt
integer               :: LayerSDL
! SDL_output(ino,NoComponente,em,pInfoSDL,maxem,LayerSDL)
    if(SDLModuleE==1)then !Electric field
		call SDL_output(nt,0,0,e,maxem,LayerSDL)    
    end if
    if(SDLModuleE==2)then !Magnetic field
		call SDL_output(nt,0,0,h,maxem,LayerSDL)    
    end if
    if(SDLField(1)>=1)then  !=1 Re(); =2 Im(); =3 ||
		call SDL_output(nt,1,1,e,maxem,LayerSDL)    
	end if
	if(SDLField(2)>=1)then	
		call SDL_output(nt,2,1,e,maxem,LayerSDL)
    end if
	if(SDLField(3)>=1)then	
		call SDL_output(nt,3,1,e,maxem,LayerSDL)
    end if  
    
    !!Warning (17jul19): H-field is not used too much so I used this variable for img(E),charge density...
	if(SDLField(4)>=1)then		  !=1 Re(); =2 Im(); =3 ||
        call SDL_output(nt,1,2,h,maxem,LayerSDL)  
    end if
	if(SDLField(5)>=1)then	    
		call SDL_output(nt,2,2,h,maxem,LayerSDL)         
    end if
	if(SDLField(6)>=1)then	
    	call SDL_output(nt,3,2,h,maxem,LayerSDL)     
    
    if(.false.)then
    call SDL_output(nt,0,-1,e,maxem,LayerSDL)!For charge density calculations    
    Print*,"Charge density (a.u.)"
    end if
    end if

return
    end subroutine setSDL_FILE

!******************************************************************!

    subroutine fromSDL_FILE(e,h,maxem,nt,LayerSDL)
use SDL_MEMBERS
use SDL_PRIVATE
implicit none

complex(prc), pointer :: e(:,:,:,:)
complex(prc), pointer :: h(:,:,:,:)
real(prc),pointer     :: maxem(:,:)
integer,intent(in)    :: nt
integer               :: LayerSDL
    if(SDLModuleE==1)then	
		call SDL_input(nt,0,0,e,maxem,LayerSDL)    
	end if
    if(SDLField(1)==1)then	
		call SDL_input(nt,1,1,e,maxem,LayerSDL)    
	end if
	if(SDLField(2)==1)then	
		call SDL_input(nt,2,1,e,maxem,LayerSDL)
    end if
	if(SDLField(3)==1)then	
		call SDL_input(nt,3,1,e,maxem,LayerSDL)
    end if
	if(SDLField(4)==1)then	
		call SDL_input(nt,1,2,h,maxem,LayerSDL)   
    end if
	if(SDLField(5)==1)then	
		call SDL_input(nt,2,2,h,maxem,LayerSDL) 
    end if
	if(SDLField(6)==1)then	
		call SDL_input(nt,3,2,h,maxem,LayerSDL) 
	end if
return
end subroutine fromSDL_FILE

!******************************************************************!

