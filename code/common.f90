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

    module COMMON_SUBROUTINES
interface

subroutine draw_field(field,index)
	use parameters
	use files	
	implicit none
    complex(prc),pointer                    :: field(:,:,:,:,:)    
    character :: index
end subroutine

end interface
contains 
! ***************************************************************************************** !

subroutine inttoASCII(ino,charno)
use files
 implicit none

 integer,intent(in)::ino
 character(5),intent(out)::charno
 integer::sCase,U,D,C,M,DM,iInter  

if (ino.LE.9) then
     sCase=0  	 
  end if

  if ((ino.GT.9).AND.(ino.LE.99)) then
	   sCase=1
  end if

  if (ino.GT.99) then
     sCase=2 	 
  end if  

  if (ino.GT.999) then
     sCase=3 	 
   end if  
   
   if (ino.GT.9999) then
     sCase=4 	 
   end if    

	select case(sCase)
	   case(0)
		    U=ino; charno=ACHAR(INT(U+48))
	   case(1)
		   D=ino/10
	       U=Mod(ino,10)
           charno=ACHAR(INT(D+48))//ACHAR(INT(U+48))
	   case(2)
		   C=ino/100
	       U=Mod(ino,10)
	  	   D=(ino-(100*C+U))/10			
           charno=ACHAR(INT(C+48))
		   charno=trim(charno)//ACHAR(INT(D+48))//ACHAR(INT((U+48)))		
	   case(3)
		   M=ino/1000
		   iInter = ino-1000*M
		   C=iInter/100
		   U=Mod(iInter,10)
		   D=(iInter-(100*C+U))/10				
		   charno=ACHAR(INT(M+48))//ACHAR(INT(C+48))//ACHAR(INT(D+48))//ACHAR(INT((U+48)))  
       case(4)
           DM=ino/10000
           iInter=ino-10000*DM
		   M=iInter/1000
           iInter=ino-10000*DM - 1000*M
		   C=iInter/100
		   U=Mod(iInter,10)
		   D=(iInter-(100*C+U))/10				
		   charno=ACHAR(INT(DM+48))//ACHAR(INT(M+48))//ACHAR(INT(C+48))//ACHAR(INT(D+48))//ACHAR(INT((U+48)))  
	   case default
		    write(logfile,*) "...Error - Subroutine:inttoASCII"
	end select   
return
end subroutine inttoASCII

!************************************************************************************************!
subroutine spectra_lossless_dielectric()
use parameters
use physconsts 
use files
	implicit none  
    integer :: iw
	real(precision) :: w,lp	
	complex(precision)  ::nII,nIII
	real(precision)  :: del
	complex(precision) ::cn,cnp	
	real(precision)  :: R,T

 ! >>> Lossless dielectric film of thickness Dh <<< !      
   do iw=nwo,nwf
		lp=lno+dlamda*(iw-1.0)	! in nm
		nII = sqrt(eps_bbox)
		nIII= sqrt(eps(3))
		del = 2*pi*nII*Dh*Q/(lp*f_abohr)
		
		cn = (dcos(del)+ci*(nIII/nII)*dsin(del) - nIII*dcos(del)-ci*nII*dsin(del))
		cn = cn/(dcos(del)+ci*(nIII/nII)*dsin(del) + nIII*dcos(del)+ci*nII*dsin(del)) 				
		R = cdabs(dconjg(cn)*cn)
        T = 1.0d0 - R       
		write(spectra,*) lp,R,T
	end do
end subroutine

end module COMMON_SUBROUTINES

!*************************************************!

subroutine draw_field(field,index)
use parameters
use files
implicit none
  complex(prc),pointer                    :: field(:,:,:,:,:)  
  real(prc),pointer::maxem(:,:)
  integer::ilayer
  integer i,j,k,r
  integer :: graphic,nolayer,il  
  logical :: flag
  character :: index
  integer       :: mxi2,mxf2,myi2,myf2,mzi2,mzf2 !To solve BUG: SERS input override by mesh_graph
  
  msg=''
  rewind(unit=inFDTD)  
  do while(msg/=msg2)
     read(inFDTD,*) msg
  end do
  
  read(inFDTD,*) graphic

 if(graphic)then
  open(unit=meshgraph,file='draw_field'//index//'.dat',status='replace',action='write')
  read(inFDTD,*) LayerMESH
  mxi2=mxi;mxf2=mxf;myi2=myi;myf2=myf;mzi2=mzi;mzf2=mzf!To solve BUG: SERS input override by mesh_graph
  read(inFDTD,*) mxi,mxf,myi,myf,mzi,mzf
  read(inFDTD,*) nolayer    

  if(nolayer/=0)then   
   il = 1
   do while(il<=nolayer)
     read(inFDTD,*) k
     write(meshgraph,*) '#',k
	   select case (LayerMESH)
	     case (1) ! >>> Plano X=cte <<< 
		   do i=myi,myf
			 write(meshgraph,*) '#'!,i
             do j=mzi,mzf			   
		       write(meshgraph,*) real(field(1,1,k,i,j))			   
			 end do 
           end do
	      case (2) ! >>> Plano Y=cte <<<	
		   !do i=1,ixmax		    
		     do i=mxi,mxf
			 write(meshgraph,*) '#'!,i
             do j=mzi,mzf
		       write(meshgraph,*) real(field(1,1,i,k,j))	
			   if((i==50).and.(j==300).and.(k==50)) print*,real(field(1,1,i,k,j))		   			   
			 end do 
           end do       
	      case (3) ! >>> Plano Z=cte <<<
		   do i=mxi,mxf
			 write(meshgraph,*) '#'!,i
             do j=myi,myf
		       write(meshgraph,*) real(field(1,1,i,j,k))			   			   
			 end do 
           end do     
	   end select	                  	   			
     il = il+1
	end do  
  end if
write(logfile,*) "...mesh-graph generated"
print*,          "...mesh-graph generated"	  
close(unit=meshgraph)
mxi=mxi2;mxf=mxf2;myi=myi2;myf=myf2;mzi=mzi2;mzf=mzf2!To solve BUG: SERS input override by mesh_graph
end if
return           	  
end subroutine draw_field


