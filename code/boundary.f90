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

module	BOUNDARY
interface

subroutine bc_xmin_bloch(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_xmin_bloch

subroutine bc_xmax_bloch(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_xmax_bloch

subroutine bc_xmin_waveguide(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_xmin_waveguide

subroutine bc_xmax_waveguide(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_xmax_waveguide

subroutine bc_ymin_waveguide(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_ymin_waveguide

subroutine bc_ymax_waveguide(e,h)
use physconsts
use parameters
implicit none
complex(prc),pointer :: e(:,:,:,:)
complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_ymax_waveguide

subroutine bc_zmin_waveguide(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_zmin_waveguide

subroutine bc_zmax_waveguide(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_zmax_waveguide

subroutine bc_ymin_bloch(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_ymin_bloch

subroutine bc_ymax_bloch(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_ymax_bloch

subroutine bc_zmin_bloch(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_zmin_bloch

subroutine bc_zmax_bloch(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_zmax_bloch


subroutine bc_xmin_metal(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_xmin_metal 

subroutine bc_xmax_metal(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_xmax_metal 

subroutine bc_ymin_metal(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_ymin_metal 

subroutine bc_ymax_metal(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_ymax_metal 

subroutine bc_zmin_metal(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_zmin_metal 

subroutine bc_zmax_metal(e,h)
use parameters
	implicit none
	complex(prc),pointer :: e(:,:,:,:)
	complex(prc),pointer :: h(:,:,:,:)
end subroutine bc_zmax_metal 

end interface
end module BOUNDARY

! -----------------------------------------------------------------------------
! Subroutine bc_xmin_bloch
! Set the boundary conditions at ix=1
! -----------------------------------------------------------------------------
subroutine bc_xmin_bloch(e,h)
use physconsts
use parameters
implicit none

complex(prc),pointer :: e(:,:,:,:)
complex(prc),pointer :: h(:,:,:,:)

integer :: iy,iz,i
complex(precision) :: fxm1

! Set the Bloch phase shift
fxm1=cdexp(-ci*akx)
! Loop over the surface at ix=1
do iy=1,iymax
   do iz=1,izmax
      do i=1,3
         e(i,0,iy,iz)=fxm1*e(i,ixmax,iy,iz)
         h(i,0,iy,iz)=fxm1*h(i,ixmax,iy,iz)
      enddo
   enddo
enddo

return
end subroutine bc_xmin_bloch

! -----------------------------------------------------------------------------
! Subroutine bc_xmax_bloch
! Set the boundary conditions at ix=ixmax
! -----------------------------------------------------------------------------
subroutine bc_xmax_bloch(e,h)
use physconsts
use parameters
implicit none

complex(prc),pointer :: e(:,:,:,:)
complex(prc),pointer :: h(:,:,:,:)

integer :: iy,iz,i
complex(precision) :: fxp1

! Set the Bloch phase shift
fxp1=cdexp(ci*akx)
! Loop over the surface at ix=ixmax
do iy=1,iymax
   do iz=1,izmax
      do i=1,3
         e(i,ixmax+1,iy,iz)=fxp1*e(i,1,iy,iz)
         h(i,ixmax+1,iy,iz)=fxp1*h(i,1,iy,iz)
      enddo
   enddo
enddo
return
    end subroutine bc_xmax_bloch

!********************!

    ! -----------------------------------------------------------------------------
! Subroutine bc_ymin_bloch
! Set the boundary conditions at iy=1
! -----------------------------------------------------------------------------
subroutine bc_ymin_bloch(e,h)
use physconsts
use parameters
implicit none

complex(prc),pointer :: e(:,:,:,:)
complex(prc),pointer :: h(:,:,:,:)

integer :: ix,iz,i
complex(precision) :: fym1

! Set the Bloch phase shift

fym1=cdexp(-ci*aky)

! Loop over the surface at iy=1

do ix=1,ixmax
   do iz=1,izmax
      do i=1,3
         e(i,ix,0,iz)=fym1*e(i,ix,iymax,iz)
         h(i,ix,0,iz)=fym1*h(i,ix,iymax,iz)
      enddo
   enddo
enddo

return
end subroutine bc_ymin_bloch

! -----------------------------------------------------------------------------
! Subroutine bc_ymax_bloch
! Set the boundary conditions at iy=iymax
! -----------------------------------------------------------------------------
subroutine bc_ymax_bloch(e,h)
use physconsts
use parameters
implicit none

complex(prc),pointer :: e(:,:,:,:)
complex(prc),pointer :: h(:,:,:,:)

integer :: ix,iz,i
complex(precision) :: fyp1

! Set the Bloch phase shift

fyp1=cdexp(ci*aky)

! Loop over the surface at iy=iymax

do ix=1,ixmax
   do iz=1,izmax
      do i=1,3
         e(i,ix,iymax+1,iz)=fyp1*e(i,ix,1,iz)
         h(i,ix,iymax+1,iz)=fyp1*h(i,ix,1,iz)
      enddo
   enddo
enddo

return
    end subroutine bc_ymax_bloch
    
!*******************************************************

subroutine bc_zmin_bloch(e,h)
use physconsts
use parameters
implicit none

complex(prc),pointer :: e(:,:,:,:)
complex(prc),pointer :: h(:,:,:,:)

integer :: ix,iy,i
complex(precision) :: fym1

! Set the Bloch phase shift

fym1=cdexp(-ci*akz)

! Loop over the surface at iz=1

do ix=1,ixmax
   do iy=1,iymax
      do i=1,3
         e(i,ix,iy,0)=fym1*e(i,ix,iy,izmax)
         h(i,ix,iy,0)=fym1*h(i,ix,iy,izmax)
      enddo
   enddo
enddo

return
end subroutine bc_zmin_bloch

! -----------------------------------------------------------------------------
! Subroutine bc_ymax_bloch
! Set the boundary conditions at iy=iymax
! -----------------------------------------------------------------------------
subroutine bc_zmax_bloch(e,h)
use physconsts
use parameters
implicit none

complex(prc),pointer :: e(:,:,:,:)
complex(prc),pointer :: h(:,:,:,:)

integer :: ix,iy,i
complex(precision) :: fyp1

! Set the Bloch phase shift

fyp1=cdexp(ci*akz)

! Loop over the surface at iy=iymax

do ix=1,ixmax
   do iy=1,iymax
      do i=1,3
         e(i,ix,iy,izmax+1)=fyp1*e(i,ix,iy,1)
         h(i,ix,iy,izmax+1)=fyp1*h(i,ix,iy,1)
      enddo
   enddo
enddo

return
end subroutine bc_zmax_bloch
! -----------------------------------------------------------------------------
! Subroutine bc_xmin_metal
! Set the boundary conditions at ix=1
! -----------------------------------------------------------------------------
subroutine bc_xmin_metal(e,h)
use parameters
implicit none

complex(prc),pointer :: e(:,:,:,:)
complex(prc),pointer :: h(:,:,:,:)

integer :: iy,iz,i

! Loop over the surface at ix=1

do iy=1,iymax
   do iz=1,izmax
      do i=1,3
         e(i,0,iy,iz)=(0.0,0.0)
         h(i,0,iy,iz)=(0.0,0.0)
      enddo
   enddo
enddo

return
end subroutine bc_xmin_metal

! -----------------------------------------------------------------------------
! Subroutine bc_xmax_metal
! Set the boundary conditions at ix=ixmax
! -----------------------------------------------------------------------------
subroutine bc_xmax_metal(e,h)
use parameters
implicit none

complex(prc),pointer :: e(:,:,:,:)
complex(prc),pointer :: h(:,:,:,:)

integer :: iy,iz,i

! Loop over the surface at ix=ixmax

do iy=1,iymax
   do iz=1,izmax
      do i=1,3
         e(i,ixmax+1,iy,iz)=(0.0,0.0)
         h(i,ixmax+1,iy,iz)=(0.0,0.0)
      enddo
   enddo
enddo

return
end subroutine bc_xmax_metal

! -----------------------------------------------------------------------------
! Subroutine bc_ymin_metal
! Set the boundary conditions at iy=1
! -----------------------------------------------------------------------------
subroutine bc_ymin_metal(e,h)
use parameters
implicit none

complex(prc),pointer :: e(:,:,:,:)
complex(prc),pointer :: h(:,:,:,:)

integer :: ix,iz,i

! Loop over the surface at iy=1

do ix=1,ixmax
   do iz=1,izmax
      do i=1,3
         e(i,ix,0,iz)=(0.0,0.0)
         h(i,ix,0,iz)=(0.0,0.0)
      enddo
   enddo
enddo

return
end subroutine bc_ymin_metal

! -----------------------------------------------------------------------------
! Subroutine bc_ymax_metal
! Set the boundary conditions at iy=iymax
! -----------------------------------------------------------------------------
subroutine bc_ymax_metal(e,h)
use parameters
implicit none

complex(prc),pointer :: e(:,:,:,:)
complex(prc),pointer :: h(:,:,:,:)

integer :: ix,iz,i

! Loop over the surface at iy=iymax

do ix=1,ixmax
   do iz=1,izmax
      do i=1,3
         e(i,ix,iymax+1,iz)=(0.0,0.0)
         h(i,ix,iymax+1,iz)=(0.0,0.0)
      enddo
   enddo
enddo

return
end subroutine bc_ymax_metal
! -----------------------------------------------------------------------------
! Subroutine bc_zmin_metal
! Set the boundary conditions at iz=1
! -----------------------------------------------------------------------------
subroutine bc_zmin_metal(e,h)
use parameters
implicit none

complex(prc),pointer :: e(:,:,:,:)
complex(prc),pointer :: h(:,:,:,:)

integer :: ix,iy,i

! Loop over the surface at iz=1

do ix=1,ixmax
   do iy=1,iymax
      do i=1,3
         e(i,ix,iy,0)=(0.0,0.0)
         h(i,ix,iy,0)=(0.0,0.0)
      enddo
   enddo
enddo

return
end subroutine bc_zmin_metal

! -----------------------------------------------------------------------------
! Subroutine bc_zmax_metal
! Set the boundary conditions at iz=izmax
! -----------------------------------------------------------------------------
subroutine bc_zmax_metal(e,h)
use parameters
implicit none

complex(prc),pointer :: e(:,:,:,:)
complex(prc),pointer :: h(:,:,:,:)

integer :: ix,iy,i

! Loop over the surface at iz=izmax

do ix=1,ixmax
   do iy=1,iymax
      do i=1,3
         e(i,ix,iy,izmax+1)=(0.0,0.0)
         h(i,ix,iy,izmax+1)=(0.0,0.0)
      enddo
   enddo
enddo
return
end subroutine bc_zmax_metal
