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

module physconsts
use parameters
! Use Atomic units
! c0=Speed of light
! pi=3.1415...
! emach=Machine accuracy
real(precision),parameter :: pi       = 3.14159265358979323846
real(precision),parameter :: c0       = 137.0360
real(precision),parameter :: eps0     = 1.0/(4.0*pi)
real(precision),parameter :: mu0      = 1.0/(eps0*(c0**2))
real(precision),parameter :: mu_inv   = 1.0
real(precision),parameter :: f_abohr  = 10/0.529177249 !Factor conversión nm -> unidades atómicas
real(precision),parameter :: x_A      = 0.0 !x_Accuracy
real(precision),parameter :: emach    = epsilon(x_A)
complex(precision),parameter :: ci    = (0.0,1.0)

real(precision),parameter :: uat      = 24.188843265E-3 !1uat -> femto seg.(1E-15 seg.)
real(precision) :: hartree              = 27.211396 !1Hartree -> 27.211396 eV
real(precision) :: heV                  = 6.626E-34       !J*s
real(precision),parameter :: ceV        = 2.99792456E8             !m/s
real(precision) :: eeV                  = 1.6021773E-19   !J/eV
real(precision) :: f_e                  = 514.2208281e9   !1 u.a. de campo eléctrico -> 514,2208281 GV/m (1G=1E9)
real(precision) :: f_h                  = 2.350518e5  !1 u.a. de campo magnético -> 2.350518e5 
real(precision) :: eVnm                 = 1239.828330d0
real(precision),parameter :: mu0SI    = 4*pi*10E-7
real(precision),parameter :: eps0SI   = 1.0/(ceV**2*mu0SI)
real(precision),parameter :: uapower  = 0.180237814  
real(precision),parameter :: charge_conv    = 1.60217733e-19 
real(precision),parameter :: velocity_conv  = 2.18769142e6 


double precision,parameter :: c0_mks       = 299792458.0d0      !m/s
double precision,parameter :: eps0_mks     = 8.854187817e-12    !F·m-1 
double precision,parameter :: mu0_mks      = 1.0/(eps0*(c0**2)) !N·A-2


!***************************************!
! 1 u.a.Energía  -->  4,35974417×10-18     J
! 1 u.a.Tiempo   -->  2,418884326505×10-17 s
! 1 u.a.Campo E  -->  514.2208281x10^9    V/m
! 1 u.a.Potencia ---> 1.80237814x10-1   Watts
! 1 u.a.<S>      ---> 6.43640931x10^19 W/m^2

! 1 m ---> 10E10/0.529177249     (u.a.l.)
! 1 s ---> 10E18/24.188843265    (u.a.t.)
! 1 V/m--> (1/514.2208281)*10^(-9) (u.a. de campo eléctrico)
! 1 V/m--> 10E-4/3               (statvolt/cm)
! 1 u.a.energía --> 4,35974417×10-18 J
! 1 u.a.tiempo  --> 2,418884326505×10-17 s
! 1 u.a.campo E --> 514.2208281)x10^9 V/m
!***************************************!
end module physconsts
