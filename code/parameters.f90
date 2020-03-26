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

module parameters


!>>> Precison=4 ->4bytes; 8 -> 8Bytes
integer,parameter ::prc = 8        ! De variables que CONSUMEN MEMORIA
integer,parameter ::precision = 8  ! De variables auxiliares que no consumen memoria
!    Precison=4 ->4bytes; 8 -> 8Bytes <<<

type  INTEGER_POINTER_ARRAY
      integer, dimension(:,:), pointer            :: p
end type INTEGER_POINTER_ARRAY

type REAL_POINTER_ARRAY
      real(prc), dimension(:,:), pointer          :: p
end type REAL_POINTER_ARRAY  

type COMPLEX_POINTER_ARRAY
      complex(prc), dimension(:,:), pointer       :: ks ! K y sigma
end type COMPLEX_POINTER_ARRAY

type COMPLEX_POINTER_ARRAY_II
      complex(prc), dimension(:,:,:,:), pointer   :: pml 
end type COMPLEX_POINTER_ARRAY_II

type COMPLEX_POINTER_ARRAY_III
      complex(prc), dimension(:,:,:,:,:), pointer :: em 
end type COMPLEX_POINTER_ARRAY_III

type COMPLEX_POINTER_ARRAY_IV
      complex(prc), dimension(:,:,:,:), pointer   :: e
end type COMPLEX_POINTER_ARRAY_IV

type COMPLEX_POINTER_ARRAY_VI
      complex(prc), dimension(:), pointer         :: fdata(:,:)
end type COMPLEX_POINTER_ARRAY_VI

type COMPLEX_POINTER_ARRAY_VII
      integer, dimension(:), pointer              :: pto(:,:)
end type COMPLEX_POINTER_ARRAY_VII


type EMR
  integer::i
  integer::j
  integer::k
  complex(precision)::ex
  complex(precision)::ey
  complex(precision)::ez    
end type

type currentFIELD  
  complex(precision)::ex
  complex(precision)::ey
  complex(precision)::hx
  complex(precision)::hy
end type

integer            :: ixmax,iymax,izmax,itmax,n_pts_store,cluster
real(precision)    :: Q,dt,dL !>>> a == periodicidad <<<
!real(precision) :: akx,aky,akz
complex(precision) :: akx,aky,akz
integer            :: iz_cur
integer            :: iy_cur
integer            :: ix_cur
integer            :: nno,nno_indc
integer, pointer   :: xd(:),yd(:),zd(:),xd_indc(:),yd_indc(:),zd_indc(:)
real(precision),    pointer   :: ois(:,:),ois_indc(:,:) !>>> os=orientation in space
integer            :: met_with_holes
integer         :: s_uno,s_dos,s_aux

real(precision) :: n(2) ! >>> Indice de refracción de la capa de medida s_uno y s_dos respectivamente <<<
integer         :: zo     
real(precision) :: GAP,GAP2

integer         :: nopto
integer         :: i_pol
integer         :: SDL_ON,SDL_ON_FOURIER,SDL_plane_on(3)
integer         :: SDL_xi1,SDL_yi2,SDL_zi3,SDL_Xextralayers,SDL_Yextralayers,SDL_Zextralayers
integer, pointer:: SDL_deltaX(:),SDL_deltaY(:),SDL_deltaZ(:)
real(precision) :: CFL
real(precision) :: lf,lno,lnf,lno_aux,lnf_aux,lf2,lno_bs !_bs = band structure (or dispersion relation)
integer         :: nwo   
integer         :: nwf  
real(precision) :: wo
real(precision) :: wf

integer         :: em_source(0:1)
integer         :: system_model !>>> 0:Single object 1:Array <<< !
logical         :: single_aperture 
integer         :: FTwindow
integer         :: spec_mode
integer         :: random_xi,random_xf,random_yi,random_yf,random_zi,random_zf,layrsSERS,izSERSshift
! >>> TF/SF
integer         :: ix_min,ix_max,iy_min,iy_max,iz_min,iz_max
integer         :: tfsf_xi,tfsf_xf,tfsf_yi,tfsf_yf,tfsf_zi,tfsf_zf
integer         :: save_simul=0
real(precision) :: RATIO_ANGLE
character(100)  :: source_path
character(100)  :: source_path_scatt
integer         :: source_time
integer         :: source_to
integer         :: source_tf
integer         :: source_tf_exp
integer         :: tfsf_slab_launcher
real(precision)            :: pump_DT,probe_DT,def_pulse,DTpump_probe
!     TF/SF
! >>> PML and UPML layers
integer            :: ABCmodel
real(precision)    :: wzPML
integer            :: nPML
integer            :: pmlxi,pmlxi_aux
integer            :: pmlxf,pmlxf_aux
integer            :: pmlyi,pmlyi_aux
integer            :: pmlyf,pmlyf_aux
integer            :: pmlzi,pmlzi_x
integer            :: pmlzf,pmlzf_x
integer            :: lcell,kcell
integer            :: pml_taflove !Taflove if it is used the conditions for conductivity (sg), and k from Taflove. pml_taflove=2 because we read pml_taflove instead of nPML in inputFDTD.dat
                                  !since 21/nov/18. nPML=2 and most of the inputs before used this value.
!    PML and UPML layers <<<

! >>> Structure parameters
integer :: Dh,hi,hf,aix,afx

! >>> Dielectric constants
real(precision)    :: eps(3)
complex(precision) :: eps_bbox
real(precision)    :: eps_in(1:3)    
real(precision) :: epsa(1:2,1:3,1:3) !>>> 1=region 1; 2=region=2
real(precision) :: epsa_inv(1:2,1:3,1:3) !>>> 1=region 1; 2=region=2
real(precision) :: mu_inv_bbox = 1.0

! >>> Parámetros modos de ejecución     
integer        :: runmode 
integer        :: metalmodel
integer,save   :: m_structure       
integer        :: no_downloads
integer        :: list_modes
integer        :: OUTPUT_TIME_STEP
integer        :: texit
real(precision) :: texit_fs

!>>> Fourier near-field
integer                        :: rwmode
integer                        :: emmode !0:Campo E y H;1:Campo E;2:Campo H
integer                        :: nof,nof_bs !nof for Near-field and nof_bs for band-structure (or dispersion relation)
real(precision),pointer        :: fw(:),fw_bs(:) !fw for Near-field and fw_bs for band-structure (or dispersion relation)
integer                        :: sxi,sxf,syi,syf,szi,szf,mxi,mxf,myi,myf,mzi,mzf
integer                        :: file_direct
real(precision)                :: dlamda
integer                        :: to_Fourier,nowav,sample_Fourier
real(precision)                :: fourier_fs
integer                        :: plane_wave
!    Fourier near-field <<<

!>>> MESH
integer :: LayerMESH
character(100)::layer_filename  
!    MESH <<<

!     Parámetros modos de ejecución <<<     

!>>> Messages List
 !>>> Message target
 character(100) msg
 character(100) msgaux
 !    Message target <<<

 !>>> Error Messages
 character(100),parameter::Err_0  = '...Error -> Unknown'
 character(100),parameter::Err_1  = '...Error -> Non-existent input file - Subroutine: setparam or setparam_grating'
 character(100),parameter::Err_2  = '...Error -> Invalid input file - Subroutine: setparam or setparam_grating'   
 character(100),parameter::Err_3  = '...Error -> Failure in energy conservation'
 character(100),parameter::Err_4  = '...Error -> Failure in div D'
 character(100),parameter::Err_5  = '...Error -> Failure in div B'
 character(100),parameter::Err_6  = '...Error -> Error in basis set generation - Subroutine: leftvector'
 character(100),parameter::Err_7  = '...Error -> No current carrying incident modes - Subroutine: leftvector'
 character(100),parameter::Err_8  = '...Error -> Unknown file path - Subroutine:'
 !    Error Messages <<<

 !>>> Structure Number and Names
 integer       ,parameter::nostr_2  = 2
 character(100),parameter::str_2    = 'Slit Array' 
 !    Structure Number and Names >>>

 !>>> Other names 
 character(100),parameter::runmode_1    = 'Time evolution options' 
 integer       ,parameter::norunmode_1  = 1
 character(100),parameter::runmode_2    = 'Fourier transformation options'  
 integer       ,parameter::norunmode_2  = 2
 
 character(100),parameter::SERS_X    = 'X=cte' 
 integer       ,parameter::noSER_X  = 1
 character(100),parameter::SERS_Y    = 'Y=cte' 
 integer       ,parameter::noSER_Y  = 2
 character(100),parameter::SERS_Z    = 'Z=cte' 
 integer       ,parameter::noSER_Z  = 3
 !    Other names >>>

 !>>> Messages
 character(100),parameter::msg0    = 'Mesh options'  
 character(100),parameter::msg1    = 'Job options'   
 character(100),parameter::msg2    = 'Structure check'
 character(100),parameter::msg3    = 'List modes'  
 character(100),parameter::msg5    = 'Dielectric constants' 
 character(100),parameter::msg7    = "Extra parameters"
 character(100),parameter::msg12    = 'Source options' 
 !    Messages >>>
!    Messages List <<<

!  >>> Para nombres de archivo en opción bucle de procesos
character(5)  :: idx
character(10) :: idx_name
!      Para nombres de archivo en opción bucle de procesos <<<

integer :: stat_alloc
real(precision)    :: delta_x
real(precision)    :: delta_y
real(precision)    :: delta_z

!>>>  NONLINEAR media
integer           :: NONLINEAR
real(precision)   :: Eoxy
character(100),parameter::msg8    = 'No Lineal'
integer           :: curr_zo,curr_zf
integer           :: cosine 
integer           :: dx,aix_slit,afx_slit,region_nl(3)
real(precision)   :: tdecay_PW,tincr_PW

!<<<  NONLINEAR media

real(precision)  :: Epump,wlpump,wlprobe
real(precision)  :: zopump,zoprobe,Gpump,Gprobe
real(precision)  :: xi_gain,ag_gain,w_gain
!<<< GAIN media

end module parameters
