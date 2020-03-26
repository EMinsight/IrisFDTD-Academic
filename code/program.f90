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

    program IrisFDTD
    use MAIN_PROGRAM
    use parameters
    implicit none
    integer :: i,io_main,if_main,i_aux
    character(100)           :: ITER
   
    open(unit=inFDTD,file='inputFDTD.dat',status='old',action='read')

    read(inFDTD,'(A100)') ITER

    if(ITER/="<<LOOP>>")then
        rewind(inFDTD)
        call main()
    else
        read(inFDTD,*) io_main
        read(inFDTD,*) if_main
        read(inFDTD,*) delta_x
        read(inFDTD,*) delta_y
        read(inFDTD,*) delta_z

        i_aux=io_main
        if(io_main/=if_main)then
            do i = io_main,if_main,1
                i_aux=i  
                call main(i_aux)
                if(i_aux==if_main) exit
            end do
        else
            i=io_main
            call main(i)
        end if
    end if

    close(unit=inFDTD)
    stop
    end program IrisFDTD
