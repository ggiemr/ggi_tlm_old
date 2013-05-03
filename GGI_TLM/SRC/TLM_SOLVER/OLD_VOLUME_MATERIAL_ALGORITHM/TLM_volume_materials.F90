!
!    GGI_TLM Time domain electromagnetic field solver based on the TLM method
!    Copyright (C) 2013  Chris Smartt
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.   
!
! SUBROUTINE set_volume_material_mesh
! SUBROUTINE calculate_volume_material_filter_coefficients
! SUBROUTINE allocate_volume_material_filter_data
! SUBROUTINE volume_material_update
!
! NAME
!     set_volume_material_mesh
!
! DESCRIPTION
!     
!     Count the number of cells of each material and construct the 
!     appropriate lists for the main solver
!     
! COMMENTS
!     Allocate a local mesh array and fill with material data. This ensures that
!     a cell doesn't get given more than one material property.
!     In overlapping volumes the last material assignment takes precidence
!
! HISTORY
!
!     started 5/09/2012 CJS
!
!
SUBROUTINE set_volume_material_mesh

USE TLM_general
USE geometry
USE TLM_volume_materials
USE mesh
USE cell_parameters

IMPLICIT NONE

! local variables

  integer	:: material_number
  integer	:: volume_number
  integer	:: cell,number_of_cells
  integer	:: total_number_of_material_cells
  integer	:: cx,cy,cz
  
  integer	:: i
  integer	:: cell_number

! START
  
  CALL write_line('CALLED: set_volume_material_mesh',0,output_to_screen_flag)

! INITIALISE VOLUME MATERIALS  
    
  CALL write_line_integer('n_volume_materials=',n_volume_materials,0,output_to_screen_flag)

  do material_number=1,n_volume_materials
  
    total_number_of_material_cells=0
    
    CALL write_line_integer('volume material number=',material_number,0,output_to_screen_flag)
    CALL write_line_integer('Number of geometric volumes=',	&
                             volume_material_list(material_number)%n_volumes,0,output_to_screen_flag)

! loop over the geometric volumes with this material type and initialise the TLM cells
    do i=1,volume_material_list(material_number)%n_volumes
    
      volume_number=volume_material_list(material_number)%volume_list(i)

      CALL write_line_integer('Geometric volume material number=',volume_number,0,output_to_screen_flag)
      
      number_of_cells=problem_volumes(volume_number)%number_of_cells

      CALL write_line_integer('Number of cells=',number_of_cells,0,output_to_screen_flag)
      
      total_number_of_material_cells=total_number_of_material_cells+number_of_cells
      
      do cell=1,number_of_cells
    
	cx=problem_volumes(volume_number)%cell_list(cell)%cell%i
	cy=problem_volumes(volume_number)%cell_list(cell)%cell%j
	cz=problem_volumes(volume_number)%cell_list(cell)%cell%k
	
        local_cell_material(cx,cy,cz)=material_number
	
      end do ! next cell
    
    end do ! next volume number

    CALL write_line_integer('Total number of material cells=',	&
                            total_number_of_material_cells,0,output_to_screen_flag)

  end do ! next volume material number
  
  CALL write_line('FINISHED: set_volume_material_mesh',0,output_to_screen_flag)

  RETURN

END SUBROUTINE set_volume_material_mesh
!
! Name calculate_volume_material_filter_coefficients
!     
!
! Description:
!      Once the timestep is found we can determine the volume material filter coefficients
!      
!
! Comments:
!      
!      
!
! History
!
!     started 04/09/12 CJS
!

SUBROUTINE calculate_volume_material_filter_coefficients()

USE TLM_general
USE TLM_volume_materials
USE filter_types
USE filter_operators
USE filter_functions
USE constants

IMPLICIT NONE

! variables passed to subroutine

! local_variables

  integer material_number
  
  type(Sfilter) :: eps_Sfilter1
  type(Zfilter) :: eps_Zfilter1,eps_Zfilter1b
  real*8	:: eps_f
  type(Zfilter) :: eps_Zfilter2
  
  type(Sfilter) :: mu_Sfilter1
  type(Zfilter) :: mu_Zfilter1,mu_Zfilter1b
  real*8	:: mu_f
  type(Zfilter) :: mu_Zfilter2
  
  type(Zfilter) :: bilinear_eps,bilinear_mu
  
  real*8	:: dt2

! function_types

! START

  write(*,*)'CALLED: calculate_volume_material_filter_coefficients'
  
  bilinear_eps%a%order=1
  bilinear_eps%b%order=1
  bilinear_eps=allocate_Zfilter(bilinear_eps%a%order,bilinear_eps%b%order)
  dt2=eps0*2d0/dt    
  bilinear_eps%a%coeff(0)= 1d0*dt2
  bilinear_eps%a%coeff(1)=-1d0*dt2
  bilinear_eps%b%coeff(0)= 1d0
  bilinear_eps%b%coeff(1)= 1d0

  bilinear_mu%a%order=1
  bilinear_mu%b%order=1
  bilinear_mu=allocate_Zfilter(bilinear_mu%a%order,bilinear_mu%b%order)
  dt2=mu0*2d0/dt    
  bilinear_mu%a%coeff(0)= 1d0*dt2
  bilinear_mu%a%coeff(1)=-1d0*dt2
  bilinear_mu%b%coeff(0)= 1d0
  bilinear_mu%b%coeff(1)= 1d0

  do material_number=1,n_volume_materials
   
! Z transform of Z parameter filter	
    
    if (volume_material_list(material_number)%type.EQ.volume_material_type_DISPERSIVE) then
    
! Z transform of permittivity filter	

!      write(*,*)'eps_S filter'
!      call write_Sfilter(volume_material_list(material_number)%eps_S,0)

      eps_Zfilter1=s_to_z(volume_material_list(material_number)%eps_S,dt) 
          
!      write(*,*)'eps_Z filter1'
!      call write_Zfilter(eps_Zfilter1,0)
    
! calculate fast permittivity and slow permittivity filter

      call Z_fast_slow_docomposition( eps_Zfilter1 ,eps_f  , eps_Zfilter1b )
          
!      write(*,*)'eps_f',eps_f
!      write(*,*)'eps_Z filter1b'
!      call write_Zfilter(eps_Zfilter1b,0)

      bilinear_eps%wnorm=eps_Zfilter1%wnorm
      bilinear_eps%T=eps_Zfilter1%T
   
      eps_Zfilter2=bilinear_eps*eps_Zfilter1b
!      write(*,*)'eps_Z filter2'
!      call write_Zfilter(eps_Zfilter2,0)
    
      volume_material_list(material_number)%eps_Z=eps_Zfilter2
      volume_material_list(material_number)%eps_f=eps_f
      
! Open circuit capacitive stub
      volume_material_list(material_number)%Ys=2d0*((volume_material_list(material_number)%eps_f/c0)*dl-2d0*dt)/dt
      	
      if (volume_material_list(material_number)%Ys.lt.0d0) then
! if we think this is a result of tolerance issuse then set to zero, else cause an error
        if (volume_material_list(material_number)%Ys.lt.-small) then
          GOTO 9000
	else
	  volume_material_list(material_number)%Ys=0d0
        end if
      end if
      
! Electric loss admittiance
      volume_material_list(material_number)%Ge=volume_material_list(material_number)%sigma_e*Z0*dl
      if (volume_material_list(material_number)%Ge.lt.0d0) then
! if we think this is a result of tolerance issuse then set to zero, else cause an error
        if (volume_material_list(material_number)%Ge.lt.-small) then
          GOTO 9010
	else
	  volume_material_list(material_number)%Ge=0d0
        end if
      end if
      
!      write(*,*)'Ys=',volume_material_list(material_number)%Ys
!      write(*,*)'Ge=',volume_material_list(material_number)%Ge
      
! Z transform of permeability filter	
    
!      call reciprocal_Sfilter(volume_material_list(material_number)%mu_S,mu_Sfilter1)      
!      mu_Zfilter1=s_to_z(mu_Sfilter1,dt) 

!      write(*,*)'mu_S filter'
!      call write_Sfilter(volume_material_list(material_number)%mu_S,0)

      mu_Zfilter1=s_to_z(volume_material_list(material_number)%mu_S,dt) 
       
!      write(*,*)'mu_Z filter1'
!      call write_Zfilter(mu_Zfilter1,0)
    
! calculate fast permeability and slow permeability filter...
    
      call Z_fast_slow_docomposition( mu_Zfilter1 ,mu_f  , mu_Zfilter1b )
          
      write(*,*)'mu_f',mu_f
!      write(*,*)'mu_Z filter1b'
!      call write_Zfilter(mu_Zfilter1b,0)

      bilinear_mu%wnorm=mu_Zfilter1%wnorm
      bilinear_mu%T=mu_Zfilter1%T
    
      mu_Zfilter2=bilinear_mu*mu_Zfilter1b
      
      write(*,*)'mu_Z filter2'
      call write_Zfilter(mu_Zfilter2,0)
    
      volume_material_list(material_number)%mu_Z=mu_Zfilter2
      volume_material_list(material_number)%mu_f=mu_f
      
! Short circuit inductive stub
      volume_material_list(material_number)%Zs=2d0*((volume_material_list(material_number)%mu_f/c0)*dl-2d0*dt)/dt	
      if (volume_material_list(material_number)%Zs.lt.0d0) then
! if we think this is a result of tolerance issuse then set to zero, else cause an error
        if (volume_material_list(material_number)%Zs.lt.-small) then
          GOTO 9020
	else
	  volume_material_list(material_number)%Zs=0d0
        end if
      end if
      
! Magnetic loss impedance
      volume_material_list(material_number)%Rm=volume_material_list(material_number)%sigma_m*dl/Z0
      if (volume_material_list(material_number)%Rm.lt.0d0) then
! if we think this is a result of tolerance issuse then set to zero, else cause an error
        if (volume_material_list(material_number)%Rm.lt.-small) then
          GOTO 9030
	else
	  volume_material_list(material_number)%Rm=0d0
        end if
      end if
      
      write(*,*)'Zs=',volume_material_list(material_number)%Zs
      write(*,*)'Rm=',volume_material_list(material_number)%Rm
		     
    end if ! material type is dispersive
    
  end do ! next material to set

  write(*,*)'FINISHED: calculate_volume_material_filter_coefficients'
  
  RETURN
  
9000 CALL write_line('Error in :calculate_volume_material_filter_coefficients',0,.TRUE.)
     CALL write_line_integer('Negative Capacitive stub, material number',material_number,0,.TRUE.)
     write(*,*)volume_material_list(material_number)%Ys
     STOP
  
9010 CALL write_line('Error in :calculate_volume_material_filter_coefficients',0,.TRUE.)
     CALL write_line_integer('Negative Electric loss, material number',material_number,0,.TRUE.)
     write(*,*)volume_material_list(material_number)%Ge
     STOP
  
9020 CALL write_line('Error in :calculate_volume_material_filter_coefficients',0,.TRUE.)
     CALL write_line_integer('Negative Inductive stub, material number',material_number,0,.TRUE.)
     write(*,*)volume_material_list(material_number)%Zs
     STOP
  
9030 CALL write_line('Error in :calculate_volume_material_filter_coefficients',0,.TRUE.)
     CALL write_line_integer('Negative Magnetic loss, material number',material_number,0,.TRUE.)
     write(*,*)volume_material_list(material_number)%Rm
     STOP

  
END SUBROUTINE calculate_volume_material_filter_coefficients
!
! NAME
!     allocate_volume_material_filter_data
!
! DESCRIPTION
!     allocate the memory required for volume material filter data
!     
! COMMENTS
!
!
! HISTORY
!
!     started 4/09/2012 CJS
!
!
SUBROUTINE allocate_volume_material_filter_data

USE TLM_general
USE TLM_volume_materials
USE filter_types
USE filter_operators
USE filter_functions
USE mesh

IMPLICIT NONE

! local variables

  integer	:: cx,cy,cz
  
  integer	:: material_number
  integer	:: cell
  integer	:: volume_filter_number
  integer       :: epsaorder,epsborder
  integer       :: muaorder,muborder

! START
  
  CALL write_line('CALLED: allocate_volume_material_filter_data',0,output_to_screen_flag)
  
  n_volume_material_cells=0
  volume_material_storage=0

! loop over special cells
  DO cell=1,n_special_cells
  
    if (cell_update_code_to_material_data(cell,1).ne.0) then
! this is a material cell add to the memory allocation list if required

      material_number=cell_update_code_to_material_data(cell,1)
      
      if ( (volume_material_list(material_number)%type.eq.volume_material_type_DISPERSIVE) ) then
	  	      
        n_volume_material_cells=n_volume_material_cells+1

        cell_update_code_to_material_data(cell,2)=volume_material_storage+1
        volume_material_storage=volume_material_storage+3   ! note 3 filters per cell, one for each field component    

      end if ! volume material at this cell cell

    end if ! is this a material cell? 

  end do ! next special cell

! allocate volume material filter data storage for each filter required
  if (volume_material_storage.ne.0) then
  
    allocate (volume_material_eps_filter_data(1:volume_material_storage)) 
    allocate (volume_material_mu_filter_data(1:volume_material_storage)) 
    
    allocate (Vs(1:volume_material_storage)) 
    allocate (Vo(1:volume_material_storage)) 
    
    Vs(1:volume_material_storage)=0d0
    Vo(1:volume_material_storage)=0d0
    
  end if
 
! loop over filters allocating memory for the particular material order
 

! loop over special cells
  DO cell=1,n_special_cells
  
    if (cell_update_code_to_material_data(cell,1).ne.0) then
! this is a material cell add to the memory allocation list if required

      material_number=cell_update_code_to_material_data(cell,1)
      
      if ( (volume_material_list(material_number)%type.eq.volume_material_type_DISPERSIVE) ) then
	  	      
	epsaorder=volume_material_list(material_number)%eps_Z%a%order
	epsborder=volume_material_list(material_number)%eps_Z%b%order
	muaorder=volume_material_list(material_number)%mu_Z%a%order
	muborder=volume_material_list(material_number)%mu_Z%b%order

! x polarisation 
        volume_filter_number=cell_update_code_to_material_data(cell,2)
	  
        volume_material_eps_filter_data(volume_filter_number)=allocate_Zfilter_response(epsaorder,epsborder)	
        volume_material_mu_filter_data(volume_filter_number)=allocate_Zfilter_response(muaorder,muborder)	

! y polarisation 
        volume_filter_number=cell_update_code_to_material_data(cell,2)+1
	  
        volume_material_eps_filter_data(volume_filter_number)=allocate_Zfilter_response(epsaorder,epsborder)	
        volume_material_mu_filter_data(volume_filter_number)=allocate_Zfilter_response(muaorder,muborder)	

! z polarisation 
        volume_filter_number=cell_update_code_to_material_data(cell,2)+2
	  
        volume_material_eps_filter_data(volume_filter_number)=allocate_Zfilter_response(epsaorder,epsborder)	
        volume_material_mu_filter_data(volume_filter_number)=allocate_Zfilter_response(muaorder,muborder)	

      end if ! volume material at this cell

    end if ! is this a material cell? 

  end do ! next special cell
 
 
  CALL write_line('FINISHED: allocate_volume_material_filter_data',0,output_to_screen_flag)


  RETURN

END SUBROUTINE allocate_volume_material_filter_data
