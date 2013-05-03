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
! SUBROUTINE scatter
! SUBROUTINE add_shunt_admittance
! SUBROUTINE add_series_impedance
!
! NAME
!     scatter
!
! DESCRIPTION
!     SCN scatter 
!     
! COMMENTS
!     New idea is to incrementally build a Thevenin equivalent at the
!     cell centre so as to include materials and cables in a unified manner
!     Field output then also becomes straightforward if recorded at this time
!
!  Note 1: cell centre sources have a factor of 0.5 applied which gives consistency
!          between cell and face excitations but leads to discrepancies in cell output at the source
!
!  Note 2: PEC and PMC cells give the correct E and H fields at the cell centres
!          but maybe not the best way to implement it as allows some field through
!
! HISTORY
!
!     started 14/08/2012 CJS
!
!
SUBROUTINE scatter

USE TLM_general
USE TLM_excitation
USE TLM_output
USE TLM_volume_materials
USE mesh
USE constants

IMPLICIT NONE

! local variables

  integer cx,cy,cz
  
  real*8 	:: Vx,Vy,Vz,Ix,Iy,Iz,V_min,V_max

! SCN link line data 
  real*8 	:: Yx,Yy,Yz,Zx,Zy,Zz
  real*8	:: VEx,VEy,VEz,VHx,VHy,VHz
  
! SCN inductive and capcaitive stub data  
  real*8 	:: Ysx,Ysy,Ysz,Zsx,Zsy,Zsz
  real*8 	:: Vsx,Vsy,Vsz,Vox,Voy,Voz
  
! SCN electric and magnetic loss data  
  real*8 	:: Gex,Gey,Gez,Rmx,Rmy,Rmz
  
  real*8	:: Jsx,Jsy,Jsz,Msx,Msy,Msz
  
  real*8	:: Z04,Y04
  
  real*8	:: const1
  real*8	:: const2

  integer	:: cell_number
  integer	:: special_cell_count
  
! Material update parameters  
  integer	:: material_number
  integer	:: material_type
  integer	:: filter_number
  
! Output parameters  
  integer	:: output_number
  integer	:: output_cell_number
  real*8	:: field(6)
  
! Excitation parameters  
  integer	:: excitation_point
  integer	:: excitation_number
  integer	:: excitation_function_number
  integer	:: field_component
  integer	:: source_type
  real*8	:: source_value
  
  integer	:: hard_source_factor
  
  logical :: op_flag


! START
  
!  CALL write_line('CALLED: scatter',0,output_to_screen_flag)

  op_flag=.FALSE.
  
  Z04=Z0*4d0
  Y04=Y0*4d0

  const1=dl*dl
  const2=dl*dl

  cell_number=0
  output_cell_number=0
  
  do cz=1,nz
    do cy=1,ny
      do cx=1,nx
      
        cell_number=cell_number+1
	
	if (cell_centre_update_code(cell_number).eq.0) then   ! free space update, no output
	
          Vx=(	  V(Vx_ymin,cx,cy,cz)+V(Vx_ymax,cx,cy,cz)	&
	  	 +V(Vx_zmin,cx,cy,cz)+V(Vx_zmax,cx,cy,cz) )/ 2d0       
	      
          Vy=(	  V(Vy_xmin,cx,cy,cz)+V(Vy_xmax,cx,cy,cz)	&
		 +V(Vy_zmin,cx,cy,cz)+V(Vy_zmax,cx,cy,cz) )/ 2d0       
	      
          Vz=(	  V(Vz_xmin,cx,cy,cz)+V(Vz_xmax,cx,cy,cz)	&
		 +V(Vz_ymin,cx,cy,cz)+V(Vz_ymax,cx,cy,cz) )/ 2d0       
	          
	  Ix=( V(Vz_ymax,cx,cy,cz)-V(Vy_zmax,cx,cy,cz)	&
	      -V(Vz_ymin,cx,cy,cz)+V(Vy_zmin,cx,cy,cz) )/ (2d0*Z0)    
	
	  Iy=( V(Vx_zmax,cx,cy,cz)-V(Vz_xmax,cx,cy,cz)	&
	      -V(Vx_zmin,cx,cy,cz)+V(Vz_xmin,cx,cy,cz) )/ (2d0*Z0)    
	
	  Iz=( V(Vy_xmax,cx,cy,cz)-V(Vx_ymax,cx,cy,cz)	&
	      -V(Vy_xmin,cx,cy,cz)+V(Vx_ymin,cx,cy,cz) )/ (2d0*Z0)   

	else
! cell_update_code points to arrays which tell us about material number and output_number

	  special_cell_count=cell_centre_update_code(cell_number)
	  
	  material_number=cell_update_code_to_material_data(special_cell_count,1)
	  if (material_number.ne.0) then
	    material_type=volume_material_list(material_number)%type
	  else
	    material_type=0
	  end if
	  output_number=cell_update_code_to_output_number(special_cell_count)
	  excitation_number=cell_update_code_to_excitation_number(special_cell_count)

! Calculate internal V and I with no stubs
! Calculate equivalent circuit parameters in x, y and z directions for Voltage calculation
          VEx=(	V(Vx_ymin,cx,cy,cz)+V(Vx_ymax,cx,cy,cz)+V(Vx_zmin,cx,cy,cz)+V(Vx_zmax,cx,cy,cz) )/2d0   	      
          VEy=(	V(Vy_xmin,cx,cy,cz)+V(Vy_xmax,cx,cy,cz)+V(Vy_zmin,cx,cy,cz)+V(Vy_zmax,cx,cy,cz) )/2d0     	      
          VEz=(	V(Vz_xmin,cx,cy,cz)+V(Vz_xmax,cx,cy,cz)+V(Vz_ymin,cx,cy,cz)+V(Vz_ymax,cx,cy,cz) )/2d0    
	  Yx=Y04
	  Yy=Y04
	  Yz=Y04
	          
! Calculate equivalent circuit parameters in x, y and z directions for Current calculation
	  VHx=2d0*( V(Vz_ymax,cx,cy,cz)-V(Vy_zmax,cx,cy,cz)-V(Vz_ymin,cx,cy,cz)+V(Vy_zmin,cx,cy,cz) )	
	  VHy=2d0*( V(Vx_zmax,cx,cy,cz)-V(Vz_xmax,cx,cy,cz)-V(Vx_zmin,cx,cy,cz)+V(Vz_xmin,cx,cy,cz) )  	
	  VHz=2d0*( V(Vy_xmax,cx,cy,cz)-V(Vx_ymax,cx,cy,cz)-V(Vy_xmin,cx,cy,cz)+V(Vx_ymin,cx,cy,cz) ) 
          Zx=Z04
	  Zy=Z04
	  Zz=Z04

          if (material_type.eq.volume_material_type_PEC) then  
! Set internal voltages to zero	  

	    VEx=0d0
            VEy=0d0
            VEz=0d0
	    
          else if (material_type.eq.volume_material_type_PMC) then
! Set internal currents to zero	  
	  
            VHx=0d0
            VHy=0d0
            VHz=0d0
	  
          else if (material_type.eq.volume_material_type_DISPERSIVE) then
	  
	    filter_number=cell_update_code_to_material_data(special_cell_count,2)

! Admittance and voltage source of capacitaive (open circuit) stubs	    
            Ysx=volume_material_list(material_number)%Ys*Y0
	    Ysy=volume_material_list(material_number)%Ys*Y0
	    Ysz=volume_material_list(material_number)%Ys*Y0
	    Vox=2d0*Vo(filter_number)
	    Voy=2d0*Vo(filter_number+1)
	    Voz=2d0*Vo(filter_number+2)
	    
	    CALL add_shunt_admittance(VEx,Yx,Vox,Ysx,VEx,Yx)
	    CALL add_shunt_admittance(VEy,Yy,Voy,Ysy,VEy,Yy)
	    CALL add_shunt_admittance(VEz,Yz,Voz,Ysz,VEz,Yz)

! Electric loss 	    
            Gex=volume_material_list(material_number)%Ge*Y0
	    Gey=volume_material_list(material_number)%Ge*Y0
	    Gez=volume_material_list(material_number)%Ge*Y0
	    
	    CALL add_shunt_admittance(VEx,Yx,0d0,Gex,VEx,Yx)
	    CALL add_shunt_admittance(VEy,Yy,0d0,Gey,VEy,Yy)
	    CALL add_shunt_admittance(VEz,Yz,0d0,Gez,VEz,Yz)
	    
! Impedance and voltage source of inductive (short circuit) stubs	    
	    Zsx=volume_material_list(material_number)%Zs*Z0
	    Zsy=volume_material_list(material_number)%Zs*Z0
	    Zsz=volume_material_list(material_number)%Zs*Z0
	    Vsx=2d0*Vs(filter_number)
	    Vsy=2d0*Vs(filter_number+1)
	    Vsz=2d0*Vs(filter_number+2)
	    
	    CALL add_series_impedance(VHx,Zx,Vsx,Zsx,VHx,Zx)
	    CALL add_series_impedance(VHy,Zy,Vsy,Zsy,VHy,Zy)
	    CALL add_series_impedance(VHz,Zz,Vsz,Zsz,VHz,Zz)

! Magnetic loss 	    
            Rmx=volume_material_list(material_number)%Rm*Z0
	    Rmy=volume_material_list(material_number)%Rm*Z0
	    Rmz=volume_material_list(material_number)%Rm*Z0
	    
	    CALL add_series_impedance(VHx,Zx,0d0,Rmx,VHx,Zx)
	    CALL add_series_impedance(VHy,Zy,0d0,Rmy,VHy,Zy)
	    CALL add_series_impedance(VHz,Zz,0d0,Rmz,VHz,Zz)  

! Calculate electric and magnetic current source terms from filters 	    
            call timeshift_Zfilter(volume_material_eps_filter_data(filter_number))
            call evaluate_Zfilter(volume_material_list(material_number)%eps_Z,	      &
	        		  volume_material_eps_filter_data(filter_number),       &
	        		  0d0)  	    
	    Jsx=volume_material_eps_filter_data(filter_number)%f    
	    
            call timeshift_Zfilter(volume_material_eps_filter_data(filter_number+1))
            call evaluate_Zfilter(volume_material_list(material_number)%eps_Z,	      &
	        		  volume_material_eps_filter_data(filter_number+1),       &
	        		  0d0)  	    
	    Jsy=volume_material_eps_filter_data(filter_number+1)%f    
	    
            call timeshift_Zfilter(volume_material_eps_filter_data(filter_number+2))
            call evaluate_Zfilter(volume_material_list(material_number)%eps_Z,	      &
	        		  volume_material_eps_filter_data(filter_number+2),       &
	        		  0d0)  	    
	    Jsz=volume_material_eps_filter_data(filter_number+2)%f    
	    
            call timeshift_Zfilter(volume_material_mu_filter_data(filter_number))
            call evaluate_Zfilter(volume_material_list(material_number)%mu_Z,	      &
	        		  volume_material_mu_filter_data(filter_number),       &
	        		  0d0)  	    
	    Msx=volume_material_mu_filter_data(filter_number)%f    
	    
            call timeshift_Zfilter(volume_material_mu_filter_data(filter_number+1))
            call evaluate_Zfilter(volume_material_list(material_number)%mu_Z,	      &
	        		  volume_material_mu_filter_data(filter_number+1),       &
	        		  0d0)  	    
	    Msy=volume_material_mu_filter_data(filter_number+1)%f    
	    
            call timeshift_Zfilter(volume_material_mu_filter_data(filter_number+2))
            call evaluate_Zfilter(volume_material_list(material_number)%mu_Z,	      &
	        		  volume_material_mu_filter_data(filter_number+2),       &
	        		  0d0)  	    
	    Msz=volume_material_mu_filter_data(filter_number+2)%f   
	     	    
	    VEx=VEx+Jsx*const1/Yx
	    VEy=VEy+Jsy*const1/Yy
	    VEz=VEz+Jsz*const1/Yz

	    VHx=VHx-Msx*const2
	    VHy=VHy-Msy*const2
	    VHz=VHz-Msz*const2
	    
	  end if  ! DISPERSIVE material

! calculate internal Voltage and current 
          Vx=VEx
          Vy=VEy
          Vz=VEz
	  
          Ix=VHx/Zx
          Iy=VHy/Zy
          Iz=VHz/Zz

! Excitation		
          if (excitation_number.ne.0) then
	  
            excitation_point=cell_update_code_to_excitation_number(special_cell_count)
	    
!	    write(*,*)'Excitation_point=',excitation_point
	    
	    excitation_function_number=excitation_points(excitation_point)%excitation_function_number
	    field_component=excitation_points(excitation_point)%field_component
	    source_type=excitation_points(excitation_point)%source_type
	    
!	    write(*,*)'excitation_function_number=',excitation_function_number
!	    write(*,*)'field_component=',field_component
!	    write(*,*)'source_type=',source_type
!	    write(*,*)'timestep=',timestep
	    
	    source_value=excitation_functions(excitation_function_number)%value(timestep)

!	    write(*,*)'source_value=',source_value
	    
	    if (source_type.EQ.source_type_soft) then
	      hard_source_factor=1
	    else
	      hard_source_factor=0
	    end if
	    
            if (field_component.eq.Ex) then
	      Vx=hard_source_factor*Vx-source_value*dl/2d0
	    else if (field_component.eq.Ey) then
	      Vy=hard_source_factor*Vy-source_value*dl/2d0
	    else if (field_component.eq.Ez) then
	      Vz=hard_source_factor*Vz-source_value*dl/2d0
	    else if (field_component.eq.Hx) then
	      Ix=hard_source_factor*Ix+source_value*dl/2d0
	    else if (field_component.eq.Hy) then
	      Iy=hard_source_factor*Iy+source_value*dl/2d0
	    else if (field_component.eq.Hz) then
	      Iz=hard_source_factor*Iz+source_value*dl/2d0
            end if
	    
!	    write(*,*)'finished excitation'
	    
	  end if ! excitation_number.ne.0
! Output		
          if (output_number.ne.0) then
	  
!	    write(*,*)'Output_point=',output_number
	    
	    output_cell_number=output_cell_number+1
	  
! calculate cell centre fields		     
	    field(Ex)=-Vx/dl
	    field(Ey)=-Vy/dl
	    field(Ez)=-Vz/dl
	
	    field(Hx)= Ix/dl
	    field(Hy)= Iy/dl
	    field(Hz)= Iz/dl
	  
            cell_output_field(output_cell_number,1:6)=field(1:6)
	  
	  end if ! output_number.ne.0

! If this is a dispersive material then update the filter functions
          if (material_type.eq.volume_material_type_DISPERSIVE) then
	  
	    filter_number=cell_update_code_to_material_data(special_cell_count,2)
	    
	    field(Ex)=-Vx/dl
	    field(Ey)=-Vy/dl
	    field(Ez)=-Vz/dl
	
	    field(Hx)= Ix/dl
	    field(Hy)= Iy/dl
	    field(Hz)= Iz/dl
	      
! permittivity filters	    	    
            call evaluate_Zfilter(volume_material_list(material_number)%eps_Z,  	      &
	    			  volume_material_eps_filter_data(filter_number),     &
	        		  field(Ex))					     
            call evaluate_Zfilter(volume_material_list(material_number)%eps_Z,  	      &
	    			  volume_material_eps_filter_data(filter_number+1),   &
	        		  field(Ey))	    
            call evaluate_Zfilter(volume_material_list(material_number)%eps_Z,  	      &
	    			  volume_material_eps_filter_data(filter_number+2),   &
	        		  field(Ez))   
				        
! permeability filters	    	    
            call evaluate_Zfilter(volume_material_list(material_number)%mu_Z,	      &
	    			  volume_material_mu_filter_data(filter_number),      &
	        		  field(Hx))					     
            call evaluate_Zfilter(volume_material_list(material_number)%mu_Z,	      &
	    			  volume_material_mu_filter_data(filter_number+1),    &
	        		  field(Hy))	    
            call evaluate_Zfilter(volume_material_list(material_number)%mu_Z,	      &
	    			  volume_material_mu_filter_data(filter_number+2),    &
	        		  field(Hz))	 

! open circuit capacitive stubs: calculate incident voltage at next timestep		      
	    Vo(filter_number)  =Vx-Vo(filter_number)
            Vo(filter_number+1)=Vy-Vo(filter_number+1)
            Vo(filter_number+2)=Vz-Vo(filter_number+2)

! short circuit inductive stubs, note - sign includes the reflection coefficient -1 at short circuit 		     
	    Vs(filter_number  )=-(Vs(filter_number  )-Ix*volume_material_list(material_number)%Zs*Z0)
	    Vs(filter_number+1)=-(Vs(filter_number+1)-Iy*volume_material_list(material_number)%Zs*Z0)
	    Vs(filter_number+2)=-(Vs(filter_number+2)-Iz*volume_material_list(material_number)%Zs*Z0)
	      
          end if ! DISPERSIVE material

        end if ! cell_update_code
	
! Calculate Scattered voltage pulses into link lines	

        V_min=V(Vx_ymin,cx,cy,cz)
        V_max=V(Vx_ymax,cx,cy,cz)     
	V(Vx_ymin,cx,cy,cz)=Vx-Z0*Iz-V_max
	V(Vx_ymax,cx,cy,cz)=Vx+Z0*Iz-V_min

        V_min=V(Vx_zmin,cx,cy,cz)
        V_max=V(Vx_zmax,cx,cy,cz)     
	V(Vx_zmin,cx,cy,cz)=Vx+Z0*Iy-V_max
	V(Vx_zmax,cx,cy,cz)=Vx-Z0*Iy-V_min
	
        V_min=V(Vy_zmin,cx,cy,cz)
        V_max=V(Vy_zmax,cx,cy,cz)     
	V(Vy_zmin,cx,cy,cz)=Vy-Z0*Ix-V_max
	V(Vy_zmax,cx,cy,cz)=Vy+Z0*Ix-V_min
	
        V_min=V(Vy_xmin,cx,cy,cz)
        V_max=V(Vy_xmax,cx,cy,cz)     
	V(Vy_xmin,cx,cy,cz)=Vy+Z0*Iz-V_max
	V(Vy_xmax,cx,cy,cz)=Vy-Z0*Iz-V_min
	
        V_min=V(Vz_xmin,cx,cy,cz)
        V_max=V(Vz_xmax,cx,cy,cz)     
	V(Vz_xmin,cx,cy,cz)=Vz-Z0*Iy-V_max
	V(Vz_xmax,cx,cy,cz)=Vz+Z0*Iy-V_min
	
        V_min=V(Vz_ymin,cx,cy,cz)
        V_max=V(Vz_ymax,cx,cy,cz)     
	V(Vz_ymin,cx,cy,cz)=Vz+Z0*Ix-V_max
	V(Vz_ymax,cx,cy,cz)=Vz-Z0*Ix-V_min
	
      end do  ! next z cell
    end do    ! next y cell
  end do      ! next x cell
  
!  CALL write_line('FINISHED: scatter',0,output_to_screen_flag)

  RETURN

END SUBROUTINE scatter
!
! NAME
!     add_shunt_admittance
!
! DESCRIPTION
!     add_shunt_admittance
!     
! COMMENTS
!     
! HISTORY
!
!     started 5/09/2012 CJS
!
!
SUBROUTINE add_shunt_admittance(V1,Y1,V2,Y2,V3,Y3)

IMPLICIT NONE

  real*8 	:: V1,Y1,V2,Y2,V3,Y3

! local variables
  
  real*8 	:: V,Y

! START

  V=(V1*Y1+V2*Y2)/(Y1+Y2)
  Y=Y1+Y2

  V3=V
  Y3=Y

  RETURN
  
END SUBROUTINE add_shunt_admittance
!
! NAME
!     add_series_impedance
!
! DESCRIPTION
!     add_series_impedance
!     
! COMMENTS
!     
! HISTORY
!
!     started 5/09/2012 CJS
!
!
SUBROUTINE add_series_impedance(V1,Z1,V2,Z2,V3,Z3)

IMPLICIT NONE

  real*8 	:: V1,Z1,V2,Z2,V3,Z3

! local variables
  
  real*8 	:: V,Z

! START

  V=V1+V2
  Z=Z1+Z2

  V3=V
  Z3=Z

  RETURN
  
END SUBROUTINE add_series_impedance
