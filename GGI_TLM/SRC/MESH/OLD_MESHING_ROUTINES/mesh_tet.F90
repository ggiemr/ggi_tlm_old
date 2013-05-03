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
!SUBROUTINE mesh_tet
!
! NAME
!     SUBROUTINE mesh_tet
!
! DESCRIPTION
!     mesh_tet:
!
!     mesh a single tet by testing whether the centre point of a cell is contained within the tet
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 30/08/2012 CJS
!
!
  SUBROUTINE mesh_tet(volume_number,tet_number,cell_count,set_mesh_flag)
  
USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE file_information
USE constants

IMPLICIT NONE

integer				:: volume_number
integer				:: tet_number
integer				:: cell_count
logical				:: set_mesh_flag

! local variables

  integer 		:: ix,iy,iz
  
  real*8		:: x1,y1,z1
  real*8		:: x2,y2,z2
  real*8		:: x3,y3,z3
  real*8		:: x4,y4,z4
 
  type(xyz)		:: min_point,max_point
  type(ijk)		:: min_cell,max_cell
  
  type(xyz)		:: point
  type(ijk)		:: cell
  
  logical		:: inside_tet
    
! START

!  CALL write_line('CALLED: mesh_tet',0,output_to_screen_flag)
  
  x1=problem_volumes(volume_number)%tet_list(tet_number)%vertex(1)%x
  y1=problem_volumes(volume_number)%tet_list(tet_number)%vertex(1)%y
  z1=problem_volumes(volume_number)%tet_list(tet_number)%vertex(1)%z
  
  x2=problem_volumes(volume_number)%tet_list(tet_number)%vertex(2)%x
  y2=problem_volumes(volume_number)%tet_list(tet_number)%vertex(2)%y
  z2=problem_volumes(volume_number)%tet_list(tet_number)%vertex(2)%z
  
  x3=problem_volumes(volume_number)%tet_list(tet_number)%vertex(3)%x
  y3=problem_volumes(volume_number)%tet_list(tet_number)%vertex(3)%y
  z3=problem_volumes(volume_number)%tet_list(tet_number)%vertex(3)%z
  
  x4=problem_volumes(volume_number)%tet_list(tet_number)%vertex(4)%x
  y4=problem_volumes(volume_number)%tet_list(tet_number)%vertex(4)%y
  z4=problem_volumes(volume_number)%tet_list(tet_number)%vertex(4)%z

! get the extent of the triangle in x,y and z
  min_point%x=min(x1,x2,x3,x4)
  min_point%y=min(y1,y2,y3,y4)
  min_point%z=min(z1,z2,z3,z4)
  
  max_point%x=max(x1,x2,x3,x4)
  max_point%y=max(y1,y2,y3,y4)
  max_point%z=max(z1,z2,z3,z4)
    
! get the extent of the triangle in cells in x,y and z, note cells can be outside the mesh here
  CALL point_to_cell_no_checks(min_point,min_cell)
  CALL point_to_cell_no_checks(max_point,max_cell)
  
! Loop over cells and check whether the centre point is inside the tet
  do ix=min_cell%i,max_cell%i
    do iy=min_cell%j,max_cell%j 
      do iz=min_cell%k,max_cell%k
      
        cell%i=ix
        cell%j=iy
        cell%k=iz
      
        CALL get_cell_centre_coordinate(cell,point)

        CALL is_point_inside_tet(point,problem_volumes(volume_number)%tet_list(tet_number),inside_tet)
	       
        if (inside_tet) then
      
	  cell_count=cell_count+1
	  
          if (set_mesh_flag) then
          
  	    problem_volumes(volume_number)%cell_list(cell_count)%cell%i=ix
  	    problem_volumes(volume_number)%cell_list(cell_count)%cell%j=iy
  	    problem_volumes(volume_number)%cell_list(cell_count)%cell%k=iz
  	    problem_volumes(volume_number)%cell_list(cell_count)%point=centre
         
          end if
	  
        end if ! cell point is inside this tet
	
      end do ! next z cell
    end do ! next y cell
  end do ! next x cell          
  
!  CALL write_line('FINISHED: mesh_tet',0,output_to_screen_flag)

  RETURN
  
  END SUBROUTINE mesh_tet
