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
!SUBROUTINE build_volume_mesh
!
! NAME
!     SUBROUTINE build_volume_mesh
!
! DESCRIPTION
!     build_volume_mesh:
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 30/08/2012 CJS
!
!
SUBROUTINE build_volume_mesh()

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: volume_number

integer	:: number_of_tets
integer	:: tet_count

type(xyz)	:: point1,point2,point3,point4
integer		:: tet_number
integer		:: cell_count
logical		:: set_mesh_flag

! START

  CALL write_line('CALLED: build_volume_mesh',0,output_to_screen_flag)

  do volume_number=1,n_volumes

    number_of_tets=problem_volumes(volume_number)%number_of_tets

! Initially only count the number of cells in the volume mesh

    set_mesh_flag=.FALSE.
    cell_count=0
    do tet_number=1,number_of_tets    
    
      CALL mesh_tet(volume_number,tet_number,cell_count,set_mesh_flag)   
      
    end do ! next tet

! Allocate memory for the mesh

    problem_volumes(volume_number)%number_of_cells=cell_count
    
    if (cell_count.gt.0) then
    
      ALLOCATE( problem_volumes(volume_number)%cell_list(1:cell_count) )
    
! Now regenerate the volume mesh and set the cell_list

      set_mesh_flag=.TRUE.
      cell_count=0
      do tet_number=1,number_of_tets    
    
        CALL mesh_tet(volume_number,tet_number,cell_count,set_mesh_flag)  
	 
      end do ! next tet
      
    end if ! number of mesh cells .gt.0
    
  CALL write_line_integer('Number of cells',cell_count,0,output_to_screen_flag)   
    
  end do ! next volume number

  CALL write_line('FINISHED: build_volume_mesh',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE build_volume_mesh
