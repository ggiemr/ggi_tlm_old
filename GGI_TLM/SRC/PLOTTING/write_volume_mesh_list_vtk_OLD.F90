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
!
! NAME
!     SUBROUTINE write_volume_mesh_list_vtk
!
! DESCRIPTION
!     write_volume_mesh_list_vtk:
!
!     Write cell list to vtk format file
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 30/08/2012 CJS
!
!
SUBROUTINE write_volume_mesh_list_vtk(file_unit,number_of_cells,cell_list)

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

integer	:: file_unit
integer	:: number_of_cells

type(cell_point)	:: cell_list(1:number_of_cells)

! local variables

integer 	:: cell_number
integer 	:: number_of_points
integer 	:: point_number
type(xyz)	:: point1,point2,point3,point4,point5,point6,point7,point8

! START

  CALL write_line('CALLED: Write_volume_mesh_list',0,output_to_screen_flag)

  number_of_points=8*number_of_cells

  write(file_unit,'(A)')'# vtk DataFile Version 2.0'
  write(file_unit,'(A)')'trim(frame_filename)'
  write(file_unit,'(A)')'ASCII'
  write(file_unit,'(A)')'DATASET POLYDATA'
  write(file_unit,'(A,I10,A)')'POINTS',number_of_points,' float'

! write point data 

    point_number=0
    
    do cell_number=1,number_of_cells
          
      CALL get_cell_corner_coordinates(cell_list(cell_number),	&
           point1,point2,point3,point4,point5,point6,point7,point8)    
	
      write(file_unit,8000)point1%x,point1%y,point1%z
      write(file_unit,8000)point2%x,point2%y,point2%z
      write(file_unit,8000)point3%x,point3%y,point3%z
      write(file_unit,8000)point4%x,point4%y,point4%z
      write(file_unit,8000)point5%x,point5%y,point5%z
      write(file_unit,8000)point6%x,point6%y,point6%z
      write(file_unit,8000)point7%x,point7%y,point7%z
      write(file_unit,8000)point8%x,point8%y,point8%z

      point_number=point_number+8
                             
    end do! next cell 
      
8000  format(3E14.5)
  
! write cell data
    write(file_unit,'(A,2I10)')'POLYGONS',number_of_cells*6,number_of_cells*30

    point_number=0
    do cell_number=1,number_of_cells
    
      write(file_unit,8010)4,point_number  ,point_number+1,point_number+2,point_number+3  ! xmin face
      write(file_unit,8010)4,point_number+4,point_number+7,point_number+6,point_number+5  ! xmax face
      write(file_unit,8010)4,point_number  ,point_number+4,point_number+5,point_number+1  ! ymin face
      write(file_unit,8010)4,point_number+3,point_number+2,point_number+6,point_number+7  ! ymax face
      write(file_unit,8010)4,point_number  ,point_number+3,point_number+7,point_number+4  ! zmin face
      write(file_unit,8010)4,point_number+1,point_number+5,point_number+6,point_number+2  ! zmax face
      point_number=point_number+8

8010  format(I3,4I8)
      
    end do ! next cell

  CALL write_line('FINISHED: Write_volume_mesh_list',0,output_to_screen_flag)
  
  RETURN
  
END SUBROUTINE write_volume_mesh_list_vtk
