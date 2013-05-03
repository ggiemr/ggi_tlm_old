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
!SUBROUTINE build_surface_mesh
!
! NAME
!     SUBROUTINE build_surface_mesh
!
! DESCRIPTION
!     build_surface_mesh:
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 8/08/2012 CJS
!
!
SUBROUTINE build_surface_mesh()

USE TLM_general
USE geometry_types
USE geometry
USE file_information
USE constants

IMPLICIT NONE

! local variables

integer	:: surface_number

integer	:: number_of_triangles
integer	:: triangle_count

integer		:: triangle_number
integer		:: face_count
logical		:: set_mesh_flag

! START

  CALL write_line('CALLED: build_surface_mesh',0,output_to_screen_flag)

  do surface_number=1,n_surfaces

    number_of_triangles=problem_surfaces(surface_number)%number_of_triangles

! Initially only count the number of faces in the surface mesh

    set_mesh_flag=.FALSE.
    face_count=0
    do triangle_number=1,number_of_triangles    
    
      CALL mesh_triangle(surface_number,triangle_number,face_count,set_mesh_flag)   
      
    end do ! next triangle

! Allocate memory for the mesh

    problem_surfaces(surface_number)%number_of_faces=face_count
    
    if (face_count.gt.0) then
    
      ALLOCATE( problem_surfaces(surface_number)%face_list(1:face_count) )
    
! Now regenerate the surface mesh and set the face_list

      set_mesh_flag=.TRUE.
      face_count=0
      do triangle_number=1,number_of_triangles    
      
        CALL mesh_triangle(surface_number,triangle_number,face_count,set_mesh_flag)   
	
      end do ! next triangle
      
    end if ! number of mesh faces .gt.0
    
  end do ! next surface number

  CALL write_line('FINISHED: build_surface_mesh',0,output_to_screen_flag)
  
  RETURN
  
  
END SUBROUTINE build_surface_mesh
