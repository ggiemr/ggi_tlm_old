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
! SUBROUTINE set_cell_update_codes
! SUBROUTINE set_face_update_codes
!
! NAME
!     set_cell_update_codes
!
! DESCRIPTION
!     set update codes on cells to indicate that the cell is in some way special
!     i.e. there is a material, excitation or output at the cell
!     Also count excitation and output cells and return the cell count in the arrays
!     local_cell_excitation and local_cell_output 
!
! COMMENTS
!
!
! HISTORY
!
!     started 15/08/2012 CJS
!
!
SUBROUTINE set_cell_update_codes

USE TLM_general
USE TLM_output
USE TLM_excitation
USE geometry
USE TLM_surface_materials
USE mesh

IMPLICIT NONE

! local variables

  integer cx,cy,cz
  
  integer cell_number
  integer special_cell_count


! START
  
  CALL write_line('CALLED: set_cell_update_codes',0,output_to_screen_flag)

  number_of_cell_centre_codes=nx*ny*nz
  
  ALLOCATE ( cell_centre_update_code(1:number_of_cell_centre_codes) )
  
  cell_centre_update_code(1:number_of_cell_centre_codes)=0

! Count special cells i.e. material, excitation or output cells  
! Also count the total number of output cells

  cell_number=0
  special_cell_count=0
  total_number_excitation_cells=0
  total_number_output_cells=0
  
  do cz=1,nz
    do cy=1,ny
      do cx=1,nx
            
        cell_number=cell_number+1
	
	if ( (local_cell_material(cx,cy,cz).NE.0).OR.	&
	     (local_cell_excitation(cx,cy,cz).NE.0).OR.	&
	     (local_cell_output(cx,cy,cz).NE.0) ) then
	     
	  special_cell_count=special_cell_count+1
	
	end if
	
	if (local_cell_excitation(cx,cy,cz).NE.0) then
	
	  total_number_excitation_cells=total_number_excitation_cells+1
	  
	end if
	
	if (local_cell_output(cx,cy,cz).NE.0) then
	
	  total_number_output_cells=total_number_output_cells+1
	  
	end if
	
      end do  ! next x cell
    end do    ! next y cell
  end do      ! next z cell

! Allocate memory for cell_code to material, excitation and output arrays 

  ALLOCATE ( cell_update_code_to_material_data(1:special_cell_count,1:2) )
  cell_update_code_to_material_data(1:special_cell_count,1:2)=0
  
  ALLOCATE ( cell_update_code_to_excitation_number(1:special_cell_count) )
  cell_update_code_to_excitation_number(1:special_cell_count)=0
  
  ALLOCATE ( cell_update_code_to_output_number(1:special_cell_count) )
  cell_update_code_to_output_number(1:special_cell_count)=0

! Allocate memory for all cell excitations and outputs 
  
  ALLOCATE( cell_excitation_field(1:total_number_excitation_cells,1:6) )
  
  ALLOCATE( cell_output_field(1:total_number_output_cells,1:6) )

! Set code lookup arrays on special cells 

  cell_number=0
  special_cell_count=0
  total_number_excitation_cells=0
  total_number_output_cells=0
  
  do cz=1,nz
    do cy=1,ny
      do cx=1,nx
            
        cell_number=cell_number+1
	
	if ( (local_cell_material(cx,cy,cz).NE.0).OR.	&
	     (local_cell_excitation(cx,cy,cz).NE.0).OR.	&
	     (local_cell_output(cx,cy,cz).NE.0) ) then
	     
	  special_cell_count=special_cell_count+1
	     
	  cell_centre_update_code(cell_number)=special_cell_count
	
	end if
	
	if (local_cell_material(cx,cy,cz).NE.0) then
	
	  cell_update_code_to_material_data(special_cell_count,1)=local_cell_material(cx,cy,cz)
	  
	end if

	
	if (local_cell_excitation(cx,cy,cz).NE.0) then
		
	  total_number_excitation_cells=total_number_excitation_cells+1
	  cell_update_code_to_excitation_number(special_cell_count)=1
	  local_cell_excitation(cx,cy,cz)=total_number_excitation_cells
	
	end if
	
	
	if (local_cell_output(cx,cy,cz).NE.0) then
	
	  total_number_output_cells=total_number_output_cells+1
	  cell_update_code_to_output_number(special_cell_count)=1
	  local_cell_output(cx,cy,cz)=total_number_output_cells
	  
	end if
		
      end do  ! next x cell
    end do    ! next y cell
  end do      ! next z cell
 
  CALL write_line_integer('Number of special cells',special_cell_count,0,output_to_screen_flag)
  CALL write_line_integer('Total number of excitation cells',total_number_excitation_cells,0,output_to_screen_flag)
  CALL write_line_integer('Total number of output cells',total_number_output_cells,0,output_to_screen_flag)
  
  n_special_cells=special_cell_count
  
  CALL write_line('FINISHED: set_cell_update_codes',0,output_to_screen_flag)


  RETURN

9000 CALL write_line('Error in set_cell_update_codes:',0,.TRUE.)
     STOP


END SUBROUTINE set_cell_update_codes
!
! NAME
!     set_face_update_codes
!
! DESCRIPTION
!     set update codes on cell faces to indicate that the face is in some way special
!     i.e. there is a material, excitation or output at the face
!     Also count excitation and output faces and return the face count in the arrays
!     local_surface_excitation and local_surface_output  
!     
! COMMENTS
!
!
! HISTORY
!
!     started 15/08/2012 CJS
!
!
SUBROUTINE set_face_update_codes

USE TLM_general
USE TLM_output
USE TLM_excitation
USE geometry
USE TLM_surface_materials
USE mesh

IMPLICIT NONE

! local variables

  integer	:: cx,cy,cz
  
  integer	:: face_number
  integer	:: special_face_count

! START
  
  CALL write_line('CALLED: set_face_update_codes',0,output_to_screen_flag)

  number_of_face_codes=(nx-1)*ny*nz+nx*(ny-1)*nz+nx*ny*(nz-1)

  ALLOCATE ( face_update_code(1:number_of_face_codes) )

! set face codes for surfaces  

  face_update_code(1:number_of_face_codes)=0  ! free space initially

! The following loops must mimic thse in the connect subroutine so as to
! get the correct face update code for each face.
! We only need to count the faces on the min side of the cell for materials because
! the opposite face gets set to -1*material_number 

! Count special cells i.e. material, excitation or output cells  
! Also count the total number of output cells

  face_number=0
  special_face_count=0
  total_number_excitation_faces=0
  total_number_output_faces=0
    
! faces normal to x i.e. the y z plane
  do cz=1,nz
    do cy=1,ny
      do cx=2,nx
            
        face_number=face_number+1	
	
	if ( (local_surface_material(cx,cy,cz,face_xmin).NE.0).OR.	&
	     (local_surface_excitation(cx,cy,cz,face_xmin).NE.0).OR.	&
	     (local_surface_output(cx,cy,cz,face_xmin).NE.0) ) then
	     
	  special_face_count=special_face_count+1
	
	end if
	
	if (local_surface_excitation(cx,cy,cz,face_xmin).NE.0) then
	
	  total_number_excitation_faces=total_number_excitation_faces+1
	  
	end if
	
	if (local_surface_output(cx,cy,cz,face_xmin).NE.0) then
	
	  total_number_output_faces=total_number_output_faces+1
	  
	end if

      end do  ! next x cell
    end do    ! next y cell
  end do      ! next z cell
  
! faces normal to y i.e. the x z plane
  do cz=1,nz
    do cy=2,ny
      do cx=1,nx
      
        face_number=face_number+1	
	
	if ( (local_surface_material(cx,cy,cz,face_ymin).NE.0).OR.	&
	     (local_surface_excitation(cx,cy,cz,face_ymin).NE.0).OR.	&
	     (local_surface_output(cx,cy,cz,face_ymin).NE.0) ) then
	     
	  special_face_count=special_face_count+1
	
	end if
	
	if (local_surface_excitation(cx,cy,cz,face_ymin).NE.0) then
	
	  total_number_excitation_faces=total_number_excitation_faces+1
	  
	end if
	
	if (local_surface_output(cx,cy,cz,face_ymin).NE.0) then
	
	  total_number_output_faces=total_number_output_faces+1
	  
	end if
	            
      end do  ! next x cell
    end do    ! next y cell
  end do      ! next z cell
  
! faces normal to z i.e. the x y plane
  do cz=2,nz
    do cy=1,ny
      do cx=1,nx
      	
        face_number=face_number+1	
	
	if ( (local_surface_material(cx,cy,cz,face_zmin).NE.0).OR.	&
	     (local_surface_excitation(cx,cy,cz,face_zmin).NE.0).OR.	&
	     (local_surface_output(cx,cy,cz,face_zmin).NE.0) ) then
	     
	  special_face_count=special_face_count+1
	
	end if
	
	if (local_surface_excitation(cx,cy,cz,face_zmin).NE.0) then
	
	  total_number_excitation_faces=total_number_excitation_faces+1
	  
	end if
	
	if (local_surface_output(cx,cy,cz,face_zmin).NE.0) then
	
	  total_number_output_faces=total_number_output_faces+1
	  
	end if
            
      end do  ! next x cell
    end do    ! next y cell
  end do      ! next z cell

! Allocate memory for face_code to material, excitation and output arrays 

  ALLOCATE ( face_update_code_to_material_data(1:special_face_count,1:2) )
  face_update_code_to_material_data(1:special_face_count,1:2)=0
  
  ALLOCATE ( face_update_code_to_excitation_number(1:special_face_count) )
  face_update_code_to_excitation_number(1:special_face_count)=0
  
  ALLOCATE ( face_update_code_to_output_number(1:special_face_count) )
  face_update_code_to_output_number(1:special_face_count)=0

! Allocate memory for all face excitations and outputs 
  
  write(*,*)'Allocating ',total_number_excitation_faces,' face excitation fields'
  ALLOCATE( face_excitation_field(1:total_number_excitation_faces,1:2,1:6) )
  
  write(*,*)'Allocating ',total_number_output_faces,' face output fields'
  ALLOCATE( face_output_field(1:total_number_output_faces,1:2,1:6) )

! Set code lookup arrays on special cells 

  face_number=0
  special_face_count=0
  total_number_excitation_faces=0
  total_number_output_faces=0
    
! faces normal to x i.e. the y z plane
  do cz=1,nz
    do cy=1,ny
      do cx=2,nx
            
        face_number=face_number+1	
	
	if ( (local_surface_material(cx,cy,cz,face_xmin).NE.0).OR.	&
	     (local_surface_excitation(cx,cy,cz,face_xmin).NE.0).OR.	&
	     (local_surface_output(cx,cy,cz,face_xmin).NE.0) ) then
	     
	  special_face_count=special_face_count+1
			     
	  face_update_code(face_number)=special_face_count
	
	end if
	
	if (local_surface_material(cx,cy,cz,face_xmin).NE.0) then
	
	  face_update_code_to_material_data(special_face_count,1)=local_surface_material(cx,cy,cz,face_xmin)
	  
	end if
	
	if (local_surface_excitation(cx,cy,cz,face_xmin).NE.0) then
	
	  total_number_excitation_faces=total_number_excitation_faces+1
	  local_surface_excitation(cx,cy,cz,face_xmin)=total_number_excitation_faces
	  face_update_code_to_excitation_number(special_face_count)=1
	  
	end if
	
	if (local_surface_output(cx,cy,cz,face_xmin).NE.0) then
	
	  total_number_output_faces=total_number_output_faces+1
	  local_surface_output(cx,cy,cz,face_xmin)=total_number_output_faces
	  face_update_code_to_output_number(special_face_count)=1
	  
	end if
		
      end do  ! next x cell
    end do    ! next y cell
  end do      ! next z cell
  
! faces normal to y i.e. the x z plane
  do cz=1,nz
    do cy=2,ny
      do cx=1,nx
      
        face_number=face_number+1	
	
	if ( (local_surface_material(cx,cy,cz,face_ymin).NE.0).OR.	&
	     (local_surface_excitation(cx,cy,cz,face_ymin).NE.0).OR.	&
	     (local_surface_output(cx,cy,cz,face_ymin).NE.0) ) then
	     
	  special_face_count=special_face_count+1
			     
	  face_update_code(face_number)=special_face_count
	
	end if
	
	if (local_surface_material(cx,cy,cz,face_ymin).NE.0) then
	
	  face_update_code_to_material_data(special_face_count,1)=local_surface_material(cx,cy,cz,face_ymin)
	  
	end if
	
	if (local_surface_excitation(cx,cy,cz,face_ymin).NE.0) then
	
	  total_number_excitation_faces=total_number_excitation_faces+1
	  local_surface_excitation(cx,cy,cz,face_ymin)=total_number_excitation_faces
	  face_update_code_to_excitation_number(special_face_count)=1
	  
	end if
	
	if (local_surface_output(cx,cy,cz,face_ymin).NE.0) then
	
	  total_number_output_faces=total_number_output_faces+1
	  local_surface_output(cx,cy,cz,face_ymin)=total_number_output_faces
	  face_update_code_to_output_number(special_face_count)=1
	  
	end if
	            
      end do  ! next x cell
    end do    ! next y cell
  end do      ! next z cell
  
! faces normal to z i.e. the x y plane
  do cz=2,nz
    do cy=1,ny
      do cx=1,nx
      	
        face_number=face_number+1	
	
	if ( (local_surface_material(cx,cy,cz,face_zmin).NE.0).OR.	&
	     (local_surface_excitation(cx,cy,cz,face_zmin).NE.0).OR.	&
	     (local_surface_output(cx,cy,cz,face_zmin).NE.0) ) then
	     
	  special_face_count=special_face_count+1
			     
	  face_update_code(face_number)=special_face_count
	
	end if
	
	if (local_surface_material(cx,cy,cz,face_zmin).NE.0) then
	
	  face_update_code_to_material_data(special_face_count,1)=local_surface_material(cx,cy,cz,face_zmin)
	  
	end if
	
	if (local_surface_excitation(cx,cy,cz,face_zmin).NE.0) then
	
	  total_number_excitation_faces=total_number_excitation_faces+1
	  local_surface_excitation(cx,cy,cz,face_zmin)=total_number_excitation_faces
	  face_update_code_to_excitation_number(special_face_count)=1
	  
	end if
	
	if (local_surface_output(cx,cy,cz,face_zmin).NE.0) then
	
	  total_number_output_faces=total_number_output_faces+1
	  local_surface_output(cx,cy,cz,face_zmin)=total_number_output_faces
	  face_update_code_to_output_number(special_face_count)=1
	  
	end if
            
      end do  ! next x cell
    end do    ! next y cell
  end do      ! next z cell
 
  CALL write_line_integer('Number of special faces',special_face_count,0,output_to_screen_flag)
  CALL write_line_integer('Total number of excitation faces',total_number_excitation_faces,0,output_to_screen_flag)
  CALL write_line_integer('Total number of output faces',total_number_output_faces,0,output_to_screen_flag)

  n_special_faces=special_face_count
  
  CALL write_line('FINISHED: set_face_update_codes',0,output_to_screen_flag)


  RETURN
  
9000 CALL write_line('Error in set_face_update_codes:',0,.TRUE.)
     CALL write_line('number_of_face_codes.NE.counted faces:',0,.TRUE.)
     CALL write_line_integer('number_of_face_codes',number_of_face_codes,0,.TRUE.)
     CALL write_line_integer('counted faces       ',face_number,0,.TRUE.)
     STOP

END SUBROUTINE set_face_update_codes
