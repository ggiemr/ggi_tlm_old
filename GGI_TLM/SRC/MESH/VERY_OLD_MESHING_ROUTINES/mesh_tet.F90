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
!     mesh a single tet by first meshing the surfaces of the 4 surface triangles then filling the 
!     enclosed cells
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

  integer 		:: ix,iy,iz,cell_face
  
  real*8		:: x1,y1,z1
  real*8		:: x2,y2,z2
  real*8		:: x3,y3,z3
  real*8		:: x4,y4,z4
    
  integer		:: mat
     
  integer,allocatable 	:: local_face_list(:,:)
  integer,allocatable 	:: local_face_list1(:,:)
  integer,allocatable 	:: local_face_list2(:,:)
  integer,allocatable 	:: local_face_list3(:,:)
  integer,allocatable 	:: local_face_list4(:,:)
  
  integer 		:: face,n_faces,max_n_faces
  integer		:: n_faces1,n_faces2,n_faces3,n_faces4
  
  integer 		:: ix_min,iy_min,iz_min,cell_face_min
  integer 		:: ix_max,iy_max,iz_max,cell_face_max
  
  integer		:: n_faces_set
  integer		:: iz_column(2),iz_face(2)
  integer		:: iswap
  
! START

!  CALL write_line('CALLED: mesh_tet',0,output_to_screen_flag)

  max_n_faces=nx*ny+nx*nz+ny*nz
  
  ALLOCATE( local_face_list(1:max_n_faces,1:4) )
  ALLOCATE( local_face_list1(1:max_n_faces,1:4) )
  ALLOCATE( local_face_list2(1:max_n_faces,1:4) )
  ALLOCATE( local_face_list3(1:max_n_faces,1:4) )
  ALLOCATE( local_face_list4(1:max_n_faces,1:4) )
  
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

! mesh surface triangles      
  CALL get_face_list(local_face_list1,max_n_faces,n_faces1,x1,y1,z1,x2,y2,z2,x3,y3,z3)
  CALL get_face_list(local_face_list2,max_n_faces,n_faces2,x1,y1,z1,x4,y4,z4,x2,y2,z2)
  CALL get_face_list(local_face_list3,max_n_faces,n_faces3,x2,y2,z2,x4,y4,z4,x3,y3,z3)
  CALL get_face_list(local_face_list4,max_n_faces,n_faces4,x1,y1,z1,x3,y3,z3,x4,y4,z4)

! combine face lists

  n_faces=0
  do face=1,n_faces1	
    n_faces=n_faces+1
    local_face_list(n_faces,1)=local_face_list1(face,1)
    local_face_list(n_faces,2)=local_face_list1(face,2)
    local_face_list(n_faces,3)=local_face_list1(face,3)
    local_face_list(n_faces,4)=local_face_list1(face,4)        
  end do ! next face in list
  do face=1,n_faces2 
    n_faces=n_faces+1
    local_face_list(n_faces,1)=local_face_list2(face,1)
    local_face_list(n_faces,2)=local_face_list2(face,2)
    local_face_list(n_faces,3)=local_face_list2(face,3)
    local_face_list(n_faces,4)=local_face_list2(face,4)        
  end do ! next face in list
  do face=1,n_faces3  
    n_faces=n_faces+1
    local_face_list(n_faces,1)=local_face_list3(face,1)
    local_face_list(n_faces,2)=local_face_list3(face,2)
    local_face_list(n_faces,3)=local_face_list3(face,3)
    local_face_list(n_faces,4)=local_face_list3(face,4)        
  end do ! next face in list
  do face=1,n_faces4   
    n_faces=n_faces+1
    local_face_list(n_faces,1)=local_face_list4(face,1)
    local_face_list(n_faces,2)=local_face_list4(face,2)
    local_face_list(n_faces,3)=local_face_list4(face,3)
    local_face_list(n_faces,4)=local_face_list4(face,4)        
  end do ! next face in list

  if (n_faces.lt.1) then
    cell_count=0
    RETURN
  end if

! get the extent of the volume containing the tet cells

  ix_min=local_face_list(1,1)
  iy_min=local_face_list(1,2)
  iz_min=local_face_list(1,3)
  ix_max=local_face_list(1,1)
  iy_max=local_face_list(1,2)
  iz_max=local_face_list(1,3)
  
  do face=2,n_faces  
    ix_min=min(local_face_list(face,1),ix_min)
    iy_min=min(local_face_list(face,2),iy_min)
    iz_min=min(local_face_list(face,3),iz_min)
    ix_max=max(local_face_list(face,1),ix_max)
    iy_max=max(local_face_list(face,2),iy_max)
    iz_max=max(local_face_list(face,3),iz_max)
  end do ! next face in list
  
! loop over xy extent of the enclosing volume

  do ix=ix_min,ix_max
    do iy=iy_min,iy_max
    
! look for two faces defining the extent in z 
      n_faces_set=0
      do face=1,n_faces  

        if ( (ix.eq.local_face_list(face,1)).AND.(iy.eq.local_face_list(face,2)) ) then 
! we have a cell with the current ix and iy coordinate


          if ( (local_face_list(face,4).EQ.face_zmin).OR.(local_face_list(face,4).EQ.face_zmax) ) then 
! face normal to z
  	    n_faces_set=n_faces_set+1
            
            if (local_face_list(face,4).eq.face_zmax) then
              iz_column(n_faces_set)=local_face_list(face,3)+1
              iz_face(n_faces_set)=face_zmin
            else
              iz_column(n_faces_set)=local_face_list(face,3)
              iz_face(n_faces_set)=face_zmin	  
            end if
            
  	  end if  ! face normal to z
	  
  	end if ! cell with the current ix and iy coordinate

      end do ! next face in list

      if (n_faces_set.eq.2) then ! set cells in this column

        if (iz_column(1).gt.iz_column(2)) then
          iswap=iz_column(1)
          iz_column(1)=iz_column(2)
          iz_column(2)=iswap
        end if
        
  	do iz=iz_column(1),iz_column(2)-1
       
          cell_count=cell_count+1
         
          if (set_mesh_flag) then
          
  	    problem_volumes(volume_number)%cell_list(cell_count)%cell%i=ix
  	    problem_volumes(volume_number)%cell_list(cell_count)%cell%j=iy
  	    problem_volumes(volume_number)%cell_list(cell_count)%cell%k=iz
  	    problem_volumes(volume_number)%cell_list(cell_count)%point=centre
         
          end if
         
        end do

      else if (n_faces_set.ne.0) then

  	write(*,*)'Error in mesh_tet'
        write(*,*)'n_faces_set=',n_faces_set
	STOP

      end if
    
    end do
  end do
    
  if ( allocated( local_face_list ) ) DEALLOCATE( local_face_list1 )
  if ( allocated( local_face_list1 ) ) DEALLOCATE( local_face_list1 )
  if ( allocated( local_face_list2 ) ) DEALLOCATE( local_face_list2 )
  if ( allocated( local_face_list3 ) ) DEALLOCATE( local_face_list3 )
  if ( allocated( local_face_list4 ) ) DEALLOCATE( local_face_list4 )
  
!  CALL write_line('FINISHED: mesh_tet',0,output_to_screen_flag)

  RETURN
  
  END SUBROUTINE mesh_tet
