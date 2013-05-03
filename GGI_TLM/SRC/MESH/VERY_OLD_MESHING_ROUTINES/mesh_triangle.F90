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
!SUBROUTINE mesh_triangle
!
! NAME
!     SUBROUTINE mesh_triangle
!
! DESCRIPTION
!     mesh_triangle:
!
!     mesh a single triangle
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 9/08/2012 CJS
!
!
  SUBROUTINE mesh_triangle(surface_number,triangle_number,face_count,set_mesh_flag)
  
USE TLM_general
USE geometry_types
USE geometry
USE cell_parameters
USE file_information
USE constants

IMPLICIT NONE

integer				:: surface_number
integer				:: triangle_number
integer				:: face_count
logical				:: set_mesh_flag

! local variables

  integer 		:: ix,iy,iz,cell_face
  
  real*8		:: x1,y1,z1
  real*8		:: x2,y2,z2
  real*8		:: x3,y3,z3
  real*8		:: xn,yn,zn
    
  integer		:: mat
     
  integer,allocatable 	:: local_face_list(:,:)
  
  integer 		:: face,n_faces,max_n_faces

! START

  max_n_faces=nx*ny+nx*nz+ny*nz
  
  ALLOCATE( local_face_list(1:max_n_faces,1:4) )
  
  x1=problem_surfaces(surface_number)%triangle_list(triangle_number)%vertex(1)%x
  y1=problem_surfaces(surface_number)%triangle_list(triangle_number)%vertex(1)%y
  z1=problem_surfaces(surface_number)%triangle_list(triangle_number)%vertex(1)%z
  
  x2=problem_surfaces(surface_number)%triangle_list(triangle_number)%vertex(2)%x
  y2=problem_surfaces(surface_number)%triangle_list(triangle_number)%vertex(2)%y
  z2=problem_surfaces(surface_number)%triangle_list(triangle_number)%vertex(2)%z
  
  x3=problem_surfaces(surface_number)%triangle_list(triangle_number)%vertex(3)%x
  y3=problem_surfaces(surface_number)%triangle_list(triangle_number)%vertex(3)%y
  z3=problem_surfaces(surface_number)%triangle_list(triangle_number)%vertex(3)%z
      
  CALL get_face_list(local_face_list,max_n_faces,n_faces,x1,y1,z1,x2,y2,z2,x3,y3,z3)
  
  if (n_faces.ne.0) then 
! there are some meshed faces
  
    CALL triangle_normal(x1,y1,z1,x2,y2,z2,x3,y3,z3,xn,yn,zn)
  
    if (set_mesh_flag) then
  
      do face=1,n_faces

! set the face accordning to the face list         
        problem_surfaces(surface_number)%face_list(face_count+face)%cell%i=local_face_list(face,1)
        problem_surfaces(surface_number)%face_list(face_count+face)%cell%j=local_face_list(face,2)
        problem_surfaces(surface_number)%face_list(face_count+face)%cell%k=local_face_list(face,3)
        problem_surfaces(surface_number)%face_list(face_count+face)%point=local_face_list(face,4)
     
! look at the triangle normal and ensure that the meshed face has the correct normal     

! FACES NORMAL TO X  
        if (xn.gt.0d0) then
! triangle normal is in +x direction, we must make sure that the meshed face is face_xmax

          if (problem_surfaces(surface_number)%face_list(face_count+face)%point.eq.face_xmin) then
! must change the face to face_xmax of the adjacent cell

            problem_surfaces(surface_number)%face_list(face_count+face)%point=face_xmax
            problem_surfaces(surface_number)%face_list(face_count+face)%cell%i=	&
	    problem_surfaces(surface_number)%face_list(face_count+face)%cell%i-1
	    
	  end if ! face=face_xmin
     
        end if ! normal in +x direction
  
        if (xn.lt.0d0) then
! triangle normal is in -x direction, we must make sure that the meshed face is face_xmin

          if (problem_surfaces(surface_number)%face_list(face_count+face)%point.eq.face_xmax) then
! must change the face to face_xmin of the adjacent cell

            problem_surfaces(surface_number)%face_list(face_count+face)%point=face_xmin
            problem_surfaces(surface_number)%face_list(face_count+face)%cell%i=	&
	    problem_surfaces(surface_number)%face_list(face_count+face)%cell%i+1
	    
	  end if ! face=face_xmax
       
        end if ! normal in -x direction
  
! FACES NORMAL TO Y
        if (yn.gt.0d0) then
! triangle normal is in +y direction, we must make sure that the meshed face is face_ymax

          if (problem_surfaces(surface_number)%face_list(face_count+face)%point.eq.face_ymin) then
! must change the face to face_ymax of the adjacent cell

            problem_surfaces(surface_number)%face_list(face_count+face)%point=face_ymax
            problem_surfaces(surface_number)%face_list(face_count+face)%cell%j=	&
	    problem_surfaces(surface_number)%face_list(face_count+face)%cell%j-1
	    
	  end if ! face=face_ymin
     
        end if ! normal in +y direction
  
        if (yn.lt.0d0) then
! triangle normal is in -y direction, we must make sure that the meshed face is face_ymin

          if (problem_surfaces(surface_number)%face_list(face_count+face)%point.eq.face_ymax) then
! must change the face to face_ymin of the adjacent cell

            problem_surfaces(surface_number)%face_list(face_count+face)%point=face_ymin
            problem_surfaces(surface_number)%face_list(face_count+face)%cell%j=	&
	    problem_surfaces(surface_number)%face_list(face_count+face)%cell%j+1
	    
	  end if ! face=face_ymax
     
        end if ! normal in -y direction
  
! FACES NORMAL TO Z
        if (zn.gt.0d0) then
! triangle normal is in +z direction, we must make sure that the meshed face is face_zmax

          if (problem_surfaces(surface_number)%face_list(face_count+face)%point.eq.face_zmin) then
! must change the face to face_zmax of the adjacent cell

            problem_surfaces(surface_number)%face_list(face_count+face)%point=face_zmax
            problem_surfaces(surface_number)%face_list(face_count+face)%cell%k=	&
	    problem_surfaces(surface_number)%face_list(face_count+face)%cell%k-1
	    
	  end if ! face=face_zmin
     
        end if ! normal in +z direction
  
        if (zn.lt.0d0) then
! triangle normal is in -z direction, we must make sure that the meshed face is face_zmin

          if (problem_surfaces(surface_number)%face_list(face_count+face)%point.eq.face_zmax) then
! must change the face to face_zmin of the adjacent cell

            problem_surfaces(surface_number)%face_list(face_count+face)%point=face_zmin
            problem_surfaces(surface_number)%face_list(face_count+face)%cell%k=	&
	    problem_surfaces(surface_number)%face_list(face_count+face)%cell%k+1
	    
	  end if ! face=face_zmax
     
        end if ! normal in -z direction
     
      end do ! next face in list

    end if  ! set_mesh_flag 
	
    face_count=face_count+n_faces
  
  end if !n_faces.ne.0 
  
  if ( allocated( local_face_list ) ) DEALLOCATE( local_face_list )
  
!  write(*,*)'FINISHED: mesh_triangle',n_faces,face_count

  RETURN
  
  END SUBROUTINE mesh_triangle
