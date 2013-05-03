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
!SUBROUTINE get_cell_centre_coordinate(cell,point)
!SUBROUTINE get_cell_point_coordinate(face,point)
!SUBROUTINE get_cell_face_corner_coordinates(face,point1,point2,point3,point4)
!SUBROUTINE get_cell_corner_coordinates(cell,point1,point2,point3,point4,point5,point6,point7,point8)
!SUBROUTINE point_to_cell(point,cell)
!SUBROUTINE get_other_side_of_face(face1,face2)
!
! Name get_cell_centre_coordinate
!     
!
! Description
!     given the i,j,k cell number return the cell centre x,y,z coordinate
!
! Comments:
!      
!
! History
!
!     started 9/08/12 CJS
!

SUBROUTINE get_cell_centre_coordinate(cell,point)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(ijk)	:: cell
type(xyz)	:: point

! local_variables

! START

  point%x=mesh_xmin+cell%i*dl-dl/2d0
  point%y=mesh_ymin+cell%j*dl-dl/2d0
  point%z=mesh_zmin+cell%k*dl-dl/2d0

  RETURN
  
END SUBROUTINE get_cell_centre_coordinate
!
! Name get_cell_point_coordinate
!     
!
! Description
!     given the cell point number return the corresponding x,y,z coordinate
!
! Comments:
!      
!
! History
!
!     started 9/08/12 CJS
!

SUBROUTINE get_cell_point_coordinate(local_cell_point,xyz_point)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(cell_point)	:: local_cell_point
type(xyz)		:: xyz_point

! local_variables

real*8	:: offset_x=0d0
real*8	:: offset_y=0d0
real*8	:: offset_z=0d0

! START

  if (local_cell_point%point.eq.face_xmin) then
    offset_x=-1.0d0
  else if (local_cell_point%point.eq.face_xmax)  then
    offset_x= 0.0d0
  else if (local_cell_point%point.eq.face_ymin)  then
    offset_y=-1.0d0
  else if (local_cell_point%point.eq.face_ymax)  then
    offset_y= 0.0d0
  else if (local_cell_point%point.eq.face_zmin)  then
    offset_z=-1.0d0
  else if (local_cell_point%point.eq.face_zmax)  then
    offset_z= 0.0d0
  else if (local_cell_point%point.eq.centre)  then
    offset_x= -0.5d0
    offset_y= -0.5d0
    offset_z= -0.5d0
  else 
    write(*,*)'Error in get_local_cell_point_coordinates'
    write(*,*)'No face defined'
    write(*,*)'Face number (local_cell_point%point)=',local_cell_point%point
    STOP
  end if

  xyz_point%x=mesh_xmin+local_cell_point%cell%i*dl+dl*offset_x
  xyz_point%y=mesh_ymin+local_cell_point%cell%j*dl+dl*offset_y
  xyz_point%z=mesh_zmin+local_cell_point%cell%k*dl+dl*offset_z

  RETURN
  
END SUBROUTINE get_cell_point_coordinate
!
! Name get_cell_face_corner_coordinates
!     
!
! Description
!     given the face number return the 4 face corner x,y,z coordinates
!     note: the order of the points returned defines the normal direction
!
! Comments:
!      
!
! History
!
!     started 9/08/12 CJS
!

SUBROUTINE get_cell_face_corner_coordinates(face,point1,point2,point3,point4)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(cell_point)	:: face
type(xyz)		:: point1,point2,point3,point4

! local_variables

! START
 
  if (face%point.eq.face_xmin) then
  
    point1%x=mesh_xmin+face%cell%i*dl-dl
    point1%y=mesh_ymin+face%cell%j*dl-dl
    point1%z=mesh_zmin+face%cell%k*dl-dl
    point2%x=mesh_xmin+face%cell%i*dl-dl
    point2%y=mesh_ymin+face%cell%j*dl-dl
    point2%z=mesh_zmin+face%cell%k*dl
    point3%x=mesh_xmin+face%cell%i*dl-dl
    point3%y=mesh_ymin+face%cell%j*dl
    point3%z=mesh_zmin+face%cell%k*dl
    point4%x=mesh_xmin+face%cell%i*dl-dl
    point4%y=mesh_ymin+face%cell%j*dl
    point4%z=mesh_zmin+face%cell%k*dl-dl

  else if (face%point.eq.face_xmax)  then
  
    point1%x=mesh_xmin+face%cell%i*dl
    point1%y=mesh_ymin+face%cell%j*dl-dl
    point1%z=mesh_zmin+face%cell%k*dl-dl
    point4%x=mesh_xmin+face%cell%i*dl
    point4%y=mesh_ymin+face%cell%j*dl-dl
    point4%z=mesh_zmin+face%cell%k*dl
    point3%x=mesh_xmin+face%cell%i*dl
    point3%y=mesh_ymin+face%cell%j*dl
    point3%z=mesh_zmin+face%cell%k*dl
    point2%x=mesh_xmin+face%cell%i*dl
    point2%y=mesh_ymin+face%cell%j*dl
    point2%z=mesh_zmin+face%cell%k*dl-dl
    
  else if (face%point.eq.face_ymin)  then
  
    point1%x=mesh_xmin+face%cell%i*dl-dl
    point1%y=mesh_ymin+face%cell%j*dl-dl
    point1%z=mesh_zmin+face%cell%k*dl-dl
    point2%x=mesh_xmin+face%cell%i*dl
    point2%y=mesh_ymin+face%cell%j*dl-dl
    point2%z=mesh_zmin+face%cell%k*dl-dl
    point3%x=mesh_xmin+face%cell%i*dl
    point3%y=mesh_ymin+face%cell%j*dl-dl
    point3%z=mesh_zmin+face%cell%k*dl
    point4%x=mesh_xmin+face%cell%i*dl-dl
    point4%y=mesh_ymin+face%cell%j*dl-dl
    point4%z=mesh_zmin+face%cell%k*dl
    
  else if (face%point.eq.face_ymax)  then
  
    point1%x=mesh_xmin+face%cell%i*dl-dl
    point1%y=mesh_ymin+face%cell%j*dl
    point1%z=mesh_zmin+face%cell%k*dl-dl
    point4%x=mesh_xmin+face%cell%i*dl
    point4%y=mesh_ymin+face%cell%j*dl
    point4%z=mesh_zmin+face%cell%k*dl-dl
    point3%x=mesh_xmin+face%cell%i*dl
    point3%y=mesh_ymin+face%cell%j*dl
    point3%z=mesh_zmin+face%cell%k*dl
    point2%x=mesh_xmin+face%cell%i*dl-dl
    point2%y=mesh_ymin+face%cell%j*dl
    point2%z=mesh_zmin+face%cell%k*dl
    
  else if (face%point.eq.face_zmin)  then
  
    point1%x=mesh_xmin+face%cell%i*dl-dl
    point1%y=mesh_ymin+face%cell%j*dl-dl
    point1%z=mesh_zmin+face%cell%k*dl-dl
    point2%x=mesh_xmin+face%cell%i*dl-dl
    point2%y=mesh_ymin+face%cell%j*dl
    point2%z=mesh_zmin+face%cell%k*dl-dl
    point3%x=mesh_xmin+face%cell%i*dl
    point3%y=mesh_ymin+face%cell%j*dl
    point3%z=mesh_zmin+face%cell%k*dl-dl
    point4%x=mesh_xmin+face%cell%i*dl
    point4%y=mesh_ymin+face%cell%j*dl-dl
    point4%z=mesh_zmin+face%cell%k*dl-dl
    
  else if (face%point.eq.face_zmax)  then
  
    point1%x=mesh_xmin+face%cell%i*dl-dl
    point1%y=mesh_ymin+face%cell%j*dl-dl
    point1%z=mesh_zmin+face%cell%k*dl
    point4%x=mesh_xmin+face%cell%i*dl-dl
    point4%y=mesh_ymin+face%cell%j*dl
    point4%z=mesh_zmin+face%cell%k*dl
    point3%x=mesh_xmin+face%cell%i*dl
    point3%y=mesh_ymin+face%cell%j*dl
    point3%z=mesh_zmin+face%cell%k*dl
    point2%x=mesh_xmin+face%cell%i*dl
    point2%y=mesh_ymin+face%cell%j*dl-dl
    point2%z=mesh_zmin+face%cell%k*dl
    
  else 
    write(*,*)'Error in get_cell_face_corner_coordinates'
    write(*,*)'No face defined'
    write(*,*)'face%point=',face%point
    STOP
  end if

  RETURN
  
END SUBROUTINE get_cell_face_corner_coordinates
!
! Name get_cell_corner_coordinates
!     
!
! Description
!     given the cell number return the 8 face corner x,y,z coordinates
!
! Comments:
!      
!
! History
!
!     started 30/08/12 CJS
!

SUBROUTINE get_cell_corner_coordinates(cell,point1,point2,point3,point4,point5,point6,point7,point8)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(ijk)	:: cell
type(xyz)	:: point1,point2,point3,point4,point5,point6,point7,point8

! local_variables

! START
 
    point1%x=mesh_xmin+cell%i*dl-dl
    point1%y=mesh_ymin+cell%j*dl-dl
    point1%z=mesh_zmin+cell%k*dl-dl
    point2%x=mesh_xmin+cell%i*dl-dl
    point2%y=mesh_ymin+cell%j*dl-dl
    point2%z=mesh_zmin+cell%k*dl
    point3%x=mesh_xmin+cell%i*dl-dl
    point3%y=mesh_ymin+cell%j*dl
    point3%z=mesh_zmin+cell%k*dl
    point4%x=mesh_xmin+cell%i*dl-dl
    point4%y=mesh_ymin+cell%j*dl
    point4%z=mesh_zmin+cell%k*dl-dl
  
    point5%x=mesh_xmin+cell%i*dl
    point5%y=mesh_ymin+cell%j*dl-dl
    point5%z=mesh_zmin+cell%k*dl-dl
    point6%x=mesh_xmin+cell%i*dl
    point6%y=mesh_ymin+cell%j*dl-dl
    point6%z=mesh_zmin+cell%k*dl
    point7%x=mesh_xmin+cell%i*dl
    point7%y=mesh_ymin+cell%j*dl
    point7%z=mesh_zmin+cell%k*dl
    point8%x=mesh_xmin+cell%i*dl
    point8%y=mesh_ymin+cell%j*dl
    point8%z=mesh_zmin+cell%k*dl-dl

  RETURN
  
END SUBROUTINE get_cell_corner_coordinates
!
! Name point_to_cell
!     
!
! Description
!      given the x,y,z coordinate return the i,j,k cell containing the point
!
! Comments:
!      
!
! History
!
!     started 10/08/12 CJS
!

SUBROUTINE point_to_cell(point,cell)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(xyz)	:: point
type(ijk)	:: cell

! local_variables

! START

  cell%i=int((point%x-mesh_xmin)/dl)+1
  cell%j=int((point%y-mesh_ymin)/dl)+1
  cell%k=int((point%z-mesh_zmin)/dl)+1
  
  if ((cell%i.lt.1).OR.(cell%i.gt.nx)) goto 9000
  if ((cell%j.lt.1).OR.(cell%j.gt.ny)) goto 9000
  if ((cell%k.lt.1).OR.(cell%k.gt.nz)) goto 9000
  
  return

9000 write(*,*)'Error in point_to_cell'
     write(*,*)'The cell is outside the defined mesh'
     write(*,*)'x =',point%x,' y =',point%y,' z =',point%z
     write(*,*)'i =',cell%i ,' j =',cell%j ,' k =',cell%k
     write(*,*)'nx=',point%x,' ny=',point%y,' nz=',point%z
     STOP

  RETURN

  RETURN
  
END SUBROUTINE point_to_cell
!
! Name get_other_side_of_face
!     
!
! Description
!     given the face number return the face centre x,y,z coordinate
!
! Comments:
!      
!
! History
!
!     started 9/08/12 CJS
!

SUBROUTINE get_other_side_of_face(face1,face2)

USE geometry_types
USE mesh

IMPLICIT NONE

! variables passed to subroutine

type(cell_point)	:: face1
type(cell_point)	:: face2

! local_variables

integer	:: ix,iy,iz

! START

  face2%cell%i=face1%cell%i
  face2%cell%j=face1%cell%j
  face2%cell%k=face1%cell%k

  if      (face1%point.eq.face_xmin) then
    face2%cell%i=face1%cell%i-1
    face2%point=face_xmax
  else if (face1%point.eq.face_xmax)  then
    face2%cell%i=face1%cell%i+1
    face2%point=face_xmin
  else if (face1%point.eq.face_ymin)  then
    face2%cell%j=face1%cell%j-1
    face2%point=face_ymax 
  else if (face1%point.eq.face_ymax)  then
    face2%cell%j=face1%cell%j+1
    face2%point=face_ymin   
  else if (face1%point.eq.face_zmin)  then
    face2%cell%k=face1%cell%k-1
    face2%point=face_zmax 
  else if (face1%point.eq.face_zmax)  then
    face2%cell%k=face1%cell%k+1
    face2%point=face_zmin   
  else 
    write(*,*)'Error in get_other_side_of_faces'
    write(*,*)'No face defined'
    write(*,*)'Face number (face1%point)=',face1%point
    STOP
  end if

  RETURN
  
END SUBROUTINE get_other_side_of_face
