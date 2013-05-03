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
MODULE model_builder

  TYPE real8_coordinates
    real*8 :: x
    real*8 :: y
    real*8 :: z
  END TYPE real8_coordinates

  TYPE integer_coordinates
    integer :: ix
    integer :: iy
    integer :: iz
  END TYPE integer_coordinates
  
  integer tot_n_edges
  
  TYPE(integer_coordinates),allocatable :: edge_list(:,:)
  
  integer MB_npx,MB_npy,MB_npz 			! number of points
  
  real*8 MB_xmin,MB_xmax,MB_ymin,MB_ymax,MB_zmin,MB_zmax  ! limits of mesh in real space
   
  real*8 MB_dl				! cell edge length
  
  real*8 MB_lx
  real*8 MB_ly
  real*8 MB_lz
 
END MODULE model_builder
!
! ________________________________________________________________
!
!
  SUBROUTINE convert_real8_integer_vertex( real_coordinate, integer_coordinate )
  
USE model_builder

IMPLICIT NONE
  
  TYPE(real8_coordinates) :: real_coordinate
  TYPE(integer_coordinates) :: integer_coordinate

! local variables

  real*8 ox,oy,oz
  integer ix,iy,iz
  
  logical error

! START

  error=.FALSE.
  
!  write(*,*)'real8 to integer'
!  write(*,*)real_coordinate%x,MB_xmin,MB_dl
!  write(*,*)real_coordinate%y,MB_ymin,MB_dl
!  write(*,*)real_coordinate%z,MB_zmin,MB_dl

  ox=real_coordinate%x-MB_xmin   	! offset from MB_xmin end of mesh
  ix=NINT(ox/MB_dl)              	! number of cells offset from MB_xmin end of mesh
  integer_coordinate%ix=ix*2	! point number
  if ( (ix.lt.0).OR.(ix.gt.MB_npx) ) then
    write(*,*)'Error in convert_real8_integer_vertex'
    write(*,*)'x point out of range'
    error=.TRUE.
  end if  

  oy=real_coordinate%y-MB_ymin   	! offset from MB_ymin end of mesh
  iy=NINT(oy/MB_dl)              	! number of cells offset from MB_ymin end of mesh
  integer_coordinate%iy=iy*2	! point number
  if ( (iy.lt.0).OR.(iy.gt.MB_npy) ) then
    write(*,*)'Error in convert_real8_integer_vertex'
    write(*,*)'y point out of range'
    error=.TRUE.
  end if  

  oz=real_coordinate%z-MB_zmin   	! offset from MB_zmin end of mesh
  iz=NINT(oz/MB_dl)              	! number of cells offset from MB_zmin end of mesh
  integer_coordinate%iz=iz*2	! point number
  if ( (iz.lt.0).OR.(iz.gt.MB_npz) ) then
    write(*,*)'Error in convert_real8_integer_vertex'
    write(*,*)'z point out of range'
    error=.TRUE.
  end if  
  
  if (.NOT.error) then
    RETURN
  else
    write(*,*)'real x coordinate',real_coordinate%x
    write(*,*)'real y coordinate',real_coordinate%y
    write(*,*)'real z coordinate',real_coordinate%z
    write(*,*)      
    write(*,*)'MB_xmin=',MB_xmin,' MB_xmax=',MB_xmax
    write(*,*)'MB_ymin=',MB_ymin,' MB_ymax=',MB_ymax
    write(*,*)'MB_zmin=',MB_zmin,' MB_zmax=',MB_zmax
    write(*,*)
    write(*,*)'integer x coordinate',integer_coordinate%ix
    write(*,*)'integer y coordinate',integer_coordinate%iy
    write(*,*)'integer z coordinate',integer_coordinate%iz
    write(*,*)      
    write(*,*)'iX npmin=',0,' iX npmax=',MB_npx
    write(*,*)'iY npmin=',0,' iY npmax=',MB_npy
    write(*,*)'iZ npmin=',0,' iZ npmax=',MB_npz
    STOP
  end if
    
  END SUBROUTINE convert_real8_integer_vertex
!
! ________________________________________________________________
!
!

  SUBROUTINE set_edges(ipoint1,ipoint2,ipoint3)
  
USE model_builder

IMPLICIT NONE
  
  TYPE(integer_coordinates) :: ipoint1,ipoint2,ipoint3

! local variables

  integer loop,edge_count
  
  logical set_edge_list_flag

! START

  do loop=1,2 ! first pass, count edges and allocate edge_list, second pass, set edge_list data
 
    if (loop.eq.1) then
      set_edge_list_flag=.FALSE.
    else
      set_edge_list_flag=.TRUE.
    end if
    
    edge_count=0
    
! set the three lines 1->2, 2->3, 3->1
    CALL set_line_edges(ipoint1,ipoint2,edge_count,set_edge_list_flag)
    CALL set_line_edges(ipoint2,ipoint3,edge_count,set_edge_list_flag)
    CALL set_line_edges(ipoint3,ipoint1,edge_count,set_edge_list_flag)
      
    tot_n_edges=edge_count
    
    if (loop.eq.1) then
      allocate( edge_list(1:tot_n_edges,1:2) )
    end if
 
  end do ! next loop
 
  RETURN
  
  END SUBROUTINE set_edges
!
! _______________________________________________________
!
!  
  SUBROUTINE set_line_edges(ipoint1_in,ipoint2_in,edge_count,set_edge_list_flag)
  
USE model_builder

IMPLICIT NONE

  integer edge_count
  
  TYPE(integer_coordinates) :: ipoint1_in,ipoint2_in

  logical set_edge_list_flag
  
! local variables

  TYPE(integer_coordinates) :: ipoint1,ipoint2

  integer iox,ioy,ioz
  integer edge,n_edges,first_edge_count,last_edge_count,local_edge_count
  integer ix,iy,iz
  integer dx,dy,dz
  TYPE(integer_coordinates) :: last_edge_point,edge_point
  TYPE(integer_coordinates) :: itrial_point1,itrial_point2,itrial_point3
  real*8		    :: dist1,dist2,dist3
  real*8		    :: lx1,ly1,lz1
  real*8		    :: lx2,ly2,lz2
  real*8 		    :: px1,py1,pz1
  real*8 		    :: px2,py2,pz2
  real*8 		    :: px3,py3,pz3
  integer		    :: min_point
  real*8		    :: min_dist
  
  integer 		    :: n_points_set
  
  logical reverse_edge
  
! function variables

  logical is_vertex, is_cell

! START

!  write(*,*)'CALLED: Set_line_edges'

! decide on an unambiguous direction in which to traverse the edge

  reverse_edge=.FALSE.
  if (ipoint1_in%ix.gt.ipoint2_in%ix) then
    ipoint1=ipoint1_in
    ipoint2=ipoint2_in
  else if (ipoint1_in%ix.lt.ipoint2_in%ix) then
    ipoint1=ipoint2_in
    ipoint2=ipoint1_in
    reverse_edge=.TRUE.
  else
  
    if (ipoint1_in%iy.gt.ipoint2_in%iy) then
      ipoint1=ipoint1_in
      ipoint2=ipoint2_in
    else if (ipoint1_in%iy.lt.ipoint2_in%iy) then
      ipoint1=ipoint2_in
      ipoint2=ipoint1_in
      reverse_edge=.TRUE.
    else
    
      if (ipoint1_in%iz.gt.ipoint2_in%iz) then
        ipoint1=ipoint1_in
        ipoint2=ipoint2_in
      else if (ipoint1_in%iz.lt.ipoint2_in%iz) then
        ipoint1=ipoint2_in
        ipoint2=ipoint1_in
        reverse_edge=.TRUE.
      else
      
!        write(*,*)'No offset between points'
!	write(*,*)ipoint1_in%ix,ipoint2_in%ix
!	write(*,*)ipoint1_in%iy,ipoint2_in%iy
!	write(*,*)ipoint1_in%iz,ipoint2_in%iz
	RETURN
  
      end if ! offset in z
  
    end if ! offset in y
    
  end if ! offset in x

! check points refer to cell vertices

  if ( (.NOT.is_vertex(ipoint1)).OR.(.NOT.is_vertex(ipoint2)) ) then
    write(*,*)'Error in set_line_edges: points are not both vertex points'
    write(*,*)'point1:',ipoint1%ix,ipoint1%iy,ipoint1%iz
    write(*,*)'point2:',ipoint2%ix,ipoint2%iy,ipoint2%iz
    STOP
  end if

! work out the integer offsets for the points
 
  iox=abs(ipoint2%ix-ipoint1%ix)/2
  ioy=abs(ipoint2%iy-ipoint1%iy)/2
  ioz=abs(ipoint2%iz-ipoint1%iz)/2
  
  dx=0
  dy=0
  dz=0
  
  if ((ipoint2%ix-ipoint1%ix).gt.0) then
    dx=2
  else
    dx=-2
  end if
  if ((ipoint2%iy-ipoint1%iy).gt.0) then
    dy=2
  else
    dy=-2
  end if
  if ((ipoint2%iz-ipoint1%iz).gt.0) then
    dz=2
  else
    dz=-2
  end if
  
  n_edges=iox+ioy+ioz
  first_edge_count=edge_count
  last_edge_count=edge_count+n_edges
  local_edge_count=0
  
  n_points_set=0
  
! set points on edges

! set point1

  ix=ipoint1%ix
  iy=ipoint1%iy
  iz=ipoint1%iz
  
  n_points_set=n_points_set+1
!  write(*,*)'Set point',ix,iy,iz
  
  last_edge_point=ipoint1

! Set points defining the edge line
  CALL get_real8_coordinates( ipoint1,lx1,ly1,lz1 )
  
  CALL get_real8_coordinates( ipoint2,lx2,ly2,lz2 )
  
!  write(*,*)'End point1',ipoint1%ix,ipoint1%iy,ipoint1%iz
!  write(*,*)'End point2',ipoint2%ix,ipoint2%iy,ipoint2%iz
!  write(*,*)'n_edges=',n_edges
!  write(*,*)'First point',ix,iy,iz
  
  do edge=1,n_edges

! set three trial points with increments of dx,dy and dz respectively
    itrial_point1=last_edge_point 
    itrial_point1%ix=itrial_point1%ix+dx
    CALL get_real8_coordinates( itrial_point1,px1,py1,pz1 )
    
    itrial_point2=last_edge_point 
    itrial_point2%iy=itrial_point2%iy+dy
    CALL get_real8_coordinates( itrial_point2,px2,py2,pz2 )
    
    itrial_point3=last_edge_point 
    itrial_point3%iz=itrial_point3%iz+dz
    CALL get_real8_coordinates( itrial_point3,px3,py3,pz3 )
  
! calculate the closest distance between the three trial points and the line

    CALL dist_line_point(lx1,ly1,lz1,lx2,ly2,lz2,px1,py1,pz1,dist1)
    CALL dist_line_point(lx1,ly1,lz1,lx2,ly2,lz2,px2,py2,pz2,dist2)
    CALL dist_line_point(lx1,ly1,lz1,lx2,ly2,lz2,px3,py3,pz3,dist3)
    
    min_point=1
    min_dist=dist1
    ix=itrial_point1%ix
    iy=itrial_point1%iy
    iz=itrial_point1%iz
    
    if (dist2.lt.min_dist) then
      min_point=2
      min_dist=dist2
      ix=itrial_point2%ix
      iy=itrial_point2%iy
      iz=itrial_point2%iz
    end if
    if (dist3.lt.min_dist) then
      min_point=3
      min_dist=dist3
      ix=itrial_point3%ix
      iy=itrial_point3%iy
      iz=itrial_point3%iz
    end if
        
    edge_point%ix=ix
    edge_point%iy=iy
    edge_point%iz=iz
!    write(*,*)'Edge point',edge_point%ix,edge_point%iy,edge_point%iz
        
    n_points_set=n_points_set+1
!    write(*,*)'Set point',ix,iy,iz
    
    edge_count=edge_count+1
    local_edge_count=local_edge_count+1
    if (set_edge_list_flag) then ! add this point to the edge list
    
      if (.NOT.reverse_edge) then
        edge_list(first_edge_count+local_edge_count,1)=last_edge_point
        edge_list(first_edge_count+local_edge_count,2)=edge_point
      else
        edge_list(last_edge_count-local_edge_count+1,1)=edge_point
        edge_list(last_edge_count-local_edge_count+1,2)=last_edge_point
      end if
    end if
    
    last_edge_point=edge_point
    
  end do ! next edge
    
!  write(*,*)'End  point',ix,iy,iz
!  write(*,*)'End point2',ipoint2%ix,ipoint2%iy,ipoint2%iz

! check endpoint coincides with point2  
  if ( (ix.ne.ipoint2%ix).OR.(iy.ne.ipoint2%iy).OR.(iz.ne.ipoint2%iz) ) then
    write(*,*)'Error in set_line_edges: end of line is not point2'
    write(*,*)'point1   =',ipoint1%ix,ipoint1%iy,ipoint1%iz
    write(*,*)'point2   =',ipoint2%ix,ipoint2%iy,ipoint2%iz
    write(*,*)'end point=',ix,iy,iz
    STOP
  end if
  
!  write(*,*)'Number of points set=',n_points_set
  
!  write(*,*)'FINISHED: Set_line_edges'
    
  RETURN
  
  END SUBROUTINE set_line_edges
!
! ___________________________________________________________
!
!
  SUBROUTINE get_face_list(local_face_list,max_n_faces,n_faces,x1,y1,z1,x2,y2,z2,x3,y3,z3)

USE model_builder
!USE geometry
!USE mesh

USE TLM_general
!USE geometry_types

IMPLICIT NONE
  
  integer 		:: max_n_faces,n_faces
  integer:: local_face_list(1:max_n_faces,1:4)
  
  real*8		:: x1,y1,z1
  real*8		:: x2,y2,z2
  real*8		:: x3,y3,z3

! local variables
  
  TYPE(real8_coordinates) :: rpoint1,rpoint2,rpoint3
  TYPE(integer_coordinates) :: ipoint1,ipoint2,ipoint3

! START

!  write(*,*)'Called: get_face_list'

! set model builder global mesh parameters
  MB_xmin=mesh_xmin
  MB_xmax=mesh_xmax
  MB_ymin=mesh_ymin
  MB_ymax=mesh_ymax
  MB_zmin=mesh_zmin
  MB_zmax=mesh_zmax
  MB_npx=nx*2
  MB_npy=ny*2
  MB_npz=nz*2
  MB_lx=nx*dl
  MB_ly=ny*dl
  MB_lz=nz*dl
  MB_dl=dl

! Put triangle points onto the integer grid

  rpoint1%x=x1
  rpoint1%y=y1
  rpoint1%z=z1
  rpoint2%x=x2
  rpoint2%y=y2
  rpoint2%z=z2
  rpoint3%x=x3
  rpoint3%y=y3
  rpoint3%z=z3

  CALL convert_real8_integer_vertex( rpoint1,ipoint1 )
  CALL convert_real8_integer_vertex( rpoint2,ipoint2 )
  CALL convert_real8_integer_vertex( rpoint3,ipoint3 )
    
! set surface triangle edges 

  CALL set_edges(ipoint1,ipoint2,ipoint3)

! set surface triangle faces 

  CALL set_faces(local_face_list,max_n_faces,n_faces,ipoint1,ipoint2,ipoint3,rpoint1,rpoint2,rpoint3)
  
  if ( allocated (edge_list) ) deallocate( edge_list )

!  write(*,*)'Finished: get_face_list'

END SUBROUTINE get_face_list
!
! _______________________________________________________
!
!  
  SUBROUTINE set_faces( local_face_list,max_n_faces,n_faces,tipoint1,tipoint2,tipoint3	&
                                                           ,trpoint1,trpoint2,trpoint3 )
USE model_builder
USE cell_parameters

IMPLICIT NONE
  
  integer 		:: max_n_faces,n_faces
  integer:: local_face_list(1:max_n_faces,1:4)
  TYPE(integer_coordinates) :: tipoint1,tipoint2,tipoint3
  TYPE(real8_coordinates) :: trpoint1,trpoint2,trpoint3

! local variables

  integer triangle,edge_count,edge
  integer,allocatable :: mesh_point(:,:,:)
  TYPE(integer_coordinates) :: ipoint1,ipoint2,ipoint3,ipoint4,ipoint5

  integer n_edges,new_n_edges,edge0,edge1,edge2,edge3,edge4
  TYPE(integer_coordinates),allocatable :: local_edge_list(:,:)
  integer,allocatable			:: edges_to_be_removed(:)
  TYPE(integer_coordinates),allocatable :: new_local_edge_list(:,:)
  integer first_edge,last_edge
  
  integer n_4point_faces
  integer n_3point_faces
  integer n_2point_faces
  integer n_doubled_edges
  integer n_edges_to_be_removed
  
  logical result
  
  integer fill_iteration
  
  integer new_corner_face
  
  integer edge_unit
  
  logical inside1
  logical is_corner1
  integer ix1,iy1,iz1
  real*8  dist1
  TYPE(integer_coordinates) :: new_ipoint1
  
  logical inside2
  logical is_corner2
  integer ix2,iy2,iz2
  real*8  dist2
  TYPE(integer_coordinates) :: new_ipoint2
  
  logical moved_point
  integer ixc,iyc,izc
  integer ox,oy,oz
  
  integer ix,iy,iz
  integer ixmin,ixmax,iymin,iymax,izmin,izmax
  integer tixmin,tixmax,tiymin,tiymax,tizmin,tizmax
  
  integer n_x_faces,n_y_faces,n_z_faces,max_faces,face
  
  integer tot_faces,n_faces_set,tot_n_faces_set
  integer n_points_set
  
  TYPE(integer_coordinates),allocatable :: x_face_list(:)
  TYPE(integer_coordinates),allocatable :: y_face_list(:)
  TYPE(integer_coordinates),allocatable :: z_face_list(:)
  
  integer n_new_faces
  
  integer n_check_points
  
  integer,save	:: xface_count=0
  integer,save	:: yface_count=0
  integer,save	:: zface_count=0
    
! function variables 
  logical :: same_point
  logical :: same_edge


! START

!  write(*,*)'Called: set_faces'
      
      tixmin=min(tipoint1%ix,tipoint2%ix,tipoint3%ix)
      tiymin=min(tipoint1%iy,tipoint2%iy,tipoint3%iy)
      tizmin=min(tipoint1%iz,tipoint2%iz,tipoint3%iz)
      
      tixmax=max(tipoint1%ix,tipoint2%ix,tipoint3%ix)
      tiymax=max(tipoint1%iy,tipoint2%iy,tipoint3%iy)
      tizmax=max(tipoint1%iz,tipoint2%iz,tipoint3%iz)
            
      allocate( mesh_point(tixmin-2:tixmax+2,tiymin-2:tiymax+2,tizmin-2:tizmax+2) )
      mesh_point(tixmin-2:tixmax+2,tiymin-2:tiymax+2,tizmin-2:tizmax+2)=0

      n_edges=tot_n_edges
      
      allocate( local_edge_list(1:n_edges,1:2) )      
      local_edge_list(1:n_edges,1:2)=edge_list(1:n_edges,1:2)

! Set point values to be 1 
      do edge=1,n_edges
      
        ipoint1=local_edge_list(edge,1)
        ipoint2=local_edge_list(edge,2)
	
	mesh_point(ipoint1%ix,ipoint1%iy,ipoint1%iz)=1
	mesh_point(ipoint2%ix,ipoint2%iy,ipoint2%iz)=1
	
      end do

! Set faces initially without regard for quality of surface  
  
      max_faces=(tixmax-tixmin)*(tiymax-tiymin)+(tiymax-tiymin)*(tizmax-tizmin)+(tizmax-tizmin)*(tixmax-tixmin)
      n_x_faces=0
      n_y_faces=0
      n_z_faces=0
      
      allocate( x_face_list(1:max_faces) )
      allocate( y_face_list(1:max_faces) )
      allocate( z_face_list(1:max_faces) )
      
! Faces normal to the x direction

      do iy=tiymin,tiymax-2,2

        izmin=tizmax
        izmax=0

        do edge=1,n_edges
      
          ipoint1=local_edge_list(edge,1)
          ipoint2=local_edge_list(edge,2)
	  
	  if (ipoint1%iy.ne.ipoint2%iy) then ! y directed edge
	
! at this y value calculate the minimum and maximum z values
     
            iymin=min(ipoint1%iy,ipoint2%iy)
	    
	    if (iymin.eq.iy) then ! calculate z range
	    
              izmin=min(ipoint1%iz,izmin)     
              izmin=min(ipoint2%iz,izmin)     
              izmax=max(ipoint1%iz,izmax)
              izmax=max(ipoint2%iz,izmax)
	
	    end if
	    
	  end if ! y directed edge
	
        end do ! next edge
	
! we now have two points which define a line of faces in the z direction

        do iz=izmin,izmax-2,2 ! loop over z working out the x value for the face
	
	  iy2=iy+1
	  iz2=iz+1
	  
          n_x_faces=n_x_faces+1
          x_face_list(n_x_faces)%ix=-1
          x_face_list(n_x_faces)%iy=iy2
          x_face_list(n_x_faces)%iz=iz2
	      
	end do ! next z point
	
      end do ! next x coordinate
      
! Faces normal to the y direction

      do iz=tizmin,tizmax-2,2

        ixmin=tixmax
        ixmax=0

        do edge=1,n_edges
      
          ipoint1=local_edge_list(edge,1)
          ipoint2=local_edge_list(edge,2)
	  
	  if (ipoint1%iz.ne.ipoint2%iz) then ! z directed edge
	
! at this z value calculate the minimum and maximum x values
     
            izmin=min(ipoint1%iz,ipoint2%iz)
	    
	    if (izmin.eq.iz) then ! calculate x range
	    
              ixmin=min(ipoint1%ix,ixmin)     
              ixmin=min(ipoint2%ix,ixmin)     
              ixmax=max(ipoint1%ix,ixmax)
              ixmax=max(ipoint2%ix,ixmax)
	
	    end if
	    
	  end if ! z directed edge
	
        end do ! next edge
	
! we now have two points which define a line of faces in the x direction

        do ix=ixmin,ixmax-2,2 ! loop over x working out the y value for the face
	
	  ix2=ix+1
	  iz2=iz+1
	  
          n_y_faces=n_y_faces+1
          y_face_list(n_y_faces)%ix=ix2
          y_face_list(n_y_faces)%iy=-1
          y_face_list(n_y_faces)%iz=iz2
	  
	end do ! next x point
	
      end do ! next z coordinate

! Faces normal to the z direction

      do ix=tixmin,tixmax-2,2

        iymin=tiymax
        iymax=0

        do edge=1,n_edges
      
          ipoint1=local_edge_list(edge,1)
          ipoint2=local_edge_list(edge,2)
	  
	  if (ipoint1%ix.ne.ipoint2%ix) then ! x directed edge
	
! at this x value calculate the minimum and maximum y values
     
            ixmin=min(ipoint1%ix,ipoint2%ix)
	    
	    if (ixmin.eq.ix) then ! calculate y range
	    	    
              iymin=min(ipoint1%iy,iymin)     
              iymin=min(ipoint2%iy,iymin)     
              iymax=max(ipoint1%iy,iymax)
              iymax=max(ipoint2%iy,iymax)
	
	    end if
	    
	  end if ! x directed edge
	
        end do ! next edge
	
! we now have two points which define a line of faces in the y direction

        do iy=iymin,iymax-2,2 ! loop over y working out the z value for the face
	
	  ix2=ix+1
	  iy2=iy+1
	  
          n_z_faces=n_z_faces+1
          z_face_list(n_z_faces)%ix=ix2
          z_face_list(n_z_faces)%iy=iy2
          z_face_list(n_z_faces)%iz=-1
	  
	end do ! next y point
	
      end do ! next x coordinate
      
      if (n_x_faces+n_y_faces+n_z_faces.eq.0) then

        if (allocated( local_edge_list )) deallocate( local_edge_list )
      
        if (allocated( x_face_list )) deallocate( x_face_list )
        if (allocated( y_face_list )) deallocate( y_face_list )
        if (allocated( z_face_list )) deallocate( z_face_list )
	
	n_faces=0
      
        return
      
      end if
      
      xface_count=xface_count+n_x_faces
      yface_count=yface_count+n_y_faces
      zface_count=zface_count+n_z_faces
            
! Set edge values to be 1 (effectively accounts for the edge on an adjacent triangle)
      do edge=1,n_edges
      
        ipoint1=local_edge_list(edge,1)
        ipoint2=local_edge_list(edge,2)
	
	ix=(ipoint1%ix+ipoint2%ix)/2
	iy=(ipoint1%iy+ipoint2%iy)/2
	iz=(ipoint1%iz+ipoint2%iz)/2
	
	mesh_point(ix,iy,iz)=1
	
      end do
      
      tot_faces=n_x_faces+n_y_faces+n_z_faces
      
      tot_n_faces_set=0

10    CONTINUE   

       n_faces_set=0
	
       n_check_points=4   ! set 4 point constraints before 3 point constraints

11      CONTINUE
	
!       set x directed faces

        do face=1,n_x_faces
	
	  ix=x_face_list(face)%ix
	  iy=x_face_list(face)%iy
	  iz=x_face_list(face)%iz
	  
	  if (ix.eq.-1) then
	  
	    do ix=tixmin,tixmax,2
	  
	      n_points_set=mesh_point(ix,iy-1,iz-1)+	&
	                   mesh_point(ix,iy-1,iz+1)+	&
	                   mesh_point(ix,iy+1,iz-1)+	&
	                   mesh_point(ix,iy+1,iz+1)
			   
	      if (n_points_set.GE.n_check_points) then ! check edge conditions 
	      
	        if ( (mesh_point(ix,iy-1,iz).NE.2) .AND. 	&
		     (mesh_point(ix,iy+1,iz).NE.2) .AND. 	&
		     (mesh_point(ix,iy,iz-1).NE.2) .AND. 	&
		     (mesh_point(ix,iy,iz+1).NE.2) ) then  ! we have not overloaded the edges
		     
	    
	          mesh_point(ix,iy-1,iz-1)=1  	! set points
	          mesh_point(ix,iy-1,iz+1)=1
	          mesh_point(ix,iy+1,iz-1)=1
	          mesh_point(ix,iy+1,iz+1)=1
		  
	          mesh_point(ix,iy,iz)=1	! set face

                  mesh_point(ix,iy-1,iz)=mesh_point(ix,iy-1,iz)+1	! set edges
                  mesh_point(ix,iy+1,iz)=mesh_point(ix,iy+1,iz)+1
                  mesh_point(ix,iy,iz-1)=mesh_point(ix,iy,iz-1)+1
                  mesh_point(ix,iy,iz+1)=mesh_point(ix,iy,iz+1)+1
		  
		  x_face_list(face)%ix=ix
	          n_faces_set=n_faces_set+1	
		  	  	      	  
!	          GOTO 1
	          GOTO 12
		  
	        end if ! not overlaoded edges 
	      
	      end if ! n_points_set.GE.n_check_points
	    
	    end do ! next x value
	  
1           CONTINUE
	  
	  end if
	  
	end do ! next face normal to x

	
!       set y directed faces

        do face=1,n_y_faces
	
	  ix=y_face_list(face)%ix
	  iy=y_face_list(face)%iy
	  iz=y_face_list(face)%iz
	  
	  if (iy.eq.-1) then
	  
	    do iy=tiymin,tiymax,2
	  
	      n_points_set=mesh_point(ix-1,iy,iz-1)+	&
	                   mesh_point(ix-1,iy,iz+1)+	&
	                   mesh_point(ix+1,iy,iz-1)+	&
	                   mesh_point(ix+1,iy,iz+1)
			   
	      if (n_points_set.GE.n_check_points) then  ! check edge conditions 
	      	      	      
	        if ( (mesh_point(ix-1,iy,iz).NE.2) .AND. 	&
		     (mesh_point(ix+1,iy,iz).NE.2) .AND. 	&
		     (mesh_point(ix,iy,iz-1).NE.2) .AND. 	&
		     (mesh_point(ix,iy,iz+1).NE.2) ) then  ! we have not overloaded the edges
		     
		     
	          mesh_point(ix-1,iy,iz-1)=1	! set points
	          mesh_point(ix-1,iy,iz+1)=1
	          mesh_point(ix+1,iy,iz-1)=1
	          mesh_point(ix+1,iy,iz+1)=1
		  
	          mesh_point(ix,iy,iz)=1      ! set face

                  mesh_point(ix-1,iy,iz)=mesh_point(ix-1,iy,iz)+1	! set edges
                  mesh_point(ix+1,iy,iz)=mesh_point(ix+1,iy,iz)+1
                  mesh_point(ix,iy,iz-1)=mesh_point(ix,iy,iz-1)+1
                  mesh_point(ix,iy,iz+1)=mesh_point(ix,iy,iz+1)+1
		  
		  y_face_list(face)%iy=iy
	          n_faces_set=n_faces_set+1
		  	  		  
!	          GOTO 2
	          GOTO 12
		  
	        end if
		
	      end if
	    
	    end do ! next y value
	  
2           CONTINUE

          end if
	  
	end do ! next face normal to y

	
!       set z directed faces

        do face=1,n_z_faces
	
	  ix=z_face_list(face)%ix
	  iy=z_face_list(face)%iy
	  iz=z_face_list(face)%iz
	  
	  if (iz.eq.-1) then
	  
	    do iz=tizmin,tizmax,2
	  
	      n_points_set=mesh_point(ix-1,iy-1,iz)+	&
	                   mesh_point(ix+1,iy-1,iz)+	&
	                   mesh_point(ix-1,iy+1,iz)+	&
	                   mesh_point(ix+1,iy+1,iz)
			 
	      if (n_points_set.GE.n_check_points) then ! check edge conditions 
	      
	        if ( (mesh_point(ix,iy-1,iz).NE.2) .AND. 	&
		     (mesh_point(ix,iy+1,iz).NE.2) .AND. 	&
		     (mesh_point(ix-1,iy,iz).NE.2) .AND. 	&
		     (mesh_point(ix+1,iy,iz).NE.2) ) then  ! we have not overloaded the edges
	    
	          mesh_point(ix-1,iy-1,iz)=1	! set points
	          mesh_point(ix+1,iy-1,iz)=1
	          mesh_point(ix-1,iy+1,iz)=1
	          mesh_point(ix+1,iy+1,iz)=1
		  
	          mesh_point(ix,iy,iz)=1	! set face

                  mesh_point(ix-1,iy,iz)=mesh_point(ix-1,iy,iz)+1	! set edges
                  mesh_point(ix+1,iy,iz)=mesh_point(ix+1,iy,iz)+1
                  mesh_point(ix,iy-1,iz)=mesh_point(ix,iy-1,iz)+1
                  mesh_point(ix,iy+1,iz)=mesh_point(ix,iy+1,iz)+1
		
		  z_face_list(face)%iz=iz
	          n_faces_set=n_faces_set+1
		  	  
!	          GOTO 3
	          GOTO 12
		  
	        end if ! 
		
	      end if !
	    
	    end do ! next z value
	  
3           CONTINUE

          end if
	  
	end do ! next face normal to z
	
12       CONTINUE

        tot_n_faces_set=tot_n_faces_set+n_faces_set
	
        if ( (tot_n_faces_set.NE.tot_faces).AND.(n_faces_set.NE.0) ) then
	  n_check_points=4
	  GOTO 10
	end if
	
        if ( (tot_n_faces_set.NE.tot_faces).AND.(n_check_points.eq.4) ) then
	  n_check_points=3
	  GOTO 11
	end if 

      if (tot_n_faces_set.NE.tot_faces) then
      
        write(*,*)'Error in set_faces'
	write(*,*)'Failed to set all of the faces'	
	write(*,*)
	write(*,*)'tot_n_faces_set=',tot_n_faces_set,' n_faces_set=',n_faces_set
	write(*,*)
	write(*,*)'Edge points'
	
        do edge=1,n_edges
      
          ipoint1=local_edge_list(edge,1)
          ipoint2=local_edge_list(edge,2)
	  write(*,8500)edge,ipoint1%ix,ipoint1%iy,ipoint1%iz,ipoint2%ix,ipoint2%iy,ipoint2%iz
8500      format(I5,6I4)		
        end do
	
	write(*,*)'X face points'
	do face=1,n_x_faces	  
          write(*,*)face,x_face_list(face)%ix,x_face_list(face)%iy,x_face_list(face)%iz	  
        end do
	
	write(*,*)'Y face points'
	do face=1,n_y_faces	  
          write(*,*)face,y_face_list(face)%ix,y_face_list(face)%iy,y_face_list(face)%iz	  
        end do
	
	write(*,*)'Z face points'
	do face=1,n_z_faces	  
          write(*,*)face,z_face_list(face)%ix,z_face_list(face)%iy,z_face_list(face)%iz	  
        end do
	
	write(*,*)'Triangle limits'
	write(*,*)'tixmin=',tixmin,' tixmax',tixmax
	write(*,*)'tiymin=',tiymin,' tiymax',tiymax
	write(*,*)'tizmin=',tizmin,' tizmax',tizmax

	STOP
	
      end if

      if (allocated( local_edge_list )) deallocate( local_edge_list )
      
      if (allocated( x_face_list )) deallocate( x_face_list )
      if (allocated( y_face_list )) deallocate( y_face_list )
      if (allocated( z_face_list )) deallocate( z_face_list )
      
! Iterative improvement of surface      

!      write(*,*)'Iterative improvement of mesh'
      
!      write(*,*)'tixmin=',tixmin,' tixmax=',tixmax
!      write(*,*)'tiymin=',tiymin,' tiymax=',tiymax
!      write(*,*)'tizmin=',tizmin,' tizmax=',tizmax
          	  
20    CONTINUE

        moved_point=.FALSE.	  
	  
	do ix=tixmin,tixmax-2,2  
	  do iy=tiymin,tiymax-2,2  
	    do iz=tizmin,tizmax-2,2  

! get cell centre coordinate	    
	      ixc=ix+1
	      iyc=iy+1
	      izc=iz+1

! check each of the 8 corners 
              do ox=-1,1,2
                do oy=-1,1,2
                  do oz=-1,1,2
		  
                    if ( (mesh_point(ixc+ox,iyc   ,izc   ).EQ.1).AND.	&
		         (mesh_point(ixc   ,iyc+oy,izc   ).EQ.1).AND.	&
		         (mesh_point(ixc   ,iyc   ,izc+oz).EQ.1)	) then
			 
! this is a corner point, check distance from triangle to this point and the opposite point,
! and swap if the second point is better

                      ipoint1%ix=ixc+ox
                      ipoint1%iy=iyc+oy
                      ipoint1%iz=izc+oz
		      CALL dist_ipoint_surface(trpoint1,trpoint2,trpoint3,ipoint1,dist1)
		      
                      ipoint2%ix=ixc-ox
                      ipoint2%iy=iyc-oy
                      ipoint2%iz=izc-oz
		      CALL dist_ipoint_surface(trpoint1,trpoint2,trpoint3,ipoint2,dist2)
		      
		      if (dist2.lt.dist1) then ! swap faces across cube
		      
		        mesh_point(ixc+ox,iyc   ,izc   )=0
		        mesh_point(ixc   ,iyc+oy,izc   )=0
		        mesh_point(ixc   ,iyc   ,izc+oz)=0
		        mesh_point(ixc-ox,iyc   ,izc   )=1
		        mesh_point(ixc   ,iyc-oy,izc   )=1
		        mesh_point(ixc   ,iyc   ,izc-oz)=1
			
                        moved_point=.TRUE.	  
		       
		      end if

	            end if			 
			 
		  end do
		end do
              end do
    
            end do ! next z
	  end do ! next y
	end do ! next x
      
      if (moved_point) goto 20 ! go for another iteration
      
!      write(*,*)'Calculate face list'

! calcualte the face list

	n_faces=0	
	  
	do ix=tixmin,tixmax,2  
	  do iy=tiymin,tiymax,2  
	    do iz=tizmin,tizmax,2  

! get cell centre coordinate	    
	      ixc=ix+1
	      iyc=iy+1
	      izc=iz+1

! check each of the 3 min faces
	      if (mesh_point(ixc-1,iyc,izc).EQ.1) then
	     	n_faces=n_faces+1
	        if (n_faces.gt.max_n_faces) goto 9000	    
	        local_face_list(n_faces,1)=(ixc+1)/2
	        local_face_list(n_faces,2)=(iyc+1)/2
	        local_face_list(n_faces,3)=(izc+1)/2
	        local_face_list(n_faces,4)=face_xmin
              end if
	      
	      if (mesh_point(ixc,iyc-1,izc).EQ.1) then
	     	n_faces=n_faces+1
	        if (n_faces.gt.max_n_faces) goto 9000	    
	        local_face_list(n_faces,1)=(ixc+1)/2
	        local_face_list(n_faces,2)=(iyc+1)/2
	        local_face_list(n_faces,3)=(izc+1)/2
	        local_face_list(n_faces,4)=face_ymin
              end if
	      
	      if (mesh_point(ixc,iyc,izc-1).EQ.1) then
	     	n_faces=n_faces+1
	        if (n_faces.gt.max_n_faces) goto 9000	    
	        local_face_list(n_faces,1)=(ixc+1)/2
	        local_face_list(n_faces,2)=(iyc+1)/2
	        local_face_list(n_faces,3)=(izc+1)/2
	        local_face_list(n_faces,4)=face_zmin
              end if
		    
            end do ! next z
	  end do ! next y
	end do ! next x
	
	if (n_faces.ne.tot_faces) then
	  write(*,*)'Discrepancy in number of faces: set_faces'
	  write(*,*)'Number of faces=',n_faces,tot_faces
	  STOP
	end if
	
      if ( allocated(mesh_point) ) deallocate( mesh_point )
      
!    write(*,*)'FINISHED: set_faces'
  
  RETURN

9000 write(*,*)'ERROR in set_faces'
     write(*,*)'Exceeded maximum number of faces in triangle'
     STOP 
  
  END SUBROUTINE set_faces
!
! ___________________________________________
!
!
  logical FUNCTION is_vertex(ipoint)

USE model_builder
  
  IMPLICIT NONE
  
  TYPE(integer_coordinates) :: ipoint
  logical result
  
  integer ix,iy,iz
  
! START

  result=.TRUE.

  ix=ipoint%ix
  iy=ipoint%iy
  iz=ipoint%iz
  
  if ( (mod(ix,2).ne.0).OR.(mod(iy,2).ne.0).OR.(mod(iz,2).ne.0) ) then
    result=.FALSE.
  end if
  
  is_vertex=result
  
  RETURN
  
  END FUNCTION is_vertex
!
! ________________________________________________________________
!
!
  SUBROUTINE get_real8_coordinates( integer_coordinate,x,y,z )
  
USE model_builder

IMPLICIT NONE
  
  TYPE(integer_coordinates) :: integer_coordinate
  real*8 x,y,z

! local variables
  
  logical error

! START

  error=.FALSE.

  if ( (integer_coordinate%ix.lt.0).OR.(integer_coordinate%ix.gt.MB_npx) ) then
    write(*,*)'Error in convert_real8_integer_cell'
    write(*,*)'x point out of range'
    error=.TRUE.
  else 
    x=MB_xmin+MB_lx*dble(integer_coordinate%ix)/dble(MB_npx)
  end if  

  if ( (integer_coordinate%iy.lt.0).OR.(integer_coordinate%iy.gt.MB_npy) ) then
    write(*,*)'Error in convert_real8_integer_cell'
    write(*,*)'y point out of range'
    error=.TRUE.
  else 
    y=MB_ymin+MB_ly*dble(integer_coordinate%iy)/dble(MB_npy)
  end if  

  if ( (integer_coordinate%iz.lt.0).OR.(integer_coordinate%iz.gt.MB_npz) ) then
    write(*,*)'Error in convert_real8_integer_cell'
    write(*,*)'z point out of range'
    error=.TRUE.
  else 
    z=MB_zmin+MB_lz*dble(integer_coordinate%iz)/dble(MB_npz)
  end if  
  
  if (.NOT.error) then
    RETURN
  else
    write(*,*)'integer x coordinate',integer_coordinate%ix
    write(*,*)'integer y coordinate',integer_coordinate%iy
    write(*,*)'integer z coordinate',integer_coordinate%iz
    write(*,*)      
    write(*,*)'iX npmin=',0,' iX npmax=',MB_npx
    write(*,*)'iY npmin=',0,' iY npmax=',MB_npy
    write(*,*)'iZ npmin=',0,' iZ npmax=',MB_npz
    write(*,*)
    STOP
  end if
    
  END SUBROUTINE get_real8_coordinates
!
! ___________________________________________
!
!
  SUBROUTINE dist_line_point(lx1,ly1,lz1,lx2,ly2,lz2,px,py,pz,dist)

! calculate the shortest distance between a point and a line
  
  real*8  :: lx1,ly1,lz1,lx2,ly2,lz2,px,py,pz,dist
  
! local variables  

  real*8 u,px2,py2,pz2,length_sqr
  
! START

  if ( (lx1.eq.lx2).AND.(ly1.eq.ly2).AND.(lz1.eq.lz2) ) then
    write(*,*)'Error: distance between line defining points is zero'
    write(*,*)'Line point 1:',lx1,ly1,lz1
    write(*,*)'Line point 2:',lx2,ly2,lz2
    write(*,*)'Point       :',px,py,pz
    STOP
  end if
  
  length_sqr=(lx2-lx1)**2+(ly2-ly1)**2+(lz2-lz1)**2
  
  u=( (px-lx1)*(lx2-lx1)+(py-ly1)*(ly2-ly1)+(pz-lz1)*(lz2-lz1) )/length_sqr
  
  px2=lx1+u*(lx2-lx1)
  py2=ly1+u*(ly2-ly1)
  pz2=lz1+u*(lz2-lz1)
  
  dist=sqrt( (px-px2)**2+(py-py2)**2+(pz-pz2)**2 )

  RETURN
  
  END
!
! _____________________________________________________________________
!
!
  logical function same_point( ipoint1,ipoint2 )
  
USE model_builder

IMPLICIT NONE
  
  TYPE(integer_coordinates) :: ipoint1,ipoint2
 
! START

  if ( (ipoint1%ix.eq.ipoint2%ix).AND.(ipoint1%iy.eq.ipoint2%iy).AND.(ipoint1%iz.eq.ipoint2%iz) ) then
    same_point=.TRUE.
  else
    same_point=.FALSE.
  end if
  
  RETURN
  
  END function same_point
  
!
! _____________________________________________________________________
!
!
  logical function same_edge( ipoint1,ipoint2,ipoint3,ipoint4 )
  
USE model_builder

IMPLICIT NONE
  
  TYPE(integer_coordinates) :: ipoint1,ipoint2,ipoint3,ipoint4

logical same_point
 
! START

  same_edge=.FALSE.
  if ( (same_point(ipoint1,ipoint3)).AND.(same_point(ipoint2,ipoint4)) ) then
    same_edge=.TRUE.
  else if ( (same_point(ipoint1,ipoint4)).AND.(same_point(ipoint2,ipoint3)) ) then
    same_edge=.TRUE.
  end if
  
  RETURN
  
  END function same_edge
!
! ______________________________________________________________________________
!
!  
  SUBROUTINE dist_ipoint_surface(rpoint1,rpoint2,rpoint3,ipoint,dist)
  
!USE solver_general
USE file_information
USE model_builder

IMPLICIT NONE

! variables passed to subroutine
   
  TYPE(real8_coordinates) :: rpoint1,rpoint2,rpoint3  
  TYPE(integer_coordinates) :: ipoint
  real*8 dist
  
! local variables
   
  TYPE(real8_coordinates) :: v1,v2,vp,normal
  real*8 mod_normal
  real*8 p
  real*8 x,y,z
  real*8 offset
  real*8 xcell,ycell,zcell
  
! START  

  CALL get_real8_coordinates( ipoint,x,y,z )

  v1%x=rpoint2%x-rpoint1%x
  v1%y=rpoint2%y-rpoint1%y
  v1%z=rpoint2%z-rpoint1%z
  
  v2%x=rpoint3%x-rpoint1%x
  v2%y=rpoint3%y-rpoint1%y
  v2%z=rpoint3%z-rpoint1%z
  
  normal%x=v1%y*v2%z-v1%z*v2%y
  normal%y=v1%z*v2%x-v1%x*v2%z
  normal%z=v1%x*v2%y-v1%y*v2%x
  
  mod_normal=sqrt(normal%x**2+normal%y**2+normal%z**2)
  
  if (mod_normal.eq.0d0) then
! error in triangle
    write(*,*)'Error in trinagle, no normal can be found'
    write(warning_file_unit,*)'point1:',rpoint1%x,rpoint1%y,rpoint1%z
    write(warning_file_unit,*)'point2:',rpoint2%x,rpoint2%y,rpoint2%z
    write(warning_file_unit,*)'point3:',rpoint3%x,rpoint3%y,rpoint3%z
    write(warning_file_unit,*)'v1:',v1%x,v1%y,v1%z
    write(warning_file_unit,*)'v2:',v2%x,v2%y,v2%z
    write(warning_file_unit,*)'normal',normal%x,normal%y,normal%z
    STOP
  end if
  
  normal%x=normal%x/mod_normal
  normal%y=normal%y/mod_normal
  normal%z=normal%z/mod_normal
  
  vp%x=x-rpoint1%x
  vp%y=y-rpoint1%y
  vp%z=z-rpoint1%z
  
  dist=abs(vp%x*normal%x+vp%y*normal%y+vp%z*normal%z)
    
  RETURN
    
  END SUBROUTINE dist_ipoint_surface
  

