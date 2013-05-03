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
! SUBROUTINE Output_cable_geometry()
! SUBROUTINE Output_cable_route()
! SUBROUTINE Output_bundle_geometry()
! SUBROUTINE Output_junction_specification()
! SUBROUTINE Output_face_junction_specification()


!
!SUBROUTINE Output_cable_geometry
!
! NAME
!     SUBROUTINE Output_cable_geometry
!
! DESCRIPTION
!      
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/11/2012 CJS
!
!
SUBROUTINE Output_cable_geometry()

USE TLM_general
USE Cables
USE File_information

IMPLICIT NONE

! local variables

  integer cable_geometry
  integer first_cable_geometry
  integer last_cable_geometry
  integer n_rows,n_cols
  integer row,col

! START

  CALL write_line('CALLED: Output_cable_geometry',0,output_to_screen_flag)

  write(*,*)'Number of cable geometries=',n_cable_geometries
  
10 continue

  write(*,*)
  write(*,*)'Enter the number of the cable geometry to output or 0 to output all of then'
  read(*,*)cable_geometry
  
  if (cable_geometry.eq.0) then
    first_cable_geometry=1
    last_cable_geometry=n_cable_geometries
  else 
  
    if ( (cable_geometry.lt.0).OR.(cable_geometry.gt.n_cable_geometries) ) then
      write(*,*)'Cable geometry is outside the available range'
      RETURN
    end if
    
    first_cable_geometry=cable_geometry
    last_cable_geometry=cable_geometry
  end if
  
  do cable_geometry=first_cable_geometry,last_cable_geometry
  
    write(*,*)
    write(*,*)'Cable geometry number		=',cable_geometry
    write(*,*)
    write(*,*)'Geometry type			=',trim(cable_geometry_list(cable_geometry)%cable_geometry_type_string)
    write(*,*)'Geometry type number		=',cable_geometry_list(cable_geometry)%cable_geometry_type
    write(*,*)'Number of conductors		=',cable_geometry_list(cable_geometry)%n_conductors
    write(*,*)'Number of Shielded conductors	=',cable_geometry_list(cable_geometry)%n_shielded_conductors
    write(*,*)'External_conductor_radius	=',cable_geometry_list(cable_geometry)%external_conductor_radius
    write(*,*)'External_dielectric_radius	=',cable_geometry_list(cable_geometry)%external_dielectric_radius
    write(*,*)'External_dielectric_permittivity	=',cable_geometry_list(cable_geometry)%external_dielectric_permittivity

! write parameters
    write(*,*)'Number of parameters		=',cable_geometry_list(cable_geometry)%n_parameters
    
    n_cols=cable_geometry_list(cable_geometry)%n_parameters
    write(*,*)'Parameter list:'
    write(*,*)(cable_geometry_list(cable_geometry)%parameters(col),	&
                                   col=1,cable_geometry_list(cable_geometry)%n_parameters)
				   
    n_rows=cable_geometry_list(cable_geometry)%n_conductors
    n_cols=n_rows
    
! write Sc
    write(*,*)'Shielded conductor flag vector (Sc):'
    write(*,8000)(cable_geometry_list(cable_geometry)%Sc(row),row=1,n_rows)
    
! write Tv
    write(*,*)'Voltage reference matrix (Tv):'
    do row=1,n_rows	   
      write(*,8000)(cable_geometry_list(cable_geometry)%Tv(row,col),col=1,n_cols)
      write(*,*)
    end do
    
! write L_internal
    write(*,*)'Internal inductance matrix (L_internal):'
    do row=1,n_rows	   
      write(*,8010)(cable_geometry_list(cable_geometry)%L_internal(row,col),col=1,n_cols)
      write(*,*)
    end do
    
! write C_internal
    write(*,*)'Internal capacitance matrix (C_internal):'
    do row=1,n_rows	   
      write(*,8010)(cable_geometry_list(cable_geometry)%C_internal(row,col),col=1,n_cols)
      write(*,*)
    end do
    
    write(*,*)'________________________________________________________'
    write(*,*)

8000 format(100I4)
8010 format(100E15.6) 
  end do ! next cable geometry
  


  CALL write_line('FINISHED: Output_cable_geometry',0,output_to_screen_flag)

    
  RETURN
    
END SUBROUTINE Output_cable_geometry


!
!SUBROUTINE Output_cable_route
!
! NAME
!     SUBROUTINE Output_cable_route
!
! DESCRIPTION
!      
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/11/2012 CJS
!
!
SUBROUTINE Output_cable_route()

USE TLM_general
USE geometry_types
USE Cables
USE Mesh
USE File_information

IMPLICIT NONE

! local variables

  integer cable
  integer first_cable
  integer last_cable
  integer segment
  integer i
  integer number_of_cell_segments

! START

  CALL write_line('CALLED: Output_cable_route',0,output_to_screen_flag)

  write(*,*)'Number of cables=',n_cables
  
10 continue

  write(*,*)
  write(*,*)'Enter the number of the cable to output or 0 to output all of then'
  read(*,*)cable
  
  if (cable.eq.0) then
    first_cable=1
    last_cable=n_cables
  else 
  
    if ( (cable.lt.0).OR.(cable.gt.n_cables) ) then
      write(*,*)'Cable is outside the available range'
      RETURN
    end if
    
    first_cable=cable
    last_cable=cable
  end if
  
  do cable=first_cable,last_cable
  
    write(*,*)
    write(*,*)'Cable number			=',cable
    write(*,*)
    write(*,*)'Cable geometry number		=',cable_list(cable)%cable_geometry_number
    write(*,*)'End 1 junction number		=',cable_list(cable)%junction_1  
    write(*,*)'End 2 junction number		=',cable_list(cable)%junction_2  
    write(*,*)'Number of lines on cable route	=',cable_list(cable)%n_lines  
    
    write(*,*)'Cable line list:'
    write(*,8000)(cable_list(cable)%line_list(i),i=1,cable_list(cable)%n_lines) 
    
    
    write(*,*)'Number of segments on cable route	=',cable_list(cable)%number_of_cable_segments  
    write(*,*)'Cable segment list:'
    write(*,*)'     Segment End 1             Segment End 2'
    write(*,*)'  i    j    k   point       i    j    k   point'
    
    do segment=1,cable_list(cable)%number_of_cable_segments
      write(*,8010)cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%i,	&
                   cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%j,	&
                   cable_list(cable)%cable_segment_list(segment)%segment_point(1)%cell%k,' ',	&
                   face_string(cable_list(cable)%cable_segment_list(segment)%segment_point(1)%point),	&
		   '    ',	&
                   cable_list(cable)%cable_segment_list(segment)%segment_point(2)%cell%i,	&
                   cable_list(cable)%cable_segment_list(segment)%segment_point(2)%cell%j,	&
                   cable_list(cable)%cable_segment_list(segment)%segment_point(2)%cell%k,' ',	&
                   face_string(cable_list(cable)%cable_segment_list(segment)%segment_point(2)%point)
                   
    end do ! next segment%i

8000 format(100I8)
8010 format(3I5,A,A,A,3I5,A,A)


    number_of_cell_segments=cable_list(cable)%number_of_cable_segments

    if (number_of_cell_segments.gt.0) then
    
! open and write line to vtk format file
      CALL open_vtk_file(cable_mesh_file_unit,cable_mesh_file_extension,cable) 
      
      CALL write_line_mesh_list_vtk(cable_mesh_file_unit,number_of_cell_segments,	&
                        cable_list(cable)%cable_segment_list)
      
      CALL close_vtk_file(cable_mesh_file_unit) 
    
    end if ! number_of_cell_segments.gt.0

  
  end do ! next cable


  CALL write_line('FINISHED: Output_cable_route',0,output_to_screen_flag)

    
  RETURN
    
END SUBROUTINE Output_cable_route

!
!SUBROUTINE Output_bundle_geometry
!
! NAME
!     SUBROUTINE Output_bundle_geometry
!
! DESCRIPTION
!      
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/11/2012 CJS
!
!
SUBROUTINE Output_bundle_geometry()

USE TLM_general
USE Cables
USE File_information

IMPLICIT NONE

! local variables

  integer bundle_geometry
  integer row
  integer col,n_cols
  integer i

! START

  CALL write_line('CALLED: Output_bundle_segment_geometry',0,output_to_screen_flag)
  
  write(*,*)'Number of cable bundle geometries=',n_bundle_segment_geometries

  write(*,*)
  write(*,*)'Enter the number of the bundle geometry to output or 0 to output a summary '
  read(*,*)bundle_geometry
  
  if (bundle_geometry.eq.0) then
    
    write(*,*)'geometry  n_cables  n_conductors	   Cable list'
    write(*,*)' number   '
      
    do bundle_geometry=1,n_bundle_segment_geometries
       
      n_cols= bundle_segment_geometry_list(bundle_geometry)%n_cables
      write(*,8000)bundle_geometry,bundle_segment_geometry_list(bundle_geometry)%n_cables,	&
                   bundle_segment_geometry_list(bundle_geometry)%n_conductors,'           ',	&
		  (bundle_segment_geometry_list(bundle_geometry)%cable_list(i),i=1,n_cols)
       
    end do
    
  end if
  
  if ( (bundle_geometry.lt.1).OR.(bundle_geometry.gt.n_bundle_segment_geometries) ) then
    write(*,*)'Bundle geometry number should greater than 0 and less than ',n_bundle_segment_geometries
    RETURN
  end if

  write(*,*)'geometry  n_cables  n_conductors	   Cable list'
  write(*,*)' number   '
  n_cols= bundle_segment_geometry_list(bundle_geometry)%n_cables
  write(*,8000)bundle_geometry,bundle_segment_geometry_list(bundle_geometry)%n_cables,      &
  	       bundle_segment_geometry_list(bundle_geometry)%n_conductors,'	      ',    &
    	      (bundle_segment_geometry_list(bundle_geometry)%cable_list(i),i=1,n_cols)

  n_cols=bundle_segment_geometry_list(bundle_geometry)%n_conductors

  write(*,*)'L'    
  do row=1,n_cols
    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%L(row,col),col=1,n_cols)
  end do ! next row
    
  write(*,*)'C'    
  do row=1,n_cols
    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%C(row,col),col=1,n_cols)
  end do ! next row
    
  write(*,*)'R'    
  do row=1,n_cols
    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%R(row,col),col=1,n_cols)	 
  end do ! next row
    
  write(*,*)'Zlink'    
  do row=1,n_cols
    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%Zlink(row,col),col=1,n_cols)	 
  end do ! next row
    
  write(*,*)'Ylink'    
  do row=1,n_cols
    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%Ylink(row,col),col=1,n_cols)	 
  end do ! next row
    
  write(*,*)'ZLstub'    
  do row=1,n_cols
    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%ZLstub(row,col),col=1,n_cols)	 
  end do ! next row
    
  write(*,*)'Yf'    
  do row=1,n_cols
    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%Yf(row,col),col=1,n_cols)	 
  end do ! next row
    
  write(*,*)'Tv'    
  do row=1,n_cols
    write(*,8020)(bundle_segment_geometry_list(bundle_geometry)%Tv(row,col),col=1,n_cols)     
  end do ! next row
    
  write(*,*)'Sc'    
  do row=1,n_cols
    write(*,*)bundle_segment_geometry_list(bundle_geometry)%SC(row)
  end do ! next row

! Not allocated yet but should be included in future    
!  write(*,*)'Zlink'    
!  do row=1,n_cols
!    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%Zlink(row,col),col=1,n_cols)	 
!  end do ! next row
!    
!  write(*,*)'Ylink'    
!  do row=1,n_cols
!    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%Ylink(row,col),col=1,n_cols)	 
!  end do ! next row
!    
!  write(*,*)'ZLstub'    
!  do row=1,n_cols
!    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%ZLstub(row,col),col=1,n_cols)	 
!  end do ! next row
!    
!  write(*,*)'Yf'    
!  do row=1,n_cols
!    write(*,8010)(bundle_segment_geometry_list(bundle_geometry)%Yf(row,col),col=1,n_cols)	 
!  end do ! next row
  
8000  format(I7,I9,I12,A,100I5)
8010  format(100E14.4)
8020  format(100I2)

  CALL write_line('FINISHED: Output_bundle_geometry',0,output_to_screen_flag)

    
  RETURN
    
END SUBROUTINE Output_bundle_geometry

!
!SUBROUTINE Output_junction_specification
!
! NAME
!     SUBROUTINE Output_junction_specification
!
! DESCRIPTION
!      
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/11/2012 CJS
!
!
SUBROUTINE Output_junction_specification()

USE TLM_general
USE cell_parameters
USE Cables
USE File_information

IMPLICIT NONE

! local variables

  integer bundle_junction
  integer ix,iy,iz,face
  integer segment
  integer n_internal,n_external
  integer row,col
  integer i
  integer filter

! START

  CALL write_line('CALLED: Output_junction_specification',0,output_to_screen_flag)

  write(*,*)'Number of bundle junctions=',n_cell_centre_junctions
  
10 continue

  write(*,*)
  write(*,*)'Enter the number of the cell centre junction to output or 0 to output a summary '
  read(*,*)bundle_junction
  
  if (bundle_junction.eq.0) then
  
    write(*,*)
    write(*,*)'junction   ix   iy   iz n_internal segment segment segment segment segment segment'
    write(*,*)'                                     xmin    ymin    zmin    xmax    ymax    zmax '
    
    do bundle_junction=1,n_cell_centre_junctions
    
     ix=cell_centre_junction_list(bundle_junction)%cell_point%cell%i
     iy=cell_centre_junction_list(bundle_junction)%cell_point%cell%j
     iz=cell_centre_junction_list(bundle_junction)%cell_point%cell%k
     n_internal=cell_centre_junction_list(bundle_junction)%n_internal_connection_nodes
    
      write(*,8000)bundle_junction,ix,iy,iz,n_internal,	&
      cell_centre_junction_list(bundle_junction)%segment_list(1),	&
      cell_centre_junction_list(bundle_junction)%segment_list(2),	&
      cell_centre_junction_list(bundle_junction)%segment_list(3),	&
      cell_centre_junction_list(bundle_junction)%segment_list(4),	&
      cell_centre_junction_list(bundle_junction)%segment_list(5),	&
      cell_centre_junction_list(bundle_junction)%segment_list(6)
8000  format(I7,I7,2I5,I11,6I8)
   
    end do
  end if
  
  if ( (bundle_junction.lt.1).OR.(bundle_junction.gt.n_cell_centre_junctions) ) then
    write(*,*)'Bundle junction number should greater than 0 and less than ',n_cell_centre_junctions
    RETURN
  end if

  ix=cell_centre_junction_list(bundle_junction)%cell_point%cell%i
  iy=cell_centre_junction_list(bundle_junction)%cell_point%cell%j
  iz=cell_centre_junction_list(bundle_junction)%cell_point%cell%k

  n_internal=cell_centre_junction_list(bundle_junction)%n_internal_connection_nodes

  write(*,*)' Bundle junction number',bundle_junction
  write(*,*)
  write(*,*)'ix=',ix,' iy=',iy,' iz=',iz
  write(*,*)
  write(*,*)'Number of internal connection nodes=',n_internal
  write(*,*)
  do face=1,6
  
    segment=cell_centre_junction_list(bundle_junction)%segment_list(face)

    write(*,*)'Face:',face_string(face),' Segment number=',segment
    if (segment.ne.0) then
    
      n_external=cell_centre_junction_list(bundle_junction)%n_external_conductors(face)
    
      write(*,*)'Number of external cables=',bundle_segment_list(segment)%n_cables
      write(*,*)'Number of external conductors=',bundle_segment_list(segment)%n_conductors,' =',n_external
      
      write(*,*)'Cable list:'
      
      do i=1,bundle_segment_list(segment)%n_cables
        write(*,*)bundle_segment_list(segment)%cable_list(i)
      end do
      
      write(*,*)'P matrix (n_internal rows x n_external columns):'
      write(*,*)
      
      do row=1,n_internal
      
        write(*,8010)(cell_centre_junction_list(bundle_junction)%P_matrix_list(face)%P(row,col),col=1,n_external)
        write(*,*)
8010	format(1000I2)
	
      end do ! next internal connection node
      
    end if ! there is a cable bundle on this face

  end do ! next face
  
  face=7
  n_external=cell_centre_junction_list(bundle_junction)%n_external_conductors(face)
     
  if (n_external.ne.0) then
    write(*,*)'Internal impedance data'	    
    write(*,*)'Number of impedance filters=',n_external	    
    write(*,*)'P matrix (n_internal rows x n_external columns):'
    write(*,*)
    do row=1,n_internal
    
      write(*,8010)(cell_centre_junction_list(bundle_junction)%P_matrix_list(face)%P(row,col),col=1,n_external)
      write(*,*)

    end do ! next internal connection node
    do filter=1,cell_centre_junction_list(bundle_junction)%n_internal_impedance_filters

      write(*,*)'! Sfilter number',filter	    
      CALL write_Sfilter(cell_centre_junction_list(bundle_junction)%Sfilter(filter),0)
      write(*,*)'! Z_f',filter	  
      write(*,*)cell_centre_junction_list(bundle_junction)%Z_f(filter)
      write(*,*)'! Zfilter number',filter	 
      CALL write_Zfilter(cell_centre_junction_list(bundle_junction)%Zfilter(filter),0)
    
    end do ! next filter
  end if
  
  CALL write_line('FINISHED: Output_junction_specification',0,output_to_screen_flag)

    
  RETURN
    
END SUBROUTINE Output_junction_specification

!
!SUBROUTINE Output_face_junction_specification
!
! NAME
!     SUBROUTINE Output_face_junction_specification
!
! DESCRIPTION
!      
!
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/11/2012 CJS
!
!
SUBROUTINE Output_face_junction_specification()

USE TLM_general
USE Cables
USE File_information

IMPLICIT NONE

! local variables

  integer face_junction
  integer ix,iy,iz,face
  integer segment
  integer n_internal,n_external
  integer row,col
  integer i

! START

  CALL write_line('CALLED: Output_face_junction_specification',0,output_to_screen_flag)

  write(*,*)'Number of face junctions=',n_face_junctions
  
10 continue

  write(*,*)
  write(*,*)'Enter the number of the face junction to output or 0 to output a summary '
  read(*,*)face_junction
  
  if (face_junction.eq.0) then
  
    write(*,*)
    write(*,*)'junction   ix   iy   iz n_internal segment1 segment2 '
    
    do face_junction=1,n_face_junctions
    
     ix=face_junction_list(face_junction)%cell_point%cell%i
     iy=face_junction_list(face_junction)%cell_point%cell%j
     iz=face_junction_list(face_junction)%cell_point%cell%k
     n_internal=face_junction_list(face_junction)%n_internal_connection_nodes
    
      write(*,8000)face_junction,ix,iy,iz,n_internal,	&
      face_junction_list(face_junction)%segment_list(1),	&
      face_junction_list(face_junction)%segment_list(2)
8000  format(I7,3I5,I11,2I9)
   
    end do
  end if
  
  if ( (face_junction.lt.1).OR.(face_junction.gt.n_face_junctions) ) then
    write(*,*)'Bundle junction number should greater than 0 and less than ',n_face_junctions
    RETURN
  end if

  ix=face_junction_list(face_junction)%cell_point%cell%i
  iy=face_junction_list(face_junction)%cell_point%cell%j
  iz=face_junction_list(face_junction)%cell_point%cell%k

  n_internal=face_junction_list(face_junction)%n_internal_connection_nodes

  write(*,*)' Bundle junction number',face_junction
  write(*,*)
  write(*,*)'ix=',ix,' iy=',iy,' iz=',iz
  write(*,*)
  write(*,*)'Number of internal connection nodes=',n_internal
  write(*,*)
  do face=1,2
  
    segment=face_junction_list(face_junction)%segment_list(face)

    write(*,*)'Face:',face,' Segment number=',segment
    if (segment.ne.0) then
    
      n_external=face_junction_list(face_junction)%n_external_conductors(face)
    
      write(*,*)'Number of external cables=',bundle_segment_list(segment)%n_cables
      write(*,*)'Number of external conductors=',bundle_segment_list(segment)%n_conductors,' =',n_external
      
      write(*,*)'Cable list:'
      
      do i=1,bundle_segment_list(segment)%n_cables
        write(*,*)bundle_segment_list(segment)%cable_list(i)
      end do
      
      write(*,*)'P matrix (n_internal rows x n_external columns):'
      write(*,*)
      
      do row=1,n_internal
      
        write(*,8010)(face_junction_list(face_junction)%P_matrix_list(face)%P(row,col),col=1,n_external)
        write(*,*)
8010	format(1000I2)
	
      end do ! next internal connection node
      
      write(*,*)'BC vector (n_internal rows ):'
      write(*,*)
      
      do row=1,n_internal
      
        write(*,8020)face_junction_list(face_junction)%BC(row)
        write(*,*)
8020	format(I4)
	
      end do ! next internal connection node
      
    end if ! there is a cable bundle on this face

  end do ! next face

  CALL write_line('FINISHED: Output_face_junction_specification',0,output_to_screen_flag)

    
  RETURN
    
END SUBROUTINE Output_face_junction_specification

