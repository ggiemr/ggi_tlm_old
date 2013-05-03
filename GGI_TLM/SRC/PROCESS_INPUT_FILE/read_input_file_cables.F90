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
! SUBROUTINE read_input_file_cables
!
! NAME
!     read input file cables
!
! DESCRIPTION
!     read cable information from the input file
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 6/08/2012 CJS
!
!
SUBROUTINE read_input_file_cables

USE TLM_general
USE file_information
USE geometry

IMPLICIT NONE

! local variables

character*256	:: input_line

! START
  
  CALL write_line('CALLED: read_input_file_cables',0,output_to_screen_flag)
    
  rewind(input_file_unit)
  
10  CONTINUE
    
! read line from input file
    read(input_file_unit,'(A)',end=100)input_line

! convert text to lower case
    CALL convert_to_lower_case(input_line,256)

    if (input_line.EQ.'cable_geometry_list') then
    
      CALL read_cable_geometry_list()
          
    else if (input_line.EQ.'cable_list') then
    
      CALL read_cable_list()
          
    else if (input_line.EQ.'cable_junction_list') then
    
      CALL read_cable_junction_list()
          
    else if (input_line.EQ.'cable_output_list') then
    
      CALL read_cable_output_list()
      
    end if
    
    GOTO 10
    
100 CONTINUE
  
  CALL write_line('FINISHED: read_input_file_cables',0,output_to_screen_flag)
  
  RETURN
  
END SUBROUTINE read_input_file_cables
