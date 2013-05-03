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
! SUBROUTINE run_TLM
!
! NAME
!     run_TLM
!
! DESCRIPTION
!     Main update loop for the TLM solution
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 14/08/2012 CJS
!
!
SUBROUTINE run_TLM

USE TLM_general
USE mesh
USE TLM_excitation
USE File_information

USE TLM_output   ! temp for testing

IMPLICIT NONE

! local variables

  character*256	:: ipline
  
  integer	:: op_time_period

! START
  
  CALL write_line('CALLED: run_TLM',0,output_to_screen_flag)
  
  timestepping_output_to_screen_flag=.FALSE.
!  timestepping_output_to_screen_flag=.TRUE.

  if (rank.eq.0) then
    write(info_file_unit,*)'' 
    write(info_file_unit,*)'Nx=',nx 
    write(info_file_unit,*)'Ny=',ny 
    write(info_file_unit,*)'Nz=',nz 
    write(info_file_unit,*)'' 
    write(info_file_unit,*)'Total number of cells=',nx*ny*nz
    write(info_file_unit,*)'' 
    write(info_file_unit,*)'Number of timesteps=',n_timesteps
    write(info_file_unit,*)'' 
    write(info_file_unit,*)'TLM solution Started:'
    call write_date_and_time(info_file_unit)
    write(info_file_unit,*)'' 
  end if
  
#if defined(MPI)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
#endif

  if (rank.eq.0) then
  
#if defined(SEQ)
    call system("ps u -C GGI_TLM_SEQ > GGI_TLM_memory_usage.txt ")
#elif defined(MPI)
    call system("ps u -C GGI_TLM_MPI > GGI_TLM_memory_usage.txt ")
#endif
    
    write(info_file_unit,'(A)')"________________________________________________________________________"
    write(info_file_unit,'(A)')""
    write(info_file_unit,'(A)')"Memory Usage:"
    write(info_file_unit,'(A)')""
!    write(info_file_unit,'(A)')"USER       PID %CPU %MEM    VSZ   RSS TTY      STAT START   TIME COMMAND"
    
    open(unit=scratch_file_unit,file='GGI_TLM_memory_usage.txt')
5   CONTINUE
    read(scratch_file_unit,'(A256)',end=6)ipline
    write(info_file_unit,'(A)')trim(ipline)
    
    GOTO 5
    
6   CONTINUE

    write(info_file_unit,'(A)')""
    write(info_file_unit,'(A)')"________________________________________________________________________"
    write(info_file_unit,'(A)')""
    
  end if

  do timestep=1,n_timesteps
  
!    CALL write_line_integer('Timestep',timestep,0,output_to_screen_flag)

    op_time_period=min(100,10**INT(log10(dble(timestep))))
    
    if (rank.eq.0) then 
      if (timestep.EQ.1) then
	write(6,'(A9,I10,A4,I10)',advance='no')'Timestep ',timestep,' of ',n_timesteps
	flush(6)
      else if (timestep.EQ.n_timesteps) then
	write(6,'(A)',advance='no')char(13)
	write(6,'(A9,I10,A4,I10)')'Timestep ',timestep,' of ',n_timesteps
	flush(6)
      else if ( mod(timestep,op_time_period).EQ.0 ) then
	write(6,'(A)',advance='no')char(13)
	write(6,'(A9,I10,A4,I10)',advance='no')'Timestep ',timestep,' of ',n_timesteps
	flush(6)    
      end if
    end if
    
    time=(timestep-1)*dt
    
    CALL excitation()
  
    CALL scatter()
  
    CALL cell_output()
    
    time=time+dt/2d0
    
    if (np.gt.1) then 
    
      CALL TLM_parallel_pass_data_1()
      
      CALL Cable_parallel_pass_data_1()
      
    end if
    
    CALL mode_stir_surfaces()
    
    CALL connect()
  
    CALL face_output()
    
    CALL cable_output()
    
    CALL outer_boundary()
    
    CALL wrap_outer_boundary()
    
    if (np.gt.1) then 
    
      CALL TLM_parallel_pass_data_2()
      
      CALL Cable_parallel_pass_data_2()
      
    end if
  
  end do ! next timestep

  if (rank.eq.0) then
    write(info_file_unit,*)'' 
    write(info_file_unit,*)'TLM solution Finished:'
    call write_date_and_time(info_file_unit)
    write(info_file_unit,*)'' 
  end if
  
  CALL write_line('FINISHED: run_TLM',0,output_to_screen_flag)

  RETURN

END SUBROUTINE run_TLM
