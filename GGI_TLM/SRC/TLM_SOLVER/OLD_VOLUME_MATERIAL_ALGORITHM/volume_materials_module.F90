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
!     MODULE TLM_volume_materials
!
! DESCRIPTION
!     volume_materials data relating to the TLM solution
!     
! COMMENTS
!     
!
! HISTORY
!
!     started 5/09/2012 CJS
!
!
MODULE TLM_volume_materials


USE filter_types

IMPLICIT NONE

TYPE::volume_material_type

  integer		:: type
  integer		:: n_volumes
  integer,allocatable	:: volume_list(:)
  
  REAL*8 		:: fmin,fmax
  
  TYPE(Sfilter)		:: eps_S
  TYPE(Sfilter)		:: mu_S  
  TYPE(Zfilter)		:: eps_Z
  TYPE(Zfilter)		:: mu_Z
  REAL*8		:: eps_f,sigma_e
  REAL*8		:: mu_f,sigma_m
  
  REAL*8	       	:: Ys
  REAL*8	       	:: Zs
  REAL*8	       	:: Ge
  REAL*8	       	:: Rm

END TYPE volume_material_type

  integer,parameter	:: volume_material_type_PEC=1

  integer,parameter	:: volume_material_type_PMC=2

  integer,parameter	:: volume_material_type_DISPERSIVE=3

  integer				  :: n_volume_materials
  type(volume_material_type),allocatable :: volume_material_list(:)
  
  integer 	:: n_volume_material_cells
  integer 	:: volume_material_storage
  
  type(Zfilter_response),allocatable	:: volume_material_eps_filter_data(:) 
  type(Zfilter_response),allocatable	:: volume_material_mu_filter_data(:) 

  real*8,allocatable			:: Vs(:)
  real*8,allocatable			:: Vo(:)
  
END MODULE TLM_volume_materials
