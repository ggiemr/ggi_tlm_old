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
!       subroutine cmatmul(A,ar,ac,B,br,bc,C,matdim)
!       subroutine ctranspose(u,nr,nc,ut,nmodes)
!       subroutine cmatvmul(A,ar,ac,B,br,C,matdim)
!       subroutine csvd(A,ar,ac,matdim,U,GAMMA,VT) 
!       subroutine csvd_invert(A,ar,ac,AI,matdim) 
!
! __________________________________________________
!
!
       subroutine ctranspose(u,nr,nc,ut,nmodes)
       
       integer nmodes
       complex*16 u(nmodes,nmodes),ut(nmodes,nmodes)
       integer nr,nc
       integer r,c
       
       do r=1,nr
         do c=1,nc
	   ut(c,r)=conjg(u(r,c))
	 end do
       end do

       return
       end
!
! __________________________________________________
!
!
       subroutine  cmatmul(A,ar,ac,B,br,bc,C,matdim)

       
       integer ar,ac,br,bc,matdim
       complex*16 A(matdim,matdim)
       complex*16 B(matdim,matdim)
       complex*16 C(matdim,matdim)
       
       complex*16 sum
       
       integer row,col,i
  
       if (ac.ne.br) then
         print*,'matrix dimension error in cmatmul'
	 print*,'ac=',ac
	 print*,'br=',br
	 stop
       end if
  
       do row=1,ar
         do col=1,bc
	   c(row,col)=0.0
           do i=1,ac
	     c(row,col)=c(row,col)+a(row,i)*b(i,col)
	   end do
	 end do
       end do
             
       return
       end
!
! __________________________________________________
!
!
       subroutine  cmatvmul(A,ar,ac,B,br,C,matdim)

! multiply complex vector by a complex matrix
       
       integer ar,ac,br,matdim
       complex*16 A(matdim,matdim)
       complex*16 B(matdim)
       complex*16 C(matdim)
       
       integer i,j
       complex*16 sum
              
       if (br.ne.ac) then
         print*,'dimension error in cmatvmul'
	 print*,'number of columns in a=',ac
	 print*,'number of rows in b   =',br
	 stop
       end if
         
       do i=1,ar
         sum=(0.0,0.0)
         do j=1,ac
	   sum=sum+A(i,j)*B(j)
	 end do
	 C(i)=sum
       end do	 
       
       return
       end
!
! __________________________________________________
!
!
       subroutine csvd(A,ar,ac,matdim,U,GAMMA,VT) 
                  
! invert matrix using SVD implemented in LAPACK

       IMPLICIT NONE
       
       integer ar,ac,matdim
       complex*16 A(matdim,matdim)
       DOUBLE PRECISION   GAMMA(matdim)
       COMPLEX*16      U(matdim,matdim),VT(matdim,matdim)
       
!      LAPACK ARGUMENTS

       CHARACTER          JOBU, JOBVT
       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
       
       DOUBLE PRECISION   RWORK(5*matdim)
       COMPLEX*16         WORK(5*matdim)
       
! LOCAL VARIABLES
       integer row,col       

! START
       
! perform singular value decomposition on A
      
      JOBU='A'
      JOBVT='A'
      M=ar
      N=ac
      LDA=matdim
      LDU=matdim
      LDVT=matdim
      LWORK=5*matdim
      
      CALL ZGESVD(JOBU,JOBVT,M,N,A,LDA,GAMMA,U,LDU,VT,LDVT,  &
                       WORK,LWORK,RWORK,INFO)
            
       if(info.ne.0) then
         if(info.lt.0) then
	   print*,'Illegal argument to zgesvd',-info
	   stop
	 else if(info.gt.0) then
	   print*,'ZGESVD did not converge, info=',info
	   stop
	 end if
       end if  
	              
      return
      end
!
! __________________________________________________
!
!
       subroutine csvd_invert(A,ar,ac,AI,matdim) 

! invert matrix using SVD implemented in LAPACK

       IMPLICIT NONE
       
       integer ar,ac,matdim
       complex*16 A(matdim,matdim)
       complex*16 AI(matdim,matdim)
                     
       DOUBLE PRECISION  GAMMA(matdim)
       COMPLEX*16      U(matdim,matdim),UT(matdim,matdim)
       DOUBLE PRECISION  GAMMA_I(matdim)
       COMPLEX*16  V(matdim,matdim),VT(matdim,matdim)
       
       COMPLEX*16      ITEST(matdim,matdim)
       
       complex*16 tm1(matdim,matdim),tm2(matdim,matdim)
       
       integer row,col

! START
     
! zero matrices 
     
       do row=1,matdim
         do col=1,matdim
	   U(row,col)=(0D0,0D0)
	   UT(row,col)=(0D0,0D0)
	   V(row,col)=(0D0,0D0)
	   VT(row,col)=(0D0,0D0)
	   tm1(row,col)=(0D0,0D0)
	   tm2(row,col)=(0D0,0D0)
	 end do
	 GAMMA(row)=0D0
	 GAMMA_I(row)=0D0
       end do
       
! generate unit matrix 
     
       do row=1,matdim
         do col=1,matdim
	   itest(row,col)=0D0
	 end do
	 itest(row,row)=1D0
       end do
              
! save matrix a for checks
       do row=1,matdim
         do col=1,matdim
	   tm1(row,col)=a(row,col)
	 end do
       end do

! perform singular value decomposition on A
       
       call csvd(A,ar,ac,matdim,U,GAMMA,VT)
               
! calculate gamma inverse
       do row=1,ac
         if (abs(gamma(row)).ne.(0D0)) then
           gamma_i(row)=1D0/gamma(row)
         else 
           gamma_i(row)=0D0
	 end if
       end do
              
! calculate transposes of U and VT
       call ctranspose(U,ar,ac,UT,matdim)
       call ctranspose(VT,ac,ac,V,matdim)       
                
! pre multiply U transpose by gamma inverse

       do row=1,ac
         do col=1,ar
	   tm2(row,col)=UT(row,col)*gamma_i(row)
	 end do
       end do
                     
! pre multiply result by V to give ai, the inverse of a
                    	
       call cmatmul(V,ac,ac,tm2,ac,ar,AI,matdim)	 
              
! restore a
       do row=1,matdim
         do col=1,matdim
	   a(row,col)=tm1(row,col)
	 end do
       end do
       
       return
       end
