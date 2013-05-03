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
! subroutine dmatmul(a,ar,ac,b,br,bc,ans,matdim)
! subroutine dmatdmul(a,ar,ac,d,dr,ans,matdim)
! subroutine dmatvmul(a,ar,ac,v,vr,ans,matdim)
! subroutine dtranspose(a,nr,nc,at,matdim)
! subroutine dmake_symmetric(a,nr,nc,at,matdim)
! subroutine dcheck_symmetry(a,nr,nc,result,matdim)
! subroutine dludcmp(a,n,indx,d,matdim)
! subroutine dlubksub(a,n,indx,b,matdim)
! subroutine dludinvert(a,n,y,matdim) 
! subroutine dmatrix_reset(nr,nc,a)
! subroutine dadd(A,ar,ac,B,br,bc,C,matdim)
! subroutine dsub(A,ar,ac,B,br,bc,C,matdim)
! subroutine dsvd(A,ar,ac,matdim,U,GAMMA,VT) 
! subroutine dsvd_invert(A,ar,ac,AI,matdim) 
! subroutine deig(A,ar,GAMMA,matdim) 
! subroutine check_difference(m1,nr,nc,m2,matdim)

!_____________________________________________________________________
!
! NAME
!    dmatmul
!
! DESCRIPTION
!    multiply real*8 matrices
!
!
! SEE ALSO
!   
!
! HISTORY
!
!     started 18/05/05 CJS
!
! COMMENTS
!   
!
  subroutine dmatmul(a,ar,ac,b,br,bc,ans,matdim)

! Modules used

IMPLICIT NONE

!Argument list variables

  integer matdim
  real*8 a(matdim,matdim),b(matdim,matdim),       &
  	     ans(matdim,matdim)
  integer ar,ac,br,bc

!Local variables

  real*8 sum
  integer row,col,i
 
! START

  if (ac.ne.br) then
    print*,'matrix dimension error in dmatmul'
    stop
  end if
  do row=1,ar
    do col=1,bc 
      sum=0d0
      do i=1,ac
  	sum=sum+a(row,i)*b(i,col)
      end do
      ans(row,col)=sum
    end do
  end do

! END
 
  return
  end
!
!_____________________________________________________________________
!
! NAME
!    dmatdmul
!
! DESCRIPTION
!    post multiply a real*8 matrix by a diagonal matrix
!
!
! SEE ALSO
!   
!
! HISTORY
!
!     started 18/05/05 CJS
!
! COMMENTS
!   
!
 subroutine dmatdmul(a,ar,ac,d,dr,ans,matdim)

! Modules used

IMPLICIT NONE

!Argument list variables

 integer matdim
 real*8 a(matdim,matdim),d(matdim),ans(matdim,matdim)
 integer ar,ac,dr

!Local variables

 integer r,c
 
! START

 if (ac.ne.dr) then
   print*,'matrix dimension error in dmatdmul'
   stop
 end if
 do r=1,ar
   do c=1,dr
     ans(r,c)=a(r,c)*d(c)
   end do
 end do

! END
 
 return
 end 
!
!_____________________________________________________________________
!
! NAME
!    dmatvmul
!
! DESCRIPTION
!    post multiply a complex vector by a real*8 matrix
!
!
! SEE ALSO
!   
!
! HISTORY
!
!     started 18/05/05 CJS
!
! COMMENTS
!   
!
 subroutine dmatvmul(a,ar,ac,v,vr,ans,matdim)

! Modules used

IMPLICIT NONE

!Argument list variables

 integer matdim
 real*8 a(matdim,matdim),v(matdim),ans(matdim)
 integer ar,ac,vr

!Local variables

 real*8 sum
 integer r,c
 
! START

 if (ac.ne.vr) then
   print*,'matrix dimension error'
   stop
 end if
 do r=1,ar
   sum=0d0
   do c=1,vr
     sum=sum+a(r,c)*v(c)
   end do
   ans(r)=sum
 end do

! END
 
 return
 end
!
!_____________________________________________________________________
!
! NAME
!    dtranspose
!
! DESCRIPTION
!    
!    calculate  transpose of a matrix
!
! SEE ALSO
!   
!
! HISTORY
!
!     started 18/05/05 CJS
!
! COMMENTS
!
!
 subroutine dtranspose(a,nr,nc,at,matdim)

! Modules used

IMPLICIT NONE

!Argument list variables

 integer matdim,nr,nc
 real*8 a(matdim,matdim),at(matdim,matdim)

!Local variables

 integer r,c
 
! START

 do r=1,nr
   do c=1,nc
    at(c,r)=a(r,c)
   end do
 end do
 
! END
 
 return
 end
!
!_____________________________________________________________________
!
! NAME
!    dmake_symmetric
!
! DESCRIPTION
!    
!    make a matrix symmetric by averaging it with its transpose
!
! SEE ALSO
!   
!
! HISTORY
!
!     started 18/05/05 CJS
!
! COMMENTS
!
!
 subroutine dmake_symmetric(a,nr,nc,at,matdim)

! Modules used

IMPLICIT NONE

!Argument list variables

 integer matdim,nr,nc
 real*8 a(matdim,matdim),at(matdim,matdim)

!Local variables

 integer r,c
 
! START

 do r=1,nr
   do c=1,nc
    at(r,c)=(a(r,c)+a(c,r))/2d0
   end do
 end do
 
! END
 
 return
 end
!
!_____________________________________________________________________
!
! NAME
!    dcheck_symmetry
!
! DESCRIPTION
!    
!   check symmetry of a matrix
!
! SEE ALSO
!   
!
! HISTORY
!
!     started 18/05/05 CJS
!
! COMMENTS
!
!
 subroutine dcheck_symmetry(a,nr,nc,result,matdim)

! Modules used

IMPLICIT NONE

!Argument list variables

 integer matdim,nr,nc
 real*8 a(matdim,matdim)
 logical result

!Local variables

 integer r,c
 
! START

 result=.TRUE.

 if (nr.ne.nc) then
   result=.FALSE.
   return
 end if

 do r=1,nr
   do c=1,nc
     if (a(c,r).NE.a(r,c)) then
       result=.FALSE.
       return
     end if
   end do
 end do
 
! END
 
 return
 end
!
!
!_____________________________________________________________________
!
! NAME
!    dludcmp
!
! DESCRIPTION
!    
!    calculate LU decomposition of a real*8 matrix
!
! SEE ALSO
!   Numerical recipes
!
! HISTORY
!
!     started 18/05/05 CJS
!
! COMMENTS
!
!

 subroutine dludcmp(a,n,indx,d,matdim)
 
! Modules used

IMPLICIT NONE

!Argument list variables

 integer matdim
 real*8 a(matdim,matdim)
 integer n,indx(matdim)
 real*8 d

!Local variables

 integer i,j,k,imax
 real*8 big,tiny,vv(matdim),temp
 real*8 sum,dum 
 
! START

 tiny=1d-20
 d=1d0
 imax=0
 do i=1,n
   big=0.0
   do j=1,n
     temp=abs(a(i,j))
     if (temp.gt.big) big=temp
   end do
   if (big.eq.0d0) then
     print*,'singular matrix in call to subroutine ludcmp'
     stop
   end if
   vv(i)=1d0/big
 end do

 do j=1,n
   do i=1,j-1
     sum=a(i,j)
     do k=1,i-1
       sum=sum-a(i,k)*a(k,j)
     end do
     a(i,j)=sum
   end do
 big=0d0
 do i=j,n
   sum=a(i,j)
   do k=1,j-1
     sum=sum-a(i,k)*a(k,j)
   end do
   a(i,j)=sum
   dum=vv(i)*abs(sum)
   if (abs(dum).ge.big) then
     big=abs(dum)
     imax=i
   end if
 end do
 if (j.ne.imax) then
   do k=1,n
     dum=a(imax,k)
     a(imax,k)=a(j,k)
     a(j,k)=dum
   end do
   d=d*(-1d0)
   vv(imax)=vv(j)
 end if
 indx(j)=imax
 if (abs(a(j,j)).eq.0d0) a(j,j)=tiny
 if (j.ne.n) then
   dum=1d0/a(j,j)
   do i=j+1,n
     a(i,j)=a(i,j)*dum
   end do
 end if
 end do

! END
 
 return
 end
!
!_____________________________________________________________________
!
! NAME
!    dlubksub
!
! DESCRIPTION
!    
!    backsubstitution through LU decomposed real*8 matrix
!
! SEE ALSO
!   Numerical recipes
!
! HISTORY
!
!     started 18/05/05 CJS
!
! COMMENTS
!
!
subroutine dlubksub(a,n,indx,b,matdim)
 
! Modules used

IMPLICIT NONE

!Argument list variables

integer matdim
real*8 a(matdim,matdim),b(matdim)
integer n,indx(matdim)

!Local variables

integer i,j,ip,ii
real*8 sum
 
! START

 ii=0
 do i=1,n
   ip=indx(i)
   sum=b(ip)
   b(ip)=b(i)
   if (ii.ne.0) then
     do j=ii,i-1
       sum=sum-a(i,j)*b(j)
     end do
   else if (sum.ne.0d0) then
     ii=i
   end if
   b(i)=sum
 end do
 do i=n,1,-1
   sum=b(i)
   do j=i+1,n
     sum=sum-a(i,j)*b(j)
   end do
   b(i)=sum/a(i,i)
 end do
 
! END
 
 return 
 end
!_______________________________________________________________
!
! NAME
!    dludinvert
!
! DESCRIPTION
!    
!    invert a real*8 matrix by LU decomposition
!
! SEE ALSO
!   Numerical recipes
!
! HISTORY
!
!     started 18/05/05 CJS
!
! COMMENTS
!
!
 subroutine dludinvert(a,n,y,matdim) 
! Modules used

IMPLICIT NONE

!Argument list variables

 integer n,matdim
 real*8 a(matdim,matdim),y(matdim,matdim)
 
!Local variables

 real*8 d
 real*8 as(matdim,matdim)
 real*8 col(1:matdim)
 integer indx(1:matdim),i,j
 
! START
 
! save matrix a
 do i=1,n
   do j=1,n
     as(i,j)=a(i,j)
   end do
 end do
!
 call dludcmp(a,n,indx,d,matdim)
 do j=1,n
   do i=1,n
     col(i)=0d0
   end do
   col(j)=1d0
   call dlubksub(a,n,indx,col,matdim)
   do i=1,n
     y(i,j)=col(i)
   end do
 end do
! restore a
 do i=1,n
   do j=1,n
     a(i,j)=as(i,j)
   end do
 end do
 
! END
 
 return
 end
!
! General Matrix manipulation routines
!

       subroutine  dmatrix_reset(nr,nc,a)

! set all elements of complex matrix a to zero
       
       integer nr,nc
       real*8 a(nr,nc)
       
       integer row,col
       
       do row=1,nr
         do col=1,nc
	   a(row,col)=0d0
	 end do
       end do
       
       return
       end
!
! __________________________________________________
!
!
!
! __________________________________________________
!
!
       subroutine  dadd(A,ar,ac,B,br,bc,C,matdim)

! multiply real*8 matrices
       
       integer ar,ac,br,bc,matdim
       real*8 A(matdim,matdim)
       real*8 B(matdim,matdim)
       real*8 C(matdim,matdim)
       
       real*8 sum
       
       integer row,col,i
  
       if (ar.ne.br) then
         print*,'matrix dimension error in dadd'
	 print*,'ar=',ar
	 print*,'br=',br
	 stop
       end if
       if (ac.ne.bc) then
         print*,'matrix dimension error in dadd'
	 print*,'ac=',ac
	 print*,'bc=',bc
	 stop
       end if
  
       do row=1,ar
         do col=1,ac
	   c(row,col)=a(row,col)+b(row,col)
	 end do
       end do
       
       return
       end
!
! __________________________________________________
!
!
       subroutine  dsub(A,ar,ac,B,br,bc,C,matdim)

! multiply real*8 matrices
       
       integer ar,ac,br,bc,matdim
       real*8 A(matdim,matdim)
       real*8 B(matdim,matdim)
       real*8 C(matdim,matdim)
       
       real*8 sum
       
       integer row,col,i
  
       if (ar.ne.br) then
         print*,'matrix dimension error in dadd'
	 print*,'ar=',ar
	 print*,'br=',br
	 stop
       end if
       if (ac.ne.bc) then
         print*,'matrix dimension error in dadd'
	 print*,'ac=',ac
	 print*,'bc=',bc
	 stop
       end if
  
       do row=1,ar
         do col=1,ac
	   c(row,col)=a(row,col)-b(row,col)
	 end do
       end do
       
       return
       end
!
! __________________________________________________
!
!
       subroutine dsvd(A,ar,ac,matdim,U,GAMMA,VT) 

! invert matrix using SVD implemented in LAPACK
       
       integer ar,ac,matdim
       real*8 A(matdim,matdim)
       DOUBLE PRECISION   GAMMA(matdim)
       DOUBLE PRECISION      U(matdim,matdim),VT(matdim,matdim)
       
!      LAPACK ARGUMENTS

       CHARACTER          JOBU, JOBVT
       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
       
       DOUBLE PRECISION   WORK(5*matdim)
       
! LOCAL VARIABLES
       integer row,col       

! START
       
! perform singular value decomposition on A
      
! these were 'S' before...
      JOBU='A'
      JOBVT='A'
      M=ar
      N=ac
      LDA=matdim
      LDU=matdim
      LDVT=matdim
      LWORK=5*matdim
                          
      CALL DGESVD(JOBU,JOBVT,M,N,A,LDA,GAMMA,U,LDU,VT,LDVT,   &
                      WORK,LWORK,INFO)
       if(info.ne.0) then
         if(info.lt.0) then
	   print*,'Illegal argument to dgesvd',-info
	   stop
	 else if(info.gt.0) then
	   print*,'DGESVD did not converge, info=',info
	   stop
	 end if
       end if     
                      
      return
      end
!
! __________________________________________________
!
!
       subroutine dsvd_invert(A,ar,ac,AI,matdim) 

! invert matrix using SVD implemented in LAPACK
       
       integer ar,ac,matdim
       real*8 A(matdim,matdim)
       real*8 AI(matdim,matdim)
                     
       DOUBLE PRECISION  GAMMA(matdim)
       DOUBLE PRECISION  GAMMA_I(matdim)
       real*8   V(matdim,matdim),VT(matdim,matdim)
       real*8   U(matdim,matdim),UT(matdim,matdim)
       
       real*8   ITEST(matdim,matdim)
       
       real*8 tm1(matdim,matdim),tm2(matdim,matdim)
       
       integer row,col

! START
     
! zero matrices 
     
       do row=1,matdim
         do col=1,matdim
	   U(row,col)=(0.0,0.0)
	   UT(row,col)=(0.0,0.0)
	   V(row,col)=(0.0,0.0)
	   VT(row,col)=(0.0,0.0)
	   tm1(row,col)=(0.0,0.0)
	   tm2(row,col)=(0.0,0.0)
	 end do
	 GAMMA(row)=0.0
	 GAMMA_I(row)=0.0
       end do
       
! generate unit matrix 
     
       do row=1,matdim
         do col=1,matdim
	   itest(row,col)=0.0
	 end do
	 itest(row,row)=1.0
       end do
              
! save matrix a for checks
       do row=1,matdim
         do col=1,matdim
	   tm1(row,col)=a(row,col)
	 end do
       end do

! perform singular value decomposition on A
       
       call dsvd(A,ar,ac,matdim,U,GAMMA,VT)
              
! calculate gamma inverse
       
       do row=1,ac
         if (abs(gamma(row)).ne.(0.0)) then
           gamma_i(row)=1.0/gamma(row)
         else 
           gamma_i(row)=0.0
	 end if
       end do
              
! calculate transposes of U and VT
       call dtranspose(U,ar,ac,UT,matdim)
       call dtranspose(VT,ac,ac,V,matdim)       
                
! pre multiply U transpose by gamma inverse

       do row=1,ac
         do col=1,ar
	   tm2(row,col)=UT(row,col)*gamma_i(row)
	 end do
       end do
                     
! pre multiply result by V to give ai, the inverse of a
                    	
       call dmatmul(V,ac,ac,tm2,ac,ar,AI,matdim)	 
              
! restore a
       do row=1,ar
         do col=1,ac
	   a(row,col)=tm1(row,col)
	 end do
       end do

       return

! checks...
      
!A=U.gamma.VT   

       print*,'A=U.gamma.VT'
       do row=1,ac
	 do col=1,ac
	  tm1(row,col)=VT(row,col)*gamma(row)
	end do
       end do	    
       call dmatmul(U,ar,ac,tm1,ac,ac,tm2,matdim)	
       call check_difference(A,ar,ac,tm2,matdim)   
	  
       print*,'UT.U=I'
       call dmatmul(UT,ac,ar,U,ar,ac,tm2,matdim)	
       call check_difference(itest,ac,ac,tm2,matdim)	   
      
       print*,'VT.V=I'
       call dmatmul(VT,ac,ac,V,ac,ac,tm2,matdim)	
      call check_difference(itest,ac,ac,tm2,matdim)	
	
       print*,'V.VT=I'
       call dmatmul(V,ac,ac,VT,ac,ac,tm2,matdim)	
       call check_difference(itest,ac,ac,tm2,matdim)
       	   
       print*,'A^-1.A=I' 
       call dmatmul(AI,ac,ar,A,ar,ac,tm1,matdim)
       call check_difference(tm1,ac,ac,itest,matdim)
       
       return
       end
!
! __________________________________________________
!
!
       subroutine deig(A,ar,GAMMA,matdim) 

! calculate the eigenvalues of A using LAPACK

       integer ar,matdim
       real*8 A(matdim,matdim)
       complex*16 GAMMA(matdim)


      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   VL( 1, matdim ), VR( 1, matdim ),           &
                        WI( matdim ), WORK( 4*matdim ), WR( matdim )
       
      
! LOCAL VARIABLES
       integer row    

! START
      
      JOBVL='N'
      JOBVR='N'
      LDA=matdim
      N=ar
      LDVL=1
      LDVR=1
      LWORK=4*matdim
      
      CALL DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,  &
                       LDVR, WORK, LWORK, INFO )
                          
       if(info.ne.0) then
         if(info.lt.0) then
	   print*,'Illegal argument to dgeev',-info
	   stop
	 else if(info.gt.0) then
	   print*,'dgeev did not converge, info=',info
	   stop
	 end if
       end if    
       
       do row=1,n
         gamma(row)=cmplx(wr(row),wi(row)) 
       end do               
		      
      return
      end
!
! __________________________________________________
!
!
       subroutine test_diagonalise_real_symmetric

     integer matdim
     parameter (matdim=10)
     
     real*8 A(matdim,matdim)
     real*8 P(matdim,matdim),D(matdim,matdim),PT(matdim,matdim)
     real*8 TM1(matdim,matdim),TM2(matdim,matdim)
     
     integer row,col,nw
     
! START

      nw=4
      
      do row=1,nw
        do col=1,nw
	  A(row,col)=max(row,col)
	end do
      end do
		
      write(*,*)     'A='		       
      do row=1,nw
        write(*,1000)(A(row,col),col=1,matdim)
      end do
		       
      call  diagonalise_real_symmetric(A,4,P,D,PT,matdim)
		
      write(*,*)     'P='		       
      do row=1,nw
        write(*,1000)(P(row,col),col=1,nw)
      end do
		
      write(*,*)     'D='		       
      do row=1,nw
        write(*,1000)(D(row,col),col=1,nw)
      end do
		
      write(*,*)     'PT='		       
      do row=1,nw
        write(*,1000)(PT(row,col),col=1,nw)
      end do

      call dmatmul(P,nw,nw,D,nw,nw,TM1,matdim)
      call dmatmul(TM1,nw,nw,PT,nw,nw,TM2,matdim)
		
      write(*,*)     'A='		       
      do row=1,nw
        write(*,1000)(TM2(row,col),col=1,nw)
      end do

      STOP

1000  format(5F10.4)

      return
      end
!
! __________________________________________________
!
!
       subroutine diagonalise_real_symmetric(A,ar,P,D,PT,matdim) 

       integer ar,matdim
       real*8 A(matdim,matdim)
       real*8 P(matdim,matdim),D(matdim,matdim),PT(matdim,matdim)
             
! LOCAL VARIABLES

     INTEGER   N,NSELECT
     INTEGER	      INFO, LDA, LDZ, LWORK, LIWORK, IL, IU, M
     DOUBLE PRECISION ABSTOL, VL, VU

     INTEGER	      ISUPPZ( 4*matdim )
     DOUBLE PRECISION,allocatable :: WORK( : )
     INTEGER,allocatable :: IWORK( : )
     DOUBLE PRECISION  W( matdim ), Z( matdim, matdim )

     integer row

!START  

     ABSTOL = -1.0
     
     LDA=matdim
     LDZ=matdim
     N=ar
     NSELECT=ar
     IL = 1
     IU = NSELECT

!     Query the workspace.
     
     LWORK=-1
     
     allocate( WORK(1:matdim) )
     allocate( IWORK(1:matdim) )
     
     CALL DSYEVR( 'V', 'A', 'Upper', N, A, LDA, VL, VU, IL,   &
		 IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,   &
		 LIWORK, INFO )
     
!     write(*,*)'Lwork=',WORK( 1 )
!     write(*,*)'ILwork=',IWORK( 1 )
     
     LWORK=int( WORK( 1 ) )
     LIWORK=int( IWORK( 1 ) )
		 
     deallocate( WORK )
     deallocate( IWORK )
     
     allocate( WORK(1:LWORK) )
     allocate( IWORK(1:LIWORK) )
		 
!     Solve eigenproblem.

     P(:,:)=A(:,:)

     CALL DSYEVR( 'V', 'A', 'Upper', N, P, LDA, VL, VU, IL,   &
    		  IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,   &
     		  LIWORK, INFO )
		 
     deallocate( WORK )
     deallocate( IWORK )

!     Check for convergence.

     IF( INFO.GT.0 ) THEN
	WRITE(*,*)'The algorithm failed to compute eigenvalues.'
	STOP
     END IF
     
     D(:,:)=0d0
     P(:,:)=Z(:,:)
     
     call dtranspose(P,N,N,PT,matdim)
     
     do row=1,N
       D(row,row)=W(row)
     end do

   return
   end

!
! __________________________________________________
!
!
       subroutine check_difference(m1,nr,nc,m2,matdim)
       
       integer matdim
       real*8 m1(matdim,matdim),m2(matdim,matdim)
       integer nr,nc
       integer r,c
       
       real*8 diff,maxdiff,max_element,max_out
       integer rmax,cmax
       
       rmax=0
       cmax=0
       maxdiff=0.0
       max_element=0.0
       
       do r=1,nr
         do c=1,nc
	   diff=abs(m1(r,c)-m2(r,c))
	   if (diff.gt.maxdiff) then
	     maxdiff=diff
	     rmax=r
	     cmax=c
	   end if
!	   if (abs(m1(r,c)).gt.max_element) max_element=abs(m1(r,c))
	   if (abs(m2(r,c)).gt.max_element) max_element=abs(m2(r,c))
	 end do
       end do
       
       if(max_element.ne.0.0) then
         max_out=maxdiff/max_element
       else
	 max_out=maxdiff
       end if
       
       print*,'Max difference= ',max_out,' at row ',rmax,' col',cmax
       print*,' '

       return
       end
!
! __________________________________________________
!
! 


