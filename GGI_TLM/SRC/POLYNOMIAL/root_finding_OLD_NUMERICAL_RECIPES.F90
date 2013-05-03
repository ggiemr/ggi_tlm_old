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
!    find_roots
!
! DESCRIPTION
!       find all the roots of a polynomial using numerical recipes routine
!
!
! SEE ALSO
!   Numerical Recipes
!
! HISTORY
!
!     started 17/05/05 CJS
!
! COMMENTS
! 
! Altered to make arrays dynamic.
! Assume the array roots is already allocated
! 

subroutine findroots(a,roots,roots_order)

! Modules used

USE polynomial_types

IMPLICIT NONE

!Argument list variables

type(polynomial), intent(in) :: a
integer     :: roots_order
complex*16  :: roots(1:roots_order)

! Local variables

integer m,np,max_order

complex*16 x,b,c
real*8 eps,test1,test2,test3
integer its,k,jj,maxm

complex*16, allocatable      :: ad(:)

! START

! set constants

maxm=100
eps=2d-10

max_order=a%order

allocate (ad(0:max_order))

m=max_order
np=max_order
do k=0,m
  ad(k)=a%coeff(k)	  
end do

do k=m,1,-1
  x=(0d0,0d0)
  
  call laguer(ad,k,x,its,np)
  
  test1=dabs(dreal((0d0,-1d0)*x))
  test2=dabs(dreal(eps*x))
  test3=dreal(x)
  if (test1.le.test2) x=dcmplx(test3)
  roots(k)=x
  
  b=ad(k)
  do jj=k-1,0,-1
    c=ad(jj)
    ad(jj)=b
    b=x*b+c
  end do
end do  

deallocate (ad)

! END
 
return
end

!
! ______________________________________________________________
!
!
! NAME
!    Laguer
!
! DESCRIPTION
!       find roots of a polynomial using Laguer's mathod
!
!
! SEE ALSO
!   Numerical Recipes
!
! HISTORY
!
!     started 17/05/05 CJS
!
! COMMENTS
!


subroutine laguer(a,m,x,its,np)

integer np
complex*16 a(0:np)
complex*16 x
integer m,its

real*8 epss
integer mr,mt,maxit

integer iter,j
real*8 abx,abp,abm,err,fmax
complex*16 dx,x1,b,d,f,g,h,sq,gp,gm,g2

! START

epss=1d-15
mr=18
mt=20
maxit=mr*mt

!print*,'Called Laguer'
!print*,'Order=',np
!print*,a(0:np)

do iter=1,maxit
!print*,iter,x
  its=iter
  b=a(m)
  err=abs(b)
  d=(0d0,0d0)
  f=(0d0,0d0)
  abx=abs(x)
  do j=m-1,0,-1
    f=x*f+d
    d=x*d+b
    b=x*b+a(j)
    err=abs(b)+abx*err
  end do
  err=err*epss
  
  if (abs(b).le.err) then
!    print*,'Returning with x='
!    print*,x
    return
  end if
  
  g=d/b
  g2=g*g
  h=g2-2d0*f/b
  sq=sqrt(dble(m-1)*((dble(m)*h)-g2))
  gp=g+sq
  gm=g-sq
  abp=abs(gp)
  abm=abs(gm)
  if (abp.lt.abm) gp=gm
  fmax=abp
  if (abm.gt.abp) fmax=abm
  if (fmax.gt.0d0) then
    dx=cmplx(m)/gp
  else
    dx=exp(log(1d0+abx))                                    &
       *( cmplx(cos(dble(iter)),sin(dble(iter))) )
  end if
  x1=x-dx
  if ((dble(x).eq.dble(x1)).and.                            &
     (dble((0d0,-1d0)*x).eq.dble((0d0,-1d0)*x1)))then
!      print*,'Returning with x='
!      print*,x
      return
  end if
  x=x1
end do

print*,'too many iterations in laguer'
stop

! END
 
return
end

!
! ______________________________________________________________
!
!
! NAME
!    root_sort
!
! DESCRIPTION
!   
!    Sort real and complex roots
!
! SEE ALSO
!  
!
! HISTORY
!
!     started 17/05/05 CJS
!
! COMMENTS
!     needs some tidying up...

 subroutine rootsort(ordert,roots,rroots,             &
                           croots,nreal,ncomplex,maxordert)

! Modules used

IMPLICIT NONE

!Argument list variables

 integer maxordert
 complex*16 roots(maxordert)
 complex*16 rroots(maxordert),croots(maxordert)
 integer ordert,nreal,ncomplex
       
!Local variables

 complex*16 swap,sum
 real*8 tangle,tanglesum
 integer i,roottype(maxordert),fpair
 integer k
 
! START

 do i=1,ordert
   if (abs(roots(i)).eq.0d0) then 
     roottype(i)=-1
   else if (abs(dble((0d0,-1d0)*roots(i))/abs(roots(i)))    &
					  .gt.1d-8) then
!  complex root found
     roottype(i)=1
   else
! real root found
     roottype(i)=-1
   end if  

 end do
!
!  reorder roots, real roots first then complex roots in pairs
!  
 nreal=0
 ncomplex=0
 do i=1,ordert
   if (roottype(i).eq.-1) then
! real root found
     nreal=nreal+1
     rroots(nreal)=roots(i)
     roottype(i)=0
   else if (roottype(i).eq.1) then
! complex root found
     ncomplex=ncomplex+1
     croots(ncomplex)=roots(i)
     roottype(i)=0
   end if
 end do 
!
 ncomplex=ncomplex/2

! pair up complex roots

 do i=1,ncomplex
   fpair=(i-1)*2+1
   if (dble(croots(fpair)).ne.0d0) then
     tangle=abs(dble((0d0,-1d0)*croots(fpair)/        &
		   dble(croots(fpair))))
   else
     tangle=1d10
   end if
   do k=fpair+1,ncomplex*2
     sum=croots(fpair)+croots(k)
     if (abs(sum).ne.0d0) then
       tanglesum=abs(dble((0d0,-1d0)*sum/	      &
			     abs(sum)))
     else
       tanglesum=0d0
     end if
     if (tanglesum/tangle.lt.1d-5) then
! conjugate root found so shift next to first root of the pair
       swap=croots(fpair+1)
       croots(fpair+1)=croots(k)
       croots(k)=swap
     end if
   end do

 end do

 do i=1,nreal
   roots(i)=rroots(i)
 end do
 do i=1,ncomplex*2
   roots(nreal+i)=croots(i)
 end do

! END

 end
