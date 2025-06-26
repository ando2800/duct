!subroutine pseudopois(D,wk,wk1,eii,m) ! size of wk1 was strange
subroutine pseudopois(m)
  !  use ctes
  use ctes,only:nk
  use fulldiagneu
  use numbers
  implicit none
  !implicit real*8(a-h,o-z)
  !      common/fulldiagneu/EN(my-2),PN(my-2,my-2),
  !     %     PNM1(my-2,my-2),APS(my-2,my-2)  ! APS --> PNT
  !      common /fulldineexv/wanted,iofzero
  integer ipiv(m-2)  
  real(nk),dimension(:,:),allocatable:: D,wk1
  real(nk),dimension(:),allocatable:: eii, wk, x  
  integer n,m,mm2,lwk,i,ii,j,jj,info
  real(nk) dummy,v1,v2,vv1,vv2,eig,eigmin
  real(nk) tch0,tch2 ! see functions below
  !      dimension D(m-2,*),wk1(m-2,*),eii(m-2)
  !      dimension wk(4*(my-2)),x(my)

  n=m-1
  mm2=m-2
  lwk=64*mm2
  
  allocate(D(m-2,m-2),wk1(m-2,m-2))
  D=0.d0; wk1=0.d0
  allocate(wk(lwk),eii(m-2),x(m))
  wk=0.d0; eii=0.d0; x=0.d0

  !pi=4.d0*datan(1.d0) ! see in numbers
  ! Chebyshev collocation nodes
  
  do j=1,m
     x(j)=-dcos(pi*dfloat(j-1)/dfloat(n))
  enddo
  
  do i=2,m-1
     ii=i-1
     do j=1,n-1
        jj=j-1
        v1=tch0(jj,x(i))
        v2=tch2(jj,x(i))
        vv1=tch0(jj+2,x(i))
        vv2=tch2(jj+2,x(i))
        PNT(ii,j)=v1-jj**2*vv1/(jj+2.d0)**2
        wk1(ii,j)=v2-jj**2*vv2/(jj+2.d0)**2
     enddo
  enddo
  call dgetrf(mm2,mm2,PNT,mm2,ipiv,info)
  !call dgetri(mm2,PNT,mm2,ipiv,wk,-1,info)
  !write(*,*) 'pseudopois: check lwk=', wk(1), lwk
  call dgetri(mm2,PNT,mm2,ipiv,wk,lwk,info)
  call dgemm('n','n',mm2,mm2,mm2,1.d0,wk1,mm2,PNT,mm2,0.d0,D,mm2)
  
  ! Solve the eigenvalue problem Neumann case
  call dgeev('n','v',mm2,D,mm2,EN,eii,dummy,1, & 
       &             PN,mm2,wk,lwk,info)

  call dcopy(mm2*mm2,PN,1,PNM1,1)
  call dgetrf(mm2,mm2,PNM1,mm2,ipiv,info)
  !call dgetri(mm2,PNM1,mm2,ipiv,wk,-1,info)
  !write(*,*) 'pseudopois: check lwk=', wk(1), lwk
  call dgetri(mm2,PNM1,mm2,ipiv,wk,lwk,info)

  ! find out the location of the zero eigenvalue

  eigmin=1.d6
  
  do i=1,mm2
     eig=dmin1(eigmin,dabs(EN(i)))
     if(eig.ne.eigmin) iofzero=i
     if(eig.ne.eigmin) eigmin=eig
  enddo
  
  wanted=0.d0

  deallocate(wk,eii,x)
  deallocate(D,wk1)

  return
 
end subroutine pseudopois

real(8) function tch0(n,x)
  implicit none 
  integer n
  real(8) x,argument
  argument=dfloat(n)*dacos(x)
  tch0=dcos(argument)
  return 
end function tch0
      
real(8) function tch2(n,x)
  
  implicit none 
  integer n
  real(8) argument,x,tch0,den
  argument=dfloat(n)*dacos(x)
  tch0=dcos(argument)
  den=1.d0-x**2
  tch2=x*dsin(argument)-dfloat(n)*tch0*dsqrt(den)
  tch2=tch2*dfloat(n)/den**(3.d0/2.d0)
  return 
end function tch2
