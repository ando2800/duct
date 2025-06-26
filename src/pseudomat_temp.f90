!subroutine pseudomaty_temp(D2T,D2,eii,dd,bound1,bound2,m)
subroutine pseudomaty_temp(bound1,bound2,m)
  use ctes,only:nk,my
  use numbers
  use chebystuff_temp
  use fulldiagtemp
  implicit none 
  !     
  ! preprocessing stuff
  ! bound either Neumann: n, or Dirichlet:d.
  !     
  integer m,n,mm2,lwk
  real(nk) y(m)
  real(nk),dimension(:,:),allocatable:: D2T,D2
  real(nk),dimension(:),allocatable:: eii,dd
  real(nk),dimension(:),allocatable:: wk
  integer ipiv(my-2)
  character*1 bound1,bound2
  integer i,j,ii,jj,icoi,l,info
  real(nk) coef,Det,ter1,ter3,dummy
  
  n=m-1
  mm2=m-2

  ! allocate for chebystuff_temp
  allocate(Dy(my,my),D2by(my,2))
  Dy=0.d0; D2by=0.d0
  ! allocate for fulldiagtemp
  allocate(ETY(my-2),PTY(my-2,my-2),PTYM1(my-2,my-2))
  ETY=0.d0; PTY=0.d0; PTYM1=0.d0

  ! allocate tmp arrays
  allocate(D2T(m-2,m-2),D2(m,m))
  D2T=0.d0; D2=0.d0
  allocate(eii(m-2),dd(m))
  eii=0.d0; dd=0.d0;

  lwk=64*(my-2)
  allocate(wk(lwk))
  wk=0.d0
  !     
  dd(1)=2.d0
  do i=2,m-1
     dd(i)=1.d0
  enddo
  dd(m)=2.d0
  
  call dcopy(m,0.d0,0,y,1)
  call dcopy(m*m,0.d0,0,Dy,1)
  !     
  ! Chebyshev collocation nodes
  !     
  do j=1,m
     y(j)=-dcos(pi*dfloat(j-1)/dfloat(n)) 
  enddo
      
  !     Assemble Chebyshev collocation derivative matrix
  !     see CHQZ p. 69      
  Dy(1,1)=-(2.d0*dfloat(n)**2+1.d0)/6.d0
  Dy(m,m)=-Dy(1,1)      
  do i=2,n
     Dy(i,i)=-0.5d0*y(i)/(1.d0-y(i)**2)
  enddo
  
  do l=1,m
     do j=1,l-1
        icoi=(l-1)+(j-1)
        coef=dd(l)/dd(j)
        Dy(l,j)=coef*(-1.d0)**icoi/(y(l)-y(j))
     enddo
     do j=l+1,m
        coef=dd(l)/dd(j);
        Dy(l,j)=coef*(-1.d0)**(l+j)/(y(l)-y(j))
     enddo
  enddo
      
  call dscal(m*m,1.d0/sy,Dy,1)
      
  ! Compute the second derivative matrix D2=D*D
  call dgemm('n','n',m,m,m,1.d0,Dy,m,Dy,m,0.d0,D2,m)
      
  call dcopy(mm2*mm2,0.d0,0,D2T,1)
  !
  do i=1,m
     D2by(i,1)=D2(i,1)
     D2by(i,2)=D2(i,m)
  enddo
  !
  if ((bound1.eq.'n').and.(bound2.eq.'n')) then
     
     Det=Dy(1,1)*Dy(m,m)-Dy(1,m)*Dy(m,1)
     
     ii=0
     do i=2,n
        ii=ii+1
        jj=0
        do j=2,n
           jj=jj+1
           ter1=D2(i,1)*(Dy(1,m)*Dy(m,j)-Dy(m,m)*Dy(1,j))/Det
           ter3=D2(i,m)*(Dy(m,1)*Dy(1,j)-Dy(1,1)*Dy(m,j))/Det
           D2T(ii,jj)=ter1+D2(i,j)+ter3
        enddo
     enddo
         
     ! find out the location of the zero eigenvalue: removed
                  
  else if ((bound1.eq.'d').and.(bound2.eq.'d')) then
         
     ii=0
     do i=2,n
        ii=ii+1
        jj=0
        do j=2,n
           jj=jj+1
           D2T(ii,jj)=D2(i,j)
        enddo
     enddo
     
  else if((bound1.eq.'d').and.(bound2.eq.'n'))then
     do i=2,m-1
        ii=i-1;
        do j=2,m-1
           jj=j-1;
           D2T(ii,jj)=D2(i,j)-D2(i,m)*Dy(m,j)/Dy(m,m);
        end do
     end do
     
  else if((bound1.eq.'n').and.(bound2.eq.'d'))then
     do i=2,m-1
        ii=i-1;
        do j=2,m-1
           jj=j-1;
           D2T(ii,jj)=D2(i,j)-D2(i,1)*Dy(1,j)/Dy(1,1);
        end do
     end do
     
  endif

  !     Solve the eigenvalue problem
  call dgeev('n','v',mm2,D2T,mm2,ETY,eii,dummy,1,PTY,mm2,wk,lwk,info)
  call dcopy(mm2*mm2,PTY,1,PTYM1,1)
  call dgetrf(mm2,mm2,PTYM1,mm2,ipiv,info)
  call dgetri(mm2,PTYM1,mm2,ipiv,wk,lwk,info)
  
  deallocate(wk)
  deallocate(eii,dd)
  deallocate(D2T,D2)
       
  return 
end subroutine pseudomaty_temp
!----|----------------------------------------------------------------|
subroutine pseudomatz_temp(bound1,bound2,m)
  use ctes,only:nk,mz
  use numbers
  use chebystuff_temp
  use fulldiagtemp 
  implicit none 
  !     
  ! preprocessing stuff
  !     
  ! bound either Neumann: n, or Dirichlet:d.
  !     
  integer m,n,mm2,lwk
  real(nk) z(m)
  real(nk),dimension(:,:),allocatable:: D2T,D2
  real(nk),dimension(:),allocatable:: eii,dd
  real(nk),dimension(:),allocatable:: wk
  integer ipiv(mz-2),ipiva(4)
  character*1 bound1,bound2
  integer i,j,ii,jj,icoi,l,info
  real(nk) coef,Det,ter1,ter3,dummy
    
  n=m-1
  mm2=m-2

  ! allocate for chebystuff_temp
  allocate(Dz(mz,mz),D2bz(mz,2))
  Dz=0.d0; D2bz=0.d0
  ! allocate for fulldiagtemp
  allocate(ETZ(mz-2),PTZ(mz-2,mz-2),PTZM1(mz-2,mz-2))
  ETZ=0.d0; PTZ=0.d0; PTZM1=0.d0

  ! allocate tmp arrays
  allocate(D2T(m-2,m-2),D2(m,m))
  D2T=0.d0; D2=0.d0
  allocate(eii(m-2),dd(m))
  eii=0.d0; dd=0.d0;  
  
  lwk=64*(mz-2)  
  allocate(wk(lwk))
  wk=0.d0
      
  dd(1)=2.d0
  do i=2,m-1
     dd(i)=1.d0
  enddo
  dd(m)=2.d0
      
  call dcopy(m,0.d0,0,z,1)
  call dcopy(m*m,0.d0,0,Dz,1)
      
  ! Chebyshev collocation nodes
      
  do j=1,m
     z(j)=-dcos(pi*dfloat(j-1)/dfloat(n))
  enddo
      
  ! Assemble Chebyshev collocation derivative matrix
  ! see CHQZ p. 69
      
  Dz(1,1)=-(2.d0*dfloat(n)**2+1.d0)/6.d0
  Dz(m,m)=-Dz(1,1)
  
  do i=2,n
     Dz(i,i)=-0.5d0*z(i)/(1.d0-z(i)**2)
  enddo
  
  do l=1,m
     do j=1,l-1
        icoi=(l-1)+(j-1)
        coef=dd(l)/dd(j)
        Dz(l,j)=coef*(-1.d0)**icoi/(z(l)-z(j))
     enddo
     do j=l+1,m
        coef=dd(l)/dd(j);
        Dz(l,j)=coef*(-1.d0)**(l+j)/(z(l)-z(j))
     enddo
  enddo
  
  call dscal(m*m,1.d0/sz,Dz,1)
      
  ! Compute the second derivative matrix D2=D*D
  call dgemm('n','n',m,m,m,1.d0,Dz,m,Dz,m,0.d0,D2,m)
  
  call dcopy(mm2*mm2,0.d0,0,D2T,1)
  
  do i=1,mz
     D2bz(i,1)=D2(i,1)
     D2bz(i,2)=D2(i,mz)
  enddo
  
  if ((bound1.eq.'n').and.(bound2.eq.'n')) then
     
     Det=Dz(1,1)*Dz(m,m)-Dz(1,m)*Dz(m,1)
     
     ii=0
     do i=2,n
        ii=ii+1
        jj=0
        do j=2,n
           jj=jj+1
           ter1=D2(i,1)*(Dz(1,m)*Dz(m,j)-Dz(m,m)*Dz(1,j))/Det
           ter3=D2(i,m)*(Dz(m,1)*Dz(1,j)-Dz(1,1)*Dz(m,j))/Det
           D2T(ii,jj)=ter1+D2(i,j)+ter3
        enddo
     enddo
         
     ! find out the location of the zero eigenvalue
         
     !c         eigmin=1.d6
     !c         
     !c         do i=1,mm2
     !c            eig=dmin1(eigmin,dabs(ETZ(i)))
     !c            if(eig.ne.eigmin) jofzero=i
     !c            if(eig.ne.eigmin) eigmin=eig
     !c         enddo
     !c         
     !c         wanted=0.d0
     !         
     !CCC   stuff for the corners: removed
     
  else if ((bound1.eq.'d').and.(bound2.eq.'d')) then
     !c         do i=1,mz
     !c            D2bz(i,1)=D2(i,1)
     !c            D2bz(i,2)=D2(i,mz)
     !c         enddo
         
     ii=0
     do i=2,n
        ii=ii+1
        jj=0
        do j=2,n
           jj=jj+1
           D2T(ii,jj)=D2(i,j)
        enddo
     enddo
     !
  else if((bound1.eq.'d').and.(bound2.eq.'n'))then
     do i=2,m-1
        ii=i-1;
        do j=2,m-1
           jj=j-1;
           D2T(ii,jj)=D2(i,j)-D2(i,m)*Dz(m,j)/Dz(m,m);
        end do
     end do
     !
  else if((bound1.eq.'n').and.(bound2.eq.'d'))then
     do i=2,m-1
        ii=i-1;
        do j=2,m-1
           jj=j-1;
           D2T(ii,jj)=D2(i,j)-D2(i,1)*Dz(1,j)/Dz(1,1);
        end do
     end do
     !
  endif
  !     Solve the eigenvalue problem
  
  call dgeev('n','v',mm2,D2T,mm2,ETZ,eii,dummy,1,PTZ,mm2,wk,lwk,info)
  call dcopy(mm2*mm2,PTZ,1,PTZM1,1)
  call dgetrf(mm2,mm2,PTZM1,mm2,ipiv,info)
  call dgetri(mm2,PTZM1,mm2,ipiv,wk,lwk,info)
       
  deallocate(wk)
  deallocate(eii,dd)
  deallocate(D2T,D2)

  !     write(6,*) 'set into:',iofzero,jofzero
  return 
end subroutine pseudomatz_temp
! -------------------------------------------------------------------- |
