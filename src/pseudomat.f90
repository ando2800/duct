#ifndef RECTANGULAR
subroutine pseudomat(bound,m)
  
  use ctes,only:nk,my,mz,ngpu
  use numbers,only:pi
  use chebystuff
  use fulldiagneu
  use fulldiagdir
#ifdef CUBLAS
  use gpudevice
#endif

  implicit none 
  !
  ! preprocessing stuff
  !
  ! bound either Neumann: n, or Dirichlet:d.
  !
  !
  !      common/fulldiagneu/EN(my-2),PN(my-2,my-2),PNM1(my-2,my-2),
  !     $     PNT(my-2,my-2),PNM1T(my-2,my-2)
  !      common /fulldineexv/AM1(4,4),wanted,iofzero
  !
  !      common/fulldiagdir/ED(my-2),PD(my-2,my-2),PDM1(my-2,my-2),
  !     $     PDT(my-2,my-2),PDM1T(my-2,my-2)
  !      common /chebystuff/ D(my,my),D2bx(my,2),D2by(my,2),x(my)
  
  real(nk),dimension(:,:),allocatable:: D2B,D2
  real(nk),dimension(:),allocatable:: dd, eii, wk, xcol
  real(nk) coef,dummy,Det,eig,eigmin,ter1,ter3
  integer ipiv(my-2),ipiva(4),info,istat
  integer m,n,mm2,lwk,i,j,k,l,icoi,ii,jj
  character bound

  n=m-1
  mm2=m-2
  lwk=64*(m-2)

  !pi=4.d0*datan(1.d0)

  ! allocate common arrays
  if (bound.eq.'n') then
     allocate(EN(my-2))
     EN=0.d0
     allocate(EN2(my-2+1,my-2+1))
     EN2=0.d0
     allocate(PN(my-2,my-2),PNM1(my-2,my-2),PNT(my-2+1,my-2+1),PNM1T(my-2+1,my-2+1))
     PN=0.d0; PNM1=0.d0; PNT=0.d0; PNM1T=0.d0
     allocate(tmprt(my-1,mz-1))
     tmprt=0.d0

  elseif (bound.eq.'d') then
     allocate(ED(my-2))
     ED=0.d0
     allocate(ED2(my-2+1,my-2+1))
     ED2=0.d0
     allocate(PD(my-2,my-2),PDM1(my-2,my-2),PDT(my-2+1,my-2+1),PDM1T(my-2+1,my-2+1))
     PD=0.d0; PDM1=0.d0; PDT=0.d0; PDM1T=0.d0

     allocate(DC(my,my),D2bx(my,2),D2by(my,2))
     DC=0.d0; D2bx=0.d0; D2by=0.d0
     !allocate(tmprt(my-1,mz-1))
     !tmprt=0.d0
  end if
  !
  ! allocate tmp arrays --> deallocate at the end of this subroutine
  allocate(D2B(m-2,m-2),D2(m,m))
  D2B=0.d0; D2=0.d0
  allocate(dd(m),eii(m-2), wk(lwk), xcol(m) )
  dd=0.d0; eii=0.d0; wk=0.d0; xcol=0.d0

  dd(1)=2.d0
  do i=2,m-1
     dd(i)=1.d0
  enddo
  dd(m)=2.d0
  
  call dcopy(m,0.d0,0,xcol,1)
  call dcopy(m*m,0.d0,0,DC,1)

  ! Chebyshev collocation nodes

  do j=1,m
     xcol(j)=-dcos(pi*dfloat(j-1)/dfloat(n))
  enddo
      
  !C Assemble Chebyshev collocation derivative matrix
  !C  see CHQZ p. 69
  
  DC(1,1)=-(2.d0*dfloat(n)**2+1.d0)/6.d0
  DC(m,m)=-DC(1,1)

  do i=2,n
     DC(i,i)=-0.5d0*xcol(i)/(1.d0-xcol(i)**2)
  enddo

  do l=1,m
     do j=1,l-1
        icoi=(l-1)+(j-1)
        coef=dd(l)/dd(j)
        DC(l,j)=coef*(-1.d0)**icoi/(xcol(l)-xcol(j))
     enddo
     do j=l+1,m
        coef=dd(l)/dd(j);
        DC(l,j)=coef*(-1.d0)**(l+j)/(xcol(l)-xcol(j))
     enddo
  enddo

  ! Compute the second derivative matrix D2=D*D
  call dgemm('n','n',m,m,m,1.d0,DC,m,DC,m,0.d0,D2,m)
  
  call dcopy(mm2*mm2,0.d0,0,D2B,1)

  if (bound.eq.'n') then
     
     Det=DC(1,1)*DC(m,m)-DC(1,m)*DC(m,1)
     
     ii=0
     do i=2,n
        ii=ii+1
        jj=0
        do j=2,n
           jj=jj+1
           ter1=D2(i,1)*(DC(1,m)*DC(m,j)-DC(m,m)*DC(1,j))/Det
           ter3=D2(i,m)*(DC(m,1)*DC(1,j)-DC(1,1)*DC(m,j))/Det
           D2B(ii,jj)=ter1+D2(i,j)+ter3
        enddo
     enddo
     
     ! Solve the eigenvalue problem Neumann case
     call dgeev('n','v',mm2,D2B,mm2,EN,eii,dummy,1, &
          &     PN,mm2,wk,lwk,info)
      
     do j=1,my-2
        do i=1,my-2
           EN2(i,j)=EN(i)+EN(j)
        enddo
     enddo  
    
     call dcopy(mm2*mm2,PN,1,PNM1,1)
     call dgetrf(mm2,mm2,PNM1,mm2,ipiv,info)
     call dgetri(mm2,PNM1,mm2,ipiv,wk,lwk,info)
     
     ! find out the location of the zero eigenvalue
     
     eigmin=1.d6
     
     do i=1,mm2
        eig=dmin1(eigmin,dabs(EN(i)))
        if(eig.ne.eigmin) iofzero=i
        if(eig.ne.eigmin) eigmin=eig
     enddo
     ! fixed wanted = 0 in iofzero+1,iofzero+1
     
     !write(*,*) 'pseudomat: iofzeo', iofzero
     wanted=0.d0
     
     AM1(1,1)=2.d0*DC(1,1)
     AM1(1,2)=DC(1,m)
     AM1(1,3)=0.d0
     AM1(1,4)=DC(1,m)
     
     AM1(2,1)=-DC(m,1)
     AM1(2,2)=DC(1,1)-DC(m,m)
     AM1(2,3)=DC(1,m)
     AM1(2,4)=0.d0
     
     AM1(3,1)=0.d0
     AM1(3,2)=DC(m,1)
     AM1(3,3)=2.d0*DC(m,m)
     AM1(3,4)=DC(m,1)
     
     AM1(4,1)=-DC(m,1)
     AM1(4,2)=0.d0
     AM1(4,3)=DC(1,m)
     AM1(4,4)=DC(1,1)-DC(m,m)
     call dgetrf(4,4,AM1,4,ipiva,info)
     call dgetri(4,AM1,4,ipiva,wk,lwk,info)
     
     ! compute and store transposed matrix Neumann case
     do i=1,my-2
        do j=1,my-2
           PNT(i,j) = PN(j,i)
           PNM1T(i,j) = PNM1(j,i)
        enddo
     enddo
     
     if (ngpu.gt.0) then
#ifdef CUBLAS
        !istat = cusparse_create(handle) ! for transpose library
        !write(*,*) 'pseudomat allocate gpudevice'
        
        istat= cublas_alloc((mm2+1)*(mm2+1), nk, dp_PNT)
        if (istat.ne.0) write(*,*) 'cublas alloc. error: PNT',istat
        istat= cublas_set_matrix( mm2+1,mm2+1,nk,PNT,mm2+1,dp_PNT,mm2+1)
        if (istat.ne.0) write(*,*) 'cublas setmat error: PNT',istat
        
        istat= cublas_alloc((mm2+1)*(mm2+1), nk, dp_PNM1T)
        if (istat.ne.0) write(*,*) 'cublas alloc. error: PNM1T',istat        
        istat= cublas_set_matrix( mm2+1,mm2+1,nk,PNM1T,mm2+1,dp_PNM1T,mm2+1)
        if (istat.ne.0) write(*,*) 'cublas setmat error: PNM1T',istat
        
        istat= cublas_alloc((mm2+1)*(mm2+1), nk, dp_EN2)
        if (istat.ne.0) write(*,*) 'cublas alloc. error: EN2',istat
        istat= cublas_set_matrix( mm2+1,mm2+1,nk,EN2,mm2+1,dp_EN2,mm2+1)
        if (istat.ne.0) write(*,*) 'cublas setmat error: EN2',istat  

        !! allocate tmp arrays 
        !istat= cublas_alloc((mm2+1)*(mm2+1), nk, devPtrF)
        !if (istat.ne.0) write(*,*) 'cublas alloc. error: F',istat  
        !istat= cublas_set_matrix(mm2+1,mm2+1,nk,0.d0,0,devPtrF,1) ! zero set
        !
        !istat= cublas_alloc((mm2+1)*(mm2+1), nk, devPtrFT)
        !if (istat.ne.0) write(*,*) 'cublas alloc. error: FT',istat   
        !istat= cublas_set_matrix(mm2+1,mm2+1,nk,0.d0,0,devPtrFT,1) ! zero set
        !
        !istat= cublas_alloc((mm2+1)*(mm2+1), nk, dp_tmp2d)
        !if (istat.ne.0) write(*,*) 'cublas alloc. error: tmp2d',istat   
        !istat= cublas_set_matrix(mm2+1,mm2+1,nk,0.d0,0,dp_tmp2d,1) ! zero set
#endif /* CUBLAS */
     end if

  else if (bound.eq.'d') then

     do i=1,m
        D2bx(i,1)=D2(i,1)
        D2bx(i,2)=D2(i,m)
        D2by(i,1)=D2(i,1)
        D2by(i,2)=D2(i,m)
     enddo
     
     ii=0
     do i=2,n
        ii=ii+1
        jj=0
        do j=2,n
           jj=jj+1
           D2B(ii,jj)=D2(i,j)
        enddo
     enddo
     
     ! Solve the eigenvalue problem Dirichlet case
     wk=0.d0
     call dgeev('n','v',mm2,D2B,mm2,ED,eii,dummy,1, &
          &     PD,mm2,wk,lwk,info)
     
     do j=1,my-2
        do i=1,my-2
           ED2(i,j)=ED(i)+ED(j)
        enddo
     enddo    

     call dcopy(mm2*mm2,PD,1,PDM1,1)
     call dgetrf(mm2,mm2,PDM1,mm2,ipiv,info)
     call dgetri(mm2,PDM1,mm2,ipiv,wk,lwk,info)
     
     ! compute and store transposed matrix Dirichlet case     
     do i=1,my-2
        do j=1,my-2
           PDT(i,j) = PD(j,i)
           PDM1T(i,j) = PDM1(j,i)
        enddo
     enddo
     
     if (ngpu.gt.0) then
#ifdef CUBLAS
        istat = cusparse_create(handle) ! for transpose library
        !write(*,*) 'pseudomat allocate gpudevice'

        istat= cublas_alloc((mm2+1)*(mm2+1), nk, dp_PDT)
        if (istat.ne.0) write(*,*) 'cublas alloc. error: PDT',istat
        istat= cublas_set_matrix( mm2+1,mm2+1,nk,PDT,mm2+1,dp_PDT,mm2+1)
        if (istat.ne.0) write(*,*) 'cublas setmat error: PDT',istat

        istat= cublas_alloc((mm2+1)*(mm2+1), nk, dp_PDM1T)
        if (istat.ne.0) write(*,*) 'cublas alloc. error: PDM1T',istat        
        istat= cublas_set_matrix( mm2+1,mm2+1,nk,PDM1T,mm2+1,dp_PDM1T,mm2+1)
        if (istat.ne.0) write(*,*) 'cublas setmat error: PDM1T',istat

        istat= cublas_alloc((mm2+1)*(mm2+1), nk, dp_ED2)
        if (istat.ne.0) write(*,*) 'cublas alloc. error: ED2',istat
        istat= cublas_set_matrix( mm2+1,mm2+1,nk,ED2,mm2+1,dp_ED2,mm2+1)
        if (istat.ne.0) write(*,*) 'cublas setmat error: ED2',istat        

        ! allocate tmp arrays 
        istat= cublas_alloc((mm2+1)*(mm2+1), nk, devPtrF)
        if (istat.ne.0) write(*,*) 'cublas alloc. error: F',istat  
        istat= cublas_set_matrix(mm2+1,mm2+1,nk,0.d0,0,devPtrF,1) ! zero set

        istat= cublas_alloc((mm2+1)*(mm2+1), nk, devPtrFT)
        if (istat.ne.0) write(*,*) 'cublas alloc. error: FT',istat   
        istat= cublas_set_matrix(mm2+1,mm2+1,nk,0.d0,0,devPtrFT,1) ! zero set

        istat= cublas_alloc((mm2+1)*(mm2+1), nk, dp_tmp2d)
        if (istat.ne.0) write(*,*) 'cublas alloc. error: tmp2d',istat   
        istat= cublas_set_matrix(mm2+1,mm2+1,nk,0.d0,0,dp_tmp2d,1) ! zero set
#endif /* CUBLAS */
     end if

  else
     
     write(6,*) 'in pseudomat error'     
     stop
     
  endif
  
  ! deallocate tmp arrays only for this routine
  deallocate(dd,eii,wk,xcol)
  deallocate(D2B,D2)
   
  return 
end subroutine pseudomat
#else /* RECTANGULAR */
subroutine pseudomat_dummy()
end subroutine pseudomat_dummy
#endif /* not RECTANGULAR */
