#define TRANS_GPU
!#undef TRANS_GPU
subroutine leg_cheb_mat(A,B,m)
  use ctes
  use legchebshen
#ifdef CUBLAS
  use gpudevice
#endif

  implicit none
  !  /********************************************************/
  !  /*                                                      */
  !  /* set LEGENDRE-CHEBYSHEV transformation matrix         */
  !  /* to common block /legcheb/                            */
  !  /*                                                      */
  !  /* A(m,m): full forward matrix  cheb --> leg            */
  !  /* B(m,m): full backward matrix  leg --> cheb           */
  !  /*         (A,B is upper triangular and sparse matrix)  */
  !  /*                                                      */
  !  /* /legcheb/ AE(mye,mye), even-even elements of A       */
  !  /*           AO(myo,myo), odd-odd elements of A         */
  !  /*           BE(mye,mye), even-even elements of B       */
  !  /*           BO(myo,myo), odd-odd elements of B         */
  !  /*         (AE,AO,BE,BO are upper triangular matrix )   */
  !  /*                                                      */
  !  /*  !only for square (my=mz) case                       */
  !  /********************************************************/
  !
  !  AE,AO,BE,BO are renamed --> LCAE, LCAO, LCBE, LCBO
  !
  real(nk) A(my,my),B(my,my)

  real(nk), dimension(:), allocatable:: wk 
  integer ipiv(my)
  integer lwk,info
  integer n,m,nevn,nodd,i0,i1
  integer k,j,jj,i,ii,istat
  real(nk) c1, c2

  n=m-1
  nevn=n/2+1 ! odd value (mye+1)
  nodd=n/2   ! even  (mye = myo+1)
  lwk=my*my

  allocate(wk(lwk))
  wk=0.d0

  if (ngpu.gt.0) then
     allocate(LCA(my+1,mz+1),LCB(my+1,mz+1))
     LCA=0.d0; LCB=0.d0;
  end if
  allocate(LCAE(mye+1,mye+1),LCAO(mye+1,mye+1), & 
       &   LCBE(mye+1,mye+1),LCBO(mye+1,mye+1))
  LCAE=0.d0; LCAO=0.d0; LCBE=0.d0; LCBO=0.d0
  
  !  /* set cheb -> leg matrix, A */      
  call dcopy((my+1)*(mz*1),    0.d0,0,A, 1)
  call dcopy(nevn*nevn,0.d0,0,LCAE,1)
  call dcopy((nodd+1)*(nodd+1),0.d0,0,LCAO,1)
  ! /* fill diagonal */
  i0=1
  i1=2 
  A(i0,i0)=1.d0
  A(i1,i1)=1.d0
  do k=2,n
     A(k+1,k+1)=2.d0*dfloat(k)/dfloat(2*k-1)*A(k,k)
  enddo

  ! /* fill first row */
  do k=2,n,2
     A(1,k+1)=-1.d0/dfloat(k*k-1)
  enddo
  ! /* fill in the rest */
  do j=2,n-1
     jj=j+1
     do i=1,j
        ii=i+1 
        c1=dfloat(2*i+2)/dfloat(2*i+3)
        c2=dfloat(2*i)/dfloat(2*i-1)
        A(ii,jj+1)=c1*A(ii+1,jj)+c2*A(ii-1,jj)-A(ii,jj-1)
     enddo
  enddo

  ! /* set Cheb. <-- Leg. matrix, B */      
  call dcopy(m*m,      0.d0,0,B, 1)
  call dcopy(nevn*nevn,0.d0,0,LCBE,1)
  call dcopy((nodd+1)*(nodd+1),0.d0,0,LCBO,1)
  ! note: B=inv(A)
  call dcopy(m*m,A,1,B,1)
  call dgetrf(m,m,B,m,ipiv,info)      ! OK
  call dgetri(m,B,m,ipiv,wk,-1,info) ! OK
  if ((wk(1)-lwk).gt.0) then
     write(*,*) 'leg_cheb_mat: check lwk=', wk(1), lwk
  end if
  call dgetri(m,B,m,ipiv,wk,lwk,info) ! OK
  !
  deallocate(wk)

  if (ngpu.gt.0) then
     LCA(1:n,1:n)=A(1:n,1:n)
     LCB(1:n,1:n)=B(1:n,1:n)
  end if

  !c     /* split into full triangular matrices */
  !c     /* AE=A(1:2:m,1:2:m); AO=A(2:2:m,2:2:m); */
  !c     /* BE=B(1:2:m,1:2:m); BO=B(2:2:m,2:2:m); */
  do j=1,nevn
     jj=(j-1)*2+1
     do i=1,nevn
        ii=(i-1)*2+1
        LCAE(i,j)=A(ii,jj)
     enddo
  enddo
  do j=1,nevn
     jj=(j-1)*2+1
     do i=1,nevn
        ii=(i-1)*2+1
        LCBE(i,j)=B(ii,jj)
     enddo
  enddo
  do j=1,nodd
     jj=j*2
     do i=1,nodd
        ii=i*2
        LCAO(i,j)=A(ii,jj)
     enddo
  enddo
  do j=1,nodd
     jj=j*2
     do i=1,nodd
        ii=i*2
        LCBO(i,j)=B(ii,jj)
     enddo
  enddo

  if (ngpu.gt.0) then

#ifdef CUBLAS
     istat= cublas_alloc((mye+1)*(mye+1), nk, dp_LCAE)
     if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_LCAE',istat
     istat= cublas_set_matrix(mye+1,mye+1,nk,LCAE,mye+1,dp_LCAE,mye+1)
     if (istat.ne.0) write(*,*) 'gpu setmat error: dp_LCAE',istat
     istat= cublas_alloc((mye+1)*(mye+1), nk, dp_LCBE)
     if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_LCBE',istat
     istat= cublas_set_matrix(mye+1,mye+1,nk,LCBE,mye+1,dp_LCBE,mye+1)
     if (istat.ne.0) write(*,*) 'gpu setmat error: dp_LCBE',istat

     istat= cublas_alloc((mye+1)*(mye+1), nk, dp_LCAO)
     if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_LCAO',istat
     istat= cublas_set_matrix(mye+1,mye+1,nk,LCAO,mye+1,dp_LCAO,mye+1)
     if (istat.ne.0) write(*,*) 'gpu setmat error: dp_LCAO',istat
     istat= cublas_alloc((mye+1)*(mye+1), nk, dp_LCBO)
     if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_LCBO',istat
     istat= cublas_set_matrix(mye+1,mye+1,nk,LCBO,mye+1,dp_LCBO,mye+1)
     if (istat.ne.0) write(*,*) 'gpu setmat error: dp_LCBO',istat

     istat= cublas_alloc((my+1)*(my+1), nk, dp_LCA)
     if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_LCA',istat
     istat= cublas_set_matrix(my+1,my+1,nk,LCA,my+1,dp_LCA,my+1)
     if (istat.ne.0) write(*,*) 'gpu setmat error: dp_LCA',istat

     istat= cublas_alloc((my+1)*(my+1), nk, dp_LCB)
     if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_LCB',istat
     istat= cublas_set_matrix(my+1,my+1,nk,LCB,my+1,dp_LCB,my+1)
     if (istat.ne.0) write(*,*) 'gpu setmat error: dp_LCB',istat
#endif

  end if

  allocate (tmpye(mye+1,mz),tmpyo(mye+1,mz),tmpze(my,mze+1),tmpzo(my,mze+1) )
  tmpye=0.d0; tmpyo=0.d0; tmpze=0.d0; tmpzo=0.d0
  if (ngpu.gt.0) then
#ifdef CUBLAS
     ! allocate tmp arrays on device
     istat= cublas_alloc((my+1)*(mz+1), nk, devPtrU)
     if (istat.ne.0) write(*,*) 'gpu alloc. error: devPtrU',istat
     istat= cublas_alloc((my+1)*(mz+1), nk, devPtrV)
     if (istat.ne.0) write(*,*) 'gpu alloc. error: devPtrV',istat

     istat= cublas_alloc((mye+1)*mz, nk, dp_tmpye)
     if (istat.ne.0) write(*,*) 'gpu alloc. error: tmpye',istat
     istat= cublas_alloc((mye+1)*mz, nk, dp_tmpyo)
     if (istat.ne.0) write(*,*) 'gpu alloc. error: tmpy0',istat
     istat= cublas_alloc(my*(mze+1), nk, dp_tmpze)
     if (istat.ne.0) write(*,*) 'gpu alloc. error: tmpze',istat
     istat= cublas_alloc(my*(mze+1), nk, dp_tmpzo)
     if (istat.ne.0) write(*,*) 'gpu alloc. error: tmpzo',istat
#endif

  end if

  return
end subroutine leg_cheb_mat
!----|----------------------------------------------------------------|
subroutine chebleg2d(uc,ul,iopt)

  use ctes, only:nk,ngpu,my,mz
  use legchebshen
#ifdef CUBLAS
  use gpudevice
#endif
  implicit none
  ! /****************************************************************/
  ! /* chebyshev legendre transformation in shen's book pp.111      */
  ! /*                     (take care for typos in this boook)      */
  ! /* uc(my,mz) : chebyshev coefficients                           */
  ! /* ul(my,mz) : legendre  coefficients                           */
  ! /*                                                              */
  ! /* iopt=+1 : cheb -> leg                                        */
  ! /*      -1 : leg  -> cheb                                       */
  ! /*                                                              */
  ! /*  ! now only for square(my=mz)                                */
  ! /****************************************************************/
  real(nk) uc(my+1,mz+1),ul(my+1,mz+1)
  integer iopt

  if (ngpu.eq.0) then
     if (iopt.gt.0) then    ! leg <-- cheb
        !write(*,*) 'chebleg'
        call chebleg_sweep2d(ul,uc,LCAE,LCAO)
     elseif (iopt.lt.0) then ! cheb <-- leg
        !write(*,*) 'legcheb'
        call chebleg_sweep2d(uc,ul,LCBE,LCBO)
     endif
  elseif (ngpu.gt.0) then
#ifdef CUBLAS
     if (iopt.gt.0) then    ! leg <-- cheb
        call chebleg_sweep2d_gpu(ul,uc,dp_LCAE,dp_LCAO) ! uses do loop
        !call chebleg_sweep2d_gpu2(ul,uc,dp_LCA) ! uses tranpose in blas
     elseif (iopt.lt.0) then ! cheb <-- leg
        call chebleg_sweep2d_gpu(uc,ul,dp_LCBE,dp_LCBO)
        !call chebleg_sweep2d_gpu2(uc,ul,dp_LCB)
     endif
#endif 
  end if

  return

end subroutine chebleg2d
!----|----------------------------------------------------------------|
subroutine chebleg_sweep2d(V,U,AE,AO)
  use ctes,only:nk,my,mye,myo,mz,mze,mzo
  ! do not use legchebshen.mod
  use legchebshen, only : tmpye,tmpyo,tmpze,tmpzo
  use timing
  !
  !! U(f) --> V(g) ! only for my=mz
  implicit none
  real(nk) V(my+1,mz+1),U(my+1,mz+1)
  ! note: mye,myo is for shen's base, M(my-2,mz-2)
  real(nk) AE(mye+1,mye+1),AO(myo+1+1,myo+1+1)
  !real*8 tmpye(mye+1,mz),tmpyo(myo+1,mz) ! allocated in leg_cheb_mat
  !real*8 tmpze(my,mze+1),tmpzo(my,mzo+1)

  integer j,jj
  !
  ! /* first sweep, A*U */
  ! /* n=my */
  ! /* g(1:2:n,1:n)=AE*f(1:2:n,1:n); */
  ! /* g(2:2:n,1:n)=AO*f(2:2:n,1:n); */ 
  call cpu_time(t1)

  do j=1,mz
     call dcopy((mye+1),U(1,j),2,tmpye(1,j),1)
  enddo
  !!!call dcopy((mye+1)*mz,U(1,1),2,tmpye(1,1),1)

  do j=1,mz
     call dcopy((myo+1),U(2,j),2,tmpyo(1,j),1)
  enddo
  !!!call dcopy((myo+1+1)*mz,U(2,1),2,tmpyo(1,1),1) ! wrapped version! (myo+1+1)

  call cpu_time(t2); prof_legcheb(2)=prof_legcheb(2) + (t2-t1)
  
  !  /*  A*U  */
  call cpu_time(t1)
  call dtrmm('l','u','n','n',mye+1,mz,1.d0,AE,mye+1,tmpye,mye+1) 
  call dtrmm('l','u','n','n',myo+1+1,mz,1.d0,AO,myo+1+1,tmpyo,myo+1+1)
  call cpu_time(t2); prof_legcheb(3)=prof_legcheb(3) + (t2-t1)

  !  /* reset full g(1:m,1:m) to V(my,mz) */
  call cpu_time(t1)
  do j=1,mz
     call dcopy((mye+1),tmpye(1,j),1,V(1,j),2)
  enddo
  !!!call dcopy((mye+1)*mz,tmpye(1,1),1,V(1,1),2)
  do j=1,mz
     call dcopy((myo+1),tmpyo(1,j),1,V(2,j),2)
  enddo
  !!!call dcopy((myo+1+1)*mz,tmpyo(1,1),1,V(2,1),2)
  ! 
  call cpu_time(t2); prof_legcheb(4)=prof_legcheb(4) + (t2-t1)
  call cpu_time(t1)
  ! V = transpose(V) ! dgeam (C = alpha*op(A) + beta*op(B)) can also do transpose 
  ! by setting alpha=1 and beta=0
  !
  ! /* second sweep, AU*A' */
  ! /* g(1:n,1:2:n)=g(1:n,1:2:n)*AE'; */
  ! /* g(1:n,2:2:n)=g(1:n,2:2:n)*AO'; */
  jj=0
  do j=1,mz,2
     jj=jj+1
     call dcopy(my,V(1,j),1,tmpze(1,jj),1)
  enddo
  jj=0
  do j=2,mz-1,2
     jj=jj+1
     call dcopy(my,V(1,j),1,tmpzo(1,jj),1)
  enddo
  call cpu_time(t2); prof_legcheb(5)=prof_legcheb(5) + (t2-t1)

  call cpu_time(t1)
  call dtrmm('r','u','t','n',my,mze+1,1.d0,AE,mye+1,tmpze,my) 
  call dtrmm('r','u','t','n',my,mzo+1+1,1.d0,AO,mzo+1+1,tmpzo,my)      
  call cpu_time(t2); prof_legcheb(6)=prof_legcheb(6) + (t2-t1)

  ! /* reset full g(1:m,1:m) to V(my,mz) */
  call cpu_time(t1)
  jj=0
  do j=1,mz,2
     jj=jj+1
     call dcopy(my,tmpze(1,jj),1,V(1,j),1)
  enddo
  jj=0
  do j=2,mz-1,2
     jj=jj+1
     call dcopy(my,tmpzo(1,jj),1,V(1,j),1)
  enddo
  call cpu_time(t2); prof_legcheb(7)=prof_legcheb(7) + (t2-t1)

  return

end subroutine chebleg_sweep2d
!----|----------------------------------------------------------------|
#ifdef CUBLAS
subroutine chebleg_sweep2d_gpu(V,U,dp_AE,dp_AO)
  use ctes
  ! do not use legchebshen.mod
  use legchebshen, only: tmpye,tmpyo,tmpze,tmpzo
  use gpudevice
  use timing
  !
  !! U(f) -->V(g) ! only for my=mz
  implicit none
  real*8 V(my+1,mz+1),U(my+1,mz+1) ! see above leg_cheb_mat
  ! note: mye,myo is for shen's base, M(my-2,mz-2)
  !real*8 AE(my/2+1,my/2+1),AO(my/2,my/2) ! should be Cptr
  integer*8 dp_AE,dp_AO
  !real*8 tmpye(mye+1,mz),tmpyo(myo+1,mz) ! legchebshen.mod
  !real*8 tmpze(my,mze+1),tmpzo(my,mzo+1)

  integer j,jj,istat
  !
  call cpu_time(t1)
  !cublas_set_matrix U-> device
  istat= cublas_set_matrix(my+1,mz+1,nk,U,my+1,devPtrU,my+1)
  call cpu_time(t2); prof_legcheb(1)=prof_legcheb(1) + (t2-t1)
  ! /* first sweep, A*U */
  ! /* n=my */
  ! /* g(1:2:n,1:n)=AE*f(1:2:n,1:n); */
  ! /* g(2:2:n,1:n)=AO*f(2:2:n,1:n); */ 
  call cpu_time(t1)
  call cublas_dcopy((mye+1)*(mz),devPtrU,2,dp_tmpye,1)
  call cublas_dcopy((mye+1)*(mz),devPtrU+nk,2,dp_tmpyo,1)
  call cpu_time(t2); prof_legcheb(2)=prof_legcheb(2) + (t2-t1)

  !  /*  A*U  */
  call cpu_time(t1)
  call cublas_dtrmm('l','u','n','n',mye+1,mz,1.d0,dp_AE,mye+1,dp_tmpye,mye+1) 
  call cublas_dtrmm('l','u','n','n',mye+1,mz,1.d0,dp_AO,mye+1,dp_tmpyo,mye+1)
  call cpu_time(t2); prof_legcheb(3)=prof_legcheb(3) + (t2-t1)

  !  /* reset full g(1:m,1:m) to V(my,mz) */
  call cpu_time(t1)
  call cublas_dcopy((mye+1)*mz,dp_tmpye,1,devPtrV,2) 
  call cublas_dcopy((mye+1)*mz,dp_tmpyo,1,devPtrV+nk,2) 
  call cpu_time(t2); prof_legcheb(4)=prof_legcheb(4) + (t2-t1)
  ! 
  call cpu_time(t1)
#ifdef TRANS_GPU 
  ! dgeam (C = alpha*op(A) + beta*op(B)) can also do transpose with alpha=1, beta=0
  !istat= cublas_get_matrix(my+1,mz+1,nk,devPtrV,my+1,V,my+1) ! device -> cpu-memory
  !write(*,'(a,5(1x,g15.8))') 'V(1:5,1)=', V(1,1),V(2,1),V(3,1),V(4,1),V(5,1)
  !write(*,'(a,5(1x,g15.8))') 'V(1:5,2)=', V(1,2),V(2,2),V(3,2),V(4,2),V(5,2)
  !write(*,'(a,5(1x,g15.8))') 'V(1:5,3)=', V(1,3),V(2,3),V(3,3),V(4,3),V(5,3)
  !write(*,'(a,5(1x,g15.8))') 'V(1:5,4)=', V(1,4),V(2,4),V(3,4),V(4,4),V(5,4)
  !write(*,'(a,5(1x,g15.8))') 'V(1:5,5)=', V(1,5),V(2,5),V(3,5),V(4,5),V(5,5)
  ! write(*,*) 'pointers', devPtrV, devPtrU
  istat = transpose_gpu(handle,int(my+1,kind=8),int(mz+1,kind=8),1.d0,devPtrV,devPtrU) ! V^T --> U 
  !call transpose_gpu(handle,my+1,mz+1,1.d0,devPtrV,devPtrU) ! V^T --> U 
  !write(*,*) 'istat =', istat
  !istat = cublas_dgeam(handle,my+1,mz+1,'t','n',1.d0,devPtrV,my+1,0.d0,devPtrV,my+1,devPtrU,my+1) ! failed 
  !istat= cublas_get_matrix(my+1,mz+1,nk,devPtrU,my+1,U,my+1) ! device -> cpu-memory
  !write(*,'(a,5(1x,g15.8))') 'U(1,1:5)=', U(1,1),U(1,2),U(1,3),U(1,4),U(1,5)
  !write(*,'(a,5(1x,g15.8))') 'U(2,1:5)=', U(2,1),U(2,2),U(2,3),U(2,4),U(2,5)
  !write(*,'(a,5(1x,g15.8))') 'U(3,1:5)=', U(3,1),U(3,2),U(3,3),U(3,4),U(3,5)
  !write(*,'(a,5(1x,g15.8))') 'U(4,1:5)=', U(4,1),U(4,2),U(4,3),U(4,4),U(4,5)
  !write(*,'(a,5(1x,g15.8))') 'U(5,1:5)=', U(5,1),U(5,2),U(5,3),U(5,4),U(5,5)
  !stop
  !
  ! quick split
  call cublas_dcopy((my)*(mze+1),devPtrU,2,dp_tmpye,1)
  call cublas_dcopy((my)*(mze+1),devPtrU+nk,2,dp_tmpyo,1) ! (mzo+1)+1 = mze+1
  call cpu_time(t2); prof_legcheb(5)=prof_legcheb(5) + (t2-t1)
  !
  ! /* second sweep, AU*A' */
  call cpu_time(t1)
  ! call cublas_dtrmm('r','u','t','n',my,mze+1,1.d0,dp_AE,mye+1,dp_tmpze,my) 
  ! call cublas_dtrmm('r','u','t','n',my,mze+1,1.d0,dp_AO,mye+1,dp_tmpzo,my)  
  call cublas_dtrmm('l','u','n','n',mze+1,my,1.d0,dp_AE,mye+1,dp_tmpye,mye+1) ! for transposed one   
  call cublas_dtrmm('l','u','n','n',mze+1,my,1.d0,dp_AO,mye+1,dp_tmpyo,mye+1) ! for transposed one
  call cpu_time(t2); prof_legcheb(6)=prof_legcheb(6)+(t2-t1)

  call cpu_time(t1)
  call cublas_dcopy((mze+1)*my,dp_tmpye,1,devPtrU,2) 
  call cublas_dcopy((mze+1)*my,dp_tmpyo,1,devPtrU+nk,2) 
  call cpu_time(t2); prof_legcheb(7)=prof_legcheb(7) + (t2-t1)
  call cpu_time(t1)
  istat = transpose_gpu(handle,int(my+1,kind=8),int(mz+1,kind=8),1.d0,devPtrU,devPtrV) ! U^T --> V 
  call cpu_time(t2); prof_legcheb(5)=prof_legcheb(5) + (t2-t1)
#else
  ! /* transpose split for the second sweep, AU*A' */
  ! /* g(1:n,1:2:n)=g(1:n,1:2:n)*AE'; */
  ! /* g(1:n,2:2:n)=g(1:n,2:2:n)*AO'; */
  jj=0
  do j=1,mz,2
     jj=jj+1
     call cublas_dcopy(my,devPtrV+(j-1)*(my+1)*nk,1,dp_tmpze+(jj-1)*my*nk,1) ! slow
  enddo
  jj=0
  do j=2,mz-1,2
     jj=jj+1
     call cublas_dcopy(my,devPtrV+(j-1)*(my+1)*nk,1,dp_tmpzo+(jj-1)*my*nk,1) ! slow
  enddo  
  call cpu_time(t2); prof_legcheb(5)=prof_legcheb(5) + (t2-t1)
  !
  ! /* second sweep, AU*A' */
  call cpu_time(t1)
  call cublas_dtrmm('r','u','t','n',my,mze+1,1.d0,dp_AE,mye+1,dp_tmpze,my) 
  call cublas_dtrmm('r','u','t','n',my,mze+1,1.d0,dp_AO,mye+1,dp_tmpzo,my)  
  call cpu_time(t2); prof_legcheb(6)=prof_legcheb(6)+(t2-t1)

  call cpu_time(t1)
  ! /* reset full g(1:m,1:m) to V(my,mz) */
  jj=0
  do j=1,mz,2
     jj=jj+1
     call cublas_dcopy(my+1,dp_tmpze+(jj-1)*my*nk,1,devPtrV+(j-1)*(my+1)*nk,1) ! slow
  enddo
  jj=0
  do j=2,mz-1,2
     jj=jj+1
     call cublas_dcopy(my+1,dp_tmpzo+(jj-1)*my*nk,1,devPtrV+(j-1)*(my+1)*nk,1) ! slow
  enddo
  call cpu_time(t2); prof_legcheb(7)=prof_legcheb(7)+(t2-t1)
#endif /* TRANS_GPU */
  call cpu_time(t1)
  ! cublas_get_matrix V <- device
  istat= cublas_get_matrix(my+1,mz+1,nk,devPtrV,my+1,V,my+1) ! device -> cpu-memory
  call cpu_time(t2); prof_legcheb(8)=prof_legcheb(8) + (t2-t1)

  return

end subroutine chebleg_sweep2d_gpu
#endif /* CUBLAS */
!----|----------------------------------------------------------------|
#ifdef CUBLAS
subroutine chebleg_sweep2d_gpu2(V,U,dp_A)
  use ctes
  ! do not use legchebshen.mod
  use legchebshen, only: tmpye,tmpyo,tmpze,tmpzo
  use gpudevice
  use timing
  !
  ! U(f) -->V(g) ! only for my=mz
  implicit none
  real*8 V(my+1,mz+1),U(my+1,mz+1)
  integer*8, value :: dp_A

  integer j,jj,istat
  !
  call cpu_time(t1)
  ! cublas_set_matrix U-> device
  istat= cublas_set_matrix(my+1,mz+1,nk,U,my+1,devPtrV,my+1)
  call cpu_time(t2); prof_legcheb(1)=prof_legcheb(1) + (t2-t1)
  ! /* first direct sweep, A*U */
  ! /* n=my */
  ! /* g(1:2:n,1:n)=AE*f(1:2:n,1:n); */
  ! /* g(2:2:n,1:n)=AO*f(2:2:n,1:n); */ 
  call cpu_time(t1)

  ! /*  A*U  */
  call cpu_time(t1)
  call cublas_dtrmm('l','u','n','n',my+1,mz+1,1.d0,dp_A,my+1,devPtrV,my+1) 
  call cpu_time(t2); prof_legcheb(3)=prof_legcheb(3) + (t2-t1)

  call cpu_time(t1)
  ! /* second sweep, AU*A' */
  ! /* g(1:n,1:2:n)=g(1:n,1:2:n)*AE'; */
  ! /* g(1:n,2:2:n)=g(1:n,2:2:n)*AO'; */
  !
  call cpu_time(t1)
  call cublas_dtrmm('r','u','t','n',my+1,mz+1,1.d0,dp_A,my+1,devPtrV,my+1) 
  call cpu_time(t2); prof_legcheb(6)=prof_legcheb(6)+(t2-t1)

  ! call cpu_time(t1)
  ! istat= cublas_get_matrix(my,mze+1,nk,dp_tmpze,my,tmpze,my) ! device -> cpu-memory
  ! istat= cublas_get_matrix(my,mzo+1,nk,dp_tmpzo,my,tmpzo,my) ! device -> cpu-memory
  ! call cpu_time(t2); prof_legcheb(8)=prof_legcheb(8) + (t2-t1)

  call cpu_time(t1)
  ! cublas_get_matrix V <- device
  istat= cublas_get_matrix(my+1,mz+1,nk,devPtrV,my+1,V,my+1) ! device -> cpu-memory
  call cpu_time(t2); prof_legcheb(8)=prof_legcheb(8) + (t2-t1)

  return

end subroutine chebleg_sweep2d_gpu2
#endif /* CUBLAS */
!----|----------------------------------------------------------------|
subroutine legshen2d(ul,us,bound,iopt)
  
  use ctes
  use legchebshen

  implicit none
  ! /****************************************************************/
  ! /* Legendre-ShenBase transformation                             */
  ! /*                                                              */
  ! /* ul(my,mz)     : Legendre coefficients                        */
  ! /* us(my-2,mz-2) : coefficients of shen's base function         */
  ! /* bound   : boundary condition                                 */
  ! /*        note : only for square homogeneous bc                 */
  ! /*                                                              */
  ! /* iopt=+1 : ul:leg -> us:shen                                  */
  ! /*      -1 : ul:leg <- us:shen                                  */
  ! /*                                                              */
  ! /****************************************************************/
  !    common /shenbkdir/ bkd(my)
  !    common /shenbkneu/ bkn(my)
  real*8 ul(my+1,mz+1),us(my-2+1,mz-2+1)
  character*1 bound
  integer iopt

  if (bound.eq.'d') then
     if (iopt.gt.0) then    ! leg --> shen
        call leg2shen_sweep2d(us,ul,bkd,my)
     elseif(iopt.lt.0) then ! leg <-- shen
        call shen2leg_sweep2d(us,ul,bkd,my)
     endif
  elseif (bound.eq.'n') then
     if (iopt.gt.0) then    ! leg --> shen
        call leg2shen_sweep2d(us,ul,bkn,my)
     elseif(iopt.lt.0) then ! leg <-- shen
        call shen2leg_sweep2d(us,ul,bkn,my)
     endif
  endif

  return
end subroutine legshen2d
! -------------------------------------------------------------------
subroutine leg2shen_sweep2d(s,g,b,m)  !g-->s

  use ctes,only:nk,my,mz
  ! use legchebshen

  real(nk) g(my+1,mz+1),s(my-2+1,mz-2+1)
  real(nk) b(my)
  integer m,jj,j,ii,i
  real(nk) ta,tb,tc,td

  ! only for square homogeneous boundary condition
  call dcopy((m-2+1)*(m-2+1),0.d0,0,s,1)
      
  do jj=1,m-2
     j=jj-1
     ta=2.d0/dfloat(2*j+1) ! MD
     tb=2.d0*b(jj)/dfloat(2*j+5)
     do ii=1,m-2
        i=ii-1
        tc=2.d0/dfloat(2*i+1)
        td=2.d0*b(ii)/dfloat(2*i+5)
        s(ii,jj)=ta*(g(ii,jj)*tc + g(ii+2,jj)*td) &
             &   + tb*(g(ii,jj+2)*tc + g(ii+2,jj+2)*td)
     enddo
  enddo

  return
end subroutine leg2shen_sweep2d
! ------------------------------------------------------------------
#ifdef CUBLAS
subroutine leg2shen_sweep2d_gpu(s,g,b,m)  !g-->s

  ! dp_tmpze, dp_tmpzo are used as the GPU device pointer 
  ! (see chebleg_sweep2d_gpu) 
  use ctes,only:nk,my,mz
  ! use legchebshen

  real(nk) g(my+1,mz+1),s(my-2+1,mz-2+1)
  real(nk) b(my)
  integer m,jj,j,ii,i
  real(nk) ta,tb,tc,td

  ! only for square homogeneous boundary condition
  call dcopy((m-2+1)*(m-2+1),0.d0,0,s,1)
      
  do jj=1,m-2
     j=jj-1
     ta=2.d0/dfloat(2*j+1) ! MD
     tb=2.d0*b(jj)/dfloat(2*j+5)
     do ii=1,m-2
        i=ii-1
        tc=2.d0/dfloat(2*i+1)
        td=2.d0*b(ii)/dfloat(2*i+5)
        s(ii,jj)=ta*(g(ii,jj)*tc + g(ii+2,jj)*td) &
             &   + tb*(g(ii,jj+2)*tc + g(ii+2,jj+2)*td)
     enddo
  enddo

  return
end subroutine leg2shen_sweep2d_gpu
#endif /* CUBLAS */
!----|----------------------------------------------------------------|
subroutine shen2leg_sweep2d(s,g,b,m)  !s-->g

  use ctes, only:nk,my,mz
  !use legchebshen

  implicit none

  real(nk) g(my+1,mz+1),s(my-2+1,mz-2+1)
  real(nk) b(my)
  real(nk) sr(my+2,mz+2)
  real(nk) bs(my)
  integer m,j,i,ii,jj
  real(nk) bi,bj

  ! only for square homogeneous boundary condition
  ! /* wrapping b and s by zero, set to br sr*/
  bs(1)=0.d0
  bs(2)=0.d0
  call dcopy(m-2,b,1,bs(3),1)
  call dcopy((m+2)*(m+2),0.d0,0,sr,1)
  do j=1,m-2
     do i=1,m-2
        sr(i+2,j+2)=s(i,j)
     enddo
  enddo
  ! /* sweep */      
  do ii=1,m
     bi=bs(ii)
     do jj=1,m
        bj=bs(jj)
        g(ii,jj)=sr(ii+2,jj+2) + sr(ii+2,jj)*bj &
             &   + sr(ii,jj+2)*bi + bi*bj*sr(ii,jj)
     enddo
  enddo

  return
end subroutine shen2leg_sweep2d
!----|----------------------------------------------------------------|
