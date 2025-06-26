#ifdef HELM_SHEN
subroutine preshen(M2D,S2D,bound)

  use ctes, only:nk,myid,master,my,mz,mye,myo,mze,mzo,ngpu
  use legchebshen
  use gpudevice

  implicit none
  !  include 'ctes3D'
  ! /*****************************************************************/
  ! /* pre-processing Shen's base 'mass' matrix and                  */
  ! /*                            'stiffness' matrix components      */
  ! /*                             (decomposed to even-odd modes)    */
  ! /* and compute generallized eigen value problem Mx=(lambda)Sx    */
  ! /*                                                               */
  ! /* bound = 'd' or 'n' note: only for homogeneous b.c.            */
  ! /*                                                               */
  ! /*****************************************************************/
  !implicit real*8(a-h,o-z)
  !    common /shenbkdir/ bkd(my)
  !    common /shenbkneu/ bkn(my)
  !    common /shenmatdir/ MED(mye),MOD(myo),
  !   &                    MPPED(mye-1),MPPOD(myo-1),
  !   &                    SED(mye),SOD(myo)
  !    real*8 MED,MOD,MPPED,MPPOD    ! this line is very important
  !    common /shenmatneu/ MEN(mye),MON(myo),
  !   &                    MPPEN(mye-1),MPPON(myo-1),
  !   &                    SEN(mye),SON(myo)
  !    real*8 MEN,MON,MPPEN,MPPON    ! this line is very important
  !    common /sheneigdir/ PD(my-2,my-2),PSM1D(my-2,my-2),ED(my-2)
  !    common /sheneigneu/ PN(my-2,my-2),PSM1N(my-2,my-2),EN(my-2)
  !    common /sheneigdir2/ PED(mye,mye),POD(myo,myo),
  !   &                     PSM1ED(mye,mye),PSM1OD(myo,myo)
  !    common /sheneigneu2/ PEN(mye,myo),PON(myo,mye),
  !   &                     PSM1EN(mye,myo),PSM1ON(myo,mye)
  real*8 M2D(my-2,my-2),S2D(my-2,my-2)
  character bound
  integer m,jj,j,istat,ip

  integer istatus, version

  m=my
  call dcopy((m-2)*(m-2),0.d0,0,M2D,1)
  call dcopy((m-2)*(m-2),0.d0,0,S2D,1)
  if (myid.eq.master) write(*,*) 'preshen bound =', bound

  if ((bound.eq.'d').or.(bound.eq.'_')) then
     if (myid.eq.master) write(*,*) 'initializing galerkin array for Direchlet '
     allocate(bkd(my))
     bkd=0.d0
     allocate(MED(mye),MDO(myo),MPPED(mye-1),MPPOD(myo-1),SED(mye),SOD(myo))
     MED=0.d0; MDO=0.d0; MPPED=0.d0; MPPOD=0.d0; SED=0.d0; SOD=0.d0
     call makeMS(MED,MDO,MPPED,MPPOD,SED,SOD,bkd,'d',my)
     if (ngpu.gt.0) then
        istat= cublas_alloc(mye, nk, dp_MED)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_MED',istat
        istat= cublas_set_vector(mye,nk,MED,1,dp_MED,1)
        !istat= cublas_alloc(mye, nk, dp_tmp1d)
        !if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_tmp1d',istat
        !istat= cublas_set_vector(mye,nk,MED,1, dp_tmp1d,1)

        !istat= cublas_alloc(mye*(my-2), nk, dp_MED)
        !if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_MED (matrix)',istat
        !istat= cublas_set_matrix(mye,1,nk,MED,1,dp_MED,1)
        !do ip=1,my-2
        !   ! /* DPE=(rp*alpha+1)*MPPED */
        !   call cublas_dcopy(mye,    dp_tmp1d ,1,dp_MED+(mye)*(ip-1)*nk,1)
        !end do


        istat= cublas_alloc(myo, nk, dp_MDO)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_MDO',istat
        istat= cublas_set_vector(myo,nk,MDO,1,dp_MDO,1)

        istat= cublas_alloc(mye-1, nk, dp_MPPED)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_MPPED',istat
        istat= cublas_set_vector(mye-1,nk,MPPED,1,dp_MPPED,1)

        istat= cublas_alloc(myo-1, nk, dp_MPPOD)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_MPPOD',istat
        istat= cublas_set_vector(myo-1,nk,MPPOD,1,dp_MPPOD,1)

        istat= cublas_alloc(mye, nk, dp_SED)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_SED',istat
        istat= cublas_set_vector(mye,nk,SED,1,dp_SED,1)

        istat= cublas_alloc(myo, nk, dp_SOD)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_SOD',istat
        istat= cublas_set_vector(myo,nk,SOD,1,dp_SOD,1)

     end if

     allocate(PD(my-2,my-2),PSM1D(my-2,my-2),ED(my-2))
     PD=0.d0; PSM1D=0.d0; ED=0.d0
     !
     call precomeigen(M2D,S2D,ED,PD,PSM1D,'d',my)


  end if

  if ((bound.eq.'n').or.(bound.eq.'_')) then
     if (myid.eq.master) write(*,*) 'initializing galerkin array for Neumann '
     allocate(bkn(my))
     bkn=0.d0;
     allocate(MEN(mye),MON(myo), MPPEN(mye-1),MPPON(myo-1),SEN(mye),SON(myo))
     MEN=0.d0; MON=0.d0; MPPEN=0.d0; MPPON=0.d0; SEN=0.d0; SON=0.d0
     call makeMS(MEN,MON,MPPEN,MPPON,SEN,SON,bkn,'n',my)
     allocate(PN(my-2,my-2),PSM1N(my-2,my-2),EN(my-2))
     PN=0.d0; PSM1N=0.d0; EN=0.d0
     call precomeigen(M2D,S2D,EN,PN,PSM1N,'n',my)
     if (ngpu.gt.0) then
        istat= cublas_alloc(my-2, nk, dp_EN)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_EN',istat
        istat= cublas_set_vector(my-2,nk,EN,1,dp_EN,1)
     end if
  endif
  ! 
  allocate(gp(my-2),vp(my-2))
  gp=0.d0; vp=0.d0
  allocate(fe(mye),fo(myo))
  fe=0.d0; fo=0.d0
  if (ngpu.gt.0) then
     istat= cublas_alloc(my-2, nk, dp_gp)
     istat= cublas_alloc(my-2, nk, dp_vp)             
     istat= cublas_alloc(mye, nk, dp_fe) ! reveted to 1d
     istat= cublas_alloc(myo, nk, dp_fo)       
  end if
  
  if (ngpu.gt.0) then
     !istat= cublas_alloc(my-2, nk, dp_ED)
     !if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_ED',istat
     !istat= cublas_set_vector(my-2,nk,ED,1,dp_ED,1)

     istat= cublas_alloc(my-2, nk, dp_rp)
     if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_rp',istat
     istat= cublas_set_vector(my-2,nk,ED,1,dp_rp,1)
     
     istat= cublas_alloc(my-2, nk, dp_ra1)
     if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_ra1',istat

  end if

  ! for solV (matrix version)
  allocate(DDEp(mye),DDOp(myo),DPEp(mye-1),DPOp(myo-1))
  !allocate(DDEp(mye,my-2),DDOp(myo,my-2),DPEp(mye-1,my-2),DPOp(myo-1,my-2))
  DDEp=0.d0; DDOp=0.d0; DPEp=0.d0; DPOp=0.d0;

  if (ngpu.gt.0) then
     istat= cublas_alloc((my-2+1)*(mz-2+1), nk, devPtrG)
     if (istat.ne.0) write(*,*) 'cublas alloc. error: devPtrG',istat
     ! for cusparse_dgtsv(me,1,dp_DPEp,dp_DDEp,dp_DPEp2,dp_fe,me,info)
     ! dl,d,du should be me-size and dl(0)=0 and du(end)=0
     
     istatus = cusparse_create(handle)
     !allocate(handle_e(my-2),handle_o(my-2))
     !do ip=1,my-2
     !   istatus = cusparse_create(handle_e(ip))
     !   istatus = cusparse_create(handle_o(ip))
     !end do
     
     if (istatus.ne.0) then
        write(*,*) "cusparse library initialization failed "
        stop
     end if
     ! Create a stream for every dcopy operation
     allocate(streams(my-2))
     do ip=1,my-2 
        istat = cuda_stream_create(streams(ip))
     enddo
     !        istatus = cusparse_set_stream(handle_e, stream_e) 
     
     istat= cublas_alloc(mye*(my-2), nk, dp_DDEp)
     if (istat.ne.0) write(*,*) 'cublas alloc. error: DDEp',istat
     istat= cublas_set_matrix(mye,(my-2),nk,0.d0,0,dp_DDEp,1) ! zero set
     
     !istat= cublas_alloc(myo*(my-2), nk, dp_DDOp)
     istat= cublas_alloc((myo+1)*(my-2), nk, dp_DDOp)
     if (istat.ne.0) write(*,*) 'cublas alloc. error: DDOp',istat
     istat= cublas_set_matrix(myo+1,(my-2),nk,0.d0,0,dp_DDOp,1) ! zero set
     
     istat= cublas_alloc(mye*(my-2), nk, dp_DPEp)
     if (istat.ne.0) write(*,*) 'cublas alloc. error: DPEp',istat
     !istat= cublas_set_vector(mye,nk,0.d0,0,dp_DPEp,1) ! zero set
     istat= cublas_set_vector(mye,(my-2),nk,0.d0,0,dp_DPEp,1) ! zero set
     
     istat= cublas_alloc((myo+1)*(my-2), nk, dp_DPOp)
     if (istat.ne.0) write(*,*) 'cublas alloc. error: DPOp',istat
     !istat= cublas_set_vector(myo,nk,0.d0,0,dp_DPOp,1) ! zero set
     istat= cublas_set_vector(myo+1,(my-2),nk,0.d0,0,dp_DPOp,1) ! zero set
     
     ! the symmetric tri-diagonal solver is not available for cusparse
     istat= cublas_alloc(mye*(my-2), nk, dp_DPEp2)
     if (istat.ne.0) write(*,*) 'cublas alloc. error: DPEp2',istat
     istat= cublas_set_vector(mye,(my-2),nk,0.d0,0,dp_DPEp2,1) ! zero set
     
     istat= cublas_alloc((myo+1)*(my-2), nk, dp_DPOp2)
     if (istat.ne.0) write(*,*) 'cublas alloc. error: DPOp2',istat
     istat= cublas_set_vector((myo+1),(my-2),nk,0.d0,0,dp_DPOp2,1) ! zero set
     
  end if

  ! /* store even-odd decomposed array */  
  if ((bound.eq.'d').or.(bound.eq.'_')) then 
     allocate(PED(mye,mye),POD(myo+1,myo+1),PSM1ED(mye,mye),PSM1OD(myo+1,myo+1))
     PED=0.d0; POD=0.d0; PSM1ED=0.d0; PSM1OD=0.d0
     !  /*   PE=P(1:2:m,1:2:m); */
     !  /*   PM1E=PM1SM1(1:2:m,1:2:m);  */
     jj=0
     do j=1,mz-2,2
        jj=jj+1
        call dcopy(mye,PD(1,j),2,PED(1,jj),1) 
        call dcopy(mye,PSM1D(1,j),2,PSM1ED(1,jj),1) 
     enddo
     if (ngpu.gt.0) then
        istat= cublas_alloc(mye*mye, nk, dp_PED)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_PED',istat
        istat= cublas_set_matrix(mye,mye,nk,PED,mye,dp_PED,mye)
        if (istat.ne.0) write(*,*) 'gpu setmat error: dp_PED',istat
        istat= cublas_alloc(mye*mye, nk, dp_PSM1ED)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_PSM1ED',istat
        istat= cublas_set_matrix(mye,mye,nk,PSM1ED,mye,dp_PSM1ED,mye)
        if (istat.ne.0) write(*,*) 'gpu setmat error: dp_PSM1ED',istat
     end if
     !
     !  /*   PO=P(2:2:m,2:2:m); */
     !  /*   PM1O=PM1SM1(2:2:m,2:2:m);  */
     jj=0
     do j=2,mz-3,2
        jj=jj+1
        call dcopy(myo,PD(2,j),2,POD(1,jj),1) 
        call dcopy(myo,PSM1D(2,j),2,PSM1OD(1,jj),1) 
     enddo
     if (ngpu.gt.0) then
        istat= cublas_alloc((myo+1)*(myo+1), nk, dp_POD)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_POD',istat
        istat= cublas_set_matrix(myo+1,myo+1,nk,POD,myo+1,dp_POD,myo+1)
        if (istat.ne.0) write(*,*) 'gpu setmat error: dp_POD',istat
        istat= cublas_alloc((myo+1)*(myo+1), nk, dp_PSM1OD)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_PSM1OD',istat
        istat= cublas_set_matrix(myo+1,myo+1,nk,PSM1OD,myo+1,dp_PSM1OD,myo+1)
        if (istat.ne.0) write(*,*) 'gpu setmat error: dp_PSM1OD',istat
     end if
     !
  end if

  if ((bound.eq.'n').or.(bound.eq.'_')) then
     allocate(PEN(mye,myo+1),PON(myo+1,mye),PSM1EN(mye,myo+1),PSM1ON(myo+1,mye))
     PEN=0.d0; PON=0.d0; PSM1EN=0.d0; PSM1ON=0.d0
     !  /*   PE=P(1:2:m,2:2:m); */
     !  /*   PM1E=PM1SM1(1:2:m,2:2:m); */
     jj=0
     do j=2,mz-3,2
        jj=jj+1
        call dcopy(mye,PN(1,j),2,PEN(1,jj),1) 
        call dcopy(mye,PSM1N(1,j),2,PSM1EN(1,jj),1) 
     enddo
     if (ngpu.gt.0) then        
        istat= cublas_alloc(mye*mye, nk, dp_PEN)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_PEN',istat
        istat= cublas_set_matrix(mye,mye,nk,PEN,mye,dp_PEN,mye)
        if (istat.ne.0) write(*,*) 'gpu setmat error: dp_PEN',istat
        istat= cublas_alloc(mye*mye, nk, dp_PSM1EN)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_PSM1EN',istat
        istat= cublas_set_matrix(mye,mye,nk,PSM1EN,mye,dp_PSM1EN,mye)
        if (istat.ne.0) write(*,*) 'gpu setmat error: dp_PSM1EN',istat
     end if
     !
     !  /*   PO=P(2:2:m,1:2:m); */
     !  /*   PM1O=PM1SM1(2:2:m,1:2:m); */
     jj=0
     do j=1,mz-2,2
        jj=jj+1
        call dcopy(myo,PN(2,j),2,PON(1,jj),1) 
        call dcopy(myo,PSM1N(2,j),2,PSM1ON(1,jj),1) 
     enddo
     if (ngpu.gt.0) then
        istat= cublas_alloc((myo+1)*(myo+1), nk, dp_PON)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_PON',istat
        istat= cublas_set_matrix(myo+1,myo+1,nk,PON,myo+1,dp_PON,myo+1)
        if (istat.ne.0) write(*,*) 'gpu setmat error: dp_PON',istat
        istat= cublas_alloc((myo+1)*(myo+1), nk, dp_PSM1ON)
        if (istat.ne.0) write(*,*) 'gpu alloc. error: dp_PSM1ON',istat
        istat= cublas_set_matrix(myo+1,myo+1,nk,PSM1ON,myo+1,dp_PSM1OD,myo+1)
        if (istat.ne.0) write(*,*) 'gpu setmat error: dp_PSM1ON',istat        
     end if

  endif

  ! wrapped version
  allocate(tmpe(mye,mz-2),tmpo(myo+1,mz-2),tmpee(mye,my-2),tmpoo(myo+1,my-2))
  tmpe=0.d0; tmpo=0.d0; tmpee=0.d0; tmpoo=0.d0

  if (ngpu.gt.0) then
     istat= cublas_alloc((my-2+1)*(mz-2+1), nk, devPtrF)
     if (istat.ne.0) write(*,*) 'cublas alloc. error: devPtrF',istat
     !cublas_set_matrix -> device
     istat= cublas_alloc(mye*(my-2), nk, dp_tmpe)
     if (istat.ne.0) write(*,*) 'cublas alloc. error: tmpe',istat
     istat= cublas_set_matrix(mye,my-2,nk,tmpe,mye,dp_tmpe,mye)
     if (istat.ne.0) write(*,*) 'cublas setmat error: tmpe',istat

     istat= cublas_alloc((myo+1)*(my-2), nk, dp_tmpo)
     if (istat.ne.0) write(*,*) 'cublas alloc. error: tmpo',istat
     istat= cublas_set_matrix(myo+1,my-2,nk,tmpo,myo+1,dp_tmpo,myo+1)
     if (istat.ne.0) write(*,*) 'cublas setmat error: tmpo',istat
     
     istat= cublas_alloc(mye*(my-2), nk, dp_tmpee)
     if (istat.ne.0) write(*,*) 'cublas alloc. error: tmpee',istat
     istat= cublas_set_matrix(mye,my-2,nk,tmpee,mye,dp_tmpee,mye)
     if (istat.ne.0) write(*,*) 'cublas setmat error: tmpee',istat

     istat= cublas_alloc((myo+1)*(my-2), nk, dp_tmpoo)
     if (istat.ne.0) write(*,*) 'cublas alloc. error: tmpoo',istat
     istat= cublas_set_matrix((myo+1),my-2,nk,tmpoo,myo+1,dp_tmpoo,myo+1)
     if (istat.ne.0) write(*,*) 'cublas setmat error: tmpoo',istat

  end if

  allocate(DDE(mye),DDO(myo))
  DDE=0.d0; DDO=0.d0;
  if (ngpu.gt.0) then
     !cublas_set_matrix -> device
     istat= cublas_alloc(mye, nk, dp_DDE)
     if (istat.ne.0) write(*,*) 'cublas alloc. error: DDE',istat
     istat= cublas_alloc(myo, nk, dp_DDO)
     if (istat.ne.0) write(*,*) 'cublas alloc. error: DDO',istat
  end if

  return
end subroutine preshen
!----|-----------------------------------------------------------------|
subroutine makeMS(Mevn,Modd,MPPevn,MPPodd, &
     &            Sevn,Sodd,b,bound,m)

  use ctes
  use legchebshen

  ! /*****************************************************************/
  ! /* compute Shen's base 'mass' matrix and                         */
  ! /*                     'stiffness' matrix components             */
  ! /*         decomposed to even-odd modes                          */
  ! /*                                                               */
  ! /* MD(my-2) : diagonal components of mass matrix                 */
  ! /* MPP(my-4) : upper upper diagonal components of mass matrix    */
  ! /* Stiff(my-2) : diagonal components of stiff matrix             */
  ! /* b(my-2) : depends on boundary condition                       */
  ! /*          (a_k is zero for 1-dir homogeneous b.c.,             */ 
  ! /*               pp.101-102 in shen's book)                      */
  ! /* bound = 'd' or 'n' note: only for homogeneous b.c.            */
  ! /*                                                               */
  ! /*****************************************************************/
  !      implicit real*8(a-h,o-z)
  real*8 Mevn(mye),Modd(myo),MPPevn(mye-1),MPPodd(myo-1),Sevn(mye),Sodd(myo)
  real*8 MD(m-2),MPP(m-4),Stiff(m-2)
  real*8 b(m)
  character bound
  integer m,n,nm2,kk,k

  n=m-1
  nm2=n-2

  if (bound.eq.'d') then
     !  /*  in this case b_k's are all -1 */
     do kk=1,m
        b(kk)=-1.d0
     enddo
  else if (bound.eq.'n') then
     do kk=1,m
        k=kk-1
        b(kk)=-dfloat((k)*(k+1))/dfloat((k+2)*(k+3))
     enddo
  else
     write(*,*) 'wrong declaration of bc'
     stop
  endif
  ! /* build matrices */
  ! /* set M, S, and also alpha*M +S.   */
  do kk=1,nm2+1
     k=kk-1
     Stiff(kk) = -b(kk)*dfloat(4*k+6)
     MD(kk) = 2.d0/dfloat(2*k+1) + 2.d0*b(kk)*b(kk)/dfloat(2*k+5)
  enddo
  call dcopy(mye,  MD(1),2, Mevn,1) ! Mass -> Mevn + Modd
  call dcopy(myo,  MD(2),2, Modd,1) !
  
  do kk=1,nm2-1
     k=kk-1
     MPP(kk)=b(kk)*2.d0/dfloat(2*k+5)
  enddo
  call dcopy(mye-1,MPP(1), 2, MPPevn,1) ! MPP -> MPPevn + MPPodd
  call dcopy(myo-1,MPP(2), 2, MPPodd,1) ! 
  call dcopy(mye, Stiff(1),2,Sevn,1)    ! Stiff -> Sevn + Sodd
  call dcopy(myo, Stiff(2),2,Sodd,1)    ! 

  return
end subroutine makeMS
!----|-----------------------------------------------------------------|
subroutine precomeigen(M2D,S2D,Eval,Pvec,PSM1,bound,m)
  use ctes
  implicit none
  ! /*****************************************************************/
  ! /* compute generalized eigen value problem, Mx=(lambda)Sx        */
  ! /*    here, M is 'mass' matrix, S is 'stiffness' matrix          */
  ! /*                                                               */
  ! /* Eval(m) : eigen values sorted by ascending order              */
  ! /* Pvec(m,m): matrix which has eigen vectors                     */
  ! /* PSM1(m,m): inv(P)inv(S)                                       */
  ! /* bound = 'd' or 'n' note: only for homogeneous b.c.            */
  ! /*                                                               */
  ! /*****************************************************************/

  real*8 M2D(my-2,my-2), S2D(my-2,my-2), St(my-2,my-2)
  real*8 Eval(my-2)
  real*8 Pvec(my-2,my-2), PSM1(my-2,my-2)
  
  ! /* tmp arrays and vectors */     
  real*8 M2D1(my-3,my-3), S2D1(my-3,my-3), St1(my-3,my-3)

  real*8 Eval1(my-3), Pvec1(my-3,my-3)
  real*8 eii(my-2) !tmp vector for dgeev
  real*8, dimension(:),allocatable::wk
  real*8 b(my)
  integer ipiv(my-2)
  character*1 bound
  integer m,n,nm2,lwk,mm,info
  integer kk,k,i,j
  
  n=m-1
  nm2=n-2
  ! note: [nm2+1=my-2]

  lwk=(m-2)*(m-2)
  allocate(wk(lwk))
  wk=0.d0

  mm=(nm2+1)*(nm2+1)
#ifdef DEBUG      
  write(*,*) 'in precomp_eigen'
#endif
  if (bound.eq.'d') then
     ! in this case b_k's are all -1
     ! S is diagonal matrix
     do kk=1,m
        b(kk)=-1.d0
     enddo
  elseif (bound.eq.'n') then
     do kk=1,m
        k=kk-1
        b(kk)=-dfloat(k*(k+1))/dfloat((k+2)*(k+3))
     enddo
  else 
     write(*,*) 'precomeigen: bound should be d or n, stop'
     stop
  endif
  ! build matrices
  !      call dcopy(mm,0.d0,0,M2D,1) !segmentation fault for big my??
  !      call dcopy(mm,0.d0,0,S2D,1)
  !      call dcopy(mm,0.d0,0,Pvec,1)
  do kk=1,nm2+1
     k=kk-1
     S2D(kk,kk)=-b(kk)*dfloat(4*k+6)
     M2D(kk,kk)=2.d0/dfloat(2*k+1)+2.d0*b(kk)*b(kk)/dfloat(2*k+5)
  enddo
  if (bound.eq.'n') then
     S2D(1,1)=2.d0
  endif
  ! set upper matrix of M
  do kk=1,nm2-1
     k=kk-1
     M2D(kk,kk+2)=b(kk)*2.d0/dfloat(2*k+5)
     M2D(kk+2,kk)=b(kk)*2.d0/dfloat(2*k+5)
  enddo
  
  if (bound.eq.'d') then
     ! /*   Solve the eigenvalue problem Dirichelt case   */
     call dcopy(mm,M2D,1,Pvec,1)
     call dcopy(mm,S2D,1,St,1)
     call dsygv(1,'V','U',nm2+1,Pvec,nm2+1,St,nm2+1, Eval, wk, lwk, info)
     !$$      do j=1,nm2+1
     !$$         do i=1,nm2+1
     !$$            M2D(i,j)=M2D(i,j)/S2D(j,j)
     !$$         enddo
     !$$      enddo
     !$$      write(*,*) 'M2D'
     !$$      write(*,*) M2D      
     !$$      call dgeev('n','v',nm2+1,M2D,nm2+1,Eval,eii,dummy,1,
     !$$     $           Pvec,nm2+1,wk,lwk,info) !!   Eval is not sorted
     !$$      write(*,*) 'info',info
     
  elseif (bound.eq.'n') then
     ! /* set M1 S1*/
     do i=1,nm2
        do j=1,nm2
           M2D1(i,j)=M2D(i+1,j+1)
           S2D1(i,j)=S2D(i+1,j+1)
        enddo
     enddo
     call dcopy(nm2*nm2,M2D1,1,Pvec1,1)
     call dcopy(nm2*nm2,S2D1,1,St1,1)
     call dsygv(1,'V','U',nm2,Pvec1,nm2, St1,nm2, Eval1, wk, lwk, info)
     ! /* set eigen value, vector matrix and S~ */ 
     Pvec(1,1)=1.d0
     S2D(1,1)=2.d0
     Eval(1)=1.d0
     do i=1,nm2
        Eval(i+1)=Eval1(i)
        do j=1,nm2
           Pvec(i+1,j+1)=Pvec1(i,j)
           S2D(i+1,j+1)=S2D1(i,j)
        enddo
     enddo
     
  endif !/* bound */
  
  call dcopy(mm,Pvec,1,PSM1,1)
  call dgetrf(nm2+1,nm2+1,PSM1,nm2+1,ipiv,info)  ! OK

  call dgetri(nm2+1,PSM1,nm2+1,ipiv,wk,-1,info)
  if ((wk(1)-lwk).gt.0) then
     write(*,*) 'precomeigen: check lwk=', wk(1), lwk
  end if
  call dgetri(nm2+1,PSM1,nm2+1,ipiv,wk,lwk,info) ! OK
  do j=1,nm2+1
     do i=1,nm2+1
        PSM1(i,j)=PSM1(i,j)/S2D(j,j)
        !
        ! note: for neumann case S2D is S~ in shen's book pp.222
        !
     enddo
  enddo

  deallocate(wk)
  
  return
end subroutine precomeigen
#else
subroutine preshen_dummy
end subroutine preshen_dummy
#endif /* HELM_SHEN */
