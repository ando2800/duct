#define TRANS_GPU
!#undef TRANS_GPU

#ifdef HELM_SHEN
subroutine solveshen(U,F,work,bound,xlambda)
  use ctes,only:my,mz,ngpu
  use timing, only: timeshen
  implicit none
#ifdef DEBUG
  
  real*8 t1,t2,t3,t4,t5,t6,t7,t8,t9

#endif
  character*1 bound
  real*8 U(my+1,mz+1),F(my-2+1,mz-2+1)
  real*8 work(my+1,mz+1)
  integer m,mm2 
  real*8 xlambda

  m=my
  mm2=my-2
  ! /* pass r.h.s. to Legendre */       ! only for square case my=mz 
#ifdef DEBUG
  !write(*,*) 'before chebleg2d'
  call cpu_time(t1)
#endif
  call chebleg2d(U,work,+1)  
  ! /* pass r.h.s  to shen's base */
#ifdef DEBUG
  call cpu_time(t2)
  timeshen(2)= timeshen(2) + (t2-t1)
#endif
  !write(*,*) 'before legshen2d'
  call legshen2d(work,F,bound,+1) ! work:Legendre --> F:Shen 
  !write(*,*)  'legshen2d', F(20,20)
#ifdef DEBUG
  call cpu_time(t3)
  timeshen(3)= timeshen(3) + t3-t2
  !write(*,*) 'before G:=inv(P)inv(S)F'
#endif
  ! /* pass to diagonalized space */
  ! /* compute G:=inv(P)inv(S)F */
  ! /* F is overwritten */
  !write(*,*) 'multevenodd'
  if (ngpu.eq.0) then
     call multevenodd(F,bound,+1) ! tmpee,tmpoo in mod.f90 are the results
  elseif (ngpu.gt.0) then 
     call multevenodd_gpu(F,bound,+1)
  end if
#ifdef DEBUG
  call cpu_time(t4)
  timeshen(4)=timeshen(4) + t4-t3
  !write(*,*) 'before solV'
#endif
  ! /* solve (N-2) linear systems  --> vs */
  if (ngpu.eq.0) then
     call solV(F,xlambda,bound) 
  elseif (ngpu.gt.0) then
     call solV_gpu(F,xlambda,bound) 
  end if
  ! F is overwritten with vs

#ifdef DEBUG
  call cpu_time(t5)
  timeshen(5)=timeshen(5) + t5-t4
  !write(*,*) 'before U=PV'
#endif
  ! /* re-set U=PV, back to shen's base */
  ! /* us=P*vs;   vs=G(my-2,my-2) --> us=F(my-2,my-2) */
  if (ngpu.eq.0) then
     call multevenodd(F,bound,-1)
  elseif (ngpu.gt.0) then
     call multevenodd_gpu(F,bound,-1)
  end if
#ifdef DEBUG
  call cpu_time(t6)
  timeshen(6)= timeshen(6) + t6-t5
  !write(*,*) 'before legshen2d'
#endif
  ! /* pass to Legendre */
  call legshen2d(work,F,bound,-1) ! work:Legendre <-- F=us:Shen
#ifdef DEBUG
  call cpu_time(t7)
  timeshen(7)= timeshen(7) + t7-t6
  !write(*,*) 'before chebleg2d'
#endif
  ! /* pass to Chebyshev */
  call chebleg2d(U,work,-1)
#ifdef DEBUG
  call cpu_time(t8)
  timeshen(8)= timeshen(8) + t8-t7
#endif
  return
end subroutine solveshen
!----|----------------------------------------------------------------|
subroutine multevenodd(F,bound,iopt)
  use legchebshen
  use ctes,only:nk,my,mz,mye,myo
  use timing
  implicit none
  !  /****************************************************************/
  !  /* F(my-2,my-2)                                                 */
  !  /* opt= +1: pass F to diagonalized space                        */
  !  /*      -1: back F to shen's base                               */
  !  /*                                                              */
  !  /*   compute C=A*B                                              */
  !  /*   here A is inv(P)inv(S) for forward (opt=+1)                */
  !  /*          is P            for backward (opt=-1)               */
  !  /* bound='d'                                                    */
  !  /*   AE=A(1:2:m,1:2:m)                                          */
  !  /*   AO=A(2:2:m,2:2:m) for Dirichlet case (supposed m is odd)   */
  !  /* bound='n'                                                    */
  !  /*   AE=A(1:2:m,2:2:m)                                          */
  !  /*   AO=A(2:2:m,1:2:m) for Neumann case (supposed m is odd)     */
  !  /*                                                              */
  !  /*   if iwhat=='dir'                                            */
  !  /*      C(1:2:n,1:n)=AE*B(1:2:n,1:n);                           */
  !  /*      C(2:2:n,1:n)=AO*B(2:2:n,1:n);                           */
  !  /*   elseif iwhat=='neu'                                        */
  !  /*      C(1:2:n,1:n)=AE*B(2:2:n,1:n);                           */
  !  /*      C(2:2:n,1:n)=AO*B(1:2:n,1:n);                           */
  !  /*      % repair constant mode of shen base spectra for neumann */
  !  /*      C(1,1:n)=fac*B(1,1:n);                                  */
  !  /*   end                                                        */
  !  /*                                                              */
  !  /****************************************************************/
  real(nk) F(my-2+1,mz-2+1)
  real(nk) bcon(my-2+1)
  character*1 bound
  integer iopt,j,k
  real(nk) fac

  call cpu_time(t1)
  !do k=1,mz-2
  !   call dcopy(mye,F(1,j),2,tmpe(1,j),1)
  !enddo
  call dcopy(mye*(mz-2),F(1,1),2,tmpe(1,1),1)
  !do k=1,mz-2
  !   call dcopy(myo,F(2,j),2,tmpo(1,j),1)
  !enddo
  call dcopy((myo+1)*(mz-2),F(2,1),2,tmpo(1,1),1)
  call cpu_time(t2); prof_shen(2)=prof_shen(2) + (t2-t1)

  if (bound.eq.'d') then
     !   /* C(1:2:n,1:n)=AE*B(1:2:n,1:n) */
     !   /* C(2:2:n,1:n)=AO*B(2:2:n,1:n) */
     call cpu_time(t1)
     if (iopt.gt.0) then
        call dgemm('n','n',mye,my-2,mye,1.d0,PSM1ED,mye, &
             &     tmpe,mye,0.d0,tmpee,mye)
        call dgemm('n','n',myo+1,my-2,myo+1,1.d0,PSM1OD,myo+1, &
             &     tmpo,myo+1,0.d0,tmpoo,myo+1)
     elseif (iopt.lt.0) then
        call dgemm('n','n',mye,my-2,mye,1.d0,PED,mye, &
             &     tmpe,mye,0.d0,tmpee,mye)
        call dgemm('n','n',myo+1,my-2,myo+1,1.d0,POD,myo+1, &
             &     tmpo,myo+1,0.d0,tmpoo,myo+1)         
     else
        write(*,*) 'error in multevenodd'
        stop
     endif
     call cpu_time(t2); prof_shen(3)=prof_shen(3) + (t2-t1)
     call cpu_time(t1)
     !do j=1,my-2
     !   call dcopy(mye,tmpee(1,j),1,F(1,j),2)
     !enddo
     call dcopy(mye*(my-2),tmpee(1,1),1,F(1,1),2)
     !do j=1,my-2
     !   call dcopy(myo,tmpoo(1,j),1,F(2,j),2)
     !enddo
     call dcopy((myo+1)*(my-2),tmpoo(1,1),1,F(2,1),2)
     call cpu_time(t2); prof_shen(4)=prof_shen(4) + (t2-t1)

  elseif (bound.eq.'n') then
     call dcopy(my-2+1,F(1,1),my-2+1,bcon,1)
     ! /* C(1:2:n,1:n)=AE*B(2:2:n,1:n); */
     ! /* C(2:2:n,1:n)=AO*B(1:2:n,1:n); */
     if (iopt.gt.0) then
        fac=0.5d0
        call dgemm('n','n',mye,my-2,myo+1,1.d0,PSM1EN,mye, &
             &     tmpo,myo+1,0.d0,tmpee,mye)
        call dgemm('n','n',myo+1,my-2,mye,1.d0,PSM1ON,myo+1, &
             &     tmpe,mye,0.d0,tmpoo,myo+1)
     elseif (iopt.lt.0) then
        fac=1.d0
        call dgemm('n','n',mye,my-2,myo+1,1.d0,PEN,mye, &
             &     tmpo,myo+1,0.d0,tmpee,mye)
        call dgemm('n','n',myo+1,my-2,mye,1.d0,PON,myo+1, &
             &     tmpe,mye,0.d0,tmpoo,myo+1)         
     else
        write(*,*) 'error in multevenodd '
        stop
     endif
     ! /* reconstruct F(m,m) */     
     !do j=1,my-2
     !   call dcopy(mye,tmpee(1,j),1,F(1,j),2)
     !enddo
     call dcopy(mye*(my-2),tmpee(1,1),1,F(1,1),2)
     !do j=1,my-2
     !   call dcopy(myo,tmpoo(1,j),1,F(2,j),2)
     !enddo
     call dcopy((myo+1)*(my-2),tmpoo(1,1),1,F(2,1),2)
     ! /* repair constant mode of shen base spectra */
     ! /* only for neumann case */
     ! /* C(1,1:n)=fac*B(1,1:n); */
     call daxpy(my-2+1,fac,bcon,1,F(1,1),my-2+1)
  endif
  return
end subroutine multevenodd
!----|----------------------------------------------------------------|
subroutine multevenodd_gpu(F,bound,iopt)
  use legchebshen
  use ctes,only:nk,my,mz,mye,myo
  use gpudevice
  use timing
  !
  implicit none
  !  /****************************************************************/
  !  /* GPU version                                                  */
  !  /****************************************************************/
  real(nk) F(my-2+1,my-2+1)
  real(nk) bcon(my-2)
  character*1 bound
  integer iopt,j,istat
  real(nk) fac
  !integer*8 devPtrF see the module gpudevice in mod.f90
  
  if (bound.eq.'d') then

     !cublas_set_matrix F-> device
     call cpu_time(t1)
     if (iopt.gt.0) istat= cublas_set_matrix(my-2,my-2,nk,F,my-2,devPtrF,my-2)
     call cpu_time(t2); prof_shen(1)=prof_shen(1) + (t2-t1)

     !write(*,*) 'cublas alloc setmat F',istat     
     call cpu_time(t1)
     call cublas_dcopy(mye*(mz-2),devPtrF,2,dp_tmpe,1)
     call cpu_time(t2); prof_shen(2)=prof_shen(2) + (t2-t1)

     call cpu_time(t1)
     call cublas_dcopy((myo+1)*(mz-2),devPtrF+nk,2,dp_tmpo,1)
     call cpu_time(t2); prof_shen(2)=prof_shen(2) + (t2-t1)

     !   /* C(1:2:n,1:n)=AE*B(1:2:n,1:n) */
     !   /* C(2:2:n,1:n)=AO*B(2:2:n,1:n) */
     
     if (iopt.gt.0) then
        call cpu_time(t1)
        call cublas_dgemm('n','n',mye,my-2,mye,1.d0,dp_PSM1ED,mye, &
             &     dp_tmpe,mye,0.d0,dp_tmpee,mye)
        call cpu_time(t2); prof_shen(3)=prof_shen(3) + (t2-t1)
        call cpu_time(t1)
        call cublas_dgemm('n','n',myo+1,my-2,myo+1,1.d0,dp_PSM1OD,myo+1, &
             &     dp_tmpo,myo+1,0.d0,dp_tmpoo,myo+1)
        call cpu_time(t2); prof_shen(3)=prof_shen(3) + (t2-t1)
     elseif (iopt.lt.0) then
        call cpu_time(t1)
        call cublas_dgemm('n','n',mye,my-2,mye,1.d0,dp_PED,mye, &
             &     dp_tmpe,mye,0.d0,dp_tmpee,mye)
        call cpu_time(t2); prof_shen(3)=prof_shen(3) + (t2-t1)
        call cpu_time(t1)
        call cublas_dgemm('n','n',myo+1,my-2,myo+1,1.d0,dp_POD,myo+1, &
             &     dp_tmpo,myo+1,0.d0,dp_tmpoo,myo+1)       
        call cpu_time(t2); prof_shen(3)=prof_shen(3) + (t2-t1)
     else
        write(*,*) 'error in multevenodd'
        stop
     endif

     call cpu_time(t1)
     call cublas_dcopy(mye*(my-2),dp_tmpee,1,devPtrF,2) 
     call cublas_dcopy((myo+1)*(my-2),dp_tmpoo,1,devPtrF+nk,2) 
     call cpu_time(t2); prof_shen(4)=prof_shen(4) + (t2-t1)
!
     call cpu_time(t1)
     !write(*,*) 'cublas_get_matrix'
     !cublas_get_matrix F <- device (very slow)
     if (iopt.lt.0) then 
        istat= cublas_get_matrix(my-2+1,mz-2+1,nk,devPtrF,my-2+1,F,my-2+1) ! device -> cpu-memory 
     end if
     call cpu_time(t2); prof_shen(5)=prof_shen(5) + (t2-t1)

  elseif (bound.eq.'n') then

     !cublas_set_matrix F-> device
     istat= cublas_set_matrix(my-2,my-2,nk,F,my-2,devPtrF,my-2)
     !write(*,*) 'cublas alloc setmat F',istat
     do j=1,my-2
        !call cublas_dcopy(mye,F(1,j),2,tmpe(1,j),1)
        call cublas_dcopy(mye,devPtrF+(j-1)*(my-2)*nk,2,dp_tmpe+(j-1)*mye*nk,1)
     enddo
     do j=1,my-2
        !call cublas_dcopy(myo,F(2,j),2,tmpo(1,j),1)
        call cublas_dcopy(myo,devPtrF+((j-1)*(my-2)+1)*nk,2,dp_tmpo+(j-1)*myo*nk,1)
     enddo

     write(*,*) 'Neumann option is not implemented yet.'
     stop

     call dcopy(my-2,F(1,1),my-2,bcon,1)
     ! /* C(1:2:n,1:n)=AE*B(2:2:n,1:n); */
     ! /* C(2:2:n,1:n)=AO*B(1:2:n,1:n); */
     if (iopt.gt.0) then
        fac=0.5d0
        call dgemm('n','n',mye,my-2,myo,1.d0,PSM1EN,mye, &
             &     tmpo,myo,0.d0,tmpee,mye)
        call dgemm('n','n',myo,my-2,mye,1.d0,PSM1ON,myo, &
             &     tmpe,mye,0.d0,tmpoo,myo)
     elseif (iopt.lt.0) then
        fac=1.d0
        call dgemm('n','n',mye,my-2,myo,1.d0,PEN,mye, &
             &     tmpo,myo,0.d0,tmpee,mye)
        call dgemm('n','n',myo,my-2,mye,1.d0,PON,myo, &
             &     tmpe,mye,0.d0,tmpoo,myo)         
     else
        write(*,*) 'error in multevenodd '
        stop
     endif
     ! /* reconstruct F(m,m) */
     do j=1,my-2
        call dcopy(mye,tmpee(1,j),1,F(1,j),2)
     enddo
     do j=1,my-2
        call dcopy(myo,tmpoo(1,j),1,F(2,j),2)
     enddo
     ! /* repair constant mode of shen base spectra */
     ! /* only for neumann case */
     ! /* C(1,1:n)=fac*B(1,1:n); */
     call daxpy(my-2,fac,bcon,1,F(1,1),my-2)
  endif
  return
end subroutine multevenodd_gpu

!----|----------------------------------------------------------------|
subroutine solV(G,alpha,bound)
  use legchebshen
  use ctes,only:nk,my,mz,mye,myo
  use timing
  implicit none
  real*8 G(my-2+1,my-2+1) ! wrapped version 
  character*1 bound
  integer m, mm2, ip
  real*8 alpha,rp,ra1

  !call cpu_time(t1)
  !G = transpose(G)
  !call cpu_time(t2); prof_solV(1)=prof_solV(1) + (t2-t1)

  m=my
  mm2=my-2

  if (bound.eq.'d') then
     ! even
     call cpu_time(t1)
     call dcopy(mye,      SED,1,DDE,1)   ! DDE <- SE
     call daxpy(mye,alpha,MED,1,DDE,1)   ! DDE <- alpha*ME+SE
     !DDE(:)=alpha*MED(:) + SED(:)
     call cpu_time(t2); prof_solV(2)=prof_solV(2) + (t2-t1)
     ! odd
     call cpu_time(t1)
     call dcopy(myo,      SOD,1,DDO,1)   ! DDO <- SO
     call daxpy(myo,alpha,MDO,1,DDO,1)   ! DDO <- alpha*MO+SO
     !DDO(:) = alpha*MDO(:) + SOD(:)
     call cpu_time(t2); prof_solV(3)=prof_solV(3) + (t2-t1)
     do ip=1,mm2
        rp=ED(ip) ! 
        ra1=(rp*alpha+1.d0)
        !write(*,*) 'rp, ra1', rp, ra1
        ! /* gp(1:mm2)=G(ip,1:mm2) */
        call cpu_time(t1)
        call dcopy(mm2,G(ip,1),mm2+1,gp,1) ! skip=my-2  --> G wrapped with my-1
        !call dcopy(mm2+1,G(1,ip),1,gp,1) ! G transposed
        call cpu_time(t2); prof_solV(4)=prof_solV(4) + (t2-t1)

        ! /* DDEp=rp*(alpha*MED+SED)+MED */
        call cpu_time(t1)
        call dcopy(mye,MED,1,DDEp,1)
        call daxpy(mye,rp,DDE,1,DDEp,1)
        !DDEp(:) = rp*DDE(:) + MED(:)
        call cpu_time(t2); prof_solV(5)=prof_solV(5) + (t2-t1)

        ! /* DDOp=rp*(alpha*MDO+SOD)+MDO */
        call cpu_time(t1)
        call dcopy(myo,MDO,1,DDOp,1)
        call daxpy(myo,rp,DDO,1,DDOp,1)
        !DDOp(:) = rp*DDO(:) + MDO(:)
        call cpu_time(t2); prof_solV(6)=prof_solV(6) + (t2-t1)

        ! /* DPE=(rp*alpha+1)*MPPED */
        call cpu_time(t1)
        call dcopy(mye-1,MPPED,1,DPEp,1)
        call dscal(mye-1,ra1,DPEp,1)
        !DPEp(:) = ra1*MPPED(:)
        call cpu_time(t2); prof_solV(7)=prof_solV(7) + (t2-t1)

        ! /* DPO=(rp*alpha+1)*MPPOD */
        call cpu_time(t1)
        call dcopy(myo-1,MPPOD,1,DPOp,1)
        call dscal(myo-1,ra1,DPOp,1)
        !DPOp(:) = ra1*MPPOD(:)
        call cpu_time(t2); prof_solV(8)=prof_solV(8) + (t2-t1)

        ! !gp --> vp
        call cpu_time(t1)
        call solevenodd(gp,DDEp,DDOp,DPEp,DPOp,vp,mye,myo)
        call cpu_time(t2); prof_solV(9)=prof_solV(9) + (t2-t1)

        ! /* G(ip,1:mm2)=vp(1:mm2) */
        call cpu_time(t1)
        call dcopy(mm2,vp,1,G(ip,1),mm2+1) ! skip=my-2  --> G with wrapping my-1
        !call dcopy(mm2+1,vp,1,G(1,ip),1)
        call cpu_time(t2); prof_solV(10)=prof_solV(10) + (t2-t1)
     enddo
     !call cpu_time(t1)
     !G = transpose(G)
     !call cpu_time(t2); prof_solV(1)=prof_solV(1) + (t2-t1)
  elseif (bound.eq.'n') then
     ip=1
     ! /* gp(1:mm2)=G(ip,1:mm2) */
     call dcopy(mm2,G(ip,1),mm2+1,gp,1) ! skip=my-2 ! --> wrapped version
     ! DDE=alpha*MEN+SEN
     call dcopy(mye,      SEN,1,DDE,1)   ! DDE <-- SE
     call daxpy(mye,alpha,MEN,1,DDE,1)   ! DDE <-- alpha*ME+SE
     ! DDO=alpha*MON+SON
     call dcopy(myo,      SON,1,DDO,1)   ! DDO <-- SO
     call daxpy(myo,alpha,MON,1,DDO,1)   ! DDO <-- alpha*MO+SO
     DDEp=DDE; DDOp=DDO; ! because solevenodd does not preserve DDE,DDO...
     ! DPE=alpha*MPPEN
     call dcopy(mye-1,MPPEN,1,DPEp,1)
     call dscal(mye-1,alpha,DPEp,1)
     ! DPO=alpha*MPPON
     call dcopy(myo-1,MPPON,1,DPOp,1)
     call dscal(myo-1,alpha,DPOp,1)
     ! G(ip,1:mm2)=vp(1:mm2)
     if (alpha.lt.1.d-10) then
        DDE(1)=1.d0
        gp(1)=0.d0
     endif
     call solevenodd(gp,DDEp,DDOp,DPEp,DPOp,vp,mye,myo)
     if (alpha.lt.1.d-10) then
        ! fix DDE(1)   ! for zero-mode           
        DDE(1)=0.d0
     endif
     call dcopy(mm2,vp,1,G(ip,1),mm2+1)
     !
     do ip=2,mm2
        rp=EN(ip)
        ra1=(rp*alpha+1.d0)
        !  gp(1:mm2)=G(ip,1:mm2)
        call dcopy(mm2,G(ip,1),mm2+1,gp,1) ! skip=my-2 --> G wrapped with my-1
        !  DDEp=rp*(alpha*MEN+SEN)+MEN
        call dcopy(mye,MEN,1,DDEp,1)
        call daxpy(mye,rp,DDE,1,DDEp,1)
        !  DDOp=rp*(alpha*MON+SON)+MON
        call dcopy(myo,MON,1,DDOp,1)
        call daxpy(myo,rp,DDO,1,DDOp,1)
        !  DPEp=(rp*alpha+1)*MPPEN
        call dcopy(mye-1,MPPEN,1,DPEp,1)
        call dscal(mye-1,ra1,DPEp,1)
        !  DPOp=(rp*alpha+1)*MPPON
        call dcopy(myo-1,MPPON,1,DPOp,1)
        call dscal(myo-1,ra1,DPOp,1)
        !gp --> vp            
        call solevenodd(gp,DDEp,DDOp,DPEp,DPOp,vp,mye,myo)
        ! G(ip,1:mm2)=vp(1:mm2)
        call dcopy(mm2,vp,1,G(ip,1),mm2+1) ! skip=my-2  --> G with wrapping my-1
     enddo
      
  endif
  return
end subroutine solV
!----|----------------------------------------------------------------|
subroutine solevenodd(f,DDE,DDO,DPE,DPO,u,me,mo)
!
!  do not use legchebshen
!  use legchebshen, only: DPEp2, DPOp2
  use legchebshen, only: fe,fo
  implicit none
  integer me,mo,info
  real*8 DDE(me),DDO(mo),DPE(me-1),DPO(mo-1)
  real*8 f(me+mo) !,fe(me),fo(mo)
  real*8 u(me+mo) !,xe(me),xo(mo)
  ! /* fe=f(1:2:m); */
  call dcopy(me,f(1),2,fe,1)
  ! /* xe=symtri(DE,DPE,fe); */
  !call symtri(DDE,DPE,fe,xe,me) 

  call DPTSV(me,1,DDE,DPE,fe,me,info)
  !DDE (N), DPE (N-1)  
  !
  !call dcopy(me-1,DPE,1,DPEp2,1)
  !call DGTSV(me,1,DPE,DDE,DPEp2,fe,me,info)

  ! /* u(1:2:m)=xe'; */
  !call dcopy(me,xe,1,u(1),2)
  call dcopy(me,fe,1,u(1),2)

  ! /* fo=f(2:2:m); */
  call dcopy(mo,f(2),2,fo,1)
  ! /* xo=symtri(DO,DPO,fo); */
  !call symtri(DDO,DPO,fo,xo,mo)
  !
  call DPTSV(mo,1,DDO,DPO,fo,mo,info)
  !
  !call dcopy(mo-1,DPO,1,DPOp2,1)
  !call DGTSV(mo,1,DPO,DDO,DPOp2,fo,mo,info)
  !
  ! /* u(2:2:m)=xo'; */
  !call dcopy(mo,xo,1,u(2),2)
  call dcopy(mo,fo,1,u(2),2)

  return
end subroutine solevenodd
!----|----------------------------------------------------------------|
subroutine symtri(b,c,r,u,n)
  
  implicit none
  ! /****************************************************************/
  ! /*  solve Au=r with A symm. tridiagonal                         */
  ! /*  input                                                       */
  ! /*  b(n)  : the diagonal of A                                   */
  ! /*  c(n-1): the upper diagonal of A                             */
  ! /*  r(n)  : r.h.s. of the system                                */
  ! /*                                                              */
  ! /*  output                                                      */
  ! /*   u(n)                                                       */
  ! /****************************************************************/
  ! this can be replaced by DPTSV in lapack.
  ! do not use DGTSV, because it modifies both DDE DPE, which affects the result.
  integer j,n
  real*8 b(n),c(n-1),r(n),u(n)
  real*8 gam(n)
  real*8 bet

  bet=b(1)
  u(1)=r(1)/bet
  do j=2,n
     gam(j)=c(j-1)/bet
     bet=b(j)-c(j-1)*gam(j)
     u(j)=(r(j)-c(j-1)*u(j-1))/bet
  enddo
  do j=n-1,1,-1
     u(j)=u(j)-gam(j+1)*u(j+1)
  enddo
  return
end subroutine symtri

!----|----------------------------------------------------------------|
subroutine solV_gpu(G,alpha,bound)
  use legchebshen
  use ctes,only:nk,my,mz,mye,myo
  use timing
  use gpudevice

  implicit none
  real*8 G(my-2+1,my-2+1) ! wrapped version 
  character*1 bound
  integer m, mm2, ip, istat, lsto
  real*8 alpha,rp,ra1
  
  integer istatus, version

  m=my
  mm2=my-2
  
  lsto=(myo+1)

  if (bound.eq.'d') then
     ! even
     !call cpu_time(t1)
     !cublas_set_matrix F-> device 
     !!!istat= cublas_set_matrix(my-1,my-1,nk,G,my-1,devPtrG,my-1)
     !call cpu_time(t2); prof_solV(1)=prof_solV(1) + (t2-t1)     
     !
     call cpu_time(t1)
     !DDE(:)=alpha*MED(:) + SED(:)
     call cublas_dcopy(mye,      dp_SED,1,dp_DDE,1) ! DDE <- SE
     call cublas_daxpy(mye,alpha,dp_MED,1,dp_DDE,1) ! DDE <- alpha*ME+SE
     !
     call cpu_time(t2); prof_solV(2)=prof_solV(2) + (t2-t1)
     ! odd
     call cpu_time(t1)
     !DDO(:) = alpha*MDO(:) + SOD(:)
     call cublas_dcopy(myo,      dp_SOD,1,dp_DDO,1) ! DDO <- SO
     call cublas_daxpy(myo,alpha,dp_MDO,1,dp_DDO,1) ! DDO <- alpha*MO+SO
     call cpu_time(t2); prof_solV(3)=prof_solV(3) + (t2-t1)

     call cpu_time(t1)
     ! transpose F --> G
#ifdef TRANS_GPU 
     istat = transpose_gpu(handle,int(mm2+1,kind=8),int(mm2+1,kind=8),1.d0,devPtrF,devPtrG) ! V^F --> G 
#else
     do ip=1,mm2
        call cublas_dcopy(mm2,devPtrF+(ip-1)*nk,mm2+1,devPtrG+(mm2+1)*(ip-1)*nk    ,1) ! slow
     end do
#endif /* TRANS_GPU */
     call cpu_time(t2); prof_solV(1)=prof_solV(1) + (t2-t1) 

     call cpu_time(t1)
     ! quick splitting G --> fe, fo
     call cublas_dcopy(mye*mm2,devPtrG,2,dp_tmpe,1)
     call cublas_dcopy(lsto*mm2,devPtrG+nk,2,dp_tmpo,1)
     call cpu_time(t2); prof_solV(4)=prof_solV(4) + (t2-t1)

     ! /* DDEp=rp*(alpha*MED+SED)+MED */
     call cpu_time(t1) 
!     call cublas_dcopy(mye*mm2,   dp_MED,1,dp_DDEp,1) 
     do ip=1,mm2
        istat = cublas_set_stream(handle, streams(ip)) ! not faster
        call cublas_dcopy(mye,   dp_MED,1,dp_DDEp+(mye)*(ip-1)*nk,1) ! slow
     end do 
     do ip=1,mm2
        ! Set CUDA stream
        istat = cublas_set_stream(handle, streams(ip))
        rp=ED(ip) ! 
        !DDEp(:) = rp*DDE(:) + MED(:)
        call cublas_daxpy(mye,rp,dp_DDE,1,dp_DDEp+(mye)*(ip-1)*nk,1)       
     end do
     call cpu_time(t2); prof_solV(5)=prof_solV(5) + (t2-t1)

     call cpu_time(t1)
     do ip=1,mm2
        ! /* DPE=(rp*alpha+1)*MPPED */
        call cublas_dcopy(mye-1,    dp_MPPED  ,1,dp_DPEp+nk+(mye)*(ip-1)*nk,1)
     end do
     do ip=1,mm2
        ra1=(ED(ip)*alpha+1.d0)        
        !DPEp(:) = ra1*MPPED(:)
        call cublas_dscal(mye-1,ra1,dp_DPEp+nk+(mye)*(ip-1)*nk,1)
     end do
     call cublas_dcopy(mye*mm2-1,    dp_DPEp+nk,1,dp_DPEp2,1)
     call cpu_time(t2); prof_solV(7)=prof_solV(7) + (t2-t1)

     call cpu_time(t1)
     istatus = cusparse_dgtsv_stridedbatch(handle,mye,dp_DPEp,dp_DDEp,dp_DPEp2,dp_tmpe,mm2,mye)
     call cpu_time(t2); prof_solV(9)=prof_solV(9) + (t2-t1)     

     call cpu_time(t1)
     do ip=1,mm2
        ! /* DDOp=rp*(alpha*MDO+SOD)+MDO */
        call cublas_dcopy(myo,   dp_MDO,1,dp_DDOp+lsto*(ip-1)*nk,1) ! slow
     end do
     do ip=1,mm2
        rp=ED(ip) ! 
        !DDOp(:) = rp*DDO(:) + MDO(:)
        call cublas_daxpy(myo,rp,dp_DDO,1,dp_DDOp+lsto*(ip-1)*nk,1) ! slow
     end do
     call cpu_time(t2); prof_solV(6)=prof_solV(6) + (t2-t1)

     call cpu_time(t1)
     do ip=1,mm2
        ! /* DPO=(rp*alpha+1)*MPPOD */
        call cublas_dcopy(myo-1,    dp_MPPOD  ,1,dp_DPOp+nk+lsto*(ip-1)*nk,1) ! slow
     end do
     do ip=1,mm2
        ra1=(ED(ip)*alpha+1.d0)     
        !DPOp(:) = ra1*MPPOD(:)
        call cublas_dscal(myo-1,ra1,dp_DPOp+nk+lsto*(ip-1)*nk,1) ! slow
     end do
     call cublas_dcopy((lsto)*mm2-1,dp_DPOp+nk,1,dp_DPOp2,1)
     call cpu_time(t2); prof_solV(8)=prof_solV(8) + (t2-t1)

     call cpu_time(t1)
     istatus = cusparse_dgtsv_stridedbatch(handle,myo,dp_DPOp,dp_DDOp,dp_DPOp2,dp_tmpo,mm2,lsto)
     call cpu_time(t2); prof_solV(9)=prof_solV(9) + (t2-t1)     

     call cpu_time(t1)
     call cublas_dcopy((mye)*mm2,dp_tmpe,1,devPtrG,2)
     call cublas_dcopy(lsto*mm2,dp_tmpo,1,devPtrG+nk,2)
     call cpu_time(t2); prof_solV(10)=prof_solV(10) + (t2-t1)

     call cpu_time(t1)
     ! transpose G --> F
#ifdef TRANS_GPU 
     istat = transpose_gpu(handle,int(mm2+1,kind=8),int(mm2+1,kind=8),1.d0,devPtrG,devPtrF) ! G^T --> F 
#else
     do ip=1,mm2
        ! /* G(ip,1:mm2)=vp(1:mm2) */
        call cublas_dcopy(mm2,devPtrG+(mm2+1)*(ip-1)*nk,1,devPtrF + (ip-1)*nk,mm2+1) 
        ! skip=my-2 --> G with wrapping my-1
     enddo
#endif
     call cpu_time(t2); prof_solV(11)=prof_solV(11) + (t2-t1)
     !
  elseif (bound.eq.'n') then
     ip=1
     write(*,*) 'not implemented solV_gpu'
     stop     
!     ! /* gp(1:mm2)=G(ip,1:mm2) */
!     call dcopy(mm2,G(ip,1),mm2,gp,1) ! skip=my-2
!     ! DDE=alpha*MEN+SEN
!     call dcopy(mye,      SEN,1,DDE,1)   ! DDE <-- SE
!     call daxpy(mye,alpha,MEN,1,DDE,1)   ! DDE <-- alpha*ME+SE
!     ! DDO=alpha*MON+SON
!     call dcopy(myo,      SON,1,DDO,1)   ! DDO <-- SO
!     call daxpy(myo,alpha,MON,1,DDO,1)   ! DDO <-- alpha*MO+SO
!     ! DPE=alpha*MPPEN
!     call dcopy(mye-1,MPPEN,1,DPEp,1)
!     call dscal(mye-1,alpha,DPEp,1)
!     ! DPO=alpha*MPPON
!     call dcopy(myo-1,MPPON,1,DPOp,1)
!     call dscal(myo-1,alpha,DPOp,1)
!     ! G(ip,1:mm2)=vp(1:mm2)
!     if (alpha.lt.1.d-10) then
!        DDE(1)=1.d0
!        gp(1)=0.d0
!     endif
!     call solevenodd(gp,DDE,DDO,DPEp,DPOp,vp,mye,myo)
!     if (alpha.lt.1.d-10) then
!        ! fix DDE(1)              
!        DDE(1)=0.d0
!     endif
!     call dcopy(mm2,vp,1,G(ip,1),mm2)
!     !
!     do ip=2,mm2
!        rp=EN(ip)
!        ra1=(rp*alpha+1.d0)
!        !  gp(1:mm2)=G(ip,1:mm2)
!        call dcopy(mm2,G(ip,1),mm2,gp,1) ! skip=my-2
!        !  DDEp=rp*(alpha*MEN+SEN)+MEN
!        call dcopy(mye,MEN,1,DDEp,1)
!        call daxpy(mye,rp,DDE,1,DDEp,1)
!        !  DDOp=rp*(alpha*MON+SON)+MON
!        call dcopy(myo,MON,1,DDOp,1)
!        call daxpy(myo,rp,DDO,1,DDOp,1)
!        !  DPEp=(rp*alpha+1)*MPPEN
!        call dcopy(mye-1,MPPEN,1,DPEp,1)
!        call dscal(mye-1,ra1,DPEp,1)
!        !  DPOp=(rp*alpha+1)*MPPON
!        call dcopy(myo-1,MPPON,1,DPOp,1)
!        call dscal(myo-1,ra1,DPOp,1)
!        !gp --> vp            
!        call solevenodd(gp,DDEp,DDOp,DPEp,DPOp,vp,mye,myo)
!        ! G(ip,1:mm2)=vp(1:mm2)
!        call dcopy(mm2,vp,1,G(ip,1),mm2)
!     enddo
      
  endif
  return
end subroutine solV_gpu
!----|----------------------------------------------------------------|
subroutine solevenodd_gpu(dp_f,dp_u,me,mo)
  !
  ! 1-d tri-diagonal solver for even-odd decomposed vector (old)
  ! with 1d-vector format matrix DPE DDE DPE2, ... 
  ! dp_DPEp,dp_DDEp,dp_DPEp2 used to be 1d-vector of the size of 'me', 'me-1'
  ! 
  use ctes, only:nk
  use gpudevice 
  !  do not use legchebshen here
  implicit none
  integer me,mo,info
  integer istatus, version

  integer*8 dp_f,dp_u

  ! /* fe=f(1:2:m); */
  call cublas_dcopy(me,dp_f,2,dp_fe,1)
  ! /* xe=symtri(DE,DPE,fe); */
  istatus = cusparse_dgtsv(handle,me,1,dp_DPEp,dp_DDEp,dp_DPEp2,dp_fe,me)
  ! /* u(1:2:m)=xe'; */
  call cublas_dcopy(me,dp_fe,1,dp_u,2)

  ! /* fo=f(2:2:m); */
  call cublas_dcopy(mo,dp_f+nk,2,dp_fo,1)
  ! /* xo=symtri(DO,DPO,fo); */
  istatus = cusparse_dgtsv(handle,mo,1,dp_DPOp,dp_DDOp,dp_DPOp2,dp_fo,mo)
  ! /* u(2:2:m)=xo'; */
  call cublas_dcopy(mo,dp_fo,1,dp_u+nk,2)

  return
end subroutine solevenodd_gpu
#endif /* HELM_SHEN */
