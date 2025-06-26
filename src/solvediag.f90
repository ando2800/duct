#ifndef RECTANGULAR
subroutine solvediag(U,F,FT,bound,xlambda)
  use ctes,only:nk,my,mz,ngpu
  use fulldiagneu
  use fulldiagdir
  use chebystuff
  use timing
#ifdef CUBLAS
  use gpudevice
#endif
  implicit none 
  !     
  !     solve the problem
  !
  real(nk) U(my,mz),F(my-2+1,mz-2+1)
  real(nk) FT(my-2+1,mz-2+1)
  real(nk) ff(4),uu(4)
  character bound,poisson
  integer m,n,mm2,i,j,k,l,istat
  real(nk) bcontv,bconth,xlambda,den,val,valore,det,sumb,sumbl,sumbr,suml,sumr,sumt,sumtl,sumtr

  m=my
  n=m-1
  mm2=m-2

  F=0.d0
  FT=0.d0
  tmprt=0.d0

  do j=2,m-1
     do i=2,m-1
        F(i-1,j-1)=U(i,j)
     enddo
  enddo

  if(bound.eq.'d') then

     !     correct rhs for boundary condition     
     do j=1,m-2
        do i=1,m-2
           bcontv=D2bx(i+1,1)*U(1,j+1)+D2bx(i+1,2)*U(my,j+1)
           bconth=D2by(j+1,1)*U(i+1,1)+D2by(j+1,2)*U(i+1,my)
           F(i,j)=F(i,j)-bcontv-bconth
        enddo
     enddo

     if (ngpu.eq.0) then
        ! Projection of rhs
        ! F=inv(P)*F*inv(P');
        !
        call dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,PDM1T,mm2+1,F, &
             &             mm2+1,0.d0,FT,mm2+1)  
        ! transpose
        do i=1,m-2+1
           do j=1,m-2+1
              tmprt(i,j)=FT(j,i)
           enddo
        enddo
        call dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,tmprt,mm2+1,PDM1T, &
             &             mm2+1,0.d0,F,mm2+1)   
        do j=1,m-2
           do i=1,m-2
              !den=(ED(i)+ED(j)-xlambda)
              FT(i,j)=F(i,j)/(ED2(i,j)-xlambda)
           enddo
        enddo
        !     
        ! reconstruct physical sol
        ! FT=P*FT*P';
        !
        call dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,PDT,mm2+1,FT, &
             &             mm2+1,0.d0,F,mm2+1) !fastest option
        ! transpose
        do i=1,m-2+1
           do j=1,m-2+1
              tmprt(i,j)=F(j,i)
           enddo
        enddo
        call dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,tmprt,mm2+1,PDT, &
             &             mm2+1,0.d0,FT,mm2+1)
        
     else
#ifdef CUBLAS
        ! GPU version
        ! Projection of rhs
        ! F=inv(P)*F*inv(P')
        call cpu_time(t1)
        istat= cublas_set_matrix(mm2+1,mm2+1,nk,F,mm2+1,devPtrF,mm2+1)
        call cpu_time(t2); prof_diag(1)=prof_diag(1) + (t2-t1)
        call cpu_time(t1)
        call cublas_dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,dp_PDM1T,mm2+1,devPtrF, & 
             &                    mm2+1,0.d0,devPtrFT,mm2+1)
        call cpu_time(t2); prof_diag(2)=prof_diag(2) + (t2-t1)
        call cpu_time(t1)
        !! transpose
        istat = transpose_gpu(handle,int(mm2+1,kind=8),int(mm2+1,kind=8),1.d0, & 
             &                devPtrFT,dp_tmp2d) ! the size should be even??
        call cpu_time(t2); prof_diag(3)=prof_diag(3) + (t2-t1)
        
        call cpu_time(t1)
        call cublas_dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,dp_tmp2d,mm2+1,dp_PDM1T, & 
             &                    mm2+1,0.d0,devPtrF,mm2+1)
        call cpu_time(t2); prof_diag(4)=prof_diag(4) + (t2-t1)
        call cpu_time(t1)
        istat = mat_invdiag_gpu(handle,int(mm2+1,kind=8),int(mm2+1,kind=8),xlambda, & 
             &                  devPtrF,dp_ED2,devPtrFT) !  FT(i,j)=F(i,j)/(ED2(i,j)-xlambda)
        call cpu_time(t2); prof_diag(5)=prof_diag(5) + (t2-t1)
        !
        ! reconstruct physical sol
        ! FT=P*FT*P'
        !
        call cpu_time(t1)
        call cublas_dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,dp_PDT,mm2+1,devPtrFT, &
             &                    mm2+1,0.d0,devPtrF,mm2+1) 
        call cpu_time(t2); prof_diag(6)=prof_diag(6) + (t2-t1)
        call cpu_time(t1)
        !! transpose
        istat = transpose_gpu(handle,int(mm2+1,kind=8),int(mm2+1,kind=8),1.d0, & 
             &                devPtrF,dp_tmp2d) ! F --> tmp 
        call cpu_time(t2); prof_diag(7)=prof_diag(7) + (t2-t1)
        call cpu_time(t1)
        call cublas_dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,dp_tmp2d,mm2+1,dp_PDT, &
             &                    mm2+1,0.d0,devPtrFT,mm2+1)
        call cpu_time(t2); prof_diag(8)=prof_diag(8) + (t2-t1)
        call cpu_time(t1)
        istat = cublas_get_matrix(mm2+1,mm2+1,nk,devPtrFT,mm2+1,FT,mm2+1)
        call cpu_time(t2); prof_diag(9)=prof_diag(9) + (t2-t1)
#endif /* CUBLAS */
     end if
     ! fill back in the solution
     do j=1,m-2
        do i=1,m-2
           U(i+1,j+1)=FT(i,j)
        enddo
     enddo
        
  else if(bound.eq.'n') then

     if ((ngpu.eq.0).or.(dabs(xlambda).lt.1.d-10)) then
        call dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,PNM1T,mm2+1,F, &
             &             mm2+1,0.d0,FT,mm2+1) ! fastest option
        ! transpose
        do i=1,m-2
           do j=1,m-2
              tmprt(i,j)=FT(j,i)
           enddo
        enddo
        call dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,tmprt,mm2+1,PNM1T, &
             &             mm2+1,0.d0,F,mm2+1)

        if(dabs(xlambda).gt.1.d-10) then
           
           !--- CASE NEUMANN LAMBDA=0
           do j=1,m-2
              do i=1,m-2
                 den=(EN(i)+EN(j)-xlambda)
                 FT(i,j)=F(i,j)/den
              enddo
           enddo
           
        else
            
           !--- CASE NEUMANN LAMBDA=0
           
           ! do till the critical column
           do j=1,iofzero-1
              do i=1,m-2
                 den=(EN(i)+EN(j)-xlambda)
                 FT(i,j)=F(i,j)/den
              enddo
           enddo
           
           ! do the iofzero column
           j=iofzero
           do i=1,iofzero-1
              den=(EN(i)+EN(j)-xlambda)
              FT(i,j)=F(i,j)/den
           enddo
           
           FT(iofzero,iofzero)=0.d0
           
           do i=iofzero+1,m-2
              den=(EN(i)+EN(j)-xlambda)
              FT(i,j)=F(i,j)/den
           enddo
           
           ! complete the other values
           do j=iofzero+1,m-2
              do i=1,m-2
                 den=(EN(i)+EN(j)-xlambda)
                 FT(i,j)=F(i,j)/den
              enddo
           enddo
           
           ! procedure to fix the solution in (iofzero,iofzero)
           valore=0.d0
           do l=1,m-2
              val=0.d0
              do k=1,m-2
                 val=val+PN(iofzero,k)*FT(l,k)
              enddo
              valore=valore+PN(iofzero,l)*val
           enddo
           FT(iofzero,iofzero)= &
                &           (wanted-valore)/(PN(iofzero,iofzero)**2)         
        endif
         
        !     solution in physical space
        
        call dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,PNT,mm2+1,FT, &
             &             mm2+1,0.d0,F,mm2+1)   !fastest option
        ! transpose
        do i=1,m-2
           do j=1,m-2
              tmprt(i,j)=F(j,i)
           enddo
        enddo
        call dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,tmprt,mm2+1,PNT, &
             &             mm2+1,0.d0,FT,mm2+1)
        
     else
#ifdef CUBLAS
        istat= cublas_set_matrix(mm2+1,mm2+1,nk,F,mm2+1,devPtrF,mm2+1)

        call cublas_dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,dp_PNM1T,mm2+1,devPtrF, & 
             &                    mm2+1,0.d0,devPtrFT,mm2+1)  
        istat = transpose_gpu(handle,int(mm2+1,kind=8),int(mm2+1,kind=8),1.d0, & 
             &                devPtrFT,dp_tmp2d) ! the size should be even??
        call cublas_dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,dp_tmp2d,mm2+1,dp_PNM1T, & 
             &                    mm2+1,0.d0,devPtrF,mm2+1)   
        istat = mat_invdiag_gpu(handle,int(mm2+1,kind=8),int(mm2+1,kind=8),xlambda, & 
             &                  devPtrF,dp_EN2,devPtrFT) !  FT(i,j)=F(i,j)/(ED2(i,j)-xlambda) 
        call cublas_dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,dp_PNT,mm2+1,devPtrFT, &
             &                    mm2+1,0.d0,devPtrF,mm2+1) 
        istat = transpose_gpu(handle,int(mm2+1,kind=8),int(mm2+1,kind=8),1.d0, & 
             &                devPtrF,dp_tmp2d) ! F --> tmp 
        call cublas_dgemm('t','n',mm2+1,mm2+1,mm2+1,1.d0,dp_tmp2d,mm2+1,dp_PNT, &
             &                    mm2+1,0.d0,devPtrFT,mm2+1) 
        istat = cublas_get_matrix(mm2+1,mm2+1,nk,devPtrFT,mm2+1,FT,mm2+1)
#endif /* CUBLAS */
     end if

     !  fill back in the solution         
     do j=1,m-2
        do i=1,m-2
           U(i+1,j+1)=FT(i,j)
        enddo
     enddo

     ! Fix Neumann values on the border
     !     bottom/top lines

     det=DC(1,1)*DC(m,m)-DC(1,m)*DC(m,1)
     do i=2,m-1
        sumb=0.d0
        sumt=0.d0
        do j=1,m-2
           sumb=sumb+DC(1,j+1)*FT(i-1,j)
           sumt=sumt+DC(m,j+1)*FT(i-1,j)
        enddo
        U(i,1)=(-DC(m,m)*sumb+DC(1,m)*sumt)/det
        U(i,m)=(DC(m,1)*sumb-DC(1,1)*sumt)/det
     enddo
     
     !  right/left lines
     
     do j=2,m-1
        suml=0.d0
        sumr=0.d0
        do i=1,m-2
           suml=suml+DC(1,i+1)*FT(i,j-1)
           sumr=sumr+DC(m,i+1)*FT(i,j-1)
        enddo
        U(1,j)=(-DC(m,m)*suml+DC(1,m)*sumr)/det
        U(m,j)=(DC(m,1)*suml-DC(1,1)*sumr)/det
     enddo
     
     !     corners 
     
     sumbl=0.d0
     do i=2,m-1
        sumbl=sumbl-DC(1,i)*(U(i,1)+U(1,i))
     enddo
     ff(1)=sumbl
     
     sumbr=0.d0
     do i=2,m-1
        sumbr=sumbr-DC(1,i)*U(m,i)+DC(m,i)*U(i,1)
     enddo
     ff(2)=sumbr
     
     sumtr=0.d0
     do i=2,m-1
        sumtr=sumtr-DC(m,i)*(U(i,m)+U(m,i))
     enddo
     ff(3)=sumtr
     
     sumtl=0.d0
     do i=2,m-1
        sumtl=sumtl-DC(1,i)*U(i,m)+DC(m,i)*U(1,i)
     enddo
     ff(4)=sumtl
     
     U(1,1)=0.d0
     do i=1,4
        U(1,1)=U(1,1)+AM1(1,i)*ff(i)
     enddo
     
     U(m,1)=0.d0
     do i=1,4
        U(m,1)=U(m,1)+AM1(2,i)*ff(i)
     enddo
     
     U(m,m)=0.d0
     do i=1,4
        U(m,m)=U(m,m)+AM1(3,i)*ff(i)
     enddo
     
     U(1,m)=0.d0
     do i=1,4
        U(1,m)=U(1,m)+AM1(4,i)*ff(i)
     enddo
     
  else
     write(6,*) 'in solvediag error in BC'
     stop
  endif
  
  return
end subroutine solvediag
#else /* RECTANGULAR */
subroutine solvediag_dummy()
end subroutine solvediag_dummy
#endif /* not RECTANGULAR */
