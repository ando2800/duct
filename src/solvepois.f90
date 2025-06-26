subroutine solvepois(U,F,FT,wk1,wk2,xlambda,m)
  use ctes,only:nk,my
  use fulldiagneu
  !
  implicit none 
  !
  real(nk) U(my,*),F(my-2,*)
  real(nk) FT(my-2,*)
  real(nk) wk1(my-2,*),wk2(my-2,*)
  real(nk) xlambda,den, val, valore,c,coei,s,snei,snej,snwi,snwj,ssei,ssej,sswi,sswj
  real(nk) valuene,valuenw,valuese,valuesw,vintne,vintnw,vintse,vintsw,vr,vl,vpol
  integer m,n,mm2,i,j,ii,jj,k,l

  n=m-1
  mm2=m-2
  do j=2,m-1
     do i=2,m-1
        F(i-1,j-1)=U(i,j)
     enddo
  enddo
  
  call dgemm('n','n',mm2,mm2,mm2,1.d0,PNM1,mm2,F,mm2,0.d0,FT,mm2)
  call dgemm('n','t',mm2,mm2,mm2,1.d0,FT,mm2,PNM1,mm2,0.d0,F,mm2)
  
  if(dabs(xlambda).gt.1.d-10) then

     ! CASE NEUMANN LAMBDA=0
     do j=1,m-2
        do i=1,m-2
           den=(EN(i)+EN(j)-xlambda)
           FT(i,j)=F(i,j)/den
        enddo
     enddo

  else

     ! CASE NEUMANN LAMBDA=0

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
     FT(iofzero,iofzero)=(wanted-valore)/(PN(iofzero,iofzero)**2)
     !         write(6,*) 'wanted=',wanted
     ! end procedure to fix the solution in (iofzero,iofzero)
     
  endif
  
  ! solution in physical space
  
  call dgemm('n','n',mm2,mm2,mm2,1.d0,PN,mm2,FT,mm2,0.d0,F,mm2)
  call dgemm('n','t',mm2,mm2,mm2,1.d0,F,mm2,PN,mm2,0.d0,FT,mm2)

  ! fill in the solution

  do j=1,m-2
     do i=1,m-2
        U(i+1,j+1)=FT(i,j)
     enddo
  enddo
  !write(6,*) 'got value', U(iofzero+1,iofzero+1)
  ! ---------- corners ----------------%
  call dgemm('n','n',mm2,mm2,mm2,1.d0,PNT,mm2,FT ,mm2,0.d0,wk1,mm2)
  call dgemm('n','t',mm2,mm2,mm2,1.d0,wk1,mm2,PNT,mm2,0.d0,wk2,mm2)

  valuesw=0.d0
  valuese=0.d0
  valuenw=0.d0
  valuene=0.d0
  do i=1,mm2
     ii=i-1
     sswi=dfloat((-1)**ii)
     ssei=1.d0
     snwi=dfloat((-1)**ii)
     snei=1.d0
     vpol=dfloat(ii+1)/dfloat((ii+2)**2)
     vintsw=0.d0
     vintse=0.d0
     vintnw=0.d0
     vintne=0.d0
     do j=1,mm2
        jj=j-1
        sswj=dfloat((-1)**jj)
        ssej=dfloat((-1)**jj)
        snwj=1.d0
        snej=1.d0
        coei=dfloat(jj+1)/dfloat((jj+2)**2)
        vintsw=vintsw+sswj*coei*wk2(i,j)
        vintse=vintse+ssej*coei*wk2(i,j)
        vintnw=vintnw+snwj*coei*wk2(i,j)
        vintne=vintne+snej*coei*wk2(i,j)
     enddo
     valuesw=valuesw+sswi*vpol*vintsw
     valuese=valuese+ssei*vpol*vintse
     valuenw=valuenw+snwi*vpol*vintnw
     valuene=valuene+snei*vpol*vintne
  enddo
  U(1,1)=16.d0*valuesw
  U(m,1)=16.d0*valuese
  U(1,m)=16.d0*valuenw
  U(m,m)=16.d0*valuene
  ! ---------- edges ----------------%
  call dgemm('n','n',mm2,mm2,mm2,1.d0,PNT,mm2,FT, &
       &             mm2,0.d0,wk1,mm2)
  call dgemm('n','t',mm2,mm2,mm2,1.d0,FT,mm2,PNT, &
       &             mm2,0.d0,wk2,mm2)
  ! left & right edge
  do j=2,m-1
     vl=0.d0
     vr=0.d0
     do i=1,m-2
        c=dfloat(i)/dfloat((i+1)**2)
        s=dfloat((-1)**(i-1))
        vl=vl+s*c*wk1(i,j-1)
        vr=vr+c*wk1(i,j-1)
     enddo
     U(1,j)=4.d0*vl
     U(m,j)=4.d0*vr
  enddo
  ! bottom & top edge
  do i=2,m-1
     vl=0.d0
     vr=0.d0
     do j=1,m-2
        c=dfloat(j)/dfloat((j+1)**2)
        s=dfloat((-1)**(j-1))
        vl=vl+s*c*wk2(i-1,j)
        vr=vr+c*wk2(i-1,j)
     enddo
     U(i,1)=4.d0*vl
     U(i,m)=4.d0*vr
  enddo
  
  return
end subroutine solvepois
