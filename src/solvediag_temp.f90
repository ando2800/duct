subroutine solvediag_temp(U,F,FT,boundy1,boundy2,boundz1, &
     &                           boundz2,xlambda,bcy,bcz)
  use ctes,only:my,mz
  use fulldiagtemp
  use chebystuff_temp
  implicit none 
  !     
  ! solve the problem
  ! bound either Neumann: n, or Dirichlet:d.
  !     
  real(nk) bcy(mz,2),bcz(my,2) 
  
  real(nk) U(my,mz),F(my-2,mz-2)
  real(nk) FT(my-2,mz-2)
  real(nk) ff(4),uu(4)
  character*1 boundy1, boundy2, boundz1, boundz2
  real(nk) bcontv(my-2,mz-2), bconth(my-2,mz-2)
  integer i,j,ny,mm2y,nz,mm2z,ii,jj
  real(nk) den,xlambda

  do j=2,mz-1
     do i=2,my-1
        F(i-1,j-1)=U(i,j)
     enddo
  enddo
      
  ny=my-1
  mm2y=my-2
  
  nz=mz-1
  mm2z=mz-2
  !---------dirichiret B.C -------------
  !     
  ! correct rhs for boundary condition
  !     
  call mod_term_FTy(boundy1,boundy2,bcontv,bcy)
  call mod_term_FTz(boundz1,boundz2,bconth,bcz)
  do j=1,mz-2
     do i=1,my-2
        F(i,j)=F(i,j)+bcontv(i,j)+bconth(i,j)
     enddo
  enddo
  !     
  ! Projection of rhs
  ! F=inv(Py)*F*inv(Pz');
  !     
  call dgemm('n','n',mm2y,mm2z,mm2y,1.d0,PTYM1,mm2y,F,mm2y,0.d0,FT,mm2y)
  call dgemm('n','t',mm2y,mm2z,mm2z,1.d0,FT,mm2y,PTZM1,mm2z,0.d0,F,mm2y)
  do j=1,mm2z
     do i=1,mm2y
        den=(ETY(i)+ETZ(j)-xlambda)
        if(abs(den).eq.1.d-12)then
           write(*,*) 'solvediag: xlamda is 0!! stop'
           stop
        endif
        FT(i,j)=F(i,j)/den
     enddo
  enddo
  !  reconstruct physical sol
  !  FT=Py*FT*Pz';
  !     
  
  call dgemm('n','n',mm2y,mm2z,mm2y,1.d0,PTY,mm2y,FT,mm2y,0.d0,F,mm2y)
  call dgemm('n','t',mm2y,mm2z,mm2z,1.d0,F,mm2y,PTZ,mm2z,0.d0,FT,mm2y)

  !  fill in the solution
  
  do j=1,mm2z
     do i=1,mm2y
        U(i+1,j+1)=FT(i,j)
     enddo
  enddo
            
  !     end procedure to fix the solution in (iofzero,jofzero)
  !$$$  write(6,*) 'wanted',wanted
  !$$$  write(6,*) 'I got',FT(iofzero+1,jofzero+1)

         
  !   Fix Neumann values on the border
  call edge_value_y(U,boundy1,boundy2,bcy)
  call edge_value_z(U,boundz1,boundz2,bcz)
  
  return
end subroutine solvediag_temp
!----+---------------------------------------------------------------+
subroutine mod_term_FTy(bc1,bc2,bcontv,bcy)
  use ctes,only:nk,my,mz
  use numbers
  use chebystuff_temp
  implicit none

  real(nk) bcy(mz,2) 

  character*1 bc1,bc2
  real(nk) bcontv(my-2,mz-2)
  integer ii,i,jj,j
  real(nk) t1,t2,Det

  !c%
  !c%build pseudo matrix with B.C.
  !c%usage::
  !c%      :input D   ; 1st derivative matrix[m*m] (chebyshev)
  !c%           % D2  ; 2nd derivative matrix[m*m] (chebyshev)
  !c%             bc1 ; B.C flag at x= -1.('d':Dirichlet,'n':Neumann)
  !c%             v1??  ; value of B.C at x= -1.
  !c%             bc2 ; B.C flag at x= +1.('d':Dirichlet,'n':Neumann)
  !c%             v2??  ; value of B.C at x= +1.
  !c%      :output 
  !c%             Tx(1:m-2); 1D modification term of r.h.s of Helmholta eq.
  !c%
  if((bc1.eq.'d').and.(bc2.eq.'d'))then
     ii=0
     do i=2,my-1
        ii=ii+1
        jj=0
        do j=2,mz-1
           jj=jj+1
           bcontv(ii,jj)=-D2by(i,1)*bcy(j,1)-D2by(i,2)*bcy(j,2)
        enddo
     enddo
     
  else if((bc1.eq.'d').and.(bc2.eq.'n'))then
     ii=0
     do i=2,my-1
        ii=ii+1
        jj=0
        do j=2,mz-1
           jj=jj+1
           bcontv(ii,jj)=-D2by(i,1)*bcy(j,1) & 
                &       - D2by(i,2)*(bcy(j,2)-DY(my,1)*bcy(j,1))/DY(my,my)
        enddo
     enddo
     
  else if((bc1.eq.'n').and.(bc2.eq.'d'))then
     ii=0
     do i=2,my-1
        ii=ii+1
        jj=0
        do j=2,mz-1
           jj=jj+1  
           bcontv(ii,jj)=-D2by(i,1)*(bcy(j,1) &
                &            - DY(1,my)*bcy(j,2))/DY(1,1) &
                &            - D2by(i,2)*bcy(j,2)
        enddo
     enddo
     
  else if((bc1.eq.'n').and.(bc2.eq.'n'))then
     ii=0
     Det=DY(1,1)*DY(my,my)-DY(1,my)*DY(my,1);
     do i=2,my-1
        ii=ii+1
        jj=0
        do j=2,mz-1
           jj=jj+1
           t1=D2by(i,1)*(DY(1,my)*bcy(j,2)-DY(my,my)*bcy(j,1))/Det
           t2=D2by(i,2)*(DY(my,1)*bcy(j,1)-DY(1,1)*bcy(j,2))/Det
           bcontv(ii,jj)=t1+t2
        enddo
     enddo
  endif
  return 
end subroutine mod_term_FTy
      
!----+---------------------------------------------------------------+
subroutine mod_term_FTz(bc1,bc2,bconth,bcz)
  use ctes,only:nk,my,mz
  use chebystuff_temp
  implicit none

  real(nk) bcz(my,2) 

  character*1 bc1,bc2
  real(nk) bconth(my-2,mz-2)
  integer i,j,ii,jj 
  real(nk) Det,t1,t2
  !c%
  !c%build pseudo matrix with B.C.
  !c%usage::
  !c%      :input D   ; 1st derivative matrix[m*m] (chebyshev)
  !c%           % D2  ; 2nd derivative matrix[m*m] (chebyshev)
  !c%             bc1 ; B.C flag at x= -1.('d':Dirichlet,'n':Neumann)
  !c%        bcz(i,1) ; value of B.C at x= -1.
  !c%             bc2 ; B.C flag at x= +1.('d':Dirichlet,'n':Neumann)
  !c%        bcz(i,2) ; value of B.C at x= +1.
  !c%      :output 
  !c%             Tx(1:m-2); 1D modification term of r.h.s of Helmholta eq.
  !c%
  if((bc1.eq.'d').and.(bc2.eq.'d'))then
     jj=0
     do j=2,mz-1
        jj=jj+1
        ii=0
        do i=2,my-1
           ii=ii+1
           bconth(ii,jj)=-D2bz(j,1)*bcz(i,1)-D2bz(j,2)*bcz(i,2)
        enddo
     enddo
     
  else if((bc1.eq.'d').and.(bc2.eq.'n'))then
     jj=0
     do j=2,mz-1
        jj=jj+1
        ii=0
        do i=2,my-1
           ii=ii+1
           bconth(ii,jj)=-D2bz(j,1)*bcz(i,1)-D2bz(j,2)*(bcz(i,2) &
                &                 - DZ(mz,1)*bcz(i,1))/DZ(mz,mz)
        enddo
     enddo
     
  else if((bc1.eq.'n').and.(bc2.eq.'d'))then
     jj=0
     do j=2,mz-1
        jj=jj+1
        ii=0
        do i=2,my-1
           ii=ii+1  
           bconth(ii,jj)=-D2bz(j,1)*(bcz(i,1) &
                &            - DZ(1,mz)*bcz(i,2))/DZ(1,1) &
                &            - D2bz(j,2)*bcz(i,2)
        enddo
     enddo
     
  else if((bc1.eq.'n').and.(bc2.eq.'n'))then
     jj=0
     Det=DZ(1,1)*DZ(mz,mz)-DZ(1,mz)*DZ(mz,1);
     do j=2,mz-1
        jj=jj+1
        ii=0
        do i=2,my-1
           ii=ii+1
           t1=D2bz(j,1)*(DZ(1,mz)*bcz(i,2)-DZ(mz,mz)*bcz(i,1))/Det
           t2=D2bz(j,2)*(DZ(mz,1)*bcz(i,1)-DZ(1,1)*bcz(i,2))/Det
           bconth(ii,jj)=t1+t2
        enddo
     enddo
  endif
  return
end subroutine mod_term_FTz
!----+--------------------------------------------------------------+
subroutine edge_value_y(U,boundy1,boundy2,bcy)
  use ctes,only:nk,my,mz
  use chebystuff_temp
  implicit none
  
  real(nk) bcy(mz,2) 
  real(nk) U(my,mz)
  character*1 boundy1,boundy2
  integer i,j
  real(nk) sum,a,b,Det
  
  if((boundy1.eq.'d').and.(boundy2.eq.'d'))then
     do j=2,mz-1
        U(1,j)=bcy(j,1)
        U(my,j)=bcy(j,2)
        ! if(abs(bcy(j,1)).gt.0.5d0) then
        !    write(*,*) 'U,bcy y=-1: ',U(1,j),bcy(j,1)
        !    write(*,*) 'U,bcy y=+1: ',U(my,j),bcy(j,2)
        ! endif
     enddo
     
  elseif((boundy1.eq.'d').and.(boundy2.eq.'n'))then
     do j=2,mz-1
        sum=0.d0
        do i=2,my-1
           sum=sum+(DY(my,i)*U(i,j))
        enddo
        
        U(1,j)=bcy(j,1)
        U(my,j)=(bcy(j,2)-DY(my,1)*U(1,j) - sum)/DY(my,my);
     enddo
  elseif((boundy1.eq.'n').and.(boundy2.eq.'d'))then
     do j=2,mz-1
        sum=0.d0
        do i=2,my-1
           sum=sum+(DY(1,i)*U(i,j))
        enddo
        
        U(my,j)=bcy(j,2)
        U(1,j)=(bcy(j,1)-DY(1,my)*U(my,j) - sum)/DY(1,1);
     enddo
  elseif((boundy1.eq.'n').and.(boundy2.eq.'n'))then
     Det=DY(1,1)*DY(my,my)-DY(1,my)*DY(my,1);
     do j=2,mz-1
        a=0.d0
        b=0.d0
        do i=2,my-1
           a=a+(DY(1,i)*U(i,j));
           b=b+(DY(my,i)*U(i,j));
        enddo
        U(1,j)=(DY(my,my)*(bcy(j,1)-a)-DY(1,my)*(bcy(j,2)-b))/Det;
        U(my,j)=(-DY(my,1)*(bcy(j,1)-a)+DY(1,1)*(bcy(j,2)-b))/Det;
     enddo
     
  endif
  return
end subroutine edge_value_y
!----+---------------------------------------------------------
subroutine edge_value_z(U,boundz1,boundz2,bcz)
  use ctes,only:nk,my,mz
  use chebystuff_temp
  implicit none
  
  real(nk) bcz(my,2) 
  real(nk) U(my,mz)
  character*1 boundz1,boundz2
  integer i,j 
  real(nk) sum,Det,a,b
    
  if((boundz1.eq.'d').and.(boundz2.eq.'d'))then
     do i=2,my-1
        U(i,1)=bcz(i,1)
        U(i,mz)=bcz(i,2)
     enddo

  elseif((boundz1.eq.'d').and.(boundz2.eq.'n'))then
     do i=2,my-1
        sum=0.d0
        do j=2,mz-1
           sum=sum+(DZ(mz,j)*U(i,j))
        enddo
        
        U(i,1)=bcz(i,1)
        U(i,mz)=(bcz(i,2)-DZ(mz,1)*U(i,1) - sum)/DZ(mz,mz);
     enddo
  elseif((boundz1.eq.'n').and.(boundz2.eq.'d'))then
     do i=2,my-1
        sum=0.d0
        do j=2,mz-1
           sum=sum+(DZ(1,j)*U(i,j))
        enddo
        
        U(i,mz)=bcz(i,2)
        U(i,1)=(bcz(i,1)-DZ(1,mz)*U(i,mz) - sum)/DZ(1,1);
     enddo
  elseif((boundz1.eq.'n').and.(boundz2.eq.'n'))then
     Det=DZ(1,1)*DZ(mz,mz)-DZ(1,mz)*DZ(mz,1);
     do i=2,my-1
        a=0.d0
        b=0.d0
        do j=2,mz-1
           a=a+(DZ(1,j)*U(i,j));
           b=b+(DZ(mz,j)*U(i,j));
        enddo
        U(i,1)=(DZ(mz,mz)*(bcz(i,1)-a)-DZ(1,mz)*(bcz(i,2)-b))/Det;
        U(i,mz)=(-DZ(mz,1)*(bcz(i,1)-a)+DZ(1,1)*(bcz(i,2)-b))/Det;
     enddo
          
  endif
  return
end subroutine edge_value_z

