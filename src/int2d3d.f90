#ifndef RECTANGULAR
subroutine int2d3d(sol,rhs,work,bound,xlambda,imode)
  use ctes,only:nk,my,mz,mye,myo
  use numbers
  implicit none
  !
  ! Interface to the 2d solver
  ! note: if bound='d' -> rhs contains boundary values in the 
  !          proper position!!
  ! 
  complex(nk) rhs(my,mz),sol(my,mz)
  real(nk) work(*)
  character bound
  real(nk) xlambda
  integer imode,iw1,iw2,iw3,iw4,iw5,iw6,iw7  

#ifdef HELM_SHEN
  ! wrapped version from LEGCHEB
  iw1=1
  iw2=(my+1)*(mz+1)+iw1
  iw3=(my-2 +1)*(mz-2 +1)+iw2
  call com2re(my,mz,rhs,my,work(iw1),my+1,'r')
  call solveshen(work(iw1),work(iw2),work(iw3), &
       &         bound,xlambda)
  call re2com(my,mz,sol,my,work(iw1),my+1,'r')
  
  if(imode.eq.0) return
  
  call com2re(my,mz,rhs,my,work(iw1),my+1,'i')
  call solveshen(work(iw1),work(iw2),work(iw3), &
       &         bound,xlambda)
  call re2com(my,mz,sol,my,work(iw1),my+1,'i')
  
#else
  
  iw1=1
  iw2=my*mz+iw1
  iw3=(my-2+1)*(mz-2+1)+iw2
  call com2re(my,mz,rhs,my,work(iw1),my,'r')
  call solvediag(work(iw1),work(iw2),work(iw3),bound,xlambda)
  call re2com(my,mz,sol,my,work(iw1),my,'r')
  
  if(imode.eq.0) return
  
  call com2re(my,mz,rhs,my,work(iw1),my,'i')
  call solvediag(work(iw1),work(iw2),work(iw3),bound,xlambda)
  call re2com(my,mz,sol,my,work(iw1),my,'i')
     
#endif /* HELM_SHEN */
        
  return

end subroutine int2d3d
!-----|---------------------------------------------------------------|
subroutine com2re(n,m,vc,lvc,vr,lvr,iwhat)
  use ctes,only:nk
  implicit none
  integer n,m,lvc,lvr,i,j
  complex(nk) vc(lvc,*)
  real(nk) vr(lvr,*)
  character*1 iwhat
  
  if (iwhat.eq.'r') then
     do j=1,m
        do i=1,n
           vr(i,j)=dreal(vc(i,j))
        enddo
     enddo
  else if (iwhat.eq.'i') then
     do j=1,m
        do i=1,n
           vr(i,j)=dimag(vc(i,j))
        enddo
     enddo
  else
     write(*,*) 'mistake in com2re'
     stop
  endif
  return
end subroutine com2re
!-----|---------------------------------------------------------------|
subroutine re2com(n,m,vc,lvc,vr,lvr,iwhat)
  use ctes,only:nk
  implicit none
  integer n,m,lvc,lvr
  complex(nk) vc(lvc,*)
  real(nk) vr(lvr,*)
  character*1 iwhat
  integer i,j

  if (iwhat.eq.'r') then
     do j=1,m
        do i=1,n
           vc(i,j)=dcmplx(vr(i,j),0.d0)
        enddo
     enddo
  else if (iwhat.eq.'i') then
     do j=1,m
        do i=1,n
           vc(i,j)=vc(i,j)+dcmplx(0.d0,vr(i,j))
        enddo
     enddo
  else
     write(*,*) 'mistake in re2com'
     stop
  endif
      
  return
end subroutine re2com
#else /* RECTANGULAR */
subroutine int2d3d_dummy()
end subroutine int2d3d_dummy
#endif /* not RECTANGULAR */
