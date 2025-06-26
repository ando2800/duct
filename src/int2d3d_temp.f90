subroutine int2d3d_temp(sol,rhs,work,boundy1,boundy2,boundz1, &
     &                  boundz2,xlambda,imode,bcy,bcz)
  use ctes,only:nk,my,mz
! --------------------------------------------------------------
!   Interface to the 2d solver
! --------------------------------------------------------------
  implicit none
  complex(nk) rhs(my,mz),sol(my,mz)
  real(nk) work(*)
  character*1 boundy1,boundy2,boundz1,boundz2
  real(nk) bcy(mz,2),bcz(my,2) 
  integer iw1,iw2,iw3,ifree,imode
  real(nk) xlambda
! --------------------------------------------------------------
  iw1=1
  iw2=my*mz+iw1
  iw3=(my-2)*(mz-2)+iw2
  ifree=(my-2)*(mz-2)+iw3
  call com2re(my,mz,rhs,my,work(iw1),my,'r')

  call solvediag_temp(work(iw1),work(iw2),work(iw3), &
       &              boundy1,boundy2,boundz1,boundz2,xlambda, &
       &              bcy,bcz)
  
  call re2com(my,mz,sol,my,work(iw1),my,'r')
      
  if(imode.eq.0)return
      
  call com2re(my,mz,rhs,my,work(iw1),my,'i')
      
  call solvediag_temp(work(iw1),work(iw2),work(iw3), &
       &              boundy1,boundy2,boundz1,boundz2,xlambda, &
       &              bcy,bcz)
      
  call re2com(my,mz,sol,my,work(iw1),my,'i')
  
  return
end subroutine int2d3d_temp


