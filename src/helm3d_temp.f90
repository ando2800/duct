#ifdef TEMPERATURE
subroutine helm3d_temp(rhs,sol,work,ipro,cvar, &
     &                       a,b,alambda,bcco)
  ! --------------------------------------------------------------!
  !   solve laplacian(u)-alambda*u = f                       
  !   if bound='d' ---> Dirichelet problem                   
  !   if bound='n' ---> Neumann problem                      
  !   cvar='u','v','w' or 't'-> indicates which equations is 
  !   being considered (only for bound='d')                  
  !   if bound='n' -> ignored                                
  !   bcco  real coefficient used for dirichlet b.c. (bound='d')
  ! --------------------------------------------------------------!
  !    computational domain: [-a,a]x[-b,b]                
  ! --------------------------------------------------------------!
  use ctes
  use numbers
  use chebystuff_temp
  implicit none 

  complex(nk) rhs(my,mz,ib:ie),sol(my,mz,ib:ie)
  ! /* work array for temperature */
  real(nk) work(*)
  integer ipro
  character*1 cvar
  !character boundy1,boundy2,boundz1,boundz2
  !real(nk) bcy(mz,2,ib:ie),bcz(my,2,ib:ie) 
  real(nk) a,b,xlambda,alambda,bcco
  integer id2b,id2,ieii,idd,ifree,imode
  
  if(ipro.le.0) then
     sy=a
     sz=b
     call pseudomaty_temp(boundy1,boundy2,my)
     call pseudomatz_temp(boundz1,boundz2,mz)
     if(ipro.eq.-1)return   !allows for pre-comp only
  endif
  
  do imode=ib,ie
     
     xlambda=alambda+xalp(imode)**2
     call int2d3d_temp(sol(1,1,imode),rhs(1,1,imode), &
          &            work,boundy1,boundy2,boundz1,boundz2, &
          &            xlambda,imode,bcy(1,1,imode),bcz(1,1,imode))
  enddo
  return
end subroutine helm3d_temp
#else
subroutine helm3d_temp_dummy()
end subroutine helm3d_temp_dummy
#endif /* TEMPERATURE */
