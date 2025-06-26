#ifdef HELM_SHEN
subroutine helm3d_shen(rhs,sol,work,ipro,bound,cvar,a,b,alambda)

  use ctes, only: nk,ib,ie,my,mz,mx1,xalp,iax
  use timing, only: timeshen
  use numbers

  implicit none
  !  /********************************************************/
  !  /*                                                      */
  !  /* solve [- laplacian(u) + alambda*u = f]               */
  !  /*    note: set f=-f to be comparable with collocation  */
  !  /*                                         method       */
  !  /* if bound='d' ---> Dirichelet problem                 */
  !  /* if bound='n' ---> Neumann problem                    */
  !  /* cvar='u','v','w' or 't'-> indicates which equations  */
  !  /* is being considered                                  */
  !  /*                                                      */
  !  /* computational domain: [-a,a]x[-b,b]                  */
  !  /*                                                      */
  !  /********************************************************/
  complex(nk) rhs(my,mz,ib:ie),sol(my,mz,ib:ie)
  ! /* work array */
  real(nk) work(*)
  integer ipro,id2b,id2,imode
  character*1 cvar
  character bound
  real(nk) a,b,alambda,xlambda
  integer i,j
  !
#ifdef DEBUG
  real*8 t0,t1,t8,t9
#endif
  !    
  if(ipro.le.0) then
     !  /* set common arrays */
     id2b=1
     id2=id2b+(my-2)**2     
#ifdef DEBUG
     write(*,*) 'preshen'
#endif
     call preshen(work(id2b),work(id2),bound)
     
     if(ipro.eq.-1) return   !allows for pre-comp only
  endif
     
  ! /* mx1 is set to 0 in ctes3D */
#ifdef DEBUG
  !write(*,*) 'cheb trans. '
  call cpu_time(t0)
#endif      
  !  rhs is OK !!!
  call transyz(rhs,'y',-1)  ! phys -> cheb
#ifdef DEBUG
  call cpu_time(t1)
  timeshen(1)= timeshen(1) + t1-t0
  call cpu_time(t0)
#endif
  call transyz(rhs,'z',-1)  ! phys -> cheb
#ifdef DEBUG
  call cpu_time(t1)
  timeshen(9)= timeshen(9) + t1-t0
#endif
  do imode=ib,ie
     xlambda=alambda+xalp(imode)**2
     call int2d3d(sol(1,1,imode),rhs(1,1,imode), &
          &       work,bound,xlambda,imode)
  enddo
#ifdef DEBUG
  ! write(*,*) 'backward cheb trans.'
  call cpu_time(t8)
#endif      
  call transyz(sol,'y',+1)  ! cheb -> phys
#ifdef DEBUG
  call cpu_time(t9)
  timeshen(1)= timeshen(1) + t9-t8
  call cpu_time(t8)
#endif  
  call transyz(sol,'z',+1)  ! cheb -> phys
#ifdef DEBUG
  call cpu_time(t9)
  timeshen(9)= timeshen(9) + t9-t8
#endif  
  return
end subroutine helm3d_shen
#else
subroutine helm3d_shen_dummy
end subroutine helm3d_shen_dummy
#endif /* HELM_SHEN */
