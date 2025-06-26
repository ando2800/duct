#ifndef RECTANGULAR      
subroutine helm3d(rhs,sol,work,ipro,bound,cvar,a,b,alambda)
  use ctes
  use numbers
  implicit none 
! -----------------------------------------------------------------!
!  solve laplacian(u)-alambda*u = f  ,X-CUT                    
!  if bound='d' ---> Dirichelet problem                        
!  if bound='n' ---> Neumann problem                           
!  cvar='u','v','w' or 't'-> indicates which equations is      
!            being considered (only for bound='d')             
!             if bound='n' -> ignored                          
! -----------------------------------------------------------------!
! NOTE THAT:                    
!        use ipro=0 for the very first call then use ipro=1  
!  M.U      ipro<=0                                  "  
!        ONLY WORKS for my=mz!!!!!
!        ONLY WORKS for homogeneous problems on the square
!                   and periodic in x direction
!        RHS as invalue; SOL as outvalue PHYS/PHYS/FOURIER both
!        a,b   are dummies (not accessed)
!              [for compatibility with RECTANGULAR version]
!
!        work(nwk) nwk1=(2*mz-2)+(mz-2)^2+mz^2
!                  nwk2=my*mz+2*(my-2)*(mz-2)
!                  nwk=max(nwk1,nwk3)
!                where: my=mz
! -----------------------------------------------------------------!
  complex(nk) rhs(my,mz,ib:ie),sol(my,mz,ib:ie)
  real(nk) work(*)
  integer ipro,imode,i,j,k,id2b,id2,ieii,idd,ifree
  integer ip1,ip2,ip3,ip4
  real(nk) a,b,xlambda,alambda
  character*1 bound,cvar
  !-----|---------------------------------------------------------------|
  if(ipro.le.0) then
     !     check my=mz, eventually stop
     if(my.ne.mz) then
        write(*,*) 'IN HELM3D: MY.NE.MZ (NOT IMPLEMENTED). ABORTING'
        stop 
     endif
     !     
     call pseudomat('d',my)
     call pseudomat('n',my)
     if(ipro.eq.-1) return   !allows for pre-comp only
  endif
  
  if(bound.eq.'d') then      !set dirichlet values in rhs
     !     /* homogeneous dirichlet in all 4 wall-planes */
     do i=ib,ie
        do j=1,my
           rhs(j,1 ,i)=zero
           rhs(j,mz,i)=zero
        enddo
        do k=1,mz
           rhs(1 ,k,i)=zero
           rhs(my,k,i)=zero
        enddo
     enddo
  endif
  
  ! start loop on the mode number
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(imode,xlambda) 
  !do imode=ib,ie
  do imode=ibf1,ief1
     !write(*,*) myid,': ompid=',ompid,': work at ',nwke*ompid+1
     xlambda=alambda+xalp(imode)**2
     ! the work array should be a private variable
     call int2d3d(sol(1,1,imode),rhs(1,1,imode), & 
          &       work(nwke*ompid+1),bound,xlambda,imode)
     
  enddo
  !$OMP END PARALLEL
  
  return
end subroutine helm3d
#else /* RECTANGULAR */
subroutine helm3d_dummy()
end subroutine helm3d_dummy
#endif /* not RECTANGULAR */
