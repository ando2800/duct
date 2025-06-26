#define DEBUG_adv0 
#undef DEBUG_adv0
#define DEBUG_adv1 
#undef DEBUG_adv1
#define DEBUG_adv6 
#undef DEBUG_adv6
#define DEBUG_advtemp1 
#undef DEBUG_advtemp1
      subroutine advection3d(u,v,w,adv1,adv2,adv3,
     $     tmp1,tmp2,tmp3,ufis,vfis,wfis,tmpzxcut,tmpx,
     $     krk,iter)      
      use ctes
      use numbers
      use running
      implicit none 
      include "mpif.h"
c     
c     /***************************************************************/
c     /* computes the advection terms of navier-stokes in diver-     */
c     /* gence formulation, for a velocity field which is re-        */
c     /* presented by fourier expansion in x and chebyshev in y,z,   */
c     /* pseudo-spectrally.                                          */
c     /*                                                             */
c     /* input                                                       */
c     /* adv1-3(nsfou) =u,v,w velocity components in fou-phys-phys   */
c     /* (Z-cut)              space; the arrays are sufficiently     */
c     /*                      dimensioned to allow for in-place      */
c     /*                      fourier transform (size: nsfis)        */
c     /* tmp1(nsfis)          BIG work array                         */
c     /* tmp2(nsfou)          SMALL work array                       */
c     /* tmp3(nsfou)          SMALL work array                       */
c     /* ufis,vfis,wfis(nsfis) BIG work arrays                       */
c     /* tmpzxcut(my*mzp*mxp) SMALLER complex*16 work array          */
c     /* tmpx(mgalx+2)        1d-in-x work array                     */
c     /* krk                  runge-kutta step-count                 */
c     /* cfl                  cfl condition                          */
c     /*                                                             */
c     /* output                                                      */
c     /* adv1,adv2,adv3(nsfou) advection terms in f-p-p space        */
c     /* dt                    admisible time step                   */
c     /*                                                             */
c     /* note: 'real*8 u(nsfour)' and 'complex*16 u(nsfou)' have     */
c     /*       equivalent size                                       */
c     /***************************************************************/
      real(nk) u(nsfour),v(nsfour),w(nsfour),
     &     adv1(nsfour),adv2(nsfour),adv3(nsfour)
      real(nk) tmp1(nsfis),tmp2(nsfour),tmp3(nsfour)
      real(nk) ufis(nsfis),vfis(nsfis),wfis(nsfis)
c     /* tmp work array for shuffling (note: complex*16)*/
      complex*16 tmpzxcut(my*mzp*mxp) !z-cut & x-cut pencil work array
c     /* local work */
      real(nk) tmpx(mgalx+2),dummy
      integer krk,iter,ierr
c
c-----|---------------------------------------------------------------|
c     /* preserve u, v, w (x-cut) */
      call dcopy(nsfour,u,1,adv1,1)  
      call dcopy(nsfour,v,1,adv2,1)
      call dcopy(nsfour,w,1,adv3,1)
c            
c     /* shuffle arrays from x-cut to z-cut */
      call chx2z(adv1,tmp2,tmpzxcut)
      call chx2z(adv2,tmp2,tmpzxcut)
      call chx2z(adv3,tmp2,tmpzxcut)
c
c     /* pass all three components to phys-phys-phys space (NOT-in-place) */
      call fourx(adv1,ufis,tmpx,+1)!fou->fis
      call fourx(adv2,vfis,tmpx,+1)!fou->fis
      call fourx(adv3,wfis,tmpx,+1)!fou->fis
c
c     /* (00) ========= check the slip velocity========================= */
      if(krk.eq.1) call slip_velocity_3d(ufis,vfis,wfis)
c
c     /* (0) ========== determine/check the time step from cfl condition */
      if(krk.eq.1) call determine_step(ufis,vfis,wfis) ! --> dt
#ifdef SET_NOSLIP
c     /* (01) ========= forced to set no slip for ufis,vfis,wfis ======*/
c     /*   by sekimoto 12/Sep/2008 */
      call set_noslip_3dfis(ufis,vfis,wfis)
      call set_noslip_3dfou(u,v,w)
#endif
c
c     /* (1) ========== compute product 'u*u' -> 'tmp1' */
c      if (myid.eq.master) write(*,*) 'advection3d: mult'
      call mult(nsfis,ufis,ufis,tmp1,'n','+')
c
c     /* pass 'u*u' back to f-p-p space -> adv1 */
c      if (myid.eq.master) write(*,*) 'advection3d: fourx'
      call fourx(adv1,tmp1,tmpx,-1)!fis->fou
c
c     /* Re-shuffle arrays from z-cut to x-cut */
      tmp2=0.d0; tmpzxcut=0.d0;
      call chz2x(adv1,tmp2,tmpzxcut)
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c
c     /* call statistics module */
      if((krk.eq.1).and.(mod(iter-1,nhist).eq.0))then
c         if (myid.eq.master) write(*,*) 'advection3d: statistics_fou2'
         call statistics_fou2(adv1,1) !u*u
      endif
c
c     /* derivative in 'x': d(u*u)/dx -> adv1 */
      call derivx(adv1,adv1,1,f1o1,'n')
c
c     /* (2) =========== compute product 'v*v' -> 'tmp1' */
      call mult(nsfis,vfis,vfis,tmp1,'n','+')
c
c     /* pass 'v*v' back to f-p-p space -> adv2 */
      call fourx(adv2,tmp1,tmpx,-1)!fis->fou
c
c     /* Re-shuffle arrays from z-cut to x-cut */
      tmp2=0.d0; tmpzxcut=0.d0;
      call chz2x(adv2,tmp2,tmpzxcut)
c
c     /* call statistics module */
      if((krk.eq.1).and.(mod(iter-1,nhist).eq.0))then
c         if (myid.eq.master) write(*,*) 'advection3d: statistics_fou2 v'
         call statistics_fou2(adv2,2) !v*v
      endif
c
c     /* pass 'v*v' to f-c-p space -> adv2 */
      call transyz(adv2,'y',-1) !fis->cheb
c
c     /* derivative in 'y': d(v*v)/dy -> tmp2 */
      call derivyz(adv2,tmp2,dummy,'y',1,f1o1,'n') !dummy??
c
c     /* pass d(v*v)/dy to f-p-p space -> tmp2 */
      call transyz(tmp2,'y',+1) !cheb->fis
c
c     /* copy d(v*v)/dy into 'adv2' */
      call dcopy(nsfour,tmp2,1,adv2,1)
c
c     /* (3) ============ compute product 'w*w' -> 'tmp1' */
      call mult(nsfis,wfis,wfis,tmp1,'n','+')
c
c     /* pass 'w*w' back to f-p-p space -> adv3 */
      call fourx(adv3,tmp1,tmpx,-1)!fis->fou
c
c     /* Re-shuffle arrays from z-cut to x-cut */
      tmp2=0.d0; tmpzxcut=0.d0;
      call chz2x(adv3,tmp2,tmpzxcut)
c
c     /* call statistics module */
      if((krk.eq.1).and.(mod(iter-1,nhist).eq.0))then
c         if (myid.eq.master) write(*,*) 'advection3d: statistics_fou2 w'
         call statistics_fou2(adv3,3) !w*w
      endif
c
c     /* pass 'w*w' to f-p-c space -> adv3 */
      call transyz(adv3,'z',-1) !fis->cheb
c
c     /* derivative in 'z': d(w*w)/dz -> tmp2 */
      call derivyz(adv3,tmp2,dummy,'z',1,f1o1,'n')
c
c     /* pass d(w*w)/dz to f-p-p space -> tmp2 */
      call transyz(tmp2,'z',+1) !cheb->fis
c
c     /* copy d(w*w)/dz into 'adv3' */
      call dcopy(nsfour,tmp2,1,adv3,1)
c
c     /* (4) ============ compute product 'u*v' -> 'tmp1' */
      call mult(nsfis,ufis,vfis,tmp1,'n','+')
c
c     /* pass 'u*v' back to f-p-p space -> tmp2 */
      call fourx(tmp2,tmp1,tmpx,-1)!fis->fou
c
c     /* Re-shuffle arrays from z-cut to x-cut */
      tmp3=0.d0; tmpzxcut=0.d0;
      call chz2x(tmp2,tmp3,tmpzxcut)
c
c     /* call statistics module */
      if((krk.eq.1).and.(mod(iter-1,nhist).eq.0))then
c         if (myid.eq.master) write(*,*) 'adv3d: statistics_fou2 uv'
         call statistics_fou2(tmp2,4) !u*v
      endif
c
c     /* derivative in 'x': d(u*v)/dx -> ++adv2 */
      call derivx(tmp2,adv2,1,f1o1,'+')
c
c     /* pass 'u*v' to f-c-p space -> tmp2 */
c      if (myid.eq.master) write(*,*) 'adv3d: transyz '
      call transyz(tmp2,'y',-1) !fis->cheb
c
c     /* derivative in 'y': d(u*v)/dy -> tmp1 */
c      if (myid.eq.master) write(*,*) 'adv3d: derivyz '
      call derivyz(tmp2,tmp1,dummy,'y',1,f1o1,'n')
c
c     /* pass d(u*v)/dy to f-p-p space -> tmp1 */
c      if (myid.eq.master) write(*,*) 'adv3d: derivyz+1 '
      call transyz(tmp1,'y',+1) !cheb->fis
c
c     /* add d(u*v)/dy into 'adv1' */
c      if (myid.eq.master) write(*,*) 'adv3d: daxpy '
      call daxpy(nsfour,+f1o1,tmp1,1,adv1,1)
c
c     /* (5) ============ compute product 'u*w' -> 'tmp1' */
c      if (myid.eq.master) write(*,*) 'adv3d: mult '
      call mult(nsfis,ufis,wfis,tmp1,'n','+')
c
c     /* pass 'u*w' back to f-p-p space -> tmp2 */
c      if (myid.eq.master) write(*,*) 'adv3d: fourx '
      call fourx(tmp2,tmp1,tmpx,-1)!fis->fou
c
c     /* Re-shuffle arrays from z-cut to x-cut */
      tmp3=0.d0; tmpzxcut=0.d0;
c      write(*,*) myid, 'adv3d: chz2x '
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call chz2x(tmp2,tmp3,tmpzxcut)
c      write(*,*) myid, 'after advection3d: chz2x'
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c 
c     /* call statistics module */
      if((krk.eq.1).and.(mod(iter-1,nhist).eq.0))then
c         if (myid.eq.master) write(*,*) 'adv3d: statistics_fou2 uw'
         call statistics_fou2(tmp2,5) !u*w
      endif
c
c     /* derivative in 'x': d(u*w)/dx -> ++adv3 */
      call derivx(tmp2,adv3,1,f1o1,'+')
c
c     /* pass 'u*w' to f-p-c space -> tmp2 */
      call transyz(tmp2,'z',-1) !fis->cheb
c
c     /* derivative in 'z': d(u*w)/dz -> tmp1 */
      call derivyz(tmp2,tmp1,dummy,'z',1,f1o1,'n')
c
c     /* pass d(u*w)/dz to f-p-p space -> tmp1 */
      call transyz(tmp1,'z',+1) !cheb->fis
c
c     /* add d(u*w)/dz into 'adv1' */
      call daxpy(nsfour,+f1o1,tmp1,1,adv1,1)
c
c     /* (6) ============ compute product 'v*w' -> 'tmp1' */
      call mult(nsfis,vfis,wfis,tmp1,'n','+')
c
c     /* pass 'v*w' back to f-p-p space -> tmp2 */
      call fourx(tmp2,tmp1,tmpx,-1)!fis->fou
c
c     /* Re-shuffle arrays from z-cut to x-cut */
      tmp3=0.d0; tmpzxcut=0.d0;
      call chz2x(tmp2,tmp3,tmpzxcut)
c
c     /* call statistics module */
      if((krk.eq.1).and.(mod(iter-1,nhist).eq.0))then
         call statistics_fou2(tmp2,6) !v*w
      endif
c
c      call zcopy(my*mz*(ie-ib+1),tmp2,1,tmp3,1)  ! tmp2=tmp3 (x-cut)
      call zcopy(nsfou,tmp2,1,tmp3,1)  ! tmp2=tmp3 (x-cut)
c
c     /* pass 'v*w' to f-c-p space -> tmp2 */
      call transyz(tmp2,'y',-1) !fis->cheb
c
c     /* derivative in 'y': d(v*w)/dy -> tmp1 */
      call derivyz(tmp2,tmp1,dummy,'y',1,f1o1,'n')
c
c     /* pass d(v*w)/dy to f-p-p space -> tmp1 */
      call transyz(tmp1,'y',+1) !cheb->fis
c
c     /* add d(v*w)/dy to 'adv3' */
      call daxpy(nsfour,+f1o1,tmp1,1,adv3,1)
c
c     /* pass 'v*w' to f-p-c space -> tmp3 */
      call transyz(tmp3,'z',-1) !fis->cheb
c
c     /* derivative in 'z': d(v*w)/dz -> tmp1 */
      call derivyz(tmp3,tmp1,dummy,'z',1,f1o1,'n')
c
c     /* pass d(v*w)/dz to f-p-p space -> tmp1 */
      call transyz(tmp1,'z',+1) !cheb->fis
c
c     /* add d(v*w)/dz to 'adv2' */
      call daxpy(nsfour,+f1o1,tmp1,1,adv2,1)
c
      return
      end subroutine advection3d
c-----|---------------------------------------------------------------|
      subroutine laplacian3d(u,flap,tmp,tmpyz,scal,what)
      use ctes
      use numbers
      implicit none
c
c     /***************************************************************/
c     /* computes the laplacian of a velocity field which is re-     */
c     /* presented by fourier expansion in x and chebyshev in y,z,   */
c     /* pseudo-spectrally.                                          */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfou)             field in fou-phys-phys space  x-cut    */
c     /* tmp(nsfour)           fourier-size 3D work array            */
c     /* tmpyz(ntmp)          1D work array: ntmp=max(my,mz)         */
c     /* scal                 scalar factor                          */
c     /* what                 'n'  - overwrite array 'flap'          */
c     /*                         flap=scal*lapl(u)                   */
c     /*                      [else] - add to array 'flap'           */
c     /*                         flap=flap+scal*lapl(u)              */
c     /*                                                             */
c     /* output                                                      */
c     /* flap(nsfou)          laplacian(u) in f-p-p space            */
c     /* u(nsfou)             field in fou-phys-phys space           */
c     /*                      [preserved]                            */
c     /*                                                             */
c     /* note: 'real*8 u(nsfour)' and 'complex*16 u(nsfou)' have     */
c     /*       equivalent size                                       */
c     /***************************************************************/
      real(nk) u(nsfour),flap(nsfour)
      real(nk) tmp(nsfour),tmpyz(my+mz)
      real(nk) scal
      character*1 what
c
c-----|---------------------------------------------------------------|
c     /* 2nd derivative in 'x': scal*d^2(u)/dx^2 -> flap   */
c     /*   OR: add to previous value of 'flap' [what]      */
      call derivx(u,flap,2,scal,what)
c------------------------------------------
c     /* pass 'u' to f-c-p space -> u */
      call transyz(u,'y',-1) !fis->cheb
c
c     /* 2nd derivative in 'y': scal*d^2(u)/dy^2 -> tmp */
      call derivyz(u,tmp,tmpyz,'y',2,scal,'n')
c
c     /* pass 'dyyu' back to f-p-p space -> tmp */
      call transyz(tmp,'y',+1) !cheb->fis
c
c     /* add to 'flap' -> ++flap */
      call daxpy(nsfour,f1o1,tmp,1,flap,1)
c
c     /* pass 'u' back to f-p-p space -> u */
      call transyz(u,'y',+1) !cheb->fis
c--------------------------------------------
c     /* pass 'u' to f-p-c space -> u */
      call transyz(u,'z',-1) !fis->cheb
c
c     /* 2nd derivative in 'z': scal*d^2(u)/dz^2 -> tmp */
      call derivyz(u,tmp,tmpyz,'z',2,scal,'n')
c
c     /* pass 'dzzu' back to f-p-p space -> tmp */
      call transyz(tmp,'z',+1) !cheb->fis
c
c     /* add to 'flap' -> ++flap */
      call daxpy(nsfour,f1o1,tmp,1,flap,1)
c
c     /* pass 'u' back to f-p-p space -> u */
      call transyz(u,'z',+1) !cheb->fis
c
      return
      end subroutine laplacian3d
c-----|---------------------------------------------------------------|
      subroutine grad3d(u,dxu,dyu,dzu,tmp,tmpyz,scal,what)
      use ctes
      use numbers
      implicit none
c
c     /***************************************************************/
c     /* computes grad(u) of a field 'u' which is re-                */ 
c     /* presented by fourier expansion in x and chebyshev in y,z,   */
c     /* pseudo-spectrally.                                          */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfou)             field in fou-phys-phys space           */
c     /* tmp(nsfour)           fourier-size 3D work array            */
c     /* tmpyz(ntmp)          1D work array: ntmp=max(my,mz)         */
c     /* scal                 scalar factor                          */
c     /* what                 'n'  - overwrite array 'flap'          */
c     /*                         flap=scal*lapl(u)                   */
c     /*                      [else] - add to array 'flap'           */
c     /*                         flap=flap+scal*lapl(u)              */
c     /*                                                             */
c     /* output                                                      */
c     /* dxu,dyu,dzu(nsfou)   components of grad(u) in f-p-p space   */
c     /* u(nsfou)             field in fou-phys-phys space           */
c     /*                      [preserved]                            */
c     /*                                                             */
c     /* note: 'real*8 u(nsfour)' and 'complex*16 u(nsfou)' have     */
c     /*       equivalent size                                       */
c     /***************************************************************/
      real(nk) u(*),dxu(*),dyu(*),dzu(*)
      real(nk) tmp(*),tmpyz(*)
      real(nk) scal
      character*1 what
c
c-----|---------------------------------------------------------------|
c     /* 1st derivative in 'x': scal*d(u)/dx -> dxu  */
      call derivx(u,dxu,1,scal,what)
c
c     /* pass 'u' to f-c-p space -> u */
      call transyz(u,'y',-1) !fis->cheb
c
c     /* 1st derivative in 'y': scal*d(u)/dy -> tmp */
      call derivyz(u,tmp,tmpyz,'y',1,scal,'n')
c     
c     /* pass 'dyu' back to f-p-p space -> tmp */
      call transyz(tmp,'y',+1)    !cheb->fis
c
      if(what.eq.'n')then
         call dcopy(nsfour,tmp,1,dyu,1)
      else
c     /* add to 'dyu' -> ++dyu */
         call daxpy(nsfour,f1o1,tmp,1,dyu,1)
      endif
c     
c     /* pass 'u' back to f-p-p space -> u */
      call transyz(u,'y',+1)    !cheb->fis
c     
c     /* pass 'u' to f-p-c space -> u */
      call transyz(u,'z',-1)    !fis->cheb
c     
c     /* 1st derivative in 'z': scal*d(u)/dz -> tmp */
      call derivyz(u,tmp,tmpyz,'z',1,scal,'n')
c     
c     /* pass 'dzu' back to f-p-p space -> tmp */
      call transyz(tmp,'z',+1)    !cheb->fis
c
      if(what.eq.'n')then
         call dcopy(nsfour,tmp,1,dzu,1)
      else
c     /* add to 'dzu' -> ++dzu */
         call daxpy(nsfour,f1o1,tmp,1,dzu,1)
      endif
c     
c     /* pass 'u' back to f-p-p space -> u */
      call transyz(u,'z',+1)    !cheb->fis
c     
      return
      end subroutine grad3d
c-----|---------------------------------------------------------------|
      subroutine div3d(u,v,w,div,tmp,tmpyz,scal)
      use ctes
      use numbers
      implicit none
c
c     /***************************************************************/
c     /* computes div.([u,v,w]) of a field 'u,v,w' which is re-      */
c     /* presented by fourier expansion in x and chebyshev in y,z,   */
c     /* pseudo-spectrally. , X-CUT                                  */
c     /*                                                             */
c     /* input                                                       */
c     /* u,v,w(nsfou)         vector field in fou-phys-phys space    */
c     /* tmp(nsfour)           fourier-size 3D work array            */
c     /* tmpyz(ntmp)          1D work array: ntmp=max(my,mz)         */
c     /* scal                 scalar factor                          */
c     /*                                                             */
c     /* output                                                      */
c     /* div(nsfou)           div.(u) in f-p-p space                 */
c     /* u(nsfou)             field in fou-phys-phys space           */
c     /*                      [preserved]                            */
c     /*                                                             */
c     /* note: 'real*8 u(nsfour)' and 'complex*16 u(nsfou)' have     */
c     /*       equivalent size                                       */
c     /***************************************************************/
      include "mpif.h"
      real(nk) u(*),v(*),w(*),div(*)
      real(nk) tmp(*),tmpyz(*)
      real(nk) scal
      integer nlocalc,nlocalr,ierr
c
c-----|---------------------------------------------------------------|
      nlocalc=my*mz*mxp !(ie-ib+1)   !size in complex
      nlocalr=2*nlocalc         !size in real
c
c     /* 1st derivative in 'x': scal*ad(u)/dx -> div  */
      call derivx(u,div,1,scal,'n')
c      write(*,*) myid,'drivx done'
c      write(*,*) myid,"after div3d div=",div(1*(my*mz)+20*my+20)
c      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c
c     /* pass 'v' to f-c-p space -> v */
      call transyz(v,'y',-1) !fis->cheb
c      write(*,*) myid,'tranyz done'
c
c     /* 1st derivative in 'y': scal*d(v)/dy -> tmp */
      call derivyz(v,tmp,tmpyz,'y',1,scal,'n') 
c      write(*,*) myid,'derivyz done'
c     
c     /* pass 'dyv' back to f-p-p space -> tmp */
      call transyz(tmp,'y',+1)    !cheb->fis ! OK
c      write(*,*) myid,'transyz done'
c
c     /* add to 'div' -> ++div */
      call daxpy(nlocalr,f1o1,tmp,1,div,1)
c     
c     /* pass 'v' back to f-p-p space -> v */
      call transyz(v,'y',+1)    !cheb->fis
c     
c     /* pass 'w' to f-p-c space -> w */
      call transyz(w,'z',-1)    !fis->cheb
c     
c     /* 1st derivative in 'z': scal*d(w)/dz -> tmp */
      call derivyz(w,tmp,tmpyz,'z',1,scal,'n')
c     
c     /* pass 'dzw' back to f-p-p space -> tmp */
      call transyz(tmp,'z',+1)    !cheb->fis
c
c     /* add to 'div' -> ++div */
      call daxpy(nlocalr,f1o1,tmp,1,div,1)
c     
c     /* pass 'w' back to f-p-p space -> w */
      call transyz(w,'z',+1)    !cheb->fis
c     
      return
      end subroutine div3d
c-----|---------------------------------------------------------------|
      subroutine update_uvwp(u,v,w,p,phi,rhsphi,tmp,tmpyz,akk,bkk)
      use ctes
      use numbers
      use running
      implicit none
c
c     /***************************************************************/
c     /* updates velocities u,v,w and pressure in a fractional-step  */
c     /* method a la rai & moin. the fields are                      */
c     /* represented by fourier expansion in x and chebyshev in y,z, */
c     /* pseudo-spectrally.                                          */
c     /*                                                             */
c     /* BCOND-PHI: also assigns the boundary conditions for the     */
c     /*            tangential velocity for the following predictor  */
c     /*            step for obtaining O(dt^2) slip-error.           */
c     /*            (common block "/phibound/")                      */
c     /*                                                             */
c     /* input                                                       */
c     /* u,v,w(nsfou)         predicted velocity field               */
c     /*                        in fou-phys-phys space               */
c     /* p(nsfou)             previous pressure field (p^{k-1}) "    */
c     /* phi(nsfou)           pseudo-pressure "                      */
c     /* rhsphi(nsfou)        r.h.s. of pseudo-pressure "            */
c     /* tmp(nsfour)           fourier-size 3D work array            */
c     /* tmpyz(ntmp)          1D work array: ntmp=max(my,mz)         */
c     /* fnu                  kinematic viscosity                    */
c     /* akk                  \alpha_k - r-k coefficient             */
c     /* bkk                  \beta_k - r-k coefficient             */
c     /* dt                   time step                              */
c     /*                                                             */
c     /* output                                                      */
c     /* u,v,w,p(nsfou)       updated velocity & pressure fields     */
c     /*                        in fou-phys-phys space               */
c     /*                                                             */
c     /* note: 'real*8 u(nsfour)' and 'complex*16 u(nsfou)' have     */
c     /*       equivalent size                                       */
c     /* note: 'phi' is destroyed; thereby we save one chebyshev     */
c     /*       transform at the end                                  */
c     /* ------------------------------------------------------------*/
c     /* MOD 21/05/2004 filter is applied to phi, not p.             */
c     /***************************************************************/
      real(nk) u(*),v(*),w(*),p(*),phi(*),rhsphi(*)
      real(nk) tmp(*),tmpyz(*)
      real(nk) scal,akk,bkk
      integer i
c     
c-----|---------------------------------------------------------------|
#ifdef FILTER
c     /* filter the spurious modes from pseudo-pressure */
      call filter_p(phi)
#endif /* FILTER */
c
c     /* update pressure: p^k=p^{k-1}+phi^{k}-alpha_k*dt*nu*rhsphi   */
      scal=(akk+bkk)*dt*fnu
      call laplacian3d(phi,rhsphi,tmp,tmpyz,scal,'n')
      do i=1,nsfour
         p(i)=p(i)+phi(i)-rhsphi(i)
      enddo
#ifdef FILTER
      call filter_p(p)
#endif /* FILTER */
c
c     /* update velocity: u^k=u^\ast-2*alpha_k*dt*grad(phi^k)        */
c      scal=-f2o1*dt*akk ! for Crank--Nicolson
      scal=-dt*(akk+bkk)
c
c     /* 1st derivative in 'x': scal*d(phi)/dx -> ++u  */
      call derivx(phi,u,1,scal,'+')
c
c     /* pass 'phi' to f-c-p space -> phi */
      call transyz(phi,'y',-1) !fis->cheb
c
c     /* 1st derivative in 'y': scal*d(phi)/dy -> tmp */
      call derivyz(phi,tmp,tmpyz,'y',1,scal,'n')
c     
c     /* pass 'tmp' back to f-p-p space -> tmp */
      call transyz(tmp,'y',+1)    !cheb->fis
c     
c     /* pass 'phi' back to f-p-p space -> phi */
      call transyz(phi,'y',+1)    !cheb->fis
c
c     /* add to 'v' -> ++v */
      call daxpy(nsfour,f1o1,tmp,1,v,1)
c     
c     /* pass 'phi' to f-p-c space -> phi */
      call transyz(phi,'z',-1)    !fis->cheb
c     
c     /* 1st derivative in 'z': scal*d(phi)/dz -> tmp */
      call derivyz(phi,tmp,tmpyz,'z',1,scal,'n')
c     
c     /* pass 'tmp' back to f-p-p space -> tmp */
      call transyz(tmp,'z',+1)    !cheb->fis
c
c     /* add to 'w' -> ++w */
      call daxpy(nsfour,f1o1,tmp,1,w,1)
c
c     /* note: the following 'comments' save one cheb transform, but */
c     /*       phi is not preserved                                  */
c$$$c     /* pass 'phi' back to f-p-p space -> phi */
c$$$      call transyz(phi,'z',+1)    !cheb->fis
c     
      return
      end subroutine update_uvwp
c-----|---------------------------------------------------------------|
      subroutine determine_step(u,v,w)
      use ctes
      use numbers
      use running
      implicit none
c
c     /***************************************************************/
c     /* computes the time step by a convective restriction  (cfl)   */
c     /*          Z-cut                                              */
c     /* input                                                       */
c     /* u,v,w(nsfis)    velocity components in p-p-p space, Z-cut   */
c     /* cfl             value of the cfl number                     */
c     /*                                                             */
c     /* output                                                      */
c     /* dt              maximum time step                           */
c     /***************************************************************/
      include "mpif.h"
      real(nk) u(mgalx+2,my,kb:ke),v(mgalx+2,my,kb:ke),
     &     w(mgalx+2,my,kb:ke)
      real(nk) reigmx,hxalp,remx,hy,hz
      character*1 chi
      integer i,j,k,ierr
c-----|---------------------------------------------------------------|
      reigmx=f0o1
c
c     /* minimum time step in fourier (x) direction */
      hxalp=alp*dfloat(mgalx)/(f2o1*pi)
      do k=kb,ke
         do j=1,my
            do i=1,mgalx
               reigmx=max(reigmx,dabs(u(i,j,k)*hxalp))
            enddo
         enddo
      enddo
c
c     /* minimum time step in chebyshev (y) direction */
      do k=kb,ke
         do j=1,my
            hy=dy(j)
            do i=1,mgalx
               reigmx=max(reigmx,dabs(v(i,j,k)/hy))
            enddo
         enddo
      enddo
c
c     /* minimum time step in chebyshev (z) direction */
      do k=kb,ke
         hz=dz(k)
         do j=1,my
            do i=1,mgalx
               reigmx=max(reigmx,dabs(w(i,j,k)/hz))
            enddo
         enddo
      enddo
      call MPI_REDUCE(reigmx,remx,1,MPI_DOUBLE_PRECISION,
     &     MPI_MAX,master,MPI_COMM_WORLD,ierr)    
c
      if(remx.le.dsmall)then
         dt=dtfixed
      else
         dt=cfl/remx
      endif
c      
      if(myid.eq.master)
     $     write(*,*)'current time step: ',dt,' (cfl=',cfl,')'
c
      return
      end subroutine determine_step
c---------------------------------------------------------------------|
      subroutine wallstressx(u,tmpyz,stress4)
      use ctes
      use numbers
      implicit none 
c
c     /***************************************************************/
c     /* computes the streamwise-mean wallstress of a velocity field */
c     /* which is represented by fourier expansion in x and che-     */
c     /* byshev in y,z,                                              */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfou)             field in fou-phys-phys space           */
c     /* tmpyz(ntmp)          2D work array: ntmp=my*mz              */
c     /*                                                             */
c     /* output                                                      */
c     /* stress4              wallstress integrated over the         */
c     /*                        each one of the 4 boundaries         */
c     /*                                                             */
c     /***************************************************************/
      real(nk) u(*)
      real(nk) tmpyz(*),stress4(4)
c
c-----|---------------------------------------------------------------|
      if(myid.ne.0) return
c     /* copy real part of fourier-0-mode to tmpyz */
      call dcopy(my*mz,u,2,tmpyz,1)
c
c     /* pass the plane to cheb-cheb space */
      call transyz_plane(tmpyz,-1) !fis -> cheb
c
c     /* evaluate the wall-stress integral */
      call cheb2_wallstress_circum(tmpyz,1,aspect,stress4,my,mz)
c
      return
      end subroutine wallstressx
c---------------------------------------------------------------------|
      subroutine wallstressx_fou(u,tmpyz,foud,dudns_int,dudnst,stress4x)
      use ctes
      use numbers
      use rfttmp
      implicit none 
c
c     /***************************************************************/
c     /* computes the streamwise wallstress of a velocity field      */
c     /* which is represented by fourier expansion in x and che-     */
c     /* byshev in y,z,                                              */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfou)             field in fou-phys-phys space (parallel)*/
c     /* tmpyz(ntmp)          2D work array: ntmp=my*mz              */
c     /* foud(mgalx+2)         tmp array for FFT                      */
c     /* output                                                      */
c     /* dudns_int(2,0:mx1,4) tmp array to compute dudnx for each    */
c     /*                      fourier mode                           */      
c     /* dudnxt(0:mx1,4)      wallstress in fourier space (          */
c     /*                      only master)                           */      
c     /* stress4x(mgalx,1:4)   wallstress integrated over the        */
c     /*                       each one of the 4 boundaries          */
c     /*                       for each plane                        */
c     /*                                                             */
c     /***************************************************************/
      include "mpif.h"
      real(nk) u(2,my,mz,ib:ie)
      real(nk) tmpyz(*)
      complex*16 foud(0:mgx)
      real(nk) dudns_int(2,0:mx1,4)
      complex*16 dudnst(0:mx1,4) ! for master reduction
      real(nk) sum4(4)
      real(nk) stress4x(mgalx+2,4)
      integer imode,iw,ierr
c-----|---------------------------------------------------------------|
c      if(myid.ne.0) ! do not return
c     /* initialize array */ 
      call dcopy(2*(mx1+1)*4,f0o1,0,dudnst,1) 
c      
      do imode=ib,ie
c     /* copy real part of fourier-mode to tmpyz */
         call dcopy(my*mz,u(1,1,1,imode),2,tmpyz,1)
c
c     /* pass the plane to cheb-cheb space */
         call transyz_plane(tmpyz,-1) !fis -> cheb
c
c     /* evaluate the wall-stress integral */
c$$$      call cheb2_wallstress(tmpyz,1,aspect,stress,my,mz)
         call cheb2_wallstress_circum(tmpyz,1,aspect,
     $                                sum4,my,mz)
c     /* restore sum4 to real part of dudns_int complex array */
         do iw=1,4
            dudns_int(1,imode,iw)=sum4(iw)
         enddo
         if (imode.ne.0) then
c     /* copy imagnary part of fourier-mode to tmpyz */
            call dcopy(my*mz,u(2,1,1,imode),2,tmpyz,1)
c
c     /* pass the plane to cheb-cheb space */
            call transyz_plane(tmpyz,-1) !fis -> cheb
c
c     /* evaluate the wall-stress integral */
            call cheb2_wallstress_circum(tmpyz,1,aspect,
     $                                   sum4,my,mz)
c     /* restore sum4 to imarginay part of dudns_int complex array */
            do iw=1,4
               dudns_int(2,imode,iw)=sum4(iw)
            enddo
         else
            do iw=1,4
               dudns_int(2,imode,iw)=f0o1
            enddo
         endif
      enddo
c
c    /* master (ib=0) gets dudns_int and back to physical space */ 
      call MPI_REDUCE(dudns_int,dudnst,(mx1+1)*4*2,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,master,MPI_COMM_WORLD,ierr)
c
c    /* inverse fft to stress4x */
      if (myid.eq.master) then
         do iw=1,4
            call dcopy(mgalx+2,f0o1,0,foud(0),1)
C     /* zero-fill & copy to contiguous -> aux */
            call zcopy(mx1+1,dudnst(0,iw),1,foud(0),1)
            call zcopy(mgx-mx1,zero,0,foud(mx1+1),1)
            !call rft(foud(0),mgalx+2,1,1)
            ! this should be checked 2020/02/25
            call dfftw_execute_dft_c2r(plan_fft_backward,
     &                                 foud(0),foud(0))
c     /* copy back into 'phys' */
            call dcopy(mgalx+2,foud(0),1,stress4x(1,iw),1)
         enddo
      end if

      return
      end subroutine wallstressx_fou
c-----|---------------------------------------------------------------|
      subroutine vorticity(u,v,w,vort,tmp,iflag)
      use ctes
      use numbers
      implicit none 
c     /***************************************************************/
c     /* input                                                       */
c     /* compute a component of vorticity in fourier space           */
c     /* u,v,w(nsfou)         velocity components in fou-phys-phys   */
c     /*                      space                                  */
c     /* iflag                '1','2','3' - component of the vort.   */
c     /*                                                             */
c     /* output                                                      */
c     /* vort(nsfou)         vorticity component in fourier space    */
c     /*                                                             */
c     /* aux:                                                        */
c     /* tmp(nsfou)          work array                              */
c     /*                                                             */
c     /***************************************************************/
      real(nk) u(*),v(*),w(*),vort(*),tmp(*)
      integer iflag
      real(nk) dummy
c-----|---------------------------------------------------------------|
      if(iflag.eq.1)then        !compute: om_x=w,y-v,z
c     /* pass 'w' to f-c-p space -> w */
         call transyz(w,'y',-1) !fis->cheb
c     
c     /* derivative in 'y': d(w)/dy -> vort */
         call derivyz(w,vort,dummy,'y',1,f1o1,'n')
c
c     /* pass 'vort' back to f-p-p space -> vort */
         call transyz(vort,'y',+1) !cheb->fis
c     
c     /* pass 'w' back to f-p-p space -> w */
         call transyz(w,'y',+1) !cheb->fis
c     
c     /* pass 'v' to f-p-c space -> v */
         call transyz(v,'z',-1) !fis->cheb
c     
c     /* derivative in 'z': d(v)/dz -> tmp */
         call derivyz(v,tmp,dummy,'z',1,f1o1,'n')
c
c     /* pass 'tmp' back to f-p-p space -> tmp */
         call transyz(tmp,'z',+1) !cheb->fis
c
c     /* add to 'vort' -> --vort */
         call daxpy(nsfour,-f1o1,tmp,1,vort,1)
c     
c     /* pass 'v' back to f-p-p space -> v */
         call transyz(v,'z',+1) !cheb->fis
c         
      else if(iflag.eq.2)then    !compute: om_y=u,z-w,x
c******* (u,z) D13
c     /* pass 'u' to f-p-c space -> u */
         call transyz(u,'z',-1) !fis->cheb
c     /* derivative in 'z': d(u)/dz (f-p-c) -> vort (f-p-p) */
         call derivyz(u,vort,dummy,'z',1,f1o1,'n')
         call transyz(vort,'z',+1) !cheb->fis
c         
c     /* pass 'u' back to f-p-p space -> u */
         call transyz(u,'z',+1) !cheb->fis         
c******* (w,x) D31    
c     /* derivative in 'x': d(w)/dx (f-p-p)  */
         call derivx(w,tmp,1,f1o1,'n')
c******* om_y=u,z-w,x
c     /* add to 'vort' ->  du/dz - (dw/dx)*/
         call daxpy(nsfour,-f1o1,tmp,1,vort,1)
c
c         write(*,*)'vorticity: oemga_y implemented by sekimoto'
c
      else if(iflag.eq.3)then    !compute: om_z=v,x-u,y
c******* (v,x) D21
c     /* derivative in 'x': d(v)/dx (f-p-p)  */
         call derivx(v,vort,1,f1o1,'n')
c******* (u,y) D12    
c     /* pass 'u' to f-c-p space -> y */
         call transyz(u,'y',-1) !fis->cheb
c     /* derivative in 'y': d(u)/dy (f-c-p) -> tmp(f-p-p) */
         call derivyz(u,tmp,dummy,'y',1,f1o1,'n')
         call transyz(tmp,'y',+1) !cheb->fis

c     /* pass 'u' back to f-p-p space -> u */
         call transyz(u,'y',+1) !cheb->fis
c******* om_z=v,x-u,y
c     /* add to 'vort' ->  dw/dy - (dv/dz)*/
         call daxpy(nsfour,-f1o1,tmp,1,vort,1)
c      
c         write(*,*)'vorticity: omega_z implemented by sekimoto'
c              
      else
         write(*,*)'vorticity: this component not implemented!',
     $        iflag
         stop
      endif
c     
      return
      end subroutine vorticity
c-----|---------------------------------------------------------------|
      subroutine energy_fou(u,energy,energy2)
      use ctes
      use numbers
      implicit none 
c
c     /***************************************************************/
c     /* computes the rms energy of a velocity field which is re-    */
c     /* presented by fourier expansion in x and chebyshev in y,z,   */
c     /* X-cut                                                       */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfou)             field in fou-phys-phys space X-cut     */
c     /*                                                             */
c     /* output                                                      */
c     /* energy(non-zero mode)  energy defined as follows:           */
c     /* energy2(with zero mode)  energy defined as follows:         */
c     /*                                                             */
c     /***************************************************************/
      include "mpif.h"
      complex*16 u(my,mz,ib:ie)
      real(nk) energy, energy2, dyy, dzz, ueprim, area, ener, ueprim2
      integer i,j,k,ierr
c-----|---------------------------------------------------------------|
      energy=f0o1
      energy2=f0o1
c     
      do k=1,mz
         dzz=dz(k)
         do j=1,my
            dyy=dy(j)
            ueprim=f0o1
c
c     /* non-zero fourier modes: */
            do i=max(1,ib),ie
               ueprim=ueprim+u(j,k,i)*dconjg(u(j,k,i))*f2o1
            enddo
c
c     /* sum the energy */
            energy=energy+ueprim*dyy*dzz
         enddo
      enddo
      if(ib.eq.0)then
         do k=1,mz
            dzz=dz(k)
            do j=1,my
               dyy=dy(j)
c     
c     /* 0-fourier mode: */
               i=0
               ueprim2=u(j,k,i)*dconjg(u(j,k,i))
               energy2=energy2+ueprim2*dyy*dzz
            enddo
         enddo
      endif
c
      call MPI_REDUCE(energy,ener,1,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,master,MPI_COMM_WORLD,ierr)     
c
      area=f4o1*aspect
      energy=dsqrt(ener/area)
      if(myid.eq.0) then
         energy2=ener+energy2 
         energy2=dsqrt(energy2/area)
      endif
c
      return
      end subroutine energy_fou
c-----|---------------------------------------------------------------|
      subroutine energy_fou_zero(u,energy)
      use ctes
      use numbers
      implicit none 
c
c     /***************************************************************/
c     /* computes the rms energy of the streamwise average of a      */
c     /* velocity field which is represented by fourier expansion    */
c     /* in x and chebyshev in y,z,                                  */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfou)             field in fou-phys-phys space           */
c     /*                                                             */
c     /* output                                                      */
c     /* energy               energy defined as follows:             */
c     /*                                                             */
c     /***************************************************************/
      complex*16 u(my,mz,ib:ie)
      real(nk) energy,dzz,dyy,ueprim,area
      integer i,j,k,i0
c-----|---------------------------------------------------------------|
      energy=f0o1
      if(myid.ne.0) return
c     
      i0=0
      do k=1,mz
         dzz=dz(k)
         do j=1,my
            dyy=dy(j)
c     /* 0-fourier mode: */
            ueprim=u(j,k,i0)*dconjg(u(j,k,i0))
c
c     /* sum the energy */
            energy=energy+ueprim*dyy*dzz
         enddo
      enddo
c
      area=f4o1*aspect
      energy=dsqrt(energy/area)
c
      return
      end subroutine energy_fou_zero
c-----|---------------------------------------------------------------|
      subroutine energy_fou_zero_sector(u,energy)
      use ctes
      use numbers
      implicit none 
c
c     /***************************************************************/
c     /* computes the rms energy of the streamwise average of a      */
c     /* velocity field which is represented by fourier expansion    */
c     /* in x and chebyshev in y,z,                                  */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfou)             field in fou-phys-phys space           */
c     /*                                                             */
c     /* output                                                      */
c     /* energy(1:5)          energy in each sector                  */
c     /*                      sector 1: y<z & y<-z                   */
c     /*                             2: y<z & y>-z                   */
c     /*                             3: y>z & y>-z                   */
c     /*                             4: y>z & y<-z                   */
c     /*                             5: else (equal sign)            */
c     /*            --------                                         */
c     /*        ^  |\  3   /|                                        */
c     /*        |  | \    / |                                        */
c     /*        |  |  \  /  |                                        */
c     /*        |  |4  \/ 2 |                                        */
c     /*        y  |   /\   |                                        */
c     /*           |  /  \  |                                        */
c     /*           | /    \ |                                        */
c     /*           |/  1   \|                                        */
c     /*            --------                                         */
c     /*            z---->                                           */
c     /*                                                             */
c     /* ATTENTION: the contribution from the diagonal is added to   */
c     /*            to some of the sectors!                          */
c     /*-------------------------------------------------------------*/
c     /* MOD: 25.9.2006 m.u. excluded contribution on diagonal; sum  */
c     /*                     to fifth component of "energy"          */
c     /***************************************************************/
      complex*16 u(my,mz,ib:ie)
      real(nk) energy(5)
      integer i,i0,j,k, icomp
      real(nk) dzz, zz, dyy, yy, ueprim, area
c
c-----|---------------------------------------------------------------|
      if(myid.ne.0) return
      do i=1,5
         energy(i)=f0o1
      enddo
c     
      i0=0
      do k=1,mz
         dzz=dz(k)
         zz=z(k)
         do j=1,my
            dyy=dy(j)
            yy=y(j)
c     /* 0-fourier mode: */
            ueprim=u(j,k,i0)*dconjg(u(j,k,i0))
c     
c     /* sum the energy for each sector */
            if(j.lt.k.and.j.lt.(mz+1-k))then
               icomp=1
            elseif(j.lt.k.and.j.gt.(mz+1-k))then
               icomp=2
            elseif(j.gt.k.and.j.gt.(mz+1-k))then
               icomp=3
            elseif(j.gt.k.and.j.lt.(mz+1-k))then
               icomp=4
            else                !j=k and/or j=(mz+1-k) => on diagonal
               icomp=5
            endif
            energy(icomp)=energy(icomp)+ueprim*dyy*dzz
         enddo
      enddo
c     
      area=f4o1*aspect
      do i=1,5
         energy(i)=(energy(i)/area)
      enddo
c     
      return
      end subroutine energy_fou_zero_sector
c-----|---------------------------------------------------------------|
      subroutine energy_fis_sector_zcut(ufou,tmpfou,tmpfis,tmpfis_ave,
     $                              tmpzxcut,tmpx,
     $                              energy,moving_ave)
      use ctes
      use numbers
      implicit none 
c     /***************************************************************/
c     /* computes the rms energy of a      */
c     /* velocity field which is represented by fourier expansion    */
c     /* in x and chebyshev in y,z,                                  */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfou)             field in fou-phys-phys space           */
c     /*                                                             */
c     /* output                                                      */
c     /* energy(1:5)          energy in each sector                  */
c     /*                      sector 1: y<z & y<-z                   */
c     /*                             2: y<z & y>-z                   */
c     /*                             3: y>z & y>-z                   */
c     /*                             4: y>z & y<-z                   */
c     /*                             5: else (equal sign)            */
c     /*            --------                                         */
c     /*        ^  |\  3   /|                                        */
c     /*        |  | \    / |                                        */
c     /*        |  |  \  /  |                                        */
c     /*        |  |4  \/ 2 |                                        */
c     /*        y  |   /\   |                                        */
c     /*           |  /  \  |                                        */
c     /*           | /    \ |                                        */
c     /*           |/  1   \|                                        */
c     /*            --------                                         */
c     /*            z---->                                           */
c     /*                                                             */
c     /* ATTENTION: the contribution from the diagonal is added to   */
c     /*            to some of the sectors!                          */
c     /***************************************************************/
      include 'mpif.h'
      complex(nk) ufou(my,mz,ib:ie)
      real(nk) tmpfou(nsfour)
      real(nk) tmpfis(mgalx+2,my,kb:ke)
      real(nk) tmpfis_ave(mgalx+2,my,kb:ke)
      real(nk) energy(mgalx,5)
      real(nk) ener_local(mgalx,5)
c
      complex(nk) tmpzxcut(my*mzp*mxp) !z-cut & x-cut pencil work array
c     /* 1d physical-size x work arrays: */
      real(nk) tmpx(mgalx+2)
      integer mspan,iix,iixs,iixe,i,k,j,icomp,isec,ierr
      real(nk) fac,dzz,zz,dyy,yy,area,moving_ave,ueprim
c
c-----|---------------------------------------------------------------|
      if(kb.eq.1) write(*,*) 'parallel-computing the energy in sectors'
      call dcopy(nsfour,ufou,1,tmpfou,1)
c     ! ufou is stored, so do not care for repair by shuffling! 
      call chx2z(tmpfou,tmpfis,tmpzxcut)
      call fourx(tmpfou,tmpfis,tmpx,+1) !fou->fis
      
      if (moving_ave.eq.1) then
         call dcopy(nsfis,f0o1,0,tmpfis_ave,1)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c         /* compute moving average for x over dLx/h=4\pi */
         mspan=81
         do k=kb,ke
            do j=1,my
               do iix=1,mgalx
                  iixs=iix-(mspan-1)/2
                  iixe=iix+(mspan-1)/2
                  do i=max(iixs,1),min(iixe,mgalx)
                     tmpfis_ave(iix,j,k) = tmpfis_ave(iix,j,k)
     $                                    + tmpfis(i,j,k)
                  end do
                  if (iixs.le.0) then
                     do i=mgalx+iixs,mgalx
                        tmpfis_ave(iix,j,k) = tmpfis_ave(iix,j,k)
     $                                       + tmpfis(i,j,k)
                     end do
                  end if
                  if (iixe.ge.mgalx+1) then
                     do i=1,mgalx-iixe
                        tmpfis_ave(iix,j,k) = tmpfis_ave(iix,j,k)
     $                                       + tmpfis(i,j,k)
                     end do
                  end if
               end do
            end do
         end do
         fac=f1o1/dfloat(mspan)
         call dscal2(nsfis,tmpfis_ave,1,fac,tmpfis,1)
         if(kb.eq.1) write(*,*) ' moving averaged '
      endif

      call dcopy(mgalx*5,f0o1,0,ener_local,1)
c     /* sum the energy for each sector */
      do k=kb,ke
         dzz=dz(k)
         zz=z(k)
         do j=1,my
            dyy=dy(j)
            yy=y(j)
            if(j.lt.k.and.j.lt.(mz+1-k))then
               icomp=1
            elseif(j.lt.k.and.j.gt.(mz+1-k))then
               icomp=2
            elseif(j.gt.k.and.j.gt.(mz+1-k))then
               icomp=3
            elseif(j.gt.k.and.j.lt.(mz+1-k))then
               icomp=4
            else                !j=k and/or j=(mz+1-k) => on diagonal
               icomp=5
            endif
            do i=1,mgalx
               ueprim=tmpfis(i,j,k)*tmpfis(i,j,k)
               ener_local(i,icomp)=ener_local(i,icomp)+ueprim*dyy*dzz
            enddo
         enddo
      enddo
c
      if (myid.eq.master) call dcopy(mgalx*5,f0o1,0,energy,1)
c    /* master (ib=0) gather ener_local to energy */ 
      call MPI_REDUCE(ener_local,energy,mgalx*5,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,master,MPI_COMM_WORLD,ierr)
cc   
      if (myid.eq.master) then
         area=f4o1*aspect
         do isec=1,5
            do i=1,mgalx
               energy(i,isec)=(energy(i,isec)/area)
            enddo
         enddo
      end if
c
      return
      end subroutine energy_fis_sector_zcut
c-----|---------------------------------------------------------------|
      subroutine energy_fis_xslice(ufis2,energy)
      use ctes
      use numbers
      implicit none 
c
c     /***************************************************************/
c     /* computes the rms energy of a velocity field which is re-    */
c     /* presented by fourier expansion in x and chebyshev in y,z,   */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfis)             field in phys-phys-phys space          */
c     /*                                                             */
c     /* output                                                      */
c     /* energy               energy defined as follows:             */
c     /*                                                             */
c     /***************************************************************/
      include 'mpif.h'
      real(nk) ufis2(my,mz,ipb:ipe)
      real(nk) energy(mgalx+2)
      real(nk) energy_local(mgalx+2)
      integer i,j,k,ierr
      real(nk) zz,dzz,yy,dyy,area
c
c-----|---------------------------------------------------------------|
c      write(*,*) ' energy_fis_xslice is now implementing!'
c     sekimoto 2010/1/20 '
c     
      if(ipb.eq.1) then
         do i=1,mgalx
            energy(i)=0.d0
         enddo
      endif
      do i=ipb,ipe
         energy_local(i)=f0o1
         do k=1,mz
            dzz=dz(k)
            do j=1,my
               dyy=dy(j)
               energy_local(i)=
     $              energy_local(i)+ufis2(j,k,i)**2*dyy*dzz
            enddo
         enddo
      enddo
c
      area=f4o1*aspect
      energy_local=dsqrt(energy_local/area)
c
c     /* send data to master processor*/
      call MPI_REDUCE(energy_local,energy,mgalx+2,
     $     MPI_DOUBLE_PRECISION,
     $     MPI_SUM,master,MPI_COMM_WORLD,ierr)
c
      return
      end subroutine energy_fis_xslice
c-----|---------------------------------------------------------------|
      subroutine energy_fis(u,energy)
      use ctes
      use numbers
      implicit none 
c
c     /***************************************************************/
c     /* computes the rms energy of a velocity field which is re-    */
c     /* presented by fourier expansion in x and chebyshev in y,z,   */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfis)             field in phys-phys-phys space          */
c     /*                                                             */
c     /* output                                                      */
c     /* energy               energy defined as follows:             */
c     /*                                                             */
c     /***************************************************************/
      real(nk) u(my,mz,mgalx+2)
      real(nk) energy,area,dzz,dyy
      integer i,j,k
c-----|---------------------------------------------------------------|
      write(*,*) 'energy_fis is not parallelized! '
      stop
      energy=f0o1
c     
      do i=1,mgalx
         do k=1,mz
            dzz=dz(k)
            do j=1,my
               dyy=dy(j)
               energy=energy+u(j,k,i)**2*dyy*dzz
            enddo
         enddo
      enddo
c
      area=f4o1*aspect*dfloat(mgalx)
      energy=dsqrt(energy/area)
c
      return
      end subroutine energy_fis
c-----|---------------------------------------------------------------|
      subroutine energy_fis_zero(u,energy)
      use ctes
      use numbers
      implicit none 
c     /***************************************************************/
c     /* computes the rms energy of the streamwise average of a      */
c     /* field in physical space                                     */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfis)             field in phys-phys-phys space          */
c     /*                                                             */
c     /* output                                                      */
c     /* energy               energy defined as follows:             */
c     /*                        \int(umean^2*dy*dz)                  */
c     /*                      where:                                 */
c     /*                        umean=\int(u*dx)/Lx                  */
c     /*                                                             */
c     /***************************************************************/
      real(nk) u(my,mz,mgalx+2)
      real(nk) energy, dzz, dyy, umean, area
      integer i,j,k
c-----|---------------------------------------------------------------|
      energy=f0o1
      write(*,*) 'energy_fis_zero is not parallelized! '
      stop
c     
      do k=1,mz
         dzz=dz(k)
         do j=1,my
            dyy=dy(j)
            umean=f0o1
            do i=1,mgalx
               umean=umean+u(j,k,i)
            enddo
            umean=umean/dfloat(mgalx)
            energy=energy+(umean)**2*dyy*dzz
         enddo
      enddo
c
      area=f4o1*aspect
      energy=dsqrt(energy/area)
c
      return
      end subroutine energy_fis_zero
c-----|---------------------------------------------------------------|
      subroutine energyabs_fis_zero_sector(u,energy)
      use ctes
      use numbers
      implicit none 
c
c     /***************************************************************/
c     /* computes the abs        of the streamwise average of a      */
c     /* field in physical space                                     */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfis)             field in phys-phys-phys space          */
c     /*                                                             */
c     /* output                                                      */
c     /* energy(1:5)          energy in each sector                  */
c     /*                        \int_Si(|umean|*dy*dz)               */
c     /*                                ^^^^^^^                      */
c     /*                      where:                                 */
c     /*                        umean=\int(u*dx)/Lx                  */
c     /*                      and where the sectors are defined as:  */
c     /*                            S1: y<z & y<-z                   */
c     /*                            S2: y<z & y>-z                   */
c     /*                            S3: y>z & y>-z                   */
c     /*                            S4: y>z & y<-z                   */
c     /*                            S5: else (equal sign)            */
c     /*            --------                                        */
c     /*        ^  |\  3   /|                                        */
c     /*        |  | \    / |                                        */
c     /*        |  |  \  /  |                                        */
c     /*        |  |4  \/ 2 |                                        */
c     /*        y  |   /\   |                                        */
c     /*           |  /  \  |                                        */
c     /*           | /    \ |                                        */
c     /*           |/  1   \|                                        */
c     /*            --------                                         */
c     /*            z---->                                           */
c     /*-------------------------------------------------------------*/
c     /* MOD: 25.9.2006 m.u. excluded contribution on diagonal; sum  */
c     /*                     to fifth component of "energy"          */
c     /***************************************************************/
      real(nk) u(my,mz,mgalx+2)
      real(nk) energy(5)
      real(nk) zz,dzz,yy,dyy,umean,area
      integer i,j,k,icomp
c
c-----|---------------------------------------------------------------|
      write(*,*) 'energyabs_fis_zero_sector is not parallelized! '
      stop
      do i=1,5
         energy(i)=f0o1
      enddo
c     
      do k=1,mz
         dzz=dz(k)
         zz=z(k)
         do j=1,my
            dyy=dy(j)
            yy=y(j)
            umean=f0o1
            do i=1,mgalx
               umean=umean+u(j,k,i)
            enddo
            umean=umean/dfloat(mgalx)
c
c     /* sum the energy for each sector */
            if(j.lt.k.and.j.lt.(mz+1-k))then
               icomp=1
            elseif(j.lt.k.and.j.gt.(mz+1-k))then
               icomp=2
            elseif(j.gt.k.and.j.gt.(mz+1-k))then
               icomp=3
            elseif(j.gt.k.and.j.lt.(mz+1-k))then
               icomp=4
            else                !j=k and/or j=(mz+1-k) => on diagonal
               icomp=5
            endif
            energy(icomp)=energy(icomp)+abs(umean)*dyy*dzz
         enddo
      enddo
c
      area=f4o1*aspect
      do i=1,5
         energy(i)=(energy(i)/area)
      enddo
c
      return
      end subroutine energyabs_fis_zero_sector
c-----|---------------------------------------------------------------|
      subroutine integral_abs_zero_sector(u,energy)
      use ctes
      use numbers
      implicit none
c
c     /***************************************************************/
c     /* computes the abs        of the streamwise average of a      */
c     /* field in fourier  space,  X_cut                             */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfou)             field in fou-phys-phys space, x-cut    */
c     /*                                                             */
c     /* output                                                      */
c     /* energy(1:5)          energy in each sector                  */
c     /*                        \int_Si(|umean|*dy*dz)               */
c     /*                                ^^^^^^^                      */
c     /*                      where:                                 */
c     /*                        umean=\int(u*dx)/Lx                  */
c     /*                      and where the sectors are defined as:  */
c     /*                            S1: y<z & y<-z                   */
c     /*                            S2: y<z & y>-z                   */
c     /*                            S3: y>z & y>-z                   */
c     /*                            S4: y>z & y<-z                   */
c     /*                            S5: else (equal sign)            */
c     /*            --------                                        */
c     /*        ^  |\  3   /|                                        */
c     /*        |  | \    / |                                        */
c     /*        |  |  \  /  |                                        */
c     /*        |  |4  \/ 2 |                                        */
c     /*        y  |   /\   |                                        */
c     /*           |  /  \  |                                        */
c     /*           | /    \ |                                        */
c     /*           |/  1   \|                                        */
c     /*            --------                                         */
c     /*            z---->                                           */
c     /*-------------------------------------------------------------*/
c     /* MOD: 25.9.2006 m.u. excluded contribution on diagonal; sum  */
c     /*                     to fifth component of "energy"          */
c     /***************************************************************/
      complex(nk) u(my,mz,ib:ie)
      real(nk) energy(5)
      real(nk) yy,dyy,zz,dzz,umean,area
      integer i,j,k,i0,icomp
c-----|---------------------------------------------------------------|
      if(myid.ne.0) return
      do i=1,5
         energy(i)=f0o1
      enddo
c     
      i0=0
      do k=1,mz
         dzz=dz(k)
         zz=z(k)
         do j=1,my
            dyy=dy(j)
            yy=y(j)
            umean=u(j,k,i0)
c     /* sum the energy for each sector */
            if(j.lt.k.and.j.lt.(mz+1-k))then
               icomp=1
            elseif(j.lt.k.and.j.gt.(mz+1-k))then
               icomp=2
            elseif(j.gt.k.and.j.gt.(mz+1-k))then
               icomp=3
            elseif(j.gt.k.and.j.lt.(mz+1-k))then
               icomp=4
            else                !j=k and/or j=(mz+1-k) => on diagonal
               icomp=5
            endif
            energy(icomp)=energy(icomp)+abs(umean)*dyy*dzz
         enddo
      enddo
c
      area=f4o1*aspect
      do i=1,5
         energy(i)=(energy(i)/area)
      enddo
c
      return
      end subroutine integral_abs_zero_sector
c-----|---------------------------------------------------------------|
      subroutine calc_energy_input_fou(u,p,energy)
      use ctes
      use numbers
      implicit none
c
c     /***************************************************************/
c     /* computes the energy input rate                              */
c     /* presented by fourier expansion in x and chebyshev in y,z,   */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfou)           u field in fou-phys-phys space (x-cut)   */
c     /* p(nsfou)           p field in fou-phys-phys space (x-cut)   */
c     /*                                                             */
c     /* output                                                      */
c     /* energy               energy input rate defined as follows:  */
c     /*                                                             */
c     /* energy input rate                                           */
c     /*  = -1/(dens*volume)\int{u \frac{dp}{dx}}dxdydz(fis)         */
c     /*  = -1/(dens*area)\int{ [u\frac{dp}{dx}]_0}dydz(fou zeromode)*/
c     /*                                                             */
c     /*  u\frac{dp}{dx}]_0} = sum_{l+m=0}{u^_l * m*p^_m}            */
c     /*           = sum_{m=1}^{mx1}{sqrt(-1)*m*dconjg(u_m)*p_m }    */
c     /*            + sum_{m=1}^{mx1}{-sqrt(-1)*m*u_m*dconjg(p_m) }  */
c     /*                                                             */
c     /***************************************************************/
      include "mpif.h"
      complex(nk) u(my,mz,ib:ie),p(my,mz,ib:ie)
      real(nk) tmp(my,mz)
      real(nk) energy,ener,area,dyy,dzz
      integer i,j,k,ierr
c
c-----|---------------------------------------------------------------|
      do k=1,mz
         do j=1,my
            tmp(j,k)=0.d0
c     /* non-zero fourier modes: */
            do i=max(1,ib),ie
               tmp(j,k)=tmp(j,k)+f2o1*alp*dfloat(i)*(
     $              - dreal(u(j,k,i))*dimag(p(j,k,i))
     $              + dimag(u(j,k,i))*dreal(p(j,k,i)))
            enddo
         enddo
      enddo
c
      energy=f0o1
c     
c     /* sum in y-z plane */
      do k=1,mz
         dzz=dz(k)
         do j=1,my
            dyy=dy(j)               
            energy=energy+tmp(j,k)*dyy*dzz
         enddo
      enddo
c
      call MPI_REDUCE(energy,ener,1,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,master,MPI_COMM_WORLD,ierr)     
c
      area=f4o1*aspect
      energy=dsqrt(ener/area)
c     /* root for the consistency with the old format */
c
      return
      end subroutine calc_energy_input_fou
c-----|---------------------------------------------------------------|
      subroutine statistics_fou(u,v,w,p,iflag)
      use ctes
      use numbers
      use running
      use statistics
      implicit none 
c
c     /***************************************************************/
c     /* accumulates the statistics. input is in fourier-physical.   */
c     /*   x-average um,vm,wm,pm. normal-stresses: uu,vv,ww,pp.      */
c     /*                                                             */
c     /* input                                                       */
c     /* u(nsfou)           u field in fou-phys-phys space,   x-cut  */
c     /* v(nsfou)           v field in fou-phys-phys space           */
c     /* w(nsfou)           w field in fou-phys-phys space           */
c     /* p(nsfou)           p field in fou-phys-phys space           */
c     /* time               present physical time                    */
c     /* iflag              '0' - initialize only                    */
c     /*                    [else] accumulate stats                  */
c     /*                                                             */
c     /* output                                                      */
c     /*                                                             */
c     /***************************************************************/
      complex*16 u(my,mz,ib:ie),v(my,mz,ib:ie),w(my,mz,ib:ie),
     $     p(my,mz,ib:ie)
      integer iflag,i,j,k,ispec
c-----|---------------------------------------------------------------|
c
      if(iflag.eq.-1) then
         ! only allocation
         allocate(um(my,mz),vm(my,mz),wm(my,mz),pm(my,mz))
         um=0.d0; vm=0.d0; wm=0.d0; pm=0.d0
         allocate(uub(my,mz),uvb(my,mz),uwb(my,mz), 
     &        vvb(my,mz),vwb(my,mz),wwb(my,mz),ppb(my,mz))
         uub=0.d0; uvb=0.d0; 
         uwb=0.d0; vvb=0.d0; 
         vwb=0.d0; wwb=0.d0; 
         ppb=0.d0

         allocate(uus(0:mx1,nspec),vvs(0:mx1,nspec),
     &        wws(0:mx1,nspec),pps(0:mx1,nspec))

         allocate(jinco(nspec),kinco(nspec))

         return

      elseif(iflag.eq.0) then
         if(myid.eq.master) write(*,*)'initializing statistics ...'
         ! allocated in initcr()
         nstat=0
         timef=time             !initial time
         call dcopy(my*mz,f0o1,0,um,1)
         call dcopy(my*mz,f0o1,0,vm,1)
         call dcopy(my*mz,f0o1,0,wm,1)
         call dcopy(my*mz,f0o1,0,pm,1)
c
         call dcopy(my*mz,f0o1,0,uub,1)
         call dcopy(my*mz,f0o1,0,uvb,1)
         call dcopy(my*mz,f0o1,0,uwb,1)
         call dcopy(my*mz,f0o1,0,vvb,1)
         call dcopy(my*mz,f0o1,0,vwb,1)
         call dcopy(my*mz,f0o1,0,wwb,1)
         call dcopy(my*mz,f0o1,0,ppb,1)
c
         call dcopy((mx1+1)*nspec,f0o1,0,uus,1)
         call dcopy((mx1+1)*nspec,f0o1,0,vvs,1)
         call dcopy((mx1+1)*nspec,f0o1,0,wws,1)
         call dcopy((mx1+1)*nspec,f0o1,0,pps,1)
         if(nspec.ne.7)then
            write(*,*)'nspec=7!'
            stop
         endif
         jinco(1)=my/16
         kinco(1)=mz/2
         jinco(2)=my/8
         kinco(2)=mz/2
         jinco(3)=my/4
         kinco(3)=mz/2
         jinco(4)=my/2
         kinco(4)=mz/2
         jinco(5)=my/16
         kinco(5)=mz/16
         jinco(6)=my/8
         kinco(6)=mz/8
         jinco(7)=my/4
         kinco(7)=mz/4
         
         return
      endif
c     /* increment the counter */
      nstat=nstat+1
c
c     /* 1st moments (x-average=0-fourier mode): um,vm,wm,pm */
      if(myid.eq.0) then
         call daxpy(my*mz,+f1o1,u,2,um,1)
         call daxpy(my*mz,+f1o1,v,2,vm,1)
         call daxpy(my*mz,+f1o1,w,2,wm,1)
         call daxpy(my*mz,+f1o1,p,2,pm,1)
      endif
c     
c     /* 2nd moments: normal stresses uu,vv,ww,pp [cross->physical space] */
c$$   /*   just only pp --> uub,vvb,wwb are computed in advection3d     */
c$$   /*   then we need not to use MPI_REDUCE in write_statistics     */
      do i=max(1,ib),ie         !non-zero fourier modes
         do k=1,mz
            do j=1,my               
c$$               uub(j,k)=uub(j,k)+u(j,k,i)*dconjg(u(j,k,i))*f2o1
c$$               vvb(j,k)=vvb(j,k)+v(j,k,i)*dconjg(v(j,k,i))*f2o1
c$$               wwb(j,k)=wwb(j,k)+w(j,k,i)*dconjg(w(j,k,i))*f2o1
               ppb(j,k)=ppb(j,k)+p(j,k,i)*dconjg(p(j,k,i))*f2o1
            enddo
         enddo
      enddo
c
c     /* zero fourier modes */
      if(myid.eq.0) then
         i=0
         do k=1,mz
            do j=1,my
c$$               uub(j,k)=uub(j,k)+u(j,k,i)*dconjg(u(j,k,i))
c$$               vvb(j,k)=vvb(j,k)+v(j,k,i)*dconjg(v(j,k,i))
c$$               wwb(j,k)=wwb(j,k)+w(j,k,i)*dconjg(w(j,k,i))
               ppb(j,k)=ppb(j,k)+p(j,k,i)*dconjg(p(j,k,i))
            enddo
         enddo
      endif
c
c     /* nspec streamwise spectra, at positions (jinco,kinco) */
c     /* DO NOT average !                                     */
      do ispec=1,nspec
         j=jinco(ispec)
         k=kinco(ispec)
         do i=max(1,ib),ie             !non-zero fourier modes
            uus(i,ispec)=uus(i,ispec)+
     $           (u(j,k,i)*dconjg(u(j,k,i)))*f2o1
            vvs(i,ispec)=vvs(i,ispec)+
     $           (v(j,k,i)*dconjg(v(j,k,i)))*f2o1
            wws(i,ispec)=wws(i,ispec)+
     $           (w(j,k,i)*dconjg(w(j,k,i)))*f2o1
            pps(i,ispec)=pps(i,ispec)+
     $           (p(j,k,i)*dconjg(p(j,k,i)))*f2o1
         enddo
      enddo
c
      if(myid.eq.0) then
         do ispec=1,nspec
            j=jinco(ispec)
            k=kinco(ispec)
            i=0                 !0-fourier mode
            uus(i,ispec)=uus(i,ispec)+(u(j,k,i)*dconjg(u(j,k,i)))
            vvs(i,ispec)=vvs(i,ispec)+(v(j,k,i)*dconjg(v(j,k,i)))
            wws(i,ispec)=wws(i,ispec)+(w(j,k,i)*dconjg(w(j,k,i)))
            pps(i,ispec)=pps(i,ispec)+(p(j,k,i)*dconjg(p(j,k,i)))
         enddo
      endif
c
      return
      end subroutine statistics_fou
c-----|---------------------------------------------------------------|
#ifdef TEMPERATURE
      subroutine statisticst_fou(t,iflag)
      use ctes
      use numbers
      use running
      use statistics
      implicit none 
c
c     /***************************************************************/
c     /* accumulates the statistics. input is in fourier-physical.   */
c     /*   x-average tm. normal-stresses: tt.                        */
c     /*                                                             */
c     /* input                                                       */
c     /* t(nsfou)           t field in fou-phys-phys space,   x-cut  */
c     /* time               present physical time                    */
c     /* iflag              '0' - initialize only                    */
c     /*                    [else] accumulate stats                  */
c     /*                                                             */
c     /* output                                                      */
c     /*                                                             */
c     /***************************************************************/
      complex*16 t(my,mz,ib:ie)
      integer iflag,i,j,k,ierr,ispec
c-----|---------------------------------------------------------------|
c
      if (iflag.eq.-1) then
         ! allocation only 
         allocate(tm(my,mz),ttb(my,mz),utb(my,mz),vtb(my,mz),wtb(my,mz)) 
         allocate(tts(0:mx1,nspec))

         return
      elseif(iflag.eq.0)then
         if(myid.eq.master) write(*,*)'initializing statistics_t...'
         call dcopy(my*mz,f0o1,0,tm,1)
c
         call dcopy(my*mz,f0o1,0,ttb,1)
         call dcopy(my*mz,f0o1,0,utb,1)
         call dcopy(my*mz,f0o1,0,vtb,1)
         call dcopy(my*mz,f0o1,0,wtb,1)
c
         call dcopy((mx1+1)*nspec,f0o1,0,tts,1)
         return
      endif
c
c     /* 1st moments (x-average=0-fourier mode): tm */
      if(myid.eq.0) then
         call daxpy(my*mz,+f1o1,t,2,tm,1)
      endif
c     
c     /* 2nd moments: normal stresses tt [cross->physical space] */
      do i=max(1,ib),ie         !non-zero fourier modes
         do k=1,mz
            do j=1,my
               ttb(j,k)=ttb(j,k)+t(j,k,i)*dconjg(t(j,k,i))*f2o1
            enddo
         enddo
      enddo
c
c     /* zero fourier modes */
      if(myid.eq.0) then
         i=0
         do k=1,mz
            do j=1,my
               ttb(j,k)=ttb(j,k)+t(j,k,i)*dconjg(t(j,k,i))
            enddo
         enddo
      endif
c
c     /* nspec streamwise spectra, at positions (jinco,kinco) */
c     /* DO NOT average !                                     */
      do ispec=1,nspec
         j=jinco(ispec)
         k=kinco(ispec)
         do i=max(1,ib),ie             !non-zero fourier modes
            tts(i,ispec)=tts(i,ispec)+
     $           (t(j,k,i)*dconjg(t(j,k,i)))*f2o1
         enddo
      enddo
c
      if(myid.eq.0) then
         do ispec=1,nspec
            j=jinco(ispec)
            k=kinco(ispec)
            i=0                 !0-fourier mode
            tts(i,ispec)=tts(i,ispec)+(t(j,k,i)*dconjg(t(j,k,i)))
         enddo
      endif
c
      return
      end subroutine statisticst_fou
#endif
c-----|---------------------------------------------------------------|
      subroutine statistics_omg2_fou(fou,iflag)
      use ctes
      use numbers
      use statistics
      implicit none 
c
c     /***************************************************************/
c     /* accumulates statistics in fou space (using only zero-mode). */
c     /* less cost version by sekimoto                               */
c     /*   2nd-moments: uu,vv,ww,uv,uw,vw; iflag=1,2,3,4,5,6         */
c     /*                                                             */
c     /* input                                                       */
c     /* fou(nsfour)   stresses field in fou-phys-phys space X-cut   */
c     /* output                                                      */
c     /*                                                             */
c     /***************************************************************/
      complex*16 fou(my,mz,ib:ie)
      integer iflag,i,j,k
c-----|---------------------------------------------------------------|
      if (iflag.eq.-1) then
         allocate(omx2(my,mz),omy2(my,mz),omz2(my,mz))
      elseif(iflag.eq.0)then
         if(myid.eq.0) then
            write(*,*)'initializing omega related statistics ...'
            call dcopy(my*mz,f0o1,0,omx2,1)
            call dcopy(my*mz,f0o1,0,omy2,1)
            call dcopy(my*mz,f0o1,0,omz2,1)
         else
            return
         endif
c     /* 2nd moments: omx,omy,omz */
      elseif(iflag.eq.1) then
c     /* zero fourier modes */
         if(myid.eq.0) then
            i=0
            do k=1,mz
               do j=1,my
                  omx2(j,k)=omx2(j,k)+
     $                 fou(j,k,i)*dconjg(fou(j,k,i))
               enddo
            enddo
         endif
         do i=max(1,ib),ie !non-zero fourier modes
            do k=1,mz
               do j=1,my
                  omx2(j,k)=omx2(j,k)+
     $                 fou(j,k,i)*dconjg(fou(j,k,i))*f2o1
               enddo
            enddo
         enddo
c     
      else if(iflag.eq.2) then
c     /* zero fourier modes */
         if(myid.eq.0) then 
            i=0
            do k=1,mz
               do j=1,my
                  omy2(j,k)=omy2(j,k)+
     $                 fou(j,k,i)*dconjg(fou(j,k,i))
               enddo
            enddo
         endif
         do i=max(1,ib),ie             !non-zero fourier modes
            do k=1,mz
               do j=1,my
                  omy2(j,k)=omy2(j,k)+
     $                 fou(j,k,i)*dconjg(fou(j,k,i))*f2o1
               enddo
            enddo
         enddo
c
      else if(iflag.eq.3) then
c     /* zero fourier modes */
         if(myid.eq.0) then
            i=0
            do k=1,mz
               do j=1,my
                  omz2(j,k)=omz2(j,k)+
     $                 fou(j,k,i)*dconjg(fou(j,k,i))
               enddo
            enddo
         endif
         do i=max(1,ib),ie             !non-zero fourier modes
            do k=1,mz
               do j=1,my
                  omz2(j,k)=omz2(j,k)+
     $                 fou(j,k,i)*dconjg(fou(j,k,i))*f2o1
               enddo
            enddo
         enddo
c     
      else 
         write(*,*) 'statistics_fou2: this flag is wrong!'
      endif
c     
      return
      end subroutine statistics_omg2_fou
c-----|---------------------------------------------------------------|
      subroutine statistics_fou2(fou,iflag)
      use ctes
      use numbers
      use statistics
      implicit none 
c
c     /***************************************************************/
c     /* accumulates statistics in fou space (using only zero-mode). */
c     /* less cost version by sekimoto                               */
c     /*   2nd-moments: uu,vv,ww,uv,uw,vw; iflag=1,2,3,4,5,6         */
c     /*                                                             */
c     /* input                                                       */
c     /* fou(nsfour)   stresses field in fou-phys-phys space X-cut   */
c     /* output                                                      */
c     /*                                                             */
c     /***************************************************************/
      complex(nk) fou(my,mz,ib:ie)
      integer i,j,k,iflag
c-----|---------------------------------------------------------------|
c
      if(myid.ne.0) return
c
c     /* 2nd moments: uu,vv,ww */
      if(iflag.eq.1) then
         call daxpy(my*mz,+f1o1,fou,2,uub,1)
c
      else if(iflag.eq.2) then
         call daxpy(my*mz,+f1o1,fou,2,vvb,1)
c
      else if(iflag.eq.3) then
         call daxpy(my*mz,+f1o1,fou,2,wwb,1)
c
c     /* 2nd moments, cross-stresses: uv,uw,vw */
      else if(iflag.eq.4) then
         call daxpy(my*mz,+f1o1,fou,2,uvb,1)
c
      else if(iflag.eq.5) then
         call daxpy(my*mz,+f1o1,fou,2,uwb,1)
c
      else if(iflag.eq.6) then
         call daxpy(my*mz,+f1o1,fou,2,vwb,1)
c
      else 
         write(*,*) 'statistics_fou2: this flag is wrong!'
      endif
c     
      return
      end subroutine statistics_fou2
c-----|---------------------------------------------------------------|
      subroutine statisticst_fou2(fou,iflag)
      use ctes
      use numbers
      use statistics
      implicit none
c
c     /***************************************************************/
c     /* accumulates statistics in fou space (using only zero-mode). */
c     /* less cost version by sekimoto                               */
c     /*   2nd-moments: ut,vt,wt; iflag=1,2,3                        */
c     /*                                                             */
c     /* input                                                       */
c     /* fou(nsfour)   stresses field in fou-phys-phys space X-cut   */
c     /* output                                                      */
c     /*                                                             */
c     /***************************************************************/
      complex*16 fou(my,mz,ib:ie)
      integer iflag,i,j,k
c-----|---------------------------------------------------------------|
c
      if(myid.ne.0) return
      !write(*,*) 'statisticst_fou2 must be checked'
c     /* 2nd moments: ut,vt,wt */
      if(iflag.eq.1) then
         call daxpy(my*mz,+f1o1,fou,2,utb,1)
c
      else if(iflag.eq.2) then
         call daxpy(my*mz,+f1o1,fou,2,vtb,1)
c
      else if(iflag.eq.3) then
         call daxpy(my*mz,+f1o1,fou,2,wtb,1)
c
      else 
         write(*,*) 'statisticst_fou2: this flag is wrong!'
      endif
c     
      return
      end subroutine statisticst_fou2
c-----|---------------------------------------------------------------|
      subroutine filter_p(p)
      use ctes
      use numbers
      implicit none
c      parallel version by sekimoto 2 Feburary, 2008    
c     /***************************************************************/
c     /* filters the spurious pressure modes of a field represented  */
c     /* by fourier & double chebyshev expansion                     */
c     /*                                                             */
c     /* according to canuto et al. (1987, pp. 233-235) there are    */
c     /* seven spurious pressure modes when using a double-chebyshev */
c     /* collocation approach and primitive variables:               */
c     /*                                                             */
c     /* I    one checkerboard mode: T_my(y)*T_mz(z)                 */
c     /*        -> set u(my,mz)=0            [chebyshev space]       */
c     /*                                                             */
c     /* II   two line modes: T_my(y)  & T_mz(z)                     */
c     /*        -> set u(1,my)=u(my,1)=0     [chebyshev space]       */
c     /*                                                             */
c     /* III  four corner modes: T'_my(y)*(1+/-y)*T'_mz(z)*(1+/-z)   */
c     /*        -> set u(my-1:my,:)=0                                */
c     /*               u(:,mz-1:mz)=0        [chebyshev space]       */
c     /*                                                             */
c     /* IV   the constant mode: T_0(y)*T_0(z)                       */
c     /*      [note: this was not mentioned by canuto et al., but    */
c     /*       e.g. in leriche & labrosse, siam j. sci. comput.,     */
c     /*       22(4), p.1386, 2000]                                  */
c     /*        -> set u(1,1)=0              [chebyshev space]       */
c     /*                                                             */
c     /* note: all of the above are performed for all fourier modes  */
c     /*-------------------------------------------------------------*/
c     /*                                                             */
c     /* input                                                       */
c     /* p(nsfou)             field in fou-phys-phys space           */
c     /*                                                             */
c     /* output                                                      */
c     /* p(nsfou)             filtered field in fou-phys-phys space  */
c     /*                                                             */
c     /***************************************************************/
      complex(nk) p(my,mz,ib:ie)
      integer i
c-----|---------------------------------------------------------------|
      if(ib.eq.0) 
     $  write(*,*)'FILTER: remove the 8 spurious modes...'
c
c     /* pass the field to f-c-c space */
      call transyz(p,'y',-1) !fis->cheb
      call transyz(p,'z',-1) !fis->cheb
c
      do i=max(1,ib),ie
c     
c     /* (0) filter the cst mode in 'y' and 'z' */
c         p(1,1,i)=zero
c     
c     /* (1) filter the line mode in 'y' */
         p(my,1,i)=zero
c     
c     /* (2) filter the line mode in 'z' */
         p(1,mz,i)=zero
c     
c     /* (3) filter the 'checkerboard' mode */
         p(my,mz,i)=zero
c
c     /* (4) corner modes: set last two coeffs to zero in both dirs */
c         do k=1,mz
c            do j=my-1,my
c               p(j,k,i)=zero
c            enddo
c         enddo
c         do k=mz-1,mz
c            do j=1,my
c               p(j,k,i)=zero
c            enddo
c         enddo
c
      enddo
c
c     /* pass the field back to f-p-p space */
      call transyz(p,'z',+1) !cheb->fis
      call transyz(p,'y',+1) !cheb->fis
c
c     /* (4) filter the corner modes: set corner values to zero */
c     /*     (in f-p-p space)                                   */
c      do i=ib,ie
c         p(1 ,1 ,i)=zero
c         p(my,1 ,i)=zero
c         p(1 ,mz,i)=zero
c         p(my,mz,i)=zero
c      enddo
c
      return
      end subroutine filter_p
c-----|---------------------------------------------------------------|
      subroutine fou_filter(u,ifou)
      use ctes
      use numbers
      implicit none 
      complex*16 u(my,mz,ib:ie)
      integer ifou,i,j,k
c-----|---------------------------------------------------------------|
      if (ib.eq.0) write(*,*)'FOU-FILTER above ',ifou
      do i=max(ib,ifou),min(ie,mx1)
         do j=1,my
            do k=1,mz
               u(j,k,i)=zero
            enddo
         enddo
      enddo
c
      return
      end 
c-----|---------------------------------------------------------------|
      subroutine fou_filter2(u,iff,ill)
      use ctes
      use numbers
      implicit none 
      complex*16 u(my,mz,ib:ie)
      integer iff,ill,i,j,k
c-----|---------------------------------------------------------------|
      if (ib.eq.0) write(*,*)'FOU-FILTER2: ',iff,ill
      do i=max(ib,iff),min(ie,ill)
         do j=1,my
            do k=1,mz
               u(j,k,i)=zero
            enddo
         enddo
      enddo
c
      return
      end 
c-----|---------------------------------------------------------------|
      subroutine cheb_filter(u,ify,ifz)
      use ctes
      use numbers
      implicit none 
      complex*16 u(0:my1,0:mz1,ib:ie) !array in c-c-f space
      integer ify,ifz,i,j,k
c-----|---------------------------------------------------------------|
      if (ib.eq.0) write(*,*)'CHEB-FILTER above ',ify,ifz
c     /* filter y-modes */
      do i=ib,ie
         do k=0,mz1
            do j=ify,my1
               u(j,k,i)=zero
            enddo
         enddo
      enddo
c     /* filter z-modes */
      do i=ib,ie
         do k=ifz,mz1            
            do j=0,my1
               u(j,k,i)=zero
            enddo
         enddo
      enddo
c     
      return
      end 
c-----|---------------------------------------------------------------|
#ifdef CMASSFLOW
      subroutine precomp_mass_flow(dtp,alphak,work,nwk)
      use ctes
      use numbers
      use massflow
      implicit none 
c     /***************************************************************/
c     /* compute the solution for the 2d problem                     */
c     /*       (dyy+dzz)u+lambda*u=1  with homogeneous dirichlet bc  */
c     /* and store it for later use in const. mass-flow-rate bc.     */
c     /*                                                             */
c     /* input                                                       */
c     /* dt (dtfixed)      time step                                 */
c     /* alphak(mrk)       runge-kutta coefficients                  */
c     /*                                                             */
c     /* work                                                        */
c     /* work(nwk)         nwk>=2*my*mz                              */
c     /***************************************************************/
      real(nk) work(nwk),alphak(mrk)
      integer nwk,ip1,ip2,ifree,krk
      integer iw1,iw2,iw3,iw4,iw5,iw6,iw7
      real(nk) dtp,flambda
c-----|---------------------------------------------------------------|
#if !defined(DT_FIXED)
      write(*,*)'CMASSFLOW only implemented with DT_FIXED!'
      stop
#endif
      ip1=1
      ip2=my*mz+ip1
      ifree=my*mz+ip2
      if(nwk.lt.ifree)then
         write(*,*)'precomp_mass_flow: work size! ',nwk,ifree
         stop
      endif
c
c     /* loop over all runge-kutta indices */
      do krk=1,mrk
         flambda=-f1o1/(alphak(krk)*dtp*fnu) ! as helco in main1.f90
c     /* set rhs=1 in interior */
         call dcopy(my*mz,f1o1,0,premassfu(1,krk),1) 
c     /* set rhs=0 on boundary */
         call dcopy(my,f0o1,0,premassfu(1,krk),1)
         call dcopy(my,f0o1,0,premassfu(my*(mz-1)+1,krk),1)
         call dcopy(mz,f0o1,0,premassfu(1,krk),my)
         call dcopy(mz,f0o1,0,premassfu(my,krk),my)
c
!#ifdef HELM_SHEN
!         ! wrapped version from LEGCHEB
!         iw1=1
!         iw2=(my+1)*(mz+1)+iw1
!         iw3=(my-2+1)*(mz-2+1)+iw2
!c
!         write(*,*) myid,'dcopy'
!         call dcopy(my*mz,premassfu(1,krk),1,work(iw1),1)
!         write(*,*) myid,'transyz'
!         call transyz_plane(work(iw1),-1) !fis -> cheb     
!         write(*,*) myid,'solveshen'
!         call solveshen(work(iw1),work(iw2),work(iw3),
!     $        'd',-flambda)
!         write(*,*) myid,'dscal2'
!         call dscal2(my*mz,work(iw1),1,-1.d0,work(iw1),1)
!         write(*,*) myid,'cheb2_int'
!         call cheb2_integral(work(iw1),1,aspect,
!     $                       premassf(krk),my,mz)
!         write(*,*) myid,'transyz back'
!         call transyz_plane(work(iw1),+1)
!         call dcopy(my*mz,work(iw1),1,premassfu(1,krk),1)
!#else
#ifdef RECTANGULAR
         call solvediagsv(premassfu(1,krk),work(ip1),work(ip2),
     $        'd',-flambda)
#else
         call solvediag(premassfu(1,krk),work(ip1),work(ip2),
     $        'd',-flambda)
#endif /* RECTANGULAR */
c     /* compute associated mass-flow-rate */
         call dcopy(my*mz,premassfu(1,krk),1,work(ip1),1)
         call transyz_plane(work(ip1),-1) !fis -> cheb
         call cheb2_integral(work(ip1),1,aspect,premassf(krk),my,mz)
!#endif/* HELM_SHEN */
!         write(*,*)'precomp_mass_flow:',krk,premassf(krk),
!     $        alphak(krk),dtp,fnu,flambda
      enddo
c
      return
      end subroutine precomp_mass_flow
c-----|---------------------------------------------------------------|
      subroutine impose_mass_flow(csum,krk,u,work,nwk,iflag)
      use ctes
      use numbers
      use massflow
      implicit none 
c     /***************************************************************/
c     /* impose a specific mass-flow-rate, using a pre-computed      */
c     /* particular solution to the 2d helmholtz problem.            */
c     /*                                                             */
c     /* input                                                       */
c     /* csum           (scalar) mass-flow-rate to be imposed        */
c     /* krk            current runge-kutta index                    */
c     /* u(nfou)        streamwise velocity field, fou-p-p space,X-CUT*/
c     /* iflag          '0' - only compute current mass-flow-rate    */
c     /*                      and return it in "csum"                */
c     /*                '1' - impose mass-flow-rate passed in "csum" */
c     /*                                                             */
c     /* work                                                        */
c     /* work(nwk)         nwk>=my*mz                                */
c     /*                                                             */
c     /* output                                                      */
c     /* u              modified streamwise velocity field.          */
c     /***************************************************************/
      complex(nk) u(my,mz,ib:ie)
      real(nk) work(nwk)
      real(nk) csum,fac,sumg
      integer krk,nwk,iflag
c
c-----|---------------------------------------------------------------|
      if(nwk.lt.my*mz)then
         write(*,*)'impose_mass_flow: work size! ',nwk,my*mz
         stop
      endif
      if(myid.ne.0) return
c
      call dcopy(my*mz,u(1,1,0),2,work,1)

      call transyz_plane(work,-1) !fis -> cheb
      call cheb2_integral(work,1,aspect,sumg,my,mz)
      if(iflag.eq.0)then
         write(*,*)'MASSFLOW: ',sumg
#ifdef RECTANGULAR
        write(*,*)' --> Reynolds =', sumg*f1o4/aspect/fnu
#else /* RECTANGULAR */
        write(*,*)' --> Reynolds =', sumg*f1o4/fnu
#endif /* RECTANGULAR */
         csum=sumg
         return
      endif
      fac=(csum-sumg)/premassf(krk)
      call daxpy(my*mz,fac,premassfu(1,krk),1,u(1,1,0),2)
c
      return
      end subroutine impose_mass_flow
#endif /* CMASSFLOW */
c-----|---------------------------------------------------------------|
      subroutine advection_temp(u,v,w,t,advt,
     $     tmp1,tmp2,tmpfou,tmpzxcut,tmpx,krk,iter)
      use ctes
      use numbers
      implicit none 
c
c     /***************************************************************/
c     /* computes the advection terms of temperature in diver-       */
c     /* gence formulation, for a velocity field which is re-        */
c     /* presented by fourier expansion in x and chebyshev in y,z,   */
c     /* pseudo-spectrally.                                          */
c     /*                                                             */
c     /* input                                                       */
c     /* u,v,w(nsfou)         velocity components in fou-phys-phys   */
c     /* t(nsfou)             temperature in fou-phys-phys           */
c     /* tmp1(nsfis)          BIG work array                         */
c     /* tmp2(nsfis)          BIG work array                         */
c     /* tmpfou(nsfou)        SMALL work array                       */
c     /* tmpx(mgalx+2)        1d-in-x work array                     */
c     /*                                                             */
c     /* output                                                      */
c     /* advt(nsfou)           advection term  in f-p-p space        */
c     /*                                                             */
c     /* local work                                                  */
c     /*                                                             */
c     /* note: 'real*8 u(nsfour)' and 'complex*16 u(nsfou)' have     */
c     /*       equivalent size                                       */
c     /***************************************************************/
      include "mpif.h"
      real(nk) u(*),v(*),w(*),t(*),advt(*)
      real(nk) tmp1(*),tmp2(*),tmpfou(*)
c     /* tmp work array for shuffling (note: complex*16)*/
      complex*16 tmpzxcut(my*mzp*mxp) !z-cut & x-cut pencil work array
      real(nk) tmpx(mgalx+2),dummy
      integer krk,iter
c
c-----|---------------------------------------------------------------|
c     /* preserve u,t (x-cut) */
      call dcopy(nsfour,u,1,tmpfou,1)  
      call dcopy(nsfour,t,1,advt,1)
c            
c     /* shuffle arrays from x-cut to z-cut */
      call chx2z(tmpfou,tmp1,tmpzxcut)
      call chx2z(advt,tmp1,tmpzxcut)
c
c     /* pass u & T to phys-phys-phys space (NOT-in-place) -> tmp1,tmp2 */
      call fourx(tmpfou,tmp1,tmpx,+1)!u: fou->fis
      call fourx(advt,tmp2,tmpx,+1)!t: fou->fis
c
c     /* (1) ========== compute product 'u*T' -> 'tmp1' */
      call mult(nsfis,tmp1,tmp2,tmp1,'n','+')
c
c     /* pass 'u*T' back to f-p-p space -> tmpfou */
      call fourx(tmpfou,tmp1,tmpx,-1)!fis->fou
c
c     /* Re-shuffle arrays from z-cut to x-cut */
      call chz2x(tmpfou,tmp1,tmpzxcut)
c
c     /* call statistics module */
      if((krk.eq.1).and.(mod(iter-1,nhist).eq.0))then
         call statisticst_fou2(tmpfou,1) !u*T
      endif
c
c     /* derivative in 'x': d(u*T)/dx -> advt */
      call derivx(tmpfou,advt,1,f1o1,'n')
c
c     /* (2) =========== compute product 'v*T' -> 'tmp1' */
c     /* preserve v (x-cut) */
      call dcopy(nsfour,v,1,tmpfou,1)  
c            
c     /* shuffle arrays from x-cut to z-cut */
      call chx2z(tmpfou,tmp1,tmpzxcut)
c
      call fourx(tmpfou,tmp1,tmpx,+1)!fou->fis
      call mult(nsfis,tmp1,tmp2,tmp1,'n','+')
      call fourx(tmpfou,tmp1,tmpx,-1)!fis->fou
c
c     /* Re-shuffle arrays from z-cut to x-cut */
      call chz2x(tmpfou,tmp1,tmpzxcut)
c
c     /* call statistics module */
      if((krk.eq.1).and.(mod(iter-1,nhist).eq.0))then
         call statisticst_fou2(tmpfou,2) !v*T
      endif
c
c     /* pass 'v*T' to f-c-p space -> tmpfou */
      call transyz(tmpfou,'y',-1) !fis->cheb
c
c     /* derivative in 'y': d(v*T)/dy -> tmp1 */
      call derivyz(tmpfou,tmp1,dummy,'y',1,f1o1,'n')
c
c     /* pass d(v*T)/dy to f-p-p space -> tmp1 */
      call transyz(tmp1,'y',+1) !cheb->fis
c
c     /* add d(v*T)/dy to 'advt' */
      call daxpy(nsfour,+f1o1,tmp1,1,advt,1)
c
c     /* (3) =========== compute product 'W*T' -> 'tmp1' */
c     /* preserve w (x-cut) */
      call dcopy(nsfour,w,1,tmpfou,1)  
c            
c     /* shuffle arrays from x-cut to z-cut */
      call chx2z(tmpfou,tmp1,tmpzxcut)
c
      call fourx(tmpfou,tmp1,tmpx,+1)!fou->fis
      call mult(nsfis,tmp1,tmp2,tmp1,'n','+')
      call fourx(tmpfou,tmp1,tmpx,-1)!fis->fou
c
c     /* Re-shuffle arrays from z-cut to x-cut */
      call chz2x(tmpfou,tmp1,tmpzxcut)
c
c     /* call statistics module */
      if((krk.eq.1).and.(mod(iter-1,nhist).eq.0))then
         call statisticst_fou2(tmpfou,3) !w*T
      endif
c
c     /* pass 'w*T' to f-p-c space -> tmpfou */
      call transyz(tmpfou,'z',-1) !fis->cheb
c
c     /* derivative in 'z': d(w*T)/dz -> tmp1 */
      call derivyz(tmpfou,tmp1,dummy,'z',1,f1o1,'n')
c
c     /* pass d(w*T)/dz to f-p-p space -> tmp1 */
      call transyz(tmp1,'z',+1) !cheb->fis
c
c     /* add d(w*T)/dz to 'advt' */
      call daxpy(nsfour,+f1o1,tmp1,1,advt,1)
c
      return
      end subroutine advection_temp
c-----|---------------------------------------------------------------|
      real*8 function get_dpdx(aspect,fnu)
      use ctes,only:nk
      use numbers
      implicit none 
c
c     /* this function determines the driving pressure gradient      */
c     /* for rectangular duct flow, using tabulated values for the   */
c     /* laminar flow solution.                                      */
c     /* particularly: -(dyy+dzz)U=dpval -->  max(U)=1.0 for each    */
c     /* aspectval.                                                  */
c     /* this means that dp/dx=-\nu*dpval for pure pressure-         */
c     /* driven, laminar flow. (assuming -1<=y<=+1)                  */
c
c     /* local variables */
      integer nval,index,i
      real(nk) aspect,fnu
      parameter(nval=19)
      real(nk) aspectval(nval),dpval(nval)
      data aspectval /1.d0,1.5d0,2.d0,2.5d0,3.d0,3.5d0,4.d0,4.5d0,5.d0,
     $     5.5d0,6.d0,6.5d0,7.d0,7.5d0,8.d0,8.5d0,9.d0,9.5d0,1.d1/
      data dpval /3.39344926982230,2.48086300035167, !1,1.5
     $     2.19545075660646,2.08474987413869, !2.0,2.5
     $     2.03778220024811,2.01705223568555, !3.0,3,5
     $     2.00773897291851,2.00352109020906, !4.0,4.5
     $     2.00160386398057,2.00073094391977, !5.0,5.5
     $     2.00033319895779,2.00015190434113, !6.0,6.5
     $     2.00006925611912,2.00003157591041, !7.0,7.5
     $     2.00001439653781,2.00000656390479, !8.0,8.5
     $     2.00000299272912,2.00000136449820, !9.0,9.5
     $     2.00000062212652/    !10.0
c-----|---------------------------------------------------------------|
      index=0
      do i=1,nval
         if(abs(aspectval(i)-aspect).le.1.d-5)then
            index=i
         endif
      enddo
      if(index.eq.0)then
         write(*,*)'this aspect ratio is not tabulated! ',aspect
         stop
      endif
      get_dpdx=-fnu*dpval(index)
c
      return
      end function get_dpdx
c-----|---------------------------------------------------------------|
      subroutine max_3din(u,divmax)
      use ctes
      use numbers
      implicit none
      include 'mpif.h'
c
c     /* max of 3d f-p-p-space field w/o boundaries, X-CUT */
c     /* 'master' holds the result                         */
c
      complex*16 u(my,mz,ib:ie)
      real(nk) fmax,divmax
      integer i,j,k,ierr
c
c-----|---------------------------------------------------------------|
      fmax=-dbig
      do i=ib,ie
         do k=2,mz-1
            do j=2,my-1
               fmax=max(fmax,abs(u(j,k,i)))
            enddo
         enddo
      enddo
      call MPI_REDUCE(fmax,divmax,1,MPI_DOUBLE_PRECISION,
     &     MPI_MAX,master,MPI_COMM_WORLD,ierr)     
c
      return
      end subroutine max_3din
c-----|---------------------------------------------------------------|
      subroutine slip_velocity_3d(u,v,w)
      use ctes
      use running, only: slipmx
      use numbers
      implicit none 

      include "mpif.h"
      real(nk) u(mgalx+2,my,kb:ke),v(mgalx+2,my,kb:ke),
     &         w(mgalx+2,my,kb:ke)
      real(nk) fmax
      integer i,j,k,ierr
c
c-----|---------------------------------------------------------------|
      fmax=-dbig
c
c     /* y=-1 & y=+1 planes */
      do k=kb,ke
         do i=1,mgalx
            fmax=max(fmax,u(i,1,k))
            fmax=max(fmax,u(i,my,k))

            fmax=max(fmax,v(i,1,k))
            fmax=max(fmax,v(i,my,k))

            fmax=max(fmax,w(i,1,k))
            fmax=max(fmax,w(i,my,k))
         enddo
      enddo
c
c     /* z=-a & z=+a planes */
      if(myid.eq.0) then
         k=1       
         do j=1,my
            do i=1,mgalx
               fmax=max(fmax,u(i,j,k))
               fmax=max(fmax,v(i,j,k))
               fmax=max(fmax,w(i,j,k))
            enddo
         enddo
      endif
      if(myid.eq.lastid) then
         k=mz         
         do j=1,my
            do i=1,mgalx
               fmax=max(fmax,u(i,j,k))
               fmax=max(fmax,v(i,j,k))
               fmax=max(fmax,w(i,j,k))
            enddo
         enddo
      endif
      call MPI_REDUCE(fmax,slipmx,1,MPI_DOUBLE_PRECISION,
     &     MPI_MAX,master,MPI_COMM_WORLD,ierr)     
c
      if(myid.eq.master) write(*,*)'SLIP-VELOCITY: ',slipmx
c
      return
      end subroutine slip_velocity_3d
c-----|---------------------------------------------------------------|
#ifdef SET_NOSLIP
c     /*   by sekimoto 12/Sep/2008 */
      subroutine set_noslip_3dfis(u,v,w)
      use ctes
      use numbers
      implicit none 
      real(nk) u(mgalx+2,my,kb:ke),v(mgalx+2,my,kb:ke)
     &         w(mgalx+2,my,kb:ke)
      integer i,j,k
c-----|---------------------------------------------------------------|
c
c     /* y=-1 & y=+1 planes */
      do k=kb,ke
         do i=1,mgalx        
            u(i,1,k)=f0o1
            u(i,my,k)=f0o1
            v(i,1,k)=f0o1
            v(i,my,k)=f0o1
            w(i,1,k)=f0o1
            w(i,my,k)=f0o1
         enddo
      enddo
c
c     /* z=-a & z=+a planes */
      if(kb.eq.1) then
         k=1
         do i=1,mgalx
            do j=1,my
               u(i,j,k)=f0o1
               v(i,j,k)=f0o1
               w(i,j,k)=f0o1
            enddo
         enddo
      endif
      if(ke.eq.mz) then
         k=mz
         do i=1,mgalx
            do j=1,my
               u(i,j,k)=f0o1
               v(i,j,k)=f0o1
               w(i,j,k)=f0o1
            enddo
         enddo
      endif
c
      return
      end subroutine set_noslip_3dfis
c-----|---------------------------------------------------------------|
c     /*   by sekimoto 12/Sep/2008 */
      subroutine set_noslip_3dfou(u,v,w)
      use ctes
      use numbers
      implicit none 
      complex(nk) u(my,mz,ib:ie),v(my,mz,ib:ie),w(my,mz,ib:ie)
      integer i,j,k
c-----|---------------------------------------------------------------|
c
c     /* y=-1 & y=+1 planes */
      do i=ib,ie
         do k=1,mz
            do j=1,my,my-1
               u(j,k,i)=zero
               v(j,k,i)=zero
               w(j,k,i)=zero
            enddo
         enddo
      enddo
c
c     /* z=-a & z=+a planes */
      do i=ib,ie         
         do k=1,mz,mz-1
            do j=1,my
               u(j,k,i)=zero
               v(j,k,i)=zero
               w(j,k,i)=zero
            enddo
         enddo
      enddo
c
      return
      end subroutine set_noslip_3dfou
#endif /* SET_NOSLIP */
c-----|---------------------------------------------------------------|
      subroutine L2_norm_fou2d_max(u,anorm_max)
      use ctes
      use numbers
      implicit none 
c
c     /* returns maximum of 2d L2 norm of u */
c
      include "mpif.h"
      complex*16 u(my,mz,ib:ie)
      real(nk) anorm,anorm_max,dzz,dyy,fmax
      integer i,j,k,ierr
c-----|---------------------------------------------------------------|
      anorm_max=-dbig
      do i=ib,ie
         anorm=f0o1
         do k=1,mz
            dzz=dz(k)
            do j=2,my-1
               dyy=dy(j)
               anorm=anorm+u(j,k,i)*dconjg(u(j,k,i))*dyy*dzz
               fmax=max(anorm_max,anorm)
            enddo
         enddo
      enddo
      call MPI_REDUCE(fmax,anorm_max,1,MPI_DOUBLE_PRECISION,
     &     MPI_MAX,master,MPI_COMM_WORLD,ierr)
c
      return
      end subroutine L2_norm_fou2d_max

c-----|---------------------------------------------------------------|
#ifdef TEMPERATURE
      subroutine set_bc_temp
      use ctes,only:nk,my,mz,mx1,ib,ie,master,myid,  
     &     boundy1,boundy2,boundz1,boundz2, bcy,bcz
      use numbers
      implicit none 

      !real(nk) bcy(mz,2,ib:ie),bcz(my,2,ib:ie)
      integer i,j,k
      !character*1 boundy1,boundy2,boundz1,boundz2

      allocate(bcy(mz,2,0:mx1),bcz(my,2,0:mx1))
      boundy1='d';boundy2='d';
      boundz1='n';boundz2='n';
c     
      if(myid.eq.master)then
         i=0
         do k=1,mz
            bcy(k,1,i)=+f1o1
            bcy(k,2,i)=-f1o1
         enddo
         do j=1,my
            bcz(j,1,i)=+f0o1
            bcz(j,2,i)=+f0o1
         enddo
      endif
        
      do i=max(ib,1),ie
         do k=1,mz
            bcy(k,1,i)=f0o1
            bcy(k,2,i)=f0o1
         enddo
         do j=1,my
            bcz(j,1,i)=f0o1
            bcz(j,2,i)=f0o1
         enddo
      enddo
      return
      end subroutine
c---------------------------------------------------------------------|
      subroutine t_corner(t)
      use ctes
      use numbers
      implicit none 

      complex*16 t(my,mz,ib:ie)
      integer i

      if(myid.eq.master)then
         i=0
         t(1,1,i)=dcmplx(f1o1,f0o1)
         t(my,1,i)=dcmplx(-f1o1,f0o1)
         t(1,mz,i)=dcmplx(f1o1,f0o1)
         t(my,mz,i)=dcmplx(-f1o1,f0o1)
      endif

      do i=max(ib,1),ie
         t(1,1,i)=zero
         t(my,1,i)=zero
         t(1,mz,i)=zero
         t(my,mz,i)=zero
      enddo

      return
      end subroutine
#endif /* TEMPERATURE */
c----|----------------------------------------------------------------|
#ifdef SHOOTING
      subroutine make_edge_state(u,v,w,beta)
      use ctes
      use numbers
      implicit none 
      parameter(iout=55)      

      character*80 fileedge
      complex*16 u(my,mz,ib:ie),v(my,mz,ib:ie),w(my,mz,ib:ie)
      real(nk) uzero(my,mz),vzero(my,mz),wzero(my,mz),pzero(my,mz)
      parameter(nparamx=15)      
      real(nk) param(nparamx),beta
c      
      if (ib.eq.0) then
c     /* read zeromode near the edge state          */
         write(*,*) 'reading zeromode of edge-state from(2D-binary):'
         read(*,*) fileedge
         nparam=7
         open(iout,file=fileedge,form="unformatted")
         read(iout)(param(ii),ii=1,nparam)
         mgalxr=int(param(1))
         myr=int(param(2))
         mzr=int(param(3))
         if(mgalxr.ne.mgalx.or.myr.ne.my.or.mzr.ne.mz)then
            write(*,*)myid,'mismatching dimensions!',
     $           mgalxr,myr,mzr
            stop
         endif
         read(iout)((uzero(j,k),j=1,my),k=1,mz)
         read(iout)((vzero(j,k),j=1,my),k=1,mz)
         read(iout)((wzero(j,k),j=1,my),k=1,mz)
         read(iout)((pzero(j,k),j=1,my),k=1,mz)
         close(iout)
c     /* compute edge state (zero-mode):            */
c     /*              u_0 = <u> + \beta (u_0 - <u>) */
         i=0
         do k=1,my
            do j=1,mz
               u(j,k,i)=dcmplx(uzero(j,k),0.d0) +
     $                  beta*(u(j,k,i)-dcmplx(uzero(j,k),0.d0))
               v(j,k,i)=dcmplx(vzero(j,k),0.d0) +
     $                  beta*(v(j,k,i)-dcmplx(vzero(j,k),0.d0))
               w(j,k,i)=dcmplx(wzero(j,k),0.d0) +
     $                  beta*(w(j,k,i)-dcmplx(wzero(j,k),0.d0))
            enddo
         enddo
      endif
      
c     /* compute edge state (non-zero-mode):            */
c     /*              u = \beta (u) */
      do i=max(ib,1),ie
         do k=1,my
            do j=1,mz
               u(j,k,i)=beta*(u(j,k,i))
               v(j,k,i)=beta*(v(j,k,i))
               w(j,k,i)=beta*(w(j,k,i))
            enddo
         enddo
      enddo
      
      return
      end subroutine
#endif /* SHOOTING */
c----|----------------------------------------------------------------|
#ifdef CUTOFF
!    cut off high fourier modes to acheive reasonable spectram, u,v,w 
!    need MAKE_PRESSURE option to reconstruct initial fields 
!    remained fourier modes.
      subroutine cutoff_fou_modes(u,v,w,icut)
      use ctes
      use numbers
      implicit none 

      complex*16 u(my,mz,ib:ie),v(my,mz,ib:ie),w(my,mz,ib:ie)
      integer i,j,k

      if(ie.gt.(mx1-icut)) then
         do i=max((mx1-icut),ib),ie
            do k=1,my
               do j=1,mz
                  u(j,k,i)=zero
                  v(j,k,i)=zero
                  w(j,k,i)=zero
               enddo
            enddo
         enddo
      endif
      return
      end subroutine
#endif /*CUTOFF*/
c----|----------------------------------------------------------------|
#ifdef FORCING_GLOBAL 
      subroutine force_global(h1,h2,h3,rfac,Am)
c    /*  add global force to h1,h2,h3 to generate turbulence */
c    /*  note: using first 2 fourier modes  */
c    /*    (zeromode and first fourier mode)*/
c    /*  Ref: Markus report on puff generation */
      use ctes
      use numbers
      implicit none 
      complex(nk) h1(my,mz,ib:ie),h2(my,mz,ib:ie),h3(my,mz,ib:ie)
      real(nk) rfac,Am,ran,r2,theta,Amthr2,xfrcng,yfrcng,zfrcng
      integer i,j,k

      ran=0.5   
!     rfac=akkc(krk)    
      do i=max(ib,0),1
         do k=1,mz
            do j=1,my
               r2=y(j)*y(j)+z(k)*z(k);
               theta=atan2(y(j),z(k));
               Amthr2=Am*(1.d0+dcos(theta))*r2
c               tf(j,k)=r2*atan2(y(j),z(k));
c               tf2(j,k)=r2*atan2(y(j),z(k))*(1+cos(theta));
               if (i.eq.0) then
                  xfrcng=Amthr2
     $                 *(1.d0-ran)*dsin(3.d0*theta)
                  h1(j,k,i)=h1(j,k,i)+dcmplx(rfac*xfrcng)
c
                  yfrcng=Amthr2
     $                 *ran*(dsin(theta)*dsin(theta));
                  h2(j,k,i)=h2(j,k,i)+dcmplx(rfac*yfrcng)
c
                  zfrcng=Amthr2
     $                 *ran*dsin(theta)*cos(theta)
                  h3(j,k,i)=h3(j,k,i)+dcmplx(rfac*zfrcng) 
c
               elseif(i.eq.1) then
                  xfrcng=Amthr2*(-0.5d0)
     $                 *(1.d0-ran)*dsin(3.d0*theta);
                  h1(j,k,i)=h1(j,k,i)+dcmplx(rfac*xfrcng)
c
                  yfrcng=Amthr2*(-0.5d0)
     $                 *ran*(dsin(theta)*dsin(theta))
                  h2(j,k,i)=h2(j,k,i)+dcmplx(rfac*yfrcng)
c
                  zfrcng=Amthr2*(-0.5d0)
     $                 *ran*dsin(theta)*cos(theta);
                  h3(j,k,i)=h3(j,k,i)+dcmplx(rfac*zfrcng)
               endif
            enddo
         enddo
      enddo      
      end subroutine force_global
#endif /* FORCING_GLOBAL */ 
c-----|---------------------------------------------------------------|
#ifdef PUFF_FORCING
      subroutine add_puff_forcing(h1,h2,h3,scal,
     $     tmpfis,tmpx,tmpfou,tmpzxcut)
      use ctes
      use numbers
      implicit none 
c     
c     /***************************************************************/
c     /* adds a localized forcing term to the rhs vector. the forcing*/
c     /* is computed in physical space, then transformed.            */
c     /*                                                             */
c     /* input                                                       */
c     /* h1,h2,h3(nsfou) r.h.s.vector in fourier-p-p space, x-cut    */
c     /* real*8 scal(1)  scalar factor (i.e. '+2*alpha_k')           */
c     /* tmpfis(nsfis)   BIG work array                              */
c     /* tmpx(mgalx+2)   1d-x work array                             */
c     /* tmpfou(nsfou)   fourier-space work array                    */
c     /* tmpzxcut(my*mzp*mxp) SMALLER complex*16 work array          */
c     /*                                                             */
c     /* output                                                      */
c     /* h1,h2,h3(nsfou) updated r.h.s.vector in f-p-p space, x-cut  */
c     /***************************************************************/

      include "mpif.h"
      real(nk) tmpfis(mgalx+2,my,kb:ke),tmpx(mgalx+2)
      complex*16 tmpfou(nsfou),h1(nsfou),h2(nsfou),h3(nsfou)
c     /* tmp work array for shuffling (note: complex*16)*/
      complex*16 tmpzxcut(my*mzp*mxp) !z-cut & x-cut pencil work array
      real(nk) scal, fampl, flambda,facx,xx,zz,yy,rr2,theta,fac
      integer myh,mzh,i,j,k
c
c-----|---------------------------------------------------------------|
      fampl=0.06d0
      flambda=0.55d0
      myh=my/2+1
      mzh=mz/2+1
      if(mod(my,2).ne.1.or.mod(mz,2).ne.1)then
         write(*,*)'add_puff_forcing: my and mz need to be odd!'
         stop
      endif
      if(myid.eq.0)write(*,*)'add_puff_forcing: ',fampl,flambda
c
c     /* x-component: start in physical space, z-cut */
      call dcopy(nsfis,f0o1,0,tmpfis,1) !initialize to zero
      do k=kb,ke
         zz=z(k)
         do j=1,my
            yy=y(j)
            rr2=yy**2+zz**2
            theta=datan2(yy,zz)
            do i=1,mgalx
               xx=x(i)
               if(xx.le.f2o1*pi)then
                  facx=(f1o1-dcos(xx))
               endif
               fac=fampl*(f1o1+dcos(theta))*facx*rr2
               tmpfis(i,j,k)=fac*dsin(f3o1*theta)*(f1o1-flambda)
            enddo
         enddo
      enddo
      if(kb.le.mzh.and.ke.ge.mzh)then
c     /* correct the origin (y=0,z=0) sinze "datan2" might return "NaN" */
         tmpfis(:,myh,mzh)=f0o1
      endif

c     /* pass 'tmpfis' to f-p-p space, z-cut -> 'tmpfou' */
      call fourx(tmpfou,tmpfis,tmpx,-1)!fis->fou
c
c     /* re-shuffle 'tmpfou' from z-cut to x-cut */
      call chz2x(tmpfou,tmpfis,tmpzxcut)
c
c     /* add 'scal*tmpfou' to r.h.s. vector 'h1' */
      call daxpy(nsfour,scal,tmpfou,1,h1,1)
c     
c     /* y-component: start in physical space, z-cut */
      call dcopy(nsfis,f0o1,0,tmpfis,1) !initialize to zero

      do k=kb,ke
         zz=z(k)
         do j=1,my
            yy=y(j)
            rr2=yy**2+zz**2
            theta=datan2(yy,zz)
            do i=1,mgalx
               xx=x(i)
               if(xx.le.f2o1*pi)then
                  facx=(f1o1-dcos(xx))
               endif
               fac=fampl*(f1o1+dcos(theta))*facx*rr2
               tmpfis(i,j,k)=fac*dsin(theta)**2*flambda
            enddo
         enddo
      enddo
      if(kb.le.mzh.and.ke.ge.mzh)then
c     /* correct the origin (y=0,z=0) since "datan2" might return "NaN" */
         tmpfis(:,myh,mzh)=f0o1
      endif

c     /* pass 'tmpfis' to f-p-p space, z-cut -> 'tmpfou' */
      call fourx(tmpfou,tmpfis,tmpx,-1)!fis->fou
c
c     /* re-shuffle 'tmpfou' from z-cut to x-cut */
      call chz2x(tmpfou,tmpfis,tmpzxcut)
c
c     /* add 'scal*tmpfou' to r.h.s. vector 'h2' */
      call daxpy(nsfour,scal,tmpfou,1,h2,1)
c     
c     /* z-component: start in physical space, z-cut */
      call dcopy(nsfis,f0o1,0,tmpfis,1) !initialize to zero
      do k=kb,ke
         zz=z(k)
         do j=1,my
            yy=y(j)
            rr2=yy**2+zz**2
            theta=datan2(yy,zz)

            do i=1,mgalx
               xx=x(i)
               if(xx.le.f2o1*pi)then
                  facx=(f1o1-dcos(xx))
               endif
               fac=fampl*(f1o1+dcos(theta))*facx*rr2
               tmpfis(j,k,i)=fac*dsin(theta)*dcos(theta)*flambda
            enddo
         enddo
      enddo
      if(kb.le.mzh.and.ke.ge.mzh)then
c     /* correct the origin (y=0,z=0) sinze "datan2" might return "NaN" */
         tmpfis(:,myh,mzh)=f0o1
      endif

c     /* pass 'tmpfis' to f-p-p space, z-cut -> 'tmpfou' */
      call fourx(tmpfou,tmpfis,tmpx,-1)!fis->fou
c
c     /* re-shuffle 'tmpfou' from z-cut to x-cut */
      call chz2x(tmpfou,tmpfis,tmpzxcut)
c
c     /* add 'scal*tmpfou' to r.h.s. vector 'h3' */
      call daxpy(nsfour,scal,tmpfou,1,h3,1)
c     
      return
      end subroutine add_puff_forcing
#endif /* PUFF_FORCING */ 
c----|---------------------------------------------------------------
      subroutine phase_shift_fou(ff,shift)
      use ctes
      use numbers
      implicit none 
c     /************************************************************/
c     /* phase shift in fourier 3d space                          */
c     /*                                                          */
c     /* ff: fourier coefficients                                 */
c     /* shift: shift for streamwise(x) direction                 */
c     /************************************************************/
      complex(nk) ff(my,mz,ib:ie)
      real(nk) shift
      integer ii,k,j
c
      do ii=ib,ie
         do k=1,mz
            do j=1,my
               ff(j,k,ii)=
     $              ff(j,k,ii)*exp(cii*dcmplx(xalp(ii)*(-shift)))
            enddo
         enddo
      enddo
c
      return
      end subroutine phase_shift_fou
c----|---------------------------------------------------------------
      subroutine make_wfou(u,v,wfou,tmp2,work2d,sol,work2db,solb)
      use ctes
      use numbers
      use chebystuff
      implicit none
c    /**********************************************************/
c    /*  reconstruct spanwise velocity using u,v (f-p-p)       */
c    /*  by following: (probably no singularity for zeromode)  */ 
c    /*  (dw/dz)= -cii*xalp(i)*u(j,k,i)-(dv/dy)                */
c    /*  for each Fourier modes solve:                         */     
c    /*   w(1:my,1:mz,i)*Dz(1:mz,1:mz)= -cii*xalp(i)*u - Dy*v  */ 
c    /*   except for the boundaries.                           */
c    /*                                                        */
c    /* this does not work, because the inverse of D is        */
c    /* really less accuracy                                   */
c    /* --> inverse of D is computed by SVD                    */
c    /*     and the boundary is set to zero                    */
c    /*                                                        */
c    /**********************************************************/
      complex(nk) u(my,mz,ib:ie),v(my,mz,ib:ie)
      complex(nk) wfou(my,mz,ib:ie),tmp2(my,mz,ib:ie)

      real(nk) work2d(my,mz),sol(my,mz)
      real(nk) work2db(my-2,mz-2),solb(my-2,mz-2)

      real(nk) dummy
      
      integer mm,lwk
      integer i,j,kk,k,jj,info,imode
c
      real(nk),dimension(:),allocatable:: wk
c
      mm=my-2
      if(isvdflag.eq.0)then
         if(ib.eq.0) then 
            write(*,*) 'compute singular value and vectors: isvdflag',
     $        isvdflag
         endif
         lwk=5*(mm) + 10*mm
         allocate(wk(lwk)) ! for SVD
c     /* compute transposed matrix of chebyshev colloc. deriv. mat.*/
         kk=0
         do k=2,mz-1
            kk=kk+1
            jj=0
            do j=2,my-1
               jj=jj+1
               DTB(jj,kk)=DC(k,j)
            enddo
         enddo
c     /* compute singular value and vectors */
         call DGESVD('a','a',mm,mm,DTB,mm,DS,DU,mm,DVT,mm,wk,lwk,info)
         isvdflag=1
c        !DTB is destroyed         
      endif
      call dcopy(nsfour,f0o1,0,wfou,1)
c
c     /* compute -(dv/dy) */
c     /* pass 'v' to f-c-p space -> v */
      call transyz(v,'y',-1)    !fis->cheb
c     
c     /* derivative in 'y': d(v)/dy -> tmp2 */
      call derivyz(v,tmp2,dummy,'y',1,-f1o1,'n')
c     
c     /* pass tmp back to f-p-p space -> tmp */
      call transyz(tmp2,'y',+1) !cheb->fis
c     
c     /* pass 'v' back to f-p-p space -> v */
      call transyz(v,'y',+1)    !cheb->fis
c
c     /* compute -(du/dx) */
      call derivx(u,tmp2,1,-f1o1,'+')
c here tmp2 is r.h.s of W*DCz' = -(dv/dy)-(du/dx)
c---------
c     /* for imode=0 */
      if(ib.eq.0) then
c     /* w(1:my,1:mz,i)*Dz(1:mz,1:mz)'=- Dy*v */
         imode=0
         call com2re(my,mz,tmp2(1,1,imode),my,work2d,my,'r')
c     /* use only inner values */
         kk=0
         do k=2,mz-1
            kk=kk+1
            jj=0
            do j=2,my-1
               jj=jj+1
               work2db(jj,kk)=work2d(j,k)
            enddo
         enddo
c     
c     /* w(:,:,0)= - Dy*v*inv(DBz) */
         call dgemm('n','n',mm,mm,mm,1.d0,work2db,mm,DCM1,
     $        mm,0.d0,solb,mm)      
c     /* boundary condition is satisfied by using DT(2:my-1,2:mz-1)*/
c     /* so it need to repair boundary  */
         kk=0
         do k=2,mz-1
            kk=kk+1
            jj=0
            do j=2,my-1
               jj=jj+1
               sol(j,k)=solb(jj,kk)
            enddo
         enddo
c     
         call re2com(my,mz,wfou(1,1,imode),my,sol,my,'r')
      endif !/*master*/
c
      do imode=max(1,ib),ie         
c     /* (1) real part */
         call com2re(my,mz,tmp2(1,1,imode),my,work2d,my,'r')
c     /* use only inner values */
         kk=0
         do k=2,mz-1
            kk=kk+1
            jj=0
            do j=2,my-1
               jj=jj+1
               work2db(jj,kk)=work2d(j,k)
            enddo
         enddo
c     /* w(:,:,0)= rhs *inv(DBz) */
         call dgemm('n','n',mm,mm,mm,1.d0,work2db,mm,DCM1,
     $        mm,0.d0,solb,mm)    
         kk=0
         do k=2,mz-1
            kk=kk+1
            jj=0
            do j=2,my-1
               jj=jj+1
               sol(j,k)=solb(jj,kk)
            enddo
         enddo
         call re2com(my,mz,wfou(1,1,imode),my,sol,my,'r')
c
c     /* (2) imaginary part */
         call com2re(my,mz,tmp2(1,1,imode),my,work2d,my,'i')
c    
c     /* use only inner values */
         kk=0
         do k=2,mz-1
            kk=kk+1
            jj=0
            do j=2,my-1
               jj=jj+1
               work2db(jj,kk)=work2d(j,k)
            enddo
         enddo
         call dgemm('n','n',mm,mm,mm,1.d0,work2db,mm,DCM1,
     $        mm,0.d0,solb,mm)    
         kk=0
         do k=2,mz-1
            kk=kk+1
            jj=0
            do j=2,my-1
               jj=jj+1
               sol(j,k)=solb(jj,kk)
            enddo
         enddo
         call re2com(my,mz,wfou(1,1,imode),my,sol,my,'i')
      enddo
c
c     /* set boundary to zero */
c     /* y=-1 & y=+1 planes   */
      do i=ib,ie         
         do k=1,mz
            do j=1,my,my-1
               wfou(j,k,i)=zero
            enddo
         enddo
      enddo
c
c     /* z=-a & z=+a planes */
      do i=ib,ie
         do k=1,mz,mz-1
            do j=1,my
               wfou(j,k,i)=zero
            enddo
         enddo
      enddo
c
      call filter_cheby(wfou,'z')
      return
      end subroutine make_wfou
c----|---------------------------------------------------------------
      subroutine filter_cheby(ff,cdir)
      use ctes
      use numbers
      implicit none 
c     /************************************************************/
c     /* apply chebyshev filter to remove (0, mz) mode            */
c     /* chebyshev transform of a 3d field in both 'y' and 'z'    */
c     /*          */
c     /* fpp <-> fcc                                              */
c     /* using chebyshev.F                                        */
c     /* call transyz(u,'z',-1) !fis->cheb                        */
c     /* call transyz(v,'z',+1) !cheb->fis                        */ 
c     /************************************************************/
      complex*16 ff(1:my,1:mz,ib:ie)
      character*1 cdir
      complex*16 const
      integer ix,j,k
      if(cdir.eq.'y') then
         write(*,*) 'no need to filter in y direction'
      elseif(cdir.eq.'z')then
c         write(*,*) 'apply cheby filter of mz mode and constant mode'
c     /* remove the last chebyshev mode */
         call transyz(ff,'z',-1)
         do ix=ib,ie
            do j=1,my
               ff(j,mz,ix)=zero
            enddo
         enddo
         call transyz(ff,'z',+1)
c     /* remove the constant mode */
         do ix=ib,ie
            do j=1,my
               const=ff(j,1,ix);
               do k=1,mz
                  ff(j,k,ix)=ff(j,k,ix)-const
               enddo
            enddo
         enddo
      endif
      return
      end subroutine filter_cheby
c----|---------------------------------------------------------------
