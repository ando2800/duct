      subroutine fourx(fou,phys,aux,iopt)
      use ctes
      use numbers
      use rfttmp
      implicit none
c     /***************************************************************/
c     /*  including fftw (inplace)version by sekimoto 11/Sep/2008    */
c     /*                                                             */
c     /*  performs real-to-complex fourier transform, Z-cut          */
c     /*     from/ to F   -Phys-Phys (array fou)                     */
c     /*     to/ from Phys-Phys-Phys (array phys)                    */
c     /*                                                             */
c     /*  iopt >=0    ===>  inverse transforms( fou ---> fis)        */
c     /*  iopt < 0    ===>  direct  transforms( fis ---> fou)        */
c     /*                                                             */
c     /*  performs full de-aliasing according to the 2/3 rule        */
c     /*    fis ---> fou : supposes high frecuency modes are not     */
c     /*                   needed (dealiasing) and throws them away  */
c     /*                                                             */
c     /* input/output                                                */
c     /* real*8   fis(my,kb:ke,mgalx+2)  physical space array        */
c     /* cmplx*16 fou(my,kb:ke,0:mx1)    fourier space array         */
c     /* cloc                         'l' - in-place transform       */
c     /*                               (fou & fis have same memory   */
c     /*                                  address)                   */
c     /*                              [else] - out-of-place          */
c     /*                                                             */
c     /* work                                                        */
c     /* real*8   aux(mgalx+2)        auxiliary buffer, contiguous   */
c     /*                                                             */
c     /* note: does NOT allow for "in-place" transform!              */
c     /***************************************************************/
      include"fftw3.f"
!      real(nk) fou(*),phys(my,kb:ke,*),aux(*)
      real(nk) fou(*),phys(mgalx+2,my,kb:ke),aux(*)
      integer iopt,k,j
c-----|---------------------------------------------------------------|
c
c     /* physical to fourier transforms  */
      if (iopt.lt.0) then
         do k=kb,ke
            do j=1,my
c     /* copy to contiguous -> aux */
               call dcopy(mgalx+2,phys(1,j,k),1,aux(1),1)
               call dfftw_execute_dft_r2c(plan_fft_forward,aux,aux)
c     /* de-alias into 'fou' */
               call pack_fftw(aux,fou,j,k) ! pack and scaled by mgalx
            enddo
         enddo
c     
c     /* fourier to physical transforms  */
      else
         do k=kb,ke
            do j=1,my
c     /* zero-fill & copy to contiguous -> aux */
               call fill(fou,aux,j,k)
               call dfftw_execute_dft_c2r(plan_fft_backward,aux,aux)
c     /* copy back into 'phys' */
               call dcopy(mgalx+2,aux,1,phys(1,j,k),1)
            enddo
         enddo
      endif
c     
      return
      end subroutine fourx
c-----|---------------------------------------------------------------|
      subroutine pack_fftw(foud,fou,j,k)
      use ctes
      use numbers
      implicit none 
c        by sekimoto 11/Sep/2008 
c          --> updated for transposed array (31/3/2020)
c     /* reduces a complex array to 2/3 size                         */
c     /* pack routine for fftw to dealias                            */
c     /* input                                                       */
c     /* foud(0:mgx)                                                 */
c     /*                                                             */
c     /* output                                                      */
c     /* fou(0:mx1,my,kb:ke)                                         */
c
      complex(nk) foud(0:mgx),fou(0:mx1,my,kb:ke)
      integer j,k,ix
      real(nk) sca
c
      sca=1.d0/dfloat(mgalx)
      do ix=0,mx1
         fou(ix,j,k)=sca*foud(ix)
      enddo
c
      return
      end subroutine pack_fftw
c-----|---------------------------------------------------------------|
      subroutine fill(fou,foud,j,k)
      use ctes
      use numbers
      implicit none 
c
c     /* inflates a complex array to 3/2 size, filling with zeros    */
c     /*                                                             */
c     /* input                                                       */
c     /* fou(0:mx1,my,mz)                                            */
c     /*                                                             */
c     /* output                                                      */
c     /* foud(0:mgx)                                                 */
c
c      include "ctes3D"
      complex(nk) foud(0:mgx),fou(0:mx1,my,kb:ke)
      integer j,k
c
c-----|---------------------------------------------------------------|     
      call zcopy(mx1+1,fou(0,j,k),1,foud(0),1)
      call zcopy(mgx-mx1,zero,0,foud(mx1+1),1)
c
      return
      end subroutine fill
c-----|---------------------------------------------------------------|
      subroutine derivx(u,du,ider,scal,what)
      use ctes
      use numbers
      implicit none
c
c     /***************************************************************/
c     /* computes 1st or 2nd derivatives in the x-fourier-direction  */
c     /*                                                             */
c     /* input                                                       */
c     /* u(0:my1,0:mz1,ib:ie)   array in c-c-f space, X-CUT          */
c     /* ider                   '1' or '2'     - order of derivative */
c     /*                        1<=ider<=2                           */
c     /* real*8 scal            REAL scalar factor                   */
c     /* what                   'n' - erase previous value of 'du'   */
c     /*                        [else] - add to previous values      */
c     /*                                                             */
c     /* output                                                      */
c     /* du(0:my1,0:mz1,ib:ie)  derivatives in c-c-f space           */
c     /*                                                             */
c     /* note: 'u' and 'du' CAN have the same memory address!        */
c     /***************************************************************/
      complex(nk) u(my,mz,ib:ie),du(my,mz,ib:ie)
      real(nk) scal
      integer ider,i,j,k
      character what      
c-----|---------------------------------------------------------------|
      if(what.eq.'n')then       !erase previous value of 'du'
         if(ider.eq.1)then
            do i=ib,ie
               do k=1,mz
                  do j=1,my
                     du(j,k,i)=cii*xalp(i)*u(j,k,i)*scal
                  enddo
               enddo
            enddo
         elseif(ider.eq.2)then
            do i=ib,ie
               do k=1,mz
                  do j=1,my
                     du(j,k,i)=-xalp(i)**2*u(j,k,i)*scal
                  enddo
               enddo
            enddo
         else
            write(*,*)'this order of deriv not implemented:',ider
            stop
         endif
      else                      !add to previous value
         if(ider.eq.1)then
            do i=ib,ie
               do k=1,mz
                  do j=1,my
                     du(j,k,i)=du(j,k,i)+cii*xalp(i)*u(j,k,i)*scal
                  enddo
               enddo
            enddo
         elseif(ider.eq.2)then
            do i=ib,ie
               do k=1,mz
                  do j=1,my
                     du(j,k,i)=du(j,k,i)-xalp(i)**2*u(j,k,i)*scal
                  enddo
               enddo
            enddo
         else
            write(*,*)'this order of deriv not implemented:',ider
            stop
         endif
      endif
c
      return
      end subroutine derivx
c-----|---------------------------------------------------------------|
      subroutine genexp
      use ctes, only:xalp,iax,mx1,mgalx,x,alp
      use numbers
      implicit none 
      integer i
c-----|---------------------------------------------------------------|
      do i=0,mx1
         iax(2*i+1) = i
         iax(2*i+2) = i
      enddo
c      
      do i=0,mx1
         xalp(i) = alp*i
      enddo
c
      do i=1,mgalx
         x(i)=dfloat(i-1)/dfloat(mgalx)*f2o1*pi/alp
      enddo
c     
      return
      end subroutine genexp
c-----|---------------------------------------------------------------|
