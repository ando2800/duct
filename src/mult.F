c     /***************************************************************/
c     /*                 computes either:                            */
c     /*                                                             */
c     /*        what = 'n':    s = (+/-) u.v                         */
c     /*        what = 'y':    s = s (+/-) u.v                       */
c     /*                                                             */
c     /*        plusm = ('+' or '-') controls the signs              */
c     /*                                                             */
c     /*  note: intended for fields in p-p-p or p-c-c space, Z-cut!  */
c     /***************************************************************/
      subroutine mult(ntot,u,v,s,what,plusm)
      implicit none
      real*8 u(*), v(*), s(*)
      integer ntot,i
      character what,plusm
      
      if (what.eq.'n') then
         
         if (plusm.eq.'+') then
            do i=1,ntot
               s(i) = u(i)*v(i)
            enddo
         else
            do i=1,ntot
               s(i) = - u(i)*v(i)
            enddo
         endif
      else

         if (plusm.eq.'+') then
            do i=1,ntot
               s(i) = s(i) + u(i)*v(i)
            enddo
         else
            do i=1,ntot
               s(i) = s(i) - u(i)*v(i)
            enddo
            
         endif
c
      endif
c
      return
      end subroutine mult
c/********************************************************************/
c/*   copies the scalar multiplication of vector a with s in array b */
c/********************************************************************/
      subroutine dscal2(n,a,nstr1,s,b,nstr2)
      implicit none
      integer n,nstr1,nstr2,i
      real*8 a(nstr1,*),b(nstr2,*),s
c
      do i=1,n
         b(1,i)=a(1,i)*s
      enddo
c
      return
      end
c-----|---------------------------------------------------------------|
