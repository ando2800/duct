      subroutine cheb_gauss_lobatto_grid(xcgl,n)
      use ctes, only: nk
      use numbers
      implicit none 
c     /***************************************************************/
c     /* sets the chebyshev-gauss-lobatto grid                       */
c     /***************************************************************/
      integer n,i
      real(nk) xcgl(0:n)

c----|----------------------------------------------------------------|
      do i=0,n
         xcgl(i)=dcos(pi*dfloat(i)/dfloat(n))
      enddo
c
      return
      end subroutine cheb_gauss_lobatto_grid
c----|----------------------------------------------------------------|
      subroutine cheb_colloc_deriv1(xcgl,deriv1,nx)
      use ctes, only: nk
      use numbers
      implicit none 
c
c     /***************************************************************/
c     /* computes the chebyshev collocation derivative matrix for    */
c     /* gauss-lobatto collocation points.                           */
c     /*                                                             */
c     /* input                                                       */
c     /* x(0:nx)       chebyshev-gauss-lobatto points                */
c     /*                                                             */
c     /* output                                                      */
c     /* deriv1(0:nx,0:nx)    collocation derivative matrix          */
c     /***************************************************************/
c
      real(nk) xcgl(0:nx)
      real(nk) deriv1(0:nx,0:nx)
      real*8 cheb_cbar,cbari,cbarj
      integer nx,i,j
c----|----------------------------------------------------------------|
c     /* use the formula provided in cqhz, p. 69, (2.4.31) */
      do j=1,nx-1
         cbarj=cheb_cbar(j,nx)
         do i=0,nx
            cbari=cheb_cbar(i,nx)
            if(i.eq.j)then
               deriv1(i,j)=-xcgl(j)/(f2o1*(f1o1-xcgl(j)**2))
            else
               deriv1(i,j)=  
     &              (-f1o1)**(i+j)*cbari/(cbarj*(xcgl(i)-xcgl(j)))
            endif
         enddo
      enddo
      j=0
      cbarj=f2o1
      do i=0,nx
         cbari=cheb_cbar(i,nx)
         if(i.eq.j)then
            deriv1(i,j)=dfloat(2*nx**2+1)/f6o1
         else
            deriv1(i,j)=(-f1o1)**(i+j)*cbari/(cbarj*(xcgl(i)-xcgl(j)))
         endif
      enddo
      j=nx
      cbarj=f2o1
      do i=0,nx
         cbari=cheb_cbar(i,nx)
         if(i.eq.j)then
            deriv1(i,j)=-dfloat(2*nx**2+1)/f6o1
         else
            deriv1(i,j)=(-f1o1)**(i+j)*cbari/(cbarj*(xcgl(i)-xcgl(j)))
         endif
      enddo
c
      return
      end subroutine cheb_colloc_deriv1
c----|----------------------------------------------------------------|
      real*8 function cheb_cbar(i,n)
      use numbers
      implicit none 
      integer i,n
c

      if(i.eq.0.or.i.eq.n)then
         cheb_cbar=f2o1
      else
         cheb_cbar=f1o1
      endif
c     
      return
      end function cheb_cbar
c----|----------------------------------------------------------------|
