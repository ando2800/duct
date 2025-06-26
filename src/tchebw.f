!-----|---------------------------------------------------------------|
      subroutine tchebi_fftw(ngaly,ngalz)
      use tcheb2tmp
      implicit none
      include 'fftw3.f'
      integer ngaly,ngalz
! /********************************************************************
! /*     initialises coefficients for the full tchebichev transform
! /*
! /*      transforms to be of size -->     ngal
! /*      assumes:
! /*           ngal-1 --> even, multiple of 2 & 3 only and >3
! /*           ngal   <=  3*nmaxt
! /********************************************************************
      if (allocated(cdumy)) then
         call dfftw_destroy_plan(planr_tcheby)
         call dfftw_destroy_plan(planc_tcheby)
         call dfftw_cleanup()
         deallocate(cdumy,cdumy2)
      end if
  
      if (allocated(cdumz)) then
         call dfftw_destroy_plan(planr_tchebz)
         call dfftw_destroy_plan(planc_tchebz)
         call dfftw_cleanup()
         deallocate(cdumz,cdumz2)
      end if

      allocate(cdumy(ngaly),cdumy2(ngaly))
      cdumy=0.d0; cdumy2=0.d0
      nny=ngaly; 
      nnyc=0
      if (nnyc.ne.0) then
         write(*,*) 'de-aliasing in y-dir, set zero:',nnyc+1,'to',ngaly
      endif
      dny1=1.d0/dfloat(ngaly-1)
      dny2=0.5d0*dny1
      call dfftw_plan_r2r_1d(planr_tcheby,nny,cdumy,cdumy, 
     $     FFTW_REDFT00, FFTW_MEASURE)
      call dfftw_plan_r2r_1d(planc_tcheby,nny,cdumy2,cdumy2, 
     $     FFTW_REDFT00, FFTW_MEASURE)
  
      allocate(cdumz(ngalz),cdumz2(ngalz))
      cdumz=0.d0; cdumz2=0.d0
      nnz=ngalz; 
      nnzc=0 
      if (nnzc.ne.0) then
         write(*,*) 'de-aliasing in z-dir, set zero:',nnzc+1,'to',ngalz
      endif
      dnz1=1.d0/dfloat(ngalz-1)
      dnz2=0.5d0*dnz1
      
      call dfftw_plan_r2r_1d(planr_tchebz,nnz,cdumz,cdumz, 
     $     FFTW_REDFT00, FFTW_MEASURE)
      call dfftw_plan_r2r_1d(planc_tchebz,nnz,cdumz2,cdumz2, 
     $     FFTW_REDFT00, FFTW_MEASURE)
   
      return
      
      end subroutine tchebi_fftw
!-----|---------------------------------------------------------------|
      subroutine tcheb_fftw(cr,ise,isa,m,cdir,iopt)
      use tcheb2tmp
      implicit none
!*********************************************************************/
!*   given  m    ==> real vectors cm(1:ngal), stored in c as
!*          ise  ==> stride between elements in each vectors
!*          isa  ==> stride between initial elements of vectors
!*
!*   each of which is a function defined at
!
!       c(j+1) = c( cos(j*pi/(ngal-1)) )
!
!*   computes chebichev expansion
!*      c(x) = sum from j=0 to (ngal-1) of c(j+1) * tp(j,x)
!*   (or viceversa)
!
!*         iopt >= 0  ==> inverse transform (x chebich --> c phys.)
!*         iopt <  0  ==> direct  transform (x physic. --> c cheb.)
!*
!*             jimenez/    may 90
!*********************************************************************/
      real(nk) cr(*)
      integer ise,isa,m,nn,nnc
      real*8 dn1,dn2
      integer i,j,k,j0,j00
      character*1 cdir
      integer iopt
!------- coefficients and savearea for chebichev expansions -----------
      if (cdir.eq.'y') then 
         nn = nny               !aux(1)
         dn1= dny1              !aux(2)
         dn2= dny2              !aux(3)
         nnc= nnyc         
!
         j00 = 1
         do k=1,m
!  *********  copy to contiguous ******
            j0 = j00
            do j=1,nn            
               cdumy(j)=cr(j0)  ! aux(3+j) = c(j0)
               j0 = j0+ise
            enddo
!-------   forward transform  (to chebichev) -------------
            if (iopt.lt.0) then
!
               call dfftw_execute_r2r(planr_tcheby,cdumy,cdumy)
               if (nnc.ne.0) then
                  cdumy((nnc+1):nn) =0.d0 ! dealiasing
               endif
! *********  copy to original ********
               j0 = j00+ise
               do j=2,nn-1
                  cr(j0) = dn1*cdumy(j) !aux(3+j)
                  j0 = j0+ise
               end do
               cr(j00) =dn2*cdumy(1) !aux(4)
               cr(j0)  =dn2*cdumy(nn) !aux(3+nn)
!----- inverse transform  (to physical) --------------
            else
               ! set fill0 1/2 modes for dealiasing
               if (nnc.ne.0) then              
                  cdumy((nnc+1):nn) = 0.d0 ! dealiasing
               endif
               cdumy(1) = 2.d0*cdumy(1) !aux(4)    = 2.d0*aux(4)
               cdumy(nn) = 2.d0*cdumy(nn) !aux(3+nn) = 2.d0*aux(3+nn)                    
               call dfftw_execute_r2r(planr_tcheby,cdumy,cdumy)
!  *********  copy to original ********
               j0 = j00
               do j=1,nn
                  cr(j0) = 0.5d0*cdumy(j) !  normalized by 0.5 for using fftw3 DCT
                  j0 = j0+ise
               end do
            endif
            j00 = j00+isa
         enddo
     
      elseif (cdir.eq.'z') then
         nn = nnz               !aux(1)
         dn1= dnz1              !aux(2)
         dn2= dnz2              !aux(3)
         nnc= nnzc
!
         j00 = 1
         do k=1,m
!  *********  copy to contiguous ******
            j0 = j00
            do j=1,nn            
               cdumz(j)=cr(j0)  ! aux(3+j) = c(j0)
               j0 = j0+ise
            enddo
!-------   forward transform  (to chebichev) -------------
            if (iopt.lt.0) then

               call dfftw_execute_r2r(planr_tchebz,cdumz,cdumz)
               if (nnc.ne.0) then               
                  cdumz((nnc+1):nn) =0.d0 ! dealiasing
               endif
! *********  copy to original ********
               j0 = j00+ise
               do j=2,nn-1
                  cr(j0) = dn1*cdumz(j) !aux(3+j)
                  j0 = j0+ise
               end do
               cr(j00) =dn2*cdumz(1) !aux(4)
               cr(j0)  =dn2*cdumz(nn) !aux(3+nn)
!----- inverse transform  (to physical) --------------
            else
               ! set fill0 1/2 modes for dealiasing
               if (nnc.ne.0) then                  
                  cdumz((nnc+1):nn) =0.d0 ! dealiasing
               endif
               cdumz(1)  = 2.d0*cdumz(1) !aux(4)    = 2.d0*aux(4)
               cdumz(nn) = 2.d0*cdumz(nn) !aux(3+nn) = 2.d0*aux(3+nn)               
               call dfftw_execute_r2r(planr_tchebz,cdumz,cdumz)
!  *********  copy to original ********
               j0 = j00
               do j=1,nn
                  cr(j0) = 0.5d0*cdumz(j) !  normalized by 0.5 for using fftw3 DCT
                  j0 = j0+ise
               end do
            endif
            j00 = j00+isa
         enddo
      else
         write(*,*) 'tcheb cdir is not set, stop' 
         stop
      end if
      
      return

      end subroutine 
!-----|---------------------------------------------------------------|
      subroutine tcheb2_fftw(cc,ise,isa,m,cdir,iopt)
      use tcheb2tmp
      implicit none
      include "fftw3.f"
!*********************************************************************/
!*   given  m    ==> real vectors cm(1:ngal), stored in c as
!*          ise  ==> stride between elements in each vectors
!*          isa  ==> stride between initial elements of vectors
!*
!*   each of which is a function defined at
!
!       c(j+1) = c( cos(j*pi/(ngal-1)) )
!
!*   computes chebichev expansion
!*      c(x) = sum from j=0 to (ngal-1) of c(j+1) * tp(j,x)
!*   (or viceversa)
!
!*         iopt >= 0  ==> inverse transform (x chebich --> c phys.)
!*         iopt <  0  ==> direct  transform (x physic. --> c cheb.)
!*
!*             jimenez/    may 90
!*********************************************************************/
      real*8 cc(2,*)
      integer ise,isa,m,nn,dn1,dn2,nnc
      integer i,j,k,j0,j00
      integer iopt
      character*1 cdir
!------- coefficients and savearea for chebichev expansions -----------
      ! this is buggy 
      if (cdir.eq.'y') then 
      nn = nny  !aux(1)
      dn1= dny1 !aux(2)
      dn2= dny2 !aux(3)
      nnc= nnyc
      !
      j00 = 1
      do k=1,m
c        !  *********  copy to contiguous ******
        j0 = j00
        do j=1,nn            
           cdumy(j)=cc(1,j0) ! aux(3+j) = c(j0)
           cdumy2(j)=cc(2,j0) ! aux(3+j) = c(j0)
           j0 = j0+ise
        enddo
c        !-------   forward transform  (to chebichev) -------------
        if (iopt.lt.0) then
           !
           call dfftw_execute_r2r(planr_tcheby,cdumy,cdumy)
           call dfftw_execute_r2r(planc_tcheby,cdumy2,cdumy2)
           if (nnc.ne.0) then
              cdumy((nnc+1):nn) =0.d0 ! dealiasing
              cdumy2((nnc+1):nn) =0.d0 ! dealiasing
           endif
c           ! *********  copy to original ********
           j0 = j00+ise
           do j=2,nn-1
              cc(1,j0) = dn1*cdumy2(j) !aux(3+j)
              cc(2,j0) = dn1*cdumy2(j) !aux(3+j)
              j0 = j0+ise
           end do
           cc(1,j00) =dn2*cdumy(1) !aux(4)
           cc(2,j00) =dn2*cdumy2(1) !aux(4)
           cc(1,j0)  =dn2*cdumy(nn) !aux(3+nn)
           cc(2,j0)  =dn2*cdumy2(nn) !aux(3+nn)
c           !----- inverse transform  (to physical) --------------
        else
c           ! set fill0 1/2 modes for dealiasing
           if (nnc.ne.0) then              
              cdumy((nnc+1):nn) = 0.d0 ! dealiasing
              cdumy2((nnc+1):nn) = 0.d0 ! dealiasing
           endif
           cdumy(1) = 2.d0*cdumy(1) !aux(4)    = 2.d0*aux(4)
           cdumy(nn) = 2.d0*cdumy(nn)!aux(3+nn) = 2.d0*aux(3+nn)                    
           cdumy2(1) = 2.d0*cdumy2(1) !aux(4)    = 2.d0*aux(4)
           cdumy2(nn) = 2.d0*cdumy2(nn) !aux(4)    = 2.d0*aux(4)

           call dfftw_execute_r2r(planr_tcheby,cdumy,cdumy)
           call dfftw_execute_r2r(planc_tcheby,cdumy2,cdumy2)
c           !  *********  copy to original ********
           j0 = j00
           do j=1,nn
              cc(1,j0) = 0.5d0*cdumy(j) !  normalized by 0.5 for using fftw3 DCT
              cc(2,j0) = 0.5d0*cdumy2(j) !  normalized by 0.5 for using fftw3 DCT   
c              ! c(j0) =aux(3+j) 
              j0 = j0+ise
           end do
        endif
         j00 = j00+isa
      enddo
     
      elseif (cdir.eq.'z') then
      nn = nnz  !aux(1)
      dn1= dnz1 !aux(2)
      dn2= dnz2 !aux(3)
      nnc= nnzc
c     !
      j00 = 1
      do k=1,m
c        !  *********  copy to contiguous ******
        j0 = j00
        do j=1,nn            
           cdumz(j)=cc(1,j0) ! aux(3+j) = c(j0)
           cdumz2(j)=cc(2,j0) ! aux(3+j) = c(j0)
           j0 = j0+ise
        enddo
c        !-------   forward transform  (to chebichev) -------------
        if (iopt.lt.0) then
           call dfftw_execute_r2r(planr_tchebz,cdumz,cdumz)
           call dfftw_execute_r2r(planc_tchebz,cdumz2,cdumz2)
           if (nnc.ne.0) then               
              cdumz((nnc+1):nn) =0.d0 ! dealiasing
              cdumz2((nnc+1):nn)=0.d0 ! dealiasing
           endif
c           ! *********  copy to original ********
           j0 = j00+ise
           do j=2,nn-1
              cc(1,j0) = dn1*cdumz(j)
              cc(2,j0) = dn1*cdumz2(j)
              j0 = j0+ise
           end do
           cc(1,j00) =dn2*cdumz(1) !aux(4)
           cc(2,j00) =dn2*cdumz2(1) !aux(4)
           cc(1,j0)  =dn2*cdumz(nn) !aux(3+nn)
           cc(2,j0)  =dn2*cdumz2(nn) !aux(3+nn)
c           !----- inverse transform  (to physical) --------------
        else
c           ! set fill0 1/2 modes for dealiasing
           if (nnc.ne.0) then
              ! dealiasing
              cdumz((nnc+1):nn) =0.d0
              cdumz2((nnc+1):nn) =0.d0
           endif
           cdumz(1)  = 2.d0*cdumz(1) !aux(4)    = 2.d0*aux(4)
           cdumz(nn) = 2.d0*cdumz(nn) !aux(3+nn) = 2.d0*aux(3+nn)
           cdumz2(1)  = 2.d0*cdumz2(1) !aux(4)    = 2.d0*aux(4)
           cdumz2(nn) = 2.d0*cdumz2(nn) !aux(3+nn) = 2.d0*aux(3+nn)

           call dfftw_execute_r2r(planr_tchebz,cdumz,cdumz)
           call dfftw_execute_r2r(planc_tchebz,cdumz2,cdumz2)
c           !  *********  copy to original ********
           j0 = j00
           do j=1,nn
              cc(1,j0) = 0.5d0*cdumz(j) !  normalized by 0.5 for using fftw3 DCT
              cc(2,j0) = 0.5d0*cdumz2(j) !  normalized by 0.5 for using fftw3 DCT
              j0 = j0+ise
           end do
        endif
        j00 = j00+isa
      enddo
      else
         write(*,*) 'tcheb2 cdir is not set, stop' 
         stop
      end if
      
      return
      
      end subroutine tcheb2_fftw

