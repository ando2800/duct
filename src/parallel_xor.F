      subroutine chx2z(xy,xz,wk2)
c     /cccccccccccccccccccccccccccccccccccccccccccccccc/
c     MODULE FOR SP2 MPI uses SENDRECV               c/
c     sends a block (my,mz,ib:ie) to a (my,kb:ke,mx) c/
c     /cccccccccccccccccccccccccccccccccccccccccccccccc/
      use ctes, only: ibeg,iend,kbeg,kend,my,mz,mx1,
     &                numerop,myid,nsfou,mzp,mxp
      use timing
      implicit none
      include 'mpif.h'
      complex*16 xz(nsfou),xy(nsfou), wk2(my*mzp*mxp)
      integer i,ip_send,ip_recv,ierr,ipnodes,ipn
c     /*cccccccccccccccccccccccccccccccccccccccccccccc*/
c     /* NOTE: wk2 dimensioned at least (mxp*mzp*my)  */
c     /*cccccccccccccccccccccccccccccccccccccccccccccc*/

      integer istat(MPI_STATUS_SIZE)

      integer mmx,mmz,ntoto,iflag,ipoxz,ireqr,nsetotr,nsetots
c
      if(numerop.eq.1)return
c      
      prof_para(1) = prof_para(1) - MPI_WTIME()

      mmx=(iend(myid)-ibeg(myid)+1)*my
      mmz=(kend(myid)-kbeg(myid)+1)*my

      if (mod(numerop,3)==0) then
         ipnodes=numerop*4/3
      else
         ipnodes=numerop
      endif      
      
c      do i=0,numerop-1
      do ipn=0,ipnodes-1
         i=ieor(myid,ipn)
         if (i<numerop) then
            ip_send = i
            ip_recv = i
            call pacy2z(wk2,xy,ip_send)
         
            if(ip_send.ne.myid)then
            
               nsetots=mmx*(kend(ip_send)-kbeg(ip_send)+1)
               nsetotr=mmz*(iend(ip_recv)-ibeg(ip_recv)+1)
            
               ipoxz=1+mmz*(ibeg(ip_recv)) !zero-based index i!!!
               call MPI_SENDRECV(wk2,nsetots,MPI_DOUBLE_COMPLEX,
     .              ip_send,0,xz(ipoxz),nsetotr,
     .              MPI_DOUBLE_COMPLEX,ip_recv,0,
     .              MPI_COMM_WORLD,istat,ierr)
            else

               ntoto=mmz*(iend(myid)-ibeg(myid)+1)
               ipoxz=1+mmz*(ibeg(myid)) 
               xz(ipoxz:ipoxz+ntoto-1)=wk2(1:ntoto)

            endif
         endif
      enddo
      ntoto=mmz*(mx1+1)
      prof_para(1) = prof_para(1) + MPI_WTIME()
c
      prof_para(3) = prof_para(3) - MPI_WTIME()
c ! transpose
c ! fou(my,kb:ke,0:mx1) [xz] --> fou(0:mx1,my,kb:ke) [xy]
      call cyzx2xyz(xy,xz)
      prof_para(3) = prof_para(3) + MPI_WTIME()
c     
      return
      end
c----------------------------------------------------------|
      subroutine chz2x(xz,xy,wk2)
c     /ccccccccccccccccccccccccccccccccccccccccccccccccccc/
c     /* MODULE FOR SP2 MPI uses SENDRECV                */
c     /* sends a block (my,kb:ke,mx) to a (my,mz,ib:ie)  */
c     /* input:                                          */
c     /* complex*16 xz(my,kb:ke,0:mx1)                   */
c     /*                                                 */
c     /* output:                                         */
c     /* complex*16 xz(my,mz,ib:ie)                      */
c     /*                                                 */
c     /* work array:                                     */
c     /* complex*16 xy(my,mz,ib:ie)                      */
c     /* complex*16 wk2(my,mzp,mxp)                      */
c     /ccccccccccccccccccccccccccccccccccccccccccccccccccc/
      use ctes, only:nsfou,my,mz,mzp,mxp,numerop,myid,
     &               kbeg,kend,ibeg,iend
      use timing
      implicit none 
      include 'mpif.h'

      complex*16 xz(nsfou),xy(nsfou), wk2(my*mzp*mxp)
      integer i,ip_send,ip_recv,ierr,ipnodes,ipn

      integer istat(MPI_STATUS_SIZE)

      integer ipoxy,mmz,mmx,nsetots,nsetotr,ireqs,iflag,ntoto
c     
      if (numerop.eq.1) return
c     
      mmz=(kend(myid)-kbeg(myid)+1)*my
      mmx=(iend(myid)-ibeg(myid)+1)*my

c ! transpose 
      prof_para(4) = prof_para(4) - MPI_WTIME()
c  fou(my,kb:ke,0:mx1) [xy] <-- fou(0:mx1,my,kb:ke) [xz]
      call cxyz2yzx(xz,xy)
      prof_para(4) = prof_para(4) + MPI_WTIME()

      if (mod(numerop,3)==0) then
         ipnodes=numerop*4/3
      else
         ipnodes=numerop
      endif      
c
      prof_para(2) = prof_para(2) - MPI_WTIME()

c      do i=0,numerop-1
      do ipn=0,ipnodes-1
         i=ieor(myid,ipn)
         if (i<numerop) then
            !ip_send = mod(myid+i,numerop)
            !ip_recv = mod(numerop+myid-i,numerop)
            ip_send = i
            ip_recv = i
            if(ip_send.ne.myid)then
               nsetots=(iend(ip_send)-ibeg(ip_send)+1)*mmz
               nsetotr=(kend(ip_recv)-kbeg(ip_recv)+1)*mmx

               ipoxy=1+mmz*(ibeg(ip_send))             

               call MPI_SENDRECV(xy(ipoxy),nsetots,
     .              MPI_DOUBLE_COMPLEX,ip_send,0,wk2,
     .              nsetotr,MPI_DOUBLE_COMPLEX,
     .              ip_recv,0,MPI_COMM_WORLD,istat,ierr)
            
            else
               ipoxy=1+mmz*(ibeg(myid)) !zero-based index i!!!
               ntoto=mmz*(iend(myid)-ibeg(myid)+1)
               wk2(1:ntoto)=xy(ipoxy:(ipoxy+ntoto-1)) 
            
            endif

            call unpacz2y(wk2,xz,ip_recv)
         endif
      enddo
      
      ntoto=mmx*mz
      prof_para(2) = prof_para(2) + MPI_WTIME()

      return
      end
c-----|---------------------------------------------------------------|
      subroutine unpacz2y(xyi,xyo,iproc)
      use ctes, only: my,mz,myid,ibeg,iend,kbeg,kend
      implicit none 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     IN xzo OUT xzi                                           c
c     unpack from kb(iproc) till ke(iproc)                     c
c     and    from ib(myid)  till ie(myid)                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer iproc,i,k,krun
      complex*16 xyi(my,kend(iproc)-kbeg(iproc)+1,*),xyo(my,mz,*)

      do i=1,iend(myid)-ibeg(myid)+1
         do k=kbeg(iproc),kend(iproc)
            krun=k-kbeg(iproc)+1
            xyo(1:my,k,i)=xyi(1:my,krun,i)
         enddo
      enddo

      return
      end
c-----|---------------------------------------------------------------|
      subroutine pacy2z(xyo,xyi,iproc)
      use ctes, only: my,mz,myid,ibeg,iend,kbeg,kend
      implicit none 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     IN xzo OUT xzi                                           c
c     unpack from kb(iproc) till ke(iproc)                     c
c     and    from ib(myid)  till ie(myid)                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer i,k,krun,iproc
      complex*16 xyi(my,mz,*),xyo(my,kend(iproc)-kbeg(iproc)+1,*)

      do i=1,iend(myid)-ibeg(myid)+1
         do k=kbeg(iproc),kend(iproc)
            
            krun=k-kbeg(iproc)+1            
            xyo(1:my,krun,i)=xyi(1:my,k,i)
         enddo
      enddo

      return
      end

c---------------------------------------------------------------------
c   for physical space
c--------------------------------------------------------------------
      subroutine chx2z_fis(xy,xz,wk2)
c     /ccccccccccccccccccccccccccccccccccccccccccccccccc/
c     MODULE FOR SP2 MPI uses SENDRECV                 c/
c     sends a block (my,mz,ipb:ipe) to a (my,kpb:kpe,mx) c/
c     /c   kpb=kp, kpe=ke                   c/
c     /ccccccccccccccccccccccccccccccccccccccccccccccccc/
      use ctes
      implicit none
      include 'mpif.h'

      real(nk) xz(*),xy(*), wk2(*)
      integer i,ip_send,ip_recv
      integer iflag,ierr,ipoxz,ireqr,mmx,mmz,nsetotr,nsetots,ntoto
c     /*cccccccccccccccccccccccccccccccccccccccccccccc*/
c     /* NOTE: wk2 dimensioned at least (mxp*mzp*my)  */
c     /*cccccccccccccccccccccccccccccccccccccccccccccc*/
#ifdef BLOCKING_COMM
      integer istat(MPI_STATUS_SIZE)
#else /* BLOCKING_COMM */            
      integer istat((numerop-1)*2*MPI_STATUS_SIZE+1)
      integer irequestr((numerop-1)+1)
      integer irequests
#endif /* BLOCKING_COMM */            
c
      if(numerop.eq.1)return
#ifndef BLOCKING_COMM
      ireqr=0
      iflag=0
#endif /* BLOCKING_COMM */            
c     
      mmx=(ipend(myid)-ipbeg(myid)+1)*my
      mmz=(kpend(myid)-kpbeg(myid)+1)*my
      
      do 10 i=0,numerop-1
         ip_send = mod(myid+i,numerop)
         ip_recv = mod(numerop+myid-i,numerop)
#ifndef BLOCKING_COMM
         if(iflag.gt.0)then
            call MPI_Wait(irequests,istat,ierr)
         endif
#endif /* BLOCKING_COMM */            
         call pacy2z_fis(wk2,xy,ip_send)
         
         if(ip_send.ne.myid)then
            
            nsetots=mmx*(kpend(ip_send)-kpbeg(ip_send)+1)
            nsetotr=mmz*(ipend(ip_recv)-ipbeg(ip_recv)+1)
            
            ipoxz=1+mmz*(ipbeg(ip_recv)) !zero-based index i!!!
#ifdef BLOCKING_COMM
            call MPI_SENDRECV(wk2,nsetots,MPI_DOUBLE_PRECISION,
     .           ip_send,0,xz(ipoxz),nsetotr,
     .           MPI_DOUBLE_COMPLEX,ip_recv,0,
     .           MPI_COMM_WORLD,istat,ierr)
#else /* BLOCKING_COMM */            
            iflag=1
            call MPI_ISEND(wk2,nsetots,MPI_DOUBLE_PRECISION,
     .           ip_send,0,MPI_COMM_WORLD,irequests,ierr)
            ireqr=ireqr+1
            call MPI_IRECV(xz(ipoxz),nsetotr,MPI_DOUBLE_PRECISION,
     $           ip_recv,0,MPI_COMM_WORLD,irequestr(ireqr),ierr)
#endif /* BLOCKING_COMM */            
         else

            ntoto=mmz*(ipend(myid)-ipbeg(myid)+1)
            ipoxz=1+mmz*(ipbeg(myid)) 
            call dcopy(ntoto,wk2,1,xz(ipoxz),1)

         endif
 10   continue
      ntoto=mmz*(mgalx+2) 
c      ntoto=mmz*(mgalx)!!?
#ifndef BLOCKING_COMM
      call MPI_Waitall(ireqr,irequestr,istat,ierr)
      call MPI_Wait(irequests,istat,ierr)
#endif /* BLOCKING_COMM */            
      call dcopy(ntoto,xz,1,xy,1)
c     
      return
      end
c----------------------------------------------------------|
      subroutine chz2x_fis(xz,xy,wk2)
c     /ccccccccccccccccccccccccccccccccccccccccccccccccccc/
c     /* MODULE FOR SP2 MPI uses SENDRECV                */
c     /* sends a block (my,kpb:kpe,mx) to a (my,mz,ipb:ipe)*/
c     /* input:                                          */
c     /* real*8 xz(my,kpb:kpe,1:mgalx)                   */
c     /*                                                 */
c     /* output:                                         */
c     /* real*8 xz(my,mz,ipb:ipe)                        */
c     /*                                                 */
c     /* work array:                                     */
c     /* real*8 xy(my,mz,ipb:ipe)                      */
c     /* real*8 wk2(my,mzp,mxp)                      */
c     /ccccccccccccccccccccccccccccccccccccccccccccccccccc/
      use ctes
      implicit none 
      include 'mpif.h'

      real(nk) xz(*),xy(*), wk2(*)
      integer i,ip_send,ip_recv
      integer ireqs,iflag,mmz,mmx,ierr,ipoxz,nsetotr,nsetots,ntoto
#ifdef BLOCKING_COMM
      integer istat(MPI_STATUS_SIZE)
#else /* BLOCKING_COMM */            
      integer istat((numerop-1)*2*MPI_STATUS_SIZE+1)
      integer irequests((numerop-1)+1)
      integer irequestr
#endif /* BLOCKING_COMM */            
c
      if(numerop.eq.1)return
#ifndef BLOCKING_COMM
      ireqs=0
      iflag=0
#endif /* BLOCKING_COMM */            
c      
      mmz=(kpend(myid)-kpbeg(myid)+1)*my
      mmx=(ipend(myid)-ipbeg(myid)+1)*my

      do 10 i=0,numerop-1
         ip_send = mod(myid+i,numerop)
         ip_recv = mod(numerop+myid-i,numerop)
         
         if(ip_send.ne.myid)then
            nsetots=(ipend(ip_send)-ipbeg(ip_send)+1)*mmz
            nsetotr=(kpend(ip_recv)-kpbeg(ip_recv)+1)*mmx

            ipoxz=1+mmz*(ipbeg(ip_send))             
#ifdef BLOCKING_COMM
            call MPI_SENDRECV(xz(ipoxz),nsetots,
     .           MPI_DOUBLE_PRECISION,ip_send,0,wk2,
     .           nsetotr,MPI_DOUBLE_PRECISION,
     .           ip_recv,0,MPI_COMM_WORLD,istat,ierr)
#else /* BLOCKING_COMM */       
            iflag=1
            ireqs=ireqs+1
            call MPI_ISEND(xz(ipoxz),nsetots,MPI_DOUBLE_PRECISION,
     $           ip_send,0,MPI_COMM_WORLD,irequests(ireqs),ierr)
            call MPI_IRECV(wk2,      nsetotr,MPI_DOUBLE_PRECISION,
     .           ip_recv,0,MPI_COMM_WORLD,irequestr,ierr)
#endif /* BLOCKING_COMM */            
            
         else
            ipoxz=1+mmz*(ipbeg(myid)) !zero-based index i!!!
            ntoto=mmz*(ipend(myid)-ipbeg(myid)+1)
            call dcopy(ntoto,xz(ipoxz),1,wk2,1)
            
         endif

#ifndef BLOCKING_COMM
         if(iflag.gt.0)then
            call MPI_Wait(irequestr,istat,ierr)
         endif
#endif /* BLOCKING_COMM */            
         call unpacz2y_fis(wk2,xy,ip_recv)
         
 10   continue
      
      ntoto=mmx*mz
#ifndef BLOCKING_COMM
      call MPI_Waitall(ireqs,irequests,istat,ierr)
#endif /* BLOCKING_COMM */            
      call dcopy(ntoto,xy,1,xz,1)

      return
      end
c-----|---------------------------------------------------------------|
      subroutine unpacz2y_fis(xyi,xyo,iproc)
      use ctes
      implicit none 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     IN xzo OUT xzi                                           c
c     unpack from kb(iproc) till ke(iproc)                     c
c     and    from ib(myid)  till ie(myid)                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer iproc
      real(nk) xyi(my,kpend(iproc)-kpbeg(iproc)+1,*),xyo(my,mz,*)
      integer i,k,kprun,j

      do i=1,ipend(myid)-ipbeg(myid)+1
         do k=kpbeg(iproc),kpend(iproc)
            kprun=k-kpbeg(iproc)+1
            do j=1,my
               xyo(j,k,i)=xyi(j,kprun,i)
            enddo
         enddo
      enddo

      return
      end


      subroutine pacy2z_fis(xyo,xyi,iproc)
      use ctes
      implicit none 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     IN xzo OUT xzi                                           c
c     unpack from kb(iproc) till ke(iproc)                     c
c     and    from ib(myid)  till ie(myid)                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer iproc,kprun
      real(nk) xyi(my,mz,*),xyo(my,kpend(iproc)-kpbeg(iproc)+1,*)
      integer i,j,k

      do i=1,ipend(myid)-ipbeg(myid)+1
         do k=kpbeg(iproc),kpend(iproc)
            
            kprun=k-kpbeg(iproc)+1
            
            do j=1,my
               xyo(j,kprun,i)=xyi(j,k,i)

            enddo
         enddo
      enddo

      return
      end

c----------------------------------------------------------|
      subroutine cxyz2yzx(phys,phys2)
      use ctes
      implicit none

      integer i,j,k,kk
      complex*16 phys(0:mx1,my,kb:ke), phys2(my,kb:ke,0:mx1)

      do k=kb,ke
         do j=1,my
            do i=0,mx1
               phys2(j,k,i) = phys(i,j,k) 
            enddo
         enddo
      enddo

c      phys2 = transpose(phys)

c      do kk=1,mgalz,blockingik2ki
c         do i=1,mgx+1
c            do k=kk,kk+blockingik2ki-1
c            phys2(k,i) = phys(i,k) 
c            enddo
c         enddo
c      enddo

      end
    
c----------------------------------------------------------|
      subroutine cyzx2xyz(phys,phys2)
      use ctes
      implicit none

      integer i,j,k,kk
      complex*16 phys(0:mx1,my,kb:ke), phys2(my,kb:ke,0:mx1)

      do k=kb,ke
         do j=1,my
            do i=0,mx1
               phys(i,j,k) = phys2(j,k,i) 
            enddo
         enddo
      enddo

c      phys = transpose(phys2)

c      do kk=1,mgalz,blockingki2ik
c         do i=1,mgx+1
c            do k=kk,kk+blockingki2ik-1
c            phys(i,k) = phys2(k,i) 
c            enddo
c         enddo
c      enddo

      end
c----------------------------------------------------------|    
