#define DEBUG_read1
#undef DEBUG_read1
#define DEBUG_read2
#undef DEBUG_read2
#define DEBUG_read3
#undef DEBUG_read3
#define DEBUG_tmp1 
#undef DEBUG_tmp1
c
      subroutine read_fft_4fields_fis(a,b,c,d,tmpfis,tmpzxcut,tmpx,
     $     rprecision,cflag)
      use ctes
      use numbers
      implicit real*8 (a-h,o-z)
c     
c     /* directly reads 4 3d fields in binary format in p-space      */
c     /* decomposited in z-dir, (my,kb:ke,mgalx)                     */
c     /*  & pass to fourier. (this saves some buffers).              */
c     /* the header line contains dimensions & physical params       */
c     /*                                                             */
c     /* rprecision:  'r8'  -  read 8byte data                       */
c     /*              'r4'   -  read 4byte data and convert to 8     */
c     /*                                                             */
c     /* cflag        't'   -  read only fifth field [TEMPERATURE]   */
c     /*                        & close unit                         */
c     /*              [else]   read first four fields  (if defined   */
c     /*                       TEMPERARTURE, leave unit open)        */
c     
      include "mpif.h"
      complex*16 a(my,kb:ke,0:mx1),b(my,kb:ke,0:mx1),c(my,kb:ke,0:mx1),
     $     d(my,kb:ke,0:mx1)    !fou-size (output)
      dimension tmpfis(nsfis) !fis-size (tmp buffer)
      complex*16 tmpzxcut(my*mzp*mxp) !z-cut & x-cut pencil work array
      dimension tmpx(mgalx+2)   !1d fis--size x work array
      character*1 cflag
      parameter(iout=27)
      character*256 filename
      character*4 ext
      character*2 rprecision
c     
      parameter(nparamx=15)
      dimension param(nparamx)
      integer istat(MPI_STATUS_SIZE)
c-----|---------------------------------------------------------------|      
!      write(ext,'(i4.4)') iinp
!#ifdef TEMPERATURE
!      call concat4b(filinp,'uvwpt_',6,ext,4,'.bin',4,filename)
!#else /* TEMPERATURE */
!      call concat4b(filinp,'uvwp_',5,ext,4,'.bin',4,filename)
!#endif /* TEMPERATURE */
      filename=filinp      
      if(myid.eq.master) then
#ifdef TEMPERATURE
         nparam=13
         if(cflag.ne.'t')then
#else
            nparam=7
#endif
            write(*,*)'for 8-byte beginning of record mark!!!!!'
            write(*,*)'file will be over written '
            if(rprecision.eq.'r8')then
               nrecl=8
               open(iout,file=filename,form="unformatted",
     $              access='direct',recl=nrecl)
               write(*,*)'reading in p-space from ',trim(filename)
               do ii=1,nparam
                  irec=1+ii
                  read(iout,rec=irec,iostat=iread)param(ii)
                  if(iread.ne.0) then
                     write(*,*)'reading error !!! iread=',iread
                     stop
                  endif
               enddo
               irec=irec+2
               mgalxr=int(param(1))
               myr=int(param(2))
               mzr=int(param(3))
               if(mgalxr.ne.mgalx.or.myr.ne.my.or.mzr.ne.mz)then
                  write(*,*)'mismatching dimensions!',
     $                 mgalxr,myr,mzr
                  stop
               endif
               time=param(6)
               write(*,*)'physical params start field: alp,aspect,fnu',
     $              param(4),param(5),param(7)
#ifdef TEMPERATURE
               write(*,*)'physical params start field: alphat,gravx,',
     $              'gravy,gravz,grashof,prandtl',
     $              (param(ii),ii=8,13)
#endif
               mmxy=mgalx*my
               do ifield=1,4
                  ifstrec=irec+(2+mmxy*mz)*(ifield-1)
                  call do_read_direct8(tmpfis,1,mgalx,1,my,kb,ke,iout,
     $                 ifstrec)
                  if(ifield.eq.1) call fourx(a,tmpfis,tmpx,kb,ke,-1) !fis->fou
                  if(ifield.eq.2) call fourx(b,tmpfis,tmpx,kb,ke,-1) !fis->fou
                  if(ifield.eq.3) call fourx(c,tmpfis,tmpx,kb,ke,-1) !fis->fou
                  if(ifield.eq.4) call fourx(d,tmpfis,tmpx,kb,ke,-1) !fis->fou
                  do id = 1, lastid
                     call do_read_direct8(tmpfis,1,mgalx,1,my,
     $                    kbeg(id),kend(id),iout,ifstrec)
                     len=mmxy*(kend(id)-kbeg(id)+1)              
                     call MPI_SEND(tmpfis,len,MPI_DOUBLE_PRECISION,
     $                    id,1,MPI_COMM_WORLD,ierr)
                  enddo 
               enddo 
            else
               write(*,*)'precision not implemented in read: ',rprecision
               stop
            endif
#ifdef TEMPERATURE
         else
            ifield=5
            irec=16
            mmxy=mgalx*my
            ifstrec=irec+(2+mmxy*mz)*(ifield-1)
            call do_read_direct8(tmpfis,1,mgalx,1,my,kb,ke,iout,
     $           ifstrec)
            call fourx(a,tmpfis,tmpx,kb,ke,-1) !fis->fou
            do id = 1, lastid
               call do_read_direct8(tmpfis,1,mgalx,1,my,
     $              kbeg(id),kend(id),iout,ifstrec)
               len=mmxy*(kend(id)-kbeg(id)+1)              
               call MPI_SEND(tmpfis,len,MPI_DOUBLE_PRECISION,
     $              id,1,MPI_COMM_WORLD,ierr)
            enddo
            close(iout)
         endif
#else 
         close(iout)
#endif /* TEMPERATURE */
c
      else /* the other id */
         len=mgalx*my*(ke-kb+1)
#ifdef TEMPERATURE
         if(cflag.ne.'t')then
#endif
            do ifield=1,4
               call MPI_RECV(tmpfis,len,MPI_DOUBLE_PRECISION,
     &              master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
               if(ifield.eq.1) call fourx(a,tmpfis,tmpx,kb,ke,-1) !fis->fou
               if(ifield.eq.2) call fourx(b,tmpfis,tmpx,kb,ke,-1) !fis->fou
               if(ifield.eq.3) call fourx(c,tmpfis,tmpx,kb,ke,-1) !fis->fou
               if(ifield.eq.4) call fourx(d,tmpfis,tmpx,kb,ke,-1) !fis->fou
            enddo
#ifdef TEMPERATURE
         else 
            call MPI_RECV(tmpfis,len,MPI_DOUBLE_PRECISION,
     &           master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
            call fourx(a,tmpfis,tmpx,kb,ke,-1) !fis->fou 
         endif
#endif
      endif
c   
#ifdef TEMPERATURE
      if(cflag.ne.'t')then
#endif
         call MPI_BCAST(time,1,MPI_DOUBLE_PRECISION,master,
     $        MPI_COMM_WORLD,ierr) 
c
#ifdef DEBUG_read1 
      write(*,*) myid,'debug_read1: writing field'
      call debugwritefouz(a,kb,ke,myid,30)
      call debugwritefouz(b,kb,ke,myid,40)
      call debugwritefouz(c,kb,ke,myid,50)
      call debugwritefouz(d,kb,ke,myid,60)
      stop
#endif
c     /* shuffle arrays from z-cut to x-cut */
         call chz2x(a,tmpfis,tmpzxcut,myid)
         call chz2x(b,tmpfis,tmpzxcut,myid)
         call chz2x(c,tmpfis,tmpzxcut,myid)
         call chz2x(d,tmpfis,tmpzxcut,myid)
c
#ifdef TEMPERATURE
      else 

#ifdef DEBUG_tmp1 
      write(*,*) myid,'writing temp field in write.F'
      call debugwritefouz(a,kb,ke,myid,50)
      write(*,*) 'read temp field successfully! '
#endif
c
         call chz2x(a,tmpfis,tmpzxcut,myid)
      endif
#endif
c
      return
      end subroutine read_fft_4fields_fis
c---------------------------------------------------------------------|
      subroutine do_read_direct8(u,iff,ill,jff,jll,kff,kll,iunit,ioff)
      use ctes,only:my,mz
      implicit real*8 (a-h,o-z)

      dimension u(jff:jll,kff:kll,iff:ill)

c---------------------------------------------------------------------|
      do i=iff,ill
         do k=kff,kll
            do j=jff,jll
               irec=(i-1)*my*mz+(k-1)*my+j+ioff
c               write(*,*)j,k,i,irec,ioff
               read(iunit,rec=irec,iostat=iread) u(j,k,i)
c               write(22,*)j,k,i, u(j,k,i)
            enddo
         enddo
      enddo
      if(iread.ne.0) then
         write(*,*)'reading error !!!',iread
         stop
      endif
c     
      return
      end subroutine do_read_direct8
c---------------------------------------------------------------------|
#ifdef TEMPERATURE
      subroutine read_5fields_fou_apart(u,v,w,p,t,filebase,nprocr, 
     &     rprecision)
#else
      subroutine read_4fields_fou_apart(u,v,w,p,filebase,nprocr, 
     &     rprecision)
#endif
      use ctes
      implicit real*8 (a-h,o-z)
c
c     /* reads local 3d fields in binary format by master processor  */
c     /* and send to each processor                                  */
c     /* the master must be 0                                        */
c     /* the filenames are per-processor and have "myid" as extension*/
c     /* the header line contains dimensions & physical params       */
c     /*                                                             */
c     /* rprecision:  'r8'  -  read 8byte data                       */
c     /*              'r4'   -  read 4byte data and convert to 8     */
c     /* cflag        't'   -  read only fifth field [TEMPERATURE]   */
c     /*                        & close unit                         */
c     /*              [else]   read first four fields  (if defined   */
c     /*                       TEMPERARTURE, leave unit open)        */
c     /*                                                             */
c     
      include "mpif.h"
      complex*16 u(nsfou),v(nsfou),w(nsfou),p(nsfou)
      complex*16,dimension(:,:,:),allocatable::ur,vr,wr,pr
#ifdef TEMPERATURE
      complex*16 t(nsfou)
      complex*16,dimension(:,:,:),allocatable::tr
#endif
      integer nprocr
      integer ibegr(0:nprocr-1),iendr(0:nprocr-1)
      character*1 cflag
      parameter(iout=27)
      character*256 filebase
      character*4 ext,ext2
      character*2 rprecision
      integer istat(MPI_STATUS_SIZE)
c-----|---------------------------------------------------------------|
      if(master.ne.0) then
         write(*,*)'read_field_fou_apart: master must be 0'
         stop
      endif

      if(myid.eq.master) then
c     /* 'filename check' */
         ispace=index(filebase,' ')
         ext=filebase(ispace-4:ispace-1)
         if(ext.eq.'.bin'.or.filebase(ispace-1:ispace-1).eq.'.') then
            write(*,*) 'please remove extension *.bin from filinp'
            stop
         endif
c     /* read file apart */
         call pointer1D(0,mx1,nprocr,ibegr,iendr)
         allocate(ur(my,mz,0:mx1),vr(my,mz,0:mx1),  
     &        wr(my,mz,0:mx1),pr(my,mz,0:mx1))
         ur=0.d0; vr=0.d0; wr=0.d0; pr=0.d0
#ifdef TEMPERATURE
         allocate(tr(my,mz,0:mx1))
         tr=0.d0
#endif         
!     master read all data
         do id=0,nprocr-1       ! lastid            
            call read_4fields_fou_local(ur(1,1,ibegr(id)),
     &           vr(1,1,ibegr(id)),wr(1,1,ibegr(id)),
     &           pr(1,1,ibegr(id)),id,
     &           ibegr(id),iendr(id),filebase,rprecision,'_')
#ifdef TEMPERATURE
            call read_4fields_fou_local(tr(1,1,ibegr(id)), 
     &           dum,dum,dum,id,
     &           ibegr(id),iendr(id),filebase,rprecision,'t')
#endif
         enddo
         ! master copy to local 
         len=my*mz*(ie-ib+1)
         call zcopy(len,ur(1,1,0),1,u(1),1) ! ib should be 0
         call zcopy(len,vr(1,1,0),1,v(1),1) ! ib should be 0
         call zcopy(len,wr(1,1,0),1,w(1),1) ! ib should be 0
         call zcopy(len,pr(1,1,0),1,p(1),1) ! ib should be 0
#ifdef TEMPERATURE
         call zcopy(len,tr(1,1,0),1,t(1),1) ! ib should be 0
#endif 
         do id=1,nproc-1

            len=my*mz*(iend(id)-ibeg(id)+1)
            call MPI_SEND(ur(1,1,ibeg(id)),len,MPI_DOUBLE_COMPLEX,
     $           id,1,MPI_COMM_WORLD,ierr)
            call MPI_SEND(vr(1,1,ibeg(id)),len,MPI_DOUBLE_COMPLEX,
     $           id,1,MPI_COMM_WORLD,ierr)
            call MPI_SEND(wr(1,1,ibeg(id)),len,MPI_DOUBLE_COMPLEX,
     $           id,1,MPI_COMM_WORLD,ierr)
            call MPI_SEND(pr(1,1,ibeg(id)),len,MPI_DOUBLE_COMPLEX,
     $           id,1,MPI_COMM_WORLD,ierr)
#ifdef TEMPERATURE
            call MPI_SEND(tr(1,1,ibeg(id)),len,MPI_DOUBLE_COMPLEX,
     $           id,1,MPI_COMM_WORLD,ierr)
#endif
         enddo
         deallocate(ur,vr,wr,pr)
#ifdef TEMPERATURE
         deallocate(tr)
#endif
      else 

         len=my*mz*(ie-ib+1)
         call MPI_RECV(u,len,MPI_DOUBLE_COMPLEX,
     &        master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(v,len,MPI_DOUBLE_COMPLEX,
     &        master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(w,len,MPI_DOUBLE_COMPLEX,
     &        master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
         call MPI_RECV(p,len,MPI_DOUBLE_COMPLEX,
     &        master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr)
#ifdef TEMPERATURE
         call MPI_RECV(t,len,MPI_DOUBLE_COMPLEX,
     &        master,MPI_ANY_TAG,MPI_COMM_WORLD,istat,ierr) 
#endif
      endif
      
      return
#ifdef TEMPERATURE
      end subroutine read_5fields_fou_apart
#else
      end subroutine read_4fields_fou_apart
#endif
c---------------------------------------------------------------------|
      subroutine read_4fields_fou_local(a,b,c,d,id,ib,ie,
     $     filebase,rprecision,cflag)
      use ctes,only:my,mz,mgalx
      use running,only:time
      implicit real*8 (a-h,o-z)
c
c     /* reads local 3d fields in binary format by master processor  */
c     /* the filenames are per-processor and have "id" as extension*/
c     /* the header line contains dimensions & physical params       */
c     /*                                                             */
c     /* rprecision:  'r8'  -  read 8byte data                       */
c     /*              'r4'   -  read 4byte data and convert to 8     */
c     /* cflag        't'   -  read only fifth field [TEMPERATURE]   */
c     /*                        & close unit                         */
c     /*              [else]   read first four fields  (if defined   */
c     /*                       TEMPERARTURE, leave unit open)        */
c     /*                                                             */
c     
      complex*16 a(my,mz,ib:ie),b(my,mz,ib:ie),c(my,mz,ib:ie),
     $     d(my,mz,ib:ie)
      character*1 cflag
      parameter(iout=27)
      character*256 filename,filebase
      character*4 ext,ext2
      character*2 rprecision
c
      parameter(nparamx=15)
      dimension param(nparamx), pllel(4)
c-----|---------------------------------------------------------------|
#ifdef TEMPERATURE
      nparam=13
      if(cflag.ne.'t')then
#else
         nparam=7
#endif
c     
c     /* 'filename check' */
         ispace=index(filebase,' ')
         ext=filebase(ispace-4:ispace-1)
         if(ext.eq.'.bin'.or.filebase(ispace-1:ispace-1).eq.'.') then
            write(*,*) 'please remove extension *.bin from filinp'
            stop
         endif
c     /* each processor read local field file uvpw_****.id    */
         write(ext2,'(i4.4)') id
         call concat3b(filebase,'.',1,ext2,4,filename)
c     
         open(iout,file=filename,form="unformatted")
         if ((myid.eq.matster).and.(ib.eq.0)) then 
            write(*,*)'reading in fou-space from ',
     $        trim(filename),' extension is id'
         endif
         if(rprecision.eq.'r8')then
            read(iout)(pllel(ii),ii=1,4)
            read(iout)(param(ii),ii=1,nparam)
            mgalxr=int(param(1))
            myr=int(param(2))
            mzr=int(param(3))
            if(mgalxr.ne.mgalx.or.myr.ne.my.or.mzr.ne.mz)then
               write(*,*) id,'mismatching dimensions!',
     $              mgalxr,myr,mzr
               stop
            endif
            numeropr=int(pllel(1))
            myidr=int(pllel(2))
            ibr=int(pllel(3))
            ier=int(pllel(4))
            !if(numeropr.ne.numerop)then
            !   write(*,*)'mismatching number of processors!',numeropr
            !   stop
            !endif
            !if(myidr.ne.myid.or.ibr.ne.ib.or.ier.ne.ie)then
            !   write(*,*)myid,'re-setting myid! ',myidr
            !   myid=myidr
            !   ib=ibr
            !   ie=ier
            !   write(*,*) 'this mode does not work, myrank.ne.myid'
            !   stop
            !endif

            if(id.eq.0) then 
               time=param(6)
               write(*,*)
     $              'physical params: alp,aspect,fnu',
     $              param(4),param(5),param(7)
            endif
#ifdef TEMPERATURE
            if(id.eq.0) write(*,*)
     $           'physical params start field: alphat,gravx,',
     $           'gravy,gravz,grashof,prandtl',
     $           (param(ii),ii=8,13)
#endif
            read(iout)(((a(j,k,i),j=1,my),k=1,mz),i=ib,ie)
            read(iout)(((b(j,k,i),j=1,my),k=1,mz),i=ib,ie)
            read(iout)(((c(j,k,i),j=1,my),k=1,mz),i=ib,ie)
            read(iout)(((d(j,k,i),j=1,my),k=1,mz),i=ib,ie)
         else
            write(*,*)'precision not implemented in read: ',rprecision
            stop
         endif
#ifdef TEMPERATURE
      else
         read(iout)(((a(j,k,i),j=1,my),k=1,mz),i=ib,ie)
         close(iout)
      endif
#else
      close(iout)
#endif
c
      return
      end subroutine read_4fields_fou_local
c---------------------------------------------------------------------|
      subroutine write_statistics(iout)
      use ctes
      use statistics
      implicit real*8 (a-h,o-z)
      include "mpif.h"
c     /* tmp array for MPI_REDUCE */
      dimension ppt(my,mz)
      dimension uust(0:mx1,nspec),vvst(0:mx1,nspec),
     $     wwst(0:mx1,nspec),ppst(0:mx1,nspec)
#ifdef TEMPERATURE
      dimension ttt(my,mz),ttst(0:mx1,nspec)
#endif
c     /* tmp arrays for MPI_reduce 2nd moment statistics of omg */
      dimension omx2t(my,mz),omy2t(my,mz),omz2t(my,mz)
      dimension rcvcnt(0:numerop-1)
c
c     /* local */
      logical lex
      character*80 filename
      character*3 ext
c-----|---------------------------------------------------------------|
      call MPI_REDUCE(ppb,ppt,my*mz,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,master,MPI_COMM_WORLD,ierr)
c
      call MPI_REDUCE(uus,uust,(mx1+1)*nspec,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,master,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(vvs,vvst,(mx1+1)*nspec,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,master,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(wws,wwst,(mx1+1)*nspec,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,master,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(pps,ppst,(mx1+1)*nspec,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,master,MPI_COMM_WORLD,ierr)
#ifdef TEMPERATURE
      call MPI_REDUCE(ttb,ttt,my*mz,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,master,MPI_COMM_WORLD,ierr)     
c
      call MPI_REDUCE(tts,ttst,(mx1+1)*nspec,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,master,MPI_COMM_WORLD,ierr)
#endif
c   /* add 2nd moment of omg by sekimoto 7th April 2009*/
      call MPI_REDUCE(omx2,omx2t,my*mz,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,master,MPI_COMM_WORLD,ierr) 
      call MPI_REDUCE(omy2,omy2t,my*mz,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,master,MPI_COMM_WORLD,ierr)     
      call MPI_REDUCE(omz2,omz2t,my*mz,MPI_DOUBLE_PRECISION,
     &     MPI_SUM,master,MPI_COMM_WORLD,ierr)         
      if(myid.eq.master) then
         icf=0
 100     continue
         icf=icf+1
         write(ext,'(i3.3)') icf
#ifdef TEMPERATURE
         call concatb(fuvwpout,'statisticst_',12,ext,3,
     $        '.bin',4,filename)
#else
         call concatb(fuvwpout,'statistics_',11,ext,3,
     $        '.bin',4,filename)
#endif
         inquire(file=filename,exist=lex)
         if(lex)goto 100
c     
         open(iout,file=filename,form="unformatted")
         write(*,*)' writing statistics (bin) to ',filename
         write(iout)timef,time,fnu,alp,aspect,dfloat(nstat),
     $        dfloat(mgalx),dfloat(my),dfloat(mz)
c$$$         fac=f1o1/dfloat(mgalx) !normalize cross-stresses 
         write(iout)
     $        ((um(j,k),j=1,my),k=1,mz),((vm(j,k),j=1,my),k=1,mz),
     $        ((wm(j,k),j=1,my),k=1,mz),((pm(j,k),j=1,my),k=1,mz),
     $        ((uub(j,k),j=1,my),k=1,mz),((uvb(j,k),j=1,my),k=1,mz),
     $        ((uwb(j,k),j=1,my),k=1,mz),((vvb(j,k),j=1,my),k=1,mz),
     $        ((vwb(j,k),j=1,my),k=1,mz),((wwb(j,k),j=1,my),k=1,mz),
     $        ((ppt(j,k),j=1,my),k=1,mz)
         write(iout)dfloat(nspec),
     $        ((uust(i,ispec),i=0,mx1),ispec=1,nspec),
     $        ((vvst(i,ispec),i=0,mx1),ispec=1,nspec),
     $        ((wwst(i,ispec),i=0,mx1),ispec=1,nspec),
     $        ((ppst(i,ispec),i=0,mx1),ispec=1,nspec),
     $        (dfloat(jinco(ispec)),ispec=1,nspec),
     $        (dfloat(kinco(ispec)),ispec=1,nspec)
#ifdef TEMPERATURE
         write(iout)alphat,gravx,gravy,gravz,fkappa,qsource
         write(iout)((tm(j,k),j=1,my),k=1,mz),
     $        ((ttt(j,k),j=1,my),k=1,mz),  
c             ttt <-- ttb by MPI      
     $        ((utb(j,k),j=1,my),k=1,mz),
     $        ((vtb(j,k),j=1,my),k=1,mz),
     $        ((wtb(j,k),j=1,my),k=1,mz)
         write(iout)
     $        ((ttst(i,ispec),i=0,mx1),ispec=1,nspec)
#endif
c   /* add 2nd moment of omg by sekimoto 7th April 2009*/
         write(iout)
     $        ((omx2t(j,k),j=1,my),k=1,mz),
     $        ((omy2t(j,k),j=1,my),k=1,mz),
     $        ((omz2t(j,k),j=1,my),k=1,mz)
         close(iout)
c     
      endif
      return
      end subroutine write_statistics
c-----|---------------------------------------------------------------|
