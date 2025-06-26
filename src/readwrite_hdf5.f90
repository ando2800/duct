!---------------------------------------------------------------------|
#ifdef TEMPERATURE
subroutine read_5fields_fou_hdf5(a,b,c,d,e,filename,rprecision)
#else
subroutine read_4fields_fou_hdf5(a,b,c,d,filename,rprecision)
#endif
  use hdf5
  use ctes
  use running,only:time
  !use timing
  implicit none 
  !
  ! /* reads local 3d fields in binary format by master processor  */
  ! /* the filenames are per-processor and have "myid" as extension*/
  ! /* the header line contains dimensions & physical params       */
  ! /* parallel hdf5 version                                       */
  ! /* rprecision:  'r8'  -  read 8byte data                       */
  ! /*              'r4'   -  read 4byte data and convert to 8     */
  ! /*              [else]   read first four fields  (if defined   */
  ! /*                       TEMPERARTURE, leave unit open)        */
  ! /*                                                             */
  !
  include 'mpif.h'
  complex(nk) a(my,mz,ib:ib+mxp-1), b(my,mz,ib:ib+mxp-1), &
       &      c(my,mz,ib:ib+mxp-1), d(my,mz,ib:ib+mxp-1)
#ifdef TEMPERATURE
  complex(nk) e(my,mz,ib:ib+mxp-1)
#endif
  integer iout,nparamx,ispace,iunit,nparam
  parameter(iout=27)
  character*256 filename,filebase,h5filename
  character*4 ext,ext2
  character*2 rprecision

  ! ------------------------- HDF5 -------------------------------
  integer(8):: bufsize    
  integer:: info
  integer(hid_t):: fid,pid
  integer:: h5err, ierr
  integer(hsize_t), dimension(4):: dims, cdims, totaldims
  integer(hsize_t), dimension(4):: offset
  integer(hid_t):: dset,attr
  integer(hid_t):: dspace,mspace,aspace      
  integer ipl,ii
  complex(8), dimension(:,:,:),allocatable:: resu 
  !---------------------------------------------------------------|
  integer mgalxr,mx1r,myr,mzr
  integer mymax,mzmax,mx1min
  real(nk) alpe,timee,t0
  
  call MPI_INFO_CREATE(info,ierr)
  fid = iunit
  ipl = ie-ib+1
  
  !resu=0.d0 ! bug fixed 2021/04
  t0=0.d0
  ! /* 'filename check' */
  filebase=trim(filename)
  ispace=index(filebase,' ')
  ext=filebase(ispace-4:ispace-1)

  if(ext.eq.'.bin'.or.filebase(ispace-1:ispace-1).eq.'.') then
     write(*,*) 'please remove extension *.bin from filinp'
     stop
  endif
  if(ext(2:4).ne.'.h5') then
     write(*,*) ' ifis=3 is for hdf5'
     stop
  endif
  !
  !   Allocate the temporarily array to read U and W.

  offset = (/0, 0, 0, 0 /)         
  
  if(myid.eq.master) then 
     write(*,*) myid,'reading in fou-space: ',trim(filename)
  end if
  if(rprecision.eq.'r8')then
     !write(*,*) myid, 'Reading from ', trim(filename) ! for debug
     ! Read uvwp
     call h5pcreate_f(H5P_FILE_ACCESS_F,pid,h5err)
     ! call h5pset_fapl_mpiposix_f(pid,commu,.true.,h5err)
     call h5pset_fapl_mpio_f(pid,MPI_COMM_WORLD,info,h5err)
     call h5fopen_f(trim(filename),H5F_ACC_RDONLY_F,fid,h5err,pid)
     if(h5err.ne. 0) then
        write(*,*) "ERROR: Problem openning file", trim(filename)
        call h5eprint_f(h5err,trim(filename))
        stop        
     endif
     call h5pclose_f(pid,h5err)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     if (myid.eq.0) then
        !write(*,*) 'reading header' 
        call readheader_hdf5(fid,mgalxr,mx1r,myr,mzr,alpe,timee)
        
        if(mgalxr.ne.mgalx.or.myr.ne.my.or.mzr.ne.mz)then
           write(*,*)myid,'mismatching dimensions!', mgalxr,myr,mzr
           write(*,*)myid,'forced to interpolate to ', mgalx,my,mz
           !stop
        endif
        if (mgalxr.eq.mgalx.and.myr.eq.my.and.mzr.eq.mz) then
           time = timee
        else
           write(*,*) 'time zero reset, timee=', timee
        end if
       
     end if
     call MPI_BCAST(mgalxr,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(mx1r,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(myr,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(mzr,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(alpe,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
     call MPI_BCAST(time,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)

     mx1min=min(mx1,mx1r)

     allocate(resu(myr,mzr,mxp))
     resu=0.d0

     if (ie.le.mx1r) then
        !write(*,*) myid,'ib, ie, mx1r',ib,ie,mx1r
        dims = (/2, myr, mzr, ipl /) ! local dimension 
     elseif (ib.le.mx1r) then
        !write(*,*) myid,'ib, ie, mx1r',ib,ie,mx1r
        dims = (/2, myr, mzr, mx1r-ib+1 /) 
     else
        dims = (/2, myr, mzr, 0 /) 
     end if

     ! Load the data to the allocated array and close the file
     call h5load_parallel(fid,"u",4,dims,myid,numerop, & 
          &               MPI_COMM_WORLD,info,resu,h5err)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
     if ((my.ne.myr).and.(mz.ne.mzr)) then
        do ii=ib,ie
           call cheb_interp(myr,mzr,resu(1,1,ii-ib+1),my,mz,a(1,1,ii)) 
        end do
     else
        a(1:my,1:mz,ib:ie) = resu(1:my,1:mz,1:ie-ib+1)
     end if

     call h5load_parallel(fid,"v",4,dims,myid,numerop,& 
          &               MPI_COMM_WORLD,info,resu,h5err)
     if ((my.ne.myr).and.(mz.ne.mzr)) then
        do ii=ib,ie
           call cheb_interp(myr,mzr,resu(1,1,ii-ib+1),my,mz,b(1,1,ii)) 
        end do
     else
        b(1:my,1:mz,ib:ie) = resu(1:my,1:mz,1:ie-ib+1)
     end if

     t0=t0-MPI_WTIME()
     call h5load_parallel(fid,"w",4,dims,myid,numerop,& 
          &               MPI_COMM_WORLD,info,resu,h5err)
     if ((my.ne.myr).and.(mz.ne.mzr)) then
        do ii=ib,ie
           call cheb_interp(myr,mzr,resu(1,1,ii-ib+1),my,mz,c(1,1,ii)) 
        end do
     else
        c(1:my,1:mz,ib:ie) = resu(1:my,1:mz,1:ie-ib+1)
     end if

     call h5load_parallel(fid,"p",4,dims,myid,numerop,& 
          &               MPI_COMM_WORLD,info,resu,h5err)
     if ((my.ne.myr).and.(mz.ne.mzr)) then
        do ii=ib,ie
           call cheb_interp(myr,mzr,resu(1,1,ii-ib+1),my,mz,d(1,1,ii)) 
        end do
     else
        d(1:my,1:mz,ib:ie) = resu(1:my,1:mz,1:ie-ib+1)
     end if

#ifdef TEMPERATURE     
     call h5load_parallel(fid,"t",4,dims,myid,numerop, &
          &               MPI_COMM_WORLD,info,resu,h5err)
     if ((my.ne.myr).and.(mz.ne.mzr)) then
        do ii=ib,ie
           write(*,*) myid, 'interpolating in y,z',ii
           call cheb_interp(myr,mzr,resu(1,1,ii-ib+1),my,mz,e(1,1,ii)) 
        end do
     else     
        e(1:my,1:mz,ib:ie) = resu(1:my,1:mz,1:ie-ib+1)
     end if
#endif     
     ! Close the file
     call h5fclose_f(fid,h5err)
  else
     write(*,*)' precision not implemented in read: ',rprecision
     stop
  endif
  return
#ifdef TEMPERATURE
end subroutine read_5fields_fou_hdf5
#else
end subroutine read_4fields_fou_hdf5
#endif
!---------------------------------------------------------------------|
subroutine h5load_parallel(fid,name,ndims,dims,rank, &
     &     size,comm,info,data,ierr)
  !
  ! double precision
  !-------------------------------------------------------------------
  use hdf5
  use ctes, only: myid
  implicit none
  include "mpif.h"
  
  integer(hid_t), intent(in):: fid
  character(len=*), intent(in):: name
  integer, intent(in):: ndims
  integer(hsize_t), dimension(ndims), intent(in):: dims
  integer, intent(in):: rank,size
  integer, intent(in):: comm,info
  real(kind = 8),intent(out):: data
  !
  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace
  integer(hid_t):: plist_id
  integer(hsize_t), dimension(ndims):: start,nooffset,totaldims
  integer, dimension(size):: lastdims
  integer:: mpierr,ierr
  
  integer:: i,lastdim
  real*8 t0


  t0=0.d0
  t0=t0-MPI_WTIME()
  
  start = 0
  nooffset = 0
  totaldims = dims
  
  lastdim = dims(ndims)     ! Don't mess with ints and longs
  
  call MPI_ALLGATHER(lastdim,1,MPI_INTEGER,lastdims,1, &
       &     MPI_INTEGER,comm,mpierr)
  
  totaldims(ndims) = sum(lastdims)
  
  !Open the global dataset and get the global dataspace
  call h5dopen_f(fid,name,dset,ierr)
  call h5dget_space_f(dset,dspace,ierr)
  
  !Create the local dataset
  call h5screate_simple_f(ndims,dims,mspace,ierr)
  call h5sselect_hyperslab_f(mspace,H5S_SELECT_SET_F,nooffset,dims,ierr)
  
  !Select the hyperslab in the global dataset
  start(ndims) = sum(lastdims(1:rank+1))-lastdims(rank+1)
  call h5sselect_hyperslab_f(dspace,H5S_SELECT_SET_F,start,dims,ierr)
  
  !     Create data transfer mode property list
  call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
  
  !Commit the memspace to the disk
  call h5dread_f(dset,H5T_NATIVE_DOUBLE,data,dims, & 
       &     ierr,mspace,dspace,plist_id)
  
  !Close property list                                                       
  call h5pclose_f(plist_id,ierr)
  
  !Close datasets and dataspaces
  call h5sclose_f(mspace,ierr)
  call h5dclose_f(dset,ierr)   
  call h5sclose_f(dspace,ierr)
  
  t0=t0+MPI_WTIME()
  if (myid.eq.0) write(*,*) '  ', trim(name), ' read: ', t0, 'sec.'

end subroutine h5load_parallel
!---------------------------------------------------------------------|
#ifdef TEMPERATURE
subroutine save_5fields_fou_hdf5(a,b,c,d,e,iimag,rprecision)
#else
subroutine save_4fields_fou_hdf5(a,b,c,d,iimag,rprecision)
#endif

  use hdf5
  use ctes

  implicit real*8 (a-h,o-z)
  !  parallel version of hdf5 write 

  !  /* the header line contains dimensions & physical params       */
  !  /*                                                             */
  !  /* rprecision:  'r8'  -  read 8byte data                       */
  !  /*              'r4'   -  read 4byte data and convert to 8     */
  !  /*                                                             */
  !     
  
  include 'mpif.h'
  complex*16 a(my,mz,ib:ie),b(my,mz,ib:ie),c(my,mz,ib:ie),d(my,mz,ib:ie)
#ifdef TEMPERATURE
  complex*16 e(my,mz,ib:ie)
#endif
  integer iout
  parameter(iout=27)
  character*256 filename
  character*4 ext,ext2
  character*2 rprecision
  
  ! ------------------------- HDF5 -------------------------------
  integer(8):: bufsize    
  integer:: info
  integer(hid_t):: fid,pid
  integer:: h5err, ierr
  integer(hsize_t), dimension(4):: dims, cdims, totaldims
  integer(hsize_t), dimension(4):: dim0, start, offset
  integer(hid_t):: dset,attr
  integer(hid_t):: dspace,mspace,aspace      
  integer ipl
  !---------------------------------------------------------------|
  call MPI_INFO_CREATE(info,ierr)
  fid = iunit
  ipl = ie-ib+1;

  dims = (/2, my, mz, ipl /) ! local data size on memory
  cdims = (/2, my, mz ,ipl /) ! local cropped date size (mgalx+2 => mgalx)
  totaldims = (/2, my, mz, mx1+1 /) ! the size of output dataset (u,v,w ...)
  start = (/0, 0, 0, 0/)
  offset = (/0, 0, 0, 0/)      

  write(ext,'(i4.4)') iimag
#ifdef TEMPERATURE
  call concat5b(fuvwpout,'uvwpt_',6,ext,4,'.',1,'h5',2,filename)
#else
  call concat5b(fuvwpout,'uvwp_',5,ext,4,'.',1,'h5',2,filename)
#endif
  if(myid.eq.0)write(*,*)' writing fields to: ',trim(filename)
     
  bufsize = 4*1024*1024

  call h5pcreate_f(H5P_FILE_ACCESS_F,pid,h5err)
  call h5pset_fapl_mpio_f(pid,MPI_COMM_WORLD,info,h5err)
  call h5pset_sieve_buf_size_f(pid, bufsize, h5err)
     
  call h5fcreate_f(filename,H5F_ACC_TRUNC_F,fid,h5err, &
       &        H5P_DEFAULT_F,pid)
  if(h5err.ne.0) then
     write(*,*) "ERROR: Problem writing file ", trim(filename)
     call h5eprint_f(h5err,trim(filename))
     stop
  endif
  call h5pclose_f(pid,h5err)
  
  if(rprecision.eq.'r8')then
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     call h5dump_parallel(fid,"u",4,dims,myid,numerop, &
          &           MPI_COMM_WORLD,info,a,h5err) 
     call h5dump_parallel(fid,"v",4,dims,myid,numerop, &
          &           MPI_COMM_WORLD,info,b,h5err) 
     call h5dump_parallel(fid,"w",4,dims,myid,numerop, &
          &           MPI_COMM_WORLD,info,c,h5err) 
     call h5dump_parallel(fid,"p",4,dims,myid,numerop, &
          &           MPI_COMM_WORLD,info,d,h5err)           
#ifdef TEMPERATURE    
     call h5dump_parallel(fid,"t",4,dims,myid,numerop, &
          &           MPI_COMM_WORLD,info,e,h5err) ! bug fixed 2023/08/23 mitani, d-->e           
#endif
     
  elseif(rprecision.eq.'r4')then
     write(*,*)' single precision is not implemented in save.' 
     stop
  else
     write(*,*)' not implemented in save, rprecision = ',rprecision
     stop
  endif
  call h5fclose_f(fid,h5err) ! here fid is for parallel
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  if (myid.eq.0) then
     call h5fopen_f(filename,H5F_ACC_RDWR_F,fid,h5err)
     call writeheader_hdf5(fid)
     call h5fclose_f(fid,h5err) ! here fid is for sequential
  end if

  return
#ifdef TEMPERATURE
end subroutine save_5fields_fou_hdf5
#else
end subroutine save_4fields_fou_hdf5
#endif
!---------------------------------------------------------------------|
subroutine write_statistics_hdf5
  use ctes
  use running
  use statistics
  use hdf5
  implicit none
  !
  include "mpif.h"
  ! /* tmp array for MPI_REDUCE */
  real(nk) ppt(my,mz)
  real(nk) uust(0:mx1,nspec),vvst(0:mx1,nspec), &
       &   wwst(0:mx1,nspec),ppst(0:mx1,nspec)
#ifdef TEMPERATURE
  real(nk) ttt(my,mz),ttst(0:mx1,nspec)
#endif
  ! /* tmp arrays for MPI_reduce 2nd moment statistics of omg */
  real(nk) omx2t(my,mz),omy2t(my,mz),omz2t(my,mz)
  real(nk) rcvcnt(0:numerop-1)
  integer icf
  !     /* local */
  logical lex
  character*256 filename
  character*4 ext
  ! ------------------------- HDF5 -------------------------------
  integer:: info
  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace
  integer(hid_t):: fid
  integer(hsize_t), dimension(2):: dimstat
  integer(hsize_t), dimension(1):: hdims
  integer:: h5err, ierr
  !---------------------------------------------------------------|
  call MPI_REDUCE(ppb,ppt,my*mz,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,master,MPI_COMM_WORLD,ierr)  
  call MPI_REDUCE(uus,uust,(mx1+1)*nspec,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,master,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(vvs,vvst,(mx1+1)*nspec,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,master,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(wws,wwst,(mx1+1)*nspec,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,master,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(pps,ppst,(mx1+1)*nspec,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,master,MPI_COMM_WORLD,ierr)
#ifdef TEMPERATURE
  call MPI_REDUCE(ttb,ttt,my*mz,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,master,MPI_COMM_WORLD,ierr)     
  
  call MPI_REDUCE(tts,ttst,(mx1+1)*nspec,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,master,MPI_COMM_WORLD,ierr)
#endif
  !   /* add 2nd moment of omg by sekimoto 7th April 2009*/
  call MPI_REDUCE(omx2,omx2t,my*mz,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,master,MPI_COMM_WORLD,ierr) 
  call MPI_REDUCE(omy2,omy2t,my*mz,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,master,MPI_COMM_WORLD,ierr)     
  call MPI_REDUCE(omz2,omz2t,my*mz,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,master,MPI_COMM_WORLD,ierr)         

  if(myid.eq.master) then
     icf=0
100  continue
     icf=icf+1
     write(ext,'(i4.4)') icf
#ifdef TEMPERATURE
     call concatb(fuvwpout,'statisticst_',12,ext,4,'.h5',3,filename)
#else
     call concatb(fuvwpout,'statistics_',11,ext,4,'.h5',3,filename)
#endif
     inquire(file=filename,exist=lex)
     if(lex)goto 100
     
     write(*,*)' writing statistics (bin) to ',trim(filename)
     call h5fcreate_f(trim(filename),H5F_ACC_TRUNC_F,fid,h5err)
     if(h5err.ne.0) then
        write(*,*) h5err, "ERROR: Problem openning ", trim(filename)
        write(*,*) "overwrite the file"
        call h5eprint_f(h5err,trim(filename))
        stop
     endif

     hdims = (/ 1 /)
     call h5write_simple_int(fid,'mgalx',mgalx,1,h5err)
     call h5write_simple_int(fid,'my',my,1,h5err)
     call h5write_simple_int(fid,'mz',mz,1,h5err)
     call h5write_simple_int(fid,'nstat',nstat,1,h5err)
     
     call h5write_simple_double(fid,'timef',timef,1,h5err)  
     call h5write_simple_double(fid,'time',time,1,h5err)
     call h5write_simple_double(fid,'fnu',fnu,1,h5err)    
     call h5write_simple_double(fid,'alp',alp,1,h5err)  
     call h5write_simple_double(fid,'aspect',aspect,1,h5err)    
     
     dimstat=(/my,mz/)
     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'um',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,um,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'vm',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,vm,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'wm',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,wm,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'pm',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,pm,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'uub',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,uub,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)  

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'uvb',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,uvb,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)  

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'uwb',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,uwb,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)  

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'vvb',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,vvb,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)  

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'vwb',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,vwb,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)  

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'wwb',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,wwb,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)  

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'ppb',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,ppt,dimstat,h5err) ! ppb ->ppt
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)  
     ! spectra
     dimstat=(/mx1+1,nspec/)
     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'uust',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,uust,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)  

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'vvst',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,vvst,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)  

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'wwst',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,wwst,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)  

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'ppst',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,ppst,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)  

     call h5write_simple_int(fid,'jinco',jinco,nspec,h5err)
     call h5write_simple_int(fid,'kinco',kinco,nspec,h5err)

#ifdef TEMPERATURE
     call h5write_simple_double(fid,'alphat',alphat,1,h5err)  
     call h5write_simple_double(fid,'gravx',gravx,1,h5err)  
     call h5write_simple_double(fid,'gravy',gravy,1,h5err)  
     call h5write_simple_double(fid,'gravz',gravz,1,h5err)  
     call h5write_simple_double(fid,'fkappa',fkappa,1,h5err)  
     call h5write_simple_double(fid,'qsource',qsource,1,h5err)  

     dimstat=(/my,mz/)
     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'tm',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,tm,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'ttb',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,ttt,dimstat,h5err)! ttb ->ttt
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'utb',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,utb,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'vtb',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,vtb,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'wtb',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,wtb,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)

     dimstat=(/mx1+1,nspec/)
     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'ttst',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,ttst,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)
#endif
     ! /* add 2nd moment of omg by sekimoto 7th April 2009*/
     dimstat=(/my,mz/)
     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'ox2',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,omx2t,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'oy2',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,omy2t,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)

     call h5screate_simple_f(2,dimstat,dspace,h5err)
     call h5dcreate_f(fid,'oz2',H5T_IEEE_F64LE,dspace,dset,h5err)
     call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,omz2t,dimstat,h5err)
     call h5dclose_f(dset,h5err)
     call h5sclose_f(dspace,h5err)

     call h5fclose_f(fid,h5err)

  endif ! master

  return

end subroutine write_statistics_hdf5
! ---------------------------------------------------------------|

! -----------HDF5 all in double precision ---
subroutine h5dump_parallel(fid,name,ndims,dims,rank,size, &
     &     comm,info,data,ierr)
  use hdf5
  use ctes, only: myid,nk
  implicit none
  include "mpif.h"

  ! only for the parallel data with the lastdim being decomposed. 
  
  integer(hid_t), intent(in):: fid
  character(len=*), intent(in):: name
  integer, intent(in):: ndims
  integer(hsize_t), dimension(ndims), intent(in):: dims
  integer, intent(in):: rank,size
  integer, intent(in):: comm,info
  
  real(kind = 8),intent(in):: data
  
  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace
  integer(hid_t):: plist_id
  integer(hsize_t), dimension(ndims):: start,nooffset,totaldims
  integer, dimension(size):: lastdims
  integer:: mpierr,ierr
  
  integer:: i,lastdim
  real(nk) t0

  t0=0.d0
  t0=t0-MPI_WTIME()

  start = 0
  nooffset = 0
  totaldims = dims
  
  lastdim = dims(ndims)     ! Don't mess with ints and longs
  
  call MPI_ALLGATHER(lastdim,1,MPI_INTEGER,lastdims,1,MPI_INTEGER,comm,mpierr)
  
  totaldims(ndims) = sum(lastdims)
  
  !     Create the global dataspace
  call h5screate_simple_f(ndims,totaldims,dspace,ierr)
  !     Create the global dataset
  call h5dcreate_f(fid,name,H5T_IEEE_F64BE,dspace,dset,ierr)
  
  !     Create the local dataset
  call h5screate_simple_f(ndims,dims,mspace,ierr)
  call h5sselect_hyperslab_f(mspace,H5S_SELECT_SET_F,nooffset,dims,ierr)
  
  ! Select the hyperslab in the global dataset
  start(ndims) = sum(lastdims(1:rank+1))-lastdims(rank+1)
  call h5sselect_hyperslab_f(dspace,H5S_SELECT_SET_F,start,dims,ierr)
      
  ! Create data transfer mode property list  
  call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)   
  
  ! Commit the memspace to the disk
  call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,data,dims,ierr,mspace,dspace,plist_id)
    
  ! Close property list                           
  call h5pclose_f(plist_id,ierr)
  
  ! Close datasets and dataspaces
  call h5sclose_f(mspace,ierr)
  call h5dclose_f(dset,ierr)   
  call h5sclose_f(dspace,ierr)

  t0=t0+MPI_WTIME()
  if (myid.eq.0) write(*,*) '  ', trim(name),' written:', t0, 'sec.'
  
end subroutine h5dump_parallel

! -----------HDF5 all in real precision ---
subroutine h5dump_real_parallel(fid,name,ndims,dims,rank,size, &
     &     comm,info,data,ierr)
  use hdf5
  use ctes, only: myid,nk
  implicit none
  include "mpif.h"

  ! only for the parallel data with the lastdim being decomposed. 
  
  integer(hid_t), intent(in):: fid
  character(len=*), intent(in):: name
  integer, intent(in):: ndims
  integer(hsize_t), dimension(ndims), intent(in):: dims
  integer, intent(in):: rank,size
  integer, intent(in):: comm,info
  
  real(kind = 4), intent(in):: data
  
  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace
  integer(hid_t):: plist_id
  integer(hsize_t), dimension(ndims):: start,nooffset,totaldims
  integer, dimension(size):: lastdims
  integer:: mpierr,ierr
  
  integer:: i,lastdim
  real(nk) t0

  t0=0.d0
  t0=t0-MPI_WTIME()

  start = 0
  nooffset = 0
  totaldims = dims
  
  lastdim = dims(ndims)     ! Don't mess with ints and longs
  
  call MPI_ALLGATHER(lastdim,1,MPI_INTEGER,lastdims,1,MPI_INTEGER,comm,mpierr)
  
  totaldims(ndims) = sum(lastdims)
  
  !     Create the global dataspace
  call h5screate_simple_f(ndims,totaldims,dspace,ierr)
  !     Create the global dataset
  call h5dcreate_f(fid,name,H5T_IEEE_F32BE,dspace,dset,ierr)
  
  !     Create the local dataset
  call h5screate_simple_f(ndims,dims,mspace,ierr)
  call h5sselect_hyperslab_f(mspace,H5S_SELECT_SET_F,nooffset,dims,ierr)
  
  ! Select the hyperslab in the global dataset
  start(ndims) = sum(lastdims(1:rank+1))-lastdims(rank+1)
  call h5sselect_hyperslab_f(dspace,H5S_SELECT_SET_F,start,dims,ierr)
      
  ! Create data transfer mode property list  
  call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
  call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)   
  
  ! Commit the memspace to the disk
  call h5dwrite_f(dset,H5T_NATIVE_REAL,data,dims,ierr,mspace,dspace,plist_id)
    
  ! Close property list                           
  call h5pclose_f(plist_id,ierr)
  
  ! Close datasets and dataspaces
  call h5sclose_f(mspace,ierr)
  call h5dclose_f(dset,ierr)   
  call h5sclose_f(dspace,ierr)

  t0=t0+MPI_WTIME()
  if (myid.eq.0) write(*,*) '  ', trim(name),' written in real:', t0, 'sec.'
  
end subroutine h5dump_real_parallel
! ---------------------------------------------------------------|
subroutine readheader_hdf5(fid,mgalxr,mx1r,myr,mzr,alpe,timee)

  use hdf5
  use ctes,only: nk,myid
  implicit none

  integer(hid_t), intent(in):: fid
  
  integer(hsize_t), dimension(1):: hdims
  integer:: h5err

  integer mgalxr,mx1r,myr,mzr,nparamse,dim
  real(nk) alpe,aspecte,timee,fnue,cfle,dtfixede,cmassflowe
#ifdef TEMPERATURE
  real(nk) alphate,gravxe,gravye,gravze,fkappae,qsourcee
#endif

  hdims = (/ 1 /)

  call h5read_simple_int(fid,'nparams',nparamse,1,h5err)
  !
  call h5read_simple_int(fid,'mgalx',mgalxr,1,h5err)  
  call h5read_simple_int(fid,'my',myr,1,h5err)  
  call h5read_simple_int(fid,'mz',mzr,1,h5err)  
  call h5read_simple_int(fid,'mx1',mx1r,1,h5err)  
  if (h5err.ne.0) then
     mx1r=2*(mgalxr/3)/2 -1
  end if
  call h5read_simple_double(fid,'alp',alpe,1,h5err)  
  call h5read_simple_double(fid,'aspect',aspecte,1,h5err)  
  call h5read_simple_double(fid,'time',timee,1,h5err)  
  call h5read_simple_double(fid,'fnu',fnue,1,h5err)  

  call h5read_simple_double(fid,'cfl',cfle,1,h5err)  
  call h5read_simple_double(fid,'dtfixed',dtfixede,1,h5err)  
  call h5read_simple_double(fid,'cmassflow',cmassflowe,1,h5err)  

#ifdef TEMPERATURE
  if (nparamse.eq.7) then
     if(myid.eq.0) write(*,*) 'this restart file does not have temperature fields, stop'
     stop
  end if
  call h5read_simple_double(fid,'alphat',alphate,1,h5err)  
  call h5read_simple_double(fid,'gravx',gravxe,1,h5err)  
  call h5read_simple_double(fid,'gravy',gravye,1,h5err)  
  call h5read_simple_double(fid,'gravz',gravze,1,h5err)  
  call h5read_simple_double(fid,'fkappa',fkappae,1,h5err)  
  call h5read_simple_double(fid,'qsource',qsourcee,1,h5err)  
#endif  
  
end subroutine readheader_hdf5

subroutine writeheader_hdf5(fid)
  use hdf5
  use ctes
  use running,only:time,cfl,dtfixed,slipmx,divmax
  use massflow
  implicit none
  integer(hid_t), intent(in):: fid
  character*4 str
  integer(hsize_t), dimension(3):: dimbc
  integer(hsize_t), dimension(1):: hdims
  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace
  integer:: h5err,nparams

#ifdef TEMPERATURE
  nparams=13
#else
  nparams=7
#endif

  hdims = (/ 1 /)
  call h5write_simple_int(fid,'nparams',nparams,1,h5err)
  call h5write_simple_int(fid,'numerop',numerop,1,h5err)
  call h5write_simple_int(fid,'mgalx',mgalx,1,h5err)  
  call h5write_simple_int(fid,'my',my,1,h5err)  
  call h5write_simple_int(fid,'mz',mz,1,h5err)  
  call h5write_simple_int(fid,'mx1',mx1,1,h5err)  
  
  call h5write_simple_double(fid,'alp',alp,1,h5err)  
  call h5write_simple_double(fid,'aspect',aspect,1,h5err)  
  call h5write_simple_double(fid,'time',time,1,h5err)  
  call h5write_simple_double(fid,'fnu',fnu,1,h5err)  

  call h5write_simple_double(fid,'cfl',cfl,1,h5err)  
  call h5write_simple_double(fid,'dtfixed',dtfixed,1,h5err)  
  call h5write_simple_double(fid,'cmassflow',cmassflow,1,h5err)  
  call h5write_simple_double(fid,'SLIP-VELOCITY',slipmx,1,h5err)    
  call h5write_simple_double(fid,'DIVERGENCE',divmax,1,h5err)    

#ifdef TEMPERATURE
  call h5write_simple_double(fid,'alphat',alphat,1,h5err)  
  call h5write_simple_double(fid,'gravx',gravx,1,h5err)  
  call h5write_simple_double(fid,'gravy',gravy,1,h5err)  
  call h5write_simple_double(fid,'gravz',gravz,1,h5err)  
  call h5write_simple_double(fid,'fkappa',fkappa,1,h5err)  
  call h5write_simple_double(fid,'qsource',qsource,1,h5err)  

  call h5write_simple_char(fid,'boundy1',boundy1,1,h5err)
  call h5write_simple_char(fid,'boundy2',boundy2,1,h5err)
  call h5write_simple_char(fid,'boundz1',boundz1,1,h5err)
  call h5write_simple_char(fid,'boundz2',boundz2,1,h5err)

  dimbc=(/mz,2,mx1+1/)
  call h5screate_simple_f(3,dimbc,dspace,h5err)
  call h5dcreate_f(fid,'bcy',H5T_IEEE_F64LE,dspace,dset,h5err)
  call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,bcy,dimbc,h5err)
  call h5dclose_f(dset,h5err)
  call h5sclose_f(dspace,h5err)
 
  dimbc=(/my,2,mx1+1/)
  call h5screate_simple_f(3,dimbc,dspace,h5err)
  call h5dcreate_f(fid,'bcz',H5T_IEEE_F64LE,dspace,dset,h5err)
  call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,bcz,dimbc,h5err)
  call h5dclose_f(dset,h5err)
  call h5sclose_f(dspace,h5err)

#endif

end subroutine writeheader_hdf5

! -- my simple read write alternative to h5lt library  --
subroutine h5write_simple_double(fid,var,array,dims,h5err)

  use hdf5
  implicit none

  character(len=*) :: var
  integer dims,h5err
  real(8) array(dims)
  ! hdf5 
  integer(hid_t):: fid
  integer(hsize_t), dimension(1):: hdims

  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace

  hdims=dims
  call h5screate_simple_f(1,hdims,dspace,h5err)
  call h5dcreate_f(fid,trim(var),H5T_IEEE_F64LE,dspace,dset,h5err)
  call h5dwrite_f(dset,H5T_NATIVE_DOUBLE,array,hdims,h5err)
  call h5dclose_f(dset,h5err)
  call h5sclose_f(dspace,h5err)

end subroutine h5write_simple_double

subroutine h5write_simple_int(fid,var,iarray,dims,h5err)

  use hdf5
  implicit none

  character(len=*) :: var
  integer dims,h5err
  integer iarray(dims)
  ! hdf5 
  integer(hid_t):: fid
  integer(hsize_t), dimension(1):: hdims

  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace

  hdims=dims
  call h5screate_simple_f(1,hdims,dspace,h5err)
  call h5dcreate_f(fid,trim(var),H5T_STD_I32LE,dspace,dset,h5err)
  call h5dwrite_f(dset,H5T_NATIVE_INTEGER,iarray,hdims,h5err)
  call h5dclose_f(dset,h5err)
  call h5sclose_f(dspace,h5err)

end subroutine h5write_simple_int

subroutine h5write_simple_char(fid,var,str,dims,h5err)

  use hdf5
  implicit none

  character(len=*) :: var
  integer dims,h5err
  character str(dims)
  ! hdf5 
  integer(hid_t):: fid
  integer(hsize_t), dimension(1):: hdims

  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace

  hdims=dims
  call h5screate_simple_f(1,hdims,dspace,h5err)
  call h5dcreate_f(fid,trim(var),H5T_FORTRAN_S1,dspace,dset,h5err)
  call h5dwrite_f(dset,H5T_NATIVE_CHARACTER,str,hdims,h5err)
  call h5dclose_f(dset,h5err)
  call h5sclose_f(dspace,h5err)

end subroutine h5write_simple_char

subroutine h5read_simple_double(fid,var,array,dims,h5err)

  use hdf5
  implicit none

  character(len=*) :: var
  integer dims,h5err
  real(8) array(dims)
  ! hdf5 
  integer(hid_t):: fid
  integer(hsize_t), dimension(1):: hdims

  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace

  hdims=dims

  call h5dopen_f(fid,var,dset,h5err)
  call h5dread_f(dset,H5T_NATIVE_DOUBLE,array,hdims,h5err)
  call h5dclose_f(dset,h5err)

end subroutine h5read_simple_double

subroutine h5read_simple_int(fid,var,array,dims,h5err)

  use hdf5
  implicit none

  character(len=*) :: var
  integer dims,h5err
  integer array(dims)
  ! hdf5 
  integer(hid_t):: fid
  integer(hsize_t), dimension(1):: hdims

  integer(hid_t):: dset
  integer(hid_t):: dspace,mspace

  hdims=dims

  call h5dopen_f(fid,var,dset,h5err)
  call h5dread_f(dset,H5T_NATIVE_INTEGER,array,hdims,h5err)
  call h5dclose_f(dset,h5err)

end subroutine h5read_simple_int
