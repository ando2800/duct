#define DEBUG_Q
#undef DEBUG_Q
#define UVW2Q
!#undef UVW2Q
#define SPECTRA
#undef SPECTRA
#define ENEROUT
!#undef ENEROUT

!--------------------------------------------------------------------
program main
  use hdf5
  use ctes
  use numbers
  use running
  implicit none
  ! /***************************************************************/
  ! /* HDF format is implemented 12/April/2012                     */
  ! /***************************************************************/
  include 'mpif.h'
  include 'fftw3.f'
  !---------------------------------------------------------------------|
  ! /* variables: fourier-space size */
  !
  complex(nk), dimension(:),allocatable:: u,v,w,p,tmpq,p_mod
#ifdef TEMPERATURE
  complex(nk), dimension(:),allocatable:: t
#endif
  !     /* set 2d array */
  !real(nk) uz(my,mz),vz(my,mz),wz(my,mz),pz(my,mz)
#ifdef TEMPERATURE
  !real(nk) tz(my,mz)
#endif
  !real(nk) psy(my,mz) ! for stream function
  ! /* set work arrays */
  ! /* 3d physical-size work arrays: */
  real(nk), dimension(:),allocatable:: tmp1,tmp4,tmp5,tmp6 ! nsfis

  real(nk), dimension(:),allocatable:: adv1,adv2,adv3 ! nsfour
  real(nk), dimension(:),allocatable:: h1,h2,h3,tmp2,tmp3
#ifdef TEMPERATURE
  real(nk),dimension(:),allocatable:: advt
#endif

  ! /* tmp work array for shuffling (note: complex*16)*/
  complex(nk),dimension(:),allocatable:: tmpzxcut
  
  ! /* helmholtz-related work arrays (roughly 2d) */
  real(nk),dimension(:),allocatable:: worke
  
  ! /* global 1d y-or-z work arrays: */
  complex(nk),dimension(:),allocatable:: tmpyz
     
  ! /* global 1d physical-size x work arrays: */
  real(nk),dimension(:),allocatable:: tmpx

  real(nk),dimension(:),allocatable:: tmpyr,tmpzr

  
  ! /* "small" arrays */
  real(nk) stress4(4),energomx4(5),energomxb4(5),stresst4(4)
  integer isave,istop,iimag,ipro,ifis 
  integer i,ist,iiend,ifield,nfields
  integer idump_hdf5
  integer iter,krk,idum
  integer iunit,nprocr
  real(nk) total_byte_l,dum,fmax,fac,total_byte
  character*256 filename,filinp0,h5filename

  parameter(iunit=34)       !some output
  character*3 csave
  parameter(csave='fou')
  character*4 ext

#ifdef SPECTRA
  dimension uus2(nsfour),vvs2(nsfour),wws2(nsfour),pps2(nsfour) 
#endif /* SPECTRA */

  ! --- hdf5 ---
  integer h5err
  ! --- MPI ---
  integer ierr
  !-------------------------------------------------------------------|
  ! /* MPI: initialize   */
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  numerop = nproc
  lastid = nproc-1
  
  ! /* hdf5 initialization */
  !if(myid.eq.master) write(*,*) 'h5open_f'
  call h5open_f(h5err)
  ngpu=0;
  
  if(myid.eq.master) then
     write(*,*) ' '
     write(*,*)'DUCT FLOW POST CODE, MPI Version by Sekimoto'
     call cpp(my,mz)
  end if
  ! /* basic initial set-up */
  call set_numeric_constants()
  !
  call initcr
  !
  allocate( u(nsfou),v(nsfou),w(nsfou),p(nsfou))
  u=0.d0; v=0.d0; w=0.d0; p=0.d0;
#ifdef TEMPERATURE
  allocate(t(nsfou)) ! complex*16
  t=0.d0
#endif /* TMEPERATURE */   

  allocate( tmp1(nsfis))
  allocate( tmp2(nsfour) )
#ifdef UVW2Q   
  allocate( tmp3(nsfour) )
  allocate( tmp4(nsfis),tmp5(nsfis),tmp6(nsfis) )
  tmp1=0.d0; tmp4=0.d0; tmp5=0.d0; tmp6=0.d0;

#endif
  !allocate( h1(nsfour),h2(nsfour),h3(nsfour) )
  allocate( h1(nsfour) )
#ifdef TEMPERATURE
  allocate(advt(nsfour))
  advt=0.d0
#endif
  ! tmp array for shuffling 
  allocate(tmpzxcut(my*mzp*mxp)) ! complex*16
  allocate(worke(nwke))
  allocate(tmpyz(my+mz)) ! complex*16
  allocate(tmpx(mgalx+2))
  allocate(tmpyr(my),tmpzr(mz))
  tmpzxcut=0.d0; worke=0.d0; tmpyz=0.d0; tmpx=0.d0; tmpyr=0.d0; tmpzr=0.d0

  if(myid.eq.master) then
     write(*,*) ' The required memory is '
     total_byte_l= 8d0*(dfloat(nsfou)*2*5+dfloat(nsfis)*4 &
          &                 +dfloat(nsfour)*8)/(2**30)
     write(*,*) total_byte_l,'G byte *',numerop, &
          &              ' = ',total_byte_l*numerop
  endif
  
  ! /* initializing for helmholtz solver*/
  ipro=-1
#ifdef HELM_SHEN
  call helm3d_shen(dum,dum,worke,ipro,'_','_',f1o1,aspect,dum)
  call leg_cheb_mat(worke,worke(my*mz+1),my)
#endif /* HELM_SHEN */
  call helm3d(dum,dum,worke,ipro,'d','_',f1o1,aspect,dum)
  ipro=1

#ifdef TEMPERATURE
  call set_bc_temp
  if(myid.eq.master) then
     write(*,*) 'TEMPERATURE ------------------------------------'
     write(*,*) 'boundary condition:',boundy1,boundy2,boundz1,boundz2
     write(*,*) '------------------------------------------------'
  endif
  ipro=-1
  call helm3d_temp(dum,dum,worke,ipro, &
       &     't',f1o1,aspect,dum,dum,dum,dum)
  ipro=1
#endif
  ! 
  idump_hdf5 = 0
  if(myid.eq.master) then
     write(*,*)'field in fis [1] or fou [2] or hdf5 [3]?'
     read(*,*) ifis
     write(*,*) 'I will read uvwp(t)_****.myid '
     write(*,*) 'set the input file folder as in runparams.dat'
     read(*,'(a)') filinp 
     write(*,*) trim(filinp)
     write(*,*) 'set the output folder as in runparams.dat'
     read(*,'(a)') fuvwpout
     write(*,*) trim(fuvwpout)
     write(*,*) 'index (istart < **** < iiend)'
     write(*,*) 'istart?'
     read(*,*) ist
     write(*,*) 'iiend?'
     read(*,*) iiend
         
     write(*,*) 'I will read ist, iiend: ',ist, iiend !,'in the path:', trim(fuvwpout)
     if (ifis.eq.2) then ! read .myid file (f-p-p)
        if (myid.eq.0) then
           write(*,*) 'ifis=2 set numproc in the file'
           read(*,*) nprocr
        end if
     end if
     if (ifis.lt.3) then
        write(*,*) 'idump_hdf5 (for ifis=3): only write f-c-c data '
        read(*,*) idump_hdf5
     end if

  endif !/* only master */
  !      call MPI_BCAST(fuvwpout,80,MPI_CHARACTER,master,
  !     $     MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ifis,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ist,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iiend,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nprocr,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(idump_hdf5,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)

  call MPI_BCAST(filinp,256,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(fuvwpout,256,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
  !     /* check total disc size */
  if(myid.eq.master) then
     nfields=iiend-ist+1
     total_byte=dfloat(mgalx)*dfloat(my)*dfloat(mz)*4d0 &
          &             *dfloat(iiend-ist+1)/dfloat(2**30)
     write(*,*) total_byte,'G byte *',nfields,' 3D fields' 
     !total_byte=dfloat(my)*dfloat(mz)*4.d0*dfloat(iiend-ist+1) & 
     !     &             /dfloat(2**30)
     !write(*,*) total_byte,'M byte *uvwp(t),psy(6) 2D fields'
  endif

  !----  start file loop  ------
  filinp0=trim(filinp)
  do iimag=ist,iiend
     ! /* set read file name */
     if (myid.eq.0) then 
        write(*,*) ' '
        write(*,*) ' - - - - ', iimag, '/', iiend, ' - - - -'
     endif
     write(ext,'(i4.4)') iimag

     filename=trim(filinp0)//ext(1:4)

     if (ifis.eq.1) then  ! read .bin file (p-p-p)
        filinp=filename
        call read_fft_4fields_fis(u,v,w,p,tmp1,tmp4,tmpx,'r8','_')
#ifdef TEMPERATURE
        call read_fft_4fields_fis(t,dum,dum,dum,tmp1,tmp4,tmpx,'r8','t')
#endif /* TEMPERATURE */ 
     elseif (ifis.eq.2) then ! read .myid file (f-p-p)
        !write(*,*) myid,'ifis=2 set numproc in the file',nprocr
#ifdef TEMPERATURE
        call read_5fields_fou_apart(u,v,w,p,t,filename,nprocr,'r8')
#else /* TEMPERATURE */
        call read_4fields_fou_apart(u,v,w,p,filename,nprocr,'r8')
#endif /* TEMPERATURE */                
     elseif (ifis.eq.3) then ! hdf5 f-p-p read
        h5filename=trim(filename)//'.h5'
        if (myid.eq.master) write(*,*) '  ifis=3: read in hdf5 format'
#ifdef TEMPERATURE
        call read_5fields_fou_hdf5(u,v,w,p,t,h5filename,'r8')
#else
        call read_4fields_fou_hdf5(u,v,w,p,h5filename,'r8')
#endif
     else
        write(*,*) ' set ifis <= 3, stop'
        stop
     end if

     if ((ifis.lt.3).and.(idump_hdf5.eq.1)) then
        if(myid.eq.master) write(*,*) 'writing f-c-c data in h5 format ...'
#ifdef TEMPERATURE 
        call save_5fields_fou_hdf5(u,v,w,p,t,iimag,'r8')
#else
        call save_4fields_fou_hdf5(u,v,w,p,iimag,'r8')
#endif      
        cycle ! skip dumping u,v,w etc, and go to next file
     end if

#ifdef SPECTRA
     if(iimag.eq.ist)then
        ! /* initialize statistics */
        if(myid.eq.master) write(*,*) 'initialize uu... spectrum'
        call spectrum(dum,dum,dum,dum,uus2,vvs2,wws2,pps2,nstat, &
             &        0)
     endif
#endif /* SPECTRA */
     !write(*,*) myid,'div3d'
     !     /* some field status output to txt file */
     call div3d(u,v,w,h1,tmp2,tmpyz,f1o1)
     call max_3din(h1,fmax)
     if(myid.eq.0)then
        write(*,*)'DIVERGENCE: ',fmax
     endif
   
     call concat4b(fuvwpout,'fields_',7,ext,4,'.ur.h5',6,filename)
     call save_hdf5_fft_3dfield_fou(u,tmp1,tmpzxcut,tmpx,&
          &        filename,'r4','ur   ',iunit)
     call concat4b(fuvwpout,'fields_',7,ext,4,'.vr.h5',6,filename)
     call save_hdf5_fft_3dfield_fou(v,tmp1,tmpzxcut,tmpx,&
          &        filename,'r4','vr   ',iunit)
     call concat4b(fuvwpout,'fields_',7,ext,4,'.wr.h5',6,filename)
     call save_hdf5_fft_3dfield_fou(w,tmp1,tmpzxcut,tmpx,&
          &        filename,'r4','wr   ',iunit)
#ifdef TEMPERATURE
     call concat4b(fuvwpout,'fields_',7,ext,4,'.tr.h5',6,filename)
     call save_hdf5_fft_3dfield_fou(t,tmp1,tmpzxcut,tmpx,&
          &        filename,'r4','tr   ',iunit)
#endif

     allocate(tmpq(nsfou))
#ifdef UVW2Q   
     ! /* compute laplacian of p from u,v,w = qcriterion, */
     ! /*  but noisy on wall*/
     ! /* (0) initializing adv1,adv2,adv3 */
     allocate( adv1(nsfour),adv2(nsfour),adv3(nsfour) )
     adv1=0.d0; adv2=0.d0; adv3=0.d0; ! initialized before advection3d
     ipro=1
     !call dcopy(nsfour,f0o1,0,adv1,1)
     !call dcopy(nsfour,f0o1,0,adv2,1)
     !call dcopy(nsfour,f0o1,0,adv3,1)      
     ! /* (1) compute advection terms, NL(u) ->adv1,adv2,adv3  */
     krk=0
     iter=0
     if (myid.eq.0) write(*,*) 'calculating Qa'
     call advection3d(u,v,w,adv1,adv2,adv3, &
          &           tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmpzxcut,tmpx, &
          &           krk,iter)
     ! /* (2) compute r.h.s.:    tmp1 = -\nabla \dot [NL(u)] = qcri  */
     call div3d(adv1,adv2,adv3,tmpq,tmp2,tmpyz,-f1o1)
     fac = f1o1/f2o1 ! bug fixed 2011/05/03 sekimoto
     ! -lap.(p) = -2 Q
     call dscal2(nsfour,tmpq,1,fac,tmpq,1) ! the same with (1/2)*lap.p
     deallocate(adv1,adv2,adv3)
#endif/* UVW2Q */ 
     ! /* (3) solve helmholtz equation  \nabla^2 (P) = ...    */
     ! /*     to compure pressure from u,v,w.                 */

#ifdef MAKE_PRESSURE
     
     allocate( p_mod(nsfou))
     p_mod=0.d0
#ifdef HELM_SHEN
     call helm3d_shen(tmpq,p_mod,worke,ipro,'n','_',f1o1,aspect,f0o1)
#else /* HELM_SHEN */
     call helm3d(tmpq,p_mod,worke,ipro,'n','_',f1o1,aspect,f0o1,dum)
#endif /* HELM_SHEN */
     deallocate(p_mod)
#endif /* MAKE_PRESSURE */
     
#ifndef UVW2Q
     fac = f1o1/f2o1 ! bug fixed 2011/05/03 sekimoto
     ! -lap.(p) = -2 Q
     call laplacian3d(p,tmpq,tmp2,tmpyz,fac,'n')
#endif      
     call repair_q_boundary(tmpq)
     
     !#ifdef DEBUG_Q
     !c$$     /* compute smoother (less noisy) qcriterion */
     !c$$        but too much cost
     !c$$         call second_invariant_fou(dfm11,dfm12,dfm13,
     !c$$     $        dfm21,dfm22,dfm23,
     !c$$     $        dfm31,dfm32,dfm33,
     !c$$     $        qcrit,tmp1,tmp4,tmpx)         
     !#endif
     
     call concat4b(fuvwpout,'fields_',7,ext,4,'.Qa.h5',6,filename)
     call save_hdf5_fft_3dfield_fou(tmpq,tmp1,tmpzxcut,tmpx,&
          &        filename,'r4','Qa   ',iunit)
     deallocate(tmpq)
#ifdef SPECTRA
     !     /* compute spectrum statistics */
     if(myid.eq.master) write(*,*) 'store spec'
     call spectrum(u,v,w,p,uus2,vvs2,wws2,pps2,nstat, &
          &        1)
#endif /* SPECTRA */
     
#ifdef SLONG
     !  /* time and x-dependent indicators for analysis of longest structures */
     if (myid.eq.master)then 
        
        ! /* save dudns and dudn (only master) */
        call concat4b(fuvwpout,'dudn_int_',9,ext,4,'.bin',4,filename)
        call save_wallstress4x(filename,stress4x,dudnst,myid,'r4')
        
        ! /* save omg^2 in sectors for indicator (only master) */
        call concat4b(fuvwpout,'energomx4x_',11,ext,4,'.bin',4,filename)
        call save_energ4x(filename,energomx4x,myid,'r4')
        
        ! /* save eneruprim (only master) */
        ! /* corresponds: E_streak */
        call concat4b(fuvwpout,'energuprim4x_',13,ext,4,'.bin',4,filename)
        call save_energ4x(filename,energuprim4x,myid,'r4')
        
        ! /* save energv (only master) */
        call concat4b(fuvwpout,'energv4x_',9,ext,4,'.bin',4,filename)
        call save_energ4x(filename,energv4x,myid,'r4')
        
        ! /* save energw (only master) */
        call concat4b(fuvwpout,'energw4x_',9,ext,4,'.bin',4,filename)
        call save_energ4x(filename,energw4x,myid,'r4')
        
        ! /* save energp (only master) */
        call concat4b(fuvwpout,'energp4x_',9,ext,4,'.bin',4,filename)
        call save_energ4x(filename,energp4x,myid,'r4')
        
     endif
#endif /* SLONG */

     !     /* compute vorticity */
     call vorticity(u,v,w,h1,tmp2,1) !compute omega_x->h1
     call concat4b(fuvwpout,'fields_',7,ext,4,'.ox.h5',6,filename)
     call save_hdf5_fft_3dfield_fou(h1,tmp1,tmpzxcut,tmpx,&
          &        filename,'r4','ox   ',iunit)

     call vorticity(u,v,w,h1,tmp2,2) !compute omega_y->h1
     call concat4b(fuvwpout,'fields_',7,ext,4,'.oy.h5',6,filename)
     call save_hdf5_fft_3dfield_fou(h1,tmp1,tmpzxcut,tmpx,&
          &        filename,'r4','oy   ',iunit)

     call vorticity(u,v,w,h1,tmp2,3) !compute omega_z->h1     
     call concat4b(fuvwpout,'fields_',7,ext,4,'.oz.h5',6,filename)
     call save_hdf5_fft_3dfield_fou(h1,tmp1,tmpzxcut,tmpx,&
          &        filename,'r4','oz   ',iunit)
     
  enddo
  !c----  end of file loop  ------
  
#ifdef SPECTRA
  !     /* save spectrum statistics */
  if(myid.eq.master) write(*,*) 'save spec',time
  call save_spectra(uus2,vvs2,wws2,pps2,nstat,&
       &     time,'to','r4')
#endif /* SPECTRA */
  
111 format('ENER: ',40(e15.8,1x))

999 if(myid.eq.master) write(*,*) 'close hdf5 file units'
  !     /* close hdf file units */
  call h5close_f(h5err)
  
!    fftw destroy
!  call dfftw_destroy_plan(plan_fft_forward)
!  call dfftw_destroy_plan(plan_fft_backward)
!  call dfftw_destroy_plan(plan_tcheb_y)
!  call dfftw_destroy_plan(plan_tcheb_z)

1000 call MPI_FINALIZE(ierr) 
  
  stop
end program main
!----|---------------------------------------------------------      
subroutine second_invariant_fou(d11,d12,d13,&
     &                          d21,d22,d23,&
     &                          d31,d32,d33,q,tmpfis,tmpfis2,tmpx)
  use ctes
  implicit real*8(a-h,o-z)
  !$$$       implemented by sekimoto (2006/April 19)      
  !     /***************************************************************/
  !     /* input      p-p-p                                            */
  !     /* compute second invalient of velocity gradient tensor        */
  !     /*                                  in fourier space           */
  !     /*                                                         .   */
  !     /*                                                             */
  !     /* output  Qcriterion : second invaliant (p-p-p :x-y-z)        */
  !$$$  /* tmp(nsfou)          temp array in fourier space             */
  !     /* dfm**               tensor D component in fourier space     */
  !     /***************************************************************/

  dimension d11(nsfour),d12(nsfour),d13(nsfour)
  dimension d21(nsfour),d22(nsfour),d23(nsfour)
  dimension d31(nsfour),d32(nsfour),d33(nsfour)
  dimension tmpfis(my,mz,mgalx+2) !fis-size (tmp buffer)
  dimension tmpfis2(my,mz,mgalx+2) !fis-size (tmp buffer)
  dimension tmpx(mgalx+2)   !1d fis--size x work array
  dimension q(my,mz,mgalx+2)
  ntotal = (mgalx+2)*mz*my  !2007/4/12

  write(*,*) 'not parallelized'
  stop

  call fourx(d11,tmpfis,tmpx,+1) !fou->fis
  do i=1,mgalx
     do k=1,mz
        do j=1,my
           q(j,k,i) = tmpfis(j,k,i)*tmpfis(j,k,i) 
        end do
     enddo
  enddo
  
  call fourx(d22,tmpfis,tmpx,+1) !fou->fis
  do i=1,mgalx
     do k=1,mz
        do j=1,my
           q(j,k,i) = q(j,k,i) + tmpfis(j,k,i)*tmpfis(j,k,i)
        end do
     enddo
  enddo
  
  call fourx(d33,tmpfis,tmpx,+1) !fou->fis
  do i=1,mgalx
     do k=1,mz
        do j=1,my
           q(j,k,i) = q(j,k,i) + tmpfis(j,k,i)*tmpfis(j,k,i)
        end do
     enddo
  enddo
  
  call fourx(d12,tmpfis,tmpx,+1) !fou->fis
  call fourx(d21,tmpfis2,tmpx,+1) !fou->fis
  do i=1,mgalx
     do k=1,mz
        do j=1,my
           q(j,k,i) = -0.5d0*q(j,k,i) - tmpfis(j,k,i)*tmpfis2(j,k,i)
        end do
     end do
  end do
  
  call fourx(d13,tmpfis,tmpx,+1) !fou->fis
  call fourx(d31,tmpfis2,tmpx,+1) !fou->fis
  do i=1,mgalx
     do k=1,mz
        do j=1,my
           q(j,k,i) = q(j,k,i) - tmpfis(j,k,i)*tmpfis2(j,k,i)
        end do
     enddo
  enddo
  
  call fourx(d23,tmpfis,tmpx,+1) !fou->fis
  call fourx(d32,tmpfis2,tmpx,+1) !fou->fis
  do i=1,mgalx
     do k=1,mz
        do j=1,my
           q(j,k,i) = q(j,k,i) - tmpfis(j,k,i)*tmpfis2(j,k,i)
        end do
     enddo
  enddo
  return
end subroutine second_invariant_fou

!---------------------------------------------------------------------|
subroutine repair_q_boundary(qcrit)
  use ctes,only:nk,my,mz,ib,ie
  use numbers
  implicit none
  complex(nk) qcrit(my,mz,ib:ie)
  integer i,k,j
  !  /* parallelized                                               */
  !  /* implemented, 2011/05/02 sekimoto                           */
  !  /* repair_boundary of Q criterion: should be zero             */
  !  /* set boundary to zero */
  !  /* y=-1 & y=+1 planes */
  do i=ib,ie
     do k=1,mz
        do j=1,my,my-1
           qcrit(j,k,i)=zero
        enddo
     enddo
  enddo

  ! /* z=-a & z=+a planes */
  do i=ib,ie
     do k=1,mz,mz-1
        do j=1,my
           qcrit(j,k,i)=zero
        enddo
     enddo
  enddo
  
  return
end subroutine repair_q_boundary
!---------------------------------------------------------------------|
#ifdef SPECTRA
subroutine spectrum(u,v,w,p,uus2,vvs2,wws2,pps2,nstat,iopt)
  use ctes
  implicit real*8(a-h,o-z)
  ! /* parallelized                                               */
  ! /* accumulate streamwise spectrum of uu, vv, ww               */
  ! /* iopt=0    : set null value to uups,...                     */
  ! /* iopt=1    : add to common data(3D)                         */
  ! /*                                                            */
  include "ctes3D" 
  complex*16 u(my,mz,ib:ie),v(my,mz,ib:ie),w(my,mz,ib:ie),&
       p(my,mz,ib:ie)
  
  dimension uus2(my,mz,ib:ie),vvs2(my,mz,ib:ie),&
       &    wws2(my,mz,ib:ie),pps2(my,mz,ib:ie) 

  if(iopt.eq.0) then
     !     /* set null */
     nstat=0
     timef=time
     
     call dcopy((ie-ib+1)*my*mz,f0o1,0,uus2,1)
     call dcopy((ie-ib+1)*my*mz,f0o1,0,vvs2,1)
     call dcopy((ie-ib+1)*my*mz,f0o1,0,wws2,1)
     call dcopy((ie-ib+1)*my*mz,f0o1,0,pps2,1)
  else
     nstat=nstat+1
     do i=max(ib,1),ie             !non-zero fourier modes
        do k=1,mz
           do j=1,my
              uus2(j,k,i)=uus2(j,k,i)+&
                   &                 (u(j,k,i)*dconjg(u(j,k,i)))*f2o1
              vvs2(j,k,i)=vvs2(j,k,i)+&
                   &                 (v(j,k,i)*dconjg(v(j,k,i)))*f2o1
              wws2(j,k,i)=wws2(j,k,i)+&
                   &                 (w(j,k,i)*dconjg(w(j,k,i)))*f2o1
              pps2(j,k,i)=pps2(j,k,i)+&
                   &                 (p(j,k,i)*dconjg(p(j,k,i)))*f2o1
           enddo
        enddo
     enddo
     if (ib.eq.0) then
        i=0                 !0-fourier mode
        do k=1,mz
           do j=1,my
              uus2(j,k,i)=uus2(j,k,i)+(u(j,k,i)*dconjg(u(j,k,i)))
              vvs2(j,k,i)=vvs2(j,k,i)+(v(j,k,i)*dconjg(v(j,k,i)))
              wws2(j,k,i)=wws2(j,k,i)+(w(j,k,i)*dconjg(w(j,k,i)))
              pps2(j,k,i)=pps2(j,k,i)+(p(j,k,i)*dconjg(p(j,k,i)))
           enddo
        enddo
     endif
  endif

  return
end subroutine spectrum
#endif /* SPECTRA */
!---------------------------------------------------------------------|

!---------------------------------------------------------------------|
!       I/O utility to make HDF5 data 
!---------------------------------------------------------------------|
subroutine save_hdf5_fft_3dfield_fou(a,tmpfis,tmpzxcut,tmpx,&
     filehdf5,rprecision,cflag,iunit)
  use hdf5
  use ctes
  implicit none 
  !  parallelized to kb:ke
  ! /* writes a 3d field in HDF5 format                            */
  ! /* note: there is no header line and append to file(iunit)     */
  ! /*       for time step data into a file                        */
  ! /* input: Z-cut data                                           */
  ! /*                                                             */
  ! /* rprecision:  'r8'  -  save 8byte data                       */
  ! /*              'r4'   -  convert to 4 byte data and           */
  ! /*                        save 4byte data                      */
  ! /*                                                             */

  include "mpif.h"
  real(nk) a(nsfour)
  real(nk) tmpfis(nsfis) !fis-size (tmp buffer)
  real(4), dimension(:,:,:), allocatable::res
  !real(4), dimension(:), allocatable::res
  real(nk), dimension(:), allocatable::atmp
  complex(nk) tmpzxcut(my*mzp*mxp) !z-cut & x-cut pencil work array
  real(nk) tmpx(mgalx+2)   !1d fis--size x work array
  character*256 filebase, filehdf5
  character*4 ext,ext2
  character*2 rprecision
  character*5 cflag
  integer h5err,kpl,leng,istat,iproc,iprocr,ierr,iwrite_in_parallel,iunit
  
  ! --- hdf5 ---
  integer(8):: bufsize    
  integer:: info_r
  integer(hid_t):: fidr, pidr, info_id
  integer(hsize_t), dimension(3):: dims, cdims, totaldims
  integer(hsize_t), dimension(3):: dim0, start, offset
  integer(hsize_t), dimension(1):: hdims
  !character(len=*):: name
  integer(hid_t):: dset,attr
  integer(hid_t):: dspace,mspace,aspace
  integer ipl
  !-----|---------------------------------------------------------------|
  call MPI_INFO_CREATE(info_r,ierr)
  fidr = iunit 
  !
  kpl=ke-kb+1;
  !
  dim0 = (/1,0,0/)
  !dims = (/ my,mzp,mgalx+2 /) ! BUG local data size on memory for tmpfis
  !dims = (/ my,mzp,mgalx /) ! local data size on memory for res
  !cdims = (/my,kpl,mgalx /) ! local cropped date size (mgalx+2 => mgalx)
  !totaldims = (/my,mz,mgalx /) ! the size of output dataset (u,v,w ...)
  !
  dims = (/ mgalx,my,mzp /) ! local data size on memory for res !mzp=maxval(kend-kbeg+1)
  cdims = (/ mgalx,my,kpl /) ! local cropped date size (mgalx+2 => mgalx)
  totaldims = (/ mgalx,my,mz /) ! the size of output dataset (u,v,w ...)
  !
  start = (/ 0, 0, 0/)
  offset = (/ 0, 0, 0/)
  !-----|---------------------------------------------------------------|
  if(kb.eq.1)  write(*,*) ' save_hdf5_fft_3dfield_fou: writing ',cflag, &
       &     ' in the precision: ',rprecision
  ! /* shuffle arrays from x-cut to z-cut */
  ! write(*,*) myid,': shuffling xcut -> zcut: a,b,c,d'
  allocate(atmp(nsfour))
  atmp = 0.d0
  ! save a: u,v,w
  atmp = a
  call chx2z(atmp,tmpfis,tmpzxcut)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  tmpfis=0.d0; tmpx=0.d0
  call fourx(atmp,tmpfis,tmpx,+1) !a: fou -> tmpfis: fis
  deallocate(atmp)

  allocate(res(mgalx,my,kpl)) ! real*4 array
  res = 0.d0
  call crop_real4(tmpfis,res,kb,ke)
  
  !write(*,*) myid, 'data croped in real4 array'

  iwrite_in_parallel=1
  if (iwrite_in_parallel.eq.1) then

     bufsize = 4*1024*1024

     call h5pcreate_f(H5P_FILE_ACCESS_F,pidr,h5err)
     call h5pset_fapl_mpio_f(pidr,MPI_COMM_WORLD,info_r,h5err)
     call h5pset_sieve_buf_size_f(pidr, bufsize, h5err)
     
     call h5fcreate_f(filehdf5,H5F_ACC_TRUNC_F,fidr,h5err, &
          &        H5P_DEFAULT_F,pidr)
     if(h5err.ne.0) then
        write(*,*) "ERROR: Problem writing file ", trim(filehdf5)
        call h5eprint_f(h5err,trim(filehdf5))
        stop
     endif
     call h5pclose_f(pidr,h5err)

     !write(*,*) myid,' write in parallel'
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     call h5dump_real_parallel(fidr,cflag,3,cdims,myid,numerop, &
          &     MPI_COMM_WORLD,info_r,res,h5err)
     call h5fclose_f(fidr,h5err) ! here fidr is for parallel
     
     if (myid.eq.0) then
        call h5fopen_f(filehdf5,H5F_ACC_RDWR_F,fidr,h5err)
        call writeheader_hdf5(fidr)
        call h5fclose_f(fidr,h5err) ! here fidr is for sequential
     end if
     
  else

     if (myid.eq.master) then

        write(ext2,'(i4.4)') myid
        
        write(*,*) myid,':writing :',trim(filehdf5)
        
        call h5fcreate_f(filehdf5,H5F_ACC_TRUNC_F,fidr,h5err)
        if(h5err.ne.0) then
           write(*,*) h5err, "ERROR: Problem openning ", trim(filehdf5)
           write(*,*) "overwrite the file"
           call h5eprint_f(h5err,trim(filehdf5))
           stop
        endif
        
        hdims = (/ 1 /)
        call h5write_simple_int(fidr,'mgalx',mgalx,1,h5err)
        call h5write_simple_int(fidr,'my',my,1,h5err)
        call h5write_simple_int(fidr,'mz',mz,1,h5err)
        
        ! main data
        call h5screate_simple_f(3,totaldims,dspace,h5err)
        call h5dcreate_f(fidr,cflag,H5T_IEEE_F32BE,dspace,dset,h5err)
        
        start(3)=0
        do iproc=0,numerop-1
           ! master receive from slaves
           if (iproc.ne.0) then
              deallocate(res)
              allocate(res(mgalx,my,kbeg(iproc):kend(iproc)))
              res=0.d0
              !write(*,*)  'master receives from ',iproc
              leng=mgalx*my*(kend(iproc)-kbeg(iproc)+1)
              iprocr = iproc
              call MPI_RECV(res,leng,MPI_REAL,iprocr, &
                   &    iprocr,MPI_COMM_WORLD,istat,ierr) 
              !write(*,*)  'master received'
           endif
           ! index should be 2 --> 3 because of updated change
           cdims(3) = kend(iproc)-kbeg(iproc)+1
           
           ! this is important line to disturb too much optimisation
           ! crop and covert to real*4 the data on memory
           !  call crop_real4(tmpfis,res,kbeg(id),kend(id))
           ! create the local data set
           call h5screate_simple_f(3,dims,mspace,h5err)
           call h5sselect_hyperslab_f(mspace,H5S_SELECT_SET_F, & 
                offset,cdims,h5err)
           
           start(3) = kbeg(iproc)-1 ! 0-start index
           !write(*,*) 'select start=',start
           call h5dget_space_f(dset,dspace,h5err)
           call h5sselect_hyperslab_f(dspace,H5S_SELECT_SET_F, & 
                start,cdims,h5err)
           !      
           !write(*,*) ' h5dwrite '
           ! Commit the memspace to the disk
           call h5dwrite_f(dset,H5T_NATIVE_REAL,res,cdims,h5err, & 
                mspace,dspace,H5P_DEFAULT_F)
           ! do not flush !
           !write(*,*) ' h5fflush'
           !call h5fflush_f(fidr,H5F_SCOPE_LOCAL_F,h5err) 
           call h5sclose_f(mspace,h5err)
           call h5sclose_f(dspace,h5err) 
           !
        enddo
        !write(*,*) ' h5sclose'
        !call h5sclose_f(mspace,h5err)        
        !Close datasets and global dataspace
        !call h5sclose_f(dspace,h5err) 
        !write(*,*) ' h5dclose'
        call h5dclose_f(dset,h5err)   
        call h5fclose_f(fidr,h5err)
        !
        !call writeheader_hdf5(filehdf5,cflag) ! fopen and fclose are inside
        
     else
        !  /* all slabes send to master */
        leng =  my*(ke-kb+1)*(mgalx)     
        call MPI_SEND(res,leng,MPI_REAL,master, &
             &        myid,MPI_COMM_WORLD,ierr)
     endif !/* non-master */
  end if ! /* i write un parallel */

  deallocate(res)

  !call chz2x(a,tmpfis,tmpzxcut) ! repair a  
  return
end subroutine save_hdf5_fft_3dfield_fou
!----|----------------------------------------------------------------|
subroutine slide(field,islide)
  use ctes,only:my,mz,mgalx
  implicit real*8(a-h,o-z)

  real*8 tmp(my,mz,mgalx+2)
  dimension field(my,mz,mgalx+2)
  do i=1,mgalx
     idash = mod(i+islide,mgalx)
     if(idash == 0) idash = mgalx;
     do k=1,mz 
        do j=1, my
           tmp(j,k,i) = field(j,k,idash)
        end do
     end do
  end do
  do i=1,mgalx
     do k=1,mz
        do j=1,my
           field(j,k,i) = tmp(j,k,i)
        end do
     end do
  end do
  
  return
end subroutine slide
!----|----------------------------------------------------------------|
subroutine save_zeromode(uz,vz,wz,pz,filename,rprecision)
  use ctes
  implicit real*8(a-h,o-z)
  !     /*                                                           */
  !     /*  store zeromode of a to uvwp_zero_****.bin                 */
  
  parameter(iout=34)       !some output
  character*256 filename
  character*2 rprecision
  dimension uz(my,mz),vz(my,mz),wz(my,mz),pz(my,mz) 
  
  open(iout,file=filename,form="unformatted")
  if(rprecision.eq.'r8')then
     write(iout)dfloat(mgalx),dfloat(my),dfloat(mz), &
          &     alp,aspect,time,fnu
     write(iout)((uz(j,k),j=1,my),k=1,mz)
     write(iout)((vz(j,k),j=1,my),k=1,mz)
     write(iout)((wz(j,k),j=1,my),k=1,mz)
     write(iout)((pz(j,k),j=1,my),k=1,mz)
  elseif(rprecision.eq.'r4')then
     write(iout)real(mgalx),real(my),real(mz),real(alp), &
          &     real(aspect),real(time),real(fnu)
     write(iout) ((real((uz(j,k))),j=1,my),k=1,mz) 
     write(iout) ((real((vz(j,k))),j=1,my),k=1,mz) 
     write(iout) ((real((wz(j,k))),j=1,my),k=1,mz) 
     write(iout) ((real((pz(j,k))),j=1,my),k=1,mz) 
  else
     write(*,*)'precision not implemented in save: ',rprecision
     stop
  endif
  close(iout)
  return
end subroutine save_zeromode
!----|----------------------------------------------------------------|
#ifdef TEMPERATURE
subroutine save_zeromode_t(uz,vz,wz,pz,tz,filename,rprecision)
  use ctes
  implicit real*8(a-h,o-z)
  !  /*                                                           */
  !  /*  store zeromode of a to uvwpt_zero_****.bin               */

  parameter(iout=34)       !some output
  character*256 filename
  character*2 rprecision
  dimension uz(my,mz),vz(my,mz),wz(my,mz),pz(my,mz),tz(my,mz) 
  !-----|---------------------------------------------------------------|
  open(iout,file=filename,form="unformatted")
  if(rprecision.eq.'r8')then
     write(iout)dfloat(mgalx),dfloat(my),dfloat(mz), &
          &         alp,aspect,time,fnu,alphat,gravx,gravy,gravz, &
          &        fkappa,qsource
     write(iout)((uz(j,k),j=1,my),k=1,mz)
     write(iout)((vz(j,k),j=1,my),k=1,mz)
     write(iout)((wz(j,k),j=1,my),k=1,mz)
     write(iout)((pz(j,k),j=1,my),k=1,mz)
     write(iout)((tz(j,k),j=1,my),k=1,mz)
  elseif(rprecision.eq.'r4')then
     write(iout)real(mgalx),real(my),real(mz),real(alp), &
          &     real(aspect),real(time),real(fnu),real(alphat), &
          &     real(gravx),real(gravy),real(gravz),real(fkappa), &
          &     real(qsource)
     write(iout) ((real((uz(j,k))),j=1,my),k=1,mz) 
     write(iout) ((real((vz(j,k))),j=1,my),k=1,mz) 
     write(iout) ((real((wz(j,k))),j=1,my),k=1,mz) 
     write(iout) ((real((pz(j,k))),j=1,my),k=1,mz) 
     write(iout) ((real((tz(j,k))),j=1,my),k=1,mz) 
  else
     write(*,*)'precision not implemented in save: ',rprecision
     stop
  endif
  close(iout)
  return
end subroutine save_zeromode_t
#endif /* TEMPETATURE */
!----|----------------------------------------------------------------|
#ifdef SLONG
subroutine save_wallstress4x(filename,stress4x,dudnst, &
     &                       rprecision)
  use ctes
  implicit real*8(a-h,o-z)
  ! /*                                                           */
  ! /*  store zeromode of a to uvwpt_zero_****.bin               */

  parameter(iout=34)       !some output
  character*256 filename
  character*2 rprecision
  !      
  dimension stress4x(mgalx+2,4)
  complex*16 dudnst(0:mx1,4) ! for master reduction
  !-----|---------------------------------------------------------------|
  open(iout,file=filename,form="unformatted")
  if(rprecision.eq.'r8')then
     write(*,*)'precision not implemented in save: ',rprecision
     stop
  elseif(rprecision.eq.'r4')then
     write(iout) ((real((stress4x(i,iw))),i=1,mgalx),iw=1,4) 
     write(iout) ((cmplx((dudnst(i,iw))),i=0,mx1),iw=1,4) 
  else
     write(*,*)'precision not implemented in save: ',rprecision
     stop
  endif
  close(iout)
  return
end subroutine save_wallstress4x
!-----|---------------------------------------------------------------|
subroutine save_energ4x(filename,energ4x,myid,rprecision)
  use ctes
  implicit real*8(a-h,o-z)
  ! /*                                         */
  ! /*                                         */
  parameter(iout=34)       !some output
  character*256 filename
  character*2 rprecision
        
  dimension energ4x(mgalx,5)
  !-----|---------------------------------------------------------------|
  open(iout,file=filename,form="unformatted")
  if(rprecision.eq.'r8')then
     write(*,*)'precision not implemented in save: ',rprecision
     stop
  elseif(rprecision.eq.'r4')then
     write(iout) ((real((energ4x(i,is))),i=1,mgalx),is=1,5) 
  else
     write(*,*)'precision not implemented in save: ',rprecision
     stop
  endif
  close(iout)
  return
end subroutine save_energ4x
#endif /* SLONG */
!----|----------------------------------------------------------------|
subroutine save_spectra(uus2,vvs2,wws2,pps2,nstat, & 
     &                  timew,cflag,rp)
  use ctes
  implicit real*8(a-h,o-z)
  ! /*                                                           */
  ! /* write streamwise spectrum of uu, vv, ww  in binary file   */
  ! /* cflag='uu','vv','ww','pp','to'                            */
  ! /*       'to' means to save all(uu,vv,ww,pp) into one file   */
  ! /*         with header                                       */
  ! /*                                                           */
  parameter(iout=34)       !some output
  character*2 cflag,rp
  character*4 ext2
  character*256 filename
  character*11 filebase
  
  dimension uus2(my,mz,ib:ie),vvs2(my,mz,ib:ie),&
       &                  wws2(my,mz,ib:ie),pps2(my,mz,ib:ie)
  ispace=index(fuvwpout,' ')
  write(ext2,'(i4.4)') myid
  filename=fuvwpout(1:ispace-1)//cflag(1:2)//'spec2.'//ext2(1:4)
  open(iout,file=filename,form='unformatted')   
  if(myid.eq.master) write(*,*)' writing spectra (bin) ' 
  !trim(filename) does not work on PRIMERGY (tatara kyushuu) 
  if(rp.eq.'r8')then
     if(cflag.eq.'uu')then
        write(iout) uus2
     elseif(cflag.eq.'vv')then
        write(iout) vvs2
     elseif(cflag.eq.'ww')then
        write(iout) wws2
     elseif(cflag.eq.'pp')then
        write(iout) pps2
     elseif(cflag.eq.'to')then
        write(iout) timef,timew,fnu,alp,aspect,dfloat(nstat),&
             &      dfloat(mgalx),dfloat(my),dfloat(mz),&
             &      dfloat(ib),dfloat(ie),dfloat(myid) 
        write(iout) (((uus2(j,k,i),j=1,my),k=1,mz),i=ib,ie),&
             &      (((vvs2(j,k,i),j=1,my),k=1,mz),i=ib,ie),&
             &      (((wws2(j,k,i),j=1,my),k=1,mz),i=ib,ie),&
             &      (((pps2(j,k,i),j=1,my),k=1,mz),i=ib,ie)
     endif
  elseif(rp.eq.'r4')then
     if(cflag.eq.'uu')then
        write(iout) real(uus2)
     elseif(cflag.eq.'vv')then
        write(iout) real(vvs2)
     elseif(cflag.eq.'ww')then
        write(iout) real(wws2)
     elseif(cflag.eq.'pp')then
        write(iout) real(pps2)
     elseif(cflag.eq.'to')then
        write(iout) real(timef),real(timew),real(fnu),real(alp),&
             &      real(aspect),real(nstat), &
             &      real(mgalx),real(my),real(mz), &
             &      real(ib),real(ie),real(myid) 
        write(iout) (((real(uus2(j,k,i)),j=1,my),k=1,mz),i=ib,ie),&
             &      (((real(vvs2(j,k,i)),j=1,my),k=1,mz),i=ib,ie),&
             &      (((real(wws2(j,k,i)),j=1,my),k=1,mz),i=ib,ie),&
             &      (((real(pps2(j,k,i)),j=1,my),k=1,mz),i=ib,ie)
     endif
  else
     write(*,*) 'wrong rprecision: in write3d_spectrum'
  endif
  close(iout)
  return
end subroutine save_spectra
!----|----------------------------------------------------------------|
!    |  for DEBUG
!----|----------------------------------------------------------------|
subroutine diff_fis(a,b)
  use ctes
  implicit real*8(a-h,o-z)

  real*8 a(my,mz,mgalx+2),b(my,mz,mgalx+2)
  a_max=0.d0
  b_max=0.d0
  d=0.d0
  do i=1,mgalx
     do k=2,mz-1
        do j=2,my-1
           a_max=max(abs(a(j,k,i)),a_max)
           b_max=max(abs(b(j,k,i)),b_max)
           d=max(abs(a(j,k,i)-b(j,k,i)),d)
        enddo
     enddo
     write(*,*) a(49,49,i),b(49,49,i)
  enddo
  write(*,*) 'max diff:',d
  write(*,*) 'max a,b:',a_max,b_max
  
  return
end subroutine diff_fis

subroutine crop_real4(tmpfis,res,kb,ke)
  use ctes, only: nk,my,mzp,mgalx
  implicit none   

  real(nk) tmpfis(mgalx+2,my,mzp) ! reshape with local kpl
  real*4 res(mgalx,my,kb:ke)
  integer kb,ke

  res(1:mgalx,1:my,kb:ke)=real(tmpfis(1:mgalx,1:my,1:(ke-kb+1)))

end subroutine crop_real4
