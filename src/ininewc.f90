subroutine initcr
  use ctes
  use running
  use statistics
  use rfttmp
  use tcheb2tmp, only:tmpyr,tmpzr
  use numbers
  use massflow
  implicit none 
  include "mpif.h"
  include"fftw3.f"
  ! /***************************************************************/
  ! /* initializes things                                          */
  ! /***************************************************************/
  ! /* MPI: tmp array for set parameter */
  real(nk) dat(20)
  integer  idat(10)
  character*80 cdat(2)
  
  integer j,k,ierr
  character*80 text
  !
  !---------------------------------------------------------------|
  ! /* reads in data MPI: only for master id          */
  if(myid.eq.master) then
     open(19,file='runparams2.dat',status='old')
     call set_parameters(19,icase,mgalx,my,mz,dat,idat,cdat)
     close(19)
     write(*,*)
     write(*,*)'runparams.dat is read by master id'
     write(*,*)
  end if

  ! /* MPI: BROADCAST parameter(only float data) to other ranks */
  call MPI_BCAST(mgalx,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(my,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(mz,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  
  call MPI_BCAST(dat,20,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(idat,10,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(cdat,2*80,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)

  call initg ! initialize grid 
  
  ! /* initializes fast fourier transforms */
  allocate(fdum(mgalx+2),bdum(mgalx+2)) 
  call dfftw_plan_dft_r2c_1d(plan_fft_forward,mgalx,fdum,fdum,FFTW_MEASURE)
  call dfftw_plan_dft_c2r_1d(plan_fft_backward,mgalx,bdum,bdum,FFTW_MEASURE)

  allocate(tmpyr(my),tmpzr(mz))
  tmpyr=0.d0; tmpzr=0.d0
  !
  ! /* initialize chebyshev transform stuff */
  ! call tchebi(my,chebsavy)
  ! call tchebi(mz,chebsavz)
  call tchebi_fftw(my,mz) ! tchebw.F

  allocate(premassfu(my*mz,mrk)) ! for cmassflow
  premassfu=0.d0

  ! /* MPI: set parameter */
  aspect=dat(1)
  alp = dat(2)
  fnu = dat(3)
  tfin = dat(4)
  dtimag = dat(5)
  cfl = dat(6)
  dtfixed = dat(7)
#ifdef TEMPERATURE
  alphat = dat(8)
  gravx = dat(9)
  gravy = dat(10)
  gravz = dat(11)
  fkappa = dat(12)
  qsource = dat(13)
#endif /* TEMPERATURE */

  icase = idat(1)
  nstep = idat(2)
  nimag = idat(3)
  nhist = idat(4)
  irestart = idat(5)
  iinp = idat(6)

  filinp = cdat(1)
  fuvwpout = cdat(2)
  
  allocate(xalp(0:mx1))
  xalp=0.d0
  allocate(iax(mx))
  iax=0
  allocate(y(my),z(mz),x(mgalx))
  y=0.d0; z=0.d0; x=0.d0
  allocate(dy(my),dz(mz))
  dy=0.d0; dz=0.d0
  ! /* compute the y- and z-coordinates */
  call dglc(-f1o1,f1o1,y,my)
  call dglc(-aspect,aspect,z,mz)
  call genexp ! setting xalp 
  
  ! /* compute the mesh size arrays in y/z for fis-space integration */
  dy(1)=(y(2)-y(1))*f1o2
  do j=2,my-1
     dy(j)=(y(j+1)-y(j-1))*f1o2
  enddo
  dy(my)=(y(my)-y(my-1))*f1o2
  
  dz(1)=(z(2)-z(1))*f1o2
  do k=2,mz-1
     dz(k)=(z(k+1)-z(k-1))*f1o2
  enddo
  dz(mz)=(z(mz)-z(mz-1))*f1o2
  
#ifdef TEMPERATURE
#ifdef CMASSFLOW
  if(abs(gravx).gt.1.d-10) then
     if(myid.eq.master) write(*,*) 'gravx should be zero with CMASSFLOW'
     call MPI_FINALIZE(ierr)
     stop
  endif
#endif /* CMASSFLOW */
#endif /* TEMPERATURE */
  ! /* write to screen */
  if(myid.eq.master) then
     write(*,'(a7,f8.5,a8,f6.3,a8,f6.3)') '  nu =',fnu,'channel'
     write(*,'(a7,f8.3,a8,f6.3,a8,f6.3)') 'alp =',alp,' aspect ratio = ',aspect
#ifdef TEMPERATURE
     write(*,*)' alphat =',alphat
     write(*,*)' gravity=',gravx,gravy,gravz
     write(*,*)' fkappa =',fkappa
     write(*,*)' qsource=',qsource
#ifdef CMASSFLOW
#ifdef RECTANGULAR
     write(*,*)' --> Reynolds = CMASS/(4*aspect)/fnu (based on h)'
     write(*,*)' --> Prandtl  =',fnu/fkappa    
     write(*,*)' --> Grashof(based on T_h - T_c and 2h)  =', &
          &  alphat*dsqrt(gravx**2+gravy**2+gravz**2)*f2o1*(2**3)/(fnu**2)
#else
     write(*,*)' --> Reynolds = CMASS/4/fnu (based on h)'
     write(*,*)' --> Prandtl  =',fnu/fkappa    
     write(*,*)' --> Grashof(based on T_h - T_c and 2h)  =', &
          &      alphat*dsqrt(gravx**2+gravy**2+gravz**2)*f2o1*(2**3)/(fnu**2)
#endif /* RECTANGULAR */
#else 
     if(dabs(qsource).lt.1.d-10) then
        write(*,*)' --> nu^-1 = f1o1/fnu'
        write(*,*)' --> Prandtl  =',fnu/fkappa
        write(*,*)' --> Grashof(based on T_h - T_c and 2h)  =', &
             &  alphat*dsqrt(gravx**2+gravy**2+gravz**2)*f2o1*(2**3)/(fnu**2)
     else
        !- inertial heating version
        write(*,*)' --> Reynolds =',f1o1/fnu
        write(*,*)' --> Prandtl  =',fnu/fkappa
        write(*,*)' --> Grashof  =', &
             & qsource*alphat*dsqrt(gravx**2+gravy**2+gravz**2)/(f2o1*fnu**2*fkappa)
        write(*,*)' --> dpdx/nu=',get_dpdx(aspect,fnu)/fnu
     endif
#endif /* CMASSFLOW */
#endif /* TEMPERATURE */
     write(*,*)
      
     write(*,'(a8,i5,a8,i5,a8,i5)') 'mgalx =',mgalx,'mz =',mz,'my =',my
     write(*,'(a8,i5,a8,i5,a8,i5)') 'mx =',mx,'mx1 =',mx1
     write(*,*)       
     write(*,'(a10,i6,a9,i6,a9,i5)') 'nstep =',nstep,'nimag =',nimag,'nhist =',nhist
     write(*,'(a10,e12.4,a8,f5.2)')  'dt   =',dtfixed,'  CFL =',CFL
     write(*,*)
     
     write(*,'(a12,i2,a9,a,a21,i5)') 'irestart =',irestart,'    from:',trim(filinp),&
          &       ' starting with no.:',iinp
     write(*,*)
  end if !/*  MPI: for master i/o */
  
  return
end subroutine initcr
!-----|---------------------------------------------------------------|
subroutine set_parameters(iunit,icase,nx,ny,nz,dat,idat,cdat)
  use ctes, only:nk
  implicit none
  !  /***************************************************************/
  !  /* reads the parameters of the problem from the input file.    */
  !  /* the file can contain lines starting with '%'or '#' which    */
  !  /* are considered as comments and are ignored.                 */
  !  /*                                                             */
  !  /* input                                                       */
  !  /* iunit          iunit for input from file                    */
  !  /*                                                             */
  !  /* output                                                      */
  !  /* [set through common blocks]                                 */
  !  /***************************************************************/
  !     /* local variable */
  character*80 text
  ! /* MPI: tmp array of parameter broad cast  */
  real(nk) dat(20)
  integer  idat(10)
  character*80 cdat(2)
  integer iunit,icase,nx,ny,nz
  
  !-----|---------------------------------------------------------------|
100 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 100
  read(text,*) nx,ny,nz

  !  /* (1) read the problem type: 1 integer                        */
                                !1   laminar profile, from 2d poiseuille
                                !2   lid-driven cavity, aspect=1
101 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 101
  !read(text,66) icase
  read(text,*) idat(1)
     
  !  /* (2) read the aspect ratio: "aspect"  - 1 float */
102 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 102
  !read(text,66) aspect
  read(text,66) dat(1)
  !     
  !  /* (3) read the streamwise wavenumber: "alp" - 1 float */
103 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 103
  !  read(text,66) alp
  read(text,66) dat(2)
  !     
  !  /* (4) read the fluid viscosity: "fnu" 1 float   */
110 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 110
  !   read(text,66) fnu
  read(text,66) dat(3)
  !     
  !  /* (5) read the max. no. of iterations: "nstep" 1 integer  */
115 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 115
  !  read(text,*) nstep
  read(text,*) idat(2)
  !     
  !  /* (6) read the iteration interval for writing images: 1 integer */
116 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 116
  !  read(text,*) nimag
  read(text,*) idat(3)
  !
  !  /* (6b) read the iteration interval for history: 1 integer */
117 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 117
  !read(text,*) nhist
  read(text,*) idat(4)
  !     
  !  /* (10) read the final time: 1 float */
118 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 118
  !  read(text,66) tfin
  read(text,66) dat(4)
  !     
  !  /* (11) read the time interval for writing images: 1 float */
119 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 119
  !  read(text,66) dtimag
  read(text,66) dat(5)
  !     
  !  /* (12) read the cfl number: 1 float */
120 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 120
  !  read(text,66) cfl
  read(text,66) dat(6)
  !     
  !  /* (13) read the fixed time-step: 1 float */
 121     read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 121
  !  read(text,66) dtfixed
  read(text,66) dat(7)
  !     
  !  /* (14) read the restart [1] from file: 1 integer             */
122 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 122
  !   read(text,*) irestart
  read(text,*) idat(5)
  !     
  !  /* (15) read the no. of first image to file: 1 integer        */
123 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 123
  ! read(text,*) iinp
  read(text,*) idat(6)
  !     
  !  /* (16) read the restart filename for uvwp-field: 1 character*80 */
130 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 130
  !  read(text,20) filinp
  read(text,20) cdat(1)
  !     
  ! /* (16b) read the output base-dir for uvwp-field: 1 character*80 */
  ! /* attention: needs trailing "/" e.g. "./" = here                */
131 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 131
  !  read(text,20) fuvwpout
  read(text,20) cdat(2)
#ifdef TEMPERATURE
  ! /* (20) read the coeff of thermal expansion: 1 float  (alphat) */
200 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 200
  !  read(text,66) alphat
  read(text,66) dat(8)
  !  /* (21) read the gravity components: 3 floats (gravx,gravy,gravz) */
201 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 201
  !  read(text,66) gravx
  read(text,66) dat(9)
202 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 202
  ! read(text,66) gravy
  read(text,66) dat(10)
203 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 203
  !  read(text,66) gravz
  read(text,66) dat(11)
  !  /* (22) read the coeff of thermal conductivity: 1 float (fkappa) */
204 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 204
  !  read(text,66) fkappa
  read(text,66) dat(12)
  !  /* (23) read the intensity of the heat source: 1 float (qsource) */
205 read(iunit,10) text
  if(text(1:1).eq.'#'.or.text(1:1).eq.'%')goto 205
  !  read(text,66) qsource
  read(text,66) dat(13)
#endif /* TEMPERATURE */
!     
66 format(10G15.3)
10 format(A)
20 format(a80)
       
  return
end subroutine set_parameters
! -------------------------------------------------------------- |


