#define DEBUG_main1 
#undef DEBUG_main1
#define DEBUG_main2 
#undef DEBUG_main2
#define DEBUG_main3 
#undef DEBUG_main3
#define DEBUG_main5 
#undef DEBUG_main5
#define DEBUG_nolap 
#undef DEBUG_nolap
#define TIMING_test 
#undef TIMING_test
#define SEKIPwritestat
#undef SEKIPwritestat

program main
  use hdf5
  use ctes
  use numbers
  use running
  use timing
  use massflow, only: cmassflow
  implicit none 
  !     /***************************************************************/
  !     /* DNS code for rectangular duct flow                          */
  !     /*                                                             */
  !     /*        +1----------------------------                       */
  !     /*           |            y            |                       */
  !     /*           |                 x       |                       */
  !     /*           |            ^            |                       */
  !     /*           |            |  /         |                       */
  !     /*           |            | /          |                       */
  !     /*           |            |/           |                       */
  !     /*         0-|            ------> z    |                       */
  !     /*           |                         |                       */
  !     /*           |                         |                       */
  !     /*           |                         |                       */
  !     /*           |                         |                       */
  !     /*           |                         |                       */
  !     /*           |                         |                       */
  !     /*        -1----------------------------                       */
  !     /*           |            |            |                       */
  !     /*        -aspect         0        +aspect                     */
  !     /*                                                             */
  !     /*                                                             */
  !     /* formulation:                                                */
  !     /* . primitive variables                                       */
  !     /*                                                             */
  !     /* time integration:                                           */
  !     /* . fractional step, semi-implicit, 3-step runge-kutta        */
  !     /* [1] rai & moin, jcp, vol 96, 15-53, 1991                    */
  !     /* [2] kim, kim, choi, jcp, vol 171, p 132ff, 2001             */
  !     /*                                                             */
  !     /* spatial discretization:                                     */
  !     /* . chebyshev collocation in (y,z)                            */
  !     /* . fourier in x                                              */
  !     /* . non-linear terms pseudo-spectrally, de-aliasing in x (2/3)*/
  !     /*                                                             */
  !     /* solver for 2d helmholtz/poisson problems (cheb-cheb):       */
  !     /* . fast matrix diagonalization technique                     */
  !     /* [3] haldenwang et al., jcp, vol 55, 115-128, 1984           */
  !     /*-------------------------------------------------------------*/
  !     /* NOTE: the inverse gauss-lobatto grid is used, i.e.          */
  !     /*   yi=-cos(i*pi/n), zj=-cos(j*pi/n)                          */
  !     /* this means for output variables in physical space that:     */
  !     /*   f(i,j)=f(yi,zj) with yi,zj in INCREASING order            */
  !     /*-------------------------------------------------------------*/
  !     /* AUTHORS                                                     */
  !     /*                                                             */
  !     /* markus uhlmann & alfredo pinelli                 march 2004 */
  !     /* atsushi sekimoto, markus uhlmann & alfredo pinelli          */
  !     /*                                     MPI version  July  2007 */
  !     /* GPU-duct, alltoall change with transpose, March 2020 by A.S */
  !     /*-------------------------------------------------------------*/
  include 'fftw3.f'
  ! /* MPI: include header file  */
  include 'mpif.h'
  ! /* variables: fourier-space size */
  ! /* MPI: nsfou is maximun size of local fou array,   */
  ! /*              cut in x-dir or in z-dir            */
  complex(nk),dimension(:),allocatable:: u,v,w,p
#ifdef TEMPERATURE 
  complex(nk),dimension(:),allocatable:: t
#endif  /* TMEPERATURE */   
  ! 
  ! /* 3d physical-size work arrays: */
  ! /* MPI: nsfis is size of local fis array (cut in z-dir) */
  real(nk),dimension(:),allocatable:: tmp1
  real(nk),dimension(:),allocatable:: tmp4, tmp5, tmp6
  ! /* 3d fourier-size work arrays: */
  real(nk),dimension(:),allocatable:: adv1,adv2,adv3
  real(nk),dimension(:),allocatable:: h1,h2,h3
  real(nk),dimension(:),allocatable:: tmp2,tmp3
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
  
  ! /* "small" arrays */
  real(nk) stress4(4),energomx4(5),energomxb4(5),stresst4(4)
  integer isave,istop,iimag,ipro,ifis 
  integer i
  integer iter,krk,idum
  integer iunit,nprocr
  real(nk) get_dpdx, dpdx

  parameter(iunit=34)       !some output
  character*3 csave
  parameter(csave='fou')

  real(nk) dum,energi,energomx,energomxb, & 
       energomy,energomyb,energomz,energomzb,energt,energtb, & 
       energu,energub,energv,energvb,energw,energwb
  real(nk) fmax,helco,anorm_max
  ! --- hdf5 ---
  integer h5err
  ! --- MPI ---
  integer ierr

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
#ifdef CUBLAS
  !if (mod(myid,32).eq.0) then
  !if (mod(myid,16).eq.0) then
  !if (mod(myid,8).eq.0) then
  !if (myid.lt.16) then
  !if (myid.eq.0) then
  ngpu=1
  call cublas_init
  !write(*,*) myid, ': use gpu id='
  !end if
#endif
  prof_shen=0.d0; prof_legcheb=0.d0; timeshen=0.d0
  prof_diag=0.d0; prof_para=0.d0; prof_total=0.d0

#ifdef PROFILO
     prof(:) = 0.d0
#endif

  istop=0
  if((myid.eq.master).and.(master.ne.0)) then
     write(*,*) 'master id must be 0: ',master
     write(*,*) 'please check ctes3D_in '
     istop = 1000 
  end if
  if(myid.eq.master) then
     write(*,*) ' '
     write(*,*)'MPI: number of processor:',nproc
     write(*,*)'DUCT FLOW CODE, MPI Version by Sekimoto & Markus '
     write(*,*)'CHEB-LEG-DUCT (testing) by A.Sekimoto '
     write(*,*)'GPU-acceleration by A.Sekimoto '
     call cpp(my,mz)
  end if
  
  call MPI_BCAST(istop,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
  if (istop.eq.1000) go to 1000
  !
  !     /* basic initial set-up */
  call set_numeric_constants()
  !
  call initcr
  !
  !write(*,*) myid, 'allocating buffers'
  ! allocating buffers
  allocate(u(nsfou),v(nsfou),w(nsfou),p(nsfou)) ! complex*16
  u=0.d0;  v=0.d0;  w=0.d0;  p=0.d0
#ifdef TEMPERATURE
  allocate(t(nsfou)) ! complex*16
  t=0.d0
#endif /* TMEPERATURE */    
  allocate(tmp1(nsfis) )
  tmp1=0.d0
  allocate(tmp4(nsfis),tmp5(nsfis),tmp6(nsfis))
  tmp4=0.d0; tmp5=0.d0; tmp6=0.d0
  
  allocate(adv1(nsfour),adv2(nsfour),adv3(nsfour))
  adv1=0.d0; adv2=0.d0; adv3=0.d0
  allocate(h1(nsfour),h2(nsfour),h3(nsfour))
  h1=0.d0; h2=0.d0; h3=0.d0
  allocate(tmp2(nsfour),tmp3(nsfour))
  tmp2=0.d0; tmp3=0.d0
#ifdef TEMPERATURE
  allocate(advt(nsfour))
  advt=0.d0
#endif
  ! tmp array for shuffling 
  allocate(tmpzxcut(my*mzp*mxp)) ! complex*16
  allocate(worke(NWKE*nthreads))
  allocate(tmpyz(my+mz)) ! complex*16
  allocate(tmpx(mgalx+2))
  tmpzxcut=0.d0; worke=0.d0; tmpyz=0.d0; tmpx=0.d0; 
  !
  !     /* initial time */
  time=f0o1
  ! ttarget=min(dtimag,tfin)
  isave=0
  iimag=0
  !
  ipro=-1
#ifdef HELM_SHEN
  call helm3d_shen(dum,dum,worke,ipro,'_','_',f1o1,aspect,dum)
  call leg_cheb_mat(worke,worke(my*mz+1),my)
#endif /* HELM_SHEN */
  call helm3d(dum,dum,worke,ipro,'d','_',f1o1,aspect,dum)
  ipro=1
  !     
  ! /* initialize the fields */
  if(irestart.eq.1) then     !read fields from file
     if(myid.eq.master) then            
        write(*,*)'field in fis [1] or fou [2] or hdf5(new) [3]?'
        !read(*,*) ifis
        ifis = 3
        write(*,*)'ifis=',ifis
     endif !/* only master id */
     call MPI_BCAST(ifis,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
     
     if (ifis.eq.1) then  ! read .bin file (p-p-p)
        call read_fft_4fields_fis(u,v,w,p, &
             &           tmp1,tmp4,tmpx,'r8','_')
#ifdef TEMPERATURE
        call read_fft_4fields_fis(t,dum,dum,dum, &
             &           tmp1,tmp4,tmpx,'r8','t')
#endif /* TEMPERATURE */        
     elseif (ifis.eq.2) then ! read .myid file (f-p-p)
        if (myid.eq.0) then
           write(*,*) 'ifis=2 set numproc in the file'
           read(*,*) nprocr
        end if
        call MPI_BCAST(nprocr,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
#ifdef TEMPERATURE
        call read_5fields_fou_apart(u,v,w,p,t,filinp,nprocr,'r8')
#else /* TEMPERATURE */
        call read_4fields_fou_apart(u,v,w,p,filinp,nprocr,'r8')
#endif /* TEMPERATURE */                
     elseif (ifis.eq.3) then ! hdf5 f-p-p read
        if (myid.eq.master) write(*,*) 'read from h5'
#ifdef TEMPERATURE
        call read_5fields_fou_hdf5(u,v,w,p,t,filinp,'r8')
#else
        !write(*,*) myid, "ifis=", ifis
        call read_4fields_fou_hdf5(u,v,w,p,filinp,'r8')
#endif
     else
        write(*,*) ' set ifis <= 3, stop'
        stop
     end if
     !
     iimag=iinp
     ! ttarget=min(time+dtimag,tfin)
  else
     write(*,*) 'init_field of TEMPARATURE is not parallelize'
  endif

  call div3d(u,v,w,h2,tmp2,tmpyz,f1o1) !f-p-p
  ! write(*,*) myid,"h2=",h2(1*(my*mz)+20*my+20) ! OK
  call max_3din(h2,divmax)
  if(myid.eq.master) then
     write(*,*)'INITIAL DIVERGENCE: ',divmax
  endif

#ifdef HELM_SHEN
  call L2_norm_fou2d_max(h2,anorm_max)
  if(myid.eq.master) then
     write(*,*)'INITIAL max L2norm of DIVERGENCE: ',anorm_max
  endif
#endif

  call statistics_fou(u,v,w,p,-1) ! for allocation (jinco,kinco are also allocated)
  ! reset statistics
  call statistics_fou(u,v,w,p,0)
  call statistics_omg2_fou(dum,-1) ! for allocation
  call statistics_omg2_fou(dum,0)
#ifdef TEMPERATURE
  call statisticst_fou(t,-1) ! for allocation
  call statisticst_fou(t,0)
#endif
  ipro=1
#ifdef SHOOTING
  ! /* recommend to use with the option, MAKE_PRESSURE*/
  if (myid.eq.0) then
     write(*,*)'shooting parameter: beta=?'
     read(*,*) beta
     write(*,*)'beta=',beta
  endif
  call MPI_BCAST(beta,1,MPI_DOUBLE_PRECISION,master, &
       &      MPI_COMM_WORLD,ierr)
  ! /* compute edge state: u = <u> + \beta (u - <u>) */
  call make_edge_state(u,v,w,beta)
#endif /* SHOOTING */
#ifdef MAKE_PRESSURE
     
  !  /* make the pressure field by solving:                          */
  !  /*        \nabla^2(p) = - \nabla \dot [NL(u)]                   */
  !  /*                                                              */
  if (myid.eq.master) then
     write(*,*) 'making pressure field from u, v, w.'
  endif
  !  /* (0) initializing adv1,adv2,adv3 */
  call dcopy(nsfour,f0o1,0,adv1,1)
  call dcopy(nsfour,f0o1,0,adv2,1)
  call dcopy(nsfour,f0o1,0,adv3,1)
  !  /* (1) compute advection terms, NL(u) ->adv1,adv2,adv3          */
  krk=0
  iter=0
  call advection3d(u,v,w,adv1,adv2,adv3, &
       &     tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmpzxcut,tmpx, &
       &     krk,iter)
  !  /* (2) compute r.h.s.:    h1 = -\nabla \dot [NL(u)]             */
  call div3d(adv1,adv2,adv3,h1,tmp2,tmpyz,-f1o1)
  !  /* (3) solve poisson equation  \nabla^2 =     */
#ifdef HELM_SHEN
  call helm3d_shen(-h1,p,worke,ipro,'n','_',f1o1,aspect,f0o1)
#else /* HELM_SHEN */
  call helm3d(h1,p,worke,ipro,'n','_',f1o1,aspect,f0o1)
#endif /* HELM_SHEN */
#endif /* MAKE_PRESSURE */

#ifdef CMASSFLOW
#ifdef RECTANGULAR
  if(myid.eq.master) then
     write(*,*) 'attention: CMASSFLOW might not work in this mode!'
     write(*,*) 'please check!!'
  endif
#endif
  call precomp_mass_flow(dtfixed,bkkc,worke,nwke)
  call impose_mass_flow(cmassflow,idum,u,worke,nwke,0)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif /* CMASSFLOW */
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
  !   /* set some unused array to zero */
  stresst4(1)=0.d0
  stresst4(2)=0.d0
  stresst4(3)=0.d0
  stresst4(4)=0.d0
  energt=0.d0
  energtb=0.d0
  !  --------------------------------------------------------------     
  !  /* ----- start the temporal loop ----------------------------- */
  !  --------------------------------------------------------------     
  do 100 iter=1,nstep
#ifdef TIMING
     t1_main=MPI_WTIME()
#endif /* TIMING */
     if(myid.eq.master) write(*,*)'STEP: ',iter,' TIME: ',time
     
     !     /* check for fixed-no. of iteration interval */
     if(mod(iter,nimag).eq.0)then
        isave=1
        if(myid.eq.master)write(*,*)' NEXT: SAVING...'
     endif
          
     do 200 krk=1,mrk       !runge-kutta step
#ifdef PROFILO
        ! start wall clock within each rk substep
        p0=MPI_WTIME()
#endif          
        !if (myid.eq.master) write(*,*) 'dscal2 '
        if(krk.gt.1)then    !rkkc(1)=0
           !  /* (i) store -\rho^k*adv(u^{k-2}) -> h1,h2,h3 */
           call dscal2(nsfour,adv1,1,-rkkc(krk),h1,1)
           call dscal2(nsfour,adv2,1,-rkkc(krk),h2,1)
           call dscal2(nsfour,adv3,1,-rkkc(krk),h3,1)
        else
           call dcopy(nsfour,f0o1,0,h1,1)
           call dcopy(nsfour,f0o1,0,h2,1)
           call dcopy(nsfour,f0o1,0,h3,1)
        endif
        !     
        ! /* (ii)  compute the non-linear terms of navier-stokes at 'k-1'   */
        ! /*      adv(u^{k-1}): -> adv1,adv2,adv3 (preserved for next step) */
        !if (myid.eq.master) write(*,*) 'advection3d '
        call advection3d(u,v,w,adv1,adv2,adv3, &
             &           tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmpzxcut,tmpx, &
             &           krk,iter)
#ifdef PROFILO
        p1=MPI_WTIME()
        prof(1)=prof(1) + (p1-p0)
        !accumulating time cost for advection3d (u,v,w)
#endif  
     
#ifdef TEMPERATURE
        !     /* (iib) add the buoyancy term (on l.h.s.) grav*alphat*T to adv   */
        !     /*     adv(u^{k-1})+grav*alphat*T^{k-1} -> adv1,adv2,adv3         */
        call daxpy(nsfour,gravx*alphat,t,1,adv1,1)
        call daxpy(nsfour,gravy*alphat,t,1,adv2,1)
        call daxpy(nsfour,gravz*alphat,t,1,adv3,1)
#endif /* TEMPERATURE */
        !     
        !     /* time-step related */
        if(krk.eq.1)then
#ifdef PROFILO
           p0=MPI_WTIME()
#endif   
#ifdef DT_FIXED
           if(myid.eq.master) then
              write(*,666) dtfixed,dtfixed/(dt/cfl)
666           format('DT_FIXED= ',e15.8,' -> CFL=',e15.8)
              if (dtfixed/(dt/cfl).gt.1.d2) then
                 write(*,*)  'Diverged'
                 istop=1000
                 stop
              endif
           endif
           dt=dtfixed
#endif /* DT_FIXED */     
           if(mod(iter-1,nhist).eq.0)then
              !  /* the normal-to-the-wall u-gradient, integrated over 4 walls */
              call wallstressx(u,tmp2,stress4)
#ifdef TEMPERATURE
              call wallstressx(t,tmp2,stresst4)
#endif /* TEMPERATURE */
     
              !  /* the rms of velocity fluctuations */
              call energy_fou(u,energu,energub)
              call energy_fou(v,energv,energvb)
              call energy_fou(w,energw,energwb)
              energi=0.d0
              ! call calc_energy_input_fou(u,p,energi)
              ! /* compute streamwise vorticity */
              ! /* compute enstrophy */
              !if (myid.eq.master) write(*,*) 'vorticity  '
              call vorticity(u,v,w,tmp2,tmp3,1) !compute omega_x->tmp2
              call statistics_omg2_fou(tmp2,1)
              call energy_fou(tmp2,energomx,energomxb) 
              !     /* intergal of (omx_xmean)^2 by sector */
              call energy_fou_zero_sector(tmp2,energomx4)
              !     /* intergal of abs(omx_xmean) by sector */
              call integral_abs_zero_sector(tmp2,energomxb4)
              call vorticity(u,v,w,tmp2,tmp3,2) !compute omega_y->tmp2
              call statistics_omg2_fou(tmp2,2)
              call energy_fou(tmp2,energomy,energomyb)
              call vorticity(u,v,w,tmp2,tmp3,3) !compute omega_z->tmp2
              call statistics_omg2_fou(tmp2,3)
              call energy_fou(tmp2,energomz,energomzb)    
#ifdef TEMPERATURE
              call energy_fou(t,energt,energtb)
#endif
              if(myid.eq.master) then
                 write(*,111)time, &
                      &   stress4(1)+stress4(2)+stress4(3)+stress4(4),&
                      &   energu+energv+energw,&
                      &   energu,energv,energw,energt,&
                      &   dsqrt(energomx4(1)+energomx4(2)+energomx4(3)&
                      &   +energomx4(4)+energomx4(5)),&
                      &   (stress4(i),i=1,4),&
                      &   (energomx4(i),i=1,4),(energomxb4(i),i=1,4),&
                      &   energomx4(5),energomxb4(5),&
                      &   stresst4(1)+stresst4(2)+stresst4(3)+stresst4(4),&
                      &   energub,energvb,energwb,energtb,& ! 24,25,26,27 are energy
                      &   (stresst4(i),i=1,4),&
                      &   energomx,energomy,energomz,&
                      &   energomxb,energomyb,energomzb,& ! 35,36,37
                      &   energi ! 38
               
              endif         !/*master ener output*/
              if (dt.lt.1.d-10) then 
                 write(*,*) 'time step: ',dt,'(cfl=',cfl,')'
                 istop=1000
              end if
              call MPI_BCAST(istop,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
              if (istop.eq.1000) go to 1000
              
              call statistics_fou(u,v,w,p,1)
#ifdef TEMPERATURE
              call statisticst_fou(t,1)
#endif
           endif
#ifdef PROFILO
           p1=MPI_WTIME()
           prof(2)= prof(2) + p1-p0
           ! accumulating time cost for energy and statistics
#endif 
        endif
     
#ifdef PROFILO     
        p0=MPI_WTIME()
#endif
        ! /* (iii) add -\gamma^k*adv(u^{k-1}) to h1,h2,h3 -> h1,h2,h3 */
        call daxpy(nsfour,-gkkc(krk),adv1,1,h1,1)
        call daxpy(nsfour,-gkkc(krk),adv2,1,h2,1)
        call daxpy(nsfour,-gkkc(krk),adv3,1,h3,1)
        !     
        ! /* (iv) compute \alpha_k\nu/2\nabla^2(u_i^{k-1}) (lapl. of vel.)  */
        ! /*     -> ++h1,++h2,++h3                                          */
#ifdef PROFILO
        p1=MPI_WTIME()
        prof(3)= prof(3) + p1-p0
        ! before laplacian3d
        p0=MPI_WTIME()
#endif    
        call laplacian3d(u,h1,tmp2,tmpyz,fnu*akkc(krk),'+')
        call laplacian3d(v,h2,tmp2,tmpyz,fnu*akkc(krk),'+')
        call laplacian3d(w,h3,tmp2,tmpyz,fnu*akkc(krk),'+')
#ifdef PROFILO
        p1=MPI_WTIME()
        prof(4) = prof(4) + p1-p0
        ! accumulating time cost for laplacian3d (u,v,w)
        p0=MPI_WTIME()
#endif  
        ! /* (v)   add u^{k-1}/dt to arrays h1,h2,h3 */
        call daxpy(nsfour,f1o1/dt,u,1,h1,1)
        call daxpy(nsfour,f1o1/dt,v,1,h2,1)
        call daxpy(nsfour,f1o1/dt,w,1,h3,1)
#ifdef PROFILO
        p1=MPI_WTIME()
        prof(3) = prof(3) + p1-p0
        p0=MPI_WTIME()
#endif       
        ! /* (vi) add -2\alpha_k*grad(p^{k-1}) to h1,h2,h3 */
        call grad3d(p,h1,h2,h3,tmp2,tmpyz, &
             &      -(akkc(krk)+bkkc(krk)),'+')
#ifdef PROFILO
        p1=MPI_WTIME()
        prof(5) = prof(5) + p1-p0
        p0=MPI_WTIME()
#endif     
        !     /* (vi-b) add driving pressure gradient to the real part  */
        !     /*         of the fourier-0-mode of the rhs of x-momentum */
        !     /*         -2\alpha_k*2dp/dx -> ++h1                      */
        !     /*         with: dp/dx=-\nu*fct(aspect) [tabulated]       */
        if(icase.eq.1) then  !gradient from 2d poiseuille profile
#if !defined(CMASSFLOW)
#if !defined(FIXED_ZERO_MODE)
           if(myid.eq.0) then
              !$$$  dpdx=-fnu*f4o1/f3o1*(f1o1+f1o1/aspect**3)
              dpdx=get_dpdx(aspect,fnu)
              write(*,*)'dpdx=',dpdx
              call dcopy(my*mz,-(akkc(krk)+bkkc(krk))*dpdx,0,tmp1,1)
              ! call dcopy(my*mz,-f2o1*akkc(krk)*dpdx,0,tmp1,1)
              call daxpy(my*mz,f1o1,tmp1,1,h1,2) !only 0-mode, real
           endif
#endif /* not FIXED_ZERO_MODE */
#endif /* not CMASSFLOW */
           
           !  call force_rolls(v,w,h2,h3,akkc(krk)) ! akkc + bkkc ??
#ifdef PUFF_FORCING
           !  /* add a localized forcing term to h1/2/3 */
           call add_puff_forcing(h1,h2,h3,-(akkc(krk)+bkkc(krk)), &
                &                tmp4,tmpx,tmp2,tmpzxcut)
#endif /* PUFF_FORCING */
#ifdef FORCING_GLOBAL
           if(myid.le.1) then
              ! /* add force terms globally to h1, h2, h3        */
              ! /* note: using zeromode and first fourier modes  */
              ! nforcing_step=int(10.d0/(1.9278334/4)/dt)            
              ! if (iter.lt.nforcing_step) then
              ss=abs(stress4(1)+stress4(2)+stress4(3)+stress4(4))
              shear_lam=1.7d0
              if (ss.lt.(shear_lam*1.1d0)) then ! switch
                 amplitude=0.06d0
                 if(myid.eq.master) write(*,*) &
                      &  'global forcing amplitude:',amplitude 
                 !     this amplitude shoud be this order 
                 !     to keep DIVERGENCE to small
                 call force_global(h1,h2,h3,-(akkc(krk)+bkkc(krk)), &
                      &   amplitude)
              endif         !/* forcing switch by dudn */
           endif            !/* using first 2 fourier modes */
#endif /* FORCING_GLOBAL */
        endif !/*icase*/
#ifdef PROFILO
        p0=MPI_WTIME()
#endif  
        ! /* (vii)  multiply h1,h2,h3 by factor -1/(\nu\alpha_k)  */
        call dscal2(nsfour,h1,1,-f1o1/(fnu*bkkc(krk)),h1,1) ! bkkc ?
        call dscal2(nsfour,h2,1,-f1o1/(fnu*bkkc(krk)),h2,1)
        call dscal2(nsfour,h3,1,-f1o1/(fnu*bkkc(krk)),h3,1)
#ifdef PROFILO
        p1=MPI_WTIME()
        prof(6) = prof(6) + p1-p0
        p0=MPI_WTIME()
#endif 
        !
        ! /* (viii) solve for predicted u^\ast -> u,v,w */
        ! /*        [hom. dirichlet conditions]         */
        !  helco=f1o1/(akkc(krk)*fnu*dt) ! for Crank--Nicolson
        helco=f1o1/(bkkc(krk)*fnu*dt)
        !write(*,*) myid,"before helm3d h1=",h1(1*(my*mz)+20*my+20)
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#ifdef HELM_SHEN
        call helm3d_shen(-h1,u,worke,ipro,'d','u',f1o1,aspect,helco) ! note: -h1
        call helm3d_shen(-h2,v,worke,ipro,'d','v',f1o1,aspect,helco) 
        call helm3d_shen(-h3,w,worke,ipro,'d','w',f1o1,aspect,helco) 
#else /* HELM_SHEN */
        call helm3d(h1,u,worke,ipro,'d','u',f1o1,aspect,helco)
        call helm3d(h2,v,worke,ipro,'d','v',f1o1,aspect,helco)
        call helm3d(h3,w,worke,ipro,'d','w',f1o1,aspect,helco)     
#endif /* HELM_SHEN */
        !write(*,*) myid,"after helm3d u=",u(1*(my*mz)+20*my+20)
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#ifdef PROFILO
        p1=MPI_WTIME()
        prof(7) = prof(7) + p1-p0
        !accumulating time cost for helm3d (u,v,w)
        p0=MPI_WTIME()
#endif       
        ! /* (ix)  compute the rhs of pseudo-pressure -> h1 */
        call div3d(u,v,w,h1,tmp2,tmpyz,f1o1/((akkc(krk)+bkkc(krk))*dt))
#ifdef PROFILO
        p1=MPI_WTIME()
        prof(8) = prof(8) + p1-p0
        p0=MPI_WTIME()
#endif        
        ! /* (x)   solve the poisson equation for \phi^{k} -> h2 */
        ! /*       [hom. neumann bcs]                            */
        !write(*,*) myid,"after div3d h1=",h1(1*(my*mz)+20*my+20)
        !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#ifdef HELM_SHEN
        call helm3d_shen(-h1,h2,worke,ipro,'n','_',f1o1,aspect,f0o1)
#else /* HELM_SHEN */
        call helm3d(h1,h2,worke,ipro,'n','_',f1o1,aspect,f0o1)
#endif /* HELM_SHEN */
!        write(*,*) myid,"after helm3d p=",h2(1*(my*mz)+20*my+20)
!        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
#ifdef PROFILO
        p1=MPI_WTIME()
        prof(7) = prof(7) + p1-p0
        !accumulating time cost for helm3d (poisson)
        p0=MPI_WTIME()
#endif        
        ! /* (xi) update velocity & pressure: --> u,v,w,p */
        call update_uvwp(u,v,w,p, &
             &           h2,h1,tmp2,tmpyz, &
             &           akkc(krk),bkkc(krk)) ! bkkc is added
#ifdef PROFILO
        p1=MPI_WTIME()
        prof(9) = prof(9) + p1-p0
        p0=MPI_WTIME()
#endif
#ifdef CMASSFLOW
        call impose_mass_flow(cmassflow,krk,u,worke,nwke,1)
#endif /* CMASSFLOW */        
        !  /* check div-free condition of velocity -> h2 */
        h2=0.d0
        tmp2=0.d0
        tmpyz=0.d0
        call div3d(u,v,w,h2,tmp2,tmpyz,f1o1)
#ifdef HELM_SHEN
        call L2_norm_fou2d_max(h2,anorm_max)
        if(myid.eq.master) then
           write(*,*)'max L2norm of DIVERGENCE: ',anorm_max
        endif
#endif /* HELM_SHEN */
        call max_3din(h2,divmax)
        if(myid.eq.master) then
           write(*,*)'DIVERGENCE: ',divmax
        endif
        
#ifdef TEMPERATURE
        ! /* solve the temperature equation with updated velocity field */
        if(krk.gt.1)then    !rkkc(1)=0
           ! /* (i) store -\rho^k*advt(u^{k-1},T^{k-2}) -> h1 */
           call dscal2(nsfour,advt,1,-rkkc(krk),h1,1)
        else
           call dcopy(nsfour,f0o1,0,h1,1)
        endif
        !     
        ! /* (ii)  compute the non-linear terms of temperature at 'k-1'   */
        ! /*      advt(u^{k},T^{k-1}): -> advt (preserved for next step)  */
        call advection_temp(u,v,w,t,advt,tmp1,tmp4,tmp2,tmpzxcut,tmpx,krk,iter)
        
        ! /* (iii) add -\gamma^k*adv(u^{k},T^{k-1}) to h1-> h1*/
        call daxpy(nsfour,-gkkc(krk),advt,1,h1,1)
        
        ! /* (iv) compute \alpha_k\kappa/2\nabla^2(T^{k-1}) -> ++h1 */
        call laplacian3d(t,h1,tmp2,tmpyz,fkappa*akkc(krk),'+')
        
        ! /* (v)   add T^{k-1}/dt to array h1 */
        call daxpy(nsfour,f1o1/dt,t,1,h1,1)
        
        ! /* (vi) add +2\alpha_k*qsource  [heat source] to h1 */
        if(myid.eq.0) then
           call dcopy(my*mz,(akkc(krk)+bkkc(krk))*qsource,0,tmp1,1)
           ! call dcopy(my*mz,f2o1*akkc(krk)*qsource,0,tmp1,1)
           call daxpy(my*mz,f1o1,tmp1,1,h1,2) !only 0-mode, real
        endif
        
        ! /* (vii)  multiply h1 by factor -1/(fkappa\alpha_k)  */
        call dscal2(nsfour,h1,1,-f1o1/(fkappa*bkkc(krk)),h1,1) ! bkkc ?
        
        ! /* (viii) solve for temperature T -> t */
        helco=f1o1/(bkkc(krk)*fkappa*dt)
        !   helco=f1o1/(akkc(krk)*fkappa*dt)
        call helm3d_temp(h1,t,worke,ipro,'t',f1o1, &
             aspect,helco,dum)
        call t_corner(t)
#endif /* TEMPERATURE */
        
        
200  enddo  !end runge-kutta sub-step
   
#ifdef TIMING
     t2_main=MPI_WTIME()
     tt=t2_main - t1_main
     
     prof_total=prof_total+tt

     if(myid.eq.master) write(*,*) ' :execution time for last step: ',tt
#endif /* TIMING */
     time=time+dt
     
     if(isave.eq.1)then
        if(myid.eq.master) write(*,*) 'writing results...'
        
        if(csave.eq.'fou')then
#ifdef TEMPERATURE 
           call save_5fields_fou_hdf5(u,v,w,p,t,iimag,'r8')
#else
           call save_4fields_fou_hdf5(u,v,w,p,iimag,'r8')
#endif
        else /* csave.ne.'fou'*/
           write(*,*) 'csave must be fou!'
           go to 1000 
        endif
        iimag=iimag+1
        !call write_statistics(iunit)
        call write_statistics_hdf5
        isave=0
     endif
     if(time.ge.(tfin-dsmall)) goto 101
100 enddo               !end full time step
101 continue        !jumped out of time loop

#ifdef PROFILO
     if(myid.eq.0) call out_profiling(iter)
#endif

! /* ---- end of temporal loop ----------------------------------*/
!     
111  format('ENER: ',40(e15.8,1x))
     
  !cc    fftw destroy
  !call dfftw_destroy_plan(plan_fft_forward)
  !call dfftw_destroy_plan(plan_fft_backward)
  !c      call dfftw_destroy_plan(plan_tcheb_y)
  !c      call dfftw_destroy_plan(plan_tcheb_z)
  ! deallocate arrays
1000 call MPI_FINALIZE(ierr)
  stop
end program main

#ifdef PROFILO 
subroutine out_profiling(iter)
  use timing 
  implicit real*8 (a-h,o-z)
  
  integer iter

  write(*,*) '   '
  write(*,*) '=== profiled results per a step in average ==='
  prof=prof/real(iter-1)
  write(*,*) 'time: advection =', prof(1)
  write(*,*) 'time: ener+stat =', prof(2)
  write(*,*) 'time: daxpy =', prof(3)
  write(*,*) 'time: laplacian3d =', prof(4)
  write(*,*) 'time: grad3d =', prof(5)
  write(*,*) 'time: dscal2 =', prof(6)
  write(*,*) 'time: helm3d =', prof(7)
  write(*,*) 'time: div3d =', prof(8)
  write(*,*) 'time: update_uvwp =', prof(9)

  write(*,*) '-----------------------------'
  write(*,*) 'total time: ', prof_total 
  write(*,*) ' ave. time: ', prof_total/real(iter-1) 

  prof_para=prof_para
  write(*,*) 'time: chx2z =', prof_para(1)
  write(*,*) 'time: chz2x =', prof_para(2)
  write(*,*) 'time: transpose (cyzx2xyz) =', prof_para(3)
  write(*,*) 'time: transpose (cxyz2yzx) =', prof_para(4)
  write(*,*) '---'
  write(*,*) 'time: comm =',(prof_para(1) + prof_para(2))/prof_total
  write(*,*) 'time: trans =',(prof_para(3) + prof_para(4))/prof_total
  write(*,*) '-----------------------------'

  return
end subroutine out_profiling
#endif  
