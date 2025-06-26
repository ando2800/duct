module ctes
  use omp_lib

  ! ****************************************************************
  !       things which are fixed for the whole run
  ! ****************************************************************
  integer,parameter:: nk=8
  ! nk=4 is not supported
  integer,parameter:: nks=8 ! for statistics
  integer,parameter:: nkr=8 ! for timestepping coefficients
  !
  ! 4, single precision; 8, double precision; 16, double-double (slow)
  integer MPI_REAL_NP
  ! ----- dimensions ---------
  !     /********* parallel version  ***********************************/
  !     /* global                                                      */
  !     /* array size in physical space:                               */
  !     /*    (my,mz,mgalx+2) (real)                                   */
  !     /* array's size in fourier space, before dealiasing            */
  !     /*    (0:mx1)*mz*my (complex)                                  */
  !     /* array's size in fourier space, extended for dealiasing      */
  !     /*    (0:mgx)*mz*my (complex)                                  */
  !     /***************************************************************/
  !     /* (my,mz are the numbers of chebyshev polys in y/z dirs)      */
  !     /* (mx must be approximatly mgalx*2/3)                         */
  !     /* warning: there is no dealiasing in y,z                      */
  !     /* note:    2d in (y,z) corresponds to 'mgalx=4'               */
  !     /***************************************************************/
  !     /* user-adjustable dimensions: */  
  integer mgalx,my,mz
  integer nproc,numerop,myid,ib,ie,kb,ke,lastid,ngpu
  integer ipb,ipe,kpb,kpe
  integer mx, mgalx1, mx1, my1, mz1, mgx
  integer mxp,mxp_fis,mzp

  !--------Stuff for doing the OMP FFTw & Cosinus transform
  integer ompid,nthreads,chunk1
  !$OMP THREADPRIVATE(ompid,nthreads,chunk1)
  integer ibf1,ief1
  !$OMP THREADPRIVATE(ibf1,ief1)

  integer nsfis,nsfoux,nsfouz,nsfou,nsfour
  integer nwke1,nwke2,nwke3
!#ifdef HELM_SHEN
  integer mye,myo,mze,mzo
  integer nwks, nwke
!#endif
  real(nk), dimension(:), allocatable:: xalp
  integer, dimension(:), allocatable:: iax 
  !------------- Parameters for SMR R-K (A.A.Wray etc.)----------! 
  real(nk) akkc(3), bkkc(3), gkkc(3), rkkc(3), tkkc(3)
  parameter (akkc =(/ 29d0/96d0, -3d0/40d0,   1d0/6d0 /)) ! alpha
  parameter (bkkc =(/ 37d0/160d0, 5d0/24d0,   1d0/6d0 /)) ! beta
  parameter (gkkc =(/ 8d0/15d0,   5d0/12d0,   3d0/4d0 /)) ! gama
  parameter (rkkc =(/ 0.d0, -17d0/60d0, -5d0/12d0/)) ! zeta
  parameter (tkkc =(/ 8d0/15d0,   2d0/3d0,    1d0/)) 
  integer mrk
  parameter(mrk=3)
  ! ---- parallelization and blocking -----
  integer master
  parameter (master=0)
  integer, dimension(:), allocatable :: ibeg, iend, kbeg, kend
  integer, dimension(:), allocatable :: ipbeg, ipend, kpbeg, kpend
  real(nk), dimension(:), allocatable :: y, z, x, dy, dz
  real(nk) fnu,alp,aspect ! /fis/
  real(nk) alphat,gravx,gravy,gravz,fkappa,qsource ! /fis2/

  integer nimag,nstep,nhist,niter ! /timacc/
  !real(nk) dtfixed,CFL,time,tfin,dtimag ! /tem/
  real(nk) tfin,dtimag ! /tem/
  integer  irestart,irun,nstart,iinp,isn
  character*256  filinp,fuvwpout
  character*256  filstt

  real(nk),dimension(:,:,:),allocatable:: bcy,bcz
  character*1 boundy1,boundy2,boundz1,boundz2

contains
  subroutine initg
    use omp_lib

    implicit none
    include 'mpif.h'
    integer id
    ! /* derived global dimensions: */
    mx =2*(mgalx/3)
    mgalx1=mgalx-1
    mx1=mx/2-1 
    my1=my-1
    mz1=mz-1
    mgx=mgalx/2

    ! /* MPI: get pointers */
    allocate(ibeg(0:nproc-1),iend(0:nproc-1),kbeg(0:nproc-1),kend(0:nproc-1))
    ibeg=0; iend=0; kbeg=0; kend=0
    allocate(ipbeg(0:nproc-1),ipend(0:nproc-1),kpbeg(0:nproc-1),kpend(0:nproc-1))
    ipbeg=0; ipend=0; kpbeg=0; kpend=0
    call pointer1D(0,mx1,nproc,ibeg,iend) !!complex fourier space size!!!
    call pointer1D(1,mz,nproc,kbeg,kend)

    call pointer1D(1,mgalx+2,nproc,ipbeg,ipend) !!physical space size!!!
    call pointer1D(1,mz,nproc,kpbeg,kpend)
#ifdef DEBUG_1
    if(myid.eq.master) then
       write(*,*) 'ibeg:',(ibeg(id),id=0,lastid)
       write(*,*) 'iend:',(iend(id),id=0,lastid)
       write(*,*) 'kbeg:',(kbeg(id),id=0,lastid)
       write(*,*) 'kend:',(kend(id),id=0,lastid)
    end if
#endif /* DEBUG_1 */

    ib = ibeg(myid) ! index starts from 0
    ie = iend(myid)
    
    kb = kbeg(myid) ! index starts from 1
    ke = kend(myid)
    
    ipb = ipbeg(myid)
    ipe = ipend(myid)

    kpb = kpbeg(myid)
    kpe = kpend(myid)

    !$OMP PARALLEL
    !nthreads=OMP_GET_NUM_THREADS()
    !ompid=OMP_GET_THREAD_NUM()
    nthreads=1 ! this is for flat MPI...
    ompid=0    ! this is for flat MPI...
    chunk1=((ie-ib+1)+nthreads-1)/nthreads
    ibf1=ompid*chunk1+ib
    ief1=min((ompid+1)*chunk1+ib-1,ie)    
    if (nthreads.gt.1) then
       write(*,'(i10,i10,a20,3i10,a1,i2)') ib,ie,'ibf1,ief1,ompid',ibf1,ief1,ompid,'/',nthreads-1
    endif
    !$OMP END PARALLEL


    ! /* local fourier space dimension in x-direction */
    mxp=maxval(iend-ibeg+1) 
    !mxp=(mx1+1)/nproc+1
    mxp_fis=maxval(ipend-ipbeg+1)   !(mgalx+2)/nproc+1
    ! /* local dimension in z-direction */
    mzp=maxval(kend-kbeg+1) 
    !mzp= mz/nproc+1
    
    !write(*,*) myid,ib,ie,'mxp',mxp,'mzp',mzp

    ! /* 3d array sizes: in physical and fourier space */
    nsfis=(mgalx+2)*my*mzp !real
    nsfoux=mxp*my*mz !complex, cut in x-dir
    nsfouz=(mx1+1)*my*mzp !complex, cut in z-dir
    
    nsfou=max(nsfoux,nsfouz) !complex, both cuts
    !nsfou=(nsfoux+nsfouz) !complex, both cuts
    !call MPI_ALLREDUCE(nsfou,nsfou_mx,1,MPI_INTEGER, &
    !     &     MPI_MAX,MPI_COMM_WORLD,ierr)     
    !nsfou=nsfou_mx
    nsfour=2*nsfou
    
    nwke1=(2*mz-2)+(mz-1)**2+(mz+1)**2
    nwke2=(2*my-2)+(my-1)**2+(my+1)**2
    nwke3=(my+1)*(mz+1)+2*(my-1)*(mz-1)
    
    mye=my1/2
    myo=my1/2-1
    mze=mz1/2
    mzo=mz1/2-1

    nwks=mye*2+myo*2+(my+1)*(mz+1)
    nwke=max(max(nwke1,nwke2),nwke3)+nwks
    !nwke=nwke1+nwke2+nwke3+nwks
    
  end subroutine initg

end module ctes

module running
  use ctes, only: nk, nkr
  real(nkr) slipmx,cfl,dt,dtfixed,time,divmax
  integer icase
end module running

module massflow ! used to be common
  use ctes, only: nk, nkr,my,mz,mrk
  real(nk) cmassflow, premassf(mrk)
  real(nk),dimension(:,:),allocatable:: premassfu 
end module massflow

module statistics

  use ctes, only: nk,nks
  real(nk) timef
  integer nstat,nspec
  !/* number of positions in cross-plane for streamwise-spectra */
  parameter(nspec=7)
  real(nks),dimension(:,:),allocatable:: um,vm,wm,pm
  real(nks),dimension(:,:),allocatable:: uub,uvb,uwb,vvb,vwb,wwb,ppb
  real(nks),dimension(:,:),allocatable:: tm,ttb,utb,vtb,wtb
  ! spec
  real(nks),dimension(:,:),allocatable:: uus,vvs,wws,pps
  real(nks),dimension(:,:),allocatable:: tts
  integer(nks),dimension(:),allocatable:: jinco,kinco
  ! 
  real(nks),dimension(:,:),allocatable:: omx2,omy2,omz2 ! (my,mz)

end module statistics

module numbers

  implicit none
  real*8 f0o1,f1o1,f2o1,f3o1,f4o1,f5o1,f6o1,f7o1,f8o1,f9o1,f1o2,f1o3,f2o3,f1o4,f3o2
  real*8 dsmall,dbig,pi
  complex*16 zero,cii,cone

contains
  subroutine set_numeric_constants()

    ! /* numbers */
    f1o1=1.d0
    f2o1=2.d0
    f3o1=3.d0
    f4o1=4.d0
    f5o1=5.d0
    f6o1=6.d0
    f7o1=7.d0
    f8o1=8.d0
    f9o1=9.d0
    f0o1=0.d0
    f1o2=.5d0
    f1o3=1.d0/3.d0
    f2o3=2.d0/3.d0
    f1o4=1.d0/4.d0
    f3o2=3.d0/2.d0

    ! /* reals */
    dsmall=1.d-8
    dbig=1.d+8
    pi=4.d0*datan(1.d0)

    ! /* complex */
    zero=dcmplx(0.d0,0.d0)
    cii=dcmplx(0.d0,1.d0)
    cone=dcmplx(1.d0,0.d0)

  end subroutine set_numeric_constants
end module numbers

module rfttmp
  ! ****************************************************************
  !   fftw tmp arrays and plans for rft 
  ! ****************************************************************  
  use ctes, only: nk
  real(nk),allocatable:: fdum(:),bdum(:)
  real(nk)  dnf,dnb
  integer   nf,nb
  integer*8 plan_fft_forward,plan_fft_backward
  !integer isuccess

end module rfttmp

module tcheb2tmp
  ! ****************************************************************
  !   fftw tmp arrays for chebyshev tranformation  
  ! ****************************************************************  
  use ctes, only: nk
  real(nk), allocatable:: cdumy(:), cdumy2(:)
  real(nk), allocatable:: cdumz(:), cdumz2(:)
  integer*8 planr_tcheby, planc_tcheby
  integer*8 planr_tchebz, planc_tchebz
  integer  nny, nnyc, nnz, nnzc
  real(nk) dny1, dny2, dnz1, dnz2

  real(nk), allocatable:: tmpyr(:), tmpzr(:)

end module tcheb2tmp

module chebystuff
  use ctes, only: nk
  ! operators, pseudomat, solvediag
  !      common /chebystuff/ DCy(my,my),D2bx(my,2),D2by(my,2),x(my)
  !      common /chebysvd/  DTB(my-2,my-2),DS(my-2),DU(my-2,my-2),
  !     $     DVT(my-2,my-2),DSP(my-2,mz-2),DCM1(my-2,my-2),isvdflag
  real(nk),dimension(:,:),allocatable:: DC, D2bx, D2by
  ! for make_wfou using SVD
  real(nk),dimension(:,:),allocatable:: DTB,DS,DU,DVT,DSP,DCM1
  integer isvdflag

end module chebystuff

!#ifdef TEMPERATURE
module chebystuff_temp
  ! helm3d_temp, pseudomat_temp, solvediag_temp
  use ctes, only: nk
  !/chebystuffy_temp/ DY(my,my),D2by(my,2),y(my),sy
  !/chebystuffz_temp/ DZ(mz,mz),D2bz(mz,2),z(mz),sz
  real(nk),dimension(:,:),allocatable:: DY,D2by
  real(nk),dimension(:,:),allocatable:: DZ,D2bz
  real(nk) sy,sz

end module chebystuff_temp

module fulldiagtemp
  use ctes, only:nk
  ! /fulldiagtempy/ ETY(my-2),PTY(my-2,my-2),PTYM1(my-2,my-2)
  real(nk),dimension(:,:),allocatable:: PTY,PTYM1
  real(nk),dimension(:),allocatable:: ETY

  ! /fulldiagtempz/ ETZ(mz-2),PTZ(mz-2,mz-2),PTZM1(mz-2,mz-2)
  real(nk),dimension(:,:),allocatable:: PTZ,PTZM1
  real(nk),dimension(:),allocatable:: ETZ

end module fulldiagtemp

!#endif 

! for helm
module fulldiagneu

  use ctes, only: nk
  real(nk),dimension(:),allocatable:: EN 
  real(nk),dimension(:,:),allocatable:: EN2
  real(nk),dimension(:,:),allocatable:: PN, PNM1, PNT, PNM1T
  ! APS --> PNT in pseudopois
  real(nk) wanted 
  real(nk) AM1(4,4)
  integer iofzero
  !EN(my-2),PN(my-2,my-2),PNM1(my-2,my-2),
  !   $     PNT(my-2,my-2),PNM1T(my-2,my-2)

end module fulldiagneu

module fulldiagdir

  use ctes, only: nk
  real(nk),dimension(:),allocatable:: ED
  real(nk),dimension(:,:),allocatable:: ED2
  real(nk),dimension(:,:),allocatable:: PD, PDM1, PDT, PDM1T
  !ED(my-2),PD(my-2,my-2),PDM1(my-2,my-2),
  !   $     PDT(my-2,my-2),PDM1T(my-2,my-2)

  real(nk),dimension(:,:),allocatable:: tmprt

end module fulldiagdir

module timing

  ! for timing and profiling
  use ctes, only: nk
  real(nk) prof_shen(10) 
  ! 1, set_mat_F; 2, dcopy; 3, dgemm; 4, dcopy; 5, get_mat_F
  real(nk) prof_legcheb(10) 
  ! 1, set_mat_U; 2, dcopy; 3, dtrmm_l; 4, dcopy; 5, dtrmm_r; 6, dcopy; 7, get_mat_V 
  real(nk) timeshen(10)
  real(nk) prof_solV(11)

  real(nk) prof_diag(10) 
  real(nk) prof(0:20)
  real*8 t1,t2,t1_main,t2_main,tt,p0,p1

  real(nk) prof_para(4) 
  ! 1, chx2z; 2, chz2x; 
  ! 3, trans (cyzx2xyz); 4, trans (cxyz2yzx)

  real(nk) prof_total 

end module timing

! ---- for helmshen_gpu  -----
module legchebshen
  
  use ctes, only: nk
  ! matrices for LEGENDRE--CHEBYSHEV transformation
  ! on device
  real(nk),dimension(:,:),allocatable:: LCAE,LCAO,LCBE,LCBO
  real(nk),dimension(:,:),allocatable:: LCA,LCB

  real(nk),dimension(:),allocatable:: ED,EN ! for preshen
  real(nk),dimension(:,:),allocatable:: PD,PN,PSM1D,PSM1N ! for preshen

  ! matrices for 2nd-order lap
  real(nk),dimension(:),allocatable:: bkd,bkn ! boundary condition
  real(nk),dimension(:),allocatable:: dld2e,dld2o,bld2e,bld2o ! even-odd
  real(nk),dimension(:),allocatable:: dln2e,dln2o,bln2e,bln2o ! even-odd
  real(nk),dimension(:),allocatable:: MED,MDO,MPPED,MPPOD,SED,SOD ! MOD=>MDO
  real(nk),dimension(:),allocatable:: MEN,MON,MPPEN,MPPON,SEN,SON

  real(nk),dimension(:),allocatable:: DDE,DDO,DPE,DPO ! tmp buffs allocated in preshen
  real(nk),dimension(:),allocatable:: DDEp,DDOp,DPEp,DPOp ! tmp arrays
  real(nk),dimension(:),allocatable:: DPEp2,DPOp2 ! tmp arrays
  real(nk),dimension(:),allocatable:: gp,vp ! 1d tmp
  real(nk),dimension(:),allocatable:: fe,fo ! 1d tmp

  real(nk),dimension(:),allocatable:: Lae,Lao,Lbe,Lbo,Dce,Dco

  ! on device
  real(nk),dimension(:,:),allocatable:: PED,POD,PSM1ED,PSM1OD
  real(nk),dimension(:,:),allocatable:: PEN,PON,PSM1EN,PSM1ON
  ! on device
  real(nk),dimension(:,:),allocatable:: tmpe,tmpee,tmpo,tmpoo
  ! real*8 tmpe(mye,mz-2),tmpo(myo+1,mz-2)
  ! real*8 tmpee(mye,my-2),tmpoo(myo+1,my-2)
  ! on device
  real(nk),dimension(:,:),allocatable:: tmpye,tmpyo,tmpze,tmpzo
  !real*8 tmpye(mye+1,mz),tmpyo(myo+1,mz)
  !real*8 tmpze(my,mze+1),tmpzo(my,mzo+1)

end module legchebshen

module gpudevice

  external cublas_init, cublas_create, cublas_set_matrix, cublas_get_matrix
  external cublas_set_vector, cublas_get_vector
  external cublas_shutdown, cublas_alloc
  external cublas_dcopy
  external cublas_dgemm, cublas_dtrmm
  external cublas_set_stream, cublas_get_stream
  external cublas_dgeam, transpose_gpu, matmul_gpu, mat_invdiag_gpu
  integer transpose_gpu, matmul_gpu, mat_invdiag_gpu
  !external cublas_stream_create, cublas_set_stream, cublas_get_stream ! no such stream functions

  integer cublas_alloc, cublas_set_matrix, cublas_get_matrix
  integer cublas_set_vector, cublas_get_vector
  integer cublas_set_stream, cublas_get_stream

  external cusparse_get_version, cusparse_create
  external cuda_set_stream
  external cuda_stream_create, cusparse_dgtsv, cusparse_dgtsv_stridedbatch

  integer cusparse_get_version, cusparse_create
  integer cuda_stream_create, cusparse_dgtsv, cusparse_dgtsv_stridedbatch

  integer*8 handle,stream_e,stream_o ! for cusparse

  integer*8,dimension(:),allocatable:: streams

  !integer :: size_of_real=4
  integer, parameter :: size_of_double=8

  !integer stream0, stream1
!#if ARCH_64
!  integer*8 dpAE,dpAO,dpBE,dpBO
  integer*8 dp_PED,dp_POD,dp_PSM1ED,dp_PSM1OD
  integer*8 dp_PEN,dp_PON,dp_PSM1EN,dp_PSM1ON
  integer*8 dp_tmpe,dp_tmpee,dp_tmpo,dp_tmpoo
  integer*8 devPtrF

  integer*8 dp_LCAE, dp_LCAO, dp_LCBE, dp_LCBO
  integer*8 dp_LCA,  dp_LCB
  integer*8 dp_tmpye,dp_tmpyo,dp_tmpze,dp_tmpzo,dp_Eu
  integer*8 devPtrU, devPtrV

  integer*8 devPtrG
  integer*8 dp_gp, dp_vp
  integer*8 dp_fe, dp_fo
  integer*8 dp_SED,dp_MED,dp_MPPED
  integer*8 dp_SOD,dp_MDO,dp_MPPOD ! note that not MDO

  integer*8 dp_ra1,dp_rp

  integer*8 dp_ED, dp_EN
  integer*8 dp_tmp1d

  integer*8 dp_DDE,  dp_DDO
  integer*8 dp_DPE,  dp_DPO
  integer*8 dp_DDEp, dp_DDOp ! tmp
  integer*8 dp_DPEp, dp_DPOp ! tmp
  integer*8 dp_DPEp2,dp_DPOp2 ! tmp
!#else
!  integer dpAE,dpAO,dpBE,dpBO
!  integer dp_PED,dp_POD,dp_PSM1ED,dp_PSM1OD
!  integer dp_tmpe,dp_tmpee,dp_tmpo,dp_tmpoo
!  integer devPtrF
!#endif

  ! for solvediag
  integer*8 devPtrFT
  integer*8 dp_PDM1
  integer*8 dp_PDM1T, dp_ED2, dp_PDT
  integer*8 dp_PNM1
  integer*8 dp_PNM1T, dp_EN2, dp_PNT
  integer*8 dp_tmp2d

end module gpudevice
