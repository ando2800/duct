subroutine derivyz(u,du,aux,cdir,ider,scal,what)

  use ctes
  use numbers

  implicit none
  !   implicit complex*16
  !   /***************************************************************/
  !   /* computes 1st or 2nd derivatives in the y- or z-chebyshev-   */
  !   /* direction                                                   */
  !   /*                                                             */
  !   /* input                                                       */
  !   /* u(0:my1,0:mz1,ib:ie)   array in c-c-f space                 */
  !   /* aux(0:my1 OR 0:mz1)    1d work array; only accessed if      */
  !   /*                        ider=2; size depends on value of     */
  !   /*                        'cdir'                               */
  !   /* cdir                   'y' or 'z' - direction of deriv      */ 
  !   /* ider                   '1' or '2'     - order of derivative */
  !   /* real*8 scal            REAL scalar factor                   */
  !   /* what                   'n'    - du=   scal*d^n(u)/dxi^n     */
  !   /*                        [else] - du=du+scal*d^n(u)/dxi^n     */
  !   /*                        1<=ider<=2                           */
  !   /*                                                             */
  !   /* output                                                      */
  !   /* du(0:my1,0:mz1,ib:ie)  derivatives in c-c-f space           */
  !   /*                                                             */
  !   /* note: 'u' and 'du' can NOT have the same memory addresses.  */
  !   /*                                                             */
  !   /* note2: the domain is supposed to be defined as: y\in[-1,1]  */
  !   /*        z\in[-aspect,+aspect] --> d/dz=...*1/aspect,         */
  !   /*                                 d^2/dz^2=...*1/aspect^2     */
  !   /***************************************************************/
  !
  complex(nk) u(0:my1,0:mz1,ib:ie),du(0:my1,0:mz1,ib:ie)
  complex(nk) aux(0:*)
  !
  character*1 cdir,what
  real(nk) scal,aspect2
  integer ider, i,j,k, ii,kk,ia
  complex(nk) aari,bri,tri,cri 
  !-----|---------------------------------------------------------------|
  !aspect=z(mz)
  aspect2=aspect**2
  !   /***************overwrite 'du'**********************************/
  if(what.eq.'n')then       !assign d^n(u)/dxi^n
     if(cdir.eq.'y')then
        if(ider.eq.1)then
           do i=ib,ie
              do k=0,mz1
                 aari=zero
                 bri=zero
                 tri=zero
                 ia=my
                 do ii=1,my1
                    kk=my-ii
                    cri = bri + 2.d0*dfloat(kk)*u(ia-1,k,i)
                    du(ia-1,k,i)=-tri*scal
                    bri = aari
                    tri = cri
                    aari = cri
                    ia = ia-1
                 enddo
                 du(ia-1,k,i)=-0.5D0*tri*scal
              enddo
           enddo
        elseif(ider.eq.2)then
           do i=ib,ie
              do k=0,mz1
                 ! /* 1st derivative -> aux */
                 aari=zero
                 bri=zero
                 tri=zero
                 ia=my
                 do ii=1,my1
                    kk=my-ii
                    cri = bri + 2.d0*dfloat(kk)*u(ia-1,k,i)
                    aux(ia-1)=-tri
                    bri = aari
                    tri = cri
                    aari = cri
                    ia = ia-1
                 enddo
                 aux(ia-1)=-0.5D0*tri
                 ! /* 2nd derivative */
                 aari=zero
                 bri=zero
                 tri=zero
                 ia=my
                 do ii=1,my1
                    kk=my-ii
                    cri = bri + 2.d0*dfloat(kk)*aux(ia-1)
                    du(ia-1,k,i)=-tri*scal
                    bri = aari
                    tri = cri
                    aari = cri
                    ia = ia-1
                 enddo
                 du(ia-1,k,i)=-0.5D0*tri*scal
              enddo
           enddo
        else
           write(*,*)'this order of deriv not implemented:',ider
           stop
        endif
     elseif(cdir.eq.'z')then
        if(ider.eq.1)then
           do i=ib,ie
              do j=0,my1
                 aari=zero
                 bri=zero
                 tri=zero
                 ia=mz
                 do ii=1,mz1
                    kk=mz-ii
                    cri = bri + 2.d0*dfloat(kk)*u(j,ia-1,i)
                    du(j,ia-1,i)=-tri*scal/aspect
                    bri = aari
                    tri = cri
                    aari = cri
                    ia = ia-1
                 enddo
                 du(j,ia-1,i)=-0.5D0*tri*scal/aspect
              enddo
           enddo
        elseif(ider.eq.2)then
           do i=ib,ie
              do j=0,my1
                 ! /* 1st derivative -> aux */
                 aari=zero
                 bri=zero
                 tri=zero
                 ia=mz
                 do ii=1,mz1
                    kk=mz-ii
                    cri = bri + 2.d0*dfloat(kk)*u(j,ia-1,i)
                    aux(ia-1)=-tri
                    bri = aari
                    tri = cri
                    aari = cri
                    ia = ia-1
                 enddo
                 aux(ia-1)=-0.5D0*tri
                 ! /* 2nd derivative */
                 aari=zero
                 bri=zero
                 tri=zero
                 ia=mz
                 do ii=1,mz1
                    kk=mz-ii
                    cri = bri + 2.d0*dfloat(kk)*aux(ia-1)
                    du(j,ia-1,i)=-tri*scal/aspect2
                    bri = aari
                    tri = cri
                    aari = cri
                    ia = ia-1
                 enddo
                 du(j,ia-1,i)=-0.5D0*tri*scal/aspect2
              enddo
           enddo
        else
           write(*,*)'this order of deriv not implemented:',ider
           stop
        endif
     else
        write(*,*)'problem in derivyz: direction ',cdir
        stop
     endif
  else   !'what'
     ! /*********************add to previous value of 'du**************/
     if(cdir.eq.'y')then
        if(ider.eq.1)then
           do i=ib,ie
              do k=0,mz1
                 aari=zero
                 bri=zero
                 tri=zero
                 ia=my
                 do ii=1,my1
                    kk=my-ii
                    cri = bri + 2.d0*dfloat(kk)*u(ia-1,k,i)
                    du(ia-1,k,i)=du(ia-1,k,i)-tri*scal
                    bri = aari
                    tri = cri
                    aari = cri
                    ia = ia-1
                 enddo
                 du(ia-1,k,i)=du(ia-1,k,i)-0.5D0*tri*scal
              enddo
           enddo
        elseif(ider.eq.2)then
           do i=ib,ie
              do k=0,mz1
                 ! /* 1st derivative -> aux */
                 aari=zero
                 bri=zero
                 tri=zero
                 ia=my
                 do ii=1,my1
                    kk=my-ii
                    cri = bri + 2.d0*dfloat(kk)*u(ia-1,k,i)
                    aux(ia-1)=-tri
                    bri = aari
                    tri = cri
                    aari = cri
                    ia = ia-1
                 enddo
                 aux(ia-1)=-0.5D0*tri
                 ! /* 2nd derivative */
                 aari=zero
                 bri=zero
                 tri=zero
                 ia=my
                 do ii=1,my1
                    kk=my-ii
                    cri = bri + 2.d0*dfloat(kk)*aux(ia-1)
                    du(ia-1,k,i)=du(ia-1,k,i)-tri*scal
                    bri = aari
                    tri = cri
                    aari = cri
                    ia = ia-1
                 enddo
                 du(ia-1,k,i)=du(ia-1,k,i)-0.5D0*tri*scal
              enddo
           enddo
        else
           write(*,*)'this order of deriv not implemented:',ider
           stop
        endif
     elseif(cdir.eq.'z')then
        if(ider.eq.1)then
           do i=ib,ie
              do j=0,my1
                 aari=zero
                 bri=zero
                 tri=zero
                 ia=mz
                 do ii=1,mz1
                    kk=mz-ii
                    cri = bri + 2.d0*dfloat(kk)*u(j,ia-1,i)
                    du(j,ia-1,i)=du(j,ia-1,i)-tri*scal/aspect
                    bri = aari
                    tri = cri
                    aari = cri
                    ia = ia-1
                 enddo
                 du(j,ia-1,i)=du(j,ia-1,i)-0.5D0*tri*scal/aspect
              enddo
           enddo
        elseif(ider.eq.2)then
           do i=ib,ie
              do j=0,my1
                 ! /* 1st derivative -> aux */
                 aari=zero
                 bri=zero
                 tri=zero
                 ia=mz
                 do ii=1,mz1
                    kk=mz-ii
                    cri = bri + 2.d0*dfloat(kk)*u(j,ia-1,i)
                    aux(ia-1)=-tri
                    bri = aari
                    tri = cri
                    aari = cri
                    ia = ia-1
                 enddo
                 aux(ia-1)=-0.5D0*tri
                 ! /* 2nd derivative */
                 aari=zero
                 bri=zero
                 tri=zero
                 ia=mz
                 do ii=1,mz1
                    kk=mz-ii
                    cri = bri + 2.d0*dfloat(kk)*aux(ia-1)
                    du(j,ia-1,i)=du(j,ia-1,i)-tri*scal/aspect2
                    bri = aari
                    tri = cri
                    aari = cri
                    ia = ia-1
                 enddo
                 du(j,ia-1,i)=du(j,ia-1,i)-0.5D0*tri*scal/aspect2
              enddo
           enddo
        else
           write(*,*)'this order of deriv not implemented:',ider
           stop
        endif
     else
        write(*,*)'problem in derivyz: direction ',cdir
        stop
     endif
  endif                     !'what'
  
  return
end subroutine derivyz
!-----|---------------------------------------------------------------|
subroutine derivyz_plane(uy,uz,duy,duz,cdir,cpdir,scal)

  use ctes
  use numbers

  implicit none
  !  /***************************************************************/
  !  /* computes 1st derivatives in the y- or z-chebyshev-dir in a  */
  !  /* plane only.                                                 */
  !  /*                                                             */
  !  /* input                                                       */
  !  /* uy,uz    are planes of a full array in fourier-c-c space    */
  !  /*                        [u(0:my1,0:mz1,0:mx1)]               */
  !  /*          where actually either 'uy' or 'uz' are accessed    */
  !  /*          depending on the value of 'cdir'. the other is then*/
  !  /*          a dummy (for reasons of dimensioning). the second  */
  !  /*          dimension of 'uy','uz' depends upon the value of   */
  !  /*          'cpdir' (indicating the 2nd direction of the res-  */
  !  /*          pective plane in which we operate).                */
  !  /*    dimension u(0:my1,0:mz1) -> cdir='y' & cpdir='z'         */
  !  /*                             OR cdir='z' & cpdir='y'         */
  !  /*    dimension u(0:my1,0:mx1) -> cdir='y' & cpdir='x'         */
  !  /*    dimension u(0:mz1,0:mx1) -> cdir='z' & cpdir='x'         */
  !  /* cdir                   'y' or 'z' - direction of deriv      */ 
  !  /* cpdir                  'y' or 'z' or 'x' - 2nd dir of plane */
  !  /*            note: only the following combinations are legal  */
  !  /*                cdir='y' & cpdir='z'                         */
  !  /*                cdir='z' & cpdir='y'                         */
  !  /*                cdir='y' & cpdir='x'                         */
  !  /*                cdir='z' & cpdir='x'                         */
  !  /* real*8 scal            REAL scalar factor                   */
  !  /*                                                             */
  !  /* output                                                      */
  !  /* duy,duz     the resulting 1st derivatives in the plane.     */
  !  /*             the same comments as for 'uy','uz' apply wrt    */
  !  /*             dimensioning and dummy variables.               */
  !  /*                                                             */
  !  /* note: 'uy'/'duy' and 'uz'/'duz' can NOT have the same       */
  !  /*            memory addresses.                                */
  !  /*                                                             */
  !  /* note2: the domain is supposed to be defined as: y\in[-1,1]  */
  !  /*        z\in[-aspect,+aspect] --> d/dz=...*1/aspect,         */
  !  /*                                 d^2/dz^2=...*1/aspect^2     */
  !  /***************************************************************/
  complex(nk) uy(0:my1,0:*),duy(0:my1,0:*)
  complex(nk) uz(0:mz1,0:*),duz(0:mz1,0:*)
  character*1 cdir,cpdir
  real(nk) scal,aspect2
  
  integer ia,i,j,k,jf,jl,kf,kl,ii,kk
  complex(nk) aari,bri,tri,cri
     
  !-----|---------------------------------------------------------------|
  write(*,*)'derivyz_plane: not yet tested, stop!'
  stop
  !aspect=z(mz)
  
  if(cdir.eq.'y')then       !derivative in y
     if(cpdir.eq.'x')then   !plane (y,x)
        kf=0
        kl=mx1
     elseif(cpdir.eq.'z')then !plane (y,z)
        kf=0
        kl=mz1
     else
        write(*,*)'derivyz_plane: ',cdir,' ',cpdir
        stop
     endif
     do k=kf,kl
        aari=zero
        bri=zero
        tri=zero
        ia=my
        do ii=1,my1
           kk=my-ii
           cri = bri + 2.d0*dfloat(kk)*uy(ia-1,k)
           duy(ia-1,k)=-tri*scal
           bri = aari
           tri = cri
           aari = cri
           ia = ia-1
        enddo
        duy(ia-1,k)=-0.5D0*tri*scal
     enddo
  elseif(cdir.eq.'z')then   !derivative in z
     if(cpdir.eq.'x')then   !plane (z,x)
        jf=0
        jl=mx1
        do j=jf,jl
           aari=zero
           bri=zero
           tri=zero
           ia=mz
           do ii=1,mz1
              kk=mz-ii
              cri = bri + 2.d0*dfloat(kk)*uz(ia-1,j)
              duz(ia-1,j)=-tri*scal/aspect
              bri = aari
              tri = cri
              aari = cri
              ia = ia-1
           enddo
           duz(ia-1,j)=-0.5D0*tri*scal/aspect
        enddo
     elseif(cpdir.eq.'y')then !plane (y,z)
        jf=0
        jl=my1
        do j=jf,jl
           aari=zero
           bri=zero
           tri=zero
           ia=mz
           do ii=1,mz1
              kk=mz-ii
              cri = bri + 2.d0*dfloat(kk)*uy(j,ia-1)
              duy(j,ia-1)=-tri*scal/aspect
              bri = aari
              tri = cri
              aari = cri
              ia = ia-1
           enddo
           duy(j,ia-1)=-0.5D0*tri*scal/aspect
        enddo
     else
        write(*,*)'derivyz_plane: ',cdir,' ',cpdir
        stop
     endif
  endif
  
  return
end subroutine derivyz_plane
!-----|---------------------------------------------------------------|
subroutine dglc(A,B,X,NP)
  ! -------------------------------------------------------------
  !          COMPUTES CHEBYSHEV NODES IN THE INTERVAL A-B
  !     A,B: end points of the interval                                   
  !     X: vector of Gauss-Lobatto nodes (dimension NP)                   
  !     NP: number of points on the interval                              
  ! -------------------------------------------------------------
  use numbers
  implicit none

  real*8 ALF,PIN,A,B
  real*8 X(NP)
  integer N,NP,I

  PI=4.0*DATAN(1.D0)
  ALF=2.0/(B-A)
  N=NP-1
  PIN=PI/dfloat(N)
  do I=1,NP
     X(I)=-DCOS((I-1)*PIN)
  end do
  do I=1,NP
     X(I)=((B-A)*X(I)+A+B)*.5D0
  end do

  return

end subroutine dglc
!-----|---------------------------------------------------------------|
subroutine transyz(cheb,cdir,iopt)

  use ctes, only: my, mz, mx1, ib, ie
  use tcheb2tmp

  implicit none
  !  /***************************************************************/
  !  /* chebyshev transform of a 3d field in either 'y' or 'z'      */
  !  /* direction                                                   */
  !  /*                                                             */
  !  /* input                                                       */
  !  /* real*8 cheb(nsize) field in either one: f-p-p,f-c-p,f-p-c,  */
  !  /*                     f-c-c space; nsize=2*my*mz*(mx1+1)      */
  !  /* cdir               'y' or 'z' - direction of transform      */
  !  /* iopt               >=0 cheb-->phys                          */
  !  /*                    < 0 phys-->cheb                          */
  !  /*                                                             */
  !  /* output                                                      */
  !  /* cheb(nsize)        transformed field                        */
  !  /*                                                             */
  !  /***************************************************************/
  !  implicit double precision (a-h,o-z)
  !include"fftw3.f"
  !real*8 cheb(2*my*mz*(mx1+1)) ! double ! used be c(*)
  !real*8 cheb(2*my*mz*(ie-ib+1)) ! double ! used be c(*)
  real*8 cheb(*) ! double ! used be c(*)
  character*1 cdir
  integer nstrid,nskip,nn,ifst,iopt
  integer i
!-----|---------------------------------------------------------------|
!     
  if(cdir.eq.'y')then       !transform in 'y'
     nstrid=2               !real & imaginary
     nskip=my*2
     !nstrid=1               !real & imaginary for tcheb2
     !nskip=my
     !nn=mz*(mx1+1) !BUG
     nn=mz*(ie-ib+1) 
     ifst=1     
     if (iopt.lt.0) then
        ! /* physical to chebi. transforms  */
        call tcheb_fftw(cheb(ifst  ),nstrid,nskip,nn,'y',-1)
        call tcheb_fftw(cheb(ifst+1),nstrid,nskip,nn,'y',-1)
        !call tcheb2_fftw(cheb(ifst),nstrid,nskip,nn,'y',-1) ! only for f-c-c 
     elseif (iopt.gt.0) then   
        ! /* chebi to physical transforms  */ 
        call tcheb_fftw(cheb(ifst  ),nstrid,nskip,nn,'y',+1)
        call tcheb_fftw(cheb(ifst+1),nstrid,nskip,nn,'y',+1)
        !call tcheb2_fftw(cheb(ifst),nstrid,nskip,nn,'y',+1) ! only for f-p-p 
     else
        write(*,*) 'transyz iopt error'
        stop
     endif          
  elseif (cdir.eq.'z')then  !transform in 'z' ! very slow for a large my
     nstrid=2*my ! for tcheb
     nskip=2
     !nstrid=my  ! for tcheb2
     !nskip=1
     nn=my
     do i=1,(ie-ib+1)         !fourier modes, real & imag part
        ifst=(i-1)*my*mz*2+1     
        if (iopt.lt.0) then
           ! /* physical to chebi. transforms  */
           call tcheb_fftw(cheb(ifst  ),nstrid,nskip,nn,'z',-1)
           call tcheb_fftw(cheb(ifst+1),nstrid,nskip,nn,'z',-1)
           !call tcheb2_fftw(cheb(ifst),nstrid,nskip,nn,'z',-1)
        elseif (iopt.gt.0) then        
           ! /* chebi to physical transforms  */         
           call tcheb_fftw(cheb(ifst  ),nstrid,nskip,nn,'z',+1)
           call tcheb_fftw(cheb(ifst+1),nstrid,nskip,nn,'z',+1)
           !call tcheb2_fftw(cheb(ifst),nstrid,nskip,nn,'z',+1)
        else
           write(*,*) 'transyz iopt error'
           stop  
        endif
     enddo
  else
     write(*,*) 'err in transyz, stop'
     stop
  endif
     
  return
end subroutine transyz
!-----|---------------------------------------------------------------|
subroutine transyz_plane(cheb,iopt)

  use ctes, only: my, mz
  use tcheb2tmp

  implicit none
  ! /***************************************************************/
  ! /* chebyshev transform of a 2d field in BOTH   'y' or 'z'      */
  ! /* direction                                                   */
  ! /*                                                             */
  ! /* input                                                       */
  ! /* real*8 cheb(nsize) field in either one: p-p-p,or p-c-c      */
  ! /*       modified comment by sekimoto                          */
  ! /*    ATT: input must be 2d physical space                     */
  ! /* %%real*8 cheb(nsize) field in either one: f-p-p,or f-c-c    */
  ! /*                         nsize=my*mz                         */
  ! /* iopt               >=0 cheb-->phys                          */
  ! /*                    < 0 phys-->cheb                          */
  ! /*                                                             */
  ! /* output                                                      */
  ! /* cheb(nsize)        transformed field                        */
  ! /*                                                             */
  ! /***************************************************************/
  !real*8 cheb(2*my*mz)
  real*8 cheb(*)
  integer iopt,nn,nstrid,nskip
  !-----|---------------------------------------------------------------|
  if (iopt.lt.0) then       !physical to chebyshev
     ! /* transform in 'y' */
     nn=mz
     nstrid=1
     nskip=my
     call tcheb_fftw(cheb,nstrid,nskip,nn,'y',-1)
     ! /* transform in 'z' */
     nn=my
     nstrid=my
     nskip=1
     call tcheb_fftw(cheb,nstrid,nskip,nn,'z',-1)
  else
     ! /* chebi to physical transforms  */         
     ! /* transform in 'y' */
     nn=mz
     nstrid=1
     nskip=my
     call tcheb_fftw(cheb,nstrid,nskip,nn,'y',+1)
     ! /* transform in 'z' */
     nn=my
     nstrid=my
     nskip=1
     call tcheb_fftw(cheb,nstrid,nskip,nn,'z',+1)
  endif
  
  return
end subroutine transyz_plane
!----|----------------------------------------------------------------|
subroutine cheb2_wallstress(u,nstri,aspect,sum,ny,nz)

  use numbers

  implicit none
  !      implicit real*8(a-h,o-z)
  !
  !   /* computes the wall-stress integrated over y,z for a field in */
  !   /* in 2d chebyshev space                                       */
  !   /*                                                             */
  !   /* input                                                       */
  !   /* u(nstri,0:ny-1,0:nz-1)   field in chebyshev space           */
  !   /* aspect                   aspect ratio: max(z)/max(y)        */
  !   /*                                                             */
  !   /* output                                                      */
  !   /* sum                      wall-stress integral               */
  !
  integer nstri,ny,nz
  real*8 u(nstri,0:ny-1,0:nz-1)
  real*8 aspect,sum
  integer k,j
  !----|----------------------------------------------------------------|
  sum=f0o1
  do k=0,nz-1,2
     do j=0,ny-1
        sum=sum+u(1,j,k) &
             &  *(-(f1o1-(-f1o1)**(j+1)) &
             &    *f2o1*aspect*dfloat(j**2)/dfloat(k**2-1) & 
             &   )
     enddo
  enddo
  do k=0,nz-1
     do j=0,ny-1,2
        sum=sum+u(1,j,k) &
             &  *(-(f1o1-(-f1o1)**(k+1)) &
             &    *f2o1/aspect*dfloat(k**2)/dfloat(j**2-1) &
             &   ) !formula (a.6) of report is only valid for even functions!!!!
     enddo
  enddo
  !     /* finally: divide the integral by the circumference */
  sum=sum/(f4o1+f4o1*aspect)

  return
end subroutine cheb2_wallstress
!----|----------------------------------------------------------------|
subroutine cheb2_wallstress_circum(u,nstri,aspect,sum,ny,nz)

  use numbers

  implicit none
!      implicit real*8(a-h,o-z)
!
!     /* computes the wall-stress as fct of circumfer.for a field in */
!     /* in 2d chebyshev space                                       */
!     /*                                                             */
!     /* input                                                       */
!     /* u(nstri,0:ny-1,0:nz-1)   field in chebyshev space           */
!     /* aspect                   aspect ratio: max(z)/max(y)        */
!     /*                                                             */
!     /* output                                                      */
!     /* sum(4)                   wall-stress summed over each of the*/
!     /*                          four walls                         */
!
  integer nstri,ny,nz
  real*8 u(nstri,0:ny-1,0:nz-1),sum(4)
  integer i,j,k
  real*8 aspect
  !----|----------------------------------------------------------------|
  call dcopy(4,f0o1,0,sum,1)
  do k=0,nz-1,2
     do j=0,ny-1
        sum(1)=sum(1)-u(1,j,k) &
             &        *f2o1*aspect*dfloat(j**2)/dfloat(k**2-1)
        sum(2)=sum(2)+u(1,j,k) & 
             &        *(-f1o1)**(j+1) &
             &        *f2o1*aspect*dfloat(j**2)/dfloat(k**2-1)
     enddo
  enddo
  do k=0,nz-1
     do j=0,ny-1,2
        sum(3)=sum(3)-u(1,j,k) & 
             &        *f2o1/aspect*dfloat(k**2)/dfloat(j**2-1)
        sum(4)=sum(4)+u(1,j,k) &
             &        *(-f1o1)**(k+1) &
             &        *f2o1/aspect*dfloat(k**2)/dfloat(j**2-1)
     enddo
  enddo
  ! /* finally: divide the integral by the circumference of the walls */
  do i=1,4
     sum(i)=sum(i)/(f4o1+f4o1*aspect)
  enddo
  
  return
end subroutine cheb2_wallstress_circum
!----|----------------------------------------------------------------|
subroutine cheb2_integral(u,nstri,aspect,sum,ny,nz)
  
  use numbers

  implicit none
  !  implicit real*8(a-h,o-z)
  !
  !  /* computes the integral over y,z for a field                  */
  !  /* in 2d chebyshev space                                       */
  !  /*                                                             */
  !  /* input                                                       */
  !  /* u(nstri,0:ny-1,0:nz-1)   field in chebyshev space           */
  !  /* aspect                   aspect ratio: max(z)/max(y)        */
  !  /*                                                             */
  !  /* output                                                      */
  !  /* sum                      integral of u over y,z             */
  !
  integer nstri,ny,nz
  real*8 u(nstri,0:ny-1,0:nz-1)
  real*8 aspect,sum
  integer j,i
  !----|----------------------------------------------------------------|
  sum=f0o1
  do j=0,nz-1,2
     do i=0,ny-1,2
        sum=sum+u(1,i,j) & 
             &  *f4o1*aspect/(dfloat(j**2-1)*dfloat(i**2-1))
     enddo
  enddo
  
  return
end subroutine cheb2_integral
!----|----------------------------------------------------------------|
subroutine cheb_interp(myr,mzr,tmp1,my,mz,tmpo1)
  use ctes,only:nk
  use numbers
  implicit none 

  ! interpolation (myr,myz) --> (my,mz) in f-p-p space
  integer myr,mzr,my,mz
  real(nk) tmp1(2,myr,mzr) 
  real(nk) tmpo1(2,my,mz) 
  real(nk) y2(my),z2(mz)
  integer i,j,k,l,kk,ll

  call tchebi_fftw(myr,mzr)

  call transyz_plane_fou(tmp1,myr,mzr,-1)

  call tchebi_fftw(my,mz)
  ! ---------------------------------------------------------------|
  ! /* compute the TARGET y- and z-coordinates */
  call dglc(-f1o1,f1o1,y2,my)
  call dglc(-f1o1,f1o1,z2,mz)
  y2=-y2;  z2=-z2

  do i=1,my                !index of y in p-space
     do j=1,mz             !index of z in p-space
        tmpo1(:,i,j)=f0o1
        do k=1,min(20,myr)           !index of y in c-space
           kk=k-1
           do l=1,min(20,mzr)        !index of z in c-space
              ll=l-1
              if ((dabs(tmp1(1,k,l))+dabs(tmp1(2,k,l))).lt.1.e-7) continue
              tmpo1(1:2,i,j)=tmpo1(1:2,i,j)+tmp1(1:2,k,l) &
                   &               *dcos(dfloat(kk)*dacos(y2(i))) &
                   &               *dcos(dfloat(ll)*dacos(z2(j)))
           enddo
        enddo
     enddo
  enddo
  
  return
end subroutine cheb_interp
      
subroutine transyz_plane_fou(cheb,my,mz,iopt)
  !use ctes, only: my, mz
  use tcheb2tmp

  implicit none
  ! /***************************************************************/
  ! /* chebyshev transform of a 2d field in BOTH   'y' or 'z'      */
  ! /* direction                                                   */
  ! /*                                                             */
  ! /* input                                                       */
  ! /* real*8 cheb(nsize) field in either one: f-p-p,or f-c-c      */
  ! /*       modified comment by sekimoto                          */
  ! /*    ATT: input must be 2d physical space                     */
  ! /*                         nsize=my*mz                         */
  ! /* iopt               >=0 cheb-->phys                          */
  ! /*                    < 0 phys-->cheb                          */
  ! /*                                                             */
  ! /* output                                                      */
  ! /* cheb(nsize)        transformed field                        */
  ! /*                                                             */
  ! /***************************************************************/
  integer my,mz,iopt
  real*8 cheb(2*my*mz)
  !real*8 cheb(*)
  integer nn,nstrid,nskip
  !-----|---------------------------------------------------------------|
  if (iopt.lt.0) then       !physical to chebyshev
     ! /* transform in 'y' */
     nn=mz
     nstrid=2
     nskip=my*2
     call tcheb_fftw(cheb(1),nstrid,nskip,nn,'y',-1)
     call tcheb_fftw(cheb(2),nstrid,nskip,nn,'y',-1)
     ! /* transform in 'z' */
     nn=my
     nstrid=2*my
     nskip=2
     call tcheb_fftw(cheb(1),nstrid,nskip,nn,'z',-1)
     call tcheb_fftw(cheb(2),nstrid,nskip,nn,'z',-1)
  else
     ! /* chebi to physical transforms  */         
     ! /* transform in 'y' */
     nn=mz
     nstrid=2
     nskip=my*2
     call tcheb_fftw(cheb(1),nstrid,nskip,nn,'y',+1)
     call tcheb_fftw(cheb(2),nstrid,nskip,nn,'y',+1)
     ! /* transform in 'z' */
     nn=my
     nstrid=2*my
     nskip=2
     call tcheb_fftw(cheb(1),nstrid,nskip,nn,'z',+1)
     call tcheb_fftw(cheb(2),nstrid,nskip,nn,'z',+1)
  endif
  
  return
end subroutine transyz_plane_fou
