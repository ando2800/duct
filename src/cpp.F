      subroutine cpp(my,mz)
      implicit none 
      integer my,mz
c
      write(*,*)'active CPP options are the following:'
#ifdef FILTER
      write(*,*)'. FILTER'
#endif
#ifdef DEBUG
      write(*,*)'. DEBUG'
#endif
#ifdef TIMING
      write(*,*)'. TIMING'
#endif
#ifdef DT_FIXED
      write(*,*)'. DT_FIXED'
#endif
#ifdef CMASSFLOW
      write(*,*)'. CMASSFLOW'
#if !defined(DT_FIXED)
      write(*,*)'CMASSFLOW only implemented with DT_FIXED!'
      stop
#endif
#if defined(TEMPERATURE)
      write(*,*)'gravx must be zero to use CMASSFLOW for TEMPERATURE!'
#endif
#else /* CMASSFLOW*/
#ifdef HELM_SHEN
      write(*,*)'constant pressure gradient may not work with SHEN'
!      stop
#endif
#endif /* CMASSFLOW*/
#ifdef TEMPERATURE
      write(*,*)'. TEMPERATURE'
#endif
#ifdef RECTANGULAR
      write(*,*)'. RECTANGULAR'
#endif

#ifdef G77
      write(*,*)'. G77'
#endif

#ifdef PROFILO
      write(*,*)'. PROFILO'
#endif

#ifdef BLOCKING_COMM
      write(*,*)'. BLOCKING_COMM'
#endif
#ifdef HELM_SHEN
      write(*,*)'. HELM_SHEN'
#endif
#ifdef FFTW
      write(*,*)'. FFTW3 for all directions'
#endif
#ifdef SET_NOSLIP
      write(*,*)'. SET_NOSLIP'
#endif
#ifdef FORCING_GLOBAL
      write(*,*)'. FORCING_GLOBAL'
#endif
#ifdef PUFF_FORCING
      write(*,*)'. PUFF_FORCING'
#endif
#ifdef NBPOISSON
      write(*,*)'. NBPOISSON'
#endif /* NBPOISSON */

#ifdef CUBLAS
      write(*,*)'. CUBLAS: use gpu device'
#else 
      write(*,*)'. CUBLAS is off: no gpu mode '
#endif /* CUBLAS*/
c
      return
      end subroutine cpp
