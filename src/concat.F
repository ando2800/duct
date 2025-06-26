c$$$      program toto
c$$$      implicit real*8 (a-h,o-z)
c$$$      character*10 cc
c$$$c
c$$$      call concat('abc',3,'def',3,'ghij',4,cc,10)
c$$$      write(*,*)cc
c$$$      stop
c$$$      end
c---------------------------------------------------------------------
      subroutine concat(a,na,b,nb,c,nc,d)
      implicit real*8 (a-h,o-z)
      character*(*) a
      character*(*) b
      character*(*) c
      character*(*) d
c
      d=a(1:na)//b(1:nb)//c(1:nc)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine concat2(a,na,b,nb,d)
      implicit real*8 (a-h,o-z)
      character*(*) a
      character*(*) b
      character*(*) d
c
      d=a(1:na)//b(1:nb)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine concatb(aa,a,na,b,nb,c,nc,d)
      implicit real*8 (a-h,o-z)
      character*(*) aa
      character*(*) a
      character*(*) b
      character*(*) c
      character*(*) d
c
      ispace=index(aa,' ')
      d=aa(1:ispace-1)//a(1:na)//b(1:nb)//c(1:nc)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine concat3b(aa,a,na,b,nb,e)
      implicit real*8 (a-h,o-z)
      character*(*) aa
      character*(*) a
      character*(*) b
      character*(*) e
c
      ispace=index(aa,' ')
      e=aa(1:ispace-1)//a(1:na)//b(1:nb)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine concat4b(aa,a,na,b,nb,c,nc,e)
      implicit real*8 (a-h,o-z)
      character*(*) aa
      character*(*) a
      character*(*) b
      character*(*) c
      character*(*) e
c
      ispace=index(aa,' ')
      e=aa(1:ispace-1)//a(1:na)//b(1:nb)//c(1:nc)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine concat5b(aa,a,na,b,nb,c,nc,d,nd,e)
      implicit real*8 (a-h,o-z)
      character*(*) aa
      character*(*) a
      character*(*) b
      character*(*) c
      character*(*) d
      character*(*) e
c
      ispace=index(aa,' ')
      e=aa(1:ispace-1)//a(1:na)//b(1:nb)//c(1:nc)//d(1:nd)
c
      return
      end
c-----------------------------------------------------------------------
