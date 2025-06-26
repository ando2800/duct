subroutine pointer1D(iff,ill,np,ibeg,iend)
  implicit none 

  integer iff,ill,np
  integer ibeg(0:np-1),iend(0:np-1)
  integer if,il,mynodes(0:np-1),irank
  integer mx,i,ishift
  !
  mx=ill-iff+1
  ! initialize mynodes(i)
  do i=0,np-1
     mynodes(i)=0
  end do
  ! set number of nodes for each ranks
  do i=1,mx
     irank=mod(i,np)
     mynodes(irank)=mynodes(irank)+1
  end do
  ibeg(0)=1
  iend(0)=mynodes(0)
  do i=1,np-1
     ibeg(i)=iend(i-1)+1
     iend(i)=ibeg(i)+mynodes(i)-1
  end do
  ishift=iff-1
  do i=0,np-1
     ibeg(i)=ibeg(i)+ishift
     iend(i)=iend(i)+ishift
  enddo
  return
end subroutine pointer1D
!----+----------------------------------------------------------------|
