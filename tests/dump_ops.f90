! writes hichan's dense D1 = A1^{-1} B1 and D1m = diag(rxi) A1u^{-1} B1u (4th order), the grid and the
! quadrature weights to ops_<grid>_<ny>.bin (stream, little-endian): ny, y(ny), wq(ny), D1(ny,ny), D1m(ny,ny)
program dump_ops
  use mod_types; use mod_compact
  implicit none
  integer, parameter :: nys(5)=[33,65,129,257,513]
  character(len=8), parameter :: names(4)=[character(len=8) :: 'uniform','tanh2','tanh3','tanh4']
  real(rp), parameter :: alphas(4)=[0._rp,2._rp,3._rp,4._rp]
  integer :: ig,k,ny,j,iu
  real(rp), allocatable :: d1(:,:),d1m(:,:),e(:)
  character(len=64) :: fn
  do ig=1,4
    do k=1,5
      ny=nys(k)
      if(allocated(yg)) deallocate(yg,a1,a2,b1,b2,t1w,t1d,t1u,t2w,t2d,t2u,c1,cn,wq,dyc,b1m,t1mw,t1md,t1mu,rxi)
      if(ig==1) then
        call init_grid(ny,'uniform',0._rp,'')
      else
        call init_grid(ny,'tanh',alphas(ig),'')
      end if
      call init_compact(ny,4)
      allocate(d1(ny,ny),d1m(ny,ny),e(ny))
      do j=1,ny
        e=0; e(j)=1
        call mat_band_apply(ny,nbw,b1,e,d1(:,j)); call ts(ny,t1w,t1d,t1u,d1(:,j))
        call mat_band_apply(ny,nbw,b1m,e,d1m(:,j)); call ts(ny,t1mw,t1md,t1mu,d1m(:,j)); d1m(:,j)=d1m(:,j)*rxi
      end do
      write(fn,'(a,a,a,i0,a)') 'ops_',trim(names(ig)),'_',ny,'.bin'
      open(newunit=iu,file=fn,access='stream',form='unformatted',status='replace')
      write(iu) ny,yg,wq,d1,d1m
      close(iu)
      deallocate(d1,d1m,e)
    end do
  end do
contains
  subroutine ts(n,w,d,u,r)
    integer :: n,i; real(rp) :: w(n),d(n),u(n),r(n)
    do i=2,n; r(i)=r(i)-w(i)*r(i-1); end do
    r(n)=r(n)*d(n)
    do i=n-1,1,-1; r(i)=(r(i)-u(i)*r(i+1))*d(i); end do
  end subroutine
end program
