! global accuracy of D1, D2 (max and wall error), D1m, quadrature and the 6-point wall derivative
program t
  use mod_types; use mod_compact
  implicit none
  integer :: io,ig,ny,k
  integer, parameter :: nys(4)=[33,65,129,257]
  real(rp) :: e1(4),e2(4),em(4),eq(4),eb1(4),eb2(4),ew(4),alpha
  character(len=8) :: gt
  real(rp), allocatable :: f(:),d(:),d2(:),ex1(:),ex2(:)
  do ig=1,2
  do io=4,6,2
    gt='tanh'; alpha=3._rp; if(ig==2) gt='uniform'
    do k=1,4
      ny=nys(k)
      if(allocated(yg)) deallocate(yg,a1,a2,b1,b2,t1w,t1d,t1u,t2w,t2d,t2u,c1,cn,wq,dyc,b1m,t1mw,t1md,t1mu,rxi)
      call init_grid(ny,gt,alpha,'')
      call init_compact(ny,io)
      allocate(f(ny),d(ny),d2(ny),ex1(ny),ex2(ny))
      f=sin(2.3_rp*yg+0.4_rp); ex1=2.3_rp*cos(2.3_rp*yg+0.4_rp); ex2=-2.3_rp**2*f
      call mat_band_apply(ny,nbw,b1,f,d); call ts(ny,t1w,t1d,t1u,d)
      call mat_band_apply(ny,nbw,b2,f,d2); call ts(ny,t2w,t2d,t2u,d2)
      e1(k)=maxval(abs(d-ex1)); e2(k)=maxval(abs(d2-ex2))
      eb1(k)=abs(d(1)-ex1(1)); eb2(k)=abs(d2(1)-ex2(1))
      ! explicit 6-point wall derivative
      ew(k)=abs(sum(ew1*f(1:6))-ex1(1))
      call mat_band_apply(ny,nbw,b1m,f,d); call ts(ny,t1mw,t1md,t1mu,d); d=d*rxi
      em(k)=maxval(abs(d-ex1))
      eq(k)=abs(sum(wq*f)-(-(cos(2.3_rp+0.4_rp)-cos(-2.3_rp+0.4_rp))/2.3_rp))
      if(k==4) print '(a,es10.2,a,i4)','  lam2max =',lam2max,'   nc1 =',nc1
      deallocate(f,d,d2,ex1,ex2)
    end do
    print '(a,a,a,i1)', 'grid=',trim(gt),' order=',io
    call rep('D1 max   ',e1); call rep('D1 wall  ',eb1); call rep('D2 max   ',e2); call rep('D2 wall  ',eb2)
    call rep('D1m max  ',em); call rep('quad     ',eq); call rep('ew1 6pt  ',ew)
  end do
  end do
contains
  subroutine rep(s,e)
    character(len=*) :: s; real(rp) :: e(4); integer :: k
    print '(a,4es10.2,a,3f6.2)', s,e,'  orders:',(log(e(k-1)/e(k))/log(2._rp),k=2,4)
  end subroutine
  subroutine ts(n,w,d,u,r)
    integer :: n,i; real(rp) :: w(n),d(n),u(n),r(n)
    do i=2,n; r(i)=r(i)-w(i)*r(i-1); end do
    r(n)=r(n)*d(n)
    do i=n-1,1,-1; r(i)=(r(i)-u(i)*r(i+1))*d(i); end do
  end subroutine
end program
