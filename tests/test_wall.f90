! wall-row coefficients of A1, A2, B2 and convergence of D2 at the wall up to ny = 1025
program t
  use mod_types; use mod_compact
  implicit none
  integer :: io,ig,ny,k,j
  integer, parameter :: nys(6)=[33,65,129,257,513,1025]
  real(rp) :: ew(6),alpha
  character(len=8) :: gt
  real(rp), allocatable :: f(:),d2(:),ex2(:)
  do ig=1,2; do io=4,6,2
    gt='tanh'; alpha=3._rp; if(ig==2) gt='uniform'
    do k=1,6
      ny=nys(k)
      if(allocated(yg)) deallocate(yg,a1,a2,b1,b2,t1w,t1d,t1u,t2w,t2d,t2u,c1,cn,wq,dyc,b1m,t1mw,t1md,t1mu,rxi)
      call init_grid(ny,gt,alpha,''); call init_compact(ny,io)
      allocate(f(ny),d2(ny),ex2(ny))
      f=sin(2.3_rp*yg+0.4_rp); ex2=-2.3_rp**2*f
      call mat_band_apply(ny,nbw,b2,f,d2); call ts(ny,t2w,t2d,t2u,d2)
      ew(k)=max(abs(d2(1)-ex2(1)),abs(d2(ny)-ex2(ny)))
      if(k==2) then
        print '(a,a,i2)','grid=',trim(gt)//' order=',io
        do j=1,3; print '(a,i1,a,3f10.5,a,f10.5)','  A2 row ',j,':',a2(:,j),'   pivot',1/t2d(j); end do
        do j=1,3; print '(a,i1,a,7f9.3)','  B2*h^2 row ',j,':',b2(:,j)*(yg(2)-yg(1))**2; end do
        do j=1,3; print '(a,i1,a,3f10.5,a,f10.5)','  A1 row ',j,':',a1(:,j),'   pivot',1/t1d(j); end do
      end if
      deallocate(f,d2,ex2)
    end do
    print '(a,6es10.2,a,5f6.2)','  D2 wall err',ew,' ord:',(log(ew(k-1)/ew(k))/log(2._rp),k=2,6)
  end do; end do
contains
  subroutine ts(n,w,d,u,r)
    integer :: n,i; real(rp) :: w(n),d(n),u(n),r(n)
    do i=2,n; r(i)=r(i)-w(i)*r(i-1); end do
    r(n)=r(n)*d(n)
    do i=n-1,1,-1; r(i)=(r(i)-u(i)*r(i+1))*d(i); end do
  end subroutine
end program
