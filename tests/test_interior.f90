! interior error (|y|<0.5), local truncation error of rows 3..ny-2, and the dense row c1 vs D1(1)
program t
  use mod_types; use mod_compact
  implicit none
  integer :: io,ig,ny,k,jm(1)
  integer, parameter :: nys(5)=[33,65,129,257,513]
  real(rp) :: e1(5),e2(5),ei1(5),ei2(5),ec(5),alpha
  character(len=8) :: gt
  real(rp), allocatable :: f(:),d(:),d2(:),ex1(:),ex2(:),r1(:),r2(:)
  do ig=1,2
  do io=4,6,2
    gt='tanh'; alpha=3._rp; if(ig==2) gt='uniform'
    do k=1,5
      ny=nys(k)
      if(allocated(yg)) deallocate(yg,a1,a2,b1,b2,t1w,t1d,t1u,t2w,t2d,t2u,c1,cn,wq,dyc,b1m,t1mw,t1md,t1mu,rxi)
      call init_grid(ny,gt,alpha,'')
      call init_compact(ny,io)
      allocate(f(ny),d(ny),d2(ny),ex1(ny),ex2(ny),r1(ny),r2(ny))
      f=sin(2.3_rp*yg+0.4_rp); ex1=2.3_rp*cos(2.3_rp*yg+0.4_rp); ex2=-2.3_rp**2*f
      ! local truncation error: residual of A f^(d) - B f with exact derivatives
      call tri(ny,a1,ex1,r1); call mat_band_apply(ny,nbw,b1,f,d); r1=r1-d
      call tri(ny,a2,ex2,r2); call mat_band_apply(ny,nbw,b2,f,d); r2=r2-d
      call mat_band_apply(ny,nbw,b1,f,d); call ts(ny,t1w,t1d,t1u,d)
      call mat_band_apply(ny,nbw,b2,f,d2); call ts(ny,t2w,t2d,t2u,d2)
      ec(k)=abs(sum(c1(1:nc1)*f(1:nc1))-d(1))
      ! error away from the walls (middle half of the channel)
      ei1(k)=maxval(abs(d-ex1),mask=abs(yg)<0.5_rp); ei2(k)=maxval(abs(d2-ex2),mask=abs(yg)<0.5_rp)
      e1(k)=maxval(abs(r1(3:ny-2))); e2(k)=maxval(abs(r2(3:ny-2)))
      if(k==5) then
        jm=maxloc(abs(d2-ex2)); print '(a,i4,a,f8.4)','  D2 max error at j=',jm(1),' y=',yg(jm(1))
        print '(a,3es10.2)','  local trunc. D2 rows 1,2,3: ',abs(r2(1:3))
      end if
      deallocate(f,d,d2,ex1,ex2,r1,r2)
    end do
    print '(a,a,a,i1)', 'grid=',trim(gt),' order=',io
    call rep('D1 |y|<.5 ',ei1); call rep('D2 |y|<.5 ',ei2)
    call rep('LTE1 3..n-2',e1); call rep('LTE2 3..n-2',e2); call rep('c1 vs D1(1)',ec)
  end do
  end do
contains
  subroutine rep(s,e)
    character(len=*) :: s; real(rp) :: e(5); integer :: k
    print '(a,5es10.2,a,4f6.2)', s,e,' ord:',(log(e(k-1)/e(k))/log(2._rp),k=2,5)
  end subroutine
  subroutine tri(n,a,x,y)
    integer :: n,i; real(rp) :: a(-1:1,n),x(n),y(n)
    do i=1,n; y(i)=a(0,i)*x(i); if(i>1) y(i)=y(i)+a(-1,i)*x(i-1); if(i<n) y(i)=y(i)+a(1,i)*x(i+1); end do
  end subroutine
  subroutine ts(n,w,d,u,r)
    integer :: n,i; real(rp) :: w(n),d(n),u(n),r(n)
    do i=2,n; r(i)=r(i)-w(i)*r(i-1); end do
    r(n)=r(n)*d(n)
    do i=n-1,1,-1; r(i)=(r(i)-u(i)*r(i+1))*d(i); end do
  end subroutine
end program
