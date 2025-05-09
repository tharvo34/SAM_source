subroutine diagnose
	
! Diagnose some useful stuff

use vars
use params
use sgs, only: sgs_diagnose
implicit none
	
integer i,j,k,kb,kc,k200,k500,k850
real(8) coef, coef1, buffer(nzm,10), buffer1(nzm,10)
real omn, omp, tmp_lwp

coef = 1./float(nx*ny)

call t_startf ('diagnose')

k200 = nzm
	
do k=1,nzm
  u0(k)=0.
  v0(k)=0.
  t01(k) = tabs0(k)
  q01(k) = q0(k)
  t0(k)=0.
  tabs0(k)=0.
  q0(k)=0.
  qcldliq0(k)=0.
  qpcpliq0(k)=0.
  qcldice0(k)=0.
  qpcpice0(k)=0.
  p0(k)=0.
  kc=min(nzm,k+1)
  kb=max(1,k-1)
  if(pres(kc).le.200..and.pres(kb).gt.200.) k200=k
  coef1 = rho(k)*dz*adz(k)*dtfactor
  do j=1,ny
    do i=1,nx
     tabs(i,j,k) = t(i,j,k)-gamaz(k)+ fac_cond * (qcl(i,j,k)+qpl(i,j,k)) +&
	                              fac_sub *(qci(i,j,k) + qpi(i,j,k))
     u0(k)=u0(k)+u(i,j,k)
     v0(k)=v0(k)+v(i,j,k)
     p0(k)=p0(k)+p(i,j,k)
     t0(k)=t0(k)+t(i,j,k)
     tabs0(k)=tabs0(k)+tabs(i,j,k)
     q0(k)=q0(k)+qv(i,j,k)+qcl(i,j,k)+qci(i,j,k)
     qcldliq0(k) = qcldliq0(k) + qcl(i,j,k)
     qpcpliq0(k) = qpcpliq0(k) + qpl(i,j,k)
     qcldice0(k) = qcldice0(k) + qci(i,j,k)
     qpcpice0(k) = qpcpice0(k) + qpi(i,j,k)

     pw_xy(i,j) = pw_xy(i,j)+qv(i,j,k)*coef1
     cw_xy(i,j) = cw_xy(i,j)+qcl(i,j,k)*coef1
     iw_xy(i,j) = iw_xy(i,j)+qci(i,j,k)*coef1

    end do
  end do
  u0(k)=u0(k)*coef
  v0(k)=v0(k)*coef
  t0(k)=t0(k)*coef
  tabs0(k)=tabs0(k)*coef
  q0(k)=q0(k)*coef
  qcldliq0(k)=qcldliq0(k)*coef
  qpcpliq0(k)=qpcpliq0(k)*coef
  qcldice0(k)=qcldice0(k)*coef
  qpcpice0(k)=qpcpice0(k)*coef
  p0(k)=p0(k)*coef

end do ! k

k500 = nzm
do k = 1,nzm
   kc=min(nzm,k+1)
   if((pres(kc).le.500.).and.(pres(k).gt.500.)) then
      if ((500.-pres(kc)).lt.(pres(k)-500.))then
         k500=kc
      else
         k500=k
      end if
   end if
end do


do j=1,ny
 do i=1,nx
  usfc_xy(i,j) = usfc_xy(i,j) + u(i,j,1)*dtfactor  
  vsfc_xy(i,j) = vsfc_xy(i,j) + v(i,j,1)*dtfactor  
  u200_xy(i,j) = u200_xy(i,j) + u(i,j,k200)*dtfactor  
  v200_xy(i,j) = v200_xy(i,j) + v(i,j,k200)*dtfactor  
  w500_xy(i,j) = w500_xy(i,j) + w(i,j,k500)*dtfactor
 end do
end do

if(dompi) then

  coef1 = 1./float(nsubdomains)
  do k=1,nzm
    buffer(k,1) = u0(k)
    buffer(k,2) = v0(k)
    buffer(k,3) = t0(k)
    buffer(k,4) = q0(k)
    buffer(k,5) = p0(k)
    buffer(k,6) = tabs0(k)
    buffer(k,7) = qcldliq0(k)
    buffer(k,8) = qpcpliq0(k)
    buffer(k,9) = qcldice0(k)
    buffer(k,10) = qpcpice0(k)
  end do
  call task_sum_real8(buffer,buffer1,nzm*10)
  do k=1,nzm
    u0(k)=buffer1(k,1)*coef1
    v0(k)=buffer1(k,2)*coef1
    t0(k)=buffer1(k,3)*coef1
    q0(k)=buffer1(k,4)*coef1
    p0(k)=buffer1(k,5)*coef1
    tabs0(k)=buffer1(k,6)*coef1
    qcldliq0(k)=buffer1(k,7)*coef1
    qpcpliq0(k)=buffer1(k,8)*coef1
    qcldice0(k)=buffer1(k,9)*coef1
    qpcpice0(k)=buffer1(k,10)*coef1
  end do

end if ! dompi

qn0 = qcldliq0 + qcldice0
qp0 = qpcpliq0 + qpcpice0
qv0 = q0 - qn0

!=====================================================
! UW ADDITIONS

! FIND VERTICAL INDICES OF 850MB, COMPUTE SWVP
k850 = 1
do k = 1,nzm
   if(pres(k).le.850.) then
      k850 = k
      EXIT
   end if
end do

do k=1,nzm
  coef1 = rho(k)*dz*adz(k)*dtfactor
  do j=1,ny
    do i=1,nx

     ! Saturated water vapor path with respect to water. Can be used
     ! with water vapor path (= pw) to compute column-average
     ! relative humidity.   
     swvp_xy(i,j) = swvp_xy(i,j)+qsatw(tabs(i,j,k),pres(k))*coef1
    end do
  end do
end do ! k

! ACCUMULATE AVERAGES OF TWO-DIMENSIONAL STATISTICS
do j=1,ny
 do i=1,nx
  psfc_xy(i,j) = psfc_xy(i,j) + (100.*pres(1) + p(i,j,1))*dtfactor  

  ! 850 mbar horizontal winds
  u850_xy(i,j) = u850_xy(i,j) + u(i,j,k850)*dtfactor  
  v850_xy(i,j) = v850_xy(i,j) + v(i,j,k850)*dtfactor  

 end do
end do

! COMPUTE CLOUD/ECHO HEIGHTS AS WELL AS CLOUD TOP TEMPERATURE
! WHERE CLOUD TOP IS DEFINED AS THE HIGHEST MODEL LEVEL WITH A
! CONDENSATE PATH OF 0.01 kg/m2 ABOVE.  ECHO TOP IS THE HIGHEST LEVEL
! WHERE THE PRECIPITATE MIXING RATIO > 0.001 G/KG.

! initially, zero out heights and set cloudtoptemp to SST
cloudtopheight = 0.
cloudtoptemp = sstxy(1:nx,1:ny) + t00
echotopheight = 0.
do j = 1,ny
   do i = 1,nx
      ! FIND CLOUD TOP HEIGHT
      tmp_lwp = 0.
      do k = nzm,1,-1
         tmp_lwp = tmp_lwp + (qcl(i,j,k)+qci(i,j,k))*rho(k)*dz*adz(k)
         if (tmp_lwp.gt.0.01) then
            cloudtopheight(i,j) = z(k)
            cloudtoptemp(i,j) = tabs(i,j,k)
            cld_xy(i,j) = cld_xy(i,j) + dtfactor
            EXIT
         end if
      end do
      ! FIND ECHO TOP HEIGHT
      do k = nzm,1,-1
         if (qpl(i,j,k)+qpi(i,j,k).gt.1.e-6) then
            echotopheight(i,j) = z(k)
            EXIT
         end if
      end do
   end do
end do

! END UW ADDITIONS
!=====================================================

!-----------------
! compute some sgs diagnostics:

call sgs_diagnose()

!-----------------

! recompute pressure levels, except at restart (saved levels are used).
if(dtfactor.ge.0.) call pressz()   ! recompute pressure levels

call t_stopf ('diagnose')

end subroutine diagnose
