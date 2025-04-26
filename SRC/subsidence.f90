
subroutine subsidence()
	
use vars
use microphysics, only: micro_field, index_water_vapor, nmicro_fields, mklsadv, mktend_vadv, flag_advect
use params, only: dotracers
use tracers, only: tracer, trlsadv
implicit none

integer i,j,k,k1,k2,n
real rdz, dq
real t_vtend, q_vtend
real t_tend(nx,ny,nzm), q_tend(nx,ny,nzm)

! Initialize large-scale vertical advective tendencies.
do k = 1,nzm
   ulsvadv(k) = 0.
   vlsvadv(k) = 0.
   qlsvadv(k) = 0.
   tlsvadv(k) = 0.
end do
mktend_vadv(:,:) = 0. ! large-scale microphysical tendencies

do k=2,nzm-1
  if(wsub(k).ge.0) then
     rdz=wsub(k)/(dz*adzw(k))	
     k1 = k
     k2 = k-1 
  else
     rdz=wsub(k)/(dz*adzw(k+1))       
     k1 = k+1
     k2 = k
  end if
  do j=1,ny
    do i=1,nx
      dudt(i,j,k,na) = dudt(i,j,k,na) - rdz*(u(i,j,k1)-u(i,j,k2)) 
      dvdt(i,j,k,na) = dvdt(i,j,k,na) - rdz*(v(i,j,k1)-v(i,j,k2)) 
      t_tend(i,j,k) =  - rdz * (t(i,j,k1)-t(i,j,k2))
      q_tend(i,j,k) =  &
       - rdz * (micro_field(i,j,k1,index_water_vapor)-micro_field(i,j,k2,index_water_vapor))
      ulsvadv(k) = ulsvadv(k) - rdz*(u(i,j,k1)-u(i,j,k2)) 
      vlsvadv(k) = vlsvadv(k) - rdz*(v(i,j,k1)-v(i,j,k2)) 
    end do
  end do

end do
do k=2,nzm-1
  t_vtend = 0.
  q_vtend = 0.
  do j=1,ny
    do i=1,nx
      t(i,j,k) = t(i,j,k) + dtn * t_tend(i,j,k)
      micro_field(i,j,k,index_water_vapor) = max(0.,micro_field(i,j,k,index_water_vapor) &
                         + dtn * q_tend(i,j,k))
      t_vtend = t_vtend + t_tend(i,j,k)
      q_vtend = q_vtend + q_tend(i,j,k)
    end do
  end do
  t_vtend = t_vtend / float(nx*ny) 
  q_vtend = q_vtend / float(nx*ny) 
  ttend(k) = ttend(k) + t_vtend
  qtend(k) = qtend(k) + q_vtend
  tlsvadv(k) = t_vtend
  qlsvadv(k) = q_vtend
end do 
	
! put qlsvadv into mklsadv(:,index_water_vapor)
mklsadv(1:nzm,index_water_vapor) = qtend(1:nzm)*float(nx*ny)

! Apply large-scale vertical advection to all microphysics fields, 
!   not just water vapor or total water.  This resolves some issues
!    when index_water_vapor refers to something other than total 
!    water (i.e., vapor+cloud).
do n = 1,nmicro_fields
  if((flag_advect(n).eq.1).AND.(n.ne.index_water_vapor)) then
    q_tend(:,:,:) = 0.

    do k=2,nzm-1
      if(wsub(k).ge.0) then
        rdz=wsub(k)/(dz*adzw(k))	
        k1 = k
        k2 = k-1 
      else
        rdz=wsub(k)/(dz*adzw(k+1))       
        k1 = k+1
        k2 = k
      end if

      do j=1,ny
        do i=1,nx
          q_tend(i,j,k) = - rdz * (micro_field(i,j,k1,n)-micro_field(i,j,k2,n))
        end do
      end do
    end do

    !bloss: add tendency separately since we cannot modify micro_field before 
    !  computing the tendency at all levels.
    do k = 2,nzm-1
      do j = 1,ny
        do i = 1,nx
          micro_field(i,j,k,n) = MAX(0., micro_field(i,j,k,n) + dtn*q_tend(i,j,k) )
          mklsadv(k,n) = mklsadv(k,n) + q_tend(i,j,k)
          mktend_vadv(k,n) = mktend_vadv(k,n) + q_tend(i,j,k)
        end do
      end do
    end do
  end if
end do

! save vadv tendencies for various water isotopes/tracers
mktend_vadv(:,:) = mktend_vadv(:,:)/float(nx*ny)
mktend_vadv(:,index_water_vapor) = qlsvadv(k)
if(dotracers) then
  !bloss: Initialize tracer large-scale tendency to zero
  !!!!NOTE: THIS SHOULD BE MOVED IF HORIZONTAL ADVECTIVE TENDENCIES ARE
  !   ALSO ADDED IN SRC/forcing.f90. !!!!
  trlsadv(:,:) = 0.

  ! Apply large-scale vertical advection to tracer fields, 
  do n = 1,ntracers
    q_tend(:,:,:) = 0.

    do k=2,nzm-1
      if(wsub(k).ge.0) then
        rdz=wsub(k)/(dz*adzw(k))	
        k1 = k
        k2 = k-1 
      else
        rdz=wsub(k)/(dz*adzw(k+1))       
        k1 = k+1
        k2 = k
      end if

      do j=1,ny
        do i=1,nx
          q_tend(i,j,k) = - rdz * (tracer(i,j,k1,n)-tracer(i,j,k2,n))
        end do
      end do
    end do

    !bloss: add tendency separately since we cannot modify tracer(:,:,:,n) before 
    !  computing the tendency at all levels.
    do k = 2,nzm-1
      do j = 1,ny
        do i = 1,nx
          tracer(i,j,k,n) = MAX(0., tracer(i,j,k,n) + dtn*q_tend(i,j,k) )
          trlsadv(k,n) = trlsadv(k,n) + q_tend(i,j,k)
        end do
      end do
    end do
  end do
end if

! normalize large-scale vertical momentum forcing
ulsvadv(:) = ulsvadv(:) / float(nx*ny) 
vlsvadv(:) = vlsvadv(:) / float(nx*ny) 

end
