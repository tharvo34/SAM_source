subroutine pressure

! call a pressure solver

use grid, only: pressure_solver, nstep, icycle, masterproc
implicit none

integer :: status = 0

call t_startf ('pressure')

if(pressure_solver.eq.0) then
   ! Use pressure_big by default
   call pressure_big(status)
   if(status.ne.0) then
      ! If pressure_big does not work because of dimensionality (2D) or
      !   domain decomposition/grid, try pressure_orig instead
      pressure_solver = 1
   end if
end if

if(pressure_solver.eq.1) then
  call pressure_orig
end if

if(nstep+icycle.eq.2) then
  if(masterproc) then
    if(pressure_solver.eq.0) write(*,*) 'Using pressure_big.f90'
    if(pressure_solver.eq.1) write(*,*) 'Using pressure_orig.f90'
  end if
end if

call t_stopf ('pressure')

end
