! Copyright (c) 2020  Sergio Alexis Paz
!
!  This file is part of GEMS. GEMS is an Extensible Molecular Simulator.
!	 .
!  GEMS is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!  .
!  GEMS is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  .
!  You should have received a copy of the GNU General Public License
!  along with GEMS.  If not, see <https://www.gnu.org/licenses/>.


module gems_fire

! FAST INERTIAL RELAXATION ENGINE
! Bitzek, E., Koskinen, P., G??hler, F., Moseler, M., & Gumbsch, P. (2006).
! Structural relaxation made simple. Physical Review Letters, 97(17), 1–4.
! http://doi.org/10.1103/PhysRevLett.97.170201

 use gems_program_types
 use gems_constants
 use gems_neighbor
 use gems_integration
 use gems_interaction
 use gems_output
 use gems_errors
 use gems_checkpoint
 use gems_clinterpreter
 
 implicit none

 contains

subroutine fire(g,steps,b_out)
type(group)           :: g
integer,intent(in)    :: steps
integer               :: ns, np
real(dp)              :: dt_ini,dt_max,alfa
logical,intent(in)    :: b_out
real(dp),parameter    :: alfa_start=0.1 
real(dp),parameter    :: fdec=0.5       
real(dp),parameter    :: finc=1.1       
real(dp),parameter    :: falf=0.99      
integer,parameter     :: nmin=5


! set all initial velocities to cero (from lammps manual)
la => g%alist
do m = 1,g%nat
  la => la%next
  o => la%o
  o%vel(:)=0._dp
enddo 

! TODO: Set all the atom masses equal:
! FIRE assumes that all degrees of freedom are comparable so all the velocities
! should be on the same scale.

! Usefull for NEB

! Save initial dt
dt_ini=dt

! Timer
call timer_start(time)

np=0

dt_max=10*dt_ini ! esto hay que ir calibrandolo

alfa=alfa_start
do ns = 1,steps
  dm_steps=dm_steps+1._dp

  ! Para dinamicas multipaso (incluyendo hyperdinamica) establezco el paso
  
  ! Integración
  call integration_stepa
  call interact(b_out)
  call integration_stepb

  ! Avance del tiempo
  time = time + dt

  ! Escribo aca para que coincida interacción y configuracion
  if (b_out) call write_out(1,dm_steps)

 ! Checkpoint
  if (b_ckp) then
    if (mod(dm_steps,real(chpeach))==0._dp)  call write_chp(ns,steps)
  endif
   
  p=0
  la => g%alist
  do m = 1,g%nat
    la => la%next
    o => la%o
    p = p + dot_product(o%vel(:)*o%force(:))
  enddo

  la => g%alist
  do m = 1,g%nat
    la => la%next
    o => la%o
    v = dot_product(o%vel(:),o%vel(:))
    facv(:) = o%force(:)/sqrt(dot_product(o%force(:),o%force(:)))
    o%vel(:)=o%vel(:)*(1._dp-alfa)+alfa*v*facv(:)
  enddo 

  if(p>0) then
    np=np+1
    if (np>nmin) then
      dt=min(dt*finc,dt_max)
      alfa=alfa*falf
    endif
  else
    np=0
    dt=dt*fec
    do m = 1,g%nat
      la => la%next
      o => la%o
      o%vel(:)=0._dp
    enddo 
    alfa=alfa_start
  endif

  ! Timming information
  if (b_out) call timer_dump(ns,time)

enddo

! Load initial dt
dt=dt_init

end subroutine
  
end module gems_fire
 
