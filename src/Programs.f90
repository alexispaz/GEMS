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


module gems_programs
 use gems_program_types
 use gems_constants
 use gems_inq_properties 
 use gems_set_properties
 use gems_neighbour

 implicit none

 ! Caracteristicas del Ensamble
 logical     :: c_lmom=.true.
 logical     :: c_amom=.true.

 real(dp)    :: wall

 contains

! DINAMICA
 
subroutine dinamic_from_xyz(steps,fname,b_out)
use gems_errors, only: wwan
use gems_interaction, only: interact
use gems_output, only: write_out
integer,intent(in)      :: steps
integer                 :: ns,u,j,k,ioflag
logical,intent(in)      :: b_out
character(*),intent(in) :: fname
character(2)            :: sym
character(200)          :: boxflag
logical                 :: ghosted

! CHECK: I think this require full ghost updates: motion jumps are arbitrary
! Sudden atom movements require fullghost
if(ghost) then
  ghosted=fullghost
  fullghost=.true.
endif

u = find_io(30)
open(u,action='read',file=fname)
boxflag=" "

do ns = 1,steps
  dm_steps=dm_steps+1._dp

  read(u,iostat=ioflag,fmt=*) k
  if(ioflag/=0) then
    call wwan("EOF?")
    return
  endif

  read(u,'(a200)') boxflag
  if(boxflag(1:3)=='box') then
    read(boxflag(4:len(boxflag)),*) tbox(1,1), tbox(2,2), tbox(3,3)
    call box_setvars()
  endif

  do j = 1, k
    read(u,*) sym, a(j)%o%pos(1:dm)
  enddo
  call pos_changed()

  call interact(b_out)
   
  ! Escribo aca para que coincida interacción y configuracion
  if (b_out) call write_out(1,dm_steps)

enddo

if(ghost) then
  fullghost=ghosted
endif

close(u)

end subroutine
  
subroutine dinamic_from_bin(steps,b_out)
use gems_output, only: write_out
! esta subrutina funciona con subsystemas
integer,intent(in)    :: steps
integer               :: ns
logical,intent(in)    :: b_out

do ns = 1,steps
  dm_steps=dm_steps+1._dp

  ! Escribo aca para que coincida interacción y configuracion
  if (b_out) call write_out(1,dm_steps)

  !call werr('Final del archivo o error en intento de lectura',read_bin())
enddo

end subroutine
    
subroutine dinamic(steps,b_out,b_time)
use gems_errors, only: timer_start, timer_dump
use gems_set_properties, only:pos_changed
use gems_integration, only: integration_stepa,integration_stepb
use gems_interaction, only: interact
use gems_output, only: write_out
use gems_checkpoint, only:b_ckp, chpeach, write_chp
use gems_input_parsing, only: load_blk, bloques, execute_block
 
integer,intent(in)    :: steps
integer               :: ns,i
logical,intent(in)    :: b_out,b_time

! Timing
call timer_start(time)
call interact(.false.)

do ns = 1,steps
  dm_steps=dm_steps+1._dp

  ! Para dinamicas multipaso
  ! de integración
  ! dt = dtaux

  ! Integración
  call integration_stepa
  call interact(b_out)
  call integration_stepb
 
  ! Escribo aca para que coincida interacción y configuracion
  if (b_out) call write_out(1,dm_steps)

 ! Checkpoint
  if (b_ckp) then
    if (mod(dm_steps,real(chpeach,dp))==0._dp)  call write_chp(ns,steps)
  endif
   
  ! Command interpreter
  if (b_load) then
    if (mod(dm_steps,real(load_each,dp))==0._dp)  call execute_block(bloques(load_blk),1,1,1,.false.)

    ! Lisent to term signal
    if (term_signal) then
      term_signal=.false.
      exit
    endif
  endif

  ! Timming information
  if (b_time) call timer_dump(ns,time,nupd_vlist)

enddo

end subroutine
 
! TOOLS

! function posdiffer(p1,p2)
!   ! Constata si pp esta en la region de pmin
!   real(dp),intent(in)    :: p1(:),p2(:)
!   logical                :: posdiffer
!   integer                :: i
!
!   posdiffer=.false.
!   do i = 1,size(p1)
!     if (abs(p1(i)-p2(i))>=dpos) then
!       posdiffer = .true.
!       return
!     endif
!   enddo
!   
! end function
!
! subroutine escape_dinamic(g,mintime,b_out)
! ! Una dinamica comun pero controla cada tanto si escape del minimo
! type(group),intent(inout)   :: g
! real(dp),intent(in)         :: mintime
! logical,intent(in)          :: b_out
! real(dp)                    :: timeaux
! real(dp)                    :: vmin
! real(dp)                    :: pmin(g%nat*dm),pmin2(g%nat*dm)
!
! kt=(kB_ui*g%fixtemp)
! beta = 1.0_dp/kt
!
! !Guardo el minimo
! call conjgrad_minimizator(g,.false.,pmin,vmin)
!
! if (b_out) call write_out_force(1)
!
! !Corro la dinamica
! timeaux=0
! do while (time<mintime)
!
!   dm_steps = dm_steps + 1._dp
!
!   ! Escribo aca para que coincida interacción y configuracion
!   if (b_out) call write_out(1,dm_steps)
!
!   call integration_stepa
!   call interact(b_out)
!
!   time = time + dt
!   if(mod(dm_steps,200._dp)==0._dp) then
!      call conjgrad_minimizator(g,.false.,pmin,vmin)
!
!
!      if (posdiffer(pmin,pmin2)) then
!        
!        pmin = pmin2 ! Guardo el nuevo minimo
!        !pos(1,:)=pos(2,:) ! Lo meto de nuevo
!        
!        print *, frame,'Cambio!',time,time-timeaux
!
!        timeaux=time
!        if (b_out) call write_out_force(1)
!      endif
!   endif
!
!   call integration_stepb
!
! enddo
!
! end subroutine escape_dinamic

   
end module gems_programs

