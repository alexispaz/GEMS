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


module gems_replicaexchange
use gems_constants, only: dp, ie1, kB_ui
use mpi_f08
  
implicit none

contains

subroutine parallel_tempering2(steps,msteps,int_id,b_out,b_time)
use gems_mpi,           only: mpi_tpc, mpi_pc, mpi_st, mpi_dp, mpi_ie1, mpi_world
use gems_program_types, only: group, group_switch_vectorial, group_switch_objeto, time
use gems_programs,      only: dinamic
use gems_interaction,   only: interact
use gems_integration,   only: integrate, its
use gems_set_properties,only: set_scal_vel
use gems_inq_properties,only: inq_pot_energy
use gems_random,        only: ranu
use gems_errors,        only: werr,timer_start, timer_dump, logunit, wlog
integer,intent(in)      :: int_id
integer,intent(in)      :: msteps,steps
logical,intent(in)      :: b_out,b_time
type(integrate),pointer :: g
real(dp)                :: par_epot,deltau,deltab,temp,r,beta,betapar
logical                 :: switched
integer(ie1)            :: doswap
integer                 :: ns,par,exchange

g => its%o(int_id)
call werr('Integration group should be in a constant T ensamble',.not.g%b_fixt)
beta=1._dp/(kB_ui*g%fixt)

par=1
exchange=0

! Timing
call timer_start(time)

! Cambio al modo de almacenamiento vectorial
call group_switch_vectorial(g,switched)

do ns = 1,steps
  
  call dinamic(msteps,b_out,.false.)

  ! Chose pairs of partners
  if (mod(ns,2)==0) then
    par=1
  else
    par=-1
  endif

  if (mod(mpi_pc,2)==0) then
    par=mpi_pc+par
  else
    par=mpi_pc-par
  endif

  if(par<0) cycle
  if(par>=mpi_tpc) cycle

  ! Compute energy
  call interact(b_out)
  call inq_pot_energy(g)

  ! Lower MPI rank will compute swapping acceptance proability
  if(par>mpi_pc)then
           
    ! Delta Beta
    call mpi_recv(temp,1,mpi_dp,par,0,mpi_world,mpi_st) 
    betapar=1._dp/(kB_ui*temp)
    deltab=beta-betapar
           
    ! Potential energy contribution
    call mpi_recv(par_epot,1,mpi_dp,par,1,mpi_world,mpi_st) 
    deltau=(par_epot-g%epot)*deltab
      
    ! Metropoli criterion
    doswap=1_ie1
    if(deltau>0) then
      r=ranu()
      if (r<exp(-deltau)) then
        doswap=1_ie1
      else
        doswap=0_ie1
      endif
    endif
         
    ! Comunicate the swaping decision
    call mpi_send(doswap,1,mpi_ie1,par,3,mpi_world) 
    deltau=(par_epot-g%epot)*deltab
                                     
  else

    ! Supply data to the lower MPI rank
    call mpi_send(g%fixt,1,mpi_dp,par,0,mpi_world) 
    call mpi_send(g%epot,1,mpi_dp,par,1,mpi_world) 

    ! Get swaping decision
    call mpi_recv(doswap,1,mpi_ie1,par,3,mpi_world,mpi_st)  

  endif

  ! If swap, exchange necesary information
  if(doswap==1_ie1) then
      
      ! Report exchange
      exchange=exchange+1
  
      ! Exchange configurations
      call mpi_sendrecv_replace(g%pp(:),size(g%pp),mpi_dp,par,0,par,0,mpi_world,mpi_st) 
      
      ! Exchange velocities
      call mpi_sendrecv_replace(g%pv(:),size(g%pv),mpi_dp,par,0,par,0,mpi_world,mpi_st) 
      
      ! Scale velocities to the new temperature
      temp=g%fixt
      call mpi_sendrecv_replace(temp,1,mpi_dp,par,0,par,0,mpi_world,mpi_st) 
      g%pv(:)=g%pv(:)*sqrt(g%fixt/temp)

  endif   

  if (b_time) call timer_dump(ns,time)

enddo

! Retomo el modo de almacenamiento anterior
if(switched) call group_switch_objeto(g)

! Report number of exchanges
call wlog('PaT'); write(logunit,'(a,i0)') "Number of Exchange: ",exchange


end subroutine parallel_tempering2

subroutine parallel_tempering(steps,msteps,int_id,b_out,b_time,b_run)
use gems_mpi,           only: mpi_tpc, mpi_pc, mpi_st, mpi_dp, mpi_ie1, mpi_world
use gems_program_types, only: group, group_switch_vectorial, group_switch_objeto, time
use gems_programs,      only: dinamic
use gems_interaction,   only: interact
use gems_integration,   only: integrate, its
use gems_set_properties,only: set_scal_vel
use gems_inq_properties,only: inq_pot_energy
use gems_random,        only: ranu
use gems_metadynamics,  only: metadynamics,dCM,bias_point_1D
use gems_errors,        only: werr,timer_start, timer_dump, logunit, wlog
integer,intent(in)      :: int_id
integer,intent(in)      :: msteps,steps
logical,intent(in)      :: b_out,b_time
type(integrate),pointer :: g
real(dp)                :: par_epot,deltau,deltab,temp,r
logical                 :: switched
integer(ie1)            :: doswap
integer                 :: ns,par,exchange
character(10)           :: b_run

real(dp)                :: potinicv1, dpot_cv1,CVi,deltabias,par_dCM
real(dp)                :: bi_i,bi_j,bj_i,bj_j,betai,betaj
integer                 :: binCV1_1,binCV1_2
real(dp),allocatable    :: interpolation_x(:)
real(dp),allocatable    :: pot_int_y(:)
real(dp),allocatable    :: for_int_y(:) 
real(dP),allocatable    :: dfor_int_y(:)




g => its%o(int_id)
call werr('Integration group should be in a constant T ensamble',.not.g%b_fixt)
temp=g%fixt

par=1
exchange=0

! Timing
call timer_start(time)

! Cambio al modo de almacenamiento vectorial
call group_switch_vectorial(g,switched)

do ns = 1,steps

  
  if (b_run == 'dmcv') then
  call metadynamics(msteps,b_out,.false.,.true.,.false.,.false.,.false.,.false.)
else if (b_run == 'wtdcm') then
  call metadynamics(msteps,b_out,.false.,.false.,.false.,.false.,.true.,.false.) 
endif    
  ! Chose pairs of partners
  if (mod(ns,2)==0) then
    par=1
  else
    par=-1
  endif

  if (mod(mpi_pc,2)==0) then
    par=mpi_pc+par
  else
    par=mpi_pc-par
  endif

  if(par<0) cycle
  if(par>=mpi_tpc) cycle
  !Calculos en el nodo en cual estamos
  ! Compute energy
  call interact(b_out)
  call inq_pot_energy(g)
  ! Compute beta
  betai=1._dp/(kB_ui*g%fixt)
  if (b_run == 'wtdcm') then
    ! Compute Bias
    bi_i=bias_point_1D(dCM)
    ! Compute CVs 
    CVi=dCM
  end if
  ! Lower MPI rank will compute swapping acceptance proability
  if(par>mpi_pc)then
    ! Delta Beta
    call mpi_recv(betaj,1,mpi_dp,par,0,mpi_world,mpi_st) 
      deltab=betai-betaj
    ! Potential energy contribution
    call mpi_recv(par_epot,1,mpi_dp,par,1,mpi_world,mpi_st) 
      deltau=(par_epot-g%epot)*deltab
    if (b_run == 'wtdcm') then
      ! Bias con la conf de j
      call mpi_recv(par_dCM,1,mpi_dp,par,4,mpi_world,mpi_st) 
      bi_j=bias_point_1D(par_dCM)
      call mpi_recv(bj_j,1,mpi_dp,par,5,mpi_world,mpi_st) 
      call mpi_send(dCM,1,mpi_dp,par,6,mpi_world) 
      call mpi_recv(bj_i,1,mpi_dp,par,7,mpi_world,mpi_st) 
      deltabias=deltau-(betai*(bi_i-bi_j))-(betaj*(bj_j-bj_i))
    endif
!write(*,*)'pc',mpi_pc,'vec',par,par_epot-g%epot,betaj-betai,betai*(bi_i-bi_j),betaj*(bj_i-bj_j)
!write(*,*)'pc',mpi_pc,'vec',deltau,deltabias
     ! Metropoli criterion
    doswap=1_ie1
      if (b_run == 'dmcv') then
        if(deltau>0) then
          r=ranu()
          if (r<exp(-deltau)) then
            doswap=1_ie1
          else
            doswap=0_ie1
          endif
        endif
      else if (b_run == 'wtdcm') then 
        if(deltabias>0) then
          r=ranu()
          if (r<exp(-deltabias)) then
            doswap=1_ie1
          else
            doswap=0_ie1
          endif
        endif
      endif        
      
      
    ! Comunicate the swaping decision
    call mpi_send(doswap,1,mpi_ie1,par,3,mpi_world) 
                                     
  else

    ! Supply data to the lower MPI rank
    call mpi_send(betai,1,mpi_dp,par,0,mpi_world) 
    call mpi_send(g%epot,1,mpi_dp,par,1,mpi_world) 
    if (b_run == 'wtdcm') then
      call mpi_send(CVi,1,mpi_dp,par,4,mpi_world) 
      call mpi_send(bi_i,1,mpi_dp,par,5,mpi_world) 
      call mpi_recv(par_dCM,1,mpi_dp,par,6,mpi_world,mpi_st) 
      bi_j=bias_point_1D(par_dCM)
      call mpi_send(bi_j,1,mpi_dp,par,7,mpi_world) 
    endif
    ! Get swaping decision
    call mpi_recv(doswap,1,mpi_ie1,par,3,mpi_world,mpi_st)  

  endif

  ! If swap, exchange necesary information
  if(doswap==1_ie1) then
      
      ! Report exchange
      exchange=exchange+1
      write(*,*) 'cambio', exchange, mpi_pc 
      ! Exchange configurations
      call mpi_sendrecv_replace(g%pp(:),size(g%pp),mpi_dp,par,0,par,0,mpi_world,mpi_st) 
      
      ! Exchange velocities
      call mpi_sendrecv_replace(g%pv(:),size(g%pv),mpi_dp,par,0,par,0,mpi_world,mpi_st) 
      
      ! Scale velocities to the new temperature
      temp=g%fixt
      call mpi_sendrecv_replace(temp,1,mpi_dp,par,0,par,0,mpi_world,mpi_st) 
      g%pv(:)=g%pv(:)*sqrt(g%fixt/temp)

  endif   

  if (b_time) call timer_dump(ns,time)
write(*,*) 'Se termino, micro=',mpi_pc,ns,'de',steps 
enddo

! Retomo el modo de almacenamiento anterior
if(switched) call group_switch_objeto(g)

! Report number of exchanges
call wlog('PaT'); write(logunit,'(a,i0)') "Number of Exchange: ",exchange


end subroutine parallel_tempering
 
end module gems_replicaexchange

