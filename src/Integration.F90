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
! 
!	References
!
! Grønbech-Jensen & Farago (2014). Constant pressure and temperature
! discrete-time Langevin molecular dynamics. JCP 141(19) 194108.
! DOI:10.1063/1.4901303
!
! Grønbech-Jensen & Farago (2012). A simple and effective Verlet-type
! algorithm for simulating Langevin dynamics.
! DOI:10.1080/00268976.2012.760055
!
! Kolb & Dünweg (1999). Optimized constant pressure stochastic
! dynamics. JCP 111(10) 4453. DOI:10.1063/1.479208  
!
! Allen & Tildesley (1987). Computer simulation of liquids.
! Oxford science publications. ISBN: 9780198556459.
!
! Snook (2006). The ermak and Generalised Langevin Apprroach to the
! Dynamics of Atomic, Polymeric and Colloidal System. Elsevier. ISBN:
! 9780444521293.
!
! Ermak & Buckholz (1980). Numerical integration of the Langevin equation:
! Monte Carlo simulation, J. Comput. Phys. 35 169.
!
! Andersen (1980). Molecular dynamics simulations at constant pressure
! and/or temperature. J. of Chem. Phys. 72(4) 2384. DOI:10.1063/1.439486 
!
! Papadopoulou et al (1993). Molecular Dynamics and Monte Carlo Simulations
! in the Grand Canonical Ensemble: Local versus Global Control. J. of Chem.
! Phys. 98(6) 4897. DOI:10.1063/1.464945. 
!
! Heffelfinger & van Swol (1998). Diffusion in Lennard‐Jones Fluids Using Dual
! Control Volume Grand Canonical Molecular Dynamics Simulation (DCV‐GCMD).
! J. of Chem. Phys. 100(10) 7548. DOI:10.1063/1.466849.
 
module gems_integration
use gems_program_types         !, nghost
use gems_groups
use gems_inq_properties
use gems_set_properties
use gems_constants, only: dp, dd, dm, kB_ui


implicit none

private

type, extends(group) :: integrate

  ! Algorithm parameters
  type(real_v)                  :: p
  type(integer_v)               :: i
  
  procedure(integrate0),pointer :: stepa=>null()  ! Integration step a
  procedure(integrate0),pointer :: stepb=>null()  ! Integration step b

  ! Constant of motion
  real(dp)                      :: fixt=0._dp       ! Fix temperature
  logical                       :: b_fixt=.false.   ! Fix temperature

  contains
  procedure   :: init_ext => integrate_construct
  procedure   :: dest_ext => integrate_destructor

end type
  
abstract interface
 subroutine integrate0(it)
  import integrate
  class(integrate)  :: it
 end subroutine
end interface
  
#define _NODE integrate_v
#define _TYPE type(integrate)
#include "vector_header.inc"

type(integrate_v),target,public :: its
                    

! Langevin piston factors
real(dp)          :: ppos(dm)=0._dp,pvel(dm)=0._dp, pistonf(dm)=0._dp

! Only 1 integrate per dimension allowed. This is to record the selection.
logical            :: piston(dm)=.false.

!TODO: If more than 1 ermak is used, the associated groups should not intersect
!because the o%pos_v and o%vel_v vectors will be overrided. A warning or error
!message should be added. So this is in some way the same limitation reflected
!by piston(dm), although that also include the change of box issue.

integer           :: fleu=0

! Condiciones de contorno
logical,public            :: cubic=.true.


logical    :: voter=.false.
real(dp)   :: beta 

public :: integrate
public :: integrate_cli
public :: polvar_integrate

public :: write_chppiston
public :: read_chppiston

public :: integration_reversea
public :: integration_stepa
public :: integration_stepb

contains

#define _NODE integrate_v
#define _TYPE type(integrate)
#include "vector_body.inc"

subroutine integrate_construct(it,g1)
class(integrate),target  :: it
type(group),intent(in)   :: g1

! Inicializo el objeto
call it%init()

! Inicializo los vectores de parametros
call it%p%init()
call it%i%init()

! Inicializo el grupo  
call it%attach(g1)

end subroutine integrate_construct
  
subroutine integrate_cli(it,w)
use gems_errors, only: werr, wref
use gems_input_parsing, only: readf, readl, readb
use gems_constants, only: atm_ui, linewidth
use gems_bias, only: biason
class(integrate),target  :: it
character(*), intent(in) :: w
character(:),allocatable :: w1
real(dp)                 :: f1,f2,f3,f4,f5
integer                  :: i, i1
logical                  :: b1
  
select case(w)
case ('voter')
  ! This is not an integration algotihm...  this just set the time to be the voter time or the regular time..  But it has to follow
  ! the temperature of an algorithm with constant T.  Thus, this should be aplied to an iteration type that allows to get a
  ! temperature
  call werr('A constant temperature algorithm is needed',.not.associated(it%stepa,target=ermak_a))
  call wwan('Hyperdinamics without a bias?',.not.biason)
  call readb(b1)
  voter=b1
  beta=1._dp/it%p%o(8)**2 
case ('v_verlet')
  it%stepa => velocity_verlet_a
  it%stepb => velocity_verlet_b
case ('pred_corr_5')
  it%stepa => nordsieck_predictor
  it%stepb => nordsieck_corrector
case('brownian')
  call readf(f1)
  call readf(f2)
  call set_cbrownian(it,f1,f2)
case('piston_lkd')
  call readf(f1)
  call readf(f2)
  call readf(f3,atm_ui)
  call readf(f4)
  call readf(f5)
  call set_lkd(it,f1,f3,f2,f4,f5)
case('piston_lgf')
  call readl(w1) 
  call readf(f1)
  call readf(f2)
  call readf(f3,atm_ui)
  call readf(f4)
  call readf(f5) 
  selectcase(w1)
  case('x')
    call set_lgf_flex(it,1,f1,f3,f2,f4,f5)
  case('y')
    call set_lgf_flex(it,2,f1,f3,f2,f4,f5)
  case('z')
    call set_lgf_flex(it,3,f1,f3,f2,f4,f5)
  case('xyz*','xzy*','zxy*','zyx*','yxz*','yzx*')
    call set_lgf(it,f1,f3,f2,f4,f5)
  endselect  
case('ermak')
  call readf(f2)
  call readf(f1)  
  call set_ermak(it,dt,f1,f2)
  call wref('ermak')
case('ermak_x','ermak_y','ermak_z')
  call readf(f2)
  call readf(f1)  
  selectcase(w(7:7))
  case('x')
    call set_ermak(it,dt,f1,f2,1)
  case('y')
    call set_ermak(it,dt,f1,f2,2)
  case('z')
    call set_ermak(it,dt,f1,f2,3)
  endselect  
case('gcmc')
  call readf(f1)  
  call readf(f2)
  call readf(f3)
  call readf(f4)
  call readi(i1)
  call readf(f5)
  call set_gcmc(it,f1,f2,f3,f4,i1,f5)
case('secfile')
  it%stepa => from_openfile_a
  it%stepb => from_openfile_b
  call reada(w1)
  fleu = find_io(30)
  open(fleu,action='read',file=w1)

  ! FIXME
  !call from_openfile_a Asumo que ya lei la primera parte del archivo al
  !setear el potencial y esas cosas
  ! call from_openfile_a(it)
  ! call interact(.false.)
  ! call from_openfile_b(it)

! Termostatos
case('scalvel') 
  call readf(f1)
  call readi(i1)
  call set_scalvel(it,i1,f1)
case('scalvel_after') 
  call readf(f1)
  call readi(i1)
  call readf(f2)
  call set_scalvel(it,i1,f1,f2)
case('andersen') 
  call readf(f1)
  f2=20; if(item<nitems) call readf(f2)
  call set_andersen(it,f1,f2)

case('clean') ! Delete all integration objects
         
  do i=1,its%size
    call its%o(i)%dest_ext
  enddo
  call its%destroy()

case default
  call werr('Integration algorithm unknown',.true.)
end select 

end subroutine integrate_cli

subroutine integrate_destructor(it)
class(integrate),target  :: it
       
it%stepa => null()
it%stepb => null()
call it%p%destroy()
call it%i%destroy()
call it%dest()

end subroutine integrate_destructor
                             
subroutine integration_stepa
use gems_neighbor, only: useghost,pbcghost_move
integer   :: i

do i=1,its%size
  if(associated(its%o(i)%stepa)) call its%o(i)%stepa()
enddo

! Debo computar de nuevo el volumen, por si hay cambios debido a algun/os
! piston/es
if(boxed) call box_setvars()
     
if(useghost)then
  call pbcghost_move
else
  ! Needed to avoid atoms outside box when doing neighboor list (on interact)
  call do_pbc(sys)
endif

call gindex_all_changed()
  
end subroutine
     
subroutine integration_stepb
use gems_neighbor, only: useghost,pbcghost_move
use gems_bias, only: biason, bias
integer   :: i
            
do i=1,its%size
  if(associated(its%o(i)%stepb)) call its%o(i)%stepb()
enddo
    
if(useghost)then
  call pbcghost_move()
else
  ! Needed to avoid atoms outside box 
  call do_pbc(sys)
endif

call gindex_all_changed()
  
! No creo que sea necesario
! do i=1,size(gr)
!   if(.not.gr(i)%agrouped) cycle
!   call followthehead_b(gr(i))
! enddo

! Avance del tiempo
if(voter) then
  if (biason) then
    time = time + dt*exp(bias*beta)
  else
    time = time + dt
  endif
else
  time = time + dt
endif

end subroutine
 
subroutine integration_reversea
use gems_constants, only: same_proc
use gems_neighbor, only: pbcghost_move
use, intrinsic :: iso_c_binding
type(integrate),pointer :: it
integer                 :: i

do i=1,its%size
  it=>its%o(i)

  if(same_proc(c_funloc(it%stepa),c_funloc(ermak_a))) call ermak_a_reverse(it) 
  ! if(associated(it%stepa,target=ermak_a)) call ermak_a_reverse(it) 
enddo
               
end subroutine

! Kolb Dünweg (@Kolb1999)
!------------------------ 

subroutine set_lkd(it,temp,press,pgama0,pgamav,pmass)
use gems_constants, only:kB_ui
class(integrate)        :: it 
real(dp),intent(in)     :: temp,press,pgama0,pgamav,pmass
real(dp)                :: fac0,facv

call werr('Two pistons in one coordinate not allowed',any(piston(:)))
piston(:) = .true.

it%stepa => lkd_stepa
it%stepb => lkd_stepb

fac0 = sqrt(kB_ui*temp*dt*pgama0)
facv = sqrt(kB_ui*temp*dt*pgamav)

call it%p%put(1,fac0)
call it%p%put(2,facv)
call it%p%put(3,press)
call it%p%put(4,pgama0)
call it%p%put(5,pgamav)
call it%p%put(6,pmass)

! Save the integration ensamble
it%b_fixt=.true.
it%fixt=temp

end subroutine

subroutine lkd_stepa(it)
use gems_random, only:rang
class(integrate)              :: it 
real(dp)                      :: r1
real(dp)                      :: pfix,pgama0,pgamav,pmass
real(dp)                      :: fac0,facv,boxold(dd),press
integer                       :: i,k
type ( atom_dclist ), pointer :: la

fac0   =it%p%o(1)
facv   =it%p%o(2)
pfix   =it%p%o(3)
pgama0 =it%p%o(4)
pgamav =it%p%o(5)
pmass  =it%p%o(6)

boxold(:)=box(:)

! Calculate first intermediate velocities
la => it%alist
do i = 1,it%nat 
  la => la%next
  do k=1,dd
    call rang(r1)
    la%o%acel(k) = la%o%force(k) * la%o%one_mass
    la%o%vel(k) = la%o%vel(k) &
             + 0.5_dp*dt*la%o%acel(k)                    &
             -0.5_dp*dt*pgama0*la%o%one_mass*la%o%vel(k) &
             + la%o%one_mass*fac0*r1
  enddo
enddo

! Calculate pressure 
call inq_pressure(it)
press = (it%pressure(1,1)+it%pressure(2,2)+it%pressure(3,3))/3

! Calculate half piston momentum
call rang(r1)
pvel(1) = pvel(1) + 0.5_dp*dt*(press-pfix) &
     - pgamav*0.5_dp*dt*pvel(1)/pmass &
     + facv*r1

! To keep this quasi-isotropic, I will keep the same proprotions
! of the sides of the box. So we have to update in this way:
!
! box(1)=(V'/V)**(1/3)*box(1)
!
! From the paper we have 
! V'=V+D with D=0.5*dt*pvel(1)/pmass(1)

! Calculate half box advance
box(1:dd) = box(1:dd)*(1._dp+0.5_dp*dt*pvel(1)/pmass/box_vol)**(1./dd)
box_vol = box(1)*box(2)*box(3)       

! TODO: Error for negative volume

! Calculate half positions
la => it%alist
do i = 1,it%nat 
  la => la%next
  la%o%pos(:) = la%o%pos(:) &
              + la%o%vel(:)*dt*boxold(:)*boxold(:)/(box(:)*box(:))
enddo

! Calculate final box size
do k=1,dd
  tbox(k,k) = box(k)*(1._dp+0.5_dp*dt*pvel(1)/pmass/box_vol)**(1./dd)
enddo
call box_setvars()
           
!Scale to final position and second intermediate velocity                                              
la => it%alist
do i = 1,it%nat 
  la => la%next
  la%o%pos(:) = la%o%pos(:)*box(:)/boxold(:)
  la%o%vel(:) = la%o%vel(:)*boxold(:)/box(:)
enddo

end subroutine
         
subroutine lkd_stepb(it)
use gems_random, only:rang
class(integrate)              :: it 
real(dp)                      :: r1
real(dp)                      :: pfix,pgama0,pgamav,pmass
real(dp)                      :: fac0,facv,press
integer                       :: i,k
type ( atom_dclist ), pointer :: la

fac0   =it%p%o(1)
facv   =it%p%o(2)
pfix   =it%p%o(3)
pgama0 =it%p%o(4)
pgamav =it%p%o(5)
pmass  =it%p%o(6)

! Calculate pressure again
call inq_pressure(it)
press = (it%pressure(1,1)+it%pressure(2,2)+it%pressure(3,3))/3

! Calculate final piston momentum
call rang(r1)
pvel(1) = pvel(1) + 0.5_dp*dt*(press-pfix) &
     - pgamav*0.5_dp*dt*pvel(1)/pmass &
     + facv*r1

! Calculate final velocities
la => it%alist
do i = 1,it%nat
  la => la%next
  do k=1,dd
    call rang(r1)
    la%o%acel(k) = la%o%force(k) * la%o%one_mass
    la%o%vel(k) = la%o%vel(k) &
             + 0.5_dp*dt*la%o%acel(k)                     &
             - 0.5_dp*dt*pgama0*la%o%one_mass*la%o%vel(k) &
             + la%o%one_mass*fac0*r1
  enddo          
enddo

end subroutine
 
! Grønbech-Jensen Farago (@Grønbech-Jensen2014) 
!----------------------------------------------

subroutine set_lgf(it,temp,press,pgama0,pgamav,pmass)
use gems_constants, only:kB_ui
class(integrate)        :: it 
real(dp),intent(in)     :: temp,press,pgama0,pgamav,pmass
real(dp)                :: fac0,facg0,kinetic
real(dp)                :: av,bv,facv,inv_pmass

call werr('Two pistons in one coordinate not allowed',any(piston(:)))
piston(:) = .true.
  
it%stepa => lgf_stepa
it%stepb => lgf_stepb

fac0 = sqrt(2._dp*pgama0*kB_ui*temp*dt)
facg0 = 0.5_dp*dt*pgama0
kinetic=it%nat*temp*kB_ui

inv_pmass=1.0_dp/pmass
bv=1._dp/(1._dp+0.5_dp*dt*inv_pmass*pgamav)
facv = sqrt(2._dp*pgamav*kB_ui*temp*dt)
av=(1._dp-0.5_dp*dt*inv_pmass*pgamav)*bv
   
call it%p%put(1,fac0)
call it%p%put(2,facg0)
call it%p%put(3,kinetic)
call it%p%put(4,press)

call it%p%put(5,inv_pmass)
call it%p%put(6,bv)
call it%p%put(7,facv)
call it%p%put(8,av)

! Thinks to comunicate between steps  
call it%p%put(9)

! Save the integration ensamble
it%b_fixt=.true.
it%fixt=temp

end subroutine
 
subroutine lgf_stepa(it)
use gems_random, only:rang
class(integrate)              :: it 
real(dp)                      :: r1
real(dp)                      :: box_old(dd),boxv_old
real(dp)                      :: b0,svirial,kinetic
real(dp)                      :: bv,facv,inv_pmass
real(dp)                      :: fac0,facg0
real(dp)                      :: rangv,pfix
logical,save                  :: first=.true.
integer                       :: i,k
type ( atom_dclist ), pointer :: la

fac0      = it%p%o(1)
facg0     = it%p%o(2)
kinetic   = it%p%o(3) 
pfix      = it%p%o(4) 
inv_pmass = it%p%o(5)
bv        = it%p%o(6)
facv      = it%p%o(7)

!FIXME: using of first makes this incompatible with more than 1 integrate
!object.. but more than one NPT implies more than one box? Not possible.
if(first) then
  first=.false.
                                                              
  ! Calculate pressure
  if(.not.b_gvirial) then
    call inq_virial(it)
    virial(:,:)=it%virial(:,:)
  endif
  svirial=(virial(1,1)+virial(2,2)+virial(3,3))/dd
  pistonf(1)=(svirial+kinetic)/box_vol-pfix
endif            

! Propagate box volume
call rang(rangv)
rangv=rangv*facv  
boxv_old=box_vol
box_old(:)=box(:)
box_vol=box_vol+bv*dt*(pvel(1)+0.5_dp*inv_pmass*(dt*pistonf(1)+rangv))

! Save rangv
it%p%o(9)=rangv

! Compute the box lengths
do k=1,dd
  tbox(k,k) = box(k)*(box_vol/boxv_old)**(1./dd)
enddo
call box_setvars()
             
! Propagate the positions
la => it%alist
do i = 1,it%nat
  la => la%next
  
  b0=1._dp/(1._dp+facg0*la%o%one_mass)
      
  do k=1,dd

    ! Tiro el random number n+1
    call rang(r1)
    la%o%pos_v(k) = la%o%one_mass*fac0*r1
                
    la%o%acel(k) = la%o%force(k)*la%o%one_mass
    la%o%pos(k) = box(k)/box_old(k)*la%o%pos(k) &
             + 2._dp*box(k)/(box_old(k)+box(k))*b0*dt*( la%o%vel(k)+0.5_dp*dt*la%o%acel(k)+0.5_dp*la%o%pos_v(k) ) 
  enddo
enddo

end subroutine
         
subroutine lgf_stepb(it)
use gems_random, only:rang
class(integrate)              :: it 
real(dp)                      :: pistonf_old
real(dp)                      :: a0,b0,svirial
real(dp)                      :: bv,av,inv_pmass
real(dp)                      :: rangv,pfix
real(dp)                      :: facg0,kinetic
integer                       :: i,k
type ( atom_dclist ), pointer :: la
 
facg0     = it%p%o(2)
kinetic   = it%p%o(3) 
pfix      = it%p%o(4) 
inv_pmass = it%p%o(5)
bv        = it%p%o(6)
av        = it%p%o(8)
rangv     = it%p%o(9)
                                                
! Calculate pressure
if(.not.b_gvirial) then
  call inq_virial(it)
  virial(:,:)=it%virial(:,:)
endif
svirial=(virial(1,1)+virial(2,2)+virial(3,3))/dd
pistonf_old=pistonf(1)
pistonf(1)=(svirial+kinetic)/box_vol-pfix

! Propagate box momentum
pvel(1) = av*pvel(1)                       & 
     + 0.5_dp*inv_pmass*dt*(av*pistonf_old+pistonf(1))  &
     + bv*inv_pmass*rangv

! Propagate velocities
la => it%alist
do i = 1,it%nat
  la => la%next
  
  b0=1._dp/(1._dp+facg0*la%o%one_mass)
  a0=(1._dp-facg0*la%o%one_mass)*b0
 
  do k = 1, dd
    la%o%vel(k) = a0*la%o%vel(k) &
                + 0.5_dp*dt*(a0*la%o%acel(k)+la%o%one_mass*la%o%force(k))&
                + b0*la%o%pos_v(k) 
    la%o%acel(k) = la%o%force(k)*la%o%one_mass
  enddo          
enddo
               
end subroutine
    

subroutine set_lgf_flex(it,dd,temp,press,pgama0,pgamav,pmass)
use gems_constants, only:kB_ui
class(integrate)        :: it 
real(dp),intent(in)     :: temp,press,pgama0
real(dp),intent(in)     :: pgamav,pmass
integer,intent(in)      :: dd
real(dp)                :: fac0,facg0,kinetic
real(dp)                :: av,bv,facv,inv_pmass

call werr('Two pistons in one coordinate not allowed',piston(dd))

it%stepa => lgf_flex_stepa
it%stepb => lgf_flex_stepb

piston(dd) = .true.
call it%i%put(1,dd)
call it%i%put(2,1) ! Boolean for first time

fac0 = sqrt(2._dp*pgama0*kB_ui*temp*dt)
facg0 = 0.5_dp*dt*pgama0
kinetic=it%nat*temp*kB_ui
 
inv_pmass=1.0_dp/pmass
bv=1._dp/(1._dp+0.5_dp*dt*inv_pmass*pgamav)
facv=sqrt(2._dp*pgamav*kB_ui*temp*dt)
av=(1._dp-0.5_dp*dt*inv_pmass*pgamav)*bv
      
call it%p%put(1,fac0)
call it%p%put(2,facg0)
call it%p%put(3,kinetic)
call it%p%put(4,press)

call it%p%put(5,inv_pmass)
call it%p%put(6,bv)
call it%p%put(7,facv)
call it%p%put(8,av)

! Thinks to comunicate between steps: gamav and boxvol  
call it%p%put(9)
         
! Save the integration ensamble
it%b_fixt=.true.
it%fixt=temp
 
end subroutine
 
subroutine lgf_flex_stepa(it)
! Following @Andersen1980, I wrote a extended Lagrangian with Lx, Ly and Lz
! as extended variable and using the reduce atomic variables px=rx/Lx and
! \dot{px}=\dot{rx}/Lx (see eq 3.2 and the discussion of eq 3.3). Then I got
! the Hamiltonian trough a Legendre transformation and find the "piston"
! force (after come back to real variables) as: fpx=V/Lx( (Wx+NKT)/V-Pext )
use gems_random, only:rang
class(integrate)              :: it 
real(dp)                      :: r1
real(dp)                      :: box_old
real(dp)                      :: b0
real(dp)                      :: bv,facv,inv_pmass
real(dp)                      :: fac0,facg0,kinetic
real(dp)                      :: rangv,pfix
integer                       :: i,dd
type(atom_dclist), pointer    :: la

dd =  it%i%o(1)

fac0    = it%p%o(1)
facg0   = it%p%o(2)
kinetic = it%p%o(3) 
pfix    = it%p%o(4) 

inv_pmass = it%p%o(5)
bv        = it%p%o(6)
facv      = it%p%o(7)

i=it%i%o(2)
if(i==1) then
  it%i%o(2)=0
     
  ! Calculate pressure
  if(.not.b_gvirial) then
    call inq_virial(it)
    virial(:,:)=it%virial(:,:)
  endif                     
  pistonf(dd)=(virial(dd,dd)+kinetic)/box_vol-pfix
  
endif            

! Propagate box length
call rang(rangv)
rangv=rangv*facv
box_old=box(dd)
box(dd)=box(dd)+bv*dt*(pvel(dd)+0.5_dp*inv_pmass*(dt*pistonf(dd)*box_vol/box(dd)+rangv))
tbox(dd,dd)=box(dd)

! WARNING: Cambio las longitudes de la caja pero no todavia su volumen. Esto es
! porque puede haber pistones en otras direcciones que contribuyan al volumen,
! y usen la ecuacion de aca arriba de nuevo, donde box_vol debe ser el mismo
! para cada direccion. Cuando todos los pistones recorran el paso a, recien
! ahi actualizo el volumen.
! call box_setvars()

! Save rangv  
call it%p%put(9,rangv)
 
! Propagate the positions
la => it%alist
do i = 1,it%nat
  la => la%next
  
  b0=1._dp/(1._dp+facg0*la%o%one_mass)
      
  ! Tiro el random number n+1
  call rang(r1)
  la%o%pos_v(dd) = la%o%one_mass*fac0*r1
              
  la%o%acel(dd) = la%o%force(dd)*la%o%one_mass
  la%o%pos(dd) = box(dd)/box_old*la%o%pos(dd) &
           + 2._dp*box(dd)/(box_old+box(dd))*b0*dt*( la%o%vel(dd)+0.5_dp*dt*la%o%acel(dd)+0.5_dp*la%o%pos_v(dd) ) 
enddo

end subroutine

subroutine lgf_flex_stepb(it)
use gems_random, only:rang
class(integrate)              :: it 
real(dp)                      :: pistonf_old(dm)
real(dp)                      :: a0,b0
real(dp)                      :: bv,av,inv_pmass
real(dp)                      :: rangv
real(dp)                      :: fac0,facg0,kinetic,pfix
integer                       :: i
type ( atom_dclist ), pointer :: la

dd =  it%i%o(1)

fac0      = it%p%o(1)
facg0     = it%p%o(2)
kinetic   = it%p%o(3) 
pfix      = it%p%o(4)

inv_pmass = it%p%o(5)
bv        = it%p%o(6)
av        = it%p%o(8)
rangv     = it%p%o(9)

! Calculate pressure
if(.not.b_gvirial) then
  call inq_virial(it)
  virial(:,:)=it%virial(:,:)
endif
pistonf_old(dd)=pistonf(dd)
pistonf(dd)=(virial(dd,dd)+kinetic)/box_vol-pfix

! Propagate box momentum
pvel(dd) = av*pvel(dd)                                                       & 
     + 0.5_dp*inv_pmass*dt*(av*pistonf_old(dd)+pistonf(dd))*box_vol/box(dd)  &
     + bv*inv_pmass*rangv

! Propagate velocities
la => it%alist
do i = 1,it%nat
  la => la%next
  
  b0=1._dp/(1._dp+facg0*la%o%one_mass)
  a0=(1._dp-facg0*la%o%one_mass)*b0
 
  la%o%vel(dd) = a0*la%o%vel(dd) &
              + 0.5_dp*dt*(a0*la%o%acel(dd)+la%o%one_mass*la%o%force(dd))&
              + b0*la%o%pos_v(dd) 
  la%o%acel(dd) = la%o%force(dd)*la%o%one_mass
 
enddo
               
end subroutine
                       
                                                                                        
! Brownian
!--------- 
subroutine set_cbrownian(it,temp,gama)
use gems_constants, only:kB_ui
class(integrate)    :: it
real(dp),intent(in) :: temp,gama
real(dp)            :: fac1,fac2

it%stepa => null()
it%stepb => cbrownian

fac1 = dt/gama
fac2 = sqrt(2._dp*kB_ui*temp*fac1)

call it%p%put(1,fac1)
call it%p%put(2,fac2)

! Save the integration ensamble
it%b_fixt=.true.
it%fixt=temp
 
end subroutine
     
subroutine cbrownian(it)
use gems_random, only:rang
class(integrate)              :: it 
real(dp)   :: fac1,fac2,r1,posold
integer    :: i,j
type (atom_dclist), pointer :: la

fac1 = it%p%o(1)
fac2 = it%p%o(2)

la => it%alist
do i = 1,it%nat

  do j = 1,dd

    call rang(r1)
    posold = la%o%pos(j)
    la%o%pos(j) = posold + la%o%one_mass*fac1*la%o%force(j) + la%o%one_sqrt_mass*r1*fac2

    ! Me parece que esta velocidad esta mal definida,
    ! ya que depende de gamma, y por ende la temperatura....
    la%o%vel(j) = (la%o%pos(j)-posold)/dt
  enddo
    
  la => la %next
enddo

end subroutine
  
! Ermak Buckholz (@Ermak1980)
!----------------------------
!Following the implementation in @Allen1987, Chap 9 "Brownian Dinamics", Pag
!263. A similar algorithm is in @Snook2006, Sec.  6.2.4 "A third first-order
!BD algorithm", Pag 118. See also the nice introduction in
!@Grønbech-Jensen2012.
     
subroutine set_ermak(it,dt,gama,temp,i)
class(integrate)    :: it
real(dp),intent(in) :: dt, gama,temp
real(dp)            :: cc0,cc1,cc2
real(dp)            :: sdr,sdv
real(dp)            :: crv1,crv2
integer,intent(in),optional  :: i

if(present(i)) then
  call it%i%put(1,i)
  it%stepa => ermak_a_x
  it%stepb => ermak_b_x
else
  it%stepa => ermak_a
  it%stepb => ermak_b
endif

!Integration constants.
cc0=exp(-dt*gama)
cc1=(1._dp -cc0)/gama       ! cc1=c1*dt en el libro
cc2=(1._dp -cc1/dt)/gama
 
! Desviaciones estandar (posicion y velocidad). Le falta el factor k*T/m
! que lo agrego en el algoritmo propiamente dicho:
sdr=sqrt( dt/gama* (2._dp-(3._dp-4._dp*cc0+cc0*cc0)/(dt*gama)) )
sdv=sqrt( 1._dp-cc0*cc0 )

!crv es el coeficiende de correlación posicion velocidad
crv1=(1._dp-cc0)*(1._dp-cc0)/(gama*sdr*sdv)
crv2=sqrt(1._dp-(crv1*crv1))

  
call it%p%put(1,cc0)
call it%p%put(2,cc1)
call it%p%put(3,cc2) 
call it%p%put(4,sdr) 
call it%p%put(5,sdv) 
call it%p%put(6,crv1) 
call it%p%put(7,crv2) 

! Factores utiles  
call it%p%put(8,sqrt(kB_ui*temp)) 

! Guardo gama y el paso del tiempo, por si hay que recalcular las constantes
! (e.g. en el caso de un multime timestep)
call it%p%put(9,dt)
call it%p%put(10,gama)

! Save the integration ensamble
it%b_fixt=.true.
it%fixt=temp
 
end subroutine

subroutine ermak_a(it)
use gems_random, only:rang
class(integrate)              :: it 
real(dp)                      :: skt
real(dp)                      :: r1 ,r2,f_lan
integer                       :: i,j
real(dp)                      :: cc1,cc2
real(dp)                      :: sdr,sdv
real(dp)                      :: crv1,crv2
type ( atom_dclist ), pointer :: la

cc1=it%p%o(2)
cc2=it%p%o(3) 
sdr=it%p%o(4) 
sdv=it%p%o(5) 
crv1=it%p%o(6) 
crv2=it%p%o(7) 
skt=it%p%o(8) 

! Si el paso de tiempo cambio, actualizo las constantes
! if(it%o(9)/=dt) call set_ermak(it,it%o(9),it%o(10))

! calculate new positions and velocities
la => it%alist
do i = 1,it%nat
  la => la%next

  ! Realizo la parte estocastica con las varianzas adecuadas
  ! ranr la almaceno en a%pos_v y ranv en a%pos_v, esto es porque
  ! quiero mantener la correlacion entre las dos
  f_lan = skt*la%o%one_sqrt_mass
  do j = 1, dd
    call rang(r1,r2)
    la%o%pos_v(j) = f_lan*sdr*r1
    la%o%vel_v(j) = f_lan*sdv*(crv1*r1+crv2*r2)
  enddo 

  ! Algoritmo sacado del libro de "Computer simulation of liquids" de Allen
  ! Chap 9, "Brownian Dinamics", Pag 263. 
  la%o%pos (1:dd) = la%o%pos(1:dd) + cc1*la%o%vel(1:dd) + cc2*dt*la%o%acel(1:dd) + la%o%pos_v(1:dd)

  ! Algoritmo realizado en la materia metodos computacionales 
  !la%o %pos = la%o %pos + cc1*la%o %vel + cc2*dt*la%o %acel + ranr
  !la%o %acel = la%o %force/mass
  !la%o %vel = cc0*la%o %vel + cc1*la%o %acel + ranv

enddo

end subroutine

subroutine ermak_b(it)
class(integrate)              :: it 
integer                       :: i
real(dp)                      :: cc0,cc1,cc2
type ( atom_dclist ), pointer :: la

cc0=it%p%o(1)
cc1=it%p%o(2)
cc2=it%p%o(3) 
 
! calculate new positions and velocities
la => it%alist
do i = 1,it%nat
  la => la%next
  la%o%vel (1:dd) = cc0*la%o%vel(1:dd) + (cc1-cc2)*la%o%acel(1:dd) + cc2*la%o%force(1:dd)* la%o%one_mass + la%o%vel_v(1:dd)
  la%o%acel(1:dd) = la%o%force(1:dd)* la%o%one_mass
enddo

end subroutine
                      
subroutine ermak_a_x(it)
use gems_random, only:rang
class(integrate)              :: it 
real(dp)                      :: skt
real(dp)                      :: r1 ,r2,f_lan
integer                       :: i
real(dp)                      :: cc1,cc2
real(dp)                      :: sdr,sdv
real(dp)                      :: crv1,crv2
type ( atom_dclist ), pointer :: la
integer                       :: x

cc1=it%p%o(2)
cc2=it%p%o(3) 
sdr=it%p%o(4) 
sdv=it%p%o(5) 
crv1=it%p%o(6) 
crv2=it%p%o(7) 
skt=it%p%o(8) 

x=it%i%o(1)  

! Si el paso de tiempo cambio, actualizo las constantes
! if(it%o(9)/=dt) call set_ermak(it,it%o(9),it%o(10))

! calculate new positions and velocities
la => it%alist
do i = 1,it%nat
  la => la%next

  ! Realizo la parte estocastica con las varianzas adecuadas
  ! ranr la almaceno en a%pos_v y ranv en a%pos_v, esto es porque
  ! quiero mantener la correlacion entre las dos
  f_lan = skt*la%o%one_sqrt_mass
  call rang(r1,r2)
  la%o%pos_v(x) = f_lan*sdr*r1
  la%o%vel_v(x) = f_lan*sdv*(crv1*r1+crv2*r2)

  ! Algoritmo sacado del libro de "Computer simulation of liquids" de Allen
  ! Chap 9, "Brownian Dinamics", Pag 263. 
  la%o%pos(x) = la%o%pos(x) + cc1*la%o%vel(x) + cc2*dt*la%o%acel(x) + la%o%pos_v(x)

  ! Algoritmo realizado en la materia metodos computacionales 
  !la%o %pos = la%o %pos + cc1*la%o %vel + cc2*dt*la%o %acel + ranr
  !la%o %acel = la%o %force/mass
  !la%o %vel = cc0*la%o %vel + cc1*la%o %acel + ranv

enddo

end subroutine

subroutine ermak_b_x(it)
class(integrate)              :: it 
integer                       :: i
real(dp)                      :: cc0,cc1,cc2
type ( atom_dclist ), pointer :: la
integer                       :: x

cc0=it%p%o(1)
cc1=it%p%o(2)
cc2=it%p%o(3) 

x=it%i%o(1)  
 
! calculate new positions and velocities
la => it%alist
do i = 1,it%nat
  la => la%next
  la%o%vel (x) = cc0*la%o%vel(x) + (cc1-cc2)*la%o%acel(x) + cc2*la%o%force(x)* la%o%one_mass + la%o%vel_v(x)
  la%o%acel(x) = la%o%force(x)* la%o%one_mass
enddo

end subroutine
                      
subroutine ermak_a_reverse(it)
! Invierte el paso de ermaka dado recientemente
class(integrate)              :: it 
integer                       :: i
real(dp)                      :: cc1,cc2
type ( atom_dclist ), pointer :: la

cc1=it%p%o(2)
cc2=it%p%o(3) 

la => it%alist
do i = 1,it%nat 
  la => la %next
  la%o%pos (1:dd) = la%o%pos(1:dd) - cc1*la%o%vel(1:dd) - cc2*dt*la%o%acel(1:dd) - la%o%pos_v(1:dd)
enddo

end subroutine ermak_a_reverse

  ! subroutine fixer
  !   integer             :: i,j,k
  !   type ( atom_dclist ), pointer :: la
  !   
  !   la => gr_fix%alist%next
  !   do  i = 1, gr_fix%nat
  !     j = la%o%idv
  !     do k = 1,dd
  !       if(fix(j+k-1)) la%o%pos(k) = fix_pos(j+k-1)
  !     enddo
  !     la => la%next
  !   enddo
  !
  ! end subroutine
 
! Explicit velocity Störmer-Verlet (Velocity-Verlet)
!---------------------------------------------------
! From @Grønbech-Jensen2014: "Verlet algorithm suffers from the problem that
! the total kinetic energy of a simulated system (which is supposed to be
! proportional to the temperature) becomes progressively depressed for
! increasing time step dt compared to the potential energy. Other
! thermodynamic observables also exhibit varia- tions with dt"
  
subroutine velocity_verlet_a(it)
class(integrate)             :: it
real(dp)                     :: acel_dx,acel_dv
integer                      :: i
type ( atom_dclist ), pointer :: la

acel_dv = dt*0.5_dp
acel_dx = acel_dv*dt

! Prediccion
la => it%alist
do  i = 1, it%nat
  la => la%next
  ! Ojo al sacar la siguiente linea (depende el neb de ella)
  la%o%acel(1:dd) = la%o%force(1:dd) * la%o%one_mass
  la%o%pos (1:dd) = la%o%pos(1:dd) + dt*la%o%vel(1:dd) + acel_dx*la%o%acel(1:dd)
  la%o%vel (1:dd) = la%o%vel(1:dd) + acel_dv*la%o%acel(1:dd) 
enddo

end subroutine
 
subroutine velocity_verlet_b(it)
class(integrate)    :: it
real(dp)            :: acel_dv
integer             :: i
type ( atom_dclist ), pointer :: la

acel_dv = dt*0.5_dp

la => it%alist
do  i = 1, it%nat 
  la => la%next
  la%o%acel(1:dd) = la%o%force(1:dd)*la%o%one_mass 
  la%o%vel(1:dd)  = la%o%vel(1:dd) + acel_dv*la%o%acel(1:dd)
enddo

end subroutine


! Predictor Corrector
!--------------------

subroutine nordsieck_predictor(it)
!nordsieck predictor corrector 5to orden algorithm to solve the equation of motion
class(integrate)             :: it
real(dp),dimension(dm)  :: p1,p2,p3,p4,p5,p6,p7,p8,p9
type ( atom_dclist ), pointer :: la
integer                :: i

! Prediccion
la => it%alist
do i = 1,it%nat
  la => la%next

  !set to integration unit
  la%o % vel(1:dd) = la%o %vel(1:dd)*dt

  !set the prediction factors
  p1(1:dd) = la%o % acel (1:dd)
  p2(1:dd) = la%o % acel2(1:dd)
  p4(1:dd) = la%o % acel3(1:dd)
  p7(1:dd) = la%o % acel4(1:dd)

  p3(1:dd) = 3*p2(1:dd)
  p5(1:dd) = 2*p4(1:dd)
  p6(1:dd) = 2*p5(1:dd)
  p8(1:dd) = 5*p7(1:dd)
  p9(1:dd) = 2*p8(1:dd)
  
  !make the prediction
  la%o%pos(1:dd)   = p7(1:dd) + p4(1:dd) + p2(1:dd) + p1(1:dd) + la%o%vel(1:dd) + la%o%pos(1:dd)
  la%o%vel(1:dd)   = p8(1:dd) + p6(1:dd) + p3(1:dd) + p1(1:dd) + p1(1:dd) + la%o%vel(1:dd)
  la%o%acel(1:dd)  = p9(1:dd) + p6(1:dd) + p5(1:dd) + p3(1:dd) + p1(1:dd)
  la%o%acel2(1:dd) = p9(1:dd) + p6(1:dd) + p2(1:dd)
  la%o%acel3(1:dd) = p8(1:dd) + p4(1:dd)

enddo

end subroutine 

subroutine nordsieck_corrector(it)
class(integrate)             :: it
!nordsieck predictor corrector 5to orden algorithm to solve the equation of motion
!real(dp),parameter :: c0 = 3.0d0/16.0d0,c1 = 251.0d0/360.0d0,c3 = 11.0d0/18.0d0, &
!                      c4 = 1.0d0/6.0d0, c5 = 1.0d0/60.0d0
real(dp),parameter :: c0 = 3.0d0/20.0d0,c1 = 251.0d0/360.0d0,c3 = 11.0d0/18.0d0, &
                      c4 = 1.0d0/6.0d0, c5 = 1.0d0/60.0d0 !del Allen 1989
real(dp)               :: half_mass,dt2_2
real(dp),dimension(dm)  :: p
type ( atom_dclist ), pointer :: la
integer                :: i

dt2_2 = dt*dt/2.0d0

! Corrección
la => it%alist
do i = 1,it%nat
  la => la%next

  half_mass = dt2_2*la%o % one_mass

  !set the corrector factor
  p(1:dd)=half_mass*la%o%force(1:dd)- la%o%acel(1:dd)

  !make the correction
  la%o%pos  (1:dd) =  la%o%pos  (1:dd) + p(1:dd)*c0
  la%o%vel  (1:dd) =  la%o%vel  (1:dd) + p(1:dd)*c1
  la%o%acel (1:dd) =  la%o%acel (1:dd) + p(1:dd)
  la%o%acel2(1:dd) =  la%o%acel2(1:dd) + p(1:dd)*c3
  la%o%acel3(1:dd) =  la%o%acel3(1:dd) + p(1:dd)*c4
  la%o%acel4(1:dd) =  la%o%acel4(1:dd) + p(1:dd)*c5

  !return to normal units
  la%o%vel(1:dd) =  la%o%vel(1:dd)/dt

enddo

end subroutine 

! Simple scaling
!---------------

subroutine set_scalvel(it,steps,temp,upto)
class(integrate)              :: it
real(dp),intent(in),optional   :: upto
integer,intent(in)  :: steps
real(dp),intent(in) :: temp

! The target temperature
call it%p%append(temp)

! This scale ensure thath temp is reached after a certain number of steps.
! A good value might be 33? Not sure where I readed that.
call it%p%append(0.5_dp/steps)

! Allow fluctuations to a certain degree
if(present(upto)) then
  it%stepa => null()
  it%stepb => scalvel_after
  call it%p%append(upto)  
else
  it%stepa => null()
  it%stepb => scalvel
endif
                
! Save the integration ensamble
it%b_fixt=.true.
it%fixt=temp
 
end subroutine

subroutine scalvel(it)
class(integrate)            :: it
real(dp)                   :: soft
real(dp)                   :: factor,newtemp
class(atom_dclist),pointer :: la
integer                    :: i

newtemp=it%p%o(1)
soft=it%p%o(2)

call inq_temperature(it)

if (it%temp==0._dp) then
  call wwan('Not scale factor, actual temp is cero',newtemp/=0._dp)
  return
endif

factor=(newtemp/it%temp)**soft

la => it%alist
do i = 1,it%nat
  la => la%next
  la%o%vel = la%o%vel * factor
enddo

it%temp = it%temp*(factor**2)
it%b_temp=.true.

end subroutine scalvel
 
subroutine scalvel_after(it)
class(integrate)            :: it
real(dp)                    :: upto,temp

temp=it%p%o(1)
upto=it%p%o(3)
call inq_temperature(it)

if(dabs(it%temp-temp)>upto*temp) call scalvel(it)

end subroutine scalvel_after

! Andersen
!---------

subroutine set_andersen(it,temp,nu)
use gems_constants, only:kB_ui
class(integrate)    :: it
real(dp),intent(in) :: temp
real(dp),intent(in) :: nu
real(dp)            :: factor

it%stepa => null()
it%stepb => andersen

! The target temperature
factor = sqrt(kB_ui*temp/((it%nat*dm)*(it%nat*dm)))!-drm)
call it%p%append(factor)

! The nu parameter
call it%p%append(nu)

! Save the integration ensamble
it%b_fixt=.true.
it%fixt=temp
 
end subroutine
 
subroutine andersen(it)
use gems_random,only: ranu,rang
class(integrate)              :: it
real(dp)                      :: r
real(dp)                      :: nudt,factor,vel_med
class(atom_dclist),pointer     :: la
integer                       :: i,j

! Factor nu
factor=it%p%o(1)
nudt=it%p%o(2)*dt

la => it%alist
do i=1,it%nat
  la => la%next

  !Test para colision
  r=ranu()
  if(r>nudt)cycle

  !Particulas que colisionan
  vel_med = factor*la%o%one_sqrt_mass
  do j = 1, dm
    call rang(r)
    la%o%vel(j) = vel_med*r
  enddo 
  
enddo

end subroutine andersen

! Gran Canonic Montecarlo (@Papadopoulou1993)
! -------------------------------------------
! Grand canonical Monte Carlo in a control volume . It fulfill the
! distribution for the given thermodynamic activity.  The control volume is
! the box slice between `z1` and `z2`. The number of attempted adjustments is
! `nadj` following @Heffelfinger1998
 
subroutine set_gcmc(g, z1, z2, act, rc, nadj, temp)
use gems_constants, only:kB_ui
class(integrate)    :: g
real(dp),intent(in) :: z1, z2, act, rc, temp
integer,intent(in)  :: nadj

! TODO: Check that all the atoms belong to the same kind (have the same name,
! etc and belong to same groups) so the atom template is any of them
 
g%stepa => null()
g%stepb => gcmc

call g%p%append(z1)
call g%p%append(z2)
call g%p%append(act)
call g%p%append(rc)
call g%p%append(temp)
call g%i%append(nadj)

end subroutine
 
subroutine gcmc(g)
! For rigid spheres and considering a particular area between two xy planes.
! TODO: Generalize with a boltzman energy
! TODO: Check overlap with a reference group (for mixtures)
use gems_neighbor, only: useghost,ghost_from_atom,ngindex,ngroup
use gems_random, only: ranu,rang
use gems_constants, only: dm, kB_ui
class(integrate)           :: g
class(ngroup),pointer      :: ng
real(dp)                   :: z1,z2,act 
real(dp)                   :: r(3), vd(3), dr, v, rc, temp, beta
type(atom_dclist), pointer :: la
type(atom),pointer         :: o, ref
integer                    :: nadj
integer                    :: i,j,n,m
 
z1=g%p%o(1)
z2=g%p%o(2)
act=g%p%o(3)
rc=g%p%o(4)
temp=g%p%o(5)
nadj=g%i%o(1)  
           
! Compute volume
v=box(1)*box(2)*(z2-z1)
     
! Count particles in the control volume
n=0
la=>g%alist
do j=1,g%nat
  la=>la%next
  if(la%o%pos(3)<z1.or.la%o%pos(3)>z2) cycle
  n=n+1
enddo
    
! Point to an atom that will work as template
! in order to add new atoms into groups.
       
! Attempted adjustments 
adj: do i=1,nadj

  ref => g%alist%next%o
  beta = sqrt(kB_ui*temp*ref%one_mass)
  call werr('No more particles',.not.associated(ref))

  ! Creation attempt
 if (ranu()<0.5) then
       
    ! Random coordinates
    r(1)=ranu()*box(1)
    r(2)=ranu()*box(2)
    r(3)=ranu()*(z2-z1)+z1

    ! Check overlap
    la=>g%alist
    do j=1,g%nat
      la=>la%next
      o => la%o

      ! ! Skip particles outside the control volume.
      ! FIXME: consider PBC
      ! if(o%pos(3)<z1-rc) cycle
      ! if(o%pos(3)>z2+rc) cycle

      ! Skip if overlapping
      vd(:) = atom_distancetopoint(o,r)
      dr = dot_product(vd,vd)
      if(dr<rc*rc) cycle adj

    enddo

    ! Metropolis acceptance
    if(act*v/(n+1)>ranu()) then
   
      ! Add particle
      n=n+1

      ! Initialize particle from template.
      allocate(o)
      call o%init()
      call atom_asign(o,ref)
      o%pos(:)=r(:)
       
      ! Give a velocity from maxwell-boltzman distribution
      ! TODO: Remove CM of added particles.
      do j = 1,dm
        call rang(dr)
        la%o%vel(j) = beta*dr
      enddo

      ! Add to the same groups of the template.
      do j=1,ref%ngr
        call ref%gr(j)%o%attach(o)
      enddo

      ! Create ghost images
      if(useghost) call ghost_from_atom(o)
             
      ! Free pointer
      o=>null()
           
    endif

  ! Destruction attempt
  else 
            
    ! Metropolis acceptance
    if(n/(v*act)>ranu()) then
                  
      ! Sort particle
      m=floor(ranu()*n)+1
      if(m>n) m=n

      la=>g%alist
      do j=1,g%nat
        la=>la%next
        o => la%o

        ! Skip particles outside the control volume.
        if(o%pos(3)<z1) cycle
        if(o%pos(3)>z2) cycle

        m=m-1  
        if(m==0) exit
      enddo
      call werr('Particle sorted does not exists',m>0)
              
      ! Remove particle
      n=n-1
      call o%dest()
      deallocate(o)
         
    endif
           
  endif 
   
enddo adj
      
           
endsubroutine

 

! Follow file
!------------

subroutine from_openfile_a(it)
! read atoms from file
use gems_elements, only: ncsym
class(integrate)             :: it
integer                    :: i,j,io
type(atom_dclist),pointer  :: la
character(ncsym)           :: sym
real(dp)                   :: acel_dx,acel_dv

acel_dv = dt*0.5_dp
acel_dx = acel_dv*dt

read(fleu,iostat=io,fmt=*) j
if(io/=0)  then
  call wstd(); write(logunit,*) 'File ended or incorrect input'
  return
endif
if(j/=it%nat) then
  call wstd(); write(logunit,*) 'Not a correct number of atoms in file'
endif
read(fleu,*)

la => it%alist
do i = 1,it%nat
  la => la%next
  la%o%acel(1:dd) = la%o%force(1:dd) * la%o%one_mass
  read(fleu,*) sym, la%o%pos
  la%o%vel (1:dd) = la%o%vel(1:dd) + acel_dv*la%o%acel(1:dd) 
enddo

end subroutine from_openfile_a
 
subroutine from_openfile_b(it)
class(integrate)             :: it
real(dp)            :: acel_dv
integer             :: i
type ( atom_dclist ), pointer :: la

acel_dv = dt*0.5_dp

la => it%alist
do  i = 1,it%nat
  la => la%next
  la%o%acel(1:dd) = la%o%force(1:dd)*la%o%one_mass 
  la%o%vel(1:dd)  = la%o%vel(1:dd) + acel_dv*la%o%acel(1:dd)
enddo

end subroutine from_openfile_b

! PBC
!----
      
! subroutine triclinic_pbc
!  !Si ahora no tengo una celda cubica, sino los vectores que indican los lados
!  !de una celda triclinica, voy a tener que la particula esta fuera de la caja
!  !en x siempre que x>tbox(1,1)+y*tbox(1,2)/tbox(1,1)+z*tbox(1,3)/tbox(1,1) donde
!  !tbox(:,:) es la matriz de los vectores de los lados organizada en columnas.
!  !Esto lo deduje partiendo de un cubo y deformandolo de a poco. En forma
!  !matricial es:
!  !
!  !|x|   |         0          tbox(1,2)/tbox(1,1) tbox(1,3)/tbox(1,1)| |x|   |tbox(1,1)|  
!  !|y| > |tbox(1,2)/tbox(2,3)          0          tbox(1,3)/tbox(2,2)| |y| + |tbox(2,2)| 
!  !|z|   |tbox(1,2)/tbox(3,3) tbox(1,3)/tbox(3,3)          0         | |z|   |tbox(3,3)|  
!  !
!  ! Entonces redefino tbox para que sea la primera matriz y box(1:3) la
!  ! diagonal. Este es un caso mas general que el cubico, pues en este la matriz
!  ! es cero, y la diagonal queda box(1:3)
!
!  integer                 :: i,k
!  real(dp)                :: aux(dm)
!
!   do i=1,natoms
!     aux=matmul(tbox,a(i)%pos)
!     do k=1,dm
!       if (a(i)%pbc(k)) then
!         if(a(i)%pos(k)>=box(k)+aux(k)) then
!           a(i)%pos(k)=a(i)%pos(k)-box(k)
!           a(i)%boxcr(k)=a(i)%boxcr(k)+1
!         elseif(a(i)%pos(k)<aux(k)) then
!           a(i)%pos(k)=a(i)%pos(k)+box(k)
!           a(i)%boxcr(k)=a(i)%boxcr(k)-1
!         endif
!       endif
!     enddo
!   enddo
! end subroutine


subroutine  write_chppiston(chpunit)
integer,intent(in)      :: chpunit
integer     :: i

write(chpunit) piston(1:dm)
do i=1,dm
  if (piston(i)) then
    write(chpunit) ppos(i)
    write(chpunit) pvel(i)
    write(chpunit) pistonf(i)
  endif
enddo
end subroutine write_chppiston
                

subroutine read_chppiston(chpunit)
integer,intent(in)      :: chpunit
integer                 :: i
logical                 :: aux_piston(dm)

read(chpunit) aux_piston(1:dm)
do i=1,dm
  if (aux_piston(i)) then
    read(chpunit) ppos(i)
    read(chpunit) pvel(i)
    read(chpunit) pistonf(i)
  endif
enddo   

end subroutine read_chppiston
              

! Variables and Labels

function polvar_integrate(var) result(g)
use gems_variables, only: polvar, polvar_find
use gems_errors, only: werr
character(*),intent(in)  :: var
type(polvar),pointer     :: pv
type(integrate),pointer      :: g


call werr('Labels should start with colon `:` symbol',var(1:1)/=':')
pv=>polvar_find(var)

g=>null()
if(.not.associated(pv)) return

call werr('Variable not linked',pv%hard)
 
! Print
select type(v=>pv%val)
type is (integrate)
  g=>v
class default
  call werr('I dont know how to return that',.true.)
end select

end function
                 
end module gems_integration

