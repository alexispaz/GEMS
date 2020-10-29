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

module gems_interaction
  use gems_output
  use gems_neighbour

  implicit none

  ! Point to the ofile element that is use to output this module
  type(outfile),pointer         :: epotparts => null()

  ! Total energy
  real(dp), public               :: tepot

contains

subroutine interact(b_out,pot,dontclean)
! calcula la fuerza y energia potencial de los atomos del systema
use gems_tb
use gems_forcefield
use gems_output, only: ptime,pnframe
use gems_cvs, only: cvs, ncvs

logical,optional,intent(in)    :: dontclean
logical,intent(in)             :: b_out
real(dp),optional,intent(out)  :: pot
integer                        :: i,j
real(dp)                       :: aux
type(boundgr_l),pointer        :: ln
type (atom_dclist), pointer    :: la
class(intergroup), pointer     :: ig
class(atom),pointer            :: ap
  
 ! hard constrain
 !do i = 1,natoms
 !  if (.not.a(i)%o%bconst) cycle
 !  rd=dot(a(i)%o%pos-a(i)%o%pconst,a(i)%o%vconst)
 !  if(rd<1.d-12) cycle
 !  if (a(i)%o%lconst) then
 !    a(i)%o%pos=a(i)%o%pconst+rd
 !  else
 !    a(i)%o%pos=a(i)%o%pos-rd
 !  endif 
 !enddo 
if(nomoreigs) then
  call test_update()
else
  ! First time here.
  ! No more new interaction groups are allowed
  nomoreigs=.true.
  ! Runing set up of ghost and interacting groups
  if(ghost) call pbcghost

  call update()
endif
    
        
tepot = 0.0_dp
if(.not.present(dontclean)) then
  virial(:,:) = 0.0_dp
  do i = 1,natoms
    ap=>a(i)%o
    ap%force = 0.0_dp
    ap%epot = 0.0_dp
  enddo
  tepot=0.0_dp
endif

call tb_preinteraction

do i = 1, igr_vop%size
  ig => igr_vop%o(i)%o

  ! Check to skip
  if (ig%disable) cycle

  call ig%interact()

  tepot=tepot+ig%epot
enddo 

! Force fields
ln => boundgrs
do while( associated(ln%next) )
  ln => ln%next
  call ln%o%interact()
  tepot=tepot+ln%o%epot
enddo 

if (b_out) call write_out(2,dm_steps)

if(present(pot)) pot=tepot

! Soft constrain
do i = 1,natoms
  ap=>a(i)%o
  if (.not.ap%bconst) cycle
  if (ap%lconst) then
    aux=dot(ap%force,ap%vconst)
    ap%force = ap%vconst*aux
    aux=dot(ap%acel,ap%vconst)
    ap%acel  = ap%vconst*aux
    aux=dot(ap%vel,ap%vconst)
    ap%vel   = ap%vconst*aux
  else
    aux=dot(ap%force,ap%vconst)
    ap%force = ap%force-aux
    aux=dot(ap%acel,ap%vconst)
    ap%acel  = ap%acel-aux
    aux=dot(ap%vel,ap%vconst)
    ap%vel   = ap%vel-aux
  endif
enddo

! Agroup the forces in the head
do i=1,ncvs
  call cvs(i)%calc()
enddo
      
end subroutine interact
     
subroutine write_eparts(op)
class(outpropa)                :: op
integer                        :: i

call werr('Error. Use epart after interaction delcarations',size(op%f)/=igr_vop%size)

do i=1,igr_vop%size
  op%f(i) = op%f(i) + igr_vop%o(i)%o%epot*ui_ev
enddo  
   
end subroutine
   
end module gems_interaction
 
