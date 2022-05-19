! Copyright (c) 2020  Sergio Alexis Paz
!
!  This file is part of GEMS. GEMS is an Extensible Molecular Simulator.
!   .
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
use gems_neighbor

implicit none

! Point to the ofile element that is use to output this module
type(outfile),pointer         :: epotparts => null()

! Total energy
real(dp), public               :: tepot

contains

subroutine interact_new(g,w1)
use gems_constants,     only: linewidth
use gems_input_parsing, only: reada, readi
use gems_bias,          only: bias_new
use gems_fields,        only: field_new
use gems_pairs,         only: pair_new
use gems_tb,            only: smatb
use gems_graphs,        only: graph_new
class(ngroup),pointer  :: g
character(*)           :: w1
character(:),allocatable   :: w2
   
call reada(w2)
selectcase(w1)
case('bias')
  call bias_new(g,w2)
case('field') 
  call field_new(g,w2)
case('graph') 
  call graph_new(g,w2) 
case('pair') 
  call pair_new(g,w2)
case('tb')
  allocate(smatb::g) 
  call reread(0)
case default
  call werr('Interaction not found',.true.)
endselect    

end subroutine interact_new

subroutine interact(b_out,pot,dontclean)
! calcula la fuerza y energia potencial de los atomos del systema
use gems_tb
use gems_forcefield
use gems_cvs, only: cvs, ncvs
use gems_bias, only: bias, cepot

logical,optional,intent(in)    :: dontclean
logical,intent(in)             :: b_out
real(dp),optional,intent(out)  :: pot
integer                        :: i
real(dp)                       :: aux
type(boundgr_l),pointer        :: ln
class(ngroup), pointer         :: g
class(atom),pointer            :: ap
type(atom_dclist),pointer :: la
  
call test_update()
  
bias = 0._dp
tepot = 0._dp
cepot = 0._dp
if(.not.present(dontclean)) then
  virial(:,:) = 0._dp
  la=>sys%alist
  do i = 1,sys%nat
    la=>la%next
    la%o%force = 0._dp
    la%o%epot = 0._dp
  enddo
  la => ghost%alist
  do i = 1,ghost%nat
    la=>la%next
    la%o%force = 0._dp
    la%o%epot = 0._dp
  enddo 
  tepot=0._dp
endif

do i = 1, ngindex%size
  g => ngindex%o(i)%o
  call g%interact()
  tepot=tepot+g%epot
  cepot=tepot
enddo 

! Force fields
ln => boundgrs
do while( associated(ln%next) )
  ln => ln%next
  call ln%o%interact()
  tepot=tepot+ln%o%epot
  cepot=tepot
enddo 

if (b_out) call write_out(2,dm_steps)

if(present(pot)) pot=tepot

! Agroup the forces in the head
do i=1,ncvs
  call cvs(i)%calc()
enddo
           
end subroutine interact
     
subroutine write_eparts(op)
class(outpropa)                :: op
integer                        :: i

call werr('Error. Use epart after interaction delcarations',size(op%f)/=ngindex%size)

do i=1,ngindex%size
  op%f(i) = op%f(i) + ngindex%o(i)%o%epot*ui_ev
enddo  
   
end subroutine

function polvar_interact(var) result(g)
use gems_variables, only: polvar, polvars
use gems_errors, only: werr
character(*),intent(in)    :: var
type(polvar),pointer       :: pv
class(ngroup),pointer  :: g


call werr('Labels should start with colon `:` symbol',var(1:1)/=':')
pv=>polvars%find(var)

g=>null()
if(.not.associated(pv)) return

call werr('Variable not linked',pv%hard)
 
! Print
select type(v=>pv%val)
class is (ngroup)
  g=>v
class default
  call werr('I dont know how to return that',.true.)
end select

end function polvar_interact
    
end module gems_interaction
 
