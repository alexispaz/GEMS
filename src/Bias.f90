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

 
module gems_bias
! Bias functions for potential energy


! OK, aca hay terminos de la funcion potencial
! que dependen de los terminos ya calculados. 
!
! Es decir, la fuerza se actualiza al estilo f=f+u(f)

use gems_interaction, only: tepot
use gems_program_types, only: group, atom_dclist
use gems_constants, only: dp
use gems_neighbour, only: intergroup, intergroup0_empty
use gems_errors, only:werr
use gems_input_parsing, only:readf

implicit none
save
private

type(intergroup), public, pointer :: igb=>null()

! The biased energy (the energy used for the bias)
real(dp), public  :: biased=0._dp, bias=0._dp

! switch to avoid boost the atom forces
logical, public   :: noboost=.false.

public :: bias_new, bias_cli, write_bias

contains

function bias_new(w,g1) result(igr)
type(group),intent(in)     :: g1
type(intergroup),pointer   :: igr
character(*),intent(in)    :: w

! Bulid the igr
allocate(igr)
call igr%init(g1=g1)

! Remember it
igb => igr

! No list
igr%lista => intergroup0_empty
   
! Select bias
select case(w)
case('mcammon')
  igr%interact => bias_mc
case('compress_below') 
  igr%interact => bias_lp
case('compress') 
  igr%interact => bias_lpe
case default
  call werr('Bias function not found')
endselect
 
end function
     
subroutine bias_cli(igr,w)
use gems_errors, only: werr
use gems_input_parsing, only: readf
use gems_constants, only: ev_ui
type(intergroup)          :: igr
character(*),intent(in)   :: w
real(dp)                  :: f1,f2
           
select case(w)
case('mcammon')
  call readf(f1,ev_ui)  ! E
  call readf(f2,ev_ui)  ! alpha
  call igr%p%put(1,f1)
  call igr%p%put(2,f2)
case('compress_below') 
  call readf(f1,ev_ui) ! E
  call readf(f2)       ! alpha
  call igr%p%put(1,f1)
  call igr%p%put(2,f2)
case('compress') 
  call readf(f1)       ! alpha
  call igr%p%put(1,f1)
case default
  call werr('Bias function not found')
endselect  
                           
end subroutine

subroutine bias_mc(ig)
! Convierte la itneraccion no boosteada en boosteada
! Debe existir fce como hypervector y necesita hd_fce 
class(intergroup),intent(inout)  :: ig 
real(dp)                         :: aux
real(dp)                         :: e,a
type(atom_dclist),pointer        :: la
integer                          :: i
               
e=ig%p%o(1)
a=ig%p%o(2)
          
biased=tepot
bias=0._dp
ig%epot=0._dp

e=e-tepot
 
if (e>0) then

  bias=e*e/(a+e)

  !corrijo la fuerza
  aux=bias/e
  aux=1._dp+aux*aux-2*aux

  if(noboost) return

  ig%epot=bias

  la => ig%a
  do i = 1,ig%n(1)
    la=>la%next
    la%o%force(:)=la%o%force(:)*aux
  enddo

endif
 
end subroutine bias_mc

subroutine bias_lpe(ig)
class(intergroup),intent(inout)  :: ig 
real(dp)                         :: aux
real(dp)                         :: a
type(atom_dclist),pointer        :: la
integer                          :: i

a=ig%p%o(1)

biased=tepot
bias=tepot*(a-1._dp)
ig%epot=0._dp
 
if(noboost) return

ig%epot=bias

la => ig%a
do i = 1,ig%n(1)
  la=>la%next                   
  la%o%force(:)=la%o%force(:)*a
enddo
 
end subroutine bias_lpe

subroutine bias_lp(ig)
class(intergroup),intent(inout)  :: ig 
real(dp)                         :: aux
real(dp)                         :: a, e
type(atom_dclist),pointer        :: la
integer                          :: i
 
e=ig%p%o(1)
a=ig%p%o(2)
e=e-tepot

biased = tepot
bias = 0._dp
ig%epot = 0._dp

if (e>0) then

  aux=a-1.0_dp
  bias=-e*aux

  if(noboost) return

  ig%epot=bias

  !corrijo la fuerza
  la => ig%a
  do i = 1,ig%n(1)
    la=>la%next
    la%o%force(:)=la%o%force(:)*a
  enddo

endif
 
end subroutine bias_lp


! TODO: Un vector como hd_fce permitiria tener mas de una funcion de bias por simulacion.
! subroutine interact_bias(b_out)
! use gems_interaction, only: interact
! logical,intent(in)   :: b_out 
! integer              :: i
!
! ! Componentes de la fuerza a boostear
! call interact(b_out,pot=mcpot,skip=skipnotbias)
!
! ! Obtengo la fuerza a boostear
! if (allocated(hd_fce)) then
!
!   i=natoms*dm
!   if(i>size(hd_fce)) then
!     deallocate(hd_fce)
!     allocate(hd_fce(i))
!   endif
!
!   do i=1,natoms
!     hd_fce((i-1)*dm+1:i*dm)=a(i)%o%force(:)
!   enddo
! endif
!
! ! Componentes de la fuerza restantes
! call interact(b_out,dontclean=.true.,skip=skipbias)
!
! end subroutine interact_bias

subroutine write_bias(op)
use gems_output
class(outpropa)     :: op

! Ya no hace falta esto si uso la variable bias
! call werr('Bias interaction no existe so it can be output.', .not.associated(igb))

! Escribo: La parte de energia a bostear |  El bias  | La suma
op%f(1) = op%f(1) + biased*ui_ev
op%f(2) = op%f(2) + bias*ui_ev
op%f(3) = op%f(3) + (biased+bias)*ui_ev
! op%f(5) = op%f(5) + btime
! op%f(6) = op%f(6) + nbtime

end subroutine
           

end module gems_bias
