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

use gems_groups, only: group, atom_dclist
use gems_constants, only: dp
use gems_neighbor, only: ngroup
use gems_errors, only:werr
use gems_input_parsing, only:readf

implicit none
save
private

! The actual potential energy, should be set externally
real(dp), public  :: cepot=0._dp

! The biased energy (the energy used for the bias)
real(dp), public  :: biased=0._dp, bias=0._dp

! switch to avoid boost the atom forces
logical, public   :: noboost=.false.
                         
! switch to avoid boost the atom forces
logical, public   :: biason=.false.
                         
public :: bias_new, bias_cli, write_bias, compress_below


type,extends(ngroup)  :: mccammon
  real(dp)    :: e,alpha
  contains
  procedure   :: interact => mc_bias
  procedure,nopass :: cli => bias_cli
end type
 
type,extends(ngroup)  :: compress_below
  real(dp)    :: e,alpha
  contains
  procedure   :: interact => lpe_bias
  procedure,nopass :: cli => bias_cli
end type
 
type,extends(ngroup)  :: compress
  real(dp)    :: alpha
  contains
  procedure   :: interact => lp_bias
  procedure,nopass :: cli => bias_cli
end type

contains

subroutine bias_new(g,w)
character(*),intent(in)   :: w
class(ngroup),pointer  :: g
              
select case(w)
case('mcammon')        ; allocate(mccammon::g)  
case('compress_below') ; allocate(compress_below::g) 
case('compress')       ; allocate(compress::g) 
case default
  call werr('Bias function not found',.true.)
endselect  
g%lista => null()
     
end subroutine bias_new
           
subroutine bias_cli(g)
use gems_errors, only: werr
use gems_input_parsing, only: readf
use gems_constants, only: ev_ui
class(ngroup),intent(inout) :: g
real(dp)                        :: f1,f2
           
biason=.true.

select type(g)
type is(mccammon)
  call readf(f1,ev_ui)  ! E
  call readf(f2,ev_ui)  ! alpha
  call mc_set(g,f1,f2)
type is(compress_below) 
  call readf(f1,ev_ui) ! E
  call readf(f2)       ! alpha
  call lpe_set(g,f1,f2)
type is(compress) 
  call readf(f1)       ! alpha
  call lp_set(g,f1)
class default
  call werr('Bias function not found',.true.)
endselect  

end subroutine

subroutine mc_set(g,f1,f2)
type(mccammon)       :: g
real(dp),intent(in)  :: f1,f2
g%e=f1
g%alpha=f2
end subroutine
             
subroutine mc_bias(g)
! Convierte la itneraccion no boosteada en boosteada
! Debe existir fce como hypervector y necesita hd_fce 
class(mccammon),intent(inout)  :: g 
real(dp)                       :: aux
real(dp)                       :: e, lbias
type(atom_dclist),pointer      :: la
integer                        :: i
               
biased=cepot
g%epot=0._dp

e=g%e-cepot
 
if (e>0) then

  lbias=e*e/(g%alpha+e)

  !corrijo la fuerza
  aux=lbias/e
  aux=1._dp+aux*aux-2*aux
 
  g%epot=lbias
  bias=bias+lbias
  
  if(noboost) return

  la => g%alist
  do i = 1,g%nat
    la=>la%next
    la%o%force(:)=la%o%force(:)*aux
  enddo

endif
 
end subroutine mc_bias
  
subroutine lpe_set(g,f1,f2)
type(compress_below)      :: g
real(dp),intent(in)       :: f1,f2
g%e=f1
g%alpha=f2
end subroutine
               
subroutine lpe_bias(g)
class(compress_below),intent(inout)  :: g 
real(dp)                             :: aux
real(dp)                             :: e, lbias
type(atom_dclist),pointer            :: la
integer                              :: i
 
e=g%e-cepot

biased = cepot
g%epot = 0._dp

if (e>0) then

  aux=g%alpha-1.0_dp
  lbias=-e*aux
 
  g%epot=lbias
  bias=bias+lbias
 
  if(noboost) return

  !corrijo la fuerza
  la => g%alist
  do i = 1,g%nat
    la=>la%next
    la%o%force(:)=la%o%force(:)*g%alpha
  enddo

endif
 
end subroutine lpe_bias
   
subroutine lp_set(g,f1)
type(compress)           :: g
real(dp),intent(in)      :: f1
g%alpha=f1
end subroutine
  
subroutine lp_bias(g)
class(compress),intent(inout)  :: g 
real(dp)                       :: lbias
type(atom_dclist),pointer      :: la
integer                        :: i

biased=cepot
lbias=cepot*(g%alpha-1._dp)
g%epot=0._dp
 
g%epot=lbias
bias=bias+lbias
       
if(noboost) return

la => g%alist
do i = 1,g%nat
  la=>la%next                   
  la%o%force(:)=la%o%force(:)*g%alpha
enddo
 
end subroutine lp_bias
                  

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
use gems_output, only: outpropa
use gems_constants, only: ui_ev
class(outpropa)     :: op

! Escribo: La parte de energia a bostear |  El bias  | La suma
op%f(1) = op%f(1) + biased*ui_ev
op%f(2) = op%f(2) + bias*ui_ev
op%f(3) = op%f(3) + (biased+bias)*ui_ev
! op%f(5) = op%f(5) + btime
! op%f(6) = op%f(6) + nbtime

end subroutine write_bias
           

end module gems_bias
