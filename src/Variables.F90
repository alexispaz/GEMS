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


module gems_variables

! Esta structura es un vector con elementos 'character' con 'len' del tamaño las
! linea de commando. Sirve para guardar bloques o cualquier linea de commandos
! que se quiera.

! Las `labels` son un par nombre/objeto polimorfico y permite
! seleccionar objetos interact e integrate.
! En el input aparecen como `:nombre` siempre que se necesiten

! Las `variables` son pares nombre/valor, donde valor es un string
! que se expande al momento de la lectura de $nombre.
! En el input aparecen como `$nombre` cuando se queire recuperar su valor y
! como `getin nombre value` cuando se establece su valor. El valor se copia textual, no se expande a menos que aparezca como 

! TODO: Incorporate all the procedures inside polvar type

use gems_errors, only: werr
use gems_constants, only: dp,linewidth
use gems_strings, only: chset_var,chset_l
private

public    :: polvars, polvar_expand, polvar_link, polvar_hard, polvar_find, polvar_save
public    :: polvar_integrate, polvar_readonly
 
type, public   :: polvar
  character(linewidth)  :: var
  class(*),pointer      :: val=>null()
  logical               :: hard=.true.
  logical               :: readonly=.false.
  contains
end type
 
#define _NODE polvar_v
#define _TYPE type(polvar)
#include "vector_header.inc"

type(polvar_v)    :: polvars

contains

#define _NODE polvar_v
#define _TYPE type(polvar)
! #define _DUMMY class(*)
#include "vector_body.inc"
        
              
subroutine polvar_link(var,val)
! Asociate the variable with a target
character(*),intent(in)  :: var
class(*),target          :: val
type(polvar),pointer     :: pv

pv=>polvar_return(var)
      
! Need to deallocate to avoid memory leak
if(associated(pv%val)) then
  if(pv%hard) deallocate(pv%val)
endif 
       
pv%val=>val
pv%hard=.false.

end subroutine  
 
subroutine polvar_readonly(var)
character(*),intent(in)  :: var 
type(polvar),pointer     :: pv
pv=>polvar_find(var)
if(.not.associated(pv)) return
pv%readonly=.true.
end subroutine  

subroutine polvar_hard(var,val)
! Set a value hard
character(*),intent(in)  :: var
class(*)                 :: val
type(polvar),pointer     :: pv

pv=>polvar_return(var)

! Need to deallocate to allow below allocation 
if(associated(pv%val)) then
  if(pv%hard) deallocate(pv%val)
endif 

! This allocates and copy the value  
allocate(pv%val,source=val)
pv%hard=.true.

end subroutine  

function polvar_expand(var) result(w)
character(*),intent(in)  :: var
character(:),allocatable :: w
type(polvar),pointer     :: pv

pv=>polvar_find(var)

if(.not.associated(pv)) then
  w='NODEF'
  return
endif 

call polvar_get(pv,w)

end function

subroutine polvar_save(var,w)
! If the variable is already there, save the value
! If not, create one hard and make it string.
character(*),intent(in)  :: var
character(*),intent(in)  :: w
type(polvar),pointer     :: pv

pv=>polvar_return(var)
if(associated(pv%val)) then
  call werr('This variable is read only',pv%readonly)
  call polvar_set(pv,w)
else 
  call polvar_hard(var,w)
endif 
       
end subroutine  

subroutine polvar_get(pv,w)
use gems_strings, only: operator(.ich.)
! Set a polvar with the value given as a character
type(polvar),pointer                  :: pv
character(:),allocatable,intent(out)  :: w

call werr('Internal error',.not.associated(pv))
call werr('Internal error',.not.associated(pv%val))

select type(v=>pv%val)
type is (integer)
  w=.ich.v
type is (real(dp))
  ! Note that v can't be an array
  w=.ich.v
  w=trim(adjustl(w))  
type is (character(*))
  w=v
class default
  call werr('I dont know how to expand that')
end select

end subroutine
 
subroutine polvar_set(pv,w)
! Set a polvar with the value given as a character
type(polvar),pointer     :: pv
character(*),intent(in)  :: w

call werr('Internal error',.not.associated(pv))
call werr('Internal error',.not.associated(pv%val))

select type(v=>pv%val)
type is (integer)
  read(w,*) v
type is (real(dp))
  ! Note that v can't be an array
  read(w,*) v
type is (character(*))
  v=w
class default
  call werr('I dont know how to set that')
end select
      
end subroutine
       
function polvar_group(var) result(g)
use gems_groups, only:group
character(*),intent(in)  :: var
type(polvar),pointer     :: pv
type(group),pointer      :: g

call werr('Labels should start with colon `:` symbol',var(1:1)/=':')
pv=>polvar_find(var)

g=>null()
if(.not.associated(pv)) return

call werr('Variable not linked',pv%hard)
   
! Print
select type(v=>pv%val)
type is (group)
  g=>v
class default
  call werr('I dont know how to return that')
end select

end function
      
function polvar_find(var) result(pv)
! Find a variable
character(*),intent(in)  :: var
type(polvar),pointer     :: pv
integer                  :: i
  
! Search for an existing variable  
do i = 1,polvars%size
  pv => polvars%o(i)
  if (trim(pv%var) == trim(var)) return
enddo
pv=>null()

end function polvar_find
 
function polvar_return(var) result(pv)
! Find or create
character(*),intent(in)  :: var
type(polvar),pointer     :: pv
  
! call werr('Variable names should not start with a number',verify(var(1:1),chset_l)/=0)
call werr('Bad character in varaible name',verify(var,chset_var)/=0)

! Search for an existing variable  
pv=>polvar_find(var)
 
if(associated(pv)) return

call polvars%append()
pv=> polvars%o(polvars%size)
pv%var=adjustl(trim(var))


end function polvar_return
!   
! function polvar_delete(var) result(pv)
! ! Find or create
! character(*),intent(in)  :: var
! type(polvar),pointer     :: pv
! integer                  :: i
!   
! ! call werr('Variable names should not start with a number',verify(var(1:1),chset_l)/=0)
! call werr('Bad character in varaible name',verify(var,chset_var)/=0)
!
! ! Search for an existing variable  
! do i = 1,polvars%size
!   pv => polvars%o(i)
!   if (trim(pv%var) == trim(var)) exit
! enddo
! call polvars%del(pv,1)
!              
! end function polvar_delete
!              
end module gems_variables
 
     
