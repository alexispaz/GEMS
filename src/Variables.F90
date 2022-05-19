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

public    :: polvar_expand, polvar_save
public    :: polvar_integrate
 
type, public   :: polvar
  character(linewidth)  :: var
  class(*),pointer      :: val=>null()
  logical               :: hard=.true.
  logical               :: readonly=.false.
  contains
end type
 
! May be linked list is better?
#define _NODE polvar_v
#define _TYPE type(polvar)
#include "vector_header.inc"

type, extends(polvar_v)  :: polvar_ev
  contains
  procedure   :: find => polvar_ev_find
  procedure   :: get => polvar_ev_return
  procedure   :: link => polvar_ev_link
  procedure   :: hard => polvar_ev_hard
  procedure   :: expand => polvar_ev_expand
  procedure   :: save => polvar_ev_save
end type

type(polvar_ev), public    :: polvars

contains

#define _NODE polvar_v
#define _TYPE type(polvar)
! #define _DUMMY class(*)
#include "vector_body.inc"
        
              
subroutine polvar_ev_link(pv,var,val,readonly)
! Link `var` in `pv` with target `val`
class(polvar_ev)         :: pv
character(*),intent(in)  :: var
class(*),target          :: val
type(polvar),pointer     :: p
logical,optional         :: readonly

p=>pv%get(var)
      
! Need to deallocate to avoid memory leak
if(associated(p%val)) then
  if(p%hard) deallocate(p%val)
endif 
       
p%val=>val
p%hard=.false.

! This is really needed?
if(present(readonly)) p%readonly=readonly

end subroutine polvar_ev_link 

subroutine polvar_ev_hard(pv,var,val)
! Set a value hard
class(polvar_ev)         :: pv
character(*),intent(in)  :: var
class(*)                 :: val
type(polvar),pointer     :: p

p=>pv%get(var)

! Need to deallocate to allow below allocation 
if(associated(p%val)) then
  if(p%hard) deallocate(p%val)
endif 

! This allocates and copy the value  
allocate(p%val,source=val)
p%hard=.true.

end subroutine polvar_ev_hard 

function polvar_ev_expand(pv,var) result(w)
class(polvar_ev)         :: pv
character(*),intent(in)  :: var
character(:),allocatable :: w
type(polvar),pointer     :: p

p=>pv%find(var)

if(.not.associated(p)) then
  w='NODEF'
  return
endif 

call polvar_get(p,w)
! allocate(character(100)::w)
! write(w,'(DT)') p
! w=trim(w)

end function polvar_ev_expand

subroutine polvar_ev_save(pv,var,w)
! If the variable is already there, save the value
! If not, create one hard and make it string.
class(polvar_ev)         :: pv
character(*),intent(in)  :: var
character(*),intent(in)  :: w
type(polvar),pointer     :: p

p=>pv%get(var)
if(associated(p%val)) then
  call werr('This variable is read only',p%readonly)
  call polvar_set(p,w)
else 
  call pv%hard(var,w)
endif 
       
end subroutine polvar_ev_save

! subroutine polvar_write(p, unit, iotype, v_list, iostat, iomsg)
! use gems_strings, only: operator(.ich.)
! class(polvar),intent(in)   :: p
! integer,intent(in)         :: unit,v_list(:)
! integer,intent(out)        :: iostat
! character(*),intent(in)    :: iotype
! character(*),intent(inout) :: iomsg  
! character(:),allocatable   :: wfmt
!    
! call werr('Internal error',.not.associated(p%val))
!
! select type(v=>p%val)
! type is (integer)
!   wfmt = '(i0)'
!   write(unit,wfmt,iostat=iostat,iomsg=iomsg)  v
! type is (real(dp))
!   wfmt = '(e25.12)'
!   write(unit,wfmt,iostat=iostat,iomsg=iomsg)  v
! type is (character(*))
!   wfmt = '(a)'
!   write(unit,wfmt,iostat=iostat,iomsg=iomsg)  v
! class default
!   call werr('I do not know how to write this',.true.)
! end select
!
! end subroutine

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
  call werr('I dont know how to expand that',.true.)
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
  call werr('I dont know how to set that',.true.)
end select
      
end subroutine
 
function polvar_ev_find(pv,var) result(p)
! Find variable `var` in polvar `pv` structure
class(polvar_ev)         :: pv
character(*),intent(in)  :: var
type(polvar),pointer     :: p
integer                  :: i
  
! Search for an existing variable  
do i = 1,pv%size
  p => pv%o(i)
  if (trim(p%var) == trim(var)) return
enddo
p=>null()

end function polvar_ev_find
 
function polvar_ev_return(pv,var) result(p)
! Find or create polvar associated with `var` in `pv` structure
class(polvar_ev)         :: pv
character(*),intent(in)  :: var
type(polvar),pointer     :: p
  
! TODO: Move this to the general variable scope
call werr('Bad character in variable name',verify(var,chset_var)/=0)

! Search for an existing variable  
p=>pv%find(var)
 
if(associated(p)) return

call pv%append()
p=> pv%o(pv%size)
p%var=adjustl(trim(var))


end function polvar_ev_return
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
 
     
