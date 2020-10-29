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

 
! linkedlist.f90 --
!     Include file for defining linked lists where each element holds
!     the same kind of data
!
!     See the example/test program for the way to use this
!
!     Note:
!     You should only use pointer variables of this type, no
!     ordinary variables, as sometimes the memory pointed to
!     will be deallocated. The subroutines and functions
!     are designed to minimize mistakes (for instance: using
!     = instead of =>)
!
!     $Id: linkedlist.f90,v 1.3 2007/01/26 09:56:43 arjenmarkus Exp $

! TOUCHED BY ALEXIS PAZ ON 2012/08/07. PRINCIPAL CHANGES ARE:
! - I introduce events inside the type to improve when two or more list are used at
!   the same time
! - I change the data variable with a precompiled variable
!   (:%s/type(LIST_DATA)/_type/g) to avoid innecesary neast in the case of
!   intrinsic types. So this change the sintax:
!      module REAL_LISTS
!      use DATA_MODULE, ONLY: LIST_DATA => REAL_DATA
!      include "linkedlist.f90"
!      end module REAL_LISTS
!   by the new sintax:
!      module REAL_LISTS
!      #define _type real(dp)
!      #include "linkedlist.f90"
!      end module REAL_LISTS
!   also it is posible to have a module that declare all the intrinsic list like
!   befora and one generic list with the 
!   On the other hand avoid the intermediate modul TWO_LIST and the use of the
!   list led in:
!      use REAL_LISTS, only: REAL_LINKED_LIST   => LINKED_LIST
!      use STRING_LISTS, only: STRING_LINKED_LIST   => LINKED_LIST

private

!
! Define the linked-list data type
!
type LINKED_LIST
    type(LINKED_LIST), pointer :: next
    _type            :: data
    contains
        procedure :: init => list_create
        procedure :: insert => list_insert
        procedure :: get => list_get_data       
        procedure :: put => list_put_data       
        procedure :: del => list_delete_element
        procedure :: dest => list_destroy           
end type LINKED_LIST 

public :: linked_list

!
! define a private (!) interface to prevent
! mistakes with ordinary assignment
!
!interface assignment(=)
!    module procedure list_assign
!end interface
!private :: list_assign

!
! Define the subroutines and functions
!
contains

! list_assign
!     Subroutine to prevent errors with assignment
! Arguments:
!     list_left   List on the left-hand side
!     list_right  List on the right-hand side
!
! NOTE:
!     This does not work because of a private/public
!     conflict
!
!subroutine list_assign( list_left, list_right )
!    type(LINKED_LIST), INTENT(OUT)  :: list_left
!    type(LINKED_LIST), INTENT(IN)   :: list_right
!   !type(LINKED_LIST), pointer      :: list_left
!   !type(LINKED_LIST), pointer      :: list_right
!
!    !
!    ! Note the order!
!    !
!    stop 'Error: ordinary assignment for lists'
!    list_left%next => null()
!end subroutine list_assign

! list_create --
!     Create and initialise a list
! Arguments:
!     list       Pointer to new linked list
!     data       The data for the first element
! Note:
!     This version assumes a shallow copy is enough
!     (that is, there are no pointers within the data
!     to be stored)
!     It also assumes the argument list does not already
!     refer to a list. Use list_destroy first to
!     destroy up an old list.
!
subroutine list_create( list, data )
    class(LINKED_LIST), target  :: list
    _type, intent(in) :: data

    list%next => null()
    list%o =  data
end subroutine list_create

! list_destroy --
!     Destroy an entire list
! Arguments:
!     list       Pointer to the list to be destroyed
! Note:
!     This version assumes that there are no
!     pointers within the data that need deallocation
!
subroutine list_destroy( list )
    class(LINKED_LIST), target  :: list

    type(LINKED_LIST), pointer  :: current
    type(LINKED_LIST), pointer  :: next

    current => list
    do while ( associated(current%next) )
        next => current%next
        deallocate( current )
        current => next
    enddo
end subroutine list_destroy

! list_count --
!     Count the number of items in the list
! Arguments:
!     list       Pointer to the list
!
integer function list_count( list )
    class(LINKED_LIST), target  :: list

    type(LINKED_LIST), pointer  :: current
    type(LINKED_LIST), pointer  :: next

    list_count = 1
    current => list
    do while ( associated(current%next) )
        current => current%next
        list_count = list_count + 1
    enddo

end function list_count

! list_insert
!     Insert a new element
! Arguments:
!     elem       Element in the linked list after
!                which to insert the new element
!     data       The data for the new element
!
subroutine list_insert( elem, data )
    class(LINKED_LIST), target  :: elem
    _type, intent(in) :: data
    type(LINKED_LIST), pointer :: next

    allocate(next)
    next%next => elem%next
    elem%next => next
    next%o =  data
end subroutine list_insert

! list_insert_head
!     Insert a new element before the first element
! Arguments:
!     list       Start of the list
!     data       The data for the new element
!
subroutine list_insert_head( list, data )
    class(LINKED_LIST), target   :: list
    _type, intent(in) :: data
    type(LINKED_LIST), pointer :: elem

    allocate(elem)
    call list%insert(list%o)
    list%o=data
end subroutine list_insert_head

! list_delete_element
!     Delete an element from the list
! Arguments:
!     list       Header of the list
!     elem       Element in the linked list to be
!                removed
!
subroutine list_delete_element( list, elem )
    class(LINKED_LIST),target    :: list
    type(LINKED_LIST), pointer  :: elem

    type(LINKED_LIST), pointer  :: current
    type(LINKED_LIST), pointer  :: prev

    if ( associated(elem,target=list) ) then
        list%o = list%next%o
        elem => list%next
    endif

    current => list
    prev    => list
    do while ( associated(current) )
        if ( associated(current,elem) ) then
            prev%next => current%next
            deallocate( current ) ! Is also "elem"
            exit
        endif
        prev    => current
        current => current%next
    enddo
!    allocate(next)
!
!    next%next => elem%next
!    elem%next => next
!    next%o =  data
end subroutine list_delete_element

! list_get_data
!     Get the data stored with a list element
! Arguments:
!     elem       Element in the linked list
!
function list_get_data( elem ) result(data)
    class(LINKED_LIST)         :: elem
    _type            :: data
    data = elem%o
end function list_get_data

! list_put_data
!     Store new data with a list element
! Arguments:
!     elem       Element in the linked list
!     data       The data to be stored
!
subroutine list_put_data( elem, data )
    class(LINKED_LIST)           :: elem
    _type, intent(in) :: data
    elem%o = data
end subroutine list_put_data

