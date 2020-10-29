#include "prep.inc"

#define _AddHardAfter    _PASTE(_NODE,_add_hard_after)
#define _AddHardcpyAfter _PASTE(_NODE,_add_hardcpy_after)
#define _AddSoftAfter    _PASTE(_NODE,_add_soft_after)
#define _AddHardHead     _PASTE(_NODE,_add_hard_head)
#define _AddSoftHead     _PASTE(_NODE,_add_soft_head)
#define _AppendHard      _PASTE(_NODE,_apend_hard)
#define _AppendSoft      _PASTE(_NODE,_apend_soft)
#define _DestroyNode     _PASTE(_NODE,_destroy_node)
#define _DestroyAll      _PASTE(_NODE,_destroy_all)
#define _Deattach        _PASTE(_NODE,_deattach)

! ## Linked list template.
!
! Invoke this template in a module like this:
!
!    module name
!
!      ... ! module declarations
!
!      #declare _NODE integer_list
!      #declare _CLASS integer
!      #include list_header.inc
!
!      ... ! module declarations
!
!    contains
!
!      ... ! procedures
!
!      #declare _NODE integer_list
!      #declare _CLASS integer
!      #include list_body.inc
!
!      ... ! procedures
!
!    end module
!
! where integer is just an example. Any other type or object can be used. A
! common way to iterate the list might be:
!
!    type(_NODE),pointer   :: list
!    class(_NODE),pointer  :: node
!
!    node => a
!    do while(associated(node))
!      ... !work with node%o ("cycle" do not advance the node)
!      node => node%next
!    enddo
!
! "cycle" fortran keyword might be wanted to advance the node. For that case it
! is possible to leave the first node of the list empty to allow this kind of
! iteration:
!
!    node => a
!    do while(associated(node%next))
!      node => node%next
!      ... !work with node%o ("cycle" advance the list)
!    enddo
!
! If the node object require deallocation prior to elimination the option
! `_ALLOCATABLE` must be defined prior to include `list_body.inc`.
! If the object is a user type defined, this type should have declared a final
! procedure.

! Declaring the node
type :: _NODE
  class(_NODE), pointer :: next => null()
  _CLASS, pointer        :: o   => null()
  logical               :: hard=.false.
  contains
    procedure :: add_hard       => _AddHardAfter
    procedure :: add_hardcpy    => _AddHardcpyAfter
    procedure :: add_soft       => _AddSoftAfter
    procedure :: add_head_hard  => _AddHardHead
    procedure :: add_head_soft  => _AddSoftHead
    procedure :: append_hard    => _AppendHard
    procedure :: append_soft    => _AppendSoft

    ! procedure :: deattach           => _Deattach
    ! Note that deattach is not used because it would consist only of 1 line:
    !   prev => node%next
    ! Where prev is the previous pointer and node is the node to deattach.

    procedure :: deattach           => _Deattach
    procedure :: destroy_node       => _DestroyNode
    procedure :: destroy_all        => _DestroyAll

    ! Bastante inestable
    !final     :: _Destroy

end type _NODE

private  _AddHardAfter
private  _AddHardcpyAfter
private  _AddSoftAfter
private  _AddHardHead
private  _AddSoftHead
private  _DestroyNode
private  _DestroyAll

#undef _NODE
#undef _CLASS
#undef _AddHardAfter
#undef _AddHardcpyAfter
#undef _AddSoftAfter
#undef _AddHardHead
#undef _AddSoftHead
#undef _AppendHard
#undef _AppendSoft
#undef _DestroyNode
#undef _DestroyAll
#undef _Deattach

