#include "prep.inc"

#define _Constructor      _PASTE(_NODE,_constructor)
#define _AddHardcpyAfter  _PASTE(_NODE,_add_hardcpy_after)
#define _AddHardAfter     _PASTE(_NODE,_add_hard_after)
#define _AddSoftAfter     _PASTE(_NODE,_add_soft_after)
#define _AddHardcpyBefore _PASTE(_NODE,_add_hardcpy_before)
#define _AddHardBefore    _PASTE(_NODE,_add_hard_before)
#define _AddSoftBefore    _PASTE(_NODE,_add_soft_before)
#define _DestroyNode      _PASTE(_NODE,_destroy_node)
#define _DestroyAll       _PASTE(_NODE,_destroy_all)
#define _Deattach         _PASTE(_NODE,_deattach)

! ## Circular double linked list template.
!
! Circular double linked lists (CDLL) might not have any particular head.
! However it is needed to keep a total count of the items to avoid infinit
! loops. Leaving the first node with a null o might be a way to avoid traking
! the number of items. In respect to DLL, CDLL simplify the procedures of adding
! and removing nodes avoiding the association testings needed for the beginning
! and final nodes of a DL. On the other hand, CDLL requires a constructor.  See
! `list_head.inc` for information in how to use the template.


! Declaring the node
type :: _NODE
  class(_NODE), pointer :: prev => null()
  class(_NODE), pointer :: next => null()
  _CLASS, pointer        :: o   => null()
  logical               :: hard=.false.
  contains
    procedure :: init               => _Constructor
    procedure :: add_hardcpy        => _AddHardcpyAfter
    procedure :: add_hard           => _AddHardAfter
    procedure :: add_soft           => _AddSoftAfter
    procedure :: add_before_hardcpy => _AddHardcpyBefore
    procedure :: add_before_hard    => _AddHardBefore
    procedure :: add_before_soft    => _AddSoftBefore
    procedure :: deattach           => _Deattach
    procedure :: destroy_node       => _DestroyNode
    procedure :: destroy_all        => _DestroyAll
    !final     :: node_destroy_all
end type _NODE

private _Constructor
private _AddHardcpyAfter
private _AddHardAfter
private _AddSoftAfter
private _AddHardcpyBefore
private _AddHardBefore
private _AddSoftBefore
private _DestroyNode
private _DestroyAll
private _Deattach

#undef _NODE
#undef _CLASS
#undef _Constructor
#undef _AddHardcpyAfter
#undef _AddHardAfter
#undef _AddSoftAfter
#undef _AddHardcpyBefore
#undef _AddHardBefore
#undef _AddSoftBefore
#undef _DestroyNode
#undef _DestroyAll
#undef _Deattach

