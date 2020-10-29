#include "prep.inc"

#define _AddHardcpyAfter  _PASTE(_NODE,_add_hardcpy_after)
#define _AddHardAfter     _PASTE(_NODE,_add_hard_after)
#define _AddSoftAfter     _PASTE(_NODE,_add_soft_after)
#define _AddHardcpyBefore _PASTE(_NODE,_add_hardcpy_before)
#define _AddHardBefore    _PASTE(_NODE,_add_hard_before)
#define _AddSoftBefore    _PASTE(_NODE,_add_soft_before)
#define _DestroyNode      _PASTE(_NODE,_destroy_node)
#define _DestroyAll       _PASTE(_NODE,_destroy_all)
#define _Deattach         _PASTE(_NODE,_deattach)
#define _AppendSoft       _PASTE(_NODE,_append_soft)

! ## Double linked list template.
!
! Double linked lists allow to insert a new node before the current
! one. See `list_head.inc` for information in how to use the template.


! Declaring the node
type :: _NODE
  class(_NODE), pointer :: prev => null()
  class(_NODE), pointer :: next => null()
  _CLASS, pointer        :: o   => null()
  logical               :: hard=.false.
  contains
    procedure :: add_hardcpy        => _AddHardcpyAfter
    procedure :: add_hard           => _AddHardAfter
    procedure :: add_soft           => _AddSoftAfter
    procedure :: add_before_hardcpy => _AddHardcpyBefore
    procedure :: add_before_hard    => _AddHardBefore
    procedure :: add_before_soft    => _AddSoftBefore
    procedure :: append_soft        => _AppendSoft
    procedure :: deattach           => _Deattach
    procedure :: destroy_node       => _DestroyNode
    procedure :: destroy_all        => _DestroyAll
    !final     :: node_destructor
end type _NODE

private _AddHardcpyAfter
private _AddHardAfter
private _AddSoftAfter
private _AddHardcpyBefore
private _AddHardBefore
private _AddSoftBefore
private _AppendSoft
private _DestroyNode 
private _DestroyAll  
private _Deattach    

#undef _NODE
#undef _CLASS
#undef _AddHardcpyAfter
#undef _AddHardAfter
#undef _AddSoftAfter
#undef _AddHardcpyBefore
#undef _AddHardBefore
#undef _AddSoftBefore
#undef _DestroyNode
#undef _DestroyAll 
#undef _Deattach   
#undef _AppendSoft
