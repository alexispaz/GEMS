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

! For more information read cdlist_header.inc.

subroutine _Constructor(node)
class(_NODE),target        :: node
node%next => node
node%prev => node
node%o  => null()
node%hard=.false.
end subroutine _Constructor

subroutine _AddSoftAfter(node, o)
!Insert a new node after `node`
class(_NODE),target        :: node
_CLASS, target, intent(in)  :: o
type(_NODE), pointer       :: new

allocate(new)
new%o   =>  o
new%next  => node%next
new%prev  => node
node%next%prev => new
node%next => new
new%hard = .false.
end subroutine _AddSoftAfter

subroutine _AddHardcpyAfter(node, o)
!Insert a new node after `node`
class(_NODE),target        :: node
_CLASS, target, intent(in)  :: o
type(_NODE), pointer       :: new

allocate(new)
allocate(new%o,source=o)
new%next  => node%next
new%prev  => node
node%next%prev => new
node%next => new
new%hard = .true.
end subroutine _AddHardcpyAfter

subroutine _AddHardAfter(node)
!Insert a new node after `node`
class(_NODE),target        :: node
type(_NODE), pointer       :: new

allocate(new)
allocate(new%o)
new%next  => node%next
new%prev  => node
node%next%prev => new
node%next => new
new%hard = .true.
end subroutine _AddHardAfter

subroutine _AddSoftBefore(node, o)
!Insert a new node bedore `node`
class(_NODE),target       :: node
_CLASS, target, intent(in) :: o
type(_NODE), pointer      :: new

allocate(new)
new%o  =>  o
new%next => node
new%prev => node%prev
node%prev%next => new
node%prev => new
new%hard = .false.
end subroutine _AddSoftBefore

subroutine _AddHardcpyBefore(node, o)
!Insert a new node bedore `node`
class(_NODE),target       :: node
_CLASS, target, intent(in) :: o
type(_NODE), pointer      :: new

allocate(new)
allocate(new%o,source=o)
new%next => node
new%prev => node%prev
node%prev%next => new
node%prev => new
new%hard = .true.
end subroutine _AddHardcpyBefore

subroutine _AddHardBefore(node)
!Insert a new node bedore `node`
class(_NODE),target       :: node
type(_NODE), pointer      :: new
allocate(new)
allocate(new%o)
new%next => node
new%prev => node%prev
node%prev%next => new
node%prev => new
new%hard = .true.
end subroutine _AddHardBefore
                     
subroutine _Deattach( node )
class(_NODE)   :: node
node%prev%next => node%next
node%next%prev => node%prev
end subroutine _Deattach

subroutine _DestroyNode( node )
! Destroy the node. You might want to use _Deattach first
! to avoid break the list integrity.
! WARNING: Do not forget to deallocate the node
class(_NODE)                     :: node

if(node%hard) then
  deallocate(node%o)
else
  node%o => null()
endif

! Set variables back to default... nah. 
! node % next => null()
! node % prev => null()
! node % hard = .false.
 
end subroutine _DestroyNode

subroutine _DestroyAll(node)
! Destroy the entire list.
! WARNING: Do not forget to deallocate the node
class(_NODE),target        :: node
class(_NODE),pointer       :: aux


do while(.not.associated(node%next,target=node))
  aux => node%next
  call aux%deattach()
  call aux%destroy_node()
  deallocate(aux)
enddo

call node%destroy_node()
             
end subroutine _DestroyAll
 

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

