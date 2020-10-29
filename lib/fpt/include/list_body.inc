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

! Hard allocation of a type migth required the allocation of its internal variables. If this allocation are slow (huge internal
! variables) it might be wanted to allocate the object and not its internal variables. May be it is possible to declare an
! allocation procedure for the type that requires minima allocation size and then to increase the size a reallocation should take
! place. However this is not considered here.
! 
! Soft and hard booleans allow for massive destroy procedures like _DestroyAll However, If _CLASS correspond to a polymorfphic
! object, the hard allocations without source will fail (i.e. `class(*)`) or only allow the allocation the parent class (e.g.
! class(animals). Thus, procedures that creates hard nodes do not allocate the object by themself unless a source is included.
! Therefore, the hard allocation without source should be something like:
!
! new=_AppendHard(node)
! allocate(node%o)
!
! Thus, hard boolean is not the only thing to check at destroy procedures. Associated should be used also if a robust deallocation is
! required.

! For more information read list_header.inc.

function _AppendSoft(node,o) result(new)
!Insert a new node after `node`
!Note that o can be null()
  class(_NODE),target        :: node
  _CLASS, target, intent(in)  :: o
  type(_NODE), pointer       :: new 
  type(_NODE),pointer       :: ln
  
  ! Go to the last
  ln => node
  do while( associated(ln%next) )
    ln => ln%next
  enddo    

  allocate(new)
  new%o   =>  o
  new%hard = .false.
  new%next  => ln%next
  ln%next => new                    

end function _AppendSoft

function _AppendHard(node) result(new)
!Insert a new node at the end. It is a function
!in order to be abble to use the node and do not
!search again
  class(_NODE),target       :: node
  type(_NODE), pointer      :: new
  type(_NODE),pointer       :: ln
  
  ! Go to the last
  ln => node
  do while(associated(ln%next))
    ln => ln%next
  enddo    

  allocate(new)
  new%hard = .true.
  new%next  => ln%next
  ln%next => new

  ! Remember to allocate the object!

end function _AppendHard
 
subroutine _AddSoftAfter(node, o)
!Insert a new node after `node`
  class(_NODE)               :: node
  _CLASS, target, intent(in)  :: o
  type(_NODE), pointer       :: new

  allocate(new)
  new%o   =>  o
  new%hard = .false.
  new%next  => node%next
  node%next => new
end subroutine _AddSoftAfter

subroutine _AddHardAfter(node)
!Insert a new node after `node`
  class(_NODE)               :: node
  type(_NODE), pointer       :: new

  allocate(new)
  new%hard = .true.
  new%next  => node%next
  node%next => new

  ! Remember to allocate the object!
  
end subroutine _AddHardAfter
          
subroutine _AddHardcpyAfter(node, o)
!Insert a new node after `node`
  class(_NODE)               :: node
  _CLASS, target, intent(in)  :: o
  type(_NODE), pointer       :: new

  allocate(new)
  allocate(new%o,source=o)
  new%hard = .true.
  new%next  => node%next
  node%next => new
end subroutine _AddHardcpyAfter

function _AddSoftHead(node, o) result(new)
!Insert a new head before node
!On output, node is the new head
  class(_NODE),target       :: node
  _CLASS, target, intent(in) :: o
  type(_NODE), pointer      :: new

  allocate(new)
  new%o  =>  o
  new%hard = .false.
  new%next => node
end function _AddSoftHead

function _AddHardHead(node, o) result(new)
!Insert a new head before node
!On output, node is the new head
  class(_NODE),target       :: node
  _CLASS, target, intent(in) :: o
  type(_NODE), pointer      :: new

  allocate(new)
  allocate(new%o,source=o)
  new%next => node
  new%hard = .true.
end function _AddHardHead
      
subroutine _Deattach( node, prev )
! WARNING: Do not forget to deallocate the node object
class(_NODE)            :: node
type(_NODE)             :: prev
prev%next => node%next
end subroutine _Deattach
 
subroutine _DestroyNode(node)
! Destroy the node. You might want deattach the node first
! to avoid break the list integrity.
! WARNING: Do not forget to deallocate the node
class(_NODE)    :: node

if(node%hard) then
  if (associated(node%o)) deallocate(node%o)
else
  node%o => null()
endif

! Set variables back to default... nah. 
! node % next => null()
! node % hard = .false.
     
end subroutine _DestroyNode
           
subroutine _DestroyAll( head )
class(_NODE),target    :: head
class(_NODE),pointer   :: aux

do while(associated(head%next))
  aux => head%next
  head%next => aux%next
  call aux%destroy_node()
  deallocate(aux)
enddo
call head%destroy_node()

end subroutine _DestroyAll

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

