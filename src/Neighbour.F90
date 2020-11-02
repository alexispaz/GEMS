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

 
module gems_neighbour

#define SMALL 0.0001

! This module compute the neigbour list of each interaction group. If the
! interaction group is `ig`, the neighboor list to construct is given by
! `ig%nn(:)` (number of neighboors of each atom) and `if%list(:,:)` (index of
!     neighboors atoms).

! There are different flavors to search and use the neigboor lists:

! self: The neighbors of atom list `a` are search in atom list `a` 
!  - cross: The neighbors of atom list `a` are search in atom list `b`
! half: the nehigboor j is in the of i, but not in the otherway
!  - full: the nehigboor j is in the of i and viceversa
! newton: when the force is computed for i, it is added also to the neighboor
! crossmol: intramolecular neighboors are not taked into account

! half usually needs newton
! full usually needs newton off
! openmp is not thread safe with self and half (and newton)
! openmp is thread safe with self and full (and newton off)
! openmp is thread safe with cross
!
! TODO: la subrrutina cross y la self son muy parecidas, se podrían unir?
!
!neigh (a,b,c) = the "atoms" index of neigh:
! a- the "atoms" index of atom
! b- the "subs" index of neighbour
! c- the count of neigh atom
!neigh (a,b,c) = the number of neigh:
! a- the "atoms" index of atom
! b- the "subs" index of neighbour
! c- 2 for neighs in core (potential cut radio), 1 for all (core + shell)
!the spirit of core-shell is update neighbour list when a shell neigh particle live the shell to the core. this event can check in the neigh loop of the force calculation.
use gems_algebra,  only: real_v, integer_v
use gems_program_types, only: boxed,box
use gems_program_types
use gems_groups
use gems_atoms
use gems_constants,     only: dp,cdm,dm,ui_ev
use gems_inq_properties, only: vdistance
use gems_errors
! use gems_algebra

implicit none

! XXX?
save
! private

public :: polvar_neighbour, polvar_intergroup

integer,parameter              :: mnb=1000             ! Maximum number of neighboor per atom


! From 1 to 13 are the upper cells. From 14 to 26 the lower. 0 is the center
integer,parameter :: map(3,0:26) = &
         reshape( [ 0, 0, 0,  &
                    1, 0, 0,   1, 1, 0,   0, 1, 0,  -1, 1, 0,&
                    1, 0,-1,   1, 1,-1,   0, 1,-1,  -1, 1,-1,&
                    1, 0, 1,   1, 1, 1,   0, 1, 1,  -1, 1, 1,&
                    0, 0, 1,  -1, 0, 0,  -1,-1, 0,   0,-1, 0,&
                    1,-1, 0,  -1, 0, 1,  -1,-1, 1,   0,-1, 1,&
                    1,-1, 1,  -1, 0,-1,  -1,-1,-1,   0,-1,-1,&
                    1,-1,-1,   0, 0,-1],[3,27])

 
! Mark the first interaction call to deny new interact commands.
! Currently, I need to perform pbcghost and update routines after all the
! interactions are set and after that any new interaction is not alowed.
! Thus, after nomoreigs is set true, no more interaction groups can be created.
! TODO: May be there is a way to change this to be incremental
logical     :: nomoreigs=.false.
 
type,public :: intergroup
  
  ! Introduzco el factor 14 por la ventana. Deberia modificar esto.
  real(dp)                      :: fac14=1.0_dp
   
  ! In some case, is better if atom indexes follows a consecutive sequence with numbers internal to the intergroup, to avoid large
  ! allocation in the internal routine matrices. Thus, if there is a change in the number of atoms in the intergroup, the numbers
  ! should be updated. Thus, each time the number of local or ghost atoms change, the internal ordering should be updated by sorting
  ! the atoms again. This sort should be an order N operation, so is not too expensive. The follwoing array of pointers introduce
  ! the described internal index to the atoms that belong to the intergroup.
  !
  ! In order to avoid distinguish between a and b list in the force field
  ! routines, the internal atom index must be unique for each atom, regardless
  ! is in a or b list. Thus if the neighbor of an atom is index 10, it is not
  ! needed to know to which list belong to compute the interaction.
  ! - 1st index: list of atoms involved in the interaction (1 or 2).
  ! - 2nd index: internal intergroup atom index: 
  !   From 1:n(1)       local atoms of list a
  !   From n(1)+1:n(2)  ghost atoms of list a
  !   From n(2)+1:n(3)  local atoms of list b
  !   From n(3)+1:n(4)  ghost atoms of list b
  type(atom_ap), pointer        :: at(:)
                      
  ! Ranges for do loops
  integer                       :: n(4)=0

  ! Neihgboor list build in this module and used in any forcefield routine.
  integer,allocatable           :: nn(:)        ! Numero de vecinos al atomo i
  integer,allocatable           :: list(:,:)    ! Vecino m-th del i-th atom of the first iteracting group (local and ghost a)

  type(atom_dclist),pointer     :: a,b,ag,bg              ! The interacting groups

  ! Index of first atom in cell for a and b list (using internal intergroup index)
  integer,allocatable           :: ahead(:,:,:),bhead(:,:,:)
   
  ! Cells
  real(dp)                      :: cmax(3),cmin(3)
  integer                       :: ncell=1
  integer                       :: ncells(3)=[1,1,1]
  integer                       :: cellaux(3)=[1,1,1]
  real(dp)                      :: cell(3)
   
  ! Index of next atom in cell (using internal intergroup index)
  integer,allocatable           :: clist(:)      
                                                    
  ! Radio de corte
  real(dp)                      :: rcut=1.e10_dp,rcut2=1.e10_dp
  
  ! Energia potencial
  real(dp)                      :: epot
  
  ! Parametros del potencial
  type(real_v)                  :: p
  type(integer_v)               :: i
  
  ! The index of the igr_vop vector. It is mailny used to
  ! reconect atoms that migrate to different cells to their corresponding
  ! intergourp. It is also used to select the interaction (for biasing, for
  !     example)
  integer :: id=0     

  ! True if is self interaction            
  logical :: self=.false.
              
  ! Control if the intramolecular interaction should be excluded.
  logical :: crossmol=.false.

  ! Add reaction forces to the neighbors. This have only sense in potentials that support cross interaction. For instance, setting
  ! this to .false. for LJ it might have sense if you do not want that `b` group "feels" `a` (non conservative). However, it is has
  ! not sense for TB, EAM, Analiticals or any other potential that does not have a `b` group or does not have the cross interaction
  ! capability.
#if _OPENMP  
  logical :: newton=.false.
#else  
  logical :: newton=.true.
#endif
                  
  ! if half, the nehigboor j is in the of i, but not in the otherway
#if _OPENMP  
  logical :: half=.false.
#else  
  logical :: half=.true.
#endif
                    
  ! This boolean allows to skip all the interaction.
  ! From the interact routine. This allows to compute
  ! forces on selected interactions
  logical :: disable=.false.
            
  ! If true, search also for neighboors of ghost atoms. Not to be confused with global variable fullghost of Program_Types which is
  ! used to define the ghost needed by pbc if not mic.
  logical :: fullghost=.false.
                   
  integer             :: nnb_cell

  procedure(intergroup0),pointer :: lista=>null()  ! funcion propiamente dicha
  procedure(intergroup0),pointer :: interact=>null()  ! funcion de interaccion

  contains
    procedure   :: ensure_alloc => intergroup_ensurealloc
    procedure   :: init => intergroup_constructor
    procedure   :: intergroup_adda_group
    procedure   :: intergroup_addb_group
    procedure   :: intergroup_adda_atom
    procedure   :: intergroup_addb_atom
    generic     :: adda => intergroup_adda_group,intergroup_adda_atom
    generic     :: addb => intergroup_addb_group,intergroup_addb_atom
    procedure   :: cleana => intergroup_cleana
    procedure   :: cleanb => intergroup_cleanb
    procedure   :: setrc => intergroup_setrc
    procedure   :: setcells => intergroup_setcells
endtype

abstract interface
 subroutine intergroup0(this)
  import intergroup
  class(intergroup),intent(inout)  :: this 
 end subroutine
end interface
     
! Double linked list for intergroups used for modules that require to keep track of their intergroups 
! (e.g. TB or EAM to perform the preinteraction)
#define _NODE intergroup_dl
#define _CLASS class(intergroup)
#include "dlist_header.inc"

#define _NODE intergroup_aop
#define _CLASS class(intergroup)
#include "arrayofptrs_header.inc"

#define _NODE intergroup_vop
#define _TYPE type(intergroup_aop)
#include "vector_header.inc"

! VOP to collect the intergroups declared inside modules
! The order of execution is important! See
! bias interaction.
type(intergroup_vop),public :: igr_vop
                 
real(dp),target     :: nb_dcut=1._dp        ! The shell length for verlet update criteria
integer             :: nupd_vlist = 0       ! Number of updates to print in the log
real(dp)            :: maxrcut=0.0_dp       ! Maximum cut ratio


contains
 
#define _NODE intergroup_dl
#define _CLASS class(intergroup)
#include "dlist_body.inc"
                         
#define _NODE intergroup_vop
#define _TYPE type(intergroup_aop)
#include "vector_body.inc"
                          
! INTERGROUP PROCEDURES

subroutine intergroup_constructor(ig,rc,g1,g2)
! Si se contiene g2 pero ademas g1, agregar sin g2 y despues incluir atomos en
! g2 con addb
class (intergroup),target         :: ig
real(dp),intent(in),optional      :: rc
type(group),intent(in)            :: g1
type(group),intent(in),optional   :: g2
integer                           :: n

! Check if the interact routine has been called, so avoid to construct new
! integroup
call werr('Command not allowed after an energy or force calculation.',nomoreigs)

! Add to igr_vop
call igr_vop%append()
n=igr_vop%size
igr_vop%o(n)%o=>ig
ig%id = n
      
! Inicializo los vectores de parametros
call ig%p%init()
call ig%i%init()
    
! Inicializo la lista con sus marcas
allocate(ig%a)
call ig%a%init()    
ig%ag=>ig%a
call ig%adda(g1)
ig%self=.true.

! radio de corte
if(present(rc)) then
  call ig%setrc(rc)
endif
 
if(present(g2)) then
  allocate(ig%b)
  call ig%b%init()    
  ig%bg=>ig%b
  call ig%addb(g2)
  ig%self=.false.
  ig%newton=.true.
endif

end subroutine intergroup_constructor
                       
subroutine intergroup_ensurealloc(ig)
! Add ghost or local atom to a group.
class ( intergroup ),target    :: ig
integer                        :: n1,n2
integer,parameter              :: spd=100

n1=ig%n(2)
n2=ig%n(4)

if(.not.allocated(ig%nn)) then
  allocate(ig%nn(n1+spd))
  allocate(ig%list(n1+spd,mnb))
  allocate(ig%at(n2+spd))
  allocate(ig%clist(n2+spd)) 
  return
endif

! TODO: set a flag when fullghost is needed to avoid excesive allocation here
! Neighboor list is constructed only for list a
if(size(ig%nn)<n1) then
  deallocate(ig%nn)
  deallocate(ig%list)
  allocate(ig%nn(n1+spd))
  allocate(ig%list(n1+spd,mnb))
endif  

! Both list a and b should be sort in cells
if(size(ig%at)<n2) then
  deallocate(ig%at)
  deallocate(ig%clist)
  allocate(ig%at(n2+spd))
  allocate(ig%clist(n2+spd)) 
endif  
 
! ! TODO: Check if neighboor of ghost atoms are to be considered
! n=g%n(1)
! if(g%fullghost) n=g%n(2)
         
end subroutine intergroup_ensurealloc
                           
subroutine intergroup_adda_group ( ig,g,ghost)
! Add ghost or local atom to a group.
class ( intergroup ),target       :: ig
type(group),intent(in)            :: g
logical,intent(in),optional       :: ghost
type(atom_dclist),pointer         :: la
integer                           :: i

la => g%alist
do i = 1,g%nat
  la => la%next

  if(present(ghost)) then
    if (ghost) then
      ! Add the atom at the end of the ghost
      call ig%a%add_before_soft(la%o)
      ig%n(2:4)=ig%n(2:4)+1
    endif
  else
    ! Add the atom just at the end of the local and the begining of the ghost
    call ig%ag%add_soft(la%o)
    ig%ag=>ig%ag%next
    ig%n(:)=ig%n(:)+1
  endif

  ! Set the integroup flags of the atom
  call la%o%addigr(ig%id,.true.)

enddo  

call ig%ensure_alloc()

end subroutine intergroup_adda_group
                 
subroutine intergroup_adda_atom ( ig,o,ghost)
! Add ghost or local atom to a group.
class ( intergroup ),target       :: ig
type(atom),pointer,intent(in)     :: o
logical,intent(in),optional       :: ghost

if(present(ghost)) then
  if (ghost) then
    ! Add the atom at the end of the ghost
    call ig%a%add_before_soft(o)
    ig%n(2:4)=ig%n(2:4)+1
  endif
else
  ! Add the atom just at the end of the local and the begining of the ghost
  call ig%ag%add_soft(o)
  ig%ag=>ig%ag%next
  ig%n(:)=ig%n(:)+1
endif

! Set the integroup flags of the atom
call o%addigr(ig%id,.true.)

call ig%ensure_alloc()

end subroutine intergroup_adda_atom
                 
subroutine intergroup_addb_group ( ig,g,ghost)
! Add ghost or local atom to a group.
class ( intergroup ),target       :: ig
type(group),intent(in)            :: g
logical,intent(in),optional       :: ghost
type(atom_dclist),pointer         :: la
integer                           :: i

la => g%alist
do i = 1,g%nat
  la => la%next

  if(present(ghost)) then
    if (ghost) then
      ! Add the atom at the end of the ghost
      call ig%b%add_before_soft(la%o)
      ig%n(4)=ig%n(4)+1
    endif
  else
    ! Add the atom just at the end of the local and the begining of the ghost
    call ig%bg%add_soft(la%o)
    ig%bg=>ig%bg%next
    ig%n(3:4)=ig%n(3:4)+1
  endif
          
  ! Set the integroup flags of the atom
  call la%o%addigr(ig%id,.false.)
              
enddo     

call ig%ensure_alloc()

end subroutine intergroup_addb_group
                 
subroutine intergroup_addb_atom ( ig,o,ghost)
! Add ghost or local atom to a group.
class ( intergroup ),target       :: ig
type(atom),pointer,intent(in)     :: o
logical,intent(in),optional       :: ghost

if(present(ghost)) then
  if (ghost) then
    ! Add the atom at the end of the ghost
    call ig%b%add_before_soft(o)
    ig%n(4)=ig%n(4)+1
  endif
else
  ! Add the atom just at the end of the local and the begining of the ghost
  call ig%bg%add_soft(o)
  ig%bg=>ig%bg%next
  ig%n(3:4)=ig%n(3:4)+1
endif

! Set the integroup flags of the atom
call o%addigr(ig%id,.true.)

call ig%ensure_alloc()

end subroutine intergroup_addb_atom
                 
subroutine intergroup_cleana (g)
! Asociate `at` in the a part and clean deleted atoms
class ( intergroup ),target       :: g
type ( atom_dclist ), pointer     :: la, prev
type ( atom ), pointer            :: o
integer                           :: i,j,d

! Acomodo en at() los atomos locales
la => g%a
do i = 1,g%n(1)
  la=>la%next
  g%at(i)%o=>la%o
enddo
  
! Acomodo en at() los ghost
d=0
do i = g%n(1)+1,g%n(2)
  la=>la%next
  o=>la%o

  if(o%delete) then
    
    d=d+1

    ! Delete the intergroup from this atom
    call o%deligr(g%id)

    ! Atempt to destroy the atom
    call del_atom(o)

    ! Delete the link
    prev => la%prev
    call la%deattach()
    call la%destroy_node() !soft

    deallocate(la)
    la=>prev
       
  else
    j=i-d
    g%at(j)%o=>la%o
  endif                

enddo

! Update index breaks
g%n(2:4)=g%n(2:4)-d

end subroutine intergroup_cleana
                 
subroutine intergroup_cleanb (g)
! Asociate `at` in the a part and clean deleted atoms
class ( intergroup ),target       :: g
type ( atom_dclist ), pointer     :: la, prev
type ( atom ), pointer            :: o
integer                           :: i,j,d
        
la => g%b
do i = g%n(2)+1,g%n(3)
  la=>la%next
  g%at(i)%o=>la%o
enddo
 
d=0
do i = g%n(3)+1,g%n(4)
  la=>la%next
  o=>la%o

  if(la%o%delete) then
    
    d=d+1

    ! Delete the intergroup from this atom
    call o%deligr(g%id)

    ! Atempt to destroy the atom
    call del_atom(o)

    ! Delete the link
    prev => la%prev
    call la%deattach()
    call la%destroy_node() !soft
    deallocate(la)
    la=>prev
                               
  else
    j=i-d
    g%at(j)%o=>la%o
  endif
enddo
  
! Update index breaks
g%n(4)=g%n(4)-d
 
end subroutine intergroup_cleanb
                 
subroutine intergroup_setrc (ig,rc)
class(intergroup)   :: ig
real(dp),intent(in) :: rc

! New cut radious
ig%rcut = rc
ig%rcut2 = rc*rc
maxrcut=max(maxrcut,rc)

end subroutine intergroup_setrc
 
subroutine intergroup_setcells(ig)
class(intergroup) :: ig
integer,parameter :: mincells=60 ! (4x4x4)

if(associated(ig%lista,target=intergroup0_empty)) return

if(.not.boxed) then
  if(ig%self) then
    if(ig%half) then
      ig%lista=>verlet_selfhalf
    else
      ig%lista=>verlet_self
    endif
  else
    ig%lista=>verlet_cross
  endif
  call wlog('NHB','Not using linked cells since box is not defined.')
  call wlog('NHB'); write(logunit,fmt="(a,f10.5)") ' -cut radious: ', ig%rcut
  return
endif 

! Redimension of ncells
ig%ncells(:)=int(box(:)/(ig%rcut+nb_dcut)) ! Number of cells in each direction 
   
! Less than 4 have no sense to use cells
if(all(ig%ncells(:)<4)) then
  if(ig%self) then
    if(ig%half) then
      ig%lista=>verlet_selfhalf
    else
      ig%lista=>verlet_self
    endif 
  else
    ig%lista=>verlet_cross
  endif
  call wlog('NHB','Not using linked cells since it would involve less than 4 cells.')
  call wlog('NHB'); write(logunit,fmt="(a,f10.5)") ' -cut radious: ', ig%rcut
  return
endif

ig%ncell = ig%ncells(1)*ig%ncells(2)*ig%ncells(3)

if(ig%self) then
  if(ig%half) then
    ig%lista=>linkedcell_selfhalf
  else
    ig%lista=>linkedcell_self
  endif 
else
  ig%lista=>linkedcell_cross
endif

! Cell size
ig%cell(:)=box(:)/ig%ncells(:)             
          
! Reallocate head array
if(allocated(ig%ahead)) deallocate(ig%ahead)
if(allocated(ig%bhead)) deallocate(ig%bhead)
if(ghost) then
  call werr('WARNING. GHOSTS AND LINKED CELLS NEVER TEST BEFORE')
  allocate(ig%ahead(0:ig%ncells(1)+1,0:ig%ncells(2)+1,0:ig%ncells(3)+1))
  allocate(ig%bhead(0:ig%ncells(1)+1,0:ig%ncells(2)+1,0:ig%ncells(3)+1)) 
else
  allocate(ig%ahead(ig%ncells(1),ig%ncells(2),ig%ncells(3)))
  allocate(ig%bhead(ig%ncells(1),ig%ncells(2),ig%ncells(3))) 
endif
 
! Used to compute cell index
ig%cellaux(1) = 1
ig%cellaux(2) = ig%ncells(1)
ig%cellaux(3) = ig%ncells(1)*ig%ncells(2)

call wlog('NHB','Using Linked Cells.')
call wlog('NHB'); write(logunit,fmt="(a,f10.5)") ' cut radious: ', ig%rcut
call wlog('NHB'); write(logunit,fmt="(a,"//cdm//"(f10.5,2x))") ' cell size: ', ig%cell(1:dm)
call wlog('NHB'); write(logunit,fmt="(a,"//cdm//"(i3,2x))") ' cell numbers: ', ig%ncells(1:dm)
 
end subroutine intergroup_setcells


! Sin lista

subroutine intergroup0_empty(g)
! FIXME???
class(intergroup),intent(inout)       :: g
end subroutine

! Lista de verlet.

subroutine verlet_self(g)
! Build neighbors list for atoms in `a` list by searching in the same `a` list.
class(intergroup),intent(inout)       :: g
type(atom_dclist),pointer            :: la,lb
integer                              :: i,j,m,n,k
real(dp)                             :: rd,vd(dm)
real(dp)                             :: rcut

! Set ceros (XXX: I think this is not needed)
g%nn(:)=0

! Cut radious
rcut=(g%rcut+nb_dcut)
rcut=rcut*rcut

! Construct array of atom pointers `at` and delete atoms flagged to be deleted. 
! XXX: Instead of make an array of atom pointers, it might by worth to create an
! array of atoms (considering a less complex atoms class with only pos property)
! to take advantage of cache when using the entire array. Furthermore, it would
! be better if atoms are also ordered in some way by proximity (e.g. in x
! coordinate). 
call g%cleana()

! Check if neighboor of ghost atoms are to be considered
n=g%n(1)
if(g%fullghost) n=g%n(2)
    
! XXX: Otra opcion para evaluar sería
! !$OMP PARALLEL DO SCHEDULE(STATIC,5) PRIVATE(m,vd,rd)
! do i = 1,g%na+g%nag-1
! m = 0
! do j = i+1,g%na+g%nag
!   vd = g%at(j)%o%pos-g%at(i)%o%pos 


!OMP: Creo varios threads
!$OMP PARALLEL 

!OMP: Un thread ejecuta, los demas quedan idle.
!$OMP SINGLE
    
la => g%a  
do i = 1,n
  la=>la%next
 
  !OMP: Se reparte la tarea entre los idle threads
  !$OMP TASK DEFAULT(NONE)    &
  !$OMP& FIRSTPRIVATE(la,i)   &
  !$OMP& SHARED(g,rcut,mic)   &
  !$OMP& PRIVATE(lb,j,vd,rd,m,k)
 
  m=0
  lb => g%a
  do j = 1 , g%n(2)
    lb=>lb%next
   
    ! Skip autointeraction
    if(i==j) cycle
   
    ! Skip intramolecular interaction.
    if(g%crossmol) then
      if(la%o%molid==lb%o%molid) cycle
    endif
  
    vd = vdistance( la%o, lb%o , mic)
    rd = dot_product(vd,vd) 

    if (rd>rcut) cycle

    ! Add j as neighboor of i.
    m=m+1
    g%list(i,m)=j

  enddo 
  g%nn(i)=m

  !$OMP END TASK
enddo  

! if(ghost) g%nn(g%n(2))=0
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL

end subroutine

subroutine verlet_selfhalf(g)
! Build neighbors list for atoms in `a` list by searching in the same `a` list.
class(intergroup),intent(inout)       :: g
type(atom_dclist),pointer            :: la,lb
integer                              :: i,j,m,n,k
real(dp)                             :: rd,vd(dm)
real(dp)                             :: rcut

! Set ceros (XXX: I think this is not needed)
g%nn(:)=0

! Cut radious
rcut=(g%rcut+nb_dcut)
rcut=rcut*rcut

call g%cleana()

! Check if neighboor of ghost atoms are to be considered
n=g%n(1)
if(g%fullghost) n=g%n(2)
    
la => g%a  
do i = 1,n
  la=>la%next
 
  m=0
  lb => la
  do j = i+1 , g%n(2)

    lb=>lb%next
   
    ! Skip intramolecular interaction.
    if(g%crossmol) then
      if(la%o%molid==lb%o%molid) cycle
    endif
  
    vd = vdistance( la%o, lb%o , mic)
    rd = dot_product(vd,vd) 

    if (rd>rcut) cycle

    ! Add j as neighboor of i.
    m=m+1
    g%list(i,m)=j

  enddo 
  g%nn(i)=m

enddo  

end subroutine

subroutine verlet_cross(g)
! Build neighbors list for atoms in `a` list by searching in `b` list.
class(intergroup),intent(inout)       :: g
type(atom_dclist),pointer            :: la,lb
integer                              :: i,j,m,k
real(dp)                             :: rd,vd(dm)
real(dp)                             :: rcut

! Set ceros (XXX: I think this is not needed)
g%nn(:)=0

! Cut radious
rcut=(g%rcut+nb_dcut)
rcut=rcut*rcut

! Construct array of atom pointers `at` and delete atoms flagged to be deleted. 
! XXX: Instead of make an array of atom pointers, it might by worth to create an
! array of atoms (considering a less complex atoms class with only pos property)
! to take advantage of cache when using the entire array. Furthermore, it would
! be better if atoms are also ordered in some way by proximity (e.g. in x
! coordinate).
call g%cleana()
call g%cleanb()
               
!OMP: Creo varios threads
!$OMP PARALLEL 

!OMP: Un thread ejecuta, los demas quedan idle.
!$OMP SINGLE
         
! Search for neighbors of local atoms of `a`.
la => g%a
do i = 1,g%n(1)
  la=>la%next
      
  !OMP: Se reparte la tarea entre los idle threads
  !$OMP TASK DEFAULT(NONE)     &
  !$OMP& FIRSTPRIVATE(la,i)    &
  !$OMP& SHARED(g,rcut,mic)    &
  !$OMP& PRIVATE(lb,j,vd,rd,m,k)
                 
  m=0
  lb => g%b
  do j = g%n(2)+1,g%n(4)
    lb => lb%next
    
    ! No need to check for inter or intramolecular interaction, since neighboors
    ! are in another list.

    vd = vdistance( la%o, lb%o , mic)
    rd =  dot_product(vd,vd) 
    
    if (rd>rcut) cycle
     
    ! Add j as neighboor of i.
    m=m+1
    g%list(i,m)=j

  enddo 
  g%nn(i)=m
     
  !$OMP END TASK
enddo     
           
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL
          
! Search for neighbors of ghost atoms of `a`. This is done in a separate
! loop in order to exclude ghost-ghost interaction
if(.not.g%fullghost) return
                 
!$OMP PARALLEL 
!$OMP SINGLE
 
do i = g%n(1)+1,g%n(2)
  la=>la%next
       
  !OMP: Se reparte la tarea entre los idle threads
  !$OMP TASK DEFAULT(NONE)        &
  !$OMP& FIRSTPRIVATE(la,i)       &
  !$OMP& PRIVATE(lb,j,vd,rd,m,k)  &
  !$OMP& SHARED(g,rcut,mic) 
     
  m=0

  lb => g%b
  do j = g%n(2)+1,g%n(3)
    lb => lb%next
  
    ! No need to check for inter or intramolecular interaction, since neighboors
    ! are in another list.
                   
    vd = vdistance( la%o, lb%o , mic)
    rd =  dot_product(vd,vd) 
    
    if (rd>rcut) cycle

    m=m+1
    g%list(i,m)=j
 
  enddo 
  g%nn(i)=m
        
  !$OMP END TASK
      
enddo     
          
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL
         
end subroutine
          
! Lista por celdas

subroutine sort_in_cells(g)
! Los atomos son distribuidos en las celdas, el primero es puesto en el vector
! g%head. Los atomos siguientes son puestos en el g%clist. Esta subrrutina
! debe llamarse siempre antes del calculo de las fuerzas Las variables
! extrañas son seteadas en la subrutina maps.
class(intergroup),intent(inout)  :: g
integer                          :: i,aux1(dm)
real(dp)                         :: one_cell(3)
type(atom),pointer               :: a

g%ahead(:,:,:) = 0
g%bhead(:,:,:) = 0
            
! Invers of cell size  
one_cell(:)=1._dp/g%cell(:)           ! Inversa de las dimensiones de cada celda
                                                               
! Atoms of list a
! !$OMP  PARALLEL DO DEFAULT(NONE) &
! !$OMP& PRIVATE(i,a,aux1) &
! !$OMP& SHARED(g,one_cell)
do i = 1,g%n(2)
  a=>g%at(i)%o
  
  ! FIXME
  ! ! Si no hago pbc podria salirse alguna fuera y dar un segfull
  ! where(a%pbc(:)) 
  !   a%pos(:)=mod(a%pos(:)+box(:),box(:))
  ! endwhere

  ! Get cell index
  aux1(:)=int(a%pos(:)*one_cell(:))+1
  ! j = 1 + dot_product(aux1,cellaux)

  ! Build cell list
  ! !$OMP CRITICAL
  g%clist(i) = g%ahead(aux1(1),aux1(2),aux1(3))
  g%ahead(aux1(1),aux1(2),aux1(3)) = i
  ! !$OMP END CRITICAL

enddo
! !$OMP END PARALLEL DO
   
! Atoms of list b
do i = g%n(2)+1,g%n(4)
  a=>g%at(i)%o

  ! FIXME
  ! ! Si no hago pbc podria salirse alguna fuera y dar un segfull
  ! where(a%pbc(:)) 
  !   a%pos(:)=mod(a%pos(:)+box(:),box(:))
  ! endwhere

  ! Get cell index
  aux1(:)=int(a%pos(:)*one_cell(:))+1
  ! j = 1 + dot_product(aux1,cellaux)
  !
  ! Build cell list
  g%clist(i) = g%bhead(aux1(1),aux1(2),aux1(3))
  g%bhead(aux1(1),aux1(2),aux1(3)) = i

enddo

endsubroutine sort_in_cells

subroutine linkedcell_self(g)
!!$ use omp_lib 
! Build neighbors list for atoms in `a` list by searching in the same `a` list.
class(intergroup),intent(inout)  :: g
integer                          :: i,j,ic
integer                          :: nabor
real(dp)                         :: rd,vd(dm)
real(dp)                         :: rcut
integer                          :: rc(dm),nc(dm)
                      
! Set ceros
g%nn(:)=0
 
! Construct array of atom pointers `at` and delete atoms flagged to be deleted. 
! XXX: Instead of make an array of atom pointers, it might by worth to create an
! array of atoms (considering a less complex atoms class with only pos property)
! to take advantage of cache when using the entire array. Furthermore, it would
! be better if atoms are also ordered in some way by proximity (e.g. in x
! coordinate).
call g%cleana()
            
! Cut radious
rcut=(g%rcut+nb_dcut)
rcut=rcut*rcut
                
call sort_in_cells(g)

! Restrict to local cells if ghost exists
! FIXME
! if(ghost) then
!   if(g%fullghost) then
!     icell(:)=0
!     fcell(:)=g%ncells(:)+1
!   endif
! else
!   icell(:)=1
!   fcell(:)=g%ncells(:)
! endif
 
!Bucle sobre las celdas
!$OMP  PARALLEL DO DEFAULT(NONE) &
!$OMP& PRIVATE(rc,ic,i,j,vd,rd,nabor,nc) &
!$OMP& SHARED(g,rcut,mic)
do ic = 1,g%ncell

  !print *, omp_get_thread_num(), ic
  rc(:)=icell2vcell(g,ic)

  !Bucle sobre las celdas vecinas superiores a rc(:)
  do nabor=0,26
    nc(:)=map(:,nabor)+rc(:)

    ! Check if cells should be considered and apply PBC
    if(.not.cell_pbc(g,nc,mic)) cycle
  
    ! Bucle sobre los atomos de g%a en la celda rc(:)
    i= g%ahead(rc(1),rc(2),rc(3))
    do while( i>0 ) 

      !Bucle sobre los atomos de la celda vecina
      j = g%ahead(nc(1),nc(2),nc(3))

      do while( j > 0 ) 
         
        ! Skip autointeraction
        if(i==j) then
          j = g%clist(j)
          cycle
        endif
               
        ! Skip intramolecular interaction.
        if(g%crossmol) then
          if(g%at(i)%o%molid==g%at(j)%o%molid) then
            j = g%clist(j)
            cycle
          endif
        endif
               
        vd = vdistance( g%at(j)%o, g%at(i)%o , mic) ! respetar el orden
        rd =  dot_product(vd,vd) 
                               
        if (rd<rcut) then
          g%nn(i)=g%nn(i)+1
          g%list(i,g%nn(i))=j
        endif
      
        j = g%clist(j)
      enddo

      i = g%clist(i) 

    enddo
  enddo
enddo
!$OMP END PARALLEL DO
     
end subroutine

subroutine linkedcell_selfhalf(g)
!!$ use omp_lib 
! Build neighbors list for atoms in `a` list by searching in the same `a` list.
class(intergroup),intent(inout)  :: g
integer                          :: i,j,ic
integer                          :: nabor
real(dp)                         :: rd,vd(dm)
real(dp)                         :: rcut
integer                          :: rc(dm),nc(dm)
                      
! Set ceros
g%nn(:)=0

call g%cleana()
            
! Cut radious
rcut=(g%rcut+nb_dcut)
rcut=rcut*rcut
                
call sort_in_cells(g)

!Bucle sobre las celdas
do ic = 1,g%ncell

  !print *, omp_get_thread_num(), ic
  rc(:)=icell2vcell(g,ic)

  !Bucle sobre las celdas vecinas superiores a rc(:)
  do nabor=0,13
    nc(:)=map(:,nabor)+rc(:)

    ! Check if cells should be considered and apply PBC
    if(.not.cell_pbc(g,nc,mic)) cycle
  
    ! Bucle sobre los atomos de g%a en la celda rc(:)
    i= g%ahead(rc(1),rc(2),rc(3))
    do while( i>0 ) 

      !Bucle sobre los atomos de la celda vecina
      if(nabor==0) then
        j=g%clist(i)
      else
        j = g%ahead(nc(1),nc(2),nc(3))
      endif

      do while( j > 0 ) 
              
        ! Skip intramolecular interaction.
        if(g%crossmol) then
          if(g%at(i)%o%molid==g%at(j)%o%molid) then
            j = g%clist(j)
            cycle
          endif
        endif
               
        vd = vdistance( g%at(j)%o, g%at(i)%o , mic) ! respetar el orden
        rd =  dot_product(vd,vd) 
                               
        if (rd<rcut) then
          g%nn(i)=g%nn(i)+1
          g%list(i,g%nn(i))=j
        endif
      
        j = g%clist(j)
      enddo

      i = g%clist(i) 

    enddo
  enddo
enddo
     
end subroutine

subroutine linkedcell_cross(g)
! Ordena en g%nn y g%list los vecinos ordenados por celdas en g%alist y g%ahead
class(intergroup),intent(inout)  :: g
integer                          :: i,j,ic
integer                          :: nabor
real(dp)                         :: rd,vd(dm)
real(dp)                         :: rcut
integer                          :: rc(dm),nc(dm)
                
! Set ceros
g%nn(:)=0
                  
! Construct array of atom pointers `at` and delete atoms flagged to be deleted. 
! XXX: Instead of make an array of atom pointers, it might by worth to create an
! array of atoms (considering a less complex atoms class with only pos property)
! to take advantage of cache when using the entire array. Furthermore, it would
! be better if atoms are also ordered in some way by proximity (e.g. in x
! coordinate).
call g%cleana()
call g%cleanb()
    
! Cut radious
rcut=(g%rcut+nb_dcut)
rcut=rcut*rcut

call sort_in_cells(g)

! Restrict to local cells if ghost exists
! FIXME
! if(ghost) then
!   if(g%fullghost) then
!     icell(:)=0
!     fcell(:)=g%ncells(:)+1
!   endif
! else
!   icell(:)=1
!   fcell(:)=g%ncells(:)
! endif
 
!Bucle sobre las celdas
do ic = 1,g%ncell

  rc(:)=icell2vcell(g,ic)
  
  !Bucle sobre las celdas vecinas a rc(:)
  do nabor=0,26
             
    nc(:)=map(:,nabor)+rc(:)

    ! Check if cells should be considered and apply PBC
    if(.not.cell_pbc(g,nc,mic)) cycle
                                     
    ! Bucle sobre los atomos de g%a en la celda rc(:)
    i= g%ahead(rc(1),rc(2),rc(3))
    do while( i>0 ) 
             
      !Bucle sobre los atomos de la celda vecina
      j = g%bhead(nc(1),nc(2),nc(3))
      do while( j > 0 ) 

        vd = vdistance( g%at(j)%o, g%at(i)%o , mic) ! respetar el orden
        rd =  dot_product(vd,vd) 
        if (rd<rcut) then
          g%nn(i)=g%nn(i)+1
          g%list(i,g%nn(i))=j
        endif
        j = g%clist(j)

      enddo
      i = g%clist(i) 

    enddo
  enddo
enddo
     
end subroutine

! Actualizacion de las listas

subroutine update()
integer                        :: i!,t1,t2
class(intergroup), pointer     :: ig
real(dp)                       :: rcut
              
! call system_clock(t1)
! print *, 'entra' 
nupd_vlist = nupd_vlist +1  

do i = 1, nlocal
  a(i)%o%pos_old = a(i)%o%pos
enddo

! Needed for NPT
if(boxed) box_old(:)=box(:)
   
! Simple loop  
! node => ilist
! do while( associated(node%next) )
!   node => node%next
!   ig => node%o
do i = 1, igr_vop%size
  ig => igr_vop%o(i)%o

  ! Check to skip
  if (ig%disable) cycle

  ! tepot(j)=node%obj%epot
  ! del sistema, cuando el sistema no esta en una caja. Esto funcionaria si
  ! modifico la caja del gems para que tenga un punto maximay un punto minimo
  ! y no que el punto minimo sea (0,0,0)
  ! if(.not.boxed) then
  !   call get_boundingbox(ig%alist,box_max,box_min,ig%n(2))
  !   call get_boundingbox(ig%blist,amax,amin,ig%n(4)-ig%n(2))
  !   box_max(:)=max(box_max(:),amax(:))
  !   box_min(:)=min(box_min(:),amin(:))
  !   set box
  ! endif
             
  if(boxed) then

    ! Redimension of cell size (needed for NPT)
    ig%cell(:)=box(:)/ig%ncells(:)

    rcut=ig%rcut+nb_dcut
      
    ! Check if the cell size drop below cut radious
    ! In that case change the number of cells, or deactivate linked cells.
    if(any(ig%cell(:)<rcut)) call ig%setcells()

    ! Check if it is possible to include more cells
    if(any(box(:)-ig%ncells(:)*rcut>rcut)) call ig%setcells()

  endif
   
  call ig%lista()
enddo 
           
! call system_clock(t2)
! print *,(t2-t1)/real(sclock_rate)
! stop
end subroutine update

subroutine test_update()
use gems_program_types, only:fullghost
!hace la lista de vecinos de ss y todos los subsistemas superiores
!notar que si ss es el ultimo, estamos seleccionando para ss
real(dp)            :: rd,dispmax1,dispmax2,vd(dm)
integer             :: i
 
! Update the ghost positions
if(ghost) call pbcghost_move
              

!vecinos con el propio sistema
dispmax1 = 1.e-16_dp
dispmax2 = 1.e-16_dp
      
do i = 1, nlocal

  vd = a(i)%o%pos - a(i)%o%pos_old

  rd = dot_product(vd,vd)

  if (rd>dispmax1) then
    dispmax2 = dispmax1
    dispmax1 = rd
  else
    if (rd>dispmax2) dispmax2 = rd
  endif

enddo

if (sqrt(dispmax1)+sqrt(dispmax2)>nb_dcut) then
  if(ghost) then
    if(fullghost) then
       call pbcfullghost()
    else
       call pbcghost()
    endif
  endif 
  call update()
endif

end subroutine test_update

#ifdef DIM3
        
function cell_pbc(g,r,mic)
! Perform PBC on cells. If cell should not be considered return .false.
integer            :: r(3)
class(intergroup),intent(in)  :: g
logical,intent(in) :: mic
logical            :: cell_pbc

cell_pbc=.true.

if (mic) then
  ! May use pbc (it will further depends on atom%pbc)
  r(:)=r(:)-1
  r(:)= mod( r(:) + g%ncells(:), g%ncells(:) )
  r(:)=r(:)+1
else
  ! Can not use PBC
  if(ghost) then
    if(any(r(:)<0)) cell_pbc=.false.
    if(any(r(:)>g%ncells(:)+1)) cell_pbc=.false.
  else
    if(any(r(:)<1)) cell_pbc=.false.
    if(any(r(:)>g%ncells(:))) cell_pbc=.false.
  endif
endif  


end function cell_pbc

function vcell2icell(g,i,j,k,mic) result(uid)
!Convert 3-index cell id into 1-index cell id.
!The cell 1,1,1 is map to the cell 1, the 2,1,1 to 2...
class(intergroup),intent(in)  :: g
integer,intent(in) :: i,j,k
integer            :: r(3)
integer            :: uid
logical,intent(in) :: mic

r(1)=i-1
r(2)=j-1
r(3)=k-1

if (mic) then
  ! May use pbc (it will further depends on atom%pbc)
  r(:) =  mod( r(:) + g%ncells(:), g%ncells(:) )
else
  ! Can not use PBC
  if(any(r(:)<1)) then
    uid=0
    return
  endif
  if(any(r(:)>g%ncells(:))) then
    uid=0
    return
  endif 
endif  

uid = 1 + dot_product(r,g%cellaux(:))

! XXX: Not sure why this (I guess never happend).
if (uid<0) uid=0

end function vcell2icell

function icell2vcell(g,ic) result(r)
!Convert 1-index cell id into 3-index cell id.
!ic shul be 1-based integer
!The cell 1,1,1 is map to the cell 1, the 2,1,1 to 2...
class(intergroup),intent(in)  :: g
integer,intent(in) :: ic
integer            :: r(3),c0
 
c0=ic-1

r(3)=int(c0/g%cellaux(3))

r(1)=mod(c0,g%cellaux(3))
r(2)=int(r(1)/g%cellaux(2))

r(1)=mod(r(1),g%cellaux(2))

r(:)=r(:)+1

end function icell2vcell
 
#endif


! GHOSTS

function islocal(r)
! Return .true. if r is inside the box
! It is equivalent to `isghost(r,0.)` 
real(dp),intent(in)  :: r(dm)
logical              :: islocal
integer              :: i

islocal=.false.
do i=1,dm
  if(r(i)<0._dp) then
    return
  elseif (r(i)>=box(i)) then
    return
  endif
enddo
islocal=.true.
            
end function islocal
  
function isghost(r,rc)
! Return .true. if r inside the box+rcut. So actually it says if some position
! belong to a ghost region but also if is a local region
real(dp),intent(in)  :: r(dm),rc
logical              :: isghost
integer              :: i

isghost=.false.
do i=1,dm
  ! s=sign(1,r(i)-box2(i))
  ! if (s*r(i)>rc+box*0.5_dp*(s+1)) return
           
  if(r(i)<-rc) then
    return
  elseif (r(i)>rc+box(i)) then
    return
  endif    
enddo
isghost=.true.
            
end function isghost
 
function isoldghost(r,rc)
! Return .true. if r inside the box+rcut. So actually it says if some position
! belong to a ghost region but also if is a local region
real(dp),intent(in)  :: r(dm),rc
logical              :: isoldghost
integer              :: i

isoldghost=.false.
do i=1,dm
  ! s=sign(1,r(i)-box2(i))
  ! if (s*r(i)>rc+box*0.5_dp*(s+1)) return
           
  if(r(i)<-rc) then
    return
  elseif (r(i)>rc+box_old(i)) then
    return
  endif    
enddo
isoldghost=.true.
            
end function isoldghost
 
function replicaisghost(n,r,rc)
! This expresion si general for any n(:) (neigh cells or not) and rc (bigger than
! box or not)
! This is the same that `call isghost(r+n*h,rc)`
real(dp),intent(in)  :: r(dm),rc
integer,intent(in)   :: n(dm)
logical              :: replicaisghost
integer              :: i,s

replicaisghost=.false.
do i=1,dm
  s=sign(1,n(i))
  if (s*r(i)>rc+box(i)*(0.5_dp*(s+1)-abs(n(i)))) return
enddo
replicaisghost=.true.

end function replicaisghost

subroutine pbcghost()
! Update ghost positions in a serial pbc simulation for only certain cut ratios.
! This should be call each time the neighboor list will be updated. This
! routine assume that the local atoms moves "slowly" between different calls. In
! other words,if there is a chance that local configurations used in two
! consecutive calls to pbcghost are uncorrelated, it is safer to use
! pbcfullghost.
use gems_atoms, only: atom_dclist
use gems_program_types, only: box,n1cells,nlocal,alocal,nghost,aghost
real(dp)                   :: rcut,r(dm),rold(dm)
type(atom_dclist), pointer :: la, prev
type(atom), pointer        :: o,o2
type(intergroup ), pointer :: ig
integer                    :: i,j,k,m
logical                    :: updatebcr
! integer                       :: ndel,nnew

rcut=maxrcut+nb_dcut

! ndel=0
! nnew=0

updatebcr=.false.
la => aghost
do i = 1,nghost
  la => la%next
  o => la%o
  
  ! Delete ghost atom from ghost list and flag it for deletion in the integroups
  if(.not.isghost(o%pos,rcut)) then
             
    ! ndel=ndel+1

    ! Attempt to destroy the atom
    call del_atom(o)
      
    ! Destroy ghost list link
    nghost=nghost-1
    prev => la%prev
    call la%deattach()
    call la%destroy_node() !soft
    deallocate(la)
    la=>prev
           
    cycle
  endif
      
  ! A ghost atom is inside the cell
  if(islocal(o%pos)) then
  
    updatebcr=.true.

    ! Fold its local image inside the box
    o2=>a(o%tag)%o
    do m=1,dm
      if(o2%pos(m)>=box(m)) then
        o2%pos(m)=o2%pos(m)-box(m)
        o2%pos_old(m)=o2%pos_old(m)-box(m)
        o2%boxcr(m)=o2%boxcr(m)+1
      elseif(o2%pos(m)<0.0_dp) then
        o2%pos(m)=o2%pos(m)+box(m)
        o2%pos_old(m)=o2%pos_old(m)+box(m)
        o2%boxcr(m)=o2%boxcr(m)-1
      endif
    enddo

    ! Pass the ghost atom to the opposite cell
    o%pos(:)=o2%pos(:)-o%boxcr(:)*box(:)

  endif
       
enddo  
     
! Find new ghost atoms
! !$OMP PARALLEL DO PRIVATE(m,i,j,la,o,k,o2,ig,r,rold)
do m =1,26
             
  la => alocal
  do i = 1,nlocal
    la => la%next
    o2 => la%o

    ! if (.not.all(la%o%pbc(:)*n1cells(m,:))) cycle

    r(:)=o2%pos(:)+n1cells(m,:)*box(:)
    if (isghost(r(:),rcut)) then

      ! Check if this was a ghost before
      rold(:)=o2%pos_old(:)+n1cells(m,:)*box_old(:)
      if (.not.isoldghost(rold(:),rcut)) then

        ! !$OMP CRITICAL

        ! nnew=nnew+1

        ! Un ghost se incorpora a la simulacion
        call new_ghostatom()
        o=>a(natoms)%o
        call atom_asign(o,o2)
        o%pos(:)=r(:)
        o%boxcr(:)=n1cells(m,:)
        o%q=o2%q
        o%tag=o2%tag

        ! Agrego el atomo a los intergroup correspondientes
        do j=1,o2%nigr
          k = o2%igr(j)
          ig=>igr_vop%o(k)%o

          if(o2%igra(j)) then
            call ig%adda(o,.true.)
          else
            call ig%addb(o,.true.)
          endif
        enddo  

        ! !$OMP END CRITICAL

      endif
    end if
  enddo
       
enddo
! !$OMP END PARALLEL DO

la => aghost
do i = 1,nghost
  la => la%next
  o => la%o

  if(o%delete) cycle
  
  ! Actualizar boxcr ya que si un local es plegado hacia adentro de la caja,
  ! muchos fantasmas van a tener el o%boxcr(:) incorrecto. Abria que ver el
  ! virial... FIXME
  if (updatebcr) o%boxcr(:)=floor(o%pos(:)*one_box(:))

  ! Save positions
  ! o%pos_old(:)=o%pos(:)
enddo  

end subroutine 
                   
subroutine pbcghost_move()
! Move ghost atoms to reflect the motion of their local images
use gems_program_types, only: box,nghost,aghost
use gems_atoms, only: atom_dclist
type ( atom_dclist ), pointer :: la
type ( atom ), pointer        :: o
integer                       :: i,j

la => aghost
do i = 1,nghost
  la => la%next
  o  => la%o
  
  j=o%tag
  o%pos(:)=a(j)%o%pos(:)+box(:)*o%boxcr(:)

enddo     

end subroutine 
                  
      
! subroutine ghost_pbc
!  use gems_program_types, only:box,atom_dclist,nlocal,nghost,aghost,alocal
!  integer                 :: ii,i,j
!  type ( atom_dclist ), pointer :: la,lb,lc
!  type ( atom ), pointer :: o
!            
!   la => aghost
!   do ii = 1,nghost
!     la => la%next
!        
!     if(all(la%o%pos(:)>0.0_dp)) then
!     if(all(la%o%pos(:)<box(:))) then
!       j = la%o%tag
!       i = la%o%id
!       lb=>a(j)%o%link
!       
!       ! print *, 'asd'
!
!       ! Ajusto alist y aghost
!       o=>la%o
!       la%o=>lb%o
!       lb%o=>o
!       
!       ! ! Ajusto a
!       ! a(i)%o=>a(j)%o
!       ! a(j)%o=>o
!       ! a(i)%id=i
!       ! a(j)%id=j
!       
!       ! Ajusto los links
!       lc=>lb%o%link
!       lb%o%link=>la%o%link
!       la%o%link=>lc
!
!       ! With half list, changeing the id is require a full update
!     endif
!     endif
!   enddo
!
! end subroutine

subroutine pbcfullghost()
! Update ghost positions in a serial pbc simulation for only certain cut ratios.
! This should be call each time the neighboor list will be updated. I think this
! routine is prefered over pbcghost when it can not be assume that the local
! atoms moves "slowly" between different calls. In other words, if there is a
! chance that local configurations used in two consecutive calls to pbcghost are
! uncorrelated, it is safer to use pbcfullghost.
use gems_atoms, only: atom_dclist
use gems_program_types, only: box,n1cells,nlocal,nghost,aghost, atom_dclist_destroyall
use gems_set_properties, only: do_pbc
real(dp)                      :: rcut,r(dm)
type ( atom_dclist ), pointer :: la
type ( atom ), pointer        :: o,o2
type ( intergroup ), pointer  :: ig
integer                       :: i,j,k,m
           
rcut=maxrcut+nb_dcut

! Destroy all ghost atoms
call atom_dclist_destroyall(aghost)
nghost=0

! Make local pbc
call do_pbc(alocal,nlocal)

! Find new ghost atoms
! !$OMP PARALLEL DO PRIVATE(m,i,j,la,o,k,o2,ig,r,rold)
do m =1,26
             
  la => alocal
  do i = 1,nlocal
    la => la%next
    o2 => la%o

    ! if (.not.all(la%o%pbc(:)*n1cells(m,:))) cycle

    r(:)=o2%pos(:)+n1cells(m,:)*box(:)
    if (isghost(r(:),rcut)) then

      ! !$OMP CRITICAL

      ! Un ghost se incorpora a la simulacion
      call new_ghostatom()
      o=>a(natoms)%o
      call atom_asign(o,o2)
      o%pos(:)=r(:)
      o%boxcr(:)=n1cells(m,:)
      o%q=o2%q
      o%tag=o2%tag

      ! Agrego el atomo a los intergroup correspondientes
      do j=1,o2%nigr
        k = o2%igr(j)
        ig=>igr_vop%o(k)%o

        if(o2%igra(j)) then
          call ig%adda(o,.true.)
        else
          call ig%addb(o,.true.)
        endif
      enddo  

      ! !$OMP END CRITICAL

    end if
  enddo
       
enddo
! !$OMP END PARALLEL DO

end subroutine 
     
function idoit(i,j,itag,jtag,vd)
! Use precomputed midpoint criterion to decide if interaction is owned.
logical                :: idoit
integer,intent(in)     :: i, j, itag, jtag
real(dp),intent(in)    :: vd(3)

idoit = .false.

if (i > nlocal) return
  
if (j < nlocal) then
  idoit = .true.
else if (itag < jtag) then
  idoit = .true.
else if (itag == jtag) then
  if (vd(3) > SMALL) then
    idoit = .true.
  else if (abs(vd(3)) < SMALL) then
    if (vd(2) > SMALL) then
      idoit = .true.
    else if (abs(vd(2)) < SMALL .and. vd(1) > SMALL) then
      idoit = .true.
    endif
  endif
endif

end function idoit

! Variables and Labels

function polvar_intergroup(var) result(g)
use gems_variables, only: polvar, polvar_find
use gems_errors, only: werr
character(*),intent(in)  :: var
type(polvar),pointer     :: pv
type(intergroup),pointer      :: g

call werr('Labels should start with colon `:` symbol',var(1:1)/=':')
pv=>polvar_find(var)

g=>null()
if(.not.associated(pv)) return

call werr('Variable not linked',pv%hard)
   
! Print
select type(v=>pv%val)
type is (intergroup)
  g=>v
class default
  call werr('I dont know how to return that')
end select

end function

subroutine polvar_neighbour()
use gems_variables, only: polvar_link
call polvar_link('nb_dcut',nb_dcut)
end subroutine
           
end module
