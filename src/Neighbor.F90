! Copyright (c) 2020  Sergio Alexis Paz
!
!  This file is part of GEMS. GEMS is an Extensible Molecular Simulator.
!   .
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


module gems_neighbor
use gems_program_types, only: boxed, box, mic
use gems_groups,        only: vdistance, sys, igroup, group, atom, atom_dclist
use gems_constants,     only: dp,cdm,dm,ui_ev
use gems_cells,         only: cgroup, map, n1cells, cell_pbc
use gems_errors

implicit none

public :: polvar_neighbor, polvar_ngroup

! Group with neighbors
! Forces are computed in group `ref` given group `b` (without newton reaction).
! ngroup gives a unique index to each atom in `ref` or `b`.
type, abstract, extends(igroup), public :: ngroup

  ! The reference group
  type(group)   :: ref

  ! The group of possible neighbors. It might be sorted in cells.
  type(cgroup)  :: b

  ! rcut
  real(dp)               :: rcut=1.e10_dp,rcut2=1.e10_dp

  ! Maximum number of neighbor per atom
  integer                :: mnb=10000

  ! Neighbor list build in this module and used in any forcefield routine.
  integer,allocatable    :: nn(:)        ! Numero de vecinos al atomo i
  integer,allocatable    :: list(:,:)    ! Vecino m-th del i-th atom of the first iteracting group (local and ghost a)

  ! Automatic switch between search algorithm
  logical :: autoswitch=.true.

  logical :: listed=.false.
  procedure(ngroup_cells),pointer :: lista => ngroup_cells
  procedure(ngroup_cells_atom),pointer :: lista_atom => ngroup_cells_atom

  contains

    ! Parent procedures that set a polymorphic pointer to an extension
    ! (see group type construct).
    procedure :: ngroup_construct
    procedure :: ngroup_attach_atom
    procedure :: ngroup_detach_atom

                                 
    procedure :: init => ngroup_construct
    procedure :: dest => ngroup_destroy
    procedure :: attach_atom => ngroup_attach_atom
    procedure :: detach_atom => ngroup_detach_atom


    ! Keep access to procedures of abstract parent
    ! see (https://fortran-lang.discourse.group/t/call-overridden-procedure-of-abstract-parent-type/590/23)
    procedure :: ngroup_init => ngroup_construct
    procedure :: ngroup_dest => ngroup_destroy

    procedure :: setrc => ngroup_setrc

    ! The calculation
    procedure(ngroup0),deferred :: interact

    ! The CLI used to read parameters
    procedure(ngroup0),nopass,deferred :: cli

endtype

abstract interface
 subroutine ngroup0(g)
  import ngroup
  class(ngroup),intent(inout)  :: g
 end subroutine
end interface

! Double linked list for ngroups used for modules that require to keep track of their ngroups
! (e.g. TB or EAM to perform the preinteraction)
#define SOFT
#define _NODE ngroup_dl
#define _CLASS class(ngroup)
#include "dlist_header.inc"

#define _NODE ngroup_aop
#define _CLASS class(ngroup)
#include "arrayofptrs_header.inc"

#define _NODE ngroup_vop
#define _TYPE type(ngroup_aop)
#include "vector_header.inc"

! VOP to collect the ngroups declared inside modules
! The order of execution is important! See
! bias interaction.
type(ngroup_vop),public :: ngindex

! Global variables
real(dp),target     :: nb_dcut=1._dp        ! The shell length for verlet update criteria
integer             :: nupd_vlist = 0       ! Number of updates to print in the log
real(dp)            :: maxrcut=0.0_dp       ! Maximum cut ratio

contains

#define SOFT
#define _NODE ngroup_dl
#define _CLASS class(ngroup)
#include "dlist_body.inc"

#define _NODE ngroup_vop
#define _TYPE type(ngroup_aop)
#include "vector_body.inc"
                 
! ngroup events
! =============

subroutine ngroup_construct(g)
class (ngroup),target             :: g

! Initialize
call g%igroup_construct()

! Ghost atoms must have its own index in order to comunicate back
! to the real image the computed properties
g%ghost=.true.

! Index the group
call ngindex%append()
ngindex%o(ngindex%size)%o=>g

! Init internal groups
call g%ref%init()
call g%b%init()
g%b%ghost=.true.

! Allocate
allocate(g%nn(g%pad))
allocate(g%list(g%pad,g%mnb))

end subroutine ngroup_construct

subroutine ngroup_destroy(g)
class (ngroup),target  :: g
integer                :: i

! Update ngindex of `g`
do i=1,ngindex%size
  if(ngindex%o(i)%o%id==g%id) exit
enddo
call ngindex%del(i,1)

! Destroy internal groups
call g%ref%dest()
call g%b%dest()

! Destroy
call g%igroup%dest()
deallocate(g%nn)
deallocate(g%list)

end subroutine ngroup_destroy

subroutine ngroup_attach_atom(g,a)
class(ngroup),target   :: g
class(atom),target     :: a
integer, allocatable   :: t_nn(:), t_list(:,:)
integer                :: n,i

! Attempt to attach
n=g%nat
call g%igroup_attach_atom(a)
if(n==g%nat) return

! Continue reallocations
! TODO: Implement a "preserve" boolean?
! TODO: Skip this if not allocated head?
n=size(g%a)
if(size(g%nn)<n) then
  allocate(t_nn(n))
  t_nn(:size(g%nn)) = g%nn(:)
  call move_alloc(to=g%nn,from=t_nn)

  allocate(t_list(n,g%mnb))
  t_list(:size(g%list,1),:) = g%list(:,:)
  call move_alloc(to=g%list,from=t_list)
endif   
 
! Sort atom into neighbor list. 
! This will only work if atom is already attached to `b` or `ref` before
! enter to this procedure. For example:
!   call g%ref%attach(a) 
!   call g%b%attach(a)  
!   call g%attach(a)  ! at last
if(g%listed) then
  i=a%gid(g)
  ! FIXME, g%listed puede ser true
  ! pero b%tesselated false...
  call ngroup_sort_atom(g,i)
endif
        
end subroutine ngroup_attach_atom

subroutine ngroup_detach_atom(g,a)
class(ngroup)          :: g
class(atom),target     :: a
integer                :: i

! Remove atom from list
if(g%listed) then
  i=a%gid(g)
  call ngroup_unsort_atom(g,i)
endif

! Detach atom
call g%ref%detach(a)
call g%b%detach(a)
call g%igroup_detach_atom(a)
 
! FIXME: QUE PASA SI el inice de B fue actualizado
! Poner que si b%no es tessellado, hay que actualizar lista.

! Clean null items
if(g%update) then
  deallocate(g%nn,g%list)
  allocate(g%nn(size(g%a)))
  allocate(g%list(size(g%a),g%mnb))
  g%listed=.false.
endif
               
end subroutine ngroup_detach_atom

subroutine ngroup_unsort_atom(g,i)
class(ngroup),intent(inout)  :: g
class(ngroup),pointer        :: ng
type(atom),pointer           :: a,aj
integer                      :: i,j,jj,ii,k
real(dp)                     :: rd,vd(dm)
type(atom_dclist),pointer    :: la
    
! Point to atom
a => g%a(i)%o
   
! Set cero 
! If atom ith is only in b, this is OK too
g%nn(i)=0

! Continue only if atom is in b.
if(a%gid(g%b)<1) return

! Remove it from list
la => g%ref%alist
do jj = 1,g%ref%nat
  la => la%next
  aj => la%o
  j = aj%gid(g)
                    
  do ii = 1,g%nn(j)
    if(g%list(j,ii)/=i) cycle

    g%nn(j)=g%nn(j)-1
    do k=ii,g%nn(j)
      g%list(j,k)=g%list(j,k+1)
    enddo
    exit

  enddo

enddo

 
end subroutine ngroup_unsort_atom

subroutine ngroup_sort_atom(g,i)
! Add atom `a` from ngroup `g` and fix neighbor list.
! Assume `a` is not already in `g` subgroup.
class(ngroup),intent(inout)  :: g
type(atom),pointer           :: a,aj
integer                      :: i,j,jj,m
real(dp)                     :: rd,vd(dm)
real(dp)                     :: rcut
type(atom_dclist),pointer    :: la

! Point to atom
g%nn(i)=0
a => g%a(i)%o

! If atom is in the ref group
if(a%gid(g%ref)==0) call g%lista_atom(i)

! If atom is in the b group
if(a%gid(g%b)>0) then
 
  ! Cut radious
  rcut=(g%rcut+nb_dcut)
  rcut=rcut*rcut
 
  la => g%ref%alist
  do jj = 1,g%ref%nat
    la => la%next
    aj => la%o
    j = aj%gid(g)
                      
    ! Skip autointeraction
    if(associated(aj,target=a)) cycle
                      
    vd = vdistance(a,aj,mic)
    rd = dot_product(vd,vd)

    if (rd>rcut) cycle

    ! Add i as neighbor of j.
    m = g%nn(j)+1
    g%list(j,m)=i
    g%nn(j)=m

  enddo
 
endif
 
end subroutine ngroup_sort_atom
     
! Set properties
! ==============

subroutine ngroup_setrc(g,rc)
! I could use the g%b%rcut directly.. but may be confuse...
class(ngroup)       :: g
real(dp),intent(in) :: rc

! New cut radious
g%rcut = rc
g%rcut2 = rc*rc
maxrcut=max(maxrcut,rc+nb_dcut)
g%b%rcut=rc+nb_dcut

end subroutine ngroup_setrc

subroutine ngroup_setlista(g)
class(ngroup) :: g

! Update cells
call g%b%tessellate()

if(.not.g%autoswitch) return

if(g%b%tessellated) then
  g%lista => ngroup_cells
  g%lista_atom => ngroup_cells_atom
else
  g%lista => ngroup_verlet
  g%lista_atom => ngroup_verlet_atom
endif

end subroutine ngroup_setlista

! Search algorithms
! =================

subroutine ngroup_verlet(g)
! Build neighbors verlet list.
class(ngroup),intent(inout)  :: g
type(atom),pointer           :: ai,aj
integer                      :: i,ii,j,m
real(dp)                     :: rd,vd(dm)
real(dp)                     :: rcut
type(atom_dclist),pointer    :: la

! Set ceros (TODO: I think this is not needed)
g%nn(:)=0

! Cut radious
rcut=(g%rcut+nb_dcut)
rcut=rcut*rcut

! TODO: Otra opcion para evaluar sería
! !$OMP PARALLEL DO SCHEDULE(STATIC,5) PRIVATE(m,vd,rd)
! do i = 1,g%na+g%nag-1
! m = 0
! do j = i+1,g%na+g%nag
!   vd = g%at(j)%o%pos-g%at(i)%o%pos

!OMP: Creo varios threads
!$OMP PARALLEL

!OMP: Un thread ejecuta, los demas quedan idle.
!$OMP SINGLE
la => g%ref%alist
do ii = 1,g%ref%nat
  la => la%next
  ai => la%o
  i = ai%gid(g)

  !OMP: Se reparte la tarea entre los idle threads
  !$OMP TASK DEFAULT(NONE)     &
  !$OMP& FIRSTPRIVATE(ai,i)    &
  !$OMP& SHARED(g,rcut,mic)    &
  !$OMP& PRIVATE(aj,j,vd,rd,m)

  ! Reset number of neighbors
  m=0

  do j = 1, g%b%amax
    aj => g%b%a(j)%o
    if(.not.associated(aj)) cycle

    ! Skip autointeraction
    if(associated(aj,target=ai)) cycle

    vd = vdistance(ai,aj,mic)
    rd = dot_product(vd,vd)

    if (rd>rcut) cycle

    ! Add j as neighbor of i.
    m=m+1
    g%list(i,m)=aj%gid(g)

  enddo
  g%nn(i)=m

  !$OMP END TASK
enddo

!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL

g%listed=.true.

end subroutine ngroup_verlet

subroutine ngroup_verlet_atom(g,i)
! Search neighbors for atom i.
class(ngroup),intent(inout)  :: g
type(atom),pointer           :: ai,aj
integer                      :: i,ii,j,m
real(dp)                     :: rd,vd(dm)
real(dp)                     :: rcut
type(atom_dclist),pointer    :: la

g%nn(i)=0
ai => g%a(i)%o

! Cut radious
rcut=(g%rcut+nb_dcut)
rcut=rcut*rcut

! Reset number of neighbors
m=0

do j = 1, g%b%amax
  aj=>g%b%a(j)%o
  if(.not.associated(aj)) cycle

  ! Skip autointeraction
  if(associated(aj,target=ai)) cycle

  vd = vdistance(ai,aj,mic)
  rd = dot_product(vd,vd)

  if (rd>rcut) cycle

  ! Add j as neighbor of i.
  m=m+1
  g%list(i,m)=aj%gid(g)

enddo
g%nn(i)=m

end subroutine ngroup_verlet_atom

subroutine ngroup_cells(g)
! Build neighbors verlet list over linked cells.
class(ngroup),intent(inout)  :: g
type(atom), pointer          :: ai, aj
integer                      :: i,ii,j,ic,k
integer                      :: nabor
real(dp)                     :: rd,vd(dm)
real(dp)                     :: rcut
integer                      :: rc(dm),nc(dm)
type(atom_dclist),pointer    :: la

! Set ceros
g%nn(:)=0
 
! Cut radious
rcut=(g%rcut+nb_dcut)
rcut=rcut*rcut

! Sort in cells
call g%b%sort()

!Bucle sobre los atomos de g
la => g%ref%alist
!$OMP PARALLEL
!$OMP SINGLE
do ii = 1,g%ref%nat
  la => la%next
  ai => la%o

  !$OMP TASK DEFAULT(NONE)             &
  !$OMP& FIRSTPRIVATE(ai)              &
  !$OMP& SHARED(g,rcut,mic)            &
  !$OMP& PRIVATE(rc,i,nabor,nc,j,aj,k,vd,rd)

  i = ai%gid(g)

  ! Get cell index
  rc(:)=int(ai%pos(:)/g%b%cell(:))+1

  !Bucle sobre las celdas vecinas a rc(:)
  do nabor=0,26
    nc(:)=map(:,nabor)+rc(:)

    ! Check if cells should be considered and apply PBC
    if(.not.cell_pbc(g%b,nc,mic)) cycle

    !Bucle sobre los atomos de la celda vecina
    j = g%b%head(nc(1),nc(2),nc(3))
    do while( j>0 )
      aj=>g%b%a(j)%o
      k = aj%gid(g)

      ! Next here to allow cycle
      j = g%b%next(j)

      ! Skip autointeraction
      if(associated(aj,target=ai)) cycle

      vd = vdistance(aj,ai,mic) ! respetar el orden
      rd =  dot_product(vd,vd)

      if (rd<rcut) then
        g%nn(i)=g%nn(i)+1
        g%list(i,g%nn(i))=k
      endif

    enddo

  enddo

  !$OMP END TASK

enddo

!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL

g%listed=.true.
 
end subroutine ngroup_cells

subroutine ngroup_cells_atom(g,i)
! Search neighbors for atom i. Asume b is already sorted in cells.
class(ngroup),intent(inout)  :: g
type(atom), pointer          :: ai, aj
integer                      :: i,ii,j,ic,k
integer                      :: nabor
real(dp)                     :: rd,vd(dm)
real(dp)                     :: rcut
integer                      :: rc(dm),nc(dm)
type(atom_dclist),pointer    :: la
  
! Set ceros
g%nn(i)=0
ai => g%a(i)%o

! Cut radious
rcut=(g%rcut+nb_dcut)
rcut=rcut*rcut

! Get cell index
rc(:)=int(ai%pos(:)/g%b%cell(:))+1

!Bucle sobre las celdas vecinas a rc(:)
do nabor=0,26
  nc(:)=map(:,nabor)+rc(:)

  ! Check if cells should be considered and apply PBC
  if(.not.cell_pbc(g%b,nc,mic)) cycle

  !Bucle sobre los atomos de la celda vecina
  j = g%b%head(nc(1),nc(2),nc(3))
  do while( j>0 )
    aj=>g%b%a(j)%o
    k = aj%gid(g)

    ! Next here to allow cycle
    j = g%b%next(j)

    ! Skip autointeraction
    if(associated(aj,target=ai)) cycle

    vd = vdistance(aj,ai,mic) ! respetar el orden
    rd =  dot_product(vd,vd)

    if (rd<rcut) then
      g%nn(i)=g%nn(i)+1
      g%list(i,g%nn(i))=k
    endif

  enddo

enddo

end subroutine ngroup_cells_atom

! Update search
! =============

subroutine update()
! Update all neighbor lists
use gems_groups, only: pbcghost, useghost
class(ngroup), pointer       :: g
integer                      :: i
type(atom_dclist),pointer    :: la

! call system_clock(t1)
nupd_vlist = nupd_vlist +1

! Add/delete ghosts
if(useghost) call pbcghost(maxrcut)
la => sys%alist
do i = 1,sys%nat
  la => la%next
  la%o%pos_old = la%o%pos
enddo

! Simple loop
do i = 1, ngindex%size
  g => ngindex%o(i)%o
  call ngroup_setlista(g)
  if(associated(g%lista)) call g%lista()
enddo

end subroutine update

subroutine inq_dispmax(g,dispmax1,dispmax2)
class(group)              :: g
type(atom_dclist),pointer :: la
integer                   :: i
real(dp),intent(out)      :: dispmax1, dispmax2
real(dp)                  :: rd,vd(dm)

!vecinos con el propio sistema
dispmax1 = 1.e-16_dp
dispmax2 = 1.e-16_dp

la => g%alist
do i = 1,g%nat
  la => la%next

  if(associated(la%o%prime)) cycle

  ! TODO: la%pos_old, should be an array inside ngroup
  ! so individual updates can be done
  vd(:) = la%o%pos(:) - la%o%pos_old(:)
  rd = dot_product(vd,vd)

  if (rd>dispmax1) then
    dispmax2 = dispmax1
    dispmax1 = rd
  else
    if (rd>dispmax2) dispmax2 = rd
  endif

enddo

end subroutine inq_dispmax

subroutine test_update()
! Check if neighbor update is needed
use gems_groups, only: pbcghost_move, useghost, do_pbc
real(dp)                   :: dispmax1,dispmax2
integer                    :: i
type(atom_dclist),pointer  :: la
class(ngroup),pointer      :: ng

if(useghost)then
  ! Update the ghost positions
  call pbcghost_move
else
  ! Needed to avoid atoms outside box when doing neighboor list (on interact)
  call do_pbc(sys)
endif
        
do i = 1, ngindex%size
  ng => ngindex%o(i)%o

  ! Capture update flag
  if (.not.ng%listed) then
    ! TODO: Individual calls to ng%lista()
    ! need old_pos to be an array of ngroup
    call update()
    exit
  endif

  ! Dispmax criteria for update
  call inq_dispmax(ng,dispmax1,dispmax2)

  if (sqrt(dispmax1)+sqrt(dispmax2)>nb_dcut) then
    call update()
    exit
  endif
enddo 

end subroutine test_update

! Variables and Labels

function polvar_ngroup(var) result(g)
use gems_variables, only: polvar, polvar_find
use gems_errors, only: werr
character(*),intent(in)      :: var
type(polvar),pointer         :: pv
class(ngroup),pointer    :: g

call werr('Labels should start with colon `:` symbol',var(1:1)/=':')
pv=>polvar_find(var)

g=>null()
if(.not.associated(pv)) return

call werr('Variable not linked',pv%hard)

! Print
select type(v=>pv%val)
class is (ngroup)
  g=>v
class default
  call werr('I dont know how to return that',.true.)
end select

end function

subroutine polvar_neighbor()
use gems_variables, only: polvar_link
call polvar_link('nb_dcut',nb_dcut)
end subroutine

end module
