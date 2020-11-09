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


module gems_groups
use gems_atoms
use gems_constants,only:sp,dp,dm,find_io
use gems_elements,only:elements,ncsym,element,inq_z

implicit none
private

type, public :: group

  ! An string helping to identify the group. Usefull in debug problems.
  character(4)                 :: nam='JDoe'

  ! Lista de atomos
  integer                     :: nat=0
  type(atom_dclist),pointer   :: alist=>null()

  ! hypervectores: vectores en el espacio de las fases
  !------------------------------------------------------------------------------
  ! Ciertas subrutinas (i.e. lbfgs) necesitan de un vector del espacio de las
  ! fases (velocidad, posicion, etc). Otras (i.e. neb, mad, claculo de
  ! correlaciones ) necesitan en una acumulacion en la memoria de estos vectores.
  ! De alli estas declaraciones. La forma de asociar estas variables con los
  ! objetos correspondientes, es recurrir a las subrutinas de asociacion para
  ! estos vectores al final de este modulo :

  real(dp),pointer            :: pp(:)=>null(),pv(:)=>null(),pf(:)=>null(),pa(:)=>null() ! vectores actuales

  ! -----Grupos allocateados
  logical                     :: empty=.true. ! Boolean for faster empty test

  !================Propiedades Mecanicas

  ! Estas propiedades deben calcularse usando el modulo de inquire_properties.
  ! Sus definiciones estan estrechamente relacionadas con este modulo


  ! -----temperatura
  real(dp)                    :: temp    =0.0_dp  ,&
                                 tempvib =0.0_dp  ,&
                                 temprot =0.0_dp  ,&
                                 mass    =0.0_dp
  ! -----energies
  real(dp)                    :: ekin    =0.0_dp  ,&
                                 epot    =0.0_dp
  ! -----center of mass
  real(dp),dimension(dm)      :: cm_pos  =0.0_dp  ,&
                                 cg_pos  =0.0_dp  ,&
                                 cm_vel  =0.0_dp
  ! -----rg
  real(dp)                    :: rg_pos  =0.0_dp
  ! -----angular properties
  real(dp)                    :: erot    =0.0_dp  ,&
                                 evib    =0.0_dp
  real(dp),dimension(3)       :: ang_mom =0.0_dp  ,&
                                 ang_vel =0.0_dp
  real(dp),dimension(3,3)     :: inercia =0.0_dp

  real(dp),dimension(3,3)     :: covar= 0.0_dp               ! Matriz de covarianza

  ! Reservado para las contribuciones a la energia potencial
  real(dp)                      :: epot1=0.0_dp, &
                                   epot2=0.0_dp, &
                                   epot3=0.0_dp, &
                                   epot4=0.0_dp, &
                                   epot5=0.0_dp


  ! ----- pressure
  real(dp),dimension(3,3)     :: virial=0.0_dp,  &
                                 pressure=0.0_dp


  ! -----geometria y morfologia
  real(dp),dimension(dm)      :: maxpos  =0.0_dp  ,&
                                 minpos  =0.0_dp
  real(sp),dimension(3)       :: mainaxis =0.0_sp

  ! Para cada propiedad le asigno un booleano que me indica si esta necesita
  ! ser calculada. Por ejemplo, me ayudaria para saber si tengo que recalcular
  ! la masa cuando estoy calculando el centro de masa o no. La idea es que el
  ! inquire sea mas eficientes

  ! -----temperatura
  logical                     :: b_temp    =.false. ,&
                                 b_tempvib =.false. ,&
                                 b_temprot =.false. ,&
                                 b_mass    =.false.
  ! -----energies
  logical                     :: b_ekin    =.false. ,&
                                 b_epot    =.false.
  ! ----- pressure
  logical                     :: b_virial   =.false. ,&
                                 b_pressure =.false.

  ! -----center of mass
  logical                     :: b_cm_pos  =.false. ,&
                                 b_cg_pos  =.false. ,&
                                 b_cm_vel  =.false.
  ! -----radio de giro
  logical                     :: b_rg_pos  =.false.

  ! -----angular properties
  logical                     :: b_erot    =.false. ,&
                                 b_evib    =.false.
  logical                     :: b_ang_mom =.false. ,&
                                 b_ang_vel =.false.
  logical                     :: b_inercia =.false. ,&
                                 b_covar   =.false.
  ! -----geometria y morfologia
  logical                     :: b_maxpos  =.false. ,&
                                 b_minpos  =.false. ,&
                                 b_mainaxis=.false.

  ! Esta indice permite saber si el grupo es usado en variables colectivas para
  ! evitar que sea borrado dejando colgado una variable colectiva.
  integer                     :: cvs=0

  contains

  procedure :: group_initialize
  generic :: init => group_initialize
  procedure :: dest => group_destroy
 !   !procedure :: del => group_atom_del
 !   procedure :: delall => group_allatom_del
  procedure :: add => group_include
  procedure :: atom => group_atombyindex ! devuelve un puntero correspondiendo a un indice desde el head

end type group

! Linked List of groups
#define _NODE group_l
#define _CLASS class(group)
#include "list_header.inc"
public :: group_l
public :: glist_add
 
! This linked list holds all allocated group.
type(group_l),public      :: glist
 
! Array of Pointers to groups
#define _NODE group_ap
#define _CLASS class(group)
#include "arrayofptrs_header.inc"
public :: group_ap
 

! Module procedures 

! public system_group_add
! public system_increase
public :: group_destroy
public :: group_atom_add
public :: dgroup_atom_add
public :: group_atom_del
public :: group_allatom_del
public :: group_switch_vectorial
public :: group_switch_objeto
         

contains

#define _NODE group_l
#define _CLASS class(group)
#include "list_body.inc"


! Constructor

! subroutine alocal2aghost(la)
! type ( atom_dclist ), pointer :: la
! integer,intent(in)            ::i
!
! ! Deattach from alocal list
! la%prev%next => la%next
! la%next%prev => la%prev
! nlocal=nlocal-1
!
! ! Attach to aghost list
! la%prev => aghost%prev
! la%next => aghost
! la%prev%next => la
! la%next%prev => la
! nghost=nghost+1
!
! end subroutine alocal2aghost

! ---------------------------- Linked List Events

! Search

!function glist_check( glist, g )
!! Esta es una lista doble circular
!  type ( group_dclist ),pointer :: glist
!  type ( group_dclist ),pointer :: l1,l2
!  type ( group ),pointer        :: g
!  logical                       :: glist_check
!
!  glist_check = .true.
!
!  l1 => glist ! Para fijar el principio
!  l2 => l1%next
!  do while (.not.associated(l1,target=l2))
!    if (associated(l2%g,target=g)) return
!    l2 => l2%next
!  enddo
!
!  glist_check = .false.
!
!end function glist_check


! ---------------------------- Group Events

! Constructor

subroutine group_initialize( g,nam)
class(group),target        :: g
character(4),optional      :: nam

! Inicializo la lista linkeada de atomos
allocate(g%alist)
call g%alist%init()
g%nat = 0

! Registro el grupo en la lista de grupos
call glist%add_after()
call glist%next%point(g)

! Initializo la cabeza del grupo para acciones agrupadas
allocate(g%alist%o)
call g%alist%o%init()

if(present(nam)) g%nam=nam

end subroutine group_initialize

! Destructor

subroutine group_destroy ( g )
class(group)         :: g

! Not needed this:
! call g%alist%destroy_all()
! If we do this:
call group_allatom_del( g )

! Borro al grupo de la lista de grupos
call group_group_del(glist,g)

! Destroy node
call g%alist%o%dest()
deallocate(g%alist%o)

! Destroy head
deallocate(g%alist)
g%nat = 0

end subroutine group_destroy

! Increase

subroutine group_atom_add ( a, g )
! Agrega a g un atomo dado por el puntero a.
! IMPORTANTE: Este debe ser un puntero
type(group),target    :: g
type(atom),target     :: a

! Check if this atom is already in group g
! FIXME: I can do this quickly if the atom
! has an integer list of group uids.
if(associated(atom_dcl_find(g%alist,a))) return

! agrego el atomo a la lista linkeada del grupo
call g%alist%add_before()
call g%alist%prev%point(a)
g%nat = g%nat + 1 ! numero de particulas

! agrego el grupo a la lista linkeada del atomo.
! call a%glist%add_soft(g) 
! FIXME: I can do this quickly if the atom
! has an integer list of group uids.

! propiedades basicas para modificar
g%mass = g%mass + a % mass ! masa

! si esta vectorial, ahora no tiene sentido
if(associated(g%pp)) then
  deallocate(g%pp)
  deallocate(g%pv)
  deallocate(g%pa)
  deallocate(g%pf)
endif

if(g%nat>0) g%empty = .false.

end subroutine group_atom_add

subroutine dgroup_atom_add( a,g )
! Agrego un hard atomo, util para los grupos draft del select create
type(group),target          :: g
type(atom),target           :: a
type(atom_dclist),pointer   :: ln 

! Create and asign the atom
call g%alist%add_before(ln)
call ln%alloc()
call ln%o%init()
call atom_asign(ln%o,a)
 
!FIXME: Le doy un id al atomo
ln%o%id = g%nat + 1
                           
g%nat = g%nat + 1 ! numero de particulas
          
! agrego el grupo a la lista linkeada del atomo.
! call a%glist%add_soft(g)
! FIXME: I can do this quickly if the atom
! has an integer list of group uids.
 
! propiedades basicas para modificar
g%mass = g%mass + a % mass ! masa

! si esta vectorial, ahora no tiene sentido
if(associated(g%pp)) then
  deallocate(g%pp)
  deallocate(g%pv)
  deallocate(g%pa)
  deallocate(g%pf)
endif

if(g%nat>0) g%empty = .false.

end subroutine dgroup_atom_add  

subroutine group_include(g2,g1)
! agrega los atomos del g1 al g2
class(group)              :: g2
type(group)               :: g1
type(atom_dclist),pointer :: la
integer                   :: i

la => g1%alist
do i = 1,g1%nat
  la => la%next
  call group_atom_add(la%o,g2)
enddo

! si esta vectorial, ahora no tiene sentido
if(associated(g2%pp)) then
  deallocate(g2%pp)
  deallocate(g2%pv)
  deallocate(g2%pa)
  deallocate(g2%pf)
endif

if(g1%nat>0) g1%empty = .false.

end subroutine group_include

subroutine glist_add(glist,g)
! Add g to glist
type(group_l),target  :: glist
type(group),target       :: g

if(group_belong(g,glist)) return
call glist%add_after()
call glist%next%point(g)
end subroutine glist_add

! Decrease

subroutine group_atom_del ( a, g )
! destruye el eslabon solo si esta asociado a un atomo identificado
type  ( group ),target         :: g
type  ( group ),pointer        :: pg
type ( atom_dclist ), pointer :: la
type ( atom ),target          :: a
integer                       :: i

pg => g

la => g%alist
do i =1,g%nat
  la => la%next

  if (la%o%id==a%id) then

    ! saco el grupo de la lista del atomo
    ! call group_group_del(a%glist,g)
    ! FIXME: I can do this quickly if the atom
    ! has an integer list of group uids.

    ! saco el atomo de la lista del grupo
    call la%deattach()
    deallocate(la)
    g%nat = g%nat - 1

    ! propiedades extras para modificar
    g%mass = g%mass - a % mass

    if(g%nat==0) g%empty = .true.

    return

  endif

enddo

end subroutine group_atom_del

subroutine group_allatom_del(g)
! destruye los eslabones... mas no la inicialicion
type  ( group ),target         :: g
type  ( group ),pointer        :: pg
type ( atom_dclist ), pointer :: la,next

pg => g

! Circulo por la lista hasta que la vacio
la => g%alist%next
do while(g%nat/=0)

  ! saco el grupo de la lista del atomo
  ! call group_group_del(la%o%glist,g)
  ! FIXME: I can do this quickly if the atom
  ! has an integer list of group uids.

  ! saco el atomo de la lista del grupo
  next => la%next
  call la%deattach()
  deallocate(la)

  g%nat = g%nat - 1

  la=>next
enddo

! propiedades extras para modificar
g%mass = 0.0_dp
g%empty = .true.

end subroutine group_allatom_del

subroutine group_group_del(glist,g)
! destruye el eslabon solo si esta asociado al grupo g
type(group),target        :: g
type(group_l),target   :: glist
type(group_l),pointer  :: prev,lg

lg => glist

do while(associated(lg%next))
  prev=>lg
  lg=>lg%next
  if (associated(lg%o,target=g)) then
    call lg%deattach(prev)
    !call lg%o%dest() I will not destroy the group.
    deallocate(lg)
    return
  endif
enddo

end subroutine group_group_del

function group_belong(g,glist) result(check)
! Check if g is in glist.
type(group),target       :: g
type(group_l),target  :: glist
type(group_l),pointer :: lg
logical                  :: check

check = .true.

lg => glist%next
do while(associated(lg))
  if (associated(lg%o,target=g)) return
  lg => lg%next
enddo

check = .false.

end function group_belong

! Seleccion

function group_atombyindex(this,i) result(at)
class(group)            :: this
type(atom),pointer      :: at
integer,intent(in)      :: i
integer                 :: j
type ( atom_dclist ), pointer :: la

if (i>this%nat .or. i<0) at => null()

la => this%alist
do j=1,i
  la => la%next
enddo
at => la%o

endfunction

! ---------------------------- Hyper vector Constructors

  ! DEL PROBLEMA FUNDAMENTAL DEL GMD

  ! El echo de que el grupo sea una lista linkeada esta indicando que es una
  ! selecciona caprichosa de atomos. Esto permite aplicar disitntas subrutinas a
  ! distintas selecciones caprichosas.

  ! En muchas subrutinas, y para muchas cosas, es util constar con vectores que
  ! representen el estado de este grupo.

  ! Fotran no puede associar a un unico puntero a una selccion
  ! caprichosa de los componentes de un vector (i.e. v(1) v(4) y v(9)), y menos
  ! a la selccion caprichosa de los componendes de distintos vectores (i.e. v(1)
  ! v(4) y r(9)).

  ! De esta manera, los punteros de posicion velocidad fuerza y aceleracion de
  ! un grupo (estos son g%pp, g%pv, gp%f y g%pa ) no pueden exisitr, ya que, una
  ! seleccion de atomos caprichosa por parte de un grupo es a priori imposible
  ! de odenar en forma vectorial.

  ! Las subrutinas que siguen a continuacion son un parche para lograr esto. Las
  ! variables atomicas, se encuentran inicialmente apuntando a vectores en el
  ! systema. Luego de la seleccion caprichos por un grupo, se puede copiar
  ! las variables atomicas a un vector y asociar los atomos a ese vector,
  ! logrando asi un ordenamiento del grupo de forma vectorial. Si este
  ! ordenamiento no existe, los punteros del grupo tendran la condicion null. Si
  ! por el contrario este ordenamiento existe, los punteros del grupo apuntaran
  ! al vector determinado. Esto involucra ciertos manejos que se deben realizar
  ! con las siguientes subrutinas:

subroutine group_switch_vectorial(g,switched)
class(group),intent(inout)     :: g
logical,optional,intent(out)   :: switched
integer                        :: i
type (atom_dclist),pointer     :: la

if(present(switched)) switched=.false.
if(associated(g%pp)) return
if(present(switched)) switched=.true.

allocate(g%pp(g%nat*dm))
allocate(g%pv(g%nat*dm))
allocate(g%pa(g%nat*dm))
allocate(g%pf(g%nat*dm))

la => g%alist
do i = 1,g%nat
  la => la%next
  g%pp((i-1)*dm+1:i*dm) =  la%o%pos
  g%pv((i-1)*dm+1:i*dm) =  la%o%vel
  g%pf((i-1)*dm+1:i*dm) =  la%o%force
  g%pa((i-1)*dm+1:i*dm) =  la%o%acel
  deallocate(la%o%pos  )
  deallocate(la%o%vel  )
  deallocate(la%o%force)
  deallocate(la%o%acel )
  la%o%pos    => g%pp((i-1)*dm+1:i*dm)
  la%o%vel    => g%pv((i-1)*dm+1:i*dm)
  la%o%force  => g%pf((i-1)*dm+1:i*dm)
  la%o%acel   => g%pa((i-1)*dm+1:i*dm)
enddo

end subroutine

subroutine group_switch_objeto(g,switched)
class(group),intent(inout)     :: g
logical,optional,intent(out)   :: switched
integer                        :: i
type (atom_dclist),pointer     :: la

if(present(switched)) switched=.false.
if(.not.associated(g%pp)) return
if(present(switched)) switched=.true.


la => g%alist
do i = 1,g%nat
  la => la%next
  allocate(la%o%pos  (dm))
  allocate(la%o%vel  (dm))
  allocate(la%o%force(dm))
  allocate(la%o%acel (dm))
  la%o%pos    = g%pp((i-1)*dm+1:i*dm)
  la%o%vel    = g%pv((i-1)*dm+1:i*dm)
  la%o%force  = g%pf((i-1)*dm+1:i*dm)
  la%o%acel   = g%pa((i-1)*dm+1:i*dm)
enddo

deallocate(g%pp); g%pp=>null()
deallocate(g%pv); g%pv=>null()
deallocate(g%pa); g%pa=>null()
deallocate(g%pf); g%pf=>null()

end subroutine

end module gems_groups

