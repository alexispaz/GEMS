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


module gems_program_types
use gems_constants,only:sp,dp,dm,find_io
use gems_elements,only:elements,ncsym,element,inq_z

!------------------------------------------------------------------------------
! this module have the gcmd variables-objects declaration :
! - atoms
! - system
! - elements
! - algorithms
! - potentials
!------------------------------------------------------------------------------

implicit none

! all the variables here declared global in the module keep value
save

! all the variables here declared are the public domain, pero no me parece que
! sea bueno que las variables se vayan heredando.
! Es interesante esta nota sobre el echo:
!One of the problems with Fortran is that when you import from
!modules you will always throw everything in the global
!namespace, as in Javascript... Every time you use USE, you are including
!everything into the global namespace.  This is one of the worst Fortran 9X
!issues. -- Stefano Borini Aug 6 '09 at 21:11
! Por ejemplo, el public aca produce que en cualquier modulo que use
! program_types, este usando las variables sp,dp,etc.. no declaradas en este
! modulo propiamente dicho
! Esto puede llevar a que los ultimos moudlos tengan un solo use pero usen
! millones de variables globales y no se lo pueda desacoplar del programa.
! la solucion: para dejar por default public un modulo, pongo public, pero corto la
! herencia con private a los use del modulo.
private

!                                                               particular vars
!------------------------------------------------------------------------------

! The program and subprograms trace this variable to get the term order
! After the term_signal is read, is set to .false. too allow subprograms end
! without quit gems.
logical,public :: term_signal=.false., cycle_signal=.false.


!                                                             program variables
!------------------------------------------------------------------------------
real(dp),target,public    :: dm_steps=0._dp
real(dp),public    :: nframe=0._dp   ! real para que se banque numeros grandes
real(dp),public    :: ptime,pnframe

!                                                           dimension variables
!------------------------------------------------------------------------------

  ! Elijo el numero de dimensiones a tiempo de compilacion con el precompilador

! Especifico numeros que establecen la dimension absolutas del programa.
real(dp),public,parameter :: dpos=0.2 ! Diferencial de posicion que es considerado cero (usado
                                      ! para minimizar y otras cosas)


!                                                           dinamical variables
!------------------------------------------------------------------------------

! Tiempo simulado
real(dp),public,target    :: time=0.0_dp

logical,public            :: binout=.true.

integer,public            :: frame=0 ! Cuando escribe: el frame actual del outunit
integer,public            :: nstep=100 ! number of integration step
real(dp),public           :: dtaux=0.002 ! integration step

real(dp),target,public    :: dt=0.002

integer,parameter,public,dimension(26,3) :: n1cells = transpose(reshape( &
                   [ 1, 0, 0, -1, 0, 0,  0, 1, 0,&
                     1, 1, 0, -1, 1, 0,  0,-1, 0,&
                     1,-1, 0, -1,-1, 0,  0, 0, 1,&
                     1, 0, 1, -1, 0, 1,  0, 1, 1,&
                     1, 1, 1, -1, 1, 1,  0,-1, 1,&
                     1,-1, 1, -1,-1, 1,  0, 0,-1,&
                     1, 0,-1, -1, 0,-1,  0, 1,-1,&
                     1, 1,-1, -1, 1,-1,  0,-1,-1,&
                     1,-1,-1, -1,-1,-1 ],[3,26]))
! integer,public,dimension(26,3) :: n1test = 0.5_dp*(sign(n1cells(:,:))+1)-abs(n1cells(:,:))

public :: box_setvars,box_expand,boxed
real(dp),public,dimension(dm,dm) :: tbox   =0.0_dp
real(dp),public,dimension(dm),target    :: box    =1.0e6_dp  
real(dp),public,dimension(dm)    :: box2   =5.0e5_dp  , &
                                    one_box=1.0e-6_dp , &
                                    one_box2=2.5e-5_dp, &
                                    box_old=0.0_dp      ! For pbcghost. Set this small to force initial ghost inclusion
real(dp),public                  :: box_vol=1.0e18_dp
logical                          :: boxed=.false.

real(dp),public             :: virial(dm,dm)=0._dp
logical,public              :: b_gvirial=.false., b_avirial=.false.
real(dp),public,allocatable :: atomvirial(:,:)
public                      :: dovirial

!                                                             dimensions arrays
!------------------------------------------------------------------------------
integer,public,parameter   :: nb_max=150 ! Vecinos maximos
integer,public             :: ncell


!                                                                    Estructuras
!------------------------------------------------------------------------------


!!!!!! ATOMOS

! TODO: Atoms deberia tener una gerarcia parecida a:
! atomo mecanico (solo r, epot y f)
! atomo dinamico (mecanico con v, a, m, etc)
! de esta manera los intergroups tendrian la posibilidad
! de tener un atomo mecanico

type :: atom

  ! Para llevar un conteo de cuantos links a un atomo hay
  ! Algo interesante aca http://stackoverflow.com/a/38363298/1342186

  ! Esta lista indica los grupos en los cuales el atomo se encuentra inmerso
  type(group_l),pointer :: glist=>null()

  ! Este es el eslabon de la lista del system donde se encuentra allocateado
  type(atom_dclist),pointer :: link=>null()

  ! ----- Propiedades relacionadas al elemento
  ! Propiedades que defino afuera de e para que se mas rapidamente accedida
  ! (en general la 1/masa esta en los cuellos de botella de los algoritmos)
  integer           :: z=119 ! The generic element
  real(dp)          :: mass=1.0_dp,one_mass=1.0_dp,one_sqrt_mass=1.0_dp
  real(dp)          :: q=0.0_dp  ! Carga
  real(dp)          :: s=1.0_dp  ! sigma
  real(dp)          :: e=0.0_dp  ! epsilon
  character(ncsym)  :: sym
  integer           :: sp=0      ! Hybridization

  ! Constrain. Si bconst=true el atomo tiene un constrain. Se colapsa la
  ! fuerza en direccion al vector vconst si lconst=T o se borra la componente
  ! de la fueza en direccion al vector si lconst=F. Idem con la velocidad. Asi
  ! la particula queda fija en un plano o en un eje. Tambien la puedo forzar
  ! directamente haciendolo con la posicion
  real(dp)              :: vconst(dm)=0.0_dp
  real(dp),allocatable  :: pconst(:) ! posicion incial del constrain
  logical               :: bconst=.false.,lconst=.false.

  !  Enlaces y moleculas.... TOFIX
  integer      :: abondid(20)=0  ! el indicie dentro de la molecula de los asociados
  integer      :: abonds=0  ! el numero de asociados
  integer      :: molid=0   ! el indice de la molecula
  integer      :: amolid=0  ! el indice dentro de la molecula

  ! ----- Propiedades mecanicas
  real(dp),pointer       :: pos(:)=>null(),   &!propieties of atom. [a][..][m/s][..]
                            force(:)=>null(), &
                            acel(:)=>null(),  & !aceleracion
                            vel(:)=>null()

  !In a local atom it has the info to unwrap coordinates. In the ghost atom,
  !it has the info of the subdomain/processor it belongs.
  integer                :: boxcr(dm)=0
  logical                :: pbc(dm)=.false. !PBC para ese atomo

  real(dp),dimension(dm) :: acel2  =0._dp,& !derivada primera de la aceleración
                            acel3  =0._dp,& !derivada segunda de la aceleración
                            acel4  =0._dp,& !derivada tercera de la aceleración
                            pos_eq =0._dp,& !para ver el desplazamiento y decidir entrar al hyperespacio
                            pos_v  =0._dp,& !posicion relativa al punto v
                            vel_v  =0._dp,& !velocidad relativa al punto v
                            vel_rot=0._dp,& !velocidad de rotacion
                            vel_vib=0._dp,& !velocidad de vibracion
                            pos_cm =0._dp,& !posicion relativa al cm del grupo
                            vel_cm =0._dp   !velocidad relativa al cm del grupo

  !para ver el desplazamiento en la lista de vecinos. Esto lo establezco bien
  !grande para forzar la primera actualizacion del verlet
  real(dp),dimension(dm) :: pos_old =1.e8_dp

  real(dp)               :: epot=0.d0,                & !energia potencial total[ev]
                            erot=0.d0,erot_ss=0.d0,   & !energia rotacional relative to system and ss [ev]
                            evib=0.d0,evib_ss=0.d0,   & !energia vibracional relative to system and ss [ev]
                            rho=0.d0,cord=0.d0,       & !densidad.. o algun otro parametro
                            border=0.d0                 !Orden de Enlace

  ! ----- Propiedades relacionadas al grupo
  ! El id hace referencia a la lista de vecinos, asique ojo.
  integer                :: id=0,idv=0,tag=0 !referrncia al systema

  ! Si un atomo aparece, se lee estas variables para decidir en que
  ! intergroup se debe ubicar
  integer                :: nigr=0   !referrncia al intergrupo
  integer,allocatable    :: igr(:)   !referrncia al intergrupo
  logical,allocatable    :: igra(:)  !T if belogn to a group or F if belong to b
  logical                :: delete=.false.

  real(dp)               :: maxdisp2=0 !desplazamiento maximo a un determinada T de grupo

  contains

  procedure :: init => atom_allocate
  procedure :: dest => atom_destroy
  procedure :: bonded => atom_inquire_bonded ! pregunta si estan enlazados
  procedure :: setbond => atom_set_bonded ! guarda el enlace entre dos atomos
  procedure :: setz => atom_setelmnt_byz
  procedure :: setsym => atom_setelmnt_bysym

  procedure :: distance_point => atom_distancetopoint
  procedure :: distance_axis => atom_distancetoaxis
  generic   :: distance => distance_axis,distance_point

  !procedure :: del => atom_atom_del
  !procedure :: delall => atom_allatom_del
  !procedure :: addatom => atom_include
  !procedure :: atom_asign
  !generic   :: assignment(=) => atom_asign

  procedure :: addigr => atom_addigr
  procedure :: deligr => atom_deligr

end type atom

#define _NODE atom_dclist
#define _CLASS class(atom)
#include "cdlist_header.inc"

#define _NODE atom_ap
#define _CLASS class(atom)
#include "arrayofptrs_header.inc"

public :: atom_ap

! Lista dura, con locales, ghost, y remotos
integer,public                            :: natoms=0
type(atom_dclist),target,public           :: atoms
type(atom_ap),target,allocatable,public   :: a(:)

! Lista blanda con atomos locales
integer,public                      :: nlocal=0
type(atom_dclist),target,public     :: alocal

! Lista blanda con atomos fantasma
logical,public                      :: mic=.true.
logical,public                      :: ghost=.false.
logical,public                      :: fullghost=.false.
integer,public                      :: nghost=0
type(atom_dclist),target,public     :: aghost

interface atom_setelmnt
  module procedure atom_setelmnt_bysym,atom_setelmnt_byz
end interface
public atom_setelmnt,atom_asign


public :: atom,atom_dclist
public :: atoms_group_add, new_ghostatom, del_atom

public :: ghost_group_add, atom_dclist_destroyall




!!!!!! GRUPOS

type :: group

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


! List and vectos of groups
#define _NODE group_l
#define _CLASS class(group)
#include "list_header.inc"

! La lista de grupos registra todos los grupos creados
type(group_l),public      :: glist
public :: group_l,group,glist_add

#define _NODE group_ap
#define _CLASS class(group)
#include "arrayofptrs_header.inc"

public :: group_ap



logical,public,target,allocatable      :: fix(:),join(:)
integer,public                         :: njoin=0

! Declaraciones de Tipos

type(group),target,public       :: sys     ! Systema total
integer,parameter,public        :: mgr=9
type(group),target,public       :: gr(mgr) ! Selecciones
type(group),target,public       :: gsel    ! grupo seleccionado
type(group),target,public       :: gnew    ! grupo creacion


logical ,public    ::   gbb=.false.,gdb=.false.


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

#define _NODE atom_dclist
#define _CLASS class(atom)
#include "cdlist_body.inc"

#define _NODE group_l
#define _CLASS class(group)
#include "list_body.inc"

! Set box size

subroutine box_setvars()
 use gems_input_parsing
 integer   :: i

 boxed=.true.

 ! For non-cubic boxes this should be box(:)=matmul(tbox(i,i),[1,1,1])
 box(1)=tbox(1,1)
 box(2)=tbox(2,2)
 box(3)=tbox(3,3)

 one_box(:) = 1.0_dp/box(:)

 ! For non-cubic boxes this should be the determinant of tbox
 box_vol = box(1)*box(2)*box(3)

end subroutine box_setvars

subroutine box_expand(a,b,c)
 real(dp)       :: a,b,c

 tbox(1,1)=tbox(1,1)*a
 tbox(2,2)=tbox(2,2)*b
 tbox(3,3)=tbox(3,3)*c

 call box_setvars()

end subroutine box_expand


! Constructor

subroutine atom_allocate(a)
use gems_errors
! inicializo los punteros y allocateables que no
! se pueden inicializar en la declaración
class(atom),intent(inout)    :: a

! crear un atomo
allocate(a%pos(dm))
allocate(a%vel(dm))
allocate(a%force(dm))
allocate(a%acel(dm))
a%acel(:)=0._dp
a%pos(:)=0._dp
a%vel(:)=0._dp
a%force(:)=0._dp
call atom_setelmnt(a,119)

! Inicializo la lista linkeada de grupos
allocate(a%glist)

end subroutine atom_allocate

! Destructor
subroutine atom_destroy(a)
class(atom)         :: a

deallocate(a%pos)
deallocate(a%vel)
deallocate(a%force)
deallocate(a%acel)
if(allocated(a%igr) ) deallocate(a%igr)
if(allocated(a%igra)) deallocate(a%igra)
call a%glist%destroy_all()
deallocate(a%glist)

end subroutine atom_destroy

! Asign

subroutine atom_asign(a1,a2)
! Copia la informacion de a2 en a1. Esto se hace sin importar como estan !
! conformados los objetos, pudiendo por ejemplo a2%pos ser un puntero slice
! o un arreglo duro, no importa.
type(atom) :: a1,a2

a1 % pos(:)   = a2 % pos(:)
a1 % vel(:)   = a2 % vel(:)
a1 % force(:) = a2 % force(:)
a1 % acel(:)  = a2 % acel(:)
call atom_setelmnt(a1,a2%z)

a1 % acel(:)    = a2 % acel(:)
a1 % acel2(:)   = a2 % acel2(:)
a1 % acel3(:)   = a2 % acel3(:)
a1 % acel4(:)   = a2 % acel4(:)
a1 % pos_old(:) = a2 % pos_old(:)
a1 % pos_v(:)   = a2 % pos_v(:)
a1 % vel_v(:)   = a2 % vel_v(:)

a1 % epot    = a2 % epot
a1 % erot    = a2 % erot
a1 % erot_ss = a2 % erot_ss
a1 % evib    = a2 % evib
a1 % evib_ss = a2 % evib_ss
a1 % rho     = a2 % rho
a1 % molid   = a2 % molid
a1 % sp      = a2 % sp
a1 % s       = a2 % s

end subroutine atom_asign

! Properties

  subroutine atom_addigr(a,i,igra)
   ! Set the atom variables related to the interaction group
   class(atom)         :: a
   integer,intent(in)  :: i
   logical,intent(in)  :: igra
   integer             :: n
   integer,parameter   :: spd=1
   integer,allocatable :: aux(:)
   logical,allocatable :: aux2(:)

   n=a%nigr

   if(.not.allocated(a%igr)) then
     allocate(a%igr(spd))
     allocate(a%igra(spd))
   endif

   if(size(a%igr)<n+1) then
     allocate( aux(n+spd) )
     aux(1:n)=a%igr(1:n)
     call move_alloc(to=a%igr,from=aux)

     allocate( aux2(n+spd) )
     aux2(1:n)=a%igra(1:n)
     call move_alloc(to=a%igra,from=aux2)
   endif

   a%igr(n+1)=i
   a%igra(n+1)=igra
   a%nigr=n+1

  end subroutine atom_addigr

subroutine atom_deligr(a,i)
! Delete an igr from atom a. This do not deallocate to speed future allocations.
class(atom)         :: a
integer,intent(in)  :: i
integer             :: j

if(.not.allocated(a%igr)) return

! Delete from array of pointers
a%nigr = a%nigr - 1
do j=i,a%nigr
  a%igr(j)=a%igr(j+1)
enddo

end subroutine atom_deligr

subroutine atom_setelmnt_byz(a,z)
  ! Establece las propiedades relacionadas al elemento en un atomo
  class(atom)               :: a
  integer                   :: z

  ! FIXME: What if is not there???

  a%z = z
  a%mass = elements%o(z)%mass
  a%sym = elements%o(z)%sym
  if(a%mass==0._dp) then
    a%one_mass = 0._dp
    a%one_sqrt_mass = 0._dp
  else
    a%one_mass = 1.0_dp/a%mass
    a%one_sqrt_mass = sqrt(a%one_mass)
  endif

end subroutine

   subroutine atom_setelmnt_bysym(a,sym)
     use gems_elements,only:add_z, inq_z
     class(atom)    :: a
     character(*)   :: sym

     ! Adding z just in case is not there
     call add_z(sym)

     call atom_setelmnt_byz(a,inq_z(sym))
   end subroutine

  function atom_inquire_bonded(this,a) result(ans)
   class(atom) :: this
   logical    :: ans
   type(atom) :: a
   integer    :: i

   ans=.false.
   if(this%molid/=a%molid) return

   ans=.true.
   do i=1,this%abonds
     if (this%abondid(i)==a%amolid) return
   enddo

   ans=.false.

  end function atom_inquire_bonded

  subroutine atom_set_bonded(this,a)
   class(atom) :: this
   type(atom) :: a
   integer    :: i

   ! Para no enlazarlo de nuevo
   do i =1,this%abonds
     if(this%abondid(i)==a%amolid) return
   enddo

   this%abonds=this%abonds+1
   this%abondid(this%abonds)=a%amolid

   a%abonds=a%abonds+1
   a%abondid(a%abonds)=this%amolid

  end subroutine atom_set_bonded

  function atom_distancetopoint(a,r) result(vd)
   !calculates the distance of two atoms also in pbc case
   class(atom),intent(in)   :: a
   real(dp),dimension(dm)  :: vd
   real(dp),intent(in)     :: r(dm)
   integer                 :: l

    ! Mas rapido usar idnint que un if
    vd=r-a%pos
    do l = 1,dm
      if (a%pbc(l)) vd(l)=vd(l)-box(l)*idnint(vd(l)*one_box(l))
    enddo

  end function atom_distancetopoint

  function atom_distancetoaxis(a,p,r) result(vd)
   ! a es el atomo, p y r dan la recta(t)=p+t*r
   !calculates the distance of two atoms also in pbc case
   class(atom),intent(in)   :: a
   real(dp),dimension(dm)  :: vd,aux
   real(dp),intent(in)    :: r(dm),p(dm)

    aux=a%distance(p)
    vd=aux-(dot_product(aux,r))*r

  end function atom_distancetoaxis

! ---------------------------- ATOMS AND GHOSTS

subroutine new_atom()
! agrega un nuevo atomo (vacio) al sistema. Luego este puede ser referenciado
! via a(natoms)%o
type ( atom ), pointer         :: o
integer                        :: j

! The atoms are allocated in the `atoms` cdlist
call atoms%add_before_hard()
o=>atoms%prev%o
call o%init()
natoms = natoms + 1

! Link the atom with its node in the `atoms` cdlist
o%link => atoms%prev
o%id  = natoms

! Indice para acomodar en vectores
j = (natoms-1)*dm+1

! The atom is registred in the array of pointers `a`
call atoms_increaseby(1)
a(natoms)%o => o

end subroutine new_atom

subroutine del_atom(o)
! Attempt to delete an atom. If an intergroup is still pointing to it, do not
! delete it but flag it as to be deleted.
! AN alternative would be just to delete it and check in the intergroups if the
! node is associated, but, I would like to have a flag system so in the future I
! can parallelize the system by regions and instead of delete it, move it to the
! other processors.
type(atom), target          :: o
type(atom_dclist), pointer  :: link
integer                     :: i

! Delete from pointer of arrays only once
if(.not.o%delete) then

  ! Flag the atom for deletion
  o%delete=.true.

  ! Delete from array of pointers
  natoms = natoms - 1
  do i=o%id,natoms
    a(i)%o => a(i+1)%o
    a(i)%o%id = i
  enddo
  a(natoms+1)%o=>null()

  ! Deattach from the hard list but keep allocated
  call o%link%deattach()
endif

! Return if there is an intergroup pointing to this atom to avoid segfault in
! that alghoritm.
if (o%nigr>0) return

! WARNING
! El atomo puede quedar en alguna otra lista que no sea un intergrup
! Como por ejemplo en la lista de fantasmas o en la lista de locales

! Deallocate
link => o%link
call o%dest()
call link%destroy_node()
deallocate(link)

end subroutine del_atom

subroutine atoms_increaseby(m)
! Arganda el systema de a pasos de spd, para asegurar m atomos puedan ser
! alojados.
integer,intent(in)             :: m
integer                        :: n,n2
type (atom_ap),allocatable     :: t_a(:)
integer,parameter              :: spd=100

n=m+natoms
n2=n+spd
if(.not.allocated(a)) then
  allocate( a(n2) )
endif

! Agrando el array of pointers si hace falta
if(n>size(a)) then
  allocate( t_a(n2) )
  t_a(1:natoms)    = a(1:natoms)
  call move_alloc(to=a,from=t_a)
endif

endsubroutine

subroutine new_localatom()
type ( atom ), pointer         :: o

! Creating a new atom
call new_atom()
o=>a(natoms)%o
o%tag=natoms

! Adding atom to the local part
call alocal%add_before_soft(o)
nlocal = nlocal + 1

call group_atom_add(o, sys) ! XGHOST

end subroutine new_localatom

subroutine new_ghostatom()
type ( atom ), pointer         :: o

! Creating a new atom
call new_atom()
o=>a(natoms)%o

! Adding atom to the ghost part
call aghost%add_before_soft(o)
nghost=nghost+1
end subroutine new_ghostatom

subroutine atom_dclist_destroyall(node)
! This subrroutine can be replaced by the instrinisc _DestroyAll when final
! procedures be implmented.
type ( atom_dclist ),target     :: node
type ( atom_dclist ), pointer   :: aux

do while(.not.associated(node%next,target=node))
  aux => node%next
  call del_atom(aux%o)
  call aux%deattach()
  call aux%destroy_node()
  deallocate(aux)
enddo

end subroutine atom_dclist_destroyall

subroutine atoms_group_add(g)
! s1 debe ser allocateado.
type(group),intent(in)         :: g
type ( atom_dclist ), pointer  :: la
type ( atom ), pointer         :: o
integer                        :: i

! Although new_atom check for space, this ensure only 1 reallocation.
call atoms_increaseby(g%nat)

la => g%alist
do i = 1,g%nat
  la => la%next

  ! Creating a new atom
  ! New atoms are always added at the end. It does not matters if this are
  ! located after ghost atoms, since the local list always can recover only the
  ! local atoms. It his really matters, then, do the a reassociations of
  ! everything below... this should not change the tag of the atoms if this are
  ! defined at the creation routine
  call new_localatom()
  o=>a(natoms)%o
  call atom_asign(o,la%o)

enddo

end subroutine atoms_group_add

subroutine ghost_group_add(g)
! s1 debe ser allocateado.
type(group),intent(in)         :: g
type ( atom_dclist ), pointer  :: la
type ( atom ), pointer         :: o
integer                        :: i

call atoms_increaseby(g%nat)

la => g%alist
do i = 1,g%nat
  la => la%next

  ! Creating a new atom
  call new_ghostatom()
  o=>a(natoms)%o
  call atom_asign(o,la%o)

enddo

end subroutine ghost_group_add

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
call glist%add_soft(g)

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
call g%alist%destroy_node()
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
if(group_belong(g,a%glist)) return

! agrego el atomo a la lista linkeada del grupo
call g%alist%add_before_soft(a)
g%nat = g%nat + 1 ! numero de particulas

! agrego el grupo a la lista linkeada del atomo.
call a%glist%add_soft(g)

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
call g%alist%add_before_hard()
ln=> g%alist%prev
call ln%o%init()
call atom_asign(ln%o,a)
 
!FIXME: Le doy un id al atomo
ln%o%id = g%nat + 1
                           
g%nat = g%nat + 1 ! numero de particulas
          
! agrego el grupo a la lista linkeada del atomo.
call a%glist%add_soft(g)
 
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
call glist%add_soft(g)
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
      call group_group_del(a%glist,g)

      ! saco el atomo de la lista del grupo
      call la%deattach()
      call la%destroy_node() !soft
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
  call group_group_del(la%o%glist,g)

  ! saco el atomo de la lista del grupo
  next => la%next
  call la%deattach()
  call la%destroy_node() !soft
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
    call lg%destroy_node()
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


subroutine dovirial(vi,virial_tmp)
use gems_constants, only:kcm_ui
real(dp),intent(in) :: virial_tmp(3,3)
integer,intent(in)  :: vi(:)
real(dp)            :: aux(6), frac, vtmp
integer             :: i, j
integer             :: n

if(.not.b_gvirial) return

n=size(vi)

aux(1) = virial_tmp(1,1)
aux(2) = virial_tmp(2,2)
aux(3) = virial_tmp(3,3)
aux(4) = virial_tmp(1,2)
aux(5) = virial_tmp(1,3)
aux(6) = virial_tmp(2,3)
virial(:,:)=virial(:,:)-virial_tmp(:,:)*kcm_ui

  ! FIXME: Creo que debería hacer esto
  ! virial(2,1)=virial(1,2)
  ! virial(3,1)=virial(1,3)
  ! virial(3,2)=virial(2,3)

if (.not.b_avirial) return

frac = 1._dp/n
do j = 1,6
  vtmp = aux(j)*frac
  do i=1,n
    atomvirial(j,vi(i)) = atomvirial(j,vi(i)) + vtmp
  end do
end do

end subroutine dovirial


end module gems_program_types

