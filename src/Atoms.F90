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


module gems_atoms
use gems_constants,only:sp,dp,dm,find_io
use gems_elements,only:elements,ncsym,element,inq_z

implicit none
private

! TODO: Atoms deberia tener una jerarquia parecida a:
! atomo mecanico (solo r, epot y f)
! atomo dinamico (mecanico con v, a, m, etc)
! de esta manera los intergroups tendrian la posibilidad
! de tener un atomo mecanico

type, public :: atom

  ! Para llevar un conteo de cuantos links a un atomo hay
  ! Algo interesante aca http://stackoverflow.com/a/38363298/1342186

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
  integer                :: id=0
  
  ! El id de el vector join... que no se que es
  integer                :: idv=0

  ! Identificador de atomo en el sistema
  integer                :: tag=0 
                               
  ! Unique atom id. See `atoms_uid` below.
  integer                :: uid=0 
                               
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

  !procedure :: del => atom_atom_del
  !procedure :: delall => atom_allatom_del
  !procedure :: addatom => atom_include
  !procedure :: atom_asign
  !generic   :: assignment(=) => atom_asign

  procedure :: addigr => atom_addigr
  procedure :: deligr => atom_deligr

end type atom

! Global counter for each atom type ever allocated.
! It allows for a uid per atom object. See `uid` in atom type.
integer :: atoms_uid=0

! Double Circular Linked List of atoms
#define _NODE atom_dclist
#define _CLASS class(atom)
#include "cdlist_header.inc"
public :: atom_dclist
public :: atom_dcl_find
                                  
! Array of Pointers to atoms
#define _NODE atom_ap
#define _CLASS class(atom)
#include "arrayofptrs_header.inc"
public :: atom_ap


! Module procedures 

interface atom_setelmnt
  module procedure atom_setelmnt_bysym,atom_setelmnt_byz
end interface
public :: atom_setelmnt,atom_asign

contains

#define _NODE atom_dclist
#define _CLASS class(atom)
#include "cdlist_body.inc"

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

! Unique id 
atoms_uid=atoms_uid+1
a%uid=atoms_uid

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

function atom_dcl_find(alist,a) result(find)
! Check if a is in alist.
type(atom_dclist),target     :: alist
class(atom),target           :: a
type(atom_dclist),pointer    :: find

find => alist
do 
  find => find%next

  if(associated(find,alist)) exit

  ! if(find%hard) then
  !   if (find%o==a) return ! XXX: Couldn't compare
  ! else
  !   if (associated(find%o,target=a)) return
  ! endif

  if (a%uid==find%o%uid) return
enddo

find => null()

end function atom_dcl_find
         


end module gems_atoms
 
