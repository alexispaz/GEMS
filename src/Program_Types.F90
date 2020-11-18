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
use gems_atoms
use gems_groups
use gems_constants,only:sp,dp,dm,find_io
use gems_elements,only:elements,ncsym,element,inq_z


implicit none
save
private

! The program and subprograms trace this variable to get the term order
! After the term_signal is read, is set to .false. too allow subprograms end
! without quit gems.
logical,public :: term_signal=.false., cycle_signal=.false.


! Program variables
real(dp),target,public    :: dm_steps=0._dp
real(dp),public    :: nframe=0._dp   ! real para que se banque numeros grandes
real(dp),public    :: ptime,pnframe

! Diferencial de posicion que es considerado cero (usado para minimizar y otras cosas)
real(dp),public,parameter :: dpos=0.2 

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

integer,public,parameter   :: nb_max=150 ! Vecinos maximos
integer,public             :: ncell


! Objetos

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


public :: new_ghostatom, del_atom

public :: ghost_group_add
 
public :: atoms_group_add, atom_distancetopoint, atom_distancetoaxis, atom_dclist_destroyall


logical,public,target,allocatable      :: fix(:),join(:)
integer,public                         :: njoin=0

! Declaraciones de Tipos

type(group),target,public       :: sys     ! Systema total
integer,parameter,public        :: mgr=9
type(group),target,public       :: gr(mgr) ! Selecciones
type(group),target,public       :: gsel    ! grupo seleccionado
type(group),target,public       :: gnew    ! grupo creacion

 
contains

subroutine box_setvars()
use gems_input_parsing

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


subroutine new_atom()
! agrega un nuevo atomo (vacio) al sistema. Luego este puede ser referenciado
! via a(natoms)%o
type ( atom ), pointer         :: o
integer                        :: j

! The atoms are allocated in the `atoms` cdlist
call atoms%add_before()
call atoms%prev%alloc()
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
call alocal%add_before()
call alocal%prev%point(o)
nlocal = nlocal + 1

call group_atom_add(o, sys) ! XGHOST

end subroutine new_localatom

subroutine new_ghostatom()
type ( atom ), pointer         :: o

! Creating a new atom
call new_atom()
o=>a(natoms)%o

! Adding atom to the ghost part
call aghost%add_before()
call aghost%prev%point(o)
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
  deallocate(aux)
enddo

end subroutine atom_dclist_destroyall
     
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


! Distancias

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

 aux=atom_distancetopoint(a,p)
 vd=aux-(dot_product(aux,r))*r

end function atom_distancetoaxis
               

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
            

subroutine atoms_group_add(g)
! Agrego atomos del grupo a la lista de atoms
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
      
end module gems_program_types
 
