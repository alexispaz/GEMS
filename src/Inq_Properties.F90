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

 
module gems_inq_properties
! conjuntos de subrutinas inq_
! energy : calculate epot and ekin and the same for each sgrp(j)
! angular_energy : 
! angular_mechanic :
! mass_center : calculate cm_pos and cm_vel and the same for each sgrp(j)
! pbc : efectuate the pbc
! disrad :

use gems_program_types, only: dt,tbox,box,one_box,mic
use gems_groups, only: group,atom,atom_dclist,vdistance
use gems_constants,     only: sp,dp,kB_ui,ui_ev,pi,pio2,dm
use gems_algebra
use gems_input_parsing
use gems_tables

implicit none

public
private :: group,atom,atom_dclist        ! program_types
private :: dt                   ! program_types
private :: dm,tbox,box,one_box  ! program_types
private :: dp,kB_ui,pio2,ui_ev,pi  ! constants

save

!All subroutines or all functions
!interface group_inq_cmpos
!  module procedure atom_dclist_inq_cmpos, group_inq_cmpos
!end interface

contains

             
! element
!--------
  
function inq_pure(alist)
  ! indica si el grupo es puro su elemento, sino devuelve 0
  type(atom_dclist),target,intent(in) :: alist
  type(atom_dclist),pointer           :: la
  integer                             :: inq_pure

  la => alist
  inq_pure=la%next%o%z
  do while (.not.associated(la%next,target=alist))
    la => la%next
    if (la%o%z/=inq_pure) then
      inq_pure=0
      exit
    endif
  enddo

end function inq_pure
             
! atoms dinstances
!-----------------

function inq_mayordistance(g1,g2) result(res)
! Busca la maxima distancia entre los atomos del grupo g1 y los del grupo g2
class(group),intent(in)    :: g1,g2
type(atom_dclist),pointer :: la,lb
real(dp)                  :: rd,res
integer                   :: i,j

res=1.e-10_dp

la => g1%alist 
do i = 1,g1%nat
 la=>la%next

 lb => g2%alist 
 do j = 1,g2%nat
   lb=>la%next
   
   rd =  rdistance2( lb%o, la%o )
   res = max(res,rd)

 enddo     
enddo     

res=sqrt(res)
    
end function inq_mayordistance

subroutine inq_pos_v(g,v)
! calcula las coordenadas referidas al punto v
class(group)               :: g
type(atom_dclist),pointer :: la
integer                   :: i
real(dp)                  :: v(dm)

la => g%alist%next
do i = 1,g%nat
  la%o%pos_v = la%o%pos - v
  la => la%next
enddo
end subroutine inq_pos_v

subroutine inq_vel_v(g,v)
! calcula las coordenadas referidas al punto v
class(group)               :: g
type(atom_dclist),pointer :: la
real(dp)                  :: v(dm)
integer                   :: i
la => g%alist%next
do i = 1,g%nat
  la%o%vel_v = la%o%vel - v
  la => la%next
enddo
end subroutine inq_vel_v

function vrs(i,j) result(vd)
!devuelve un versor que va del atomo i al j
real(dp),dimension(dm)  :: vd
type(atom),intent(in)   :: i,j

vd = vdistance(i,j, mic)
vd = vd/sqrt(dot_product(vd,vd))

end function vrs

function rdistance(i,j) result(rd)
!calculates the distance of two atoms also in pbc case
real(dp)                 :: rd
real(dp),dimension(dm)   :: vd
type(atom),intent(in)   :: i,j

vd=vdistance(i,j, mic)
rd=sqrt(dot_product(vd,vd))

end function rdistance

function rdistance2(i,j) result(rd)
!calculates the distance of two atoms also in pbc case
real(dp)                 :: rd
real(dp),dimension(dm)   :: vd
type(atom),intent(in)   :: i,j

vd=vdistance(i,j, mic)
rd=dot_product(vd,vd)

end function rdistance2

! De nuevo
! --------

subroutine bond_order_groups(g1,g2,r1,r2)
! Compute the bond order number of atoms in "g1" respect to atoms in group
! "g2". "r1" and "r2" are the cut radious of the "bond_order" function.
use gems_algebra, only: bond_order
class(group),intent(in)       :: g2
class(group),intent(inout)    :: g1
real(dp),intent(in)          :: r1,r2
real(dp)                     :: vd(dm),rd
type (atom_dclist),pointer   :: la,lb
integer                      :: i,j

la => g1%alist
do i = 1,g1%nat
  la => la%next

  la%o%border=0.0_dp

  lb => g2%alist
  do j = 1,g2%nat
    lb => lb%next

    if(associated(lb%o,target=la%o)) cycle

    vd = vdistance(lb%o,la%o, mic)
    rd = sqrt(dot_product(vd,vd))
    la%o%border=la%o%border+bond_order(rd,r1,r2)
  enddo

enddo

end subroutine bond_order_groups
                
subroutine bond_order_list(a,b,r1,r2)
! Compute the bond order number of atoms in "g1" respect to atoms in group
! "g2". "r1" and "r2" are the cut radious of the "bond_order" function.
use gems_algebra, only: bond_order
real(dp),intent(in)          :: r1,r2
real(dp)                     :: vd(dm),rd
type(atom_dclist),target     :: a,b
type (atom_dclist),pointer   :: la,lb


la => a%next
do while(.not.associated(la,target=a))
  la%o%border=0.0_dp

  lb => b%next
  do while(.not.associated(lb,target=b))

    if(associated(lb%o,target=la%o)) cycle

    vd = vdistance(lb%o,la%o, mic)
    rd = sqrt(dot_product(vd,vd))
    la%o%border = la%o%border+bond_order(rd,r1,r2)
    lb => lb%next
  enddo
  la => la%next
enddo

end subroutine bond_order_list
                
! atoms angles
!-------------

#ifdef DIM3
function inq_dihedral(i,j,k,l) result(phi)
!Calcula el angulo dihedro entre 4 atomos.
real(dp)                :: phi
real(dp),dimension(3)   :: vbxc
real(dp),dimension(3)   :: da,db,dc
real(dp)                :: ab,bc,ac,at,bt,ct,axb,bxc,rnum,den,signo,cosphi
type(atom),intent(in)   :: i,j,k,l
integer                 :: m
logical                 :: pbc(3)

da=vdistance(j,i, mic)
db=vdistance(k,j, mic)
dc=vdistance(l,k, mic)

! Producto cruz b con c
vbxc=cross_product(db,dc)
pbc=i%pbc.or.j%pbc
do m = 1,3
  if (pbc(m)) vbxc(m)=vbxc(m)-box(m)*idnint(vbxc(m)*one_box(m))
enddo

! Productos escalares cruzados
ab=dot_product(da,db)
bc=dot_product(db,dc)
ac=dot_product(da,dc)

! Modulo cuadrado de cada distancia
at=dot_product(da,da)
bt=dot_product(db,db)
ct=dot_product(dc,dc)

! Productors vectoriales al cuadrado entre distancias
axb=at*bt-ab*ab
bxc=bt*ct-bc*bc

! 
rnum=ab*bc-ac*bt ! El producto cuadrado de los vectoriales?
den=axb*bxc      ! El producto cuadrado de los vectoriales

!chequeo de que átomos no están en la misma línea
if(den>1.d-12)then
  den=dsqrt(den) ! El producto de los vectoriales

  !coseno del ángulo
  cosphi=rnum/den
  cosphi=dmin1(cosphi,1.0_dp)
  cosphi=dmax1(cosphi,-1.0_dp)

  !signo del ángulo
  signo=dot_product(da,vbxc)

  !definición del ángulo
  phi=dsign(dacos(cosphi),signo) 
endif

end function inq_dihedral
#endif

! energies
!---------

subroutine inq_pot_energy(g)
class(group)               :: g
type(atom_dclist),pointer :: la
integer                   :: i

if(g%b_epot) return

g%epot = 0.0_dp

la => g%alist
do i = 1,g%nat
  la => la%next
  g%epot = g%epot + la%o%epot
enddo

g%b_epot=.true.

end subroutine inq_pot_energy
     
subroutine inq_abskin_energy(g)
class(group)         :: g
type(atom_dclist),pointer   :: la
integer              :: i

if(g%b_ekin) return

g%ekin = 0.0_dp


! Prerequisitos
la => g%alist
do i = 1,g%nat
  la => la%next
  g%ekin = g%ekin + dot_product(la%o%vel,la%o%vel)*la%o%mass
enddo

g%ekin = g%ekin*0.5_dp

g%b_ekin=.true.
end subroutine inq_abskin_energy
    
subroutine inq_kin_energy(g)
class(group)         :: g
type(atom_dclist),pointer   :: la
integer              :: i

if(g%b_ekin) return

g%ekin = 0.0_dp


! Prerequisitos
call inq_cm_vel(g)
call inq_vel_v(g,g%cm_vel)

la => g%alist
do i = 1,g%nat
  la => la%next
  g%ekin = g%ekin + dot_product(la%o%vel_v,la%o%vel_v)*la%o%mass
enddo

g%ekin = g%ekin*0.5_dp

g%b_ekin=.true.
end subroutine inq_kin_energy

subroutine inq_temperature(g)
class(group)         :: g

if(g%b_temp) return

! Prerequisitos
call inq_kin_energy(g)

g%temp  = 2.0_dp * g%ekin / (kB_ui*(dm*g%nat))!-dm)) ! El -dm es cuando la energia interna esta referida al centro de masa y la dinamica conservativa????
g%b_temp=.true.
end subroutine inq_temperature

! pressure
!---------

subroutine inq_virial(g)
use gems_program_types, only:b_gvirial,virial
use gems_groups, only: ghost
class(group)              :: g
type(atom_dclist),pointer :: la
type(atom),pointer        :: o
integer                   :: i,m

! Using equation 27 of 
! Thompson, A. P., Plimpton, S. J., & Mattson, W. (2009). General formulation of
! pressure and stress tensor for arbitrary many-body interaction potentials under
! periodic boundary conditions. Journal of Chemical Physics, 131(15).
! http://doi.org/10.1063/1.3245303

if(b_gvirial) return
if(g%b_virial) return
      
g%virial(:,:) = 0.0_dp

la => g%alist
do m = 1,g%nat
  la => la%next
  o => la%o
  do i = 1,3
    g%virial(i,:) = g%virial(i,:) + o%pos(i)*o%force(:)
  enddo
enddo
virial(:,:) = g%virial(:,:)

! Avoiding repeted terms in the virial:
!
! https://sourceforge.net/p/lammps/mailman/message/35720312/
!
! Consider a 1 D system with a diatomic molecule centered just in the limit of
! the box with lenght L. So the system including ghost atoms looks like:
!
!            0                           L 
! -------2---|---1------------------2---|---1--------
!        \_____/ group
!
!
! Where I also indicated the "group associated with the local cell" of this
! system. Using equation 25 of the Thompshon paper virial will be:
!
! W=x1*f1 + (x2-L)*f2
!
! I did not include the term x2*f2 and (x1+L)*f1 since atoms at x2 and x1+L do
! not interact with any group associated with the local cell.
!
! This terms will not be there with the typical lammps setup (i.e. newton's 3rd
! law enabled for interactions crossing subdomain boundaries). In this case,
! each subdomain crossing interaction is listed only once in the neighbor
! lists. When newton's 3rd law is disabled, no forces are accumulated on
! ghost atoms and the reverse communication will not be needed.

! Entonces, esta froma de calcular el virial supone que la fuerza de los atomos
! ghost solo es almacenada en ellos cuando j>i (cuando el tag del ghost sea
! mayor que el local). De esta manera solo la mitad de las interacciones con
! los atomos ghost es tenida en cuenta.

la => ghost%alist
do m = 1,ghost%nat
  la=>la%next
  o => la%o
  do i = 1,3
    ! if(o%pos(i)>box(i)) g%virial(i,j) = g%virial(i,j) + box(i)*o%force(j)
    ! g%virial(i,:) = g%virial(i,:) + box(i)*floor(o%pos(i)*one_box(i))*o%force(:)
    g%virial(i,:) = g%virial(i,:) + box(i)*o%boxcr(i)*o%force(:)
  enddo
enddo

call wwan('Calculation of virial with MIC is not fully implemented in all cases.',mic)
! call wwan('Calculation of virial with MIC is not implemented yet. If 
!     no periodic atom is included in the virial requested, then
!     the calculation is fine and ignore this message. If there is periodic
!     atoms but the box is large enough to avoid interaction with the replicas,
!     then the virial calculation is also fine and ignore this message')

! Note that the virial as defined by Thmopson does not include the volume as
! denominator
! g%virial(:,:) = g%virial(:,:)/box_vol

g%b_virial=.true.

end subroutine inq_virial
 
subroutine inq_pressure(g)
use gems_program_types, only:box_vol,b_gvirial,virial
class(group)              :: g
type(atom_dclist),pointer :: la
type(atom),pointer        :: a
integer                   :: i,j,m

if(.not.b_gvirial) then
  call inq_virial(g)
  virial(:,:)=g%virial(:,:)
endif

! call inq_kin_energy(g)

if(g%b_pressure) return

g%pressure(:,:) = 0.0_dp

la => g%alist
do m = 1,g%nat
  la => la%next
  a => la%o
  do i = 1,3
    do j = i,3
      g%pressure(i,j) =  g%pressure(i,j) + a%vel(i)*a%vel(j)*a%mass
    enddo
  enddo
enddo

do i = 2,3
  do j = 1,i-1
    g%pressure(i,j) = g%pressure(j,i)
  enddo
enddo

do i = 1,3
  do j = 1,3
    g%pressure(i,j) = (g%pressure(i,j) + virial(i,j)) /box_vol
  enddo
enddo

g%b_pressure=.true.
         

end subroutine inq_pressure
 

! center of mass
!---------------
 
subroutine inq_mass(g)
class(group)         :: g
type(atom_dclist),pointer   :: la
integer              :: i

if(g%b_mass) return

g%mass = 0.0_dp

la => g%alist%next
do i = 1,g%nat
  g%mass = g%mass + la%o%mass
  la => la%next
enddo

g%b_mass=.true.
end subroutine inq_mass

subroutine group_inq_cmpos(g)
class(group)         :: g
type(atom_dclist),pointer   :: la
integer              :: i

if(g%b_cm_pos) return

g%cm_pos(1:dm) = 0.0_dp

! Prerequisitos
call inq_mass(g)

la => g%alist
do i = 1,g%nat
  la => la%next
  g%cm_pos(1:dm) = g%cm_pos(1:dm) + la%o % mass * la%o%pos(1:dm)
enddo
g%cm_pos(1:dm) = g%cm_pos(1:dm)/g%mass

g%b_cm_pos=.true.
end subroutine


subroutine group_inq_rg(g)
class(group)         :: g
type(atom_dclist),pointer   :: la
integer              :: i
real(dp)             :: sumaRg,vd(dm)

if(g%b_rg_pos) return

g%rg_pos = 0.0_dp
sumaRg= 0.0_dp
! Prerequisitos
call group_inq_cmpos(g)

la => g%alist
do i = 1,g%nat
  la => la%next
  vd(1:dm)=la%o%pos(1:dm)-g%cm_pos(1:dm)
   sumaRg= sumaRg + dot_product(vd,vd)
enddo
g%rg_pos =sqrt(sumaRg / real (g%nat,dp))

g%b_rg_pos=.true.
end subroutine




function atom_dclist_inq_cmpos(alist,mass) result(rcm)
type(atom_dclist),target,intent(in)  :: alist
real(dp),intent(out),optional        :: mass
real(dp)                             :: rcm(dm),m
type(atom_dclist),pointer            :: la

rcm=0._dp
m=0._dp

la => alist%next
do while(.not.associated(la,target=alist))
  rcm = rcm + la%o%mass * la%o%pos
  m=m+la%o%mass
  la => la%next
enddo    

rcm = rcm/m

if(present(mass)) mass=m

end function

subroutine inq_cm_vel(g)
class(group)         :: g
type(atom_dclist),pointer   :: la
integer              :: i

if(g%b_cm_vel) return

! Prerequisitos
call inq_mass(g)

g%cm_vel(1:dm) = 0.0_dp
la => g%alist
do i = 1,g%nat
  la => la%next
  g%cm_vel(1:dm) = g%cm_vel(1:dm) + la%o%mass * la%o%vel(1:dm)
enddo
g%cm_vel(1:dm) = g%cm_vel(1:dm)/g%mass

g%b_cm_vel=.true.
end subroutine


! geometry
!--------- 
  
subroutine inq_cg(g)
class(group)         :: g
type(atom_dclist),pointer   :: la
integer              :: i

if (g%b_cg_pos) return

g%cg_pos = 0.0_dp
la => g%alist%next
do i = 1,g%nat
  g%cg_pos = g%cg_pos + la%o%pos
  la => la%next
enddo
g%cg_pos = g%cg_pos/g%nat

g%b_cg_pos=.true.
end subroutine
           
subroutine inq_refers_cg(g)
class(group)                 :: g
type(atom_dclist),pointer   :: la
integer                     :: i

! Prerequisitos
call inq_cg(g)

la => g%alist%next
do i = 1,g%nat
  la%o%pos_v = la%o%pos - g%cg_pos
  la => la%next
enddo
end subroutine inq_refers_cg
           
subroutine inq_boundingbox(g)
class(group)                 :: g

if (g%b_maxpos.and.g%b_minpos) return

call get_boundingbox(g%alist,g%maxpos,g%minpos,g%nat)

g%b_maxpos=.true.
g%b_minpos=.true.
end subroutine inq_boundingbox

subroutine get_boundingbox(dclist,maxpos,minpos,nat)
! Me pregunto cual de los dos metodos de loop de lista
! es mas rapido. TODO: Escribir un ejemplo en el otro proyecto de estructuras
class(atom_dclist),target   :: dclist
real(dp),intent(out)        :: maxpos(dm),minpos(dm)
integer,optional,intent(in) :: nat
type(atom_dclist),pointer   :: la
integer                     :: i,k

maxpos(:)=1.0e-8_dp
minpos(:)=1.0e8_dp

if(present(nat)) then
  la => dclist
  do i = 1,nat
    la => la%next
    do k = 1,dm
      if (la%o%pos(k)>maxpos(k)) maxpos(k) = la%o%pos(k)
      if (la%o%pos(k)<minpos(k)) minpos(k) = la%o%pos(k)
    enddo
  enddo
  return
endif

la => dclist
do i = 1,nat
  la => la%next
  if(associated(la,target=dclist)) exit
  do k = 1,dm
    if (la%o%pos(k)>maxpos(k)) maxpos(k) = la%o%pos(k)
    if (la%o%pos(k)<minpos(k)) minpos(k) = la%o%pos(k)
  enddo
enddo

end subroutine get_boundingbox 
! Bond-Order             
! ----------

subroutine inq_bondorder(g,r1,r2)
class(group),intent(in)      :: g
real(dp)                    :: vd(dm),rd,aux
real(dp)                    :: r1,r2
type(atom_dclist),pointer   :: la,lb
integer                     :: i,j


if (r1>r2) then
  aux=r2
  r2=r1
  r1=aux
endif

la => g%alist
do i = 1,g%nat
  la => la%next

  aux=0.0_dp

  lb => g%alist
  do j = 1,g%nat
    lb => lb%next

    if(associated(lb%o,target=la%o)) cycle

    vd = vdistance(lb%o,la%o, mic)
    rd = sqrt(dot_product(vd,vd))
    aux=aux+bond_order(rd,r1,r2)

  enddo

 la%o%border=aux

end do

end subroutine  
! angular properties
! ------------------
subroutine inq_covariance(g)
class(group)                :: g
type(atom_dclist),pointer  :: la
integer                    :: k,l,i

if(g%b_covar) return

g%covar = 0.0_dp
la => g%alist%next
do i = 1,g%nat
  do l=1,dm
    do k=1,dm
      g%covar(l,k) = g%covar(l,k) + la%o%pos(l)*la%o%pos(k)
    enddo                 
  enddo
  la => la%next
enddo
g%covar = g%covar/g%nat

g%b_covar = .true.

end subroutine inq_covariance

! CUIDADO QUE PUEDEN TENER BUGS

subroutine inq_inerciaxy(g)
class(group)               :: g
type(atom_dclist),pointer :: la
integer                   :: k,l,i
real(dp)                  :: aux(2)

if(g%b_inercia) return

g%inercia = 0.0_dp
la => g%alist%next
do i = 1,g%nat
  do l=1,3
    aux=la%o%pos(1:2)
    g%inercia(l,l)=g%inercia(l,l) + dot_product(aux,aux)*la%o%mass
    do k=1,2
      g%inercia(l,k)=g%inercia(l,k) - la%o%pos(k)*la%o%pos(l)*la%o%mass
    enddo
  enddo
  la => la%next
enddo

g%b_inercia = .true.
end subroutine
 
subroutine inq_inercia(g)
class(group)               :: g
type(atom_dclist),pointer :: la
integer                   :: k,l,i

if(g%b_inercia) return

g%inercia = 0.0_dp
la => g%alist%next
do i = 1,g%nat
  do l=1,dm
    g%inercia(l,l)=g%inercia(l,l) + dot_product(la%o%pos,la%o%pos)*la%o%mass
    do k=1,dm
      g%inercia(l,k)=g%inercia(l,k) - la%o%pos(k)*la%o%pos(l)*la%o%mass
    enddo
  enddo
  la => la%next
enddo

g%b_inercia = .true.
end subroutine

subroutine inq_angular_mom(g)
use gems_algebra, only:cross_product
class(group)               :: g
type(atom_dclist),pointer :: la
integer                   :: i

if(g%b_ang_mom) return

! Prerequisitos
call group_inq_cmpos(g)
call inq_cm_vel(g)
call inq_pos_v(g,g%cm_pos)
call inq_vel_v(g,g%cm_vel)

g%ang_mom = 0.0_dp
la => g%alist%next
do i = 1,g%nat
  g%ang_mom = g%ang_mom + cross_product(la%o%pos_v,la%o%vel_v)*la%o%mass
  la => la%next
enddo

g%b_ang_mom = .true.
end subroutine

subroutine inq_angular_vel(g)
class(group)               :: g
real(sp),dimension(3,3)   :: inercia_aux
real(sp),dimension(3,1)   :: ang_vel

if(g%b_ang_mom) return

! Prerequisitos
call inq_angular_mom(g)
call inq_inercia(g)

g%ang_vel = g%ang_mom
inercia_aux=real(g%inercia,sp)
ang_vel(:,1)=real(g%ang_vel,sp)

! FIXME: Necesito una subrrutina equivalente al NR de Gaussj
! call gaussj(inercia_aux,ang_vel)
g%inercia = inercia_aux
g%ang_vel = ang_vel(:,1)
!g%ang_vel = cramers_3x3(g%inercia,g%ang_mom) ! otra opcion

g%b_ang_mom = .true.
end subroutine

! SEGURO TIENE BUGS

subroutine inq_angular_energy(g)
class(group)               :: g
type(atom_dclist),pointer :: la
integer                   :: i
real(dp)                  :: aux(3)

if(g%b_evib.and.g%b_erot) return

! Prerequisitos
call inq_angular_vel(g)
call group_inq_cmpos(g)
call inq_cm_vel(g)
call inq_pos_v(g,g%cm_pos)
call inq_vel_v(g,g%cm_vel)

g%erot = 0.0_dp
g%evib = 0.0_dp
la => g%alist%next
do i = 1,g%nat
  aux = cross_product(g%ang_vel,la%o%pos_v)
  la%o%vel_rot = aux(1:dm)
  la%o%vel_vib = la%o%vel_v-la%o%vel_rot
  la%o%erot_ss = dot_product(la%o%vel_rot,la%o%vel_rot)*la%o%mass
  la%o%evib_ss = dot_product(la%o%vel_vib,la%o%vel_vib)*la%o%mass
  g%evib = g%evib + la%o%evib_ss
  g%erot = g%erot + la%o%erot_ss
  la => la%next
enddo
g%evib = g%evib*0.5_dp
g%erot = g%erot*0.5_dp

g%b_evib = .true.
g%b_erot = .true.
end subroutine

! morfology
!----------

subroutine inq_principal_geometric_axis(g,trns)
!previus calcs needed: g%covar_gc
real(sp)                      :: aux(3,3),aux2(3,3)
real(dp),optional,intent(out) :: trns(3,3)  ! Matriz de transformacion (use rotate)
! integer                       :: nrot
class(group)                   :: g

if (g%b_mainaxis) return

! Prerequisitos
call inq_covariance(g)
                  
aux=real(g%covar,sp)
g%mainaxis=0.d0
! FIXME: Necesito una subrrutina equivalente a jacobi y eigsrt del numerical recipes
! call jacobi(aux,g%mainaxis,aux2,nrot)
! call eigsrt(g%mainaxis,aux2)
if(present(trns)) trns=transpose(aux2)
g%mainaxis = sqrt(g%mainaxis)

!triax_param = princ_axis(2)**2/(princ_axis(1)*princ_axis(dm))

!    prolat_____
!   |          /|o
!  b/c       /  |b
!   |   sphere  |l
!   !    /      |a
!   |  /        |t
!   [/__ a/b____|

g%b_mainaxis = .true.

end subroutine inq_principal_geometric_axis

subroutine inq_principal_mass_axis(g,trns)
!previus calcs needed: g%covar_gc
real(sp)                      :: aux(3,3),aux2(3,3)
real(dp),optional,intent(out) :: trns(3,3)  ! Matriz de transformacion (use rotate)
! integer                       :: nrot
class(group)                   :: g

if (g%b_mainaxis) return

! Prerequisitos
call inq_inercia(g)
                  
aux=real(g%inercia,sp)
g%mainaxis=0.d0
! FIXME: Necesito una subrrutina equivalente a jacobi y eigsrt del numerical recipes
! call jacobi(aux,g%mainaxis,aux2,nrot)
! call eigsrt(g%mainaxis,aux2)
if(present(trns)) trns=transpose(aux2)
g%mainaxis = sqrt(g%mainaxis)

!triax_param = princ_axis(2)**2/(princ_axis(1)*princ_axis(dm))

!    prolat_____
!   |          /|o
!  b/c       /  |b
!   |   sphere  |l
!   !    /      |a
!   |  /        |t
!   [/__ a/b____|

g%b_mainaxis = .true.

end subroutine inq_principal_mass_axis
 
subroutine inq_principal_mass_xyaxis(g,trns)
! use nr,only: jacobi,eigsrt
!previus calcs needed: g%covar_gc
real(sp)                      :: aux(3,3),aux2(3,3)
real(dp),optional,intent(out) :: trns(3,3)  ! Matriz de transformacion (use rotate)
! integer                       :: nrot
class(group)                   :: g

if (g%b_mainaxis) return

! Prerequisitos
call inq_inerciaxy(g)
                  
aux=real(g%inercia,sp)
g%mainaxis=0.d0
! FIXME: Necesito una subrrutina equivalente a jacobi y eigsrt del numerical recipes
! call jacobi(aux,g%mainaxis,aux2,nrot)
g%mainaxis(3)=-1.0e8_dp
! FIXME: Necesito una subrrutina equivalente a jacobi y eigsrt del numerical recipes
! call eigsrt(g%mainaxis,aux2)
if(present(trns)) trns=transpose(aux2)
g%mainaxis = sqrt(g%mainaxis)

!triax_param = princ_axis(2)**2/(princ_axis(1)*princ_axis(dm))

!    prolat_____
!   |          /|o
!  b/c       /  |b
!   |   sphere  |l
!   !    /      |a
!   |  /        |t
!   [/__ a/b____|

g%b_mainaxis = .true.

end subroutine inq_principal_mass_xyaxis

function inq_triaxial_param(g)
real(dp)          :: inq_triaxial_param
class(group)       :: g
 
call inq_principal_geometric_axis(g)
                   
#ifdef DIM3
inq_triaxial_param = g%mainaxis(3)*g%mainaxis(1)/(g%mainaxis(2)*g%mainaxis(2))
#endif
#ifdef DIM2
inq_triaxial_param = g%mainaxis(1)/g%mainaxis(2)
#endif
#ifdef DIM1
inq_triaxial_param = 1.0_dp
#endif
inq_triaxial_param = log(inq_triaxial_param)

end function inq_triaxial_param

function inq_disrad(npoints,long,g)
integer                      :: i,j,npoints,counts(npoints),dist
real(dp)                     :: volfac,dr,dr2,long,long2,rd
real(dp),allocatable         :: inq_disrad(:,:)
real(dp)                     :: vd(dm)
class(group)                  :: g
type(atom_dclist),pointer    :: la,lb

allocate(inq_disrad(2,npoints))

!para no tener que sacar la raiz en cada iteracion
long2=long*long

dr = long/npoints
dr2 = long2/npoints
counts=0

!make the histogram
la => g%alist%next
do i = 1, g%nat - 1

  lb => la%next
  do j = i+1, g%nat
    vd = vdistance(la%o,lb%o, mic)
    rd = dot_product(vd,vd)
    if (rd.lt.long2) then
      dist = int(rd/dr2)+1
      counts(dist) = counts(dist) + 1
    endif
    lb => lb%next
  enddo

  la => la%next
enddo

volfac = 4.0_dp*pi/dm

!normalize the histogram
do i = 1,npoints
  inq_disrad(1,i) = dr*(i-0.5_dp)
  inq_disrad(2,i) = counts(i)/(volfac*(i*dr)**dm)
enddo

inq_disrad(2,:) = inq_disrad(2,:) / maxval(inq_disrad)

end function
      
function inq_pdisrad(p,g) result(pdis)
use gems_program_types, only: distance
real(dp),intent(in)          :: p(dm)
class(group),intent(in)       :: g
real(dp)                     :: pdis(g%nat)
real(dp)                     :: vd(dm)
type(atom_dclist),pointer    :: la
integer                      :: i

la => g%alist
do i = 1, g%nat
  la => la%next
  vd = distance(la%o%pos,p,la%o%pbc)
  pdis(i)=sqrt(dot_product(vd,vd))
enddo

end function inq_pdisrad

!                                              configuration
!----------------------------------------------------------

subroutine inq_max_change(g,nsteps,temp)
class(group)         :: g
integer,intent(in)  :: nsteps
real(dp),intent(in) :: temp
integer             :: i
real(dp)            :: fact
type(atom_dclist),pointer :: la

! Esto tiene sentido en una simulación canonica
! Sabiendo que para una particula:
!   __     ____
!   Ec=m_i*v**2/2=3kT/2
!
! Podemos escribir que
!   _____
!   dx**2=dm*k*T*dt**2/m
!
! El maximo lo encontramos si el atomo en cada paso de
! dinamica es desplazada hacia la misma dirección:
!   ____
!   dx**2=dm*k*T*(dt*nsteps)**2/m
!
! Dado que la masa esta dividiendo esto es una propiedad intrinseca de cada
! atomo y no del grupo (aunque depende de la temepratura de este)..... salvo
! que el grupo sea de un elemento puro, con lo cual se puede definir un
! desplaciamiento maximo en general.

fact = dm*kB_ui*temp*dt*nsteps*dt*nsteps

la => g%alist%next
do i = 1,g%nat
  la%o % maxdisp2 = fact*la%o%one_mass
  la => la%next
enddo

end subroutine inq_max_change

!                                                   select
!----------------------------------------------------------



function inq_far(ctr,g) result(at)
! selecciona el atomo mas cercano al punto
type(atom),pointer         :: at
class(group),intent(in)     :: g
real(dp)                   :: rd,rad,ctr(dm)
type (atom_dclist),pointer :: la
integer                    :: i

at=>null()

call inq_pos_v(g,ctr)

rad = 1.0e-8_dp

la => g%alist%next
do i = 1,g%nat
  rd = dot_product(la%o % pos_v,la%o % pos_v)
  if (rd > rad ) then
    rad=rd
    at=>la%o
  endif
  la => la%next
enddo

end function inq_far
                   


end module gems_inq_properties

!                                          molecular volume
!----------------------------------------------------------

!  function volume(g)
!   class(group)       :: g
!   type(atom_dclist) :: la
!   real(dp)               :: volume
!   logical,allocatable    :: vatom(:,:,:)
!   integer,allocatable    :: vmolec(:,:,:)
!   integer                :: i,j,k,tri(3),diam,resol,rad

!   volume=0.d0
!   resol = 100  ! resolucion de la grilla

   ! allocateo la matrix que representa el espacio donde esta un atomo
   ! esto abria que generalizarlo para que cada tipo de atomo tenga su matrix
   ! incluso se podria hacer mas veloz si tomara solo un cuarto del espacio
   ! y luego lo repitiera con las reflexiones correspondientes
   ! diam es un numero par de elementos y siempre divisible por dos
!   rad = int( ss%e%metalradi*resol )
!   diam = rad*2 
!   allocate( vatom(diam,diam,diam) )
!   vatom = .false.
!   do i = -rad,rad
!     do j =  -rad,rad
!       do k =  -rad,rad
!         if( (i**2+j**2+k**2) <= rad**2 ) vatom(i+rad,j+rad,k+rad) = .true.
!       enddo
!     enddo
!   enddo

   ! allocateo la matrix que representa el espacio donde esta la molecula
!   call max_points(g)
!   tri = int( s%g%maxpos-s%g%minpos*resol )
!   allocate( vmolec(tri(1),tri(2),tri(3)) )

  
   !positionen the atoms in the molecule matrix. can use unpack or merge
!   vmolec = .false.
!   do i = 1,s%nat
!     tri = int( s%a(i)%pos*resol )
!     tri = tri - rad
!     vmolec(tri(1):tri(1)+diam,tri(2):tri(2)+diam,tri(3):tri(3)+diam)=                             & 
!     unpack([1],mask=vatom,field=vmolec(tri(1):tri(1)+diam,tri(2):tri(2)+diam,tri(3):tri(3)+diam))
!   enddo

!   volume=sum(vmolec)*resol  ! [amstrong**3]

!  end function


!

!Widely used definitions of atomic radius include:
! - Van der Waals radius: in principle, half the minimum distance between the
!  nuclei of two atoms of the element that are not bound to the same molecule.
! - Ionic radius: the nominal radius of the ions of an element in a specific
!  ionization state, deduced from the spacing of atomic nuclei in cyrstalline salts of
!  that include that ion. In principle, the spacing between two adjacent oppositely
!  charged ions (the length of the ionic bond between them) should equal the sum of
!  their ionic radii.[3] 
! - Covalent radius: the nominal radius of the atoms of an element when covalently
!  bound to other atoms, as deduced the separation between the atomic nuclei in
!  molecules. In principle, the distance between two atoms that are bound to each
!  other in a molecule (the length of that covalent bond) should equal the sum of
!  their covalent radii.[3]
! - Metallic radius: the nominal!  radius of atoms of an element when joined to
!  other atoms by metallic bonds.
! - Bohr radius: the radius of the lowest-energy electron orbit predicted by Bohr
!  model of the atom (1913). It is only applicable to atoms and ions with a single
!  electron, such as hydrogen, singly-ionized helium, and positronium. Although
!  the model itself is now obsolete, the Bohr radius for the hydrogen atom is still
!  regarded as an important physical constant.
 


 



