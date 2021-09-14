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



!SETS DE PARAMETROS HOMOATOMICOS:
!   (  Co-Co; Cu-Cu;  Pd-Pd; Ag-Ag; Pt-Pt; Au-Au   )

!Para Interacciones HOMO, se tiene que :

!r0=Parametro de red / sqrt(2)  Primeros Vecinos
!rcut_1=r0*sqrt(2)              Segundos Vecinos (Donde se da el Empalme con el polinomio soft-cutt)
!rcut_2=r0*sqrt(3)              Terceros Vecinos (Donde el polinomio soft-cutt se hace cero)

!SETS DE PARAMETROS HETEROATOMICOS:
!   (  Co-Ag; Ni-Ag; Cu-Ag; Pd-Ag;  )
!   (  Co-Au; Pt-Au                 )

!Para interacciones HETERO, se tiene que :

! r0=[r0(A-A)+r0(B-B)] / 2
! rcut_1= Parametro de red (A)
! rcut_2= Parametro de red (B)
! Si se cumple que : Parametro de red (A) < Parametro de red (B)


module gems_tb
use gems_program_types
use gems_constants
use gems_constants
use gems_algebra
use gems_tables
use gems_neighbor
!$ use omp_lib
implicit none

private

! NOTA: Si se desea declarar un tipo derivado del integroup que contenga los
! parametros de la interaccion, se tiene el problema de que la subrrutina
! interact tiene que tener como passed argument la clase ngroup, forzando a
! usar select type dentro de ella. Se podría declarar al ngroup como tipo
! abstract, permitiendo hacer que interact sea un deferred procedure, pero en
! este caso todas las interacciones deben tener un tipo derivados, ya que los
! abstract types no pueden ser instanciados.
! ver http://stackoverflow.com/a/25413107/1342186
! El problema de esto es que cada procedimiento de interact fuerza a definir un
! nuevo tipo derivado ya que no se puede reapuntar el interact puesto que no es
! puntero

type,public,extends(ngroup) :: smatb

  ! The number of internal atom types
  integer     :: nz = 0

  ! The global-internal maping of atom types
  integer     :: z(118) = 0

  ! Enable/disable interaction for each pair of atom types
  logical,allocatable  ::  prm(:,:)

  ! One set of parameters for each pair of atom types
  real(dp),allocatable ::  ca(:,:) ,&
                           eps(:,:),&
                           p(:,:)  ,&
                           q(:,:)  ,&
                           r0(:,:)

  ! Smooth function for each pair of atom types
  real(dp),allocatable ::  rci(:,:),&
                           rce(:,:),&
                           b1(:,:) ,&
                           b2(:,:) ,&
                           b3(:,:) ,&
                           r1(:,:) ,&
                           r2(:,:) ,&
                           r3(:,:)

  ! NOTE: 2 smatb object does not shares their density
  ! This should be sys size, so can be references by a%id
  real(dp),allocatable :: eband(:), band(:)

  contains

  procedure :: attach_atom => smatb_attach

  procedure :: interact => smatb_interact
  procedure,nopass :: cli => smatb_cli

end type

contains

subroutine smatb_attach(g,a)
class(smatb),target        :: g
class(atom),target         :: a
type(atom_dclist),pointer  :: la
integer                    :: i,n,k  

! Save current atom number
n=g%nat
 
! Attempt to attach
call g%ngroup_attach_atom(a)

! Return if atom was already in the group
if(n==g%nat) return
                                  
! Search new internal atom types
la=>g%alist
do i = 1,g%nat
  la=>la%next
  k=la%o%z

  if(g%z(k)/=0) cycle
  g%nz=g%nz+1
  g%z(k)=g%nz
enddo

n=g%nz
  
! XXX: Just using more memory to avoid this operations
! n=0.5_dp*n*(n+1)

if(allocated(g%prm)) then
  if(size(g%prm,1)<n) then
    deallocate(g%prm)
    deallocate(g%ca, g%eps)
    deallocate(g%p, g%q, g%r0)
    deallocate(g%rci, g%rce)
    deallocate(g%b1, g%b2)
    deallocate(g%b3, g%r1)
    deallocate(g%r2, g%r3)
  endif
endif

if(.not.allocated(g%prm)) then
  allocate(g%prm(n,n))
  g%prm(:,:)=.false.

  allocate(g%ca(n,n))
  allocate(g%eps(n,n))
  allocate(g%p(n,n))
  allocate(g%q(n,n))
  allocate(g%r0(n,n))

  allocate(g%rci(n,n))
  allocate(g%rce(n,n))

  allocate(g%b1(n,n))
  allocate(g%b2(n,n))
  allocate(g%b3(n,n))
  allocate(g%r1(n,n))
  allocate(g%r2(n,n))
  allocate(g%r3(n,n))
endif
                        
! Reallocate if needed
if(allocated(g%band)) then
  n=size(g%band)
  if(g%nat<n) return
  deallocate(g%band)
  deallocate(g%eband)
endif
  
n=g%nat+g%pad
allocate(g%band(n))
allocate(g%eband(n))
 
end subroutine
 
subroutine smatb_set(g,i,j,a,eps,p,q,r0,rci,rce)
type(smatb)         :: g
integer,intent(in)  :: i,j
real(dp),intent(in) :: a,eps,p,q,r0,rce,rci
real(dp)            :: dcut,ab,bb,cb,aux

g%prm(i,j)  = .true.

g%ca(i,j)  = a
g%eps(i,j) = eps
g%p(i,j)   = p
g%q(i,j)   = q
g%r0(i,j)  = r0

!Polinomio de suavizado.
!Logsdail etal. RSC Advances 2(2012)5863. 10.1039/c2ra20309j (en el supporting esta el polinomio)
!Logsdail dice que este polinomio machea el suavisado de Baletto etal. (2002). JCP, 116(9), 3856–3863. 10.1063/1.1448484. En este
!paper Balleto menciona que usa entre segundos y terceros vecinos.
!En Rapallo etal. JCP 122(2005)194308. 10.1063/1.1898223 (en este ultimo dice que rext es el 3er vecino)
g%rci(i,j) = rci
g%rce(i,j) = rce

! Aleaciones vieja tirada: Como estaba antes
! rci = 4.0729350596D0 !rci del Au
! rce = 4.4271218640D0

! Polinomial coefficient for band energy (ver supporting info of Logsdail et al.).
! Notese aqui que aux es la raiz de la energía de banda (ver supporting info of Logsdail et al.).
dcut=rce-rci
ab=-   1._dp/dcut**3
bb=-(q/r0)   /dcut**2
cb=-(q/r0)**2/dcut
aux=eps*dexp(-q*(rci/r0-1._dp))
g%b1(i,j)=aux*(12._dp*ab-6._dp*bb+cb)/(2._dp*dcut**2)
g%b2(i,j)=aux*(15._dp*ab-7._dp*bb+cb)/dcut
g%b3(i,j)=aux*(20._dp*ab-8._dp*bb+cb)/2._dp

! Polynomia coeffciens for repulstion energy (ver supporting info of Logsdail et al.).
aux=a*dexp(-p*(rci/r0-1._dp))
bb=-(p/r0)   /dcut**2
cb=-(p/r0)**2/dcut
g%r1(i,j)=aux*(12._dp*ab-6._dp*bb+cb)/(2._dp*dcut**2)
g%r2(i,j)=aux*(15._dp*ab-7._dp*bb+cb)/dcut
g%r3(i,j)=aux*(20._dp*ab-8._dp*bb+cb)/2._dp

! In case of alloys
g%prm(j,i) = g%prm(i,j)

g%ca(j,i)  = g%ca(i,j)
g%eps(j,i) = g%eps(i,j)
g%p(j,i)   = g%p(i,j)
g%q(j,i)   = g%q(i,j)
g%r0(j,i)  = g%r0(i,j)

g%rci(j,i) = g%rci(i,j)
g%rce(j,i) = g%rce(i,j)

g%b1(j,i)  = g%b1(i,j)
g%b2(j,i)  = g%b2(i,j)
g%b3(j,i)  = g%b3(i,j)
g%r1(j,i)  = g%r1(i,j)
g%r2(j,i)  = g%r2(i,j)
g%r3(j,i)  = g%r3(i,j)

end subroutine
    
subroutine smatb_cli(g)
use gems_input_parsing
class(ngroup),intent(inout)  :: g
real(dp)                  :: a,eps,p,q,r0,rci,rce
integer                   :: i,j
character(:),allocatable  :: w1

select type(g)
type is(smatb)  
  ! Set Parameters
  call reada(w1)
  select case(w1)
  case('read')
    call reada(w1)
    call smatb_readprm(g,w1)
  case default
    call reread(0)
    call readf(a)
    call readf(eps)
    call readf(p)
    call readf(q)
    call readf(r0)  ! Angstroms
    call readf(rci) ! Angstroms
    call readf(rce) ! Angstroms
    do i=1,g%nz-1
      do j=i+1,g%nz
        call smatb_set(g,i,j,a,eps,p,q,r0,rci,rce)
      enddo
    enddo
  end select

  ! Radio de corte
  call g%setrc(maxval(g%rce))

class default
  call werr('Interaction type mismatch. Expected smatb type.')
end select  
 
end subroutine smatb_cli

subroutine smatb_readprm(g,prmfile)
! Subrroutine to read the parameter file
use gems_input_parsing, only:opts, input_options, gems_iopts
use gems_elements, only:inq_z
use gems_errors, only:wlog
use gems_algebra, only:sort_int
type(smatb)                 :: g
character(*),intent(in)     :: prmfile
type(input_options), target :: iopts
character(:),allocatable    :: clase,w1,w2
real(dp)                    :: a,eps,p,q,r0,rci,rce
integer                     :: u,i,j

! Redirigiendo el input parsing a las opciones locales
u = find_io(30)
open(u,action='read',file=prmfile)
call iopts%init(skip_blank_lines=.true.,in=u,echo_lines=.false.,comments='!')
opts=>iopts

do
  call read_line(eof)
  if (eof) exit

  call reada(w1)
  if (w1(1:1)=='!') cycle

  ! Selecciono la clase
  selectcase(w1)
  case('GUPTA','SMATB')
    clase=w1
    cycle
  endselect

  selectcase(clase)
  case('GUPTA','SMATB')

    call reada(w2) !nomb2
    call readf(a)
    call readf(eps)
    call readf(p)
    call readf(q)
    call readf(r0)  ! Angstroms
    call readf(rci) ! Angstroms
    call readf(rce) ! Angstroms

    ! Getting internal atom types
    i=g%z(inq_z(w1))
    j=g%z(inq_z(w2))

    ! Cycle if this entree is not relevant for g
    if(i*j==0) cycle

    ! XXX: Just using more memory to avoid this operations
    ! ! Solving symmetry in bonds
    ! k=min(i,j)
    ! j=max(i,j)
    ! i=k
    !
    ! ! Getting pair index
    ! n=j+(i-1)*g%nz-i*(i-1)/2

    ! Add parameters to the corresponding bondgr
    call smatb_set(g,i,j,a,eps,p,q,r0,rci,rce)

  endselect
enddo

! Devuelvo el input parsing al gems
opts=>gems_iopts
close(u)

! FIXME: This warning is always true if the interaction is between two groups
call wwan('No parameters found for some atom types',.not.any(g%prm))

end subroutine smatb_readprm

subroutine smatb_preinteraction(g)
use gems_groups, only: atom, atom_dclist
class(smatb)                :: g
real(dp)                    :: dr,vd(dm),banda,r0
real(dp)                    :: drm,drm2,drm3,drm4,drm5
integer                     :: z1,z2
integer                     :: i,ii,j,l
type(atom),pointer          :: o1,o2
type(atom_dclist),pointer   :: la

! Set zeros
g%band(:)=0.0_dp
g%eband(:)=0.0_dp

! Interaccion

! Sobre los atomos
la => g%ref%alist
do ii = 1,g%ref%nat
  la => la%next
  o1 => la%o
  i = o1%gid(g)
  z1 = g%z(o1%z)
 
  ! sobre los vecinos
  do l = 1, g%nn(i)

    j  = g%list(i,l)
    o2 => g%a(j)%o
    z2 = g%z(o2%z)

    vd = vdistance( o2, o1, mic)
    dr =  sqrt(dot_product(vd,vd))

    if(dr>g%rce(z1,z2)) cycle

    if(dr<g%rci(z1,z2))then

      r0=g%r0(z1,z2)

      banda=g%eps(z1,z2)*g%eps(z1,z2)*exp(-2.0_dp*g%q(z1,z2)*(dr/r0-1.0_dp))

    else

      drm=dr-g%rce(z1,z2)
      drm2=drm *drm
      drm3=drm2*drm
      drm4=drm3*drm
      drm5=drm4*drm

      banda=g%b1(z1,z2)*drm5+g%b2(z1,z2)*drm4+g%b3(z1,z2)*drm3
      banda=banda*banda

    end if

    !OMP: Un thread puede estar accediendo al mismo sector del indice i a
    !partir del indice j mas abajo. Me pregunto si es necesario ponerlo aca
    !o no
    !!$OMP ATOMIC
    g%band(i)=g%band(i)+banda

  enddo

enddo

! Set the band for the ghost
la => g%b%alist
do l=1,g%b%nat
  la => la%next

  o1=>la%o
  o2=>o1%ghost
  if(.not.associated(o2)) cycle

  i = o1%gid(g)
  j = o2%gid(g)
  g%band(i) = g%band(j)

enddo

! eband=band**(-0.5_dp)
! !!$OMP WORKSHARE
g%eband(:)=1.0_dp/sqrt(g%band(:))
! !!$OMP END WORKSHARE

end subroutine smatb_preinteraction
 
subroutine smatb_interact(g)
class(smatb),intent(inout)  :: g
real(dp)                    :: dr,vd(dm),r0,aux1,aux2,lev_ui,repul
real(dp)                    :: drm,drm2,drm3,drm4,drm5,factor,factor2(dm),epot
integer                     :: z1,z2
integer                     :: i,ii,j,l,n
type(atom),pointer          :: o1,o2
type(atom_dclist),pointer   :: la

call smatb_preinteraction(g)

! Make ev_ui inside the lexical extent of OMP
lev_ui=ev_ui

epot=0._dp

la => g%ref%alist
do ii = 1,g%ref%nat
  la => la%next
  o1 => la%o
  i = o1%gid(g)
  z1 = g%z(o1%z)

  do l = 1, g%nn(i)  ! sobre los vecinos

    j  = g%list(i,l)
    o2 => g%a(j)%o
    z2 = g%z(o2%z)

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr =  sqrt(dot_product(vd,vd))

    factor=0.0_dp

    if(dr>g%rce(z1,z2)) cycle

    if(dr<=g%rci(z1,z2))then

      r0=g%r0(z1,z2)

      ! Energia de repulsion
      repul=g%ca(z1,z2)*exp(-g%p(z1,z2)*(dr/r0-1.0_dp))

      ! Fuerza de repulsion y de banda
      factor=2.0_dp*repul*(g%p(z1,z2)/r0)                 &
            -g%eps(z1,z2)*g%eps(z1,z2)*1.0_dp/2.0_dp    &
             *exp(-2.0_dp*g%q(z1,z2)*(dr/r0-1.0_dp))      &
             *(2.0_dp*g%q(z1,z2)/r0)*(g%eband(i)+g%eband(j))

    else if(dr<=g%rce(z1,z2))then

      drm=dr-g%rce(z1,z2)
      drm2=drm *drm
      drm3=drm2*drm
      drm4=drm3*drm
      drm5=drm4*drm

      ! Fuerza de banda
      factor= (  g%b1(z1,z2)*drm5 &
                +g%b2(z1,z2)*drm4 &
                +g%b3(z1,z2)*drm3)&
             *(  5.0_dp*g%b1(z1,z2)*drm4 &
                +4.0_dp*g%b2(z1,z2)*drm3 &
                +3.0_dp*g%b3(z1,z2)*drm2)&
             *(g%eband(i)+g%eband(j))

      ! Fuerza de repulsion
      factor=factor-2.0_dp*(5.0_dp*g%r1(z1,z2)*drm4+4.0_dp*g%r2(z1,z2)*drm3+3.0_dp*g%r3(z1,z2)*drm2)

      ! Energia de repulsion
      repul=g%r1(z1,z2)*drm5+g%r2(z1,z2)*drm4+g%r3(z1,z2)*drm3
    end if

    factor2(1:dm) = factor *lev_ui * vd(1:dm)/dr

    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm)
    o1%epot = o1%epot + repul*lev_ui
    epot = epot + repul*lev_ui


    if (associated(o2%ghost)) then
      if (o2%gid(g)>o1%gid(g)) then

        o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
        o2%epot = o2%epot + repul*lev_ui
        epot = epot + repul*lev_ui

        if(b_gvirial) then
          do n = 1,3
            virial(n,:) = virial(n,:) + box(n)*floor(o2%pos(n)*one_box(n))*factor2(:)
          enddo
        endif

      endif
    endif

    !!$OMP END CRITICAL

  enddo

  aux2=sqrt(g%band(i))*lev_ui

  ! Energia de coohesion por atomo
  o1%epot = o1%epot - aux2

  ! Energia total
  epot = epot - aux2


enddo

g%epot=epot

end subroutine

end module gems_tb
