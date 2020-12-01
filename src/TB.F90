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

!r0=Parametro de red / dsqrt(2)  Primeros Vecinos
!rcut_1=r0*dsqrt(2)              Segundos Vecinos (Donde se da el Empalme con el polinomio soft-cutt)
!rcut_2=r0*dsqrt(3)              Terceros Vecinos (Donde el polinomio soft-cutt se hace cero)

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
use gems_neighbour
!$ use omp_lib
implicit none

private
public  :: smatb_new

! NOTA: Si se desea declarar un tipo derivado del integroup que contenga los
! parametros de la interaccion, se tiene el problema de que la subrrutina
! interact tiene que tener como passed argument la clase intergroup, forzando a
! usar select type dentro de ella. Se podría declarar al intergroup como tipo
! abstract, permitiendo hacer que interact sea un deferred procedure, pero en
! este caso todas las interacciones deben tener un tipo derivados, ya que los
! abstract types no pueden ser instanciados.
! ver http://stackoverflow.com/a/25413107/1342186
! El problema de esto es que cada procedimiento de interact fuerza a definir un
! nuevo tipo derivado ya que no se puede reapuntar el interact puesto que no es
! puntero

type,extends(intergroup) :: smatb

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
  real(dp),allocatable :: eband(:), band(:)

  contains

  procedure :: interact => smatb_interact
  procedure,nopass :: cli => smatb_cli

end type

contains

subroutine smatb_set(ig,i,j,a,eps,p,q,r0,rci,rce)
type(smatb)         :: ig
integer,intent(in)  :: i,j
real(dp),intent(in) :: a,eps,p,q,r0,rce,rci
real(dp)            :: dcut,ab,bb,cb,aux

ig%prm(i,j)  = .true.

ig%ca(i,j)  = a
ig%eps(i,j) = eps
ig%p(i,j)   = p
ig%q(i,j)   = q
ig%r0(i,j)  = r0

!Polinomio de suavizado.
!Logsdail etal. RSC Advances 2(2012)5863. 10.1039/c2ra20309j (en el supporting esta el polinomio)
!Logsdail dice que este polinomio machea el suavisado de Baletto etal. (2002). JCP, 116(9), 3856–3863. 10.1063/1.1448484. En este
!paper Balleto menciona que usa entre segundos y terceros vecinos.
!En Rapallo etal. JCP 122(2005)194308. 10.1063/1.1898223 (en este ultimo dice que rext es el 3er vecino)
ig%rci(i,j) = rci
ig%rce(i,j) = rce

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
ig%b1(i,j)=aux*(12._dp*ab-6._dp*bb+cb)/(2._dp*dcut**2)
ig%b2(i,j)=aux*(15._dp*ab-7._dp*bb+cb)/dcut
ig%b3(i,j)=aux*(20._dp*ab-8._dp*bb+cb)/2._dp

! Polynomia coeffciens for repulstion energy (ver supporting info of Logsdail et al.).
aux=a*dexp(-p*(rci/r0-1._dp))
bb=-(p/r0)   /dcut**2
cb=-(p/r0)**2/dcut
ig%r1(i,j)=aux*(12._dp*ab-6._dp*bb+cb)/(2._dp*dcut**2)
ig%r2(i,j)=aux*(15._dp*ab-7._dp*bb+cb)/dcut
ig%r3(i,j)=aux*(20._dp*ab-8._dp*bb+cb)/2._dp

! In case of alloys
ig%prm(j,i) = ig%prm(i,j)

ig%ca(j,i)  = ig%ca(i,j)
ig%eps(j,i) = ig%eps(i,j)
ig%p(j,i)   = ig%p(i,j)
ig%q(j,i)   = ig%q(i,j)
ig%r0(j,i)  = ig%r0(i,j)

ig%rci(j,i) = ig%rci(i,j)
ig%rce(j,i) = ig%rce(i,j)

ig%b1(j,i)  = ig%b1(i,j)
ig%b2(j,i)  = ig%b2(i,j)
ig%b3(j,i)  = ig%b3(i,j)
ig%r1(j,i)  = ig%r1(i,j)
ig%r2(j,i)  = ig%r2(i,j)
ig%r3(j,i)  = ig%r3(i,j)

end subroutine
 
subroutine smatb_new(pg,g1,g2)
use gems_inq_properties, only: inq_pure
use gems_groups, only: group
use gems_atoms, only: atom_dclist
type(smatb),pointer             :: ig
class(intergroup),pointer       :: pg
type(group),intent(in)          :: g1
type(group),intent(in),optional :: g2
type(atom_dclist),pointer       :: la
integer                         :: n,k,i

! Return a intergroup class pointer
allocate(ig)
pg=>ig

! Initialize the integroup
if(present(g2)) then
  call ig%init(g1=g1,g2=g2)
else
  call ig%init(g1=g1)
endif

! Set internal types from group 1
n=0
la => g1%alist
do i = 1,g1%nat
  la => la%next
  k=la%o%z
  if(ig%z(k)==0) n=n+1
  ig%z(k)=n
enddo

! Set internal types from group 2
if(present(g2)) then
  la => g2%alist
  do i = 1,g2%nat
    la => la%next
    k=la%o%z
    if(ig%z(k)==0) n=n+1
    ig%z(k)=n
  enddo
endif

ig%nz=n

! XXX: Just using more memory to avoid this operations
! n=0.5_dp*n*(n+1)

allocate(ig%prm(n,n))
ig%prm(:,:)=.false.

allocate(ig%ca(n,n))
allocate(ig%eps(n,n))
allocate(ig%p(n,n))
allocate(ig%q(n,n))
allocate(ig%r0(n,n))

allocate(ig%rci(n,n))
allocate(ig%rce(n,n))

allocate(ig%b1(n,n))
allocate(ig%b2(n,n))
allocate(ig%b3(n,n))
allocate(ig%r1(n,n))
allocate(ig%r2(n,n))
allocate(ig%r3(n,n))

! Alloc the band, eband
allocate(ig%eband(ig%n(4)))
allocate(ig%band(ig%n(4)))

! if(allocated(band)) then
!   call wwan('Note that TB energy of a system can not be expresed as the sum of two subsystems TB energies.'  )
!   return
! endif

end subroutine smatb_new
     
subroutine smatb_cli(ig)
use gems_input_parsing
class(intergroup),intent(inout)  :: ig
real(dp)                  :: a,eps,p,q,r0,rci,rce
integer                   :: i,j
character(len=linewidth)  :: w1

select type(ig)
type is(smatb)  

  ! Parameters default
  call reada(w1)
  select case(w1)
  case('read')
    call reada(w1)
    call smatb_readprm(ig,w1)
  case default
    call reread(0)
    call readf(a)
    call readf(eps)
    call readf(p)
    call readf(q)
    call readf(r0)  ! Angstroms
    call readf(rci) ! Angstroms
    call readf(rce) ! Angstroms
    do i=1,ig%nz-1
      do j=i+1,ig%nz
        call smatb_set(ig,i,j,a,eps,p,q,r0,rci,rce)
      enddo
    enddo
  end select

  ! Radio de corte
  call ig%setrc(maxval(ig%rce))

class default
  call werr('Interaction type mismatch. Expected smatb type.')
end select  
 
end subroutine smatb_cli

subroutine smatb_readprm(ig,prmfile)
! Subrroutine to read the parameter file
use gems_input_parsing, only:opts, input_options, gems_iopts
use gems_elements, only:inq_z
use gems_errors, only:wlog
use gems_algebra, only:sort_int
type(smatb)                 :: ig
character(*),intent(in)     :: prmfile
type(input_options), target :: iopts
character(90)               :: clase,w1,w2
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
    i=ig%z(inq_z(w1))
    j=ig%z(inq_z(w2))

    ! Cycle if this entree is not relevant for ig
    if(i*j==0) cycle

    ! XXX: Just using more memory to avoid this operations
    ! ! Solving symmetry in bonds
    ! k=min(i,j)
    ! j=max(i,j)
    ! i=k
    !
    ! ! Getting pair index
    ! n=j+(i-1)*ig%nz-i*(i-1)/2

    ! Add parameters to the corresponding bondgr
    call smatb_set(ig,i,j,a,eps,p,q,r0,rci,rce)

  endselect
enddo

! Devuelvo el input parsing al gems
opts=>gems_iopts
close(u)

! FIXME: This warning is always true if the interaction is between two groups
call wwan('No parameters found for some atom types',any(ig%prm))

end subroutine smatb_readprm

subroutine smatb_preinteraction(ig)
class(smatb)              :: ig
real(dp)                  :: dr,vd(dm),banda,r0
real(dp)                  :: drm,drm2,drm3,drm4,drm5
integer                   :: z1,z2
integer                   :: i,j,l
type(atom),pointer        :: o1,o2


! TODO MPI: Redimensionar band y eband si aumenta el numero de atomos

! Set zeros
ig%band(:)=0.0_dp
ig%eband(:)=0.0_dp

! Interaccion

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP& FIRSTPRIVATE(ig)    &
!$OMP& PRIVATE(l,z1,j,z2,dr,vd,banda,drm,drm2,drm3,drm4,drm5,r0,o1,o2) &
!$OMP& SHARED(mic,nlocal)

! Sobre los atomos
do i = 1,ig%n(1)

  o1 => ig%at(i)%o
  z1 = ig%z(o1%z)

  ! sobre los vecinos
  do l = 1, ig%nn(i)

    j  = ig%list(i,l)
    o2 => ig%at(j)%o
    z2 = ig%z(o2%z)

    vd = vdistance( o2, o1, mic)
    dr =  dsqrt(dot_product(vd,vd))

    if(dr>ig%rce(z1,z2)) cycle

    if(dr<ig%rci(z1,z2))then

      r0=ig%r0(z1,z2)

      banda=ig%eps(z1,z2)*ig%eps(z1,z2)*exp(-2.0_dp*ig%q(z1,z2)*(dr/r0-1.0_dp))

    else

      drm=dr-ig%rce(z1,z2)
      drm2=drm *drm
      drm3=drm2*drm
      drm4=drm3*drm
      drm5=drm4*drm

      banda=ig%b1(z1,z2)*drm5+ig%b2(z1,z2)*drm4+ig%b3(z1,z2)*drm3
      banda=banda*banda

    end if

    !OMP: Un thread puede estar accediendo al mismo sector del indice i a
    !partir del indice j mas abajo. Me pregunto si es necesario ponerlo aca
    !o no
    !!$OMP ATOMIC
    ig%band(i)=ig%band(i)+banda

    if(.not.ig%newton) cycle
    if(j>nlocal) cycle

    !!$OMP ATOMIC
    ig%band(j)=ig%band(j)+banda

  enddo

enddo

!$OMP END PARALLEL DO

! band for the ghost
do i=1+ig%n(1),ig%n(2)
  j = ig%at(i)%o%tag
  ig%band(i) = ig%band(j)
enddo

! eband=band**(-0.5_dp)
! !!$OMP WORKSHARE
ig%eband(:)=1.0_dp/sqrt(ig%band(:))
! !!$OMP END WORKSHARE

end subroutine smatb_preinteraction

subroutine smatb_interact(ig)
class(smatb),intent(inout)   :: ig
real(dp)                     :: dr,vd(dm),r0,aux1,aux2,lev_ui,repul
real(dp)                     :: drm,drm2,drm3,drm4,drm5,factor,factor2(dm),epot
integer                      :: z1,z2
integer                      :: i,j,l,n
type(atom),pointer           :: o1,o2

call smatb_preinteraction(ig)

! Make ev_ui inside the lexical extent of OMP
lev_ui=ev_ui

epot=0._dp

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP& PRIVATE(i,o1,z1,l,j,o2,z2,vd,dr,n,aux1,aux2,repul)  &
!$OMP& PRIVATE(factor,factor2,r0,drm,drm2,drm3,drm4,drm5)  &
!$OMP& SHARED(ig,lev_ui,b_gvirial,virial)       &
!$OMP& SHARED(mic,nlocal,one_box,box)      &
!$OMP& REDUCTION(+:epot)

do i = 1,ig%n(1)

  o1 => ig%at(i)%o
  z1 = ig%z(o1%z)

  do l = 1, ig%nn(i)  ! sobre los vecinos

    j  = ig%list(i,l)
    o2 => ig%at(j)%o
    z2 = ig%z(o2%z)

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr =  dsqrt(dot_product(vd,vd))

    factor=0.0_dp

    if(dr>ig%rce(z1,z2)) cycle

    if(dr<=ig%rci(z1,z2))then

      r0=ig%r0(z1,z2)

      ! Energia de repulsion
      repul=ig%ca(z1,z2)*exp(-ig%p(z1,z2)*(dr/r0-1.0_dp))

      ! Fuerza de repulsion y de banda
      factor=2.0_dp*repul*(ig%p(z1,z2)/r0)                 &
            -ig%eps(z1,z2)*ig%eps(z1,z2)*1.0_dp/2.0_dp    &
             *exp(-2.0_dp*ig%q(z1,z2)*(dr/r0-1.0_dp))      &
             *(2.0_dp*ig%q(z1,z2)/r0)*(ig%eband(i)+ig%eband(j))

    else if(dr<=ig%rce(z1,z2))then

      drm=dr-ig%rce(z1,z2)
      drm2=drm *drm
      drm3=drm2*drm
      drm4=drm3*drm
      drm5=drm4*drm

      ! Fuerza de banda
      factor= (  ig%b1(z1,z2)*drm5 &
                +ig%b2(z1,z2)*drm4 &
                +ig%b3(z1,z2)*drm3)&
             *(  5.0_dp*ig%b1(z1,z2)*drm4 &
                +4.0_dp*ig%b2(z1,z2)*drm3 &
                +3.0_dp*ig%b3(z1,z2)*drm2)&
             *(ig%eband(i)+ig%eband(j))

      ! Fuerza de repulsion
      factor=factor-2.0_dp*(5.0_dp*ig%r1(z1,z2)*drm4+4.0_dp*ig%r2(z1,z2)*drm3+3.0_dp*ig%r3(z1,z2)*drm2)

      ! Energia de repulsion
      repul=ig%r1(z1,z2)*drm5+ig%r2(z1,z2)*drm4+ig%r3(z1,z2)*drm3
    end if

    factor2(1:dm) = factor *lev_ui * vd(1:dm)/dr

    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm)
    o1%epot = o1%epot + repul*lev_ui
    epot = epot + repul*lev_ui

    if(j<=nlocal) then

      if(.not.ig%newton) cycle

      o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
      o2%epot = o2%epot + repul*lev_ui
      epot = epot + repul*lev_ui

      if(b_gvirial) then
        do n = 1,3
          virial(n,:) = virial(n,:) + o2%pos(n)*factor2(:)
        enddo
      endif

    else if (o2%tag>o1%tag) then

      o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
      o2%epot = o2%epot + repul*lev_ui
      epot = epot + repul*lev_ui

      if(b_gvirial) then
        do n = 1,3
          virial(n,:) = virial(n,:) + box(n)*floor(o2%pos(n)*one_box(n))*factor2(:)
        enddo
      endif

    endif

    !!$OMP END CRITICAL

  enddo

  aux2=dsqrt(ig%band(i))*lev_ui

  ! Energia de coohesion por atomo
  o1%epot = o1%epot - aux2

  ! Energia total
  epot = epot - aux2

enddo

!$OMP END PARALLEL DO

ig%epot=epot

end subroutine

end module gems_tb
