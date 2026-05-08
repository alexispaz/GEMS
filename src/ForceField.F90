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


module gems_forcefield

! FIXME: This module was build without considering that sys index may have
! null atoms due to detaching.

! Existen los boundgr que son grupos de objetos bound. Los bound, son objetos
! que contienen 2,3,4 atomos unidos. Cada boundgr tiene asignada una funcion y
! parametros que va a aplicar a los objetos bounds que contenga. Por ejemplo, un
! bound puede ser un dihedro, que tiene 4 atomos, y abrá un boundgr que contenga
! el dihedro y que tendra la informacion de como usarlo en el computo de energía

use gems_constants
use gems_program_types
use gems_algebra
use gems_neighbor
use gems_output

implicit none

private

public   :: write_ebend, write_estretch, write_etors

public   :: boundgrs,read_psf,read_prm,boundgr_l

! Estas variables son para escribir
! logical,private   :: noadv=.true.
! real(dp),private  :: var

! Por cada juego de parametros leido debe haber al menos 1 enlace o 1 angulo,
! etc. Entonces conviene construir como objetos los parametros y que estos
! posean adentro una lista de enlaces o angulos que cumplan con estos
! parametros. Entonces los parametros son los equivalentes a los atom groups o
! mas bien a los intergrups pero no son estos porque no requieren lista de
! vecinos.

type :: bound
  ! Group of atoms considered togheter like in a bond, angle, etc.
  type(atom_ap),allocatable  :: a(:)
end type bound

#define _NODE bound_dcl
#define _CLASS class(bound)
#include "cdlist_header.inc"

type :: boundgr
  ! Group of bounds with the same chemical species and therefore interaction
  ! parameters.

  ! Identidades atomicas
  integer,allocatable     :: z(:)

   ! Lista de objetos bound of size size(z)
  type(bound_dcl),pointer :: list
  integer                 :: n=0

  ! Counter for pointers to the list. This is used in case that two different
  ! interact procedures are needed for the same bounds with the same chemical
  ! species.
  integer                 :: links=0

  ! Parametros reales del potencial
  type(real_v) :: p

  ! Parametros enteros del potencial
  type(integer_v) :: i

  ! Energia Potencial
  real(dp)                :: epot=0._dp

  ! Tipo de interaccion
  procedure(boundgr0),pointer :: interact=>null()  ! funcion de interaccion

  contains
    procedure   :: attach => boundgr_add
    procedure   :: init => boundgr_init
    procedure   :: destroy => boundgr_destroy
    ! final       :: boundgr_destroy
end type

abstract interface
  subroutine boundgr0(this)
   import boundgr
   class(boundgr),intent(inout)  :: this
  end subroutine
end interface

#define _NODE boundgr_l
#define _CLASS type(boundgr)
#define _ALLOCATABLE
#include "list_header.inc"

type(boundgr_l),target  :: boundgrs

contains

#define _NODE bound_dcl
#define _CLASS class(bound)
#include "cdlist_body.inc"

#define _NODE boundgr_l
#define _CLASS type(boundgr)
#include "list_body.inc"

subroutine charm_stretching(bg)
  !The functional form can be found here:
  !  MacKerell, A. D., Banavali, N., & Foloppe, N. (2000). Development and current
  !  status of the CHARMM force field for nucleic acids. Biopolymers, 56(4),257–265.
  !This also works for OPLS, see:
  !  Jorgensen, W. L., Maxwell, D. S., & Tirado-Rives, J. (1996). Development
  !  and testing of the OPLS all-atom force field on conformational energetics and
  !  properties of organic liquids. J. Am. Chem. Soc., 118, 11225–11236.
  use gems_inq_properties, only:vdistance
  class(boundgr),intent(inout)  :: bg
  type(atom),pointer            :: i,j
  type(bound_dcl),pointer       :: ln
  integer                       :: m
  real(dp)                      :: factor(dm),pstret,fstret,difr
  real(dp)                      :: ke,de
  real(dp)                      :: vd(dm),d

  bg%epot=0._dp
  ke=bg%p%o(1)
  de=bg%p%o(2)

  ln => bg%list
  do m = 1,bg%n
    ln=>ln%next

    ! Los atomos
    i => ln%o%a(1)%o
    j => ln%o%a(2)%o

    ! Distancias
    call vdistance(vd,i,j, mic)
    d  = dsqrt(dot_product(vd,vd))

    difr=d-de
    fstret=-2_dp*ke*difr
    pstret=ke*difr*difr

    i%epot=i%epot+pstret*0.5_dp
    j%epot=j%epot+pstret*0.5_dp
    bg%epot = bg%epot + pstret

    factor(1:dm) = fstret*vd(1:dm)/d
    i%force(1:dm) = i%force(1:dm) + factor(1:dm)
    j%force(1:dm) = j%force(1:dm) - factor(1:dm)
 
  end do

end subroutine charm_stretching

subroutine charm_bending(bg)
  !The functional form can be found here:
  !  MacKerell, A. D., Banavali, N., & Foloppe, N. (2000). Development and current
  !  status of the CHARMM force field for nucleic acids. Biopolymers, 56(4),257–265.
  !
  !This also works for OPLS, see:
  !  Jorgensen, W. L., Maxwell, D. S., & Tirado-Rives, J. (1996). Development
  !  and testing of the OPLS all-atom force field on conformational energetics and
  !  properties of organic liquids. J. Am. Chem. Soc., 118, 11225--11236.
  !
  use gems_inq_properties, only:vdistance
  use gems_elements
  class(boundgr),intent(inout)  :: bg
  type(atom),pointer            :: i,j,k
  type(bound_dcl),pointer       :: ln
  real(dp)                      :: dcosthda(dm),dcosthdb(dm)
  real(dp)                      :: tita,pbending,fbending
  real(dp)                      :: diftita,costh,sinth
  integer                       :: m
  real(dp)                      :: ke,te
  real(dp)                      :: a(3),am,b(3),bm

  bg%epot=0._dp
  ke=bg%p%o(1)
  te=bg%p%o(2)

  ln => bg%list
  do m = 1,bg%n
    ln=>ln%next

    ! los atomos
    i => ln%o%a(1)%o
    j => ln%o%a(2)%o
    k => ln%o%a(3)%o

    ! Distancia ij
    call vdistance(a,i,j, mic)
    am =  dsqrt(dot_product(a,a))

    ! Distancia jk
    call vdistance(b,k,j, mic)
    bm =  dsqrt(dot_product(b,b))

    ! Angulo
    costh=dot_product(a,b)/(am*bm)
    costh=min(costh,1.0_dp)
    costh=max(costh,-1.0_dp)
    tita=acos(costh)

    ! La referencia de arriba tiene la ecuacion esta multiplicada por 2
    diftita=tita-te
    pbending=ke*diftita*diftita

    i%epot = i%epot+pbending/3.0_dp
    j%epot = j%epot+pbending/3.0_dp
    k%epot = k%epot+pbending/3.0_dp
    bg%epot = bg%epot + pbending

    ! sinth = sin(theta)
    sinth = sqrt(1._dp-costh*costh)

    ! pi angle

    if (sinth<1.e-6_dp) then
 
      ! FIXME
      ! Acer lo mismo que en torsion (es decir usar las dos formulas)
      print *, 'FIXME! pi angle in bending'

      ! Cuando theta0 es 0 o pi, entonces (theta-theta0)/sin(theta) es un
      ! limite notable y tiende a 1, por lo que se puede 
      if (diftita < 0._dp) then
        fbending=2*ke
      else
        fbending= -2*ke
      endif
      
      ! TODO: Deducir el angulo con el producto vectorial? el angulo entre los
      ! dos vectores sería deducido como un seno y la derivada sería del
      ! coseno. Algo parecido a lo que se hizo en torsion
    else
 
      ! Using dcos(theta)/di=-sin(theta) dtheta/di we have:
      !
      !   dtheta/di=-1/sin(theta)*dcos(theta)/di  ; sin(theta)/=0
      !

      fbending=-2*ke*diftita/sinth

      ! Using the chain rule for theta(A(i)):
      !
      !   dcos(theta)/di= (dA/di)^T dcos(theta)/dA
      !
      ! Where dA/di is the Jacobian J_A(i) so it is important to keep the order of
      ! the terms.
      !
      ! Then, using costh=a.b/|a||b| (see Blondel, A., & Karplus, M. (1996).)
      !
      ! dcosthda = b/|a||b| + a.b/|b| d(1/|a|)/da
      !            = b/|a||b| + a.b/|b| (-1/|a|^2) d(|a|)/da
      !            = b/|a||b| + a.b/|b| (-1/|a|^2) a/|a|
      !            = b/|a||b| - costh /|a|^2
      !
      ! And da/i=I , the identity

      dcosthda=b/(am*bm)-costh*a/(am*am)
      dcosthdb=a/(am*bm)-costh*b/(bm*bm)

      dcosthda=fbending*dcosthda
      dcosthdb=fbending*dcosthdb

    endif

    i%force(1:dm) = i%force(1:dm) -dcosthda
    j%force(1:dm) = j%force(1:dm) +dcosthda+dcosthdb
    k%force(1:dm) = k%force(1:dm) -dcosthdb

    ! print *, 'ASD', m-1, elements%o(bg%z(1))%sym, elements%o(bg%z(2))%sym, elements%o(bg%z(3))%sym
    ! print *, bg%p%o(1)*ui_kcm,bg%p%o(2)*rad_grad
  end do

end subroutine charm_bending

subroutine charm_torsion(bg)
  use gems_elements
!
!The functional form can be found here:
!  MacKerell, A. D., Banavali, N., & Foloppe, N. (2000). Development and current
!  status of the CHARMM force field for nucleic acids. Biopolymers, 56(4),257--265.
!
!The tecnical details of the implementation are here:
!  Blondel, A., & Karplus, M. (1996). New formulation for derivatives of torsion
!  angles and improper torsion angles in molecular mechanics: Elimination of
!  singularities. Journal of Computational Chemistry, 17(9), 1132--1141.
!
! La forma funcional del charm es:
!
!   ke*cos(n*phi-delta)
!
! Para ciertos valores de delta, esta formula se puede convertir usando la formula del angulo multiple
! (multiple-angle formula):
!
! cos(nx)=n!/((n-k)!k!)cos^k(x)sin^(n-k)(x)cos((n-k)pi/2)
!
! Si bien la mayoria de los delta son 0 o 180 en principio este puede tener cualquier
! valor. Por eso, charm_torsion y opls_torsion son diferentes.
!
class(boundgr),intent(inout)   :: bg
type(bound_dcl),pointer        :: ln
type(atom),pointer             :: i,j,k,l
integer                        :: m,n,t
real(dp)                       :: dphidf(dm),dphidg(dm),dphidh(dm),dcosda(dm),dcosdb(dm),dsindb(dm),dsindv(dm)
real(dp)                       :: cphi,sphi,ptors,ftors,phi
real(dp)                       :: f(dm),g(dm),h(dm),a(dm),b(dm),v(dm)
real(dp)                       :: ke,delta,invm_a,invm_b,invm_v

bg%epot = 0.0_dp

ln => bg%list
do m = 1,bg%n
  ln=>ln%next

  ! los atomos
  i => ln%o%a(1)%o
  j => ln%o%a(2)%o
  k => ln%o%a(3)%o
  l => ln%o%a(4)%o

  ! Distancias y nomenclatura siguiendo Blondel, A., & Karplus, M. (1996).
  call vdistance(f,i,j, mic)
  call vdistance(g,j,k, mic)
  call vdistance(h,k,l, mic)  ! With a minus

  ! Productos cruz.
  a=cross_product(f,g)
  b=cross_product(g,h)
  v=cross_product(g,a)

  ! Modulos
  invm_a=1._dp/sqrt(dot_product(a,a))
  invm_b=1._dp/sqrt(dot_product(b,b))
  invm_v=1._dp/sqrt(dot_product(v,v))

  ! The vectors A,G,and V=AxG form an orthogonal frame. Vector a is in the
  ! plane (V,G) of this orthogonal frame. The projection of H in the plane (V,A)
  ! are the coordinates of H in this orthogonal frame, because H origin is in G.
  ! So the angle between (V,G) plane and H coordinates in the plane (V,A) is the
  ! angle phi. Coordinates of H in V and A axis are the cos and sin of phi
  ! respectively. Phi can be found also from G since is just h rotated 90
  ! degrees in the (V,G,A) frame.
  ! Proyecciones de G me dan sphi y cphi. Me pregunto porque no proyecto
  ! directamente H.
  cphi = dot_product(b,a)*(invm_a*invm_b);
  sphi = dot_product(b,v)*(invm_v*invm_b); ! With a minus
  phi = -atan2(sphi,cphi)

  ! Multipe dihedral terms
  ptors=0._dp
  ftors=0._dp
  do t=1,bg%i%o(1)
    ke    = bg%p%o(2*(t-1)+1)
    delta = bg%p%o(2*(t-1)+2)
    n     = bg%i%o(t+1)

    ! if(n==0) then
    !   ! TODO: reference of this
    !   diff=(phi-delta)
    !   if (diff<-pi) then
    !     diff = diff + twopi
    !   else if (diff>pi) then
    !     diff = diff - twopi
    !   endif
    !
    !   ptors=ptors+ke*diff*diff
    !   ftors=ftors+2*ke*diff
    !
    !   print *, 'FIXME: Ver charm_torsion'
    ! else
      ptors=ptors+ke*(1+cos(n*phi-delta))
      ftors=ftors-ke*n*sin(n*phi-delta)
    ! endif
  enddo

  i%epot = i%epot + ptors/4
  j%epot = j%epot + ptors/4
  k%epot = k%epot + ptors/4
  l%epot = l%epot + ptors/4
  bg%epot = bg%epot + ptors

  ! Normalizo b
  b=b*invm_b

  ! So, I need the dphi/dr1 (which is a 3d vector):
  if (abs(sphi)>0.1) then

    ! Using dcos(phi)/di=-sin(phi) dphi/di we have:
    !
    !   dphi/di=-1/sin(phi)*dcos(phi)/di  ; sin(phi)/=0
    !
    ! Using the chain rule for phi(F(i)):
    !
    !   dcos(phi)/di= (dF/di)^T dcos(phi)/dF
    !
    ! Where dF/di is the Jacobian J_F(i) so it is important to keep the order of
    ! the terms. Making phi(A(F)) we have
    !
    !   dcos(phi)/dF= (dA/dF)^T dcos(phi)/dA  ; Where dA/dF is the Jacobian J_A(F)
    !
    ! and using this in the precious equation to find the chain rule for phi(A(F(i)))
    !
    !   dcos(phi)/di= (dF/di)^T (dA/dF)^T dcos(phi)/dA
    !
    ! Then dF/di=I with I the identity matrix and, since A=FxG, dA/dF can be
    ! wrote as the operator IxG as (IxG)P=PxG. The transpose (dA/dF)^T gives the
    ! operator -IxG since it is antisymmetric, so:
    !
    !   dcos(phi)/di= (dA/dF)^T dcos(phi)/dA
    !               = -dcos(phi)/dA x G

    ! Normalizo a
    a=a*invm_a

    ! Eq 13 of Blondel, A., & Karplus, M. (1996) and symmetric equivalent
    dcosda = invm_a*(cphi*a-b)  ! With a minus
    dcosdb = invm_b*(cphi*b-a)  ! With a minus

    ! From Eq 9 of Blondel, A., & Karplus, M. (1996).
    ftors = ftors/sphi

    ! Minus Eq 20 of Blondel, A., & Karplus, M. (1996)
    dphidf = ftors*cross_product(g,dcosda) ! With a minus

    ! Symmetric equivalente of Eq 20 of Blondel, A., & Karplus, M. (1996)
    dphidh = ftors*cross_product(dcosdb,g) ! With a minus

    ! Eq 23 of Blondel, A., & Karplus, M. (1996)
    ! Note that here dB/dG=IxH so (dB/dG)^T dphi/dB = -ftors*(-dcosdb)xH
    dphidg = ftors*(cross_product(dcosda,f)+cross_product(h,dcosdb)) ! With a minus

  else

    ! Improtante, notar que aca dvecotr/dvector denota jacobiano mientras que
    ! dscalar/dvector denota gradiente (vector columna), es decir, la notacion
    ! no es consistente (es mejor trabajar con vectores filas, como resulta de
    ! usar jacobianos).
    !
    ! Using dsin(phi)/di=cos(phi) dphi/di  we have:
    !
    !   dphi/di=1/cos(phi)*dsin(phi)/di  ; cos(phi)/=0
    !
    ! Then using sin(phi)=B.V/(|V||B|)
    !
    !   dsin(phi)/dF= (dA/dF)^T (dV/dA)^T dsin(phi)/dV ; (dA/dF)^T=-IxG ; (dV/dA)^T=-GxI
    !               = (G x dsindv) x G
    !
    !   dsin(phi)/dG= (dV/dG)^T dsin(phi)/dV + (dB/dG)^T dsin(phi)/dB
    !               = ( IxA +(GxI)(FxI) )^T dsin(phi)/dV + (IxH)^T dsin(phi)/dB
    !               = (-IxA +(FxI)(GxI) ) dsin(phi)/dV - (IxH) dsin(phi)/dB
    !               = (-dsindv x A + Fx(Gx dsindv) - dsindb x H 
    !               = -dsindv x (FxG) +Fx(G x dsindv) - dsindb x H 
    !
    !   dsin(phi)/dH= (dB/dH)^T dsin(phi)/dB
    !

    ! Normalizo V=GxA
    v=v*invm_v

    ! As Eq 13 of Blondel, A., & Karplus, M. (1996) from sin(phi)=B.V/(|V||B|)
    dsindv = invm_v*(sphi*v-b) ! With a minus
    dsindb = invm_b*(sphi*b-v) ! With a minus

    ftors = -ftors/cphi  ! With a minus

    ! TODO: If there is a more efficient way, declare a 3 cross product function
    dphidf = cross_product(g,dsindv)
    dphidf = cross_product(dphidf,g)
    dphidf = ftors*dphidf ! With a minus from ftors

    ! (dA/dG)^T=-FxI ; (dV/dG)^T=-IxA ; (dB/dG)^T=-IxH 
    dphidg = cross_product(g,dsindv)
    dphidg = cross_product(f,dphidg)
    dphidg = dphidg+cross_product(a,dsindv)
    dphidg = dphidg+cross_product(h,dsindb)
    dphidg = ftors*dphidg ! With a minus from ftors

    ! (dB/dH)^T=-GxI
    dphidh = ftors*cross_product(dsindb,g); ! With a minus from ftors

  endif
 
  ! Using eq 26 of Blondel, A., & Karplus, M. (1996)
  i%force = i%force + dphidf
  j%force = j%force + dphidg - dphidf
  k%force = k%force + dphidh - dphidg
  l%force = l%force - dphidh
             
end do

! print *, 'ASD', m-1, elements%o(bg%z(1))%sym, elements%o(bg%z(2))%sym, elements%o(bg%z(3))%sym, elements%o(bg%z(4))%sym
! print *, bg%p%o(1)*ui_kcm,bg%i%o(1),bg%i%o(2),bg%p%o(2)*rad_grad

end subroutine charm_torsion

subroutine opls_torsion(bg)
! OPLS
!Jorgensen, W. L., Maxwell, D. S., & Tirado-rives, J. (1996). Development and
!Testing of the OPLS All-Atom Force Field on Conformational Energetics and
!Properties of Organic Liquids, 118(15), 11225-11236.
class(boundgr),intent(inout)   :: bg
type(bound_dcl),pointer        :: ln
type(atom),pointer             :: i,j,k,l
integer                        :: m
real(dp)                       :: dud1(dm),dud2(dm),dud3(dm),dud4(dm)
real(dp)                       :: dnum(dm),dden(dm)
real(dp)                       :: bxc(3),ab,bc,aa,bb,cc,axb_mod,bxc_mod
real(dp)                       :: ac,cosphi,cphi,cphi2,ptors,ftors,de1,den
real(dp)                       :: a(dm),b(dm),c(dm)

bg%epot = 0.0_dp

ln => bg%list
do m = 1,bg%n
  ln=>ln%next

  ! los atomos
  i => ln%o%a(1)%o
  j => ln%o%a(2)%o
  k => ln%o%a(3)%o
  l => ln%o%a(4)%o

  ! Distancias y nomenclatura
  call vdistance(a,j,i, mic)
  call vdistance(b,k,j, mic)
  call vdistance(c,l,k, mic)

  ! Producto cruz b con c
  bxc=cross_product(b,c)
  bxc(1)=bxc(1)-box(1)*idnint(bxc(1)/box(1))
  bxc(2)=bxc(2)-box(2)*idnint(bxc(2)/box(2))
  bxc(3)=bxc(3)-box(3)*idnint(bxc(3)/box(3))

  ! Productos escalares cruzados
  ab=dot_product(a,b)
  bc=dot_product(b,c)
  ac=dot_product(a,c)

  ! Modulo cuadrado de cada distancia
  aa=dot_product(a,a)
  bb=dot_product(b,b)
  cc=dot_product(c,c)

  ! Productors vectoriales al cuadrado entre distancias
  axb_mod=aa*bb-ab*ab
  bxc_mod=bb*cc-bc*bc

  ! Denominador para calcular el coseno de phi
  den=axb_mod*bxc_mod

  !chequeo de que átomos no están en la misma línea
  if(den<1.e-12_dp) cycle

  den=dsqrt(den) ! El producto de los vectoriales

  ! Coseno del ángulo
  ! Se puede calcular a partir del producto escalar de los vectores normales a
  ! los planos que forman el angulo dihedro. Es decir:
  !
  !   cosphi=-dot_product(axb,bxc)/(axb_mod*bxc_mod)
  !
  ! Ahora el numerador se puede calcular facilmente usando la identidad
  !
  !   dot_product(axb,cxd)=(b^T(c^Ta)I- ca^T)d=ac*bd-ad*bc
  !
  ! Donde ^T indica transpuesta.
  !
  cosphi=(ab*bc-ac*bb)/den
  cosphi=min(cosphi,1.0_dp)
  cosphi=max(cosphi,-1.0_dp)

  ! Ryckaert−Bellemans (RB) potential function
  cphi=-cosphi  ! cos(phi-pi)
  cphi2=cphi*cphi

  ptors= bg%p%o(0)                   &
        +bg%p%o(1)*cphi              &
        +bg%p%o(2)*cphi2             &
        +bg%p%o(3)*cphi2*cphi        &
        +bg%p%o(4)*cphi2*cphi2       &
        +bg%p%o(5)*cphi2*cphi2*cphi
  ptors=ptors*0.25_dp

  ftors= bg%p%o(1)                           &
        +bg%p%o(2)*2.0_dp*cphi               &
        +bg%p%o(3)*3.0_dp*cphi2              &
        +bg%p%o(4)*4.0_dp*cphi2*cphi         &
        +bg%p%o(5)*5.0_dp*cphi2*cphi2

  ! ! debug para testear topologias. Fija una determinada torsion.
  ! if(m==tfix) then
  !   ptors=0.1_dp*(dabs(phi)-dabs(phifix))**2
  !   ftors=-0.2_dp*dabs((dabs(phi)-dabs(phifix)))
  ! endif

  i%epot = i%epot + ptors
  j%epot = j%epot + ptors
  k%epot = k%epot + ptors
  l%epot = l%epot + ptors
  bg%epot = bg%epot + 4*ptors

  de1=ftors/den
  axb_mod=axb_mod/den*cosphi
  bxc_mod=bxc_mod/den*cosphi

  dnum=c*bb-b*bc
  dden=(ab*b-a*bb)*bxc_mod
  dud1=(dnum-dden)*de1

  dnum=((b-a)*bc-ab*c)+(2.0_dp*ac*b-c*bb)
  dden=axb_mod*(bc*c-b*cc)+(a*bb-aa*b-ab*(b-a))*bxc_mod
  dud2=(dnum-dden)*de1

  dnum=ab*b-a*bb
  dden=axb_mod*(bb*c-bc*b)
  dud4=(dnum-dden)*de1

  dud3=-(dud1+dud2+dud4)

  i%force = i%force + dud1
  j%force = j%force + dud2
  k%force = k%force + dud3
  l%force = l%force + dud4


end do

end subroutine opls_torsion

!  subroutine impropers(f,epot)
!   ! OPLS improper
!   use gems_inq_properties, only:vdistance
!   use gems_input_parsing, only:ioprefix
!   class(ffunction),intent(in) :: f
!   type(ffnode),pointer        :: fn
!   type(atom),pointer          :: i,j,k,l
!   integer                     :: m
!   real(dp),intent(out)        :: epot
!   real(dp)                   :: dud1(dm),dud2(dm),dud3(dm),dud4(dm)
!   real(dp)                   :: vij(dm),vjk(dm),vkl(dm),dnum(dm),dden(dm)
!   real(dp)                   :: bxc(3)
!   real(dp)                   :: ab,bc,at,bt,ct,axb_mod,bxc_mod,rnum
!   real(dp)                   :: ac,cosphi,cphi,cphi2,sphi,signo,ptors,ftors,de1,den,phi
!   real(dp)                   :: k0,n0,phis
!
!   epot = 0.0_dp
!
!   fn => f%head
!   do m=1,f%nmatch
!     fn => fn%next
!
!     ! los atomos
!     i => fn%a
!     j => fn%b
!     k => fn%c
!     l => fn%d
!
!     ! los parametros
!     phis = f%param(fn%id,1)
!     k0 = f%param(fn%id,2)
!     n0 = f%param(fn%id,3)
!
!     ! Distancias
!     vij=vdistance(j,i, mic)
!     vjk=vdistance(k,j, mic)
!     vkl=vdistance(l,k, mic)
!
!     ! Producto cruz b con c
!     bxc(1)=vjk(2)*vkl(3)-vjk(3)*vkl(2)
!     bxc(2)=vjk(3)*vkl(1)-vjk(1)*vkl(3)
!     bxc(3)=vjk(1)*vkl(2)-vjk(2)*vkl(1)
!     bxc(:)=bxc(:)-box(:)*idnint(bxc(:)*one_box(:))
!
!     ! Productos escalares cruzados
!     ab=dot(vij,vjk)
!     bc=dot(vjk,vkl)
!     ac=dot(vij,vkl)
!
!     ! Modulo cuadrado de cada distancia
!     at=dot(vij,vij)
!     bt=dot(vjk,vjk)
!     ct=dot(vkl,vkl)
!
!     ! Productors vectoriales al cuadrado entre distancias
!     axb_mod=at*bt-ab*ab
!     bxc_mod=bt*ct-bc*bc
!
!     !
!     rnum=ab*bc-ac*bt ! El producto cuadrado de los vectoriales?
!     den=axb_mod*bxc_mod      ! El producto cuadrado de los vectoriales
!
!     !chequeo de que átomos no están en la misma línea
!
!     if(den<1.d-12) cycle
!
!     den=dsqrt(den) ! El producto de los vectoriales
!
!     !coseno del ángulo
!     cosphi=rnum/den
!     cosphi=dmin1(cosphi,1.0_dp)
!     cosphi=dmax1(cosphi,-1.0_dp)
!
!     !signo del ángulo
!     signo=dot_product(vij,bxc)
!
!     !definición del ángulo
!     phi=dsign(dacos(cosphi),signo)
!
!     !if (impropia) phi = n*phi -phi_s
!     !phis(i)=phi
!
!     ! Ryckaert−Bellemans (RB) potential function
!     cphi=dcos(phi-pi)
!     sphi=dsin(phi-pi)
!     cphi2=cphi*cphi
!
!     ptors=k0*(1+dcos(n0*phi-phis))
!
!     ftors=n0*k0*dsin(n0*phi-phis)
!
!
!     ! debug para testear topologias. Fija una determinada torsion.
!     if(m==tfix) then
!       ptors=0.1_dp*(dabs(phi)-dabs(phifix))**2
!       ftors=-0.2_dp*dabs((dabs(phi)-dabs(phifix)))
!     endif
!
!     i%epot = i%epot + ptors*ev_ui*0.25_dp
!     j%epot = j%epot + ptors*ev_ui*0.25_dp
!     k%epot = k%epot + ptors*ev_ui*0.25_dp
!     l%epot = l%epot + ptors*ev_ui*0.25_dp
!
!     epot = epot + ptors*ev_ui
!
!     if(dabs(sphi)<1.d-12) sphi=dsign(1.d-12,sphi)
!     de1=ftors/den/sphi
!     axb_mod=axb_mod/den*cosphi
!     bxc_mod=bxc_mod/den*cosphi
!
!     dnum=vkl*bt-vjk*bc
!     dden=(ab*vjk-vij*bt)*bxc_mod
!     dud1=(dnum-dden)*de1
!
!     dnum=((vjk-vij)*bc-ab*vkl)+(2.0_dp*ac*vjk-vkl*bt)
!     dden=axb_mod*(bc*vkl-vjk*ct)+(vij*bt-at*vjk-ab*(vjk-vij))*bxc_mod
!     dud2=(dnum-dden)*de1
!
!     dnum=ab*vjk-vij*bt
!     dden=axb_mod*(bt*vkl-bc*vjk)
!     dud4=(dnum-dden)*de1
!
!     dud3=-(dud1+dud2+dud4)
!
!
!     i%force = i%force + dud1*ev_ui
!     j%force = j%force + dud2*ev_ui
!     k%force = k%force + dud3*ev_ui
!     l%force = l%force + dud4*ev_ui
!
!
!  end do
!
!  end subroutine impropers
!
!  ! Si los parametros estan dados para la forma de Fourier, cambio a los
!  ! parametros adecuados para poder calcularla con la funcion RB que
!  ! esta programada y es mas eficiente. Para esta conversion consultar:
!  !
!  !  Apol, E.; Apostolov, R.; Berendsen, H. J. C.; van Buuren, A.; Bjelkmar, P.;
!  !  van Drunen, R.; Feenstra, A.; Groenhof, G.; Kasson, P.; Larsson, P.;
!  !  Meulenhoff, P.; Murtola, T.; Pall, S.; Pronk, S.; Schulz, R.; Shirts, M.;
!  !  Sijbers, A.; Tieleman, P.; Hess, B.; van der Spoel, D.; Lindahl, E. GROMACS
!  !  User Manual, version 4.5.4; Royal Institute of Technology and Uppsala
!  !  University: Stockholm and Uppsala, Sweden, 2010. Pagina 76
!
!  f1=this%f(i)%param(j,1)
!  f2=this%f(i)%param(j,2)
!  f3=this%f(i)%param(j,3)
!  f4=this%f(i)%param(j,4)
!
!  this%f(i)%param(j,1)=f2+0.5_dp*(f1+f3)
!  this%f(i)%param(j,2)=0.5_dp*(-f1+3.0_dp*f3)
!  this%f(i)%param(j,3)=-f2+4.0_dp*f4
!  this%f(i)%param(j,4)=-2.0_dp*f3
!  !this%f(i)%param(j,5)=-4.0_dp*f4
!  !this%f(i)%param(j,6)=-0.0_dp

subroutine read_prm(prmfile)
  use gems_input_parsing
  ! use gems_elements, only:inq_z, set_z
  use gems_elements
  use gems_errors, only:wlog
  use gems_algebra, only:sort_int
  character(*),intent(in)  :: prmfile
  type(boundgr_l),pointer   :: ln,lp
  type(input_options), target   :: iopts
  character(:),allocatable :: clase,w1,w2,w3,w4
  character(10)            :: zx(4)
  real(dp)                 :: f1,f2
  integer                  :: ix(4),u,i1,i2,i

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
    case('ATOMS','BONDS','ANGLES','DIHEDRALS','IMPROPER','CMAP')
      clase=w1
      cycle
    endselect

    selectcase(clase)
    case('ATOMS')

      if(adjustl(w1)=='MASS') then
        call readi(ix(1))
        call reada(w1) ! nomb
        call readf(f2) ! mass
        i1=inq_z(w1)
        if(i1/=0) call set_z(i1,mass=f2)
      endif

    case('BONDS')

      !V(bond) = Kb(b - b0)**2
      !Kb: kcal/mole/A**2
      !b0: A

      call reada(w2) !nomb2
      call readf(f1); f1=f1*kcm_ui
      call readf(f2)

      ! Solving symmetry in bonds
      ix(1)=inq_z(w1)
      ix(2)=inq_z(w2)
      call sort_int(ix(1:2))

      ! Add parameters to the corresponding bondgr
      ln => boundgrs
      do while(associated(ln%next))
        ln=>ln%next
        if (size(ln%o%z)/=2) cycle
        if (any(ix(1:2)/=ln%o%z(:))) cycle

        ! TODO: If manyterms is true here then it is needed to replicate the
        ! boundgr and assoicate it with charm. This should be needed if there is
        ! the same chemical species in two force fields, for example in OPLS and
        ! Charmm, and both are actives.
        if (.not.associated(ln%o%interact)) ln%o%interact=>charm_stretching

        call ln%o%p%put(1,f1)
        call ln%o%p%put(2,f2)
        exit
      enddo

    case('ANGLES')

      !V(angle) = Ktheta(Theta - Theta0)**2
      !V(Urey-Bradley) = Kub(S - S0)**2
      !Ktheta: kcal/mole/rad**2
      !Theta0: degrees
      !Kub: kcal/mole/A**2 (Urey-Bradley)
      !S0: A

      call reada(w2) !nomb2
      call reada(w3) !nomb3
      call readf(f1); f1=f1*kcm_ui
      call readf(f2); f2=f2*grad_rad

      ! Sort before compare
      ix(1)=inq_z(w1)
      ix(2)=inq_z(w3)
      call sort_int(ix(1:2))
      ix(3)=ix(2)
      ix(2)=inq_z(w2)

      ! Add parameters to the corresponding anglegr
      ln => boundgrs
      do while(associated(ln%next))
        ln=>ln%next
        if (size(ln%o%z)/=3) cycle
        if (any(ix(1:3)/=ln%o%z(:))) cycle

        ! TODO: If manyterms is true here then it is needed to replicate the
        ! boundgr and assoicate it with charm. This should be needed if there is
        ! the same chemical species in two force fields, for example in OPLS and
        ! Charmm, and both are actives.
        if (.not.associated(ln%o%interact)) ln%o%interact=>charm_bending

        call ln%o%p%put(1,f1)
        call ln%o%p%put(2,f2)
        exit
      enddo

    case('DIHEDRALS')

			!V(dihedral) = Kchi(1 + cos(n(chi) - delta))
			!Kchi: kcal/mole
			!n: multiplicity
			!delta: degrees

      call reada(w2) !nomb2
      call reada(w3) !nomb3
      call reada(w4) !nomb4
      call readf(f1); f1=f1*kcm_ui
      call readi(i1) ! Multiplicity
      call readf(f2); f2=f2*grad_rad

      ! Trough isoforms
      do i=1,2

        select case(i)
        case(1)
          ix(1)=inq_z(w1)
          ix(2)=inq_z(w2)
          ix(3)=inq_z(w3)
          ix(4)=inq_z(w4)
          zx(1)=trim(adjustl(w1))
          zx(2)=trim(adjustl(w2))
          zx(3)=trim(adjustl(w3))
          zx(4)=trim(adjustl(w4))
        case(2)
          ix(1)=inq_z(w4)
          ix(2)=inq_z(w3)
          ix(3)=inq_z(w2)
          ix(4)=inq_z(w1)
          zx(1)=trim(adjustl(w4))
          zx(2)=trim(adjustl(w3))
          zx(3)=trim(adjustl(w2))
          zx(4)=trim(adjustl(w1))
        end select
         
        ln => boundgrs
        do while(associated(ln%next))
          lp=>ln
          ln=>ln%next

          if (size(ln%o%z)/=4) cycle
          if (any(ix(1:4)/=ln%o%z(1:4))) cycle

          ! TODO: If manyterms is true here then it is needed to replicate the
          ! boundgr and assoicate it with charm. This should be needed if there is
          ! the same chemical species in two force fields, for example in OPLS and
          ! Charmm, and both are actives.
          if (.not.associated(ln%o%interact)) ln%o%interact=>charm_torsion

          if (f1==0._dp) then
            call wlog ('-PRM->'); write(logunit,'(a,4(a5))') 'No dihedral for:  ', zx(1:4)
            lp%next=>ln%next
            deallocate(ln)
            exit
          endif

          ! Adding a new term (for dihedral multiple entrees)
          i2=ln%o%i%o(1)
          call ln%o%p%put(2*i2+1,f1)
          call ln%o%p%put(2*i2+2,f2)
          i2=i2+1
          call ln%o%i%put(1,i2)
          call ln%o%i%put(i2+1,i1)
                     
          exit
        enddo
      enddo

    case('IMPROPER')
    case('CMAP')
    endselect

  enddo

  ! Read again the file in search for dihedrals without parameters that may use the wildcard atom type X
  rewind(u)
  do
    call read_line(eof)
    if (eof) exit

    call reada(w1)
    if (w1(1:1)=='!') cycle

    ! Selecciono la clase
    selectcase(w1)
    case('ATOMS','BONDS','ANGLES','DIHEDRALS','IMPROPER','CMAP')
      clase=w1
      cycle
    endselect

    selectcase(clase)
    case('DIHEDRALS')
      call reada(w2) !nomb2
      call reada(w3) !nomb3
      call reada(w4) !nomb4
      call readf(f1); f1=f1*kcm_ui
      call readi(i1) ! Multiplicity
      call readf(f2); f2=f2*grad_rad
     
      ! Trough isoforms
      do i=1,2

        select case(i)
        case(1)
          ix(1)=inq_z(w1)
          ix(2)=inq_z(w2)
          ix(3)=inq_z(w3)
          ix(4)=inq_z(w4)
          zx(1)=trim(adjustl(w1))
          zx(2)=trim(adjustl(w2))
          zx(3)=trim(adjustl(w3))
          zx(4)=trim(adjustl(w4))
        case(2)
          ix(1)=inq_z(w4)
          ix(2)=inq_z(w3)
          ix(3)=inq_z(w2)
          ix(4)=inq_z(w1)
          zx(1)=trim(adjustl(w4))
          zx(2)=trim(adjustl(w3))
          zx(3)=trim(adjustl(w2))
          zx(4)=trim(adjustl(w1))
        end select

        ln => boundgrs
        do while(associated(ln%next))
          lp=>ln
          ln=>ln%next

          if (size(ln%o%z)/=4) cycle
          if (associated(ln%o%interact)) cycle
          if (ix(1)/=ln%o%z(1) .and. zx(1)/='X') cycle
          if (ix(2)/=ln%o%z(2) .and. zx(2)/='X') cycle
          if (ix(3)/=ln%o%z(3) .and. zx(3)/='X') cycle
          if (ix(4)/=ln%o%z(4) .and. zx(4)/='X') cycle

           call wlog ('-PRM->'); write(logunit,'(2(a,4(a5)))') &
             'Taking wildcard ', zx(1:4), &
             ' for dihedral ', elements%o(ln%o%z(1))%sym, elements%o(ln%o%z(2))%sym,&
                               elements%o(ln%o%z(3))%sym, elements%o(ln%o%z(4))%sym
          ln%o%interact => charm_torsion
          if (f1==0._dp) then
            call wlog ('-PRM->'); write(logunit,'(a,4(a5))') 'No dihedral for:  ', zx(1:4)
            lp%next=>ln%next
            deallocate(ln)
            exit
          endif

          ! Adding a new term (for dihedral multiple entrees)
          i2=ln%o%i%o(1)
          call ln%o%p%put(2*i2+1,f1)
          call ln%o%p%put(2*i2+2,f2)
          i2=i2+1
          call ln%o%i%put(1,i2)
          call ln%o%i%put(i2+1,i1)

        enddo

      enddo

    endselect

  enddo

  ! TODO: Chequeo que todos los atomos tengan al menos la masa por defecto
 
  ! Chequeo que todos los boundgrs tengan parametros. Asumo que si no los tienen
  ! se pueden borrar. Cuidado, cuando lea muchos archivos .prm esto no sirve
  ! TODO: Dar error si hay un dihedro sin parametro
  ln => boundgrs
  do while(associated(ln%next))
    lp => ln
    ln => ln%next
    if (.not.associated(ln%o%interact)) then
      call wlog ('-PRM->'); write(logunit,'(a,4(a5))') &
           'Ignoring dihedral types ', elements%o(ln%o%z(1))%sym, elements%o(ln%o%z(2))%sym,&
                                       elements%o(ln%o%z(3))%sym, elements%o(ln%o%z(4))%sym
      lp%next=>ln%next
      deallocate(ln)
      ln=>lp
    endif
  enddo
     

  ! Devuelvo el input parsing al gems
  opts=>gems_iopts
  close(u)
end subroutine read_prm

subroutine read_psf(topfile)
  use gems_elements, only: add_z, inq_z
  use gems_groups, only: atom_setelmnt
  character(*),intent(in)  :: topfile
  character(9)             :: clase
  integer                  :: i,j,k,l,u,n,natoms,m,ioflag,ix(4),iy(4)
  character(5)             :: seg,resn,nomb,tipo
  real(dp)                 :: carga, masa
  type(boundgr),pointer     :: bg

  u = find_io(30)
  open(u,action='read',file=topfile)

  read(u,fmt=*) ! PSF CMAP
  read(u,fmt=*)
  read(u,iostat=ioflag,fmt=*) n, clase

  ioflag=0
  do while(ioflag==0)
    ! if(ioflag/=0) exit

    selectcase(clase)
    case('!NTITLE')
      do m=1,n
        read(u,*)
      enddo

    case('!NATOM')

      ! Asumo que ya estan allocateado los atomos
      ! leyendo por ejemplo el pdb
      ! FIXME: sys index may have null atoms
      natoms=n
      do m=1,n
        read(u,*) i,seg,j,resn,nomb,tipo,carga,masa
        call add_z(tipo,masa,carga)
        call atom_setelmnt(sys%a(i)%o,inq_z(tipo))
      enddo

    case('!NBOND:')

      do m=1,n

        read(u,'(2(i8))',advance='no') ix(1:2)
        if(mod(m,4)==0) read(u,*)

        ! Solving symmetry in bonds. Need to do it in ix???
        iy(1)=sys%a(ix(1))%o%z
        iy(2)=sys%a(ix(2))%o%z
        call sort_int(iy(1:2))

        bg=>boundgrs_include(iy(1:2))
        call bg%attach(ix(1:2))

      enddo

    case('!NTHETA:')

      do m=1,n
        read(u,'(3(i8))',advance='no') ix(1:3)
        if(mod(m,3)==0) read(u,*)

        ! Solving symmetry in angles.
        iy(1)=sys%a(ix(1))%o%z
        iy(2)=sys%a(ix(3))%o%z
        call sort_int(iy(1:2))
        iy(3)=iy(2)
        iy(2)=sys%a(ix(2))%o%z

        bg=>boundgrs_include(iy(1:3))
        call bg%attach(ix(1:3))

      enddo

    case('!NPHI:')

      do m=1,n
        read(u,'(4(i8))',advance='no') ix(1:4)
        if(mod(m,2)==0) read(u,*)

        ! Solving symmetry in torsion.
        iy(1)=sys%a(ix(1))%o%z
        iy(2)=sys%a(ix(2))%o%z
        iy(3)=sys%a(ix(3))%o%z
        iy(4)=sys%a(ix(4))%o%z
        if(iy(1)>iy(4)) then
          iy(1)=iy(4)
          iy(2)=iy(3)
          iy(3)=sys%a(ix(2))%o%z
          iy(4)=sys%a(ix(1))%o%z
        endif

        bg=>boundgrs_include(iy(1:4))
        call bg%attach(ix(1:4))

      enddo
    case('!NIMPHI:')
      do m=1,n
        read(u,'(4(i8))',advance='no') i,j,k,l
        if(mod(m,2)==0) read(u,*)
      enddo
    case('!NDON:')
      ! FIXME
      read(u,*)
    case('!NACC:')
      ! FIXME
      read(u,*)
    case('!NNB:')
      ! FIXME
      read(u,*)
      do m=1,natoms
        read(u,'(i8)',advance='no') i
        if(mod(m,8)==0) read(u,*)
      enddo
    case('!NCRTERM:')
      do m=1,n
        read(u,'(8(i8))') i,j,k,l,i,j,k,l
      enddo
    endselect

    read(u,*)
    read(u,iostat=ioflag,fmt=*) n, clase

  enddo

  close(u)
end subroutine read_psf

function boundgrs_include(z) result(bg)
  use gems_algebra, only:sort_int
  ! Returns the boundgr in boundgrs that hold z. If there is not boundgr it
  ! creates one. Note that z is not sorted in this subrroutine, so it will
  ! consider a bond C-H different that a bond H-C.
  type(boundgr_l),pointer   :: ln
  type(boundgr),pointer     :: bg
  integer,intent(inout)     :: z(:)

  ! Compruebo si ya existe un boundgr para alojar este bound
  ln => boundgrs
  do while(associated(ln%next))
    ln=>ln%next

    if (size(ln%o%z)/=size(z)) cycle
    if (any(z(:)/=ln%o%z(:))) cycle

    bg=>ln%o
    return
  enddo

  ! Si no existe, creo un nuevo boundgr
  call boundgrs%add_after()
  call boundgrs%next%alloc()
  bg=>boundgrs%next%o
  call bg%init(z)

end function boundgrs_include

subroutine boundgr_init(bg,z)
  class(boundgr)   :: bg
  integer,intent(inout)     :: z(:)

  allocate(bg%list)
  call bg%list%init()

  allocate(bg%z,source=z)
  call bg%i%init()
  call bg%p%init()

  ! Set cero the first integer in case it is needed to use it as a counter
  bg%i%o(1)=0
   
end subroutine boundgr_init

subroutine boundgr_destroy(bg)
  class(boundgr)   :: bg
  ! type(boundgr)   :: bg
 
  call bg%list%destroy_all()
  deallocate(bg%list)

  deallocate(bg%z)
  call bg%i%destroy()
  call bg%p%destroy()

end subroutine boundgr_destroy

subroutine boundgr_add(bg,id)
  use gems_algebra, only:sort_int
  use gems_errors, only:werr
  ! Add a bound to the corresponding boundgroup. If there is not a boundgroup for
  ! this bound it creates onw.
  class(boundgr)            :: bg
  type(bound),pointer       :: b
  integer,intent(inout)     :: id(:)
  integer                   :: i

  ! TODO: delete when ready
  call werr('Trying to add different bound sizes',(size(bg%z)/=size(id)))

  ! XXX: Quizas tendria que relizar una busqueda para ver si ya se encuentra
  ! allocateado el enlace y en ese caso añadirlo soft. Capaz se podría hacer si
  ! se agregan de manera ordenada siguiente un uid que sea el hash de los dos
  ! atoms de manera de buscar de forma eficiente
  call bg%list%add_after()
  call bg%list%next%alloc()
  bg%n = bg%n + 1
  b=>bg%list%next%o
  allocate(b%a(size(id)))

  do i=1,size(id)
    b%a(i)%o => sys%a(id(i))%o
  enddo

end subroutine boundgr_add


! subroutine impropers(f,epot)
! ! OPLS improper
! use gems_inq_properties, only:vdistance
! use gems_input_parsing, only:ioprefix
! class(ffunction),intent(in) :: f
! type(ffnode),pointer        :: fn
! type(atom),pointer          :: i,j,k,l
! integer                     :: m
! real(dp),intent(out)        :: epot
! real(dp)                   :: dud1(dm),dud2(dm),dud3(dm),dud4(dm)
! real(dp)                   :: vij(dm),vjk(dm),vkl(dm),dnum(dm),dden(dm)
! real(dp)                   :: bxc(3)
! real(dp)                   :: ab,bc,at,bt,ct,axb_mod,bxc_mod,rnum
! real(dp)                   :: ac,cosphi,cphi,cphi2,sphi,signo,ptors,ftors,de1,den,phi
! real(dp)                   :: k0,n0,phis
!
! epot = 0.0_dp
!
! fn => f%head
! do m=1,f%nmatch
!   fn => fn%next
!
!   ! los atomos
!   i => fn%a
!   j => fn%b
!   k => fn%c
!   l => fn%d
!
!   ! los parametros
!   phis = f%param(fn%id,1)
!   k0 = f%param(fn%id,2)
!   n0 = f%param(fn%id,3)
!
!   ! Distancias
!   vij=vdistance(j,i, mic)
!   vjk=vdistance(k,j, mic)
!   vkl=vdistance(l,k, mic)
!
!   ! Producto cruz b con c
!   bxc(1)=vjk(2)*vkl(3)-vjk(3)*vkl(2)
!   bxc(2)=vjk(3)*vkl(1)-vjk(1)*vkl(3)
!   bxc(3)=vjk(1)*vkl(2)-vjk(2)*vkl(1)
!   bxc(:)=bxc(:)-box(:)*idnint(bxc(:)*one_box(:))
!
!   ! Productos escalares cruzados
!   ab=dot(vij,vjk)
!   bc=dot(vjk,vkl)
!   ac=dot(vij,vkl)
!
!   ! Modulo cuadrado de cada distancia
!   at=dot(vij,vij)
!   bt=dot(vjk,vjk)
!   ct=dot(vkl,vkl)
!
!   ! Productors vectoriales al cuadrado entre distancias
!   axb_mod=at*bt-ab*ab
!   bxc_mod=bt*ct-bc*bc
!
!   !
!   rnum=ab*bc-ac*bt ! El producto cuadrado de los vectoriales?
!   den=axb_mod*bxc_mod      ! El producto cuadrado de los vectoriales
!
!   !chequeo de que átomos no están en la misma línea
!
!   if(den<1.d-12) cycle
!
!   den=dsqrt(den) ! El producto de los vectoriales
!
!   !coseno del ángulo
!   cosphi=rnum/den
!   cosphi=dmin1(cosphi,1.0_dp)
!   cosphi=dmax1(cosphi,-1.0_dp)
!
!   !signo del ángulo
!   signo=dot_product(vij,bxc)
!
!   !definición del ángulo
!   phi=dsign(dacos(cosphi),signo)
!
!   !if (impropia) phi = n*phi -phi_s
!   !phis(i)=phi
!
!   ! Ryckaert−Bellemans (RB) potential function
!   cphi=dcos(phi-pi)
!   sphi=dsin(phi-pi)
!   cphi2=cphi*cphi
!
!   ptors=k0*(1+dcos(n0*phi-phis))
!
!   ftors=n0*k0*dsin(n0*phi-phis)
!
!
!   ! debug para testear topologias. Fija una determinada torsion.
!   if(m==tfix) then
!     ptors=0.1_dp*(dabs(phi)-dabs(phifix))**2
!     ftors=-0.2_dp*dabs((dabs(phi)-dabs(phifix)))
!   endif
!
!   i%epot = i%epot + ptors*ev_ui*0.25_dp
!   j%epot = j%epot + ptors*ev_ui*0.25_dp
!   k%epot = k%epot + ptors*ev_ui*0.25_dp
!   l%epot = l%epot + ptors*ev_ui*0.25_dp
!
!   epot = epot + ptors*ev_ui
!
!   if(dabs(sphi)<1.d-12) sphi=dsign(1.d-12,sphi)
!   de1=ftors/den/sphi
!   axb_mod=axb_mod/den*cosphi
!   bxc_mod=bxc_mod/den*cosphi
!
!   dnum=vkl*bt-vjk*bc
!   dden=(ab*vjk-vij*bt)*bxc_mod
!   dud1=(dnum-dden)*de1
!
!   dnum=((vjk-vij)*bc-ab*vkl)+(2.0_dp*ac*vjk-vkl*bt)
!   dden=axb_mod*(bc*vkl-vjk*ct)+(vij*bt-at*vjk-ab*(vjk-vij))*bxc_mod
!   dud2=(dnum-dden)*de1
!
!   dnum=ab*vjk-vij*bt
!   dden=axb_mod*(bt*vkl-bc*vjk)
!   dud4=(dnum-dden)*de1
!
!   dud3=-(dud1+dud2+dud4)
!
!
!   i%force = i%force + dud1*ev_ui
!   j%force = j%force + dud2*ev_ui
!   k%force = k%force + dud3*ev_ui
!   l%force = l%force + dud4*ev_ui
!
!
!  end do
!
!  end subroutine impropers

subroutine write_ebend(op)
  class(outpropa)   :: op
  type(boundgr_l),pointer         :: ln
  
  ln => boundgrs
  do while( associated(ln%next) )
    ln => ln%next
    if(size(ln%o%z)/=3) cycle
    op%f(1) = op%f(1) + ln%o%epot*ui_kcm
  enddo 
     
end subroutine

subroutine write_etors(op)
  class(outpropa)   :: op
  type(boundgr_l),pointer         :: ln
  
  ln => boundgrs
  do while( associated(ln%next) )
    ln => ln%next
    if(size(ln%o%z)/=4) cycle
    op%f(1) = op%f(1) + ln%o%epot*ui_kcm
  enddo 
     
end subroutine
                
subroutine write_estretch(op)
  class(outpropa)   :: op
  type(boundgr_l),pointer         :: ln
  
  ln => boundgrs
  do while( associated(ln%next) )
    ln => ln%next
    if(size(ln%o%z)/=2) cycle
    op%f(1) = op%f(1) + ln%o%epot*ui_kcm
  enddo 
     
end subroutine
                       
end module gems_forcefield



! References for the all-atom additive CHARMM force fields

!
! Empirical Force Field Reviews
!
! MacKerell, Jr., A.D. "Empirical Force Fields for Biological Macromolecules: Overview and Issues," Journal of Computational Chemistry 25: 1584-1604, 2004, [DOI]
! Lopes, P.E.M., Guvench, O., and MacKerell, A.D., Jr., Current Status of Protein Force Fields for Molecular Dynamics, In Molecular Modeling of Proteins, 2nd edition, A. Kukol, Editor, Humana Press. Chapter 3, pp. 47 -72, 2014, [DOI], (Methods in Molecular Biology, 1215: 47-72, 2014).
! Guvench, O. and MacKerell, Jr. A.D. "Comparison of protein force fields for molecular dynamics simulations," In Molecular Modeling of Proteins, A. Kukol, Ed. Humana Press, Humana Press. 2008 (Methods Mol Biol. 2008;443:63-88. PMID: 18446282)

!
! CHARMM General FF (CGenFF)
!
! K. Vanommeslaeghe, E. Hatcher, C. Acharya, S. Kundu, S. Zhong, J. Shim, E. Darian, O. Guvench, P. Lopes, I. Vorobyov, A. D. MacKerell Jr., "CHARMM General Force Field (CGenFF): A force field for drug-like molecules compatible with the CHARMM all-atom additive biological force fields," Journal of Computational Chemistry 31: 671-690, 2010.
!
! W. Yu, X. He, K. Vanommeslaeghe, and A. D. MacKerell Jr., "Extension of the CHARMM general force field to sulfonyl-containing compounds and its utility in biomolecular simulations," Journal of Computational Chemistry, 33: 2451-2468, 2012.
!
! Before downloading please read this warning!
!
! toppar_all36a_cgenff.tgz
!
! The force field can automatically be applied to an arbitrary organic molecule using the CGenFF program, which can be conveniently be accessed through the cgenff.paramchem.org web interface. Click here for usage information. The resulting parameters and charges are accompanied by penalty scores. If these penalty scores are high, it is recommended to re-optimize the parameters, as described in the above reference and the tutorial.
!
! Frequent users of the CGenFF program may wish to obtain a binary license. The procedure for obtaining a free-of-charge not-for-profit license is initiated by e-mailing us; it may take up to a few weeks and will require someone with signature authority at your institution to sign a license agreement that needs to be sent back to us.
!
! For-profit users may obtain the CGenFF program from SilcsBio, LLC.
!
!
!
! Proteins
!
! Best, R.B., Zhu, X., Shim, J., Lopes, P.E.M., Mittal, J., Feig, M., and MacKerell Jr., A.D. "Optimization of the additive CHARMM all-atom protein force field targeting improved sampling of the backbone phi, psi and side-chain chi1 and chi2 dihedral angles," Journal of Chemical Theory and Computation, 8: 3257-3273, 2012, [DOI] PMC3549273
! MacKerell, Jr., A. D., Bashford, D., Bellott, M., Dunbrack Jr., R.L., Evanseck, J.D., Field, M.J., Fischer, S., Gao, J., Guo, H., Ha, S., Joseph-McCarthy, D., Kuchnir, L., Kuczera, K., Lau, F.T.K., Mattos, C., Michnick, S., Ngo, T., Nguyen, D.T., Prodhom, B., Reiher, III, W.E., Roux, B., Schlenkrich, M., Smith, J.C., Stote, R., Straub, J., Watanabe, M., Wiorkiewicz-Kuczera, J., Yin, D., ad Karplus, M. "All-atom empirical potential for molecular modeling and dynamics studies of proteins," Journal of Physical Chemistry B 102: 3586-3616, 1998, [DOI]
! Detailed presentation of the CMAP procedure and its use to treat the conformational properties of the protein backbone:
!
! MacKerell, Jr., A.D., Feig, M., and Brooks, III, C.L. "Extending the treatment of backbone energetics in protein force fields: limitations of gas-phase quantum mechanics in reproducing protein conformational distributions in molecular dynamics simulations," Journal of Computational Chemistry, 25: 1400-1415, 2004. [DOI]
!
!
! Nucleic Acids
!
! C36 DNA
!
! Hart, K., Foloppe, N., Baker, C.M., Denning, E.J., Nilsson, L. and MacKerell, A.D., Jr., "Optimization of the CHARMM additive force field for DNA: Improved treatment of the BI/BII conformational equilibrium," Journal of Chemical Theory and Computation, 8: 348-362, 2012, [DOI] PMC3285246
! C36 RNA
!
! Denning, E.J. Priyakumar, U.D. Nilsson, L. and MacKerell, Jr. A.D. "Impact of 2'-hydroxyl sampling on the conformational properties of RNA: Update of the CHARMM all-atom additive force field for RNA," Journal of Computational Chemistry 32: 1929-1943, 2011, [DOI]
! C27 RNA and DNA
!
! Foloppe, N. and MacKerell, Jr., A.D. "All-Atom Empirical Force Field for Nucleic Acids: 1) Parameter Optimization Based on Small Molecule and Condensed Phase Macromolecular Target Data," Journal of Computational Chemistry 21: 86-104, 2000, [DOI]
! The supplemental material of Foloppe and MacKerell is actually the full, unabridged version and can be obtained from http://journals.wiley.com/jcc/. As the published manuscript is an abbreviated version, it is strongly suggested that the full version be obtained in order to get all the details of the parameter optimization procedure.
!
! MacKerell, Jr., A.D. and Banavali, N. "All-Atom Empirical Force Field for Nucleic Acids: 2) Application to Molecular Dynamics Simulations of DNA and RNA in Solution," Journal of Computational Chemistry 21: 105-120, 2000. [DOI] Also see http://journals.wiley.com/jcc/ for the supplemental material or view versions of Supplemental Material Figure 6a and Figure 6b here.
!
!
! Lipids
!
! C36 lipids
!
! Klauda, J.B., Venable, R.M., Freites, J.A., O'Connor, J.W., Tobias, D.J., Mondragon-Ramirez, C., Vorobyov, I., MacKerell, Jr., A.D., and Pastor, R.W. "Update of the CHARMM All-Atom Additive Force Field for Lipids: Validation on Six Lipid Types," Journal of Physical Chemistry B, 114: 7830-7843, 2010, [DOI]
! C22 and C27 lipids
!
! Schlenkrich, M., Brickmann, J. MacKerell, Jr., A.D. and Karplus, M. "An Empirical Potential Energy Function for Phospholipids: Criteria for Parameter Optimization and Applications," in Biological Membranes: A Molecular Perspective from Computation and Experiment K.M. Merz, Jr. and B. Roux, Eds. Birkhauser, Boston, 1996.
! Feller., S.E., Yin, D., Pastor, R.W. and MacKerell, Jr., A.D. "Molecular Dynamics Simulation of Unsaturated Lipids at Low Hydration: Parametrization and Comparison with Diffraction Studies," Biophysical Journal 73: 2269-2279, 1997.
!
!
! Carbohydrates
!
! Guvench, O., Mallajosyula, S.S. Raman, E.P., Hatcher, E. Vanommeslaeghe, K., Foster, T.J., Jamison II, F.W., and MacKerell, A.D., Jr. "CHARMM additive all-atom force field for carbohydrate derivatives and their utility in polysaccharide and carbohydrate-protein modeling," Journal of Chemical Theory and Computing, 7: 3162-3180, 2011, [DOI], PMC3224046
! Mallajosyula, S.S., Guvench, O., Hatcher, E., MacKerell, A.D., Jr., "CHARMM Additive All-Atom Force Field for Phosphate and Sulfate Linked to Carbohydrates," Journal of Chemical Theory and Computation, 8: 759-776, 2012, [DOI], PMC3367516
! Raman, P. Guvench, O. MacKerell, Jr. A.D. "CHARMM Additive All-Atom Force Field for Glycosidic Linkages in Carbohydrates Involving Furanoses," Journal of Physical Chemistry B 114: 12981-12994, 2010, [DOI], PMC2958709
! Guvench, O. Hatcher, E.R. Venable, R.M. Pastor, R.W. and MacKerell Jr., A.D. "Additive Empirical CHARMM Force Field for glycosyl linked hexopyranoses," Journal of Chemical Theory and Computation 5: 2353-2370, 2009, [DOI], PMC2757763
!
