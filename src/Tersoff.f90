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

 
module gems_tersoff
! Sacado del paper original.
! Utilizo el promedio de los Bon Angle Term para hacer bucle con i>j
! Como introduccion ver el codigo bond order del GULP.
! Para enquilombarse la vida ver Bitacora 09/03/011 y sucesivas (ahi tengo
! algunas sugerencias que quizas sirvan para mejorar)

use gems_program_types
use gems_constants
use gems_tables
use gems_inq_properties
use gems_output
use gems_neighbour

implicit none

private

public  :: tsf_new
!public gradb_k

type,extends(intergroup) :: tsf
  real(dp),allocatable  :: gradb_k(:,:)
  contains
  procedure :: interact => tsf_interact
  procedure,nopass :: cli => intergroup0_empty
end type


! Estas variables son para escribir
logical,private   :: noadv=.true.
real(dp),private  :: var
public            :: write_tsfangs

!-----------------

!! Parametros Silicio
!real(dp),parameter       :: ja=1.8308e+3, B=4.7118e+2, lambda1=2.4799
!                            lambda2=1.7322, bt=1.1000e-6, np=7.8734e-1,  
!                            c=1.0039e+5, d=1.6217e+1, h=-5.9825e-1, 
!                            r1=2.7, r2=3.0, chi=1.0

! Parametros Carbono
!Lindsay, L. 2010. Optimized Tersoff and Brenner empirical potential parameters
!for lattice dynamics and phonon thermal transport in carbon nanotubes and
!graphene. Physical Review B 81, no. 20 (May): 1-6.
!doi:10.1103/PhysRevB.81.205441.
!http://link.aps.org/doi/10.1103/PhysRevB.81.205441.
!real(dp),parameter       :: ja=1393.6_dp*ev_ui  ,&
!                            B = 346.74_dp*ev_ui ,& 
!                            lambda1 = 3.4879_dp   ,&
!                            lambda2 = 2.2119_dp      ,&
!                            np = 0.72751_dp     ,& 
!                            bt = 1.5724e-7_dp   ,&
!                            c = 38049.0_dp      ,& !the strength of the angular effect
!                            d = 4.6484_dp       ,& !how sharp the dependence on angle is
!                            h = -0.57058_dp     ,& !g has a minilambda2m for h
!                            r1 = 1.8_dp         ,& !R=1.95, D=0.15
!                            r2 = 2.0_dp         ,&
!                            chi=1.0_dp         

!Paper original tersoff
real(dp),parameter       :: ja=1393.6_dp*ev_ui  ,&
                            B = 346.74_dp*ev_ui ,& 
                            lambda1 = 3.4879_dp   ,&
                            lambda2 = 2.2119_dp      ,&
                            np = 0.72751_dp     ,& 
                            bt = 1.5724e-7_dp   ,&
                            c = 38049.0_dp      ,& !the strength of the angular effect
                            d = 4.3484_dp       ,& !how sharp the dependence on angle is
                            h = -0.57058_dp     ,& !g has a minilambda2m for h
                            r1 = 1.8_dp         ,& !R=1.95, D=0.15
                            r2 = 2.1_dp         ,&
                            chi=1.0_dp             

real(dp),public,parameter       :: tsf_rcut=r2
  
real(dp),parameter       :: btnp = bt**np   ,& 
                            d2 = d*d        ,& 
                            c2 = c*c  

real(dp),dimension(dm)             :: gradb_i,gradb_j ! Terminos de la fuerza por el gradiente de b
real(dp),dimension(dm),allocatable :: gradb_k(:,:) 


logical,public   :: get_tsfangs

contains

                     
subroutine tsf_new(pg,g)
class(intergroup),pointer  :: pg
type(tsf),pointer          :: ig
type(group),intent(in)     :: g

! Return a intergroup class pointer
allocate(ig)
pg=>ig

! Init the internal tersoff intergroup and request neighbor list with
call ig%init(g1=g)

! Solo carbono
! do i = 1,zmax
!   if(i==6) cycle
!   call error(comp1(i),1020)
! enddo

allocate(ig%gradb_k(ig%n(4),dm))

end subroutine tsf_new

subroutine tsf_interact(ig)
! No es posible escribir una interaccion del tipo s1-s2 en un potencial como
! el Tersoff.
! La notación y el desarrollo trata de seguir la descripcion echa en:
! Tersoff, J. 1988. New empirical approach for the structure and energy of
! covalent systems. Physical Review B 37, no. 12: 6991–7000.
! http://link.aps.org/doi/10.1103/PhysRevB.37.6991.
class(tsf),intent(inout)  :: ig
real(dp)                  :: fc,dfc                  ! Funcion de corte y derivada
real(dp)                  :: bij                     ! Bond Angle Term, es lo que hace covalente al Tersoff
real(dp),dimension(dm)    :: v                       ! Versores
real(dp),dimension(dm)    :: vd                      ! Distancia Vectorial
real(dp)                  :: rd,factor               ! Distancia
real(dp)                  :: va,vr,dva,dvr           ! Pot de a pares y derivadas
real(dp)                  :: dv1,dv2         
integer                   :: i,j,l,m,n,k
type(atom_dclist),pointer :: la


ig%epot = 0.0_dp

la => ig%a     ! Bucle sobre los atomos (=> i)
do l = 1,ig%n(1)
  la=>la%next

  i  = la%o%tag

  do m = 1, ig%nn(i)  ! Bucle sobre los vecinos de i (=> j/=i)

    j  = ig%list(i,m)
    j  = ig%at(j)%o%tag


    ! Calculo las distancias
    vd = vdistance( a(i)%o, a(j)%o , mic) 
    rd =  dsqrt(dot_product(vd,vd)) 

    ! Excluye por radio de corte
    if(rd>r2) cycle

    ! Potencial de a pares atractivo  y derivada
    va=-b*exp(-lambda2*rd)  
    dva=-lambda2*va

    ! Potencial de a pares repulsivo y derivada
    vr=ja*exp(-lambda1*rd) 
    dvr=-lambda1*vr

    ! Funcion de corte y derivada
    fc=fcut_dfcut(rd,r1,r2,dfc) 

    ! Calculo del bond angle term simetrizado (promedio)
    v=vd/rd
    call bond_angle_term(ig,i,j,v,rd,bij)

    ! Energia
    factor = 0.5_dp*fc*(vr+bij*va)

    !! Si bien bij no es bji, cuadno termine el do cada atomo va a tener una
    !! energia 0.5*fc*( vr+va( (bij+bji)*0.5 ) ), es decir un b promedio simetrico.
    !a(j)%o%epot = a(j)%o%epot + factor*0.5_dp

    a(i)%o%epot = a(i)%o%epot + factor!*0.5_dp
    ig%epot = ig%epot + factor

    ! Factor de fuerza direccion i<-j
    dv1=0.5_dp*(dfc*(vr+bij*va)+fc*(dvr+bij*dva))
    
    ! Factor de fuerza en direccion del gradiente del angulo
    dv2=0.5_dp*fc*va 

    ! Fuerza
    a(i)%o%force = a(i)%o%force - dv1*v + dv2*gradb_i

    ! Nuevamente cuando el do termine el 0.5 sirve para (bij+bji)/2
    a(j)%o%force = a(j)%o%force + dv2*gradb_j
    !a(j)%o%force = a(j)%o%force + dv1*v + dv2*gradb_j*0.5_dp

    do n = 1, ig%nn(i)  ! Bucle sobre todos los otros vecinos de i(=> k/=j/=i)
      ! Esto es solo sobre los vecinos de i porque el potencial tiene la funcion
      ! de corte en relacion a i
      k  = ig%list(i,n)

      !Excluye por autointeraccion
      if(k==j) cycle

      !Distancia i<-k
      vd = vdistance( a(i)%o, a(k)%o , mic) 
      rd = dsqrt(dot_product(vd,vd)) 
   
      !excluye por radio de corte
      if(rd>r2) cycle

      a(k)%o%force = a(k)%o%force + dv2*ig%gradb_k(k,:)

    enddo
    
  enddo
enddo

end subroutine

subroutine bond_angle_term(ig,i,j,vij,rij,bat_ij)
class(tsf)                   :: ig
integer,intent(in)           :: i,j
real(dp),intent(in)          :: vij(dm),rij
real(dp),intent(out)         :: bat_ij
real(dp),dimension(dm)       :: v,vd,gradct_i,gradct_j,gradct_k
real(dp)                     :: rd,fc,zij,dfc,ct,aux,g,dg,fzij
integer                      :: n,k

! Devuelve el bond angle term (bat_ij) y su gradiente

! Inicializo sumatorias
zij=0.0_dp
gradb_i=0.0_dp
gradb_j=0.0_dp
gradb_k=0.0_dp
 
 
noadv=.true. ! Para escribir en la linea

do n = 1, ig%nn(i)  ! Bucle sobre todos los otros vecinos de i(=> k/=j)
  ! Esto es solo sobre los vecinos de i porque el potencial tiene la funcion
  ! de corte en relacion a i
  k = ig%list(i,n)

  !Excluye por autointeraccion
  if(k==j) cycle

  !Distancia i<-k
  vd = vdistance( a(i)%o, a(k)%o , mic) 
  rd = dsqrt(dot_product(vd,vd)) 

  !excluye por radio de corte
  if(rd>r2) cycle

  !Versor i<-k
  v=vd/rd

  !Coseno del angulo kij (el angulo sobre i)
  ct=dot_product(vij,v)  
  ct=max(ct,-1.0_dp)
  ct=min(ct,1.0_dp)
  !if(ct<-1.0_dp) ct=-1.0_dp
  !if(ct>+1.0_dp) ct=+1.0_dp

  if (get_tsfangs) then
    var= acos(ct)
    call write_out(4,dm_steps)
    if(n==ig%nn(i)) noadv=.false.
  endif

  !Funcion de corte y derivada respecto a componente rik
  fc=fcut_dfcut(rd,r1,r2,dfc)

  !Funcion G(ct) 
  aux=1.0_dp/(d2+(h-ct)**2) 
  g=1.0_dp+c2/d2-c2*aux                       

  !Funcion zeta_ij 
  zij=zij+fc*g  

  !Derivada de G(ct) respecto a ct
  dg=-2.0_dp*c2*(h-ct)*aux*aux

  !expterm   = EXP((lambda3*(drij-drik))**3)
  !expterm_d = (3.0_dp)*(lambda3**3)*((drij-drik)**2)*expterm

  !Menos gradiente de ct_ijk respecto a j
  gradct_j=(v-ct*vij)/rij

  !Menos gradiente de ct_ijk respecto a k
  gradct_k=(vij-ct*v)/rd

  !Gradiente de ct_ijk respecto a i
  gradct_i=-(gradct_k+gradct_j)

  !Termino del gradiente de zij a sumarle directamente a i 
  gradb_i=gradb_i+dfc*g*v+fc*dg*gradct_i
  !+ f_C * gterm * expterm_d * (-rij_hat+rik_hat)

  !Termino del Gradiente de zij respecto a j a sumarle a j 
  ! (Adelanto en la sumatoria para completar la fuerza en i cuando se
  ! inviertan las etiquetas)
  gradb_j=gradb_j+fc*dg*gradct_j
  !+ f_C   * gterm   * expterm_d * (rij_hat)

  ! Termino del Gradiente de zij respecto a k a sumarle a k 
  ! (Adelanto en la sumatoria para completar la fuerza en i cuando se
  ! inviertan las etiquetas)
  ig%gradb_k(k,:)=-dfc*g*v+fc*dg*gradct_k
  !+ f_C   * gterm   * expterm_d * (-rik_hat)

enddo ! fin atomo k

if(zij>1.0e-15_dp) then 
  aux=1.0_dp/(1.0_dp+btnp*zij**np)
  bat_ij=chi*aux**(0.5_dp/np)
  fzij=-0.5_dp*bat_ij*aux*btnp*zij**(np-1.0_dp)
  gradb_i=fzij*gradb_i
  gradb_j=fzij*gradb_j
  ig%gradb_k(:,:)=fzij*gradb_k(:,:)
else
  bat_ij=1.0_dp
  gradb_i=0.0_dp
  gradb_j=0.0_dp
  ig%gradb_k(:,:)=0.0_dp
endif

end subroutine

      ! ESCRITURA

subroutine write_tsfangs(op,un)
use gems_output
class(outpropa)     :: op
integer,intent(in)  :: un

write(un,fmt='(e20.12)',advance='no') var
!op%adv=noadv

end subroutine

!  subroutine bo_cm_interaction(s1,s2) 
!   ! Esta subrutina asume que s1 es carbono y s2 metal
!   integer                 :: i,j,l,is1,is2,k
!   real(dp)                :: rd,nc,vr,va
!   real(dp)                :: factor
!   real(dp),dimension(dm)  :: factor2
!   type(subsystem)         :: s1,s2
!    
!   is1=s1%id
!   is2=s2%id
!
!   do i = 1, s1%nat        ! sobre los atomos
!
!     do l = 1 , s1%nnb(is2,i)  ! sobre los vecinos
!       j  = s1%nb(is2,i,l)
!       k = z_boi( s1%a(i)%z,s2%a(j)%z )
!
!       rd = s1%nbrd(is2,i,l)
! 
!       if(rd<e2s(k)%fin) then
!
!          nc=s1%a(i)%cord 
!
!          vr=linear_inter(e2s(k),rd)*linear_inter(f2s(k),nc)
!          va=linear_inter(es2(k),rd)*linear_inter(fs2(k),nc)
!
!          ! force calculation (restricted)
!          factor=vr*sqrt_s2(k)-va*sqrt_2s(k)
!          factor2(1:dm) = factor*s1%vd(1:dm)/rd
!          s1%a(i)%force(1:dm) = s1%a(i)%force(1:dm) + factor2(1:dm) 
!          s2%a(j)%force(1:dm) = s2%a(j)%force(1:dm) - factor2(1:dm)
!           
!          ! energy calculation
!          factor=(va-vr)*0.5
!          s1%a(i)%epot = s1%a(i)%epot + factor 
!          s2%a(j)%epot = s2%a(j)%epot + factor 
!        
!        endif
!      enddo
!    enddo
!
!  end subroutine
!
!  subroutine bo_mm_interaction(s1,s2) 
!   ! Esta subrutina asume que s1 y s2 son metales
!   integer                 :: i,j,k,l,is1,is2
!   real(dp)                :: rd,nc,vr,va
!   real(dp)                :: factor
!   real(dp),dimension(dm)  :: factor2
!   type(subsystem)         :: s1,s2
!    
!   is1=s1%id
!   is2=s2%id
!
!
!  end subroutine
!

end module gems_tersoff
