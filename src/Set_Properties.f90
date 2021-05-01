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


module gems_set_properties
 use gems_program_types, only: dt,a,natoms,gsel,sys
 use gems_groups, only: group, group_l,glist
 use gems_atoms, only: atom,atom_dclist
 use gems_constants,only:dp,pi,dm,find_io
 use gems_algebra
 use gems_inq_properties

  implicit none

  public
  private :: group,atom,atom_dclist,gsel,sys   ! program_types
  private :: dt,dm                             ! program_types
  private :: dp,pi                             ! constants

  save

  real(dp)     :: dc,tc

  interface do_pbc
    module procedure    :: do_pbc_1,do_pbc_2
  end interface

 contains

!                                     control de booleanos
!----------------------------------------------------------

subroutine mass_changed
class(group_l),pointer :: lg
lg => glist%next
do while(associated(lg))
  lg%o%b_mass    =.false.
  lg => lg %next
enddo
end subroutine mass_changed

subroutine vel_changed
class(group_l),pointer :: lg

lg => glist%next
do while(associated(lg))
  lg%o%b_ekin    =.false.
  lg%o%b_temp    =.false.
  lg%o%b_tempvib =.false.
  lg%o%b_temprot =.false.
  lg%o%b_erot    =.false.
  lg%o%b_evib    =.false.
  lg%o%b_cm_vel  =.false.
  lg%o%b_ang_mom =.false.
  lg%o%b_ang_vel =.false.
  lg%o%b_pressure =.false.
  lg => lg %next
enddo
end subroutine vel_changed

subroutine pos_changed
class(group_l),pointer :: lg

lg => glist%next
do while(associated(lg))
  lg%o%b_erot     =.false.
  lg%o%b_evib     =.false.
  lg%o%b_cm_vel   =.false.
  lg%o%b_ang_mom  =.false.
  lg%o%b_ang_vel  =.false.

  lg%o%b_maxpos   =.false.
  lg%o%b_minpos   =.false.
  lg%o%b_mainaxis =.false.
  lg%o%b_cm_pos   =.false.
  lg%o%b_rg_pos   =.false.
  lg%o%b_covar    =.false.
  lg%o%b_inercia  =.false.

  lg%o%b_virial =.false.
  lg%o%b_pressure =.false.

  lg=>lg%next
enddo

call epot_changed

end subroutine pos_changed

subroutine epot_changed
class (group_l),pointer :: lg

lg => glist%next
do while(associated(lg))
  lg%o%b_epot  =.false.
  lg%o%b_virial =.false.
  lg=>lg%next
enddo

end subroutine epot_changed

subroutine posvel_changed
class(group_l),pointer :: lg

lg => glist%next
do while(associated(lg))
  lg%o%b_ekin    =.false.
  lg%o%b_temp    =.false.
  lg%o%b_tempvib =.false.
  lg%o%b_temprot =.false.
  lg%o%b_erot    =.false.
  lg%o%b_evib    =.false.
  lg%o%b_cm_vel  =.false.
  lg%o%b_rg_pos   =.false.
  lg%o%b_ang_mom =.false.
  lg%o%b_ang_vel =.false.

  lg%o%b_maxpos   =.false.
  lg%o%b_minpos   =.false.
  lg%o%b_mainaxis =.false.
  lg%o%b_cm_pos   =.false.
  lg%o%b_covar    =.false.
  lg%o%b_inercia  = .false.

  lg%o%b_virial =.false.
  lg%o%b_pressure =.false.

  lg => lg%next
enddo

call epot_changed

end subroutine posvel_changed

!                                                 elementos
!----------------------------------------------------------

subroutine set_z(g,z)
type(group)               :: g
class(atom_dclist),pointer :: la
integer,intent(in)        :: z
integer                   :: i

la => g%alist%next
do i = 1,g%nat
  call la%o%setz(z)
  la => la%next
enddo

call mass_changed()

end subroutine set_z

subroutine set_mass(g,f)
type(group)                :: g
class(atom_dclist),pointer :: la
real(dp),intent(in)        :: f
integer                    :: i

la => g%alist%next
do i = 1,g%nat
  la%o%mass = f

  if(la%o%mass==0._dp) then
    la%o%one_mass = 0._dp
    la%o%one_sqrt_mass = 0._dp
  else
    la%o%one_mass = 1.0_dp/la%o%mass
    la%o%one_sqrt_mass = sqrt(la%o%one_mass)
  endif

  la => la%next
enddo

call mass_changed()

end subroutine set_mass

subroutine set_molid(g,id)
type(group)               :: g
class(atom_dclist),pointer :: la
integer,intent(in)        :: id
integer                   :: i

la => g%alist%next
do i = 1,g%nat
  !if(la%o%molid/=0) then
  !  werr 'un atomo en mas de una molecula?'
  !endif
  ! FIXME: Pero si aplico dos topoligas distintas a una molecula eso no es
  ! cierto
  la%o%molid=id
  la%o%amolid=i
  la => la%next
enddo
end subroutine set_molid

subroutine set_sigma(g,s)
type(group)                :: g
class(atom_dclist),pointer :: la
real(dp),intent(in)        :: s
integer                    :: i

la => g%alist%next
do i = 1,g%nat
  la%o%s=s
  la => la%next
enddo
end subroutine set_sigma

!                                               temperatura
!----------------------------------------------------------

! NOTE that if a lattice is started with all of the atoms perfectly place, wich
! is the lowest possivle energy position for the lattice, about half of the
! kinetic energy added goes into the potential energy needed to place the atom
! off their perfect positions. After a short period of time (say 200 time step),
! the temperature of the lattice will be half of the temperature when the
! initial temp command was first issued. (extracted from xmd manual)

  function ptemp_gexp(amp,lon,step,tfin,nc)
   real(dp)            :: ptemp_gexp,amp,tfin
   integer             :: step,nc,c,lon
    !devuelve la temperatura correspondiente a cada instante par un programa de
    !temperaturas de "nc" ciclos de decaimientos exponenciales con amplitudes
    !moduladas por una exponencial global. cada ciclo tiene una duracion igual a
    !"lon". ojo, esta subrutina se llama recursivamente

    !    |
    !    | |i
    !    | | \         |l
    !    | |  i        | \
    !    | |   \       |  l     |l
    !    | |    `_     |   `_   | \
    !    | |      - _  |     -  |  `_
    !    |________________________________

    !encuentro la amplitud de la exitacion para el c actual
    c = mod(step,lon)
    tc = amp*(tfin/amp)**(dfloat(c)/nc)

    !exponencial en cada ciclo
    if (step>lon) step = step - lon
    ptemp_gexp = ptemp_exp(tc,dc,step)

  end function ptemp_gexp

  function ptemp_exp(amp,lon,step)
   real(dp)            :: ptemp_exp,amp,lon
   integer             :: step
    !devuelve la temperatura correspondiente a cada instante de un programa de
    !temperaturas con un decaimiento exponenciales que cumple con t(0)=amp,
    !t(lon/2)=amp/3, t(lon)=0. el tiempo se mide en pasos.  ojo, esta subrutina
    !se llama recursivamente

    !  t | |i
    !    | | \
    !    | |  i
    !    | |   \
    !    | |    `_
    !    | |      - _
    !    |_____________  t

    ptemp_exp = amp/3.0_dp*(4.0_dp*4.0_dp**(-step/lon)-1)

  end function ptemp_exp


!                               distribucion de velocidades
!----------------------------------------------------------
  subroutine set_add_cmvel(g,f)
    type(group)                :: g
    class(atom_dclist),pointer  :: la
    real(dp),intent(in)        :: f(dm)
    integer                    :: j

    call group_inq_cmpos(sys)
    la => g%alist%next
    do j =1, g%nat
      la%o%vel = la%o%vel + f
      la => la%next
    enddo

    call vel_changed()
  end subroutine set_add_cmvel

subroutine set_gdist(g,newtemp)
use gems_random, only:rang
use gems_constants
real(dp)            :: factor,newtemp,vel_med,f(dm)
real(dp)            :: r1
class(group)        :: g
class(atom_dclist),pointer   :: la
integer              :: i,j

! Esto se usa para sacar la traslacion del CM
call inq_cm_vel(g)

! Sacando la traslacion del CM, los grados de libertad son (dm*g%nat-dm), entocnes:
!   Ec=(dm*n-dm)/2*k*T
! N particulas de energia cinetcia 1/2mv**2 tiene una Ec:
!   Ec=n/2*m*vel2_med
! De estas dos despejo la velocidad cuadratica media
!   vel2_med=(dm*n-dm)*k*T/(n*m)
! Sabiendo que vel2_med=vx2_med+vy2_med+vz2_med y como no tengo razones para
! suponer que las componentes difieren tengo que vx2_med=vel2_med/dm y queda:
!   vx2_med=(dm*n-dm)*k*T/(n*m*dm)
!   vx_med=sqrt((dm*n-dm)*k*T/(n*m*dm))
! Para obtener la distribucion de maxwell, cada componente debe seguir una
! distribucion gaussiana, cuya norma sea la vx_med

factor = kB_ui*newtemp!/(g%nat*dm)*(dm*g%nat)-drm)

la => g%alist%next
do i = 1,g%nat
  !take the sqrt of the average square velocity
  vel_med = dsqrt(factor* la%o%one_mass)
  do j = 1,dm
    call rang(r1)
    la%o%vel(j) = vel_med*r1+g%cm_vel(j)
  enddo
  la => la%next
enddo

call vel_changed()

! Para compenzar la pobre distribución
call set_scal_vel(g,newtemp)

! Para compenzar desplazamiento del centro de masa
call inq_cm_vel(g)
f=-g%cm_vel
call set_add_cmvel(g,f)

end subroutine set_gdist

subroutine set_scal_vel(g,newtemp)
class(group)               :: g
real(dp)                   :: factor,newtemp
class(atom_dclist),pointer :: la
integer                    :: i

call inq_temperature(g)

if (g%temp==0._dp) then
  call wwan('Not scale factor, actual temp is cero',newtemp/=0._dp)
  return
endif

factor=sqrt(newtemp/g%temp)

la => g%alist
do i = 1,g%nat
  la => la%next
  la%o%vel = la%o%vel * factor
enddo

call vel_changed()

g%temp = g%temp*(factor**2)
g%b_temp=.true.

end subroutine set_scal_vel

subroutine set_vel_from_file(g,archivo,frame)
! read atoms from file
use gems_elements,only:ncsym
integer                     :: i,j,u
character(*),intent(in)     :: archivo
type(group)                 :: g
class(atom_dclist),pointer  :: la
character(ncsym)            :: sym
integer,intent(in)          :: frame
character(3)                :: ext

i=len(trim(adjustl(archivo)))
ext=archivo(i-2:i)

u = find_io(30)
open(u,action='read',file=trim(adjustl(archivo)))

select case(ext)
case ('xyz')

  ! In case a specific frame wants to be loaded
  do i=1,frame-1
    read(u,*) j
    read(u,*)
    do j = 1, j
      read(u,*)
    enddo
  enddo

  ! Reading the frame
  read(u,*) j
  read(u,*)
  la => g%alist%next
  do i = 1,g%nat
    read(u,*) sym, la%o%vel
    la => la%next
  enddo

case default
  call werr(); write(logunit,*) 'File format ',trim(ext),' unknown'
end select

close(u)

call vel_changed()
end subroutine set_vel_from_file

!                                      posición y velocidad
!----------------------------------------------------------

subroutine set_pbc(g,pbc)
! rota un grupo un angulo r sobre el punto v
type(group),intent(inout)  :: g
class(atom_dclist),pointer  :: la
logical,intent(in)         :: pbc(dm)
integer                    :: i

la => g%alist
do i = 1,g%nat
  la => la%next
  la%o%pbc = pbc
enddo

end subroutine set_pbc

subroutine do_pbc_1(dclist,nat)
! Put the atoms inside the box using pbc
use gems_program_types, only: box
use gems_atoms, only: atom_dclist
class(atom_dclist),target  :: dclist
integer,intent(in)         :: nat
class(atom_dclist),pointer :: la
integer                    :: i,k

la => dclist
do i = 1,nat
  la => la%next

  do k=1,dm
    if (la%o%pbc(k)) then
      if(la%o%pos(k)>=box(k)) then
        la%o%pos(k)=la%o%pos(k)-box(k)
        la%o%boxcr(k)=la%o%boxcr(k)+1
      elseif(la%o%pos(k)<0.0_dp) then
        la%o%pos(k)=la%o%pos(k)+box(k)
        la%o%boxcr(k)=la%o%boxcr(k)-1
      endif
    endif
  enddo
enddo

end subroutine do_pbc_1

subroutine do_pbc_2(g)
class(group),intent(inout)  :: g

call do_pbc_1(g%alist,g%nat)

end subroutine do_pbc_2

   subroutine move(g,r)
   ! rota un grupo un angulo r sobre el punto v
   real(dp)                   :: r(dm)
   type(group)                :: g
   class(atom_dclist),pointer  :: la
   integer                    :: i

    la => g%alist%next
    do i = 1,g%nat
      la%o%pos = la%o%pos + r
      la => la%next
    enddo

    call pos_changed()

  end subroutine move

  subroutine rotate_axis_ref(g,tita,p,v)
   ! rotacion de un angulo tita usando un eje paralelo a algun eje cartesiano
   ! (p) que pasa por algun punto v.  El eje de rotacion viene dado por p:
   ! [1,2]=z; [1,3]=y; etc..
   real(dp)   ,intent(in)     :: tita
   integer    ,intent(in)     :: p(2)
   type(group),intent(in)     :: g
   real(dp)                   :: c,s,aux(2),t,v(dm)
   integer                    :: i
   class(atom_dclist),pointer  :: la

    ! supongo tita en grados
    call inq_pos_v(g,v)

    t=tita*pi/180.0_dp

    c=cos(t)
    s=sin(t)

    la => g%alist%next
    do i = 1,g%nat
      aux(1) = c*(la%o%pos_v(p(1))) + s*(la%o%pos_v(p(2)))
      aux(2) = c*(la%o%pos_v(p(2))) - s*(la%o%pos_v(p(1)))
      la%o%pos(p(1)) = aux(1)+v(p(1))
      la%o%pos(p(2)) = aux(2)+v(p(2))
      la => la%next
    enddo

    call pos_changed()

  end subroutine rotate_axis_ref

  subroutine rotate_rodrigues(g,tita,v,p)
   ! rotacion de un angulo tita (rad) usando el eje v que pasa por algun punto p.
   real(dp),intent(in)        :: p(dm),v(dm)
   real(dp),intent(in)        :: tita
   type(group),intent(in)     :: g
   real(dp)                   :: aux(dm)
   integer                    :: i
   class(atom_dclist),pointer  :: la

    la => g%alist%next
    do i = 1,g%nat
      aux = la%o%pos - p
      la%o%pos = rodrigues_rotation(tita,v,aux) + p
      la => la%next
    enddo

    call pos_changed()

  end subroutine rotate_rodrigues

  subroutine rotate(g,matrix,center)
   ! rotacion usando una matriz dada y un centro
   type(group),intent(in)     :: g
   real(dp),intent(in)        :: center(dm)
   real(dp)                   :: matrix(dm,dm)
   integer                    :: i
   class(atom_dclist),pointer  :: la

    call inq_pos_v(g,center)

    la => g%alist
    do i = 1,g%nat
      la => la%next
      la%o%pos = matmul(matrix,la%o%pos_v)+center
    enddo

    call pos_changed()

  end subroutine rotate

  subroutine givens_rotation(g,tita,p)
   real(dp)   ,intent(in)     :: tita
   integer    ,intent(in)     :: p(2)
   type(group),intent(in)     :: g
   real(dp)                   :: c,s,aux(2),t
   integer                    :: i
   class(atom_dclist),pointer  :: la

    ! supongo tita en grados

    t=tita*pi/180.0_dp

    c=cos(t)
    s=sin(t)

    la => g%alist%next
    do i = 1,g%nat
      aux(1) = c*(la%o%pos(p(1))) + s*(la%o%pos(p(2)))
      aux(2) = c*(la%o%pos(p(2))) - s*(la%o%pos(p(1)))
      la%o%pos(p(1)) = aux(1)
      la%o%pos(p(2)) = aux(2)
      la => la%next
    enddo

    call pos_changed()

  end subroutine givens_rotation

  subroutine align(g)
   ! orienta el objeto en los ejes principales de simetria
   type(group),intent(inout)  :: g
   real(dp)                   :: aux3(3,3) ,aux2(dm)

    call inq_principal_geometric_axis(g,aux3)
    if(dabs(inq_triaxial_param(g))<1.0e-8_dp) then
       call wwan('Ajustado por una esfera: No se puede orientar')
      return
    endif

    aux2=0.0_dp
    call rotate(g,aux3,aux2)

  end subroutine align

  subroutine align_mass(g)
   ! orienta el objeto en los ejes principales de simetria teniendo en cuenta la
   ! masa de cada particula
   type(group),intent(inout)  :: g
   real(dp)                   :: aux3(3,3) ,aux2(dm)

    call inq_principal_mass_axis(g,aux3)
    if(dabs(inq_triaxial_param(g))<1.0e-8_dp) then
       call wwan('Ajustado por una esfera: No se puede orientar')
       return
    endif

    aux2=0.0_dp
    call rotate(g,aux3,aux2)

  end subroutine align_mass

  subroutine alignxy_mass(g)
   ! orienta el objeto en los ejes principales de simetria teniendo en cuenta la
   ! masa de cada particula
   type(group),intent(inout)  :: g
   real(dp)                   :: aux3(3,3) ,aux2(dm)

    call inq_principal_mass_xyaxis(g,aux3)

    if(dabs(inq_triaxial_param(g))<1.0e-8_dp) then
       call wwan('Ajustado por una esfera: No se puede orientar')
       return
    endif

    aux2=0.0_dp
    call rotate(g,aux3,aux2)

  end subroutine alignxy_mass

  subroutine align_tutor(g,p,eje)
  ! orienta el objeto de manera que el punto p quede en el eje eje
   type(group),intent(inout)  :: g
   real(dp)                   :: aux(3),phi,cphi,aux2(3)
   real(dp),intent(in)        :: p(dm),eje(dm)

    ! Coseno del angulo de rotacion (p.dot.ejez)
    cphi = dot_product(p,eje)/sqrt(dot_product(p,p)*dot_product(eje,eje))
    cphi = dmin1(cphi,1.0_dp)
    cphi = dmax1(cphi,-1.0_dp)

    !definición del ángulo
    phi=dacos(cphi)
    if(phi<1.0e-8_dp) return


    if(dabs(phi-pi)<1e-8) then
      aux=ortogonal(eje)
    elseif(dabs(phi)<1e-8) then
      aux=ortogonal(eje)
    else
      ! Busco el eje de rotacion: (p)x(eje)
      aux(1)=p(2)*eje(3)-p(3)*eje(2)
      aux(2)=p(3)*eje(1)-p(1)*eje(3)
      aux(3)=p(1)*eje(2)-p(2)*eje(1)
    endif
    aux=aux/sqrt(dot_product(aux,aux))

    aux2=0.0_dp
    call rotate_rodrigues(g,phi,aux,aux2)

  end subroutine align_tutor

  subroutine bialign(g1,g2)
   ! orienta un grupo en funcion de los ejes principales de otro. Util para
   ! orientar sistemas esfericos como un octahedro truncado
   type(group),intent(inout)  :: g1,g2
   real(dp)                   :: aux3(3,3)

    call inq_principal_geometric_axis(g2,aux3)
    if(dabs(inq_triaxial_param(g2))<1.0e-8_dp) then
       call wwan('Ajustado por una esfera: No se puede orientar')
      return
    endif
    call rotate(g1,aux3,g1%cm_pos)

  end subroutine bialign

  subroutine bialignxy_mass(g1,g2)
   ! orienta un grupo en funcion de los ejes principales de otro. Util para
   ! orientar sistemas esfericos como un octahedro truncado
   type(group),intent(inout)  :: g1,g2
   real(dp)                   :: aux3(3,3),aux2(dm,dm)

    call inq_principal_mass_xyaxis(g2,aux3)
    if(dabs(inq_triaxial_param(g2))<1.0e-8_dp) then
       call wwan('Ajustado por una esfera: No se puede orientar')
      return
    endif
    aux2=0.0_dp
    call rotate(g1,aux3,aux2)

  end subroutine bialignxy_mass
!
!
!
!  subroutine rotate_2d(g,r,v)
!   ! rota un grupo un angulo r sobre el punto v
!   real(dp),intent(in)        :: v(2),r
!   real(dp)                   :: m(2,2)
!   type(group)                :: g
!   class(atom_dclist),pointer  :: la
!   integer                    :: i,j
!
!    m = rotmatrix_2x2(r)
!
!    call inq_pos_v(g,v)
!
!
!    la => g%alist%next
!    do i = 1,g%nat
!      la % o % pos_v = matmul( m, la % o % pos_v )
!      la%o%pos = la%o%pos_v + v
!      la => la%next
!    enddo
!
!    call pos_changed()
!
!  end subroutine rotate_2d
!
!  subroutine rotate_3d(g,r,v)
!   ! rota un grupo un angulo r sobre el punto v
!   real(dp),intent(in)        :: v(dm),r(:)
!   real(dp)                   :: m(dm,dm)
!   type(group)                :: g
!   class(atom_dclist),pointer  :: la
!   integer                    :: i,j
!
!    m = rotmatrix(r)
!
!    call inq_pos_v(g,v)
!
!
!    la => g%alist%next
!    do i = 1,g%nat
!      la % o % pos_v = matmul( m, la % o % pos_v )
!      la%o%pos = la%o%pos_v + v
!      la => la%next
!    enddo
!
!    call pos_changed()
!
!  end subroutine rotate_3d
!

  subroutine minterdist(g,rm)
   ! Normaliza las distancias interatomicas con la minima
   real(dp)                   :: rm,r(dm),vd(dm),rd,m
   type(group)                :: g
   class(atom_dclist),pointer  :: la,lb
   integer                    :: i,j

    la => g%alist%next
    m = 1e9_dp
    do i = 1,g%nat-1
      lb => la%next
      do j = i+1,g%nat
        vd = vdistance(la%o,lb%o, mic)
        rd = dot_product(vd,vd)
        if (rd<m) m = rd
        lb => lb%next
      enddo
      la => la%next
    enddo

    r(:)=rm/sqrt(m)

    call group_inq_cmpos(g)
    call expand(g,r,g%cm_pos)

  end subroutine

  subroutine expand(g,r,v)
   real(dp)                   :: v(dm),r(dm)
   type(group)                :: g
   class(atom_dclist),pointer  :: la
   integer                    :: i,j

    call inq_pos_v(g,v)

    la => g%alist%next
    do i = 1,g%nat
      do j = 1,dm
        la % o % pos_v(j) = la % o % pos_v(j) * r(j)
      enddo
      la%o%pos = la%o%pos_v + v
      la => la%next
    enddo

    call pos_changed()

  end subroutine expand

  subroutine set_sp(g,sp)
   type(group)                 :: g
   integer,intent(in)          :: sp
   integer                     :: i
   type(atom_dclist),pointer   :: la

    la => g%alist
    do i = 1,g%nat
      la => la%next
      la%o % sp = sp
    enddo

  end subroutine set_sp

  subroutine set_maxpos(g,r,k)
   real(dp)                    :: r
   type(group)                 :: g
   class(atom_dclist),pointer   :: la
   integer                     :: i,k

    call inq_boundingbox(g)

    r = r - g % maxpos(k)

    la => g%alist%next
    do i = 1,g%nat
      la % o % pos(k) = la % o % pos(k) + r
      la => la%next
    enddo

    call pos_changed()

  end subroutine set_maxpos

  subroutine set_minpos(g,r,k)
   real(dp)                    :: r
   type(group)                 :: g
   class(atom_dclist),pointer   :: la
   integer                     :: i,k

    call inq_boundingbox(g)

    r = r - g % minpos(k)

    la => g%alist%next
    do i = 1,g%nat
      la % o % pos(k) = la % o % pos(k) + r
      la => la%next
    enddo

    call pos_changed()

  end subroutine set_minpos

  subroutine set_cm_pos(g,v)
   real(dp)                    :: v(dm)
   type(group)                 :: g
   class(atom_dclist),pointer   :: la
   integer                     :: i

    call group_inq_cmpos(g)

    v = v - g % cm_pos
    la => g%alist%next
    do i = 1,g%nat
      la % o % pos = la % o % pos + v
      la => la%next
    enddo

    call pos_changed()

  end subroutine set_cm_pos

subroutine set_element(g,i1)
use gems_atoms, only: atom_dclist, atom_setelmnt
type(group)                 :: g
class(atom_dclist),pointer   :: la
integer                     :: i1,i

la => g%alist
do i = 1,g%nat
  la => la%next
  call atom_setelmnt(la%o,i1)
enddo

call mass_changed()

end subroutine set_element

  subroutine set_cg_pos(g,v)
   real(dp)                    :: v(dm)
   type(group)                 :: g
   class(atom_dclist),pointer   :: la
   integer                     :: i

    call inq_cg(g)

    v = v - g % cg_pos
    la => g%alist%next
    do i = 1,g%nat
      la % o % pos = la % o % pos + v
      la => la%next
    enddo

    call pos_changed()

  end subroutine set_cg_pos

  subroutine set_pos_from_file(g,archivo,frame)
    ! read atoms from file
    use gems_elements,only:ncsym
    integer                     :: i,j,u
    character(*),intent(in)     :: archivo
    type(group)                 :: g
    class(atom_dclist),pointer  :: la
    character(ncsym)            :: sym
    integer,intent(in)          :: frame
    character(3)                :: ext

    i=len(trim(adjustl(archivo)))
    ext=archivo(i-2:i)

    u = find_io(30)
    open(u,action='read',file=trim(adjustl(archivo)))

    select case(ext)
    case ('xyz')

      ! In case a specific frame wants to be loaded
      do i=1,frame-1
        read(u,*) j
        read(u,*)
        do j = 1, j
          read(u,*)
        enddo
      enddo

      ! Reading the frame
      read(u,*) j
      read(u,*)
      la => g%alist%next
      do i = 1,g%nat
        read(u,*) sym, la%o%pos
        la => la%next
      enddo

    case default
      call werr(); write(logunit,*) 'File format ',trim(ext),' unknown'
    end select

    close(u)

    call pos_changed()
  end subroutine set_pos_from_file

  subroutine set_clean_acel(g)
   type(group)                 :: g
   class(atom_dclist),pointer   :: la
   integer                     :: i

    la => g%alist%next
    do i = 1,g%nat
      la % o % acel  = 0.0_dp
      la % o % acel2 = 0.0_dp
      la % o % acel3 = 0.0_dp
      la % o % acel4 = 0.0_dp
      la => la%next
    enddo

  end subroutine set_clean_acel

  subroutine set_cm_vel(g,v)
   real(dp)                    :: v(dm)
   type(group)                 :: g
   class(atom_dclist),pointer   :: la
   integer                     :: i

    call inq_cm_vel(g)

    v = v - g % cm_vel
    la => g%alist%next
    do i = 1,g%nat
      la % o % vel = la % o % vel + v
      la => la%next
    enddo

    call vel_changed()

  end subroutine set_cm_vel

end module gems_set_properties

