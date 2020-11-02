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

 
module gems_analitic_pot
use gems_program_types
use gems_groups
use gems_atoms
use gems_constants
use gems_algebra 
use gems_neighbour, only:intergroup,intergroup0_empty

implicit none
save
private

real(dp)     :: pa1d_a(2),pa1d_a0(2),pa1d_c(2)
real(dp)     :: vt_d1,vt_d2,vt_d3,vt_d4
real(dp)     :: le_a,le_b,le_alpha,le_beta
real(dp)     :: sho_k
real(dp)     :: prof1,prof2,prof3,anch1,anch2,anch3,cent1,cent2,cent3

public :: sho_line_set
public :: sho_plane_init,sho_plane_set
public :: analitical_new,analitical_cli,write_halfsho
real(dp)  :: numf=0, sumf=0

public :: lucas1d_set

contains

!--------------------------------- interacción 

subroutine pa1d(ig)
!Pozo armonico de una dimension
class(intergroup),intent(inout)   :: ig 
integer                           :: i,m
real(dp)                          :: x,aux,f,e
type(atom_dclist),pointer         :: la
                
ig%epot=0._dp

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next
                      
  m=1
  if(la%o%pos(1)>0) m=2

  x = la%o%pos(1)
  aux = pa1d_a(m)*exp(2.0_dp*(1.0_dp-x/pa1d_a0(m)))

  e = aux*(x/pa1d_a0(m))**2
  f = (1-x/pa1d_a0(m)) * 2*aux*x/(pa1d_a0(m)*pa1d_a0(m))

  if(abs(x)>abs(pa1d_c(m)))then
    f = f*((x-pa1d_c(m))**4+1)+e*4.0_dp*(x-pa1d_c(m))**3
    e = e*((x-pa1d_c(m))**4+1)
  endif

  la%o%epot = la%o%epot+e 
  la%o%force = la%o%force-f 
  ig%epot = ig%epot + e

enddo

end subroutine

subroutine voter2d(ig)
!"A method for accelerating the molecular dynamics simulation of infrequent
!events". Arthur F. Voter. J. Chem- Phys 106 (11) 1997
class(intergroup),intent(inout)   :: ig 
integer                   :: i
real(dp)                  :: x,y,e
type(atom_dclist),pointer :: la
                
ig%epot=0._dp

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next
    

  x=la%o%pos(1)
  y=la%o%pos(2)
  e = cos(twopi*x)*(1+vt_d1*y)*ev_ui+2.0_dp*vt_d2*(pi*y)**2+vt_d3*cos(twopi*x/vt_d4)
  e = e +1.203*ev_ui !corrigo para que sea negativo (al menos en la parte donde hay boost)
  la%o%force(1) = twopi*(1+vt_d1*y)*sin(twopi*x)*ev_ui+vt_d3*twopi/vt_d4*sin(twopi/vt_d4*x)
  la%o%force(2) =-vt_d1*cos(twopi*x)*ev_ui -4*pi**2*vt_d2*y
  !la%o%epot     = 0.5_dp*la%o%epot     
  !la%o%force(1) = 0.5_dp*la%o%force(1) 
  !la%o%force(2) = 0.5_dp*la%o%force(2) 

  la%o%epot = la%o%epot+e
  ig%epot = ig%epot + e
  
enddo
end subroutine

subroutine leiva1d(ig)
!"A method for accelerating the molecular dynamics simulation of infrequent
!events". Arthur F. Voter. J. Chem- Phys 106 (11) 1997
class(intergroup),intent(inout) :: ig 
integer                         :: i
real(dp)                        :: x,e
type(atom_dclist),pointer       :: la
               
ig%epot=0._dp

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next
    
  x=la%o%pos(1)
  if (x>0) then
    e = le_B*x**2*exp(-le_beta*x)
    la%o % force(1) = -le_B*x*exp(-le_beta*x)*(2-x*le_beta)
  else
    e = le_A*x**2*exp(le_alpha*x)
    la%o % force(1) = -le_A*x*exp(le_alpha*x)*(2+x*le_alpha)
  endif
  la%o%epot =  la%o%epot + e
  ig%epot = ig%epot + e
enddo
end subroutine

subroutine oscar2d(ig)
! Sacada del examen del curso de metodos computacionales 2011
class(intergroup),intent(inout)   :: ig 
integer                 :: i
real(dp)                :: g,f,    &
                           h,h1,h2,&
                           v,v1,v2,&
                           x,y,b1,b2
real(dp)                :: e
type(atom_dclist),pointer :: la
               
ig%epot=0._dp

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next
    
  x=la%o%pos(1)
  y=la%o%pos(2)
  g = -3.0_dp*exp(-(x**2))
  h1 = exp(-((y-(5.0_dp/3.0_dp))**2))
  h2 = -exp(-((y-(1.0_dp/3.0_dp))**2))
  f = -5.0_dp*exp(-(y**2))
  v1 = exp(-((x-1.33_dp)**2))
  v2 = exp(-((x+1.0_dp)**2))
  b1 = cos(2.0_dp*pi*x)*(1.0_dp+0.2_dp*y)
  b2 = (1.0e-3_dp*pi*y**2)


  h=h1+h2 
  v=v1+v2
  e = g*h + f*v + b1 + b2
  la%o%force(1) = -2.0_dp*x*g*h-f*2.0_dp*((x-1.33_dp)*v1-(x+1.0_dp)*v2) &
                  -2.0_dp*pi*sin(2.0_dp*pi*x)*(1.0_dp+0.2_dp*y)
  la%o%force(2) = -g*2.0_dp*((y-(5.0_dp/3.0_dp))*h1+(y-(1.0_dp/3.0_dp)*h2))-2.0_dp*y*f*v &
                  + 0.2_dp*cos(2.0_dp*pi*x) - 2.0_dp*pi*3.0_dp*y

  la%o%epot = la%o%epot + e
  ig%epot = ig%epot + e
  
enddo 

end subroutine
              
subroutine sho2d(ig)
class(intergroup),intent(inout)   :: ig 
integer                 :: i
real(dp)                :: x,y,e
type(atom_dclist),pointer :: la
               
ig%epot=0._dp

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next
    
  x=la%o%pos(1) 
  y=la%o%pos(2)
  la%o%force(1) = -sho_k*x
  la%o%force(2) = -sho_k*y
  e = sho_k*(x**2+y**2)*0.5_dp
  !e = e -2*ev_ui !corrigo para que sea negativo
  la%o%epot = la%o%epot + e
  ig%epot = ig%epot + e

enddo

end subroutine

subroutine sho1d_set(g1,p1)
  type(group),intent(in)     :: g1
  real(dp),intent(in)        :: p1
  type(intergroup),pointer   :: igr

  allocate(igr)
  call igr%init(g1=g1)
  igr%lista => intergroup0_empty
  igr%interact => sho1d
  
  ! parameters
  sho_k=p1
      
end subroutine

subroutine sho1d(ig)
class(intergroup),intent(inout)   :: ig 
integer                   :: i
real(dp)                  :: x,e
type(atom_dclist),pointer :: la
                 
ig%epot=0._dp

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next
    
  x=la%o%pos(1) 
  la%o%force(1) = -sho_k*x

  e= sho_k*x**2*0.5_dp

  la%o%epot = la%o%epot + e
  ig%epot = ig%epot + e

enddo
end subroutine
   


subroutine sho_line_set(g1,p1,p,v)
type(group),intent(in)     :: g1
real(dp),intent(in)        :: p1
type(intergroup),pointer   :: igr
real(dp),intent(in)        :: p(dm)
real(dp)                   :: v(dm)
integer                    :: i,j


! Initialize the integroup
allocate(igr)
call igr%init(g1=g1)
igr%lista => intergroup0_empty
igr%interact => sho_line ! Deberia elegir aca

! parametros
j=1
call igr%p%put(j,p1)
v(:)=v(:)/sqrt(dot_product(v,v))
do i=1,dm
  j=j+1
  call igr%p%put(j,v(i))
enddo  
do i=1,dm
  j=j+1
  call igr%p%put(j,p(i))
enddo  

end subroutine

subroutine sho_line(ig)
class(intergroup),intent(inout)   :: ig 
type(atom_dclist),pointer         :: la
real(dp)                          :: x(dm),v(dm),p(dm),k,e
integer                           :: i
                 
k=ig%p%o(1)
v(1:dm)=ig%p%o(2:dm+1)
p(1:dm)=ig%p%o(dm+2:2*dm+1)
      
ig%epot=0._dp

la => ig%a
do i = 1,ig%n(1)
  la=>la%next
    
  ! Calcular el vector desplazamiento al eje
  x(:)=la%o%pos(1:dm)-p(1:dm)
  x(:)=x(:)-dot_product(x,v)*v(:)

  la%o%force(:) = la%o%force(:)-k*x(:)
  e=k*dot_product(x,x)*0.5_dp

  la%o%epot = la%o%epot + e
  ig%epot = ig%epot + e

enddo
end subroutine


           
function analitical_new(w,g1) result(igr)
! TODO: Generalize for any plane
type(group),intent(in)     :: g1
type(intergroup),pointer   :: igr
character(*),intent(in)    :: w

! Bulid the igr
allocate(igr)
call igr%init(g1=g1)
          
! No list
igr%lista => intergroup0_empty

select case(w)
case('sho')
  igr%interact => halfsho_plane
case('sho2d')
  igr%interact => sho2d
case('voter2d')
  igr%interact => voter2d
case('oscar2d')
  igr%interact => oscar2d
case('leiva1d')
  igr%interact => leiva1d
case('pozoa1d')
  igr%interact => pa1d
case('halfsho_plane')
  igr%interact => halfsho_plane
case default
  call werr('Analitical interaction not found')
endselect
 
end function

subroutine analitical_cli(igr,w)
use gems_input_parsing
type(intergroup)          :: igr
character(*),intent(in)   :: w
integer                   :: i1!,i2,i3,i4
real(dp)                  :: f1,f2,f3,f4
           
select case(w)
case('halfsho_plane') 

  call readf(f1,ev_ui)
  call readf(f2)
  call readi(i1)

  call werr("Dimension out of bound.",abs(i1)>dm)

  call igr%p%put(1,f1)
  call igr%p%put(2,f2)
  call igr%i%put(1,i1)
            
case('sho2d') 
  call readf(f1,ev_ui)
  sho_k=f1
  !call sho2d_set(g1,p1)
case('voter2d')
  call readf(f1)
  call readf(f2,ev_ui)
  call readf(f3,ev_ui)
  call readf(f4)
  vt_d1=f1
  vt_d2=f2
  vt_d3=f3
  vt_d4=f4
       
case('oscar2d')

case('leiva1d')
   
  call read_energy(f2) ! Maximo izquierda
  call readf(f1,ev_ui) ! Kfuerza  izquierda
  call read_energy(f4) ! Maximo derecha  
  call readf(f3,ev_ui) ! Kfuerza  derecha
  
  le_a=f1
  le_alpha=-sqrt(f2*ev_ui*exp(2.0_dp)/le_a)
  le_alpha=-2.0_dp/le_alpha
  le_b=f3
  le_beta=sqrt(f4*ev_ui*exp(2.0_dp)/le_b)
  le_beta=2.0_dp/le_beta
            
case('pozoa1d')
 
  call readf(f1)
  call readf(f2,ev_ui)
  
  pa1d_a0(:) =f1
  pa1d_a(:)  =f2
  pa1d_c(:)  =pa1d_a0*4
  
case default
  call werr('SHO type not found')
endselect

                           
end subroutine
                             
   
subroutine halfsho_plane(ig)
class(intergroup),intent(inout)   :: ig 
type(atom_dclist),pointer         :: la
real(dp)                          :: x,r,k,e
integer                           :: i,j,l
                 
k=ig%p%o(1)
r=ig%p%o(2)
l=ig%i%o(1)
j=abs(l)
      
ig%epot=0._dp

la => ig%a
do i = 1,ig%n(1)
  la=>la%next
    
  ! Calcular el vector desplazamiento al eje
  x=la%o%pos(j)-r

  ! Choose half of the harmonic potential
  if (x*l<0) cycle

  sumf=sumf-k*x
  numf=numf+1

  la%o%force(j) = la%o%force(j)-k*x
  e=k*x**2*0.5_dp

  la%o%epot = la%o%epot + e
  ig%epot = ig%epot + e

enddo

end subroutine
 
subroutine write_halfsho(op)
  use gems_output
  class(outpropa)   :: op
  op%f(1)  = op%f(1)  + sumf
  op%f(2)  = op%f(2)  + numf
  sumf=0
  numf=0
end subroutine write_halfsho



           
subroutine sho_plane_init(g1,p1,p2,k)
! TODO: Generalize for larger planes  
type(group),intent(in)     :: g1
real(dp),intent(in)        :: p1,p2
type(intergroup),pointer   :: igr
integer                    :: k

! Initialize the integroup
allocate(igr)
call igr%init(g1=g1)
igr%lista => intergroup0_empty
igr%interact => sho_plane ! Deberia elegir aca

! parametros
call sho_plane_set(igr,p1,p2,k)

end subroutine
            
subroutine sho_plane_set(igr,p1,p2,k)
class(intergroup)          :: igr
real(dp),intent(in)        :: p1,p2
integer                    :: k

! parametros
call igr%p%put(1,p1)
call igr%p%put(2,p2)
call igr%i%put(1,k)

end subroutine

subroutine sho_plane(ig)
class(intergroup),intent(inout)   :: ig 
type(atom_dclist),pointer         :: la
real(dp)                          :: x,r,k,e
integer                           :: i,j
                 
k=ig%p%o(1)
r=ig%p%o(2)
j=ig%i%o(1)
      
ig%epot=0._dp

la => ig%a
do i = 1,ig%n(1)
  la=>la%next
    
  ! Calcular el vector desplazamiento al eje
  x=la%o%pos(j)-r

  la%o%force(j) = la%o%force(j)-k*x
  e=k*x**2*0.5_dp

  la%o%epot = la%o%epot + e
  ig%epot = ig%epot + e

enddo

end subroutine
           
 
subroutine lucas1d_set(g1,a1,a2,a3,b1,b2,b3,c1,c2,c3)
  type(group),intent(in)     :: g1
  real(dp),intent(in)        :: a1,a2,a3,b1,b2,b3,c1,c2,c3
  type(intergroup),pointer   :: igr

  ! Add new node to the internal intergroup list
  allocate(igr)
  call igr%init(g1=g1)
  igr%lista => intergroup0_empty
  igr%interact => lucas1d
  
  ! parameters
  prof1=a1
  prof2=a2
  prof3=a3
  anch1=b1
  anch2=b2
  anch3=b3
  cent1=c1
  cent2=c2
  cent3=c3
      
end subroutine
                  
subroutine lucas1d(ig)
class(intergroup),intent(inout)   :: ig 
integer                 :: i
real(dp)                :: x,e,w
type(atom_dclist),pointer :: la
               
ig%epot=0._dp

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next
    
  x=la%o%pos(1) 
  la%o%force(1) = la%o%force(1) -2.0_dp*prof1*anch1*(x-cent1)*exp(-anch1*(x-cent1)*(x-cent1)) &
                  -2.0_dp*prof2*anch2*(x-cent2)*exp(-anch2*(x-cent2)*(x-cent2)) &
                  -2.0_dp*prof3*anch3*(x-cent3)*exp(-anch3*(x-cent3)*(x-cent3)) 

  e =  -prof1*exp(-anch1*(x-cent1)*(x-cent1))-prof2*exp(-anch2*(x-cent2)*(x-cent2))-prof3*exp(-anch3*(x-cent3)*(x-cent3)) 
  !e = e -2*ev_ui !corrigo para que sea negativo
  la%o%epot = la%o%epot + e
  ig%epot = ig%epot + e
  w=-2.0_dp*prof1*anch1*(x-cent1)*exp(-anch1*(x-cent1)*(x-cent1)) &
                  -2.0_dp*prof2*anch2*(x-cent2)*exp(-anch2*(x-cent2)*(x-cent2)) &
                  -2.0_dp*prof3*anch3*(x-cent3)*exp(-anch3*(x-cent3)*(x-cent3))

enddo

end subroutine


end module gems_analitic_pot
