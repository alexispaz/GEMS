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

 
module gems_fields
use gems_program_types
use gems_groups
use gems_atoms
use gems_constants
use gems_algebra 
use gems_neighbour, only:intergroup,intergroup0_empty

implicit none
save
private

public    :: fields_cli,fields_new,write_halfsho
real(dp)  :: numf=0, sumf=0

type,extends(intergroup) :: sho1d
  real(dp)  :: k
  contains
  procedure :: interact => sho1d_interact
  procedure,nopass :: cli => fields_cli
end type
              
type,extends(intergroup) :: sho2d
  real(dp)  :: kx, ky
  contains
  procedure :: interact => sho2d_interact
  procedure,nopass :: cli => fields_cli
end type
                            
type,extends(intergroup) :: shoplane
  real(dp)  :: k, r
  integer   :: j
  contains
  procedure :: interact => shoplane_interact
  procedure,nopass :: cli => fields_cli
end type
                             
type,extends(intergroup) :: hshoplane
  real(dp)  :: k, r
  integer   :: l
  contains
  procedure :: interact => hshoplane_interact
  procedure,nopass :: cli => fields_cli
end type
           
type,extends(intergroup) :: voter2d
  real(dp)  :: d1,d2,d3,d4
  contains
  procedure :: interact => voter2d_interact
  procedure,nopass :: cli => fields_cli
end type
                    
type,extends(intergroup) :: oscar2d
  contains
  procedure :: interact => oscar2d_interact
  procedure,nopass :: cli => fields_cli
end type
                   
type,extends(intergroup) :: leiva1d
  real(dp)  :: ca, cb
  real(dp)  :: alpha, beta
  contains
  procedure :: interact => leiva1d_interact
  procedure,nopass :: cli => fields_cli
end type
                 
type,extends(intergroup) :: pa1d
  real(dp),dimension(2)  :: a0,c,a1
  contains
  procedure :: interact => pa1d_interact
  procedure,nopass :: cli => fields_cli
end type
                    
type,extends(intergroup) :: sholine
  real(dp)  :: k
  real(dp)  :: v(dm)
  real(dp)  :: p(dm)
  contains
  procedure :: interact => sholine_interact
  procedure,nopass :: cli => fields_cli
end type
                    
type,extends(intergroup) :: lucas1d
  real(dp)     :: prof1,prof2,prof3
  real(dp)     :: anch1,anch2,anch3
  real(dp)     :: cent1,cent2,cent3
  contains
  procedure :: interact => lucas1d_interact
  procedure,nopass :: cli => fields_cli
end type
       

contains

!--------------------------------- interacción 


subroutine fields_new(ig,g1,w)
character(*),intent(in)   :: w
type(group)               :: g1
class(intergroup),pointer :: ig

select case(w)
case('halfsho_plane') 
  allocate(hshoplane::ig)
case('sho2d') 
  allocate(sho2d::ig)
case('voter2d')
  allocate(voter2d::ig)
case('oscar2d')
  allocate(oscar2d::ig)
case('leiva1d')
  allocate(leiva1d::ig)
case('pozoa1d')
  allocate(pa1d::ig)
case('sho_plane') 
  allocate(shoplane::ig)
case('sho_line') 
  allocate(sholine::ig)
case('lucas1d')
  allocate(lucas1d::ig)
case default
  call werr('Analitical potential not found')
endselect
   
call ig%init(g1=g1) 
ig%lista => intergroup0_empty

end subroutine fields_new

subroutine fields_cli(ig)
use gems_input_parsing
integer                   :: i1!,i2,i3,i4
real(dp)                  :: f1,f2,f3,f4,f5,f6,f7,f8,f9
real(dp)                  :: fv(dm),fv2(dm)
class(intergroup),intent(inout) :: ig
           
select type(ig)
type is(hshoplane) 
  call readf(f1,ev_ui)
  call readf(f2)
  call readi(i1)
  call werr("Dimension out of bound.",abs(i1)>dm)
  call hshoplane_set(ig,f1,f2,i1)
type is(sho1d) 
  call readf(f1,ev_ui)
  call sho1d_set(ig,f1) 
type is(sho2d) 
  call readf(f1,ev_ui)
  call readf(f2,ev_ui)
  call sho2d_set(ig,f1,f2)
type is(voter2d)
  call readf(f1)
  call readf(f2,ev_ui)
  call readf(f3,ev_ui)
  call readf(f4)
  call voter2d_set(ig,f1,f2,f3,f4)
type is(oscar2d)
type is(leiva1d)
  call read_energy(f2) ! Maximo izquierda
  call readf(f1,ev_ui) ! Kfuerza  izquierda
  call read_energy(f4) ! Maximo derecha  
  call readf(f3,ev_ui) ! Kfuerza  derecha

  f2=-sqrt(f2*ev_ui*exp(2.0_dp)/f1)
  f4=sqrt(f4*ev_ui*exp(2.0_dp)/f3)
  call leiva1d_set(ig,f1,f3,-2.0_dp/f2,2.0_dp/f4)
type is(pa1d)
  call readf(f1)
  call readf(f2,ev_ui)
  fv(:)=f2
  fv2(:)=f1
  call pa1d_set(ig,fv(1:2),fv2(1:2),fv2(1:2)*4._dp)
type is(shoplane) 
  call readf(f1,ev_ui)
  call readf(f2)
  call readi(i1)
  call werr("Dimension out of bound.",i1>dm)
  call shoplane_set(ig,f1,f2,i1)
type is(sholine) 
  call readf(f1,ev_ui)
  call readf(fv) ! axis_p
  call readf(fv2) ! axis_v
  call sholine_set(ig,f1,fv,fv2)
type is(lucas1d)
  ! Add new node to the internal intergroup list
  call readf(f1,kjm_ui)
  call readf(f2,kjm_ui)
  call readf(f3,kjm_ui)
  call readf(f4)
  call readf(f5)
  call readf(f6)
  call readf(f7)
  call readf(f8)
  call readf(f9)
  call lucas1d_set(ig,f1,f2,f3,f4,f5,f6,f7,f8,f9)
class default
  call werr('Analitical potential not found')
endselect

                           
end subroutine
                             
      
subroutine pa1d_set(ig,a,a0,c)
type(pa1d)     :: ig
real(dp)       :: a(2),a0(2),c(2)
ig%a0(:)=a0(:)
ig%a1(:)=a(:)
ig%c(:)=c(:)
end subroutine
            
subroutine pa1d_interact(ig)
!Pozo armonico de una dimension
class(pa1d),intent(inout)   :: ig 
integer                    :: i,m
real(dp)                   :: x,aux,f,e
type(atom_dclist),pointer  :: la
                
ig%epot=0._dp

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next
                      
  m=1
  if(la%o%pos(1)>0) m=2

  x = la%o%pos(1)
  aux = ig%a1(m)*exp(2.0_dp*(1.0_dp-x/ig%a0(m)))

  e = aux*(x/ig%a0(m))**2
  f = (1-x/ig%a0(m)) * 2*aux*x/(ig%a0(m)*ig%a0(m))

  if(abs(x)>abs(ig%c(m)))then
    f = f*((x-ig%c(m))**4+1)+e*4.0_dp*(x-ig%c(m))**3
    e = e*((x-ig%c(m))**4+1)
  endif

  la%o%epot = la%o%epot+e 
  la%o%force = la%o%force-f 
  ig%epot = ig%epot + e

enddo

end subroutine
     
subroutine voter2d_set(ig,d1,d2,d3,d4)
type(voter2d)  :: ig
real(dp)       :: d1,d2,d3,d4
ig%d1=d1
ig%d2=d2
ig%d3=d3
ig%d4=d4
end subroutine
      
subroutine voter2d_interact(ig)
!"A method for accelerating the molecular dynamics simulation of infrequent
!events". Arthur F. Voter. J. Chem- Phys 106 (11) 1997
class(voter2d),intent(inout) :: ig 
real(dp)                    :: x,y,e
type(atom_dclist),pointer   :: la
integer                     :: i
                
ig%epot=0._dp

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next
    

  x=la%o%pos(1)
  y=la%o%pos(2)
  e = cos(twopi*x)*(1+ig%d1*y)*ev_ui+2.0_dp*ig%d2*(pi*y)**2+ig%d3*cos(twopi*x/ig%d4)
  e = e +1.203*ev_ui !corrigo para que sea negativo (al menos en la parte donde hay boost)
  la%o%force(1) = twopi*(1+ig%d1*y)*sin(twopi*x)*ev_ui+ig%d3*twopi/ig%d4*sin(twopi/ig%d4*x)
  la%o%force(2) =-ig%d1*cos(twopi*x)*ev_ui -4*pi**2*ig%d2*y
  !la%o%epot     = 0.5_dp*la%o%epot     
  !la%o%force(1) = 0.5_dp*la%o%force(1) 
  !la%o%force(2) = 0.5_dp*la%o%force(2) 

  la%o%epot = la%o%epot+e
  ig%epot = ig%epot + e
  
enddo
end subroutine
        
subroutine leiva1d_set(ig,a,b,alpha,beta)
type(leiva1d)  :: ig
real(dp)       :: a,b,alpha,beta
ig%alpha=alpha
ig%cb=b
ig%beta=beta
ig%ca=a
end subroutine
      
subroutine leiva1d_interact(ig)
class(leiva1d),intent(inout) :: ig 
integer                     :: i
real(dp)                    :: x,e
type(atom_dclist),pointer   :: la
               
ig%epot=0._dp

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next
    
  x=la%o%pos(1)
  if (x>0) then
    e = ig%cb*x**2*exp(-ig%beta*x)
    la%o%force(1) = -ig%cb*x*exp(-ig%beta*x)*(2-x*ig%beta)
  else
    e = ig%ca*x**2*exp(ig%alpha*x)
    la%o%force(1) = -ig%ca*x*exp(ig%alpha*x)*(2+x*ig%alpha)
  endif
  la%o%epot =  la%o%epot + e
  ig%epot = ig%epot + e
enddo
end subroutine
            
subroutine oscar2d_interact(ig)
! Sacada del examen del curso de metodos computacionales 2011
class(oscar2d),intent(inout)   :: ig 
integer                       :: i
real(dp)                      :: g,f,    &
                                 h,h1,h2,&
                                 v,v1,v2,&
                                 x,y,b1,b2
real(dp)                      :: e
type(atom_dclist),pointer     :: la
               
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
      
subroutine sho2d_set(ig,kx,ky)
real(dp),intent(in) :: kx,ky
type(sho2d)         :: ig
ig%kx=kx
ig%ky=ky
end subroutine
                      
subroutine sho2d_interact(ig)
class(sho2d),intent(inout) :: ig 
integer                   :: i
real(dp)                  :: x,y,e
type(atom_dclist),pointer :: la
               
ig%epot=0._dp

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next
    
  x=la%o%pos(1) 
  y=la%o%pos(2)
  la%o%force(1) = -ig%kx*x
  la%o%force(2) = -ig%ky*y
  e = (ig%kx*x**2+ig%ky*y**2)*0.5_dp
  !e = e -2*ev_ui !corrigo para que sea negativo
  la%o%epot = la%o%epot + e
  ig%epot = ig%epot + e

enddo

end subroutine

subroutine sho1d_set(ig,p1)
real(dp),intent(in)   :: p1
type(sho1d)           :: ig
ig%k=p1
end subroutine

subroutine sho1d_interact(ig)
class(sho1d),intent(inout)   :: ig 
integer                     :: i
real(dp)                    :: x,e
type(atom_dclist),pointer   :: la
                 
ig%epot=0._dp

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next
    
  x=la%o%pos(1) 
  la%o%force(1) = -ig%k*x

  e= ig%k*x**2*0.5_dp

  la%o%epot = la%o%epot + e
  ig%epot = ig%epot + e

enddo
end subroutine sho1d_interact
   


subroutine sholine_set(ig,k,p,v)
real(dp),intent(in)  :: k
type(sholine)        :: ig
real(dp),intent(in)  :: p(dm)
real(dp),intent(in)  :: v(dm)

! parametros
ig%k=k
ig%p(:)=p(:)
ig%v(:)=v(:)
ig%v(:)=ig%v(:)/sqrt(dot_product(ig%v(:),ig%v(:)))

end subroutine sholine_set

subroutine sholine_interact(ig)
class(sholine),intent(inout)   :: ig 
type(atom_dclist),pointer         :: la
real(dp)                          :: x(dm),e
integer                           :: i
                 
ig%epot=0._dp

la => ig%a
do i = 1,ig%n(1)
  la=>la%next
    
  ! Calcular el vector desplazamiento al eje
  x(:)=la%o%pos(1:dm)-ig%p(1:dm)
  x(:)=x(:)-dot_product(x,ig%v)*ig%v(:)

  la%o%force(:) = la%o%force(:)-ig%k*x(:)
  e=ig%k*dot_product(x,x)*0.5_dp

  la%o%epot = la%o%epot + e
  ig%epot = ig%epot + e

enddo

end subroutine sholine_interact


                

subroutine hshoplane_set(ig,k,r,l)
type(hshoplane)      :: ig
real(dp),intent(in)  :: k, r
integer,intent(in)   :: l
ig%k=k
ig%r=r
ig%l=l
end subroutine hshoplane_set
                     
subroutine hshoplane_interact(ig)
class(hshoplane),intent(inout)   :: ig 
type(atom_dclist),pointer       :: la
real(dp)                        :: x,e
integer                         :: i,j
                 
j=abs(ig%l)
      
ig%epot=0._dp

la => ig%a
do i = 1,ig%n(1)
  la=>la%next
    
  ! Calcular el vector desplazamiento al eje
  x=la%o%pos(j)-ig%r

  ! Choose half of the harmonic potential
  if (x*ig%l<0) cycle

  sumf=sumf-ig%k*x
  numf=numf+1

  la%o%force(j) = la%o%force(j)-ig%k*x
  e=ig%k*x**2*0.5_dp

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


subroutine shoplane_set(ig,k,r,j)
type(shoplane)      :: ig
real(dp),intent(in) :: k,r
integer,intent(in)  :: j
ig%k=k
ig%r=r
ig%j=j
end subroutine

subroutine shoplane_interact(ig)
class(shoplane),intent(inout)   :: ig 
type(atom_dclist),pointer         :: la
real(dp)                          :: x,e
integer                           :: i,j
                 
j=ig%j
      
ig%epot=0._dp

la => ig%a
do i = 1,ig%n(1)
  la=>la%next
    
  ! Calcular el vector desplazamiento al eje
  x=la%o%pos(j)-ig%r

  la%o%force(j) = la%o%force(j)-ig%k*x
  e=ig%k*x**2*0.5_dp

  la%o%epot = la%o%epot + e
  ig%epot = ig%epot + e

enddo

end subroutine
           
 
subroutine lucas1d_set(ig,a1,a2,a3,b1,b2,b3,c1,c2,c3)
real(dp),intent(in)  :: a1,a2,a3,b1,b2,b3,c1,c2,c3
type(lucas1d)        :: ig
ig%prof1=a1
ig%prof2=a2
ig%prof3=a3
ig%anch1=b1
ig%anch2=b2
ig%anch3=b3
ig%cent1=c1
ig%cent2=c2
ig%cent3=c3
end subroutine
                  
subroutine lucas1d_interact(ig)
class(lucas1d),intent(inout) :: ig 
real(dp)                    :: x,e
real(dp)                    :: a1,a2,a3,b1,b2,b3,c1,c2,c3
type(atom_dclist),pointer   :: la
integer                     :: i
               
ig%epot=0._dp

! parameters     
a1=ig%prof1 
a2=ig%prof2 
a3=ig%prof3 
                 
b1=ig%anch1 
b2=ig%anch2 
b3=ig%anch3 
                 
c1=ig%cent1 
c2=ig%cent2 
c3=ig%cent3 

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next
    
  x=la%o%pos(1) 
  la%o%force(1) = la%o%force(1) &
                - 2.0_dp*a1*b1*(x-c1)*exp(-b1*(x-c1)*(x-c1)) &
                - 2.0_dp*a2*b2*(x-c2)*exp(-b2*(x-c2)*(x-c2)) &
                - 2.0_dp*a3*b3*(x-c3)*exp(-b3*(x-c3)*(x-c3)) 

  e = -a1*exp(-b1*(x-c1)*(x-c1)) &
      -a2*exp(-b2*(x-c2)*(x-c2)) &
      -a3*exp(-b3*(x-c3)*(x-c3)) 

  !e = e -2*ev_ui !corrijo para que sea negativo
  la%o%epot = la%o%epot + e
  ig%epot = ig%epot + e

enddo

end subroutine


end module gems_fields
