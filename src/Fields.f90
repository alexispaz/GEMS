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
use gems_constants
use gems_algebra 
use gems_neighbor, only:ngroup

implicit none
save
private

public    :: field_cli,field_new,write_halfsho
real(dp)  :: numf=0, sumf=0

! This is like a sho0d
! type,extends(ngroup) :: shofix
!   real(dp)              :: k
!   real(dp),allocatable  :: r(:,:)
!   contains
!   procedure :: interact => shofix_interact
!   procedure,nopass :: cli => pair_cli
! end type
                         
type,extends(ngroup) :: sho1d
  real(dp)  :: k
  contains
  procedure :: interact => sho1d_interact
  procedure,nopass :: cli => field_cli
end type
              
type,extends(ngroup) :: sho2d
  real(dp)  :: kx, ky
  contains
  procedure :: interact => sho2d_interact
  procedure,nopass :: cli => field_cli
end type
                            
type,extends(ngroup) :: shoplane
  real(dp)  :: k, r
  integer   :: j
  contains
  procedure :: interact => shoplane_interact
  procedure,nopass :: cli => field_cli
end type
                             
type,extends(ngroup) :: hshoplane
  real(dp)  :: k, r
  integer   :: l
  contains
  procedure :: interact => hshoplane_interact
  procedure,nopass :: cli => field_cli
end type
           
type,extends(ngroup) :: voter2d
  real(dp)  :: d1,d2,d3,d4
  contains
  procedure :: interact => voter2d_interact
  procedure,nopass :: cli => field_cli
end type
                    
type,extends(ngroup) :: oscar2d
  contains
  procedure :: interact => oscar2d_interact
  procedure,nopass :: cli => field_cli
end type
                   
type,extends(ngroup) :: leiva1d
  real(dp)  :: ca, cb
  real(dp)  :: alpha, beta
  contains
  procedure :: interact => leiva1d_interact
  procedure,nopass :: cli => field_cli
end type
                 
type,extends(ngroup) :: pa1d
  real(dp),dimension(2)  :: a0,c,a1
  contains
  procedure :: interact => pa1d_interact
  procedure,nopass :: cli => field_cli
end type
                    
type,extends(ngroup) :: sholine
  real(dp)  :: k
  real(dp)  :: v(dm)
  real(dp)  :: p(dm)
  contains
  procedure :: interact => sholine_interact
  procedure,nopass :: cli => field_cli
end type
                    
type,extends(ngroup) :: lucas1d
  real(dp)     :: prof1,prof2,prof3
  real(dp)     :: anch1,anch2,anch3
  real(dp)     :: cent1,cent2,cent3
  contains
  procedure :: interact => lucas1d_interact
  procedure,nopass :: cli => field_cli
end type
       

contains

!--------------------------------- interacción 


subroutine field_new(g,w)
character(*),intent(in)   :: w
class(ngroup),pointer :: g

select case(w)
case('halfsho_plane') 
  allocate(hshoplane::g)
case('sho2d') 
  allocate(sho2d::g)
case('voter2d')
  allocate(voter2d::g)
case('oscar2d')
  allocate(oscar2d::g)
case('leiva1d')
  allocate(leiva1d::g)
case('pozoa1d')
  allocate(pa1d::g)
case('sho_plane') 
  allocate(shoplane::g)
case('sho_line') 
  allocate(sholine::g)
case('lucas1d')
  allocate(lucas1d::g)
case default
  call werr('Field type not found')
endselect
g%lista => null()

end subroutine field_new

subroutine field_cli(g)
use gems_input_parsing
integer                   :: i1!,i2,i3,i4
real(dp)                  :: f1,f2,f3,f4,f5,f6,f7,f8,f9
real(dp)                  :: fv(dm),fv2(dm)
class(ngroup),intent(inout) :: g
           
select type(g)
type is(hshoplane) 
  call readf(f1,ev_ui)
  call readf(f2)
  call readi(i1)
  call werr("Dimension out of bound.",abs(i1)>dm)
  call hshoplane_set(g,f1,f2,i1)
type is(sho1d) 
  call readf(f1,ev_ui)
  call sho1d_set(g,f1) 
type is(sho2d) 
  call readf(f1,ev_ui)
  call readf(f2,ev_ui)
  call sho2d_set(g,f1,f2)
type is(voter2d)
  call readf(f1)
  call readf(f2,ev_ui)
  call readf(f3,ev_ui)
  call readf(f4)
  call voter2d_set(g,f1,f2,f3,f4)
type is(oscar2d)
type is(leiva1d)
  call read_energy(f2) ! Maximo izquierda
  call readf(f1,ev_ui) ! Kfuerza  izquierda
  call read_energy(f4) ! Maximo derecha  
  call readf(f3,ev_ui) ! Kfuerza  derecha

  f2=-sqrt(f2*ev_ui*exp(2.0_dp)/f1)
  f4=sqrt(f4*ev_ui*exp(2.0_dp)/f3)
  call leiva1d_set(g,f1,f3,-2.0_dp/f2,2.0_dp/f4)
type is(pa1d)
  call readf(f1)
  call readf(f2,ev_ui)
  fv(:)=f2
  fv2(:)=f1
  call pa1d_set(g,fv(1:2),fv2(1:2),fv2(1:2)*4._dp)
type is(shoplane) 
  call readf(f1,ev_ui)
  call readf(f2)
  call readi(i1)
  call werr("Dimension out of bound.",i1>dm)
  call shoplane_set(g,f1,f2,i1)
type is(sholine) 
  call readf(f1,ev_ui)
  call readf(fv) ! axis_p
  call readf(fv2) ! axis_v
  call sholine_set(g,f1,fv,fv2)
type is(lucas1d)
  ! Add new node to the internal ngroup list
  call readf(f1,kjm_ui)
  call readf(f2,kjm_ui)
  call readf(f3,kjm_ui)
  call readf(f4)
  call readf(f5)
  call readf(f6)
  call readf(f7)
  call readf(f8)
  call readf(f9)
  call lucas1d_set(g,f1,f2,f3,f4,f5,f6,f7,f8,f9)
class default
  call werr('Analitical potential not found')
endselect

                           
end subroutine
                             
      
subroutine pa1d_set(g,a,a0,c)
type(pa1d)     :: g
real(dp)       :: a(2),a0(2),c(2)
g%a0(:)=a0(:)
g%a1(:)=a(:)
g%c(:)=c(:)
end subroutine
            
subroutine pa1d_interact(g)
!Pozo armonico de una dimension
class(pa1d),intent(inout)   :: g 
integer                    :: i,m
real(dp)                   :: x,aux,f,e
type(atom_dclist),pointer  :: la
                
g%epot=0._dp

la => g%alist      ! sobre los atomos
do i = 1,g%nat
  la=>la%next
                      
  m=1
  if(la%o%pos(1)>0) m=2

  x = la%o%pos(1)
  aux = g%a1(m)*exp(2.0_dp*(1.0_dp-x/g%a0(m)))

  e = aux*(x/g%a0(m))**2
  f = (1-x/g%a0(m)) * 2*aux*x/(g%a0(m)*g%a0(m))

  if(abs(x)>abs(g%c(m)))then
    f = f*((x-g%c(m))**4+1)+e*4.0_dp*(x-g%c(m))**3
    e = e*((x-g%c(m))**4+1)
  endif

  la%o%epot = la%o%epot+e 
  la%o%force = la%o%force-f 
  g%epot = g%epot + e

enddo

end subroutine
     
subroutine voter2d_set(g,d1,d2,d3,d4)
type(voter2d)  :: g
real(dp)       :: d1,d2,d3,d4
g%d1=d1
g%d2=d2
g%d3=d3
g%d4=d4
end subroutine
      
subroutine voter2d_interact(g)
!"A method for accelerating the molecular dynamics simulation of infrequent
!events". Arthur F. Voter. J. Chem- Phys 106 (11) 1997
class(voter2d),intent(inout) :: g 
real(dp)                    :: x,y,e
type(atom_dclist),pointer   :: la
integer                     :: i
                
g%epot=0._dp

la => g%alist      ! sobre los atomos
do i = 1,g%nat
  la=>la%next
    

  x=la%o%pos(1)
  y=la%o%pos(2)
  e = cos(twopi*x)*(1+g%d1*y)*ev_ui+2.0_dp*g%d2*(pi*y)**2+g%d3*cos(twopi*x/g%d4)
  e = e +1.203*ev_ui !corrigo para que sea negativo (al menos en la parte donde hay boost)
  la%o%force(1) = twopi*(1+g%d1*y)*sin(twopi*x)*ev_ui+g%d3*twopi/g%d4*sin(twopi/g%d4*x)
  la%o%force(2) =-g%d1*cos(twopi*x)*ev_ui -4*pi**2*g%d2*y
  !la%o%epot     = 0.5_dp*la%o%epot     
  !la%o%force(1) = 0.5_dp*la%o%force(1) 
  !la%o%force(2) = 0.5_dp*la%o%force(2) 

  la%o%epot = la%o%epot+e
  g%epot = g%epot + e
  
enddo
end subroutine
        
subroutine leiva1d_set(g,a,b,alpha,beta)
type(leiva1d)  :: g
real(dp)       :: a,b,alpha,beta
g%alpha=alpha
g%cb=b
g%beta=beta
g%ca=a
end subroutine
      
subroutine leiva1d_interact(g)
class(leiva1d),intent(inout) :: g 
integer                     :: i
real(dp)                    :: x,e
type(atom_dclist),pointer   :: la
               
g%epot=0._dp

la => g%alist      ! sobre los atomos
do i = 1,g%nat
  la=>la%next
    
  x=la%o%pos(1)
  if (x>0) then
    e = g%cb*x**2*exp(-g%beta*x)
    la%o%force(1) = -g%cb*x*exp(-g%beta*x)*(2-x*g%beta)
  else
    e = g%ca*x**2*exp(g%alpha*x)
    la%o%force(1) = -g%ca*x*exp(g%alpha*x)*(2+x*g%alpha)
  endif
  la%o%epot =  la%o%epot + e
  g%epot = g%epot + e
enddo
end subroutine
            
subroutine oscar2d_interact(g)
! Sacada del examen del curso de metodos computacionales 2011
class(oscar2d),intent(inout)   :: g 
integer                       :: i
real(dp)                      :: c,f,    &
                                 h,h1,h2,&
                                 v,v1,v2,&
                                 x,y,b1,b2
real(dp)                      :: e
type(atom_dclist),pointer     :: la
               
g%epot=0._dp

la => g%alist      ! sobre los atomos
do i = 1,g%nat
  la=>la%next
    
  x=la%o%pos(1)
  y=la%o%pos(2)
  c = -3.0_dp*exp(-(x**2))
  h1 = exp(-((y-(5.0_dp/3.0_dp))**2))
  h2 = -exp(-((y-(1.0_dp/3.0_dp))**2))
  f = -5.0_dp*exp(-(y**2))
  v1 = exp(-((x-1.33_dp)**2))
  v2 = exp(-((x+1.0_dp)**2))
  b1 = cos(2.0_dp*pi*x)*(1.0_dp+0.2_dp*y)
  b2 = (1.0e-3_dp*pi*y**2)


  h=h1+h2 
  v=v1+v2
  e = c*h + f*v + b1 + b2
  la%o%force(1) = -2.0_dp*x*c*h-f*2.0_dp*((x-1.33_dp)*v1-(x+1.0_dp)*v2) &
                  -2.0_dp*pi*sin(2.0_dp*pi*x)*(1.0_dp+0.2_dp*y)
  la%o%force(2) = -c*2.0_dp*((y-(5.0_dp/3.0_dp))*h1+(y-(1.0_dp/3.0_dp)*h2))-2.0_dp*y*f*v &
                  + 0.2_dp*cos(2.0_dp*pi*x) - 2.0_dp*pi*3.0_dp*y

  la%o%epot = la%o%epot + e
  g%epot = g%epot + e
  
enddo 

end subroutine
      
subroutine sho2d_set(g,kx,ky)
real(dp),intent(in) :: kx,ky
type(sho2d)         :: g
g%kx=kx
g%ky=ky
end subroutine
                      
subroutine sho2d_interact(g)
class(sho2d),intent(inout) :: g 
integer                   :: i
real(dp)                  :: x,y,e
type(atom_dclist),pointer :: la
               
g%epot=0._dp

la => g%alist      ! sobre los atomos
do i = 1,g%nat
  la=>la%next
    
  x=la%o%pos(1) 
  y=la%o%pos(2)
  la%o%force(1) = -g%kx*x
  la%o%force(2) = -g%ky*y
  e = (g%kx*x**2+g%ky*y**2)*0.5_dp
  !e = e -2*ev_ui !corrigo para que sea negativo
  la%o%epot = la%o%epot + e
  g%epot = g%epot + e

enddo

end subroutine

subroutine sho1d_set(g,p1)
real(dp),intent(in)   :: p1
type(sho1d)           :: g
g%k=p1
end subroutine

subroutine sho1d_interact(g)
class(sho1d),intent(inout)   :: g 
integer                     :: i
real(dp)                    :: x,e
type(atom_dclist),pointer   :: la
                 
g%epot=0._dp

la => g%alist      ! sobre los atomos
do i = 1,g%nat
  la=>la%next
    
  x=la%o%pos(1) 
  la%o%force(1) = -g%k*x

  e= g%k*x**2*0.5_dp

  la%o%epot = la%o%epot + e
  g%epot = g%epot + e

enddo
end subroutine sho1d_interact
   


subroutine sholine_set(g,k,p,v)
real(dp),intent(in)  :: k
type(sholine)        :: g
real(dp),intent(in)  :: p(dm)
real(dp),intent(in)  :: v(dm)

! parametros
g%k=k
g%p(:)=p(:)
g%v(:)=v(:)
g%v(:)=g%v(:)/sqrt(dot_product(g%v(:),g%v(:)))

end subroutine sholine_set

subroutine sholine_interact(g)
class(sholine),intent(inout)   :: g 
type(atom_dclist),pointer         :: la
real(dp)                          :: x(dm),e
integer                           :: i
                 
g%epot=0._dp

la => g%alist
do i = 1,g%nat
  la=>la%next
    
  ! Calcular el vector desplazamiento al eje
  x(:)=la%o%pos(1:dm)-g%p(1:dm)
  x(:)=x(:)-dot_product(x,g%v)*g%v(:)

  la%o%force(:) = la%o%force(:)-g%k*x(:)
  e=g%k*dot_product(x,x)*0.5_dp

  la%o%epot = la%o%epot + e
  g%epot = g%epot + e

enddo

end subroutine sholine_interact


! subroutine shofix_set(g,k)
! real(dp),intent(in)       :: k
! integer                   :: i
! type(atom_dclist),pointer :: la
! type(shofix)              :: g
!    
! g%k=k
!
! allocate(g%r(g%nat,3))
!
! la => g%alist      ! sobre los atomos
! do i = 1,g%nat
!   la=>la%next
!   g%r(i,:)=la%o%pos(:)
! enddo     
!              
! end subroutine shofix_set
!   
! subroutine shofix_interact(g)
! use gems_program_types, only: atom_distancetopoint
! class(shofix),intent(inout)   :: g
! integer                       :: i
! real(dp)                      :: vd(dm),dr,factor2(dm)
! real(dp)                      :: p,f
! type(atom_dclist),pointer :: la
!  
! g%epot = 0._dp
!
! la => g%ref%alist      ! sobre los atomos
! do i = 1,g%nat
!   la=>la%next
!
!   vd = atom_distancetopoint(la%o,g%r(i,:))  ! vector a-->p
!   dr = dot_product(vd,vd) 
!
!   p=-0.5_dp*g%k*dr*dr
!      
!   la%o%epot = la%o%epot + p*0.5_dp
!   g%epot = g%epot + p
!
!   f=-g%k*dr
!
!   factor2(1:dm)=f*vd(1:dm)   ! el dividido dr esta en f
!
!   la%o%force(1:dm) = la%o%force(1:dm) - factor2(1:dm) 
!   
! enddo     
!
! end subroutine  
              
                

subroutine hshoplane_set(g,k,r,l)
type(hshoplane)      :: g
real(dp),intent(in)  :: k, r
integer,intent(in)   :: l
g%k=k
g%r=r
g%l=l
end subroutine hshoplane_set
                     
subroutine hshoplane_interact(g)
class(hshoplane),intent(inout)   :: g 
type(atom_dclist),pointer       :: la
real(dp)                        :: x,e
integer                         :: i,j
                 
j=abs(g%l)
      
g%epot=0._dp

la => g%alist
do i = 1,g%nat
  la=>la%next
    
  ! Calcular el vector desplazamiento al eje
  x=la%o%pos(j)-g%r

  ! Choose half of the harmonic potential
  if (x*g%l<0) cycle

  sumf=sumf-g%k*x
  numf=numf+1

  la%o%force(j) = la%o%force(j)-g%k*x
  e=g%k*x**2*0.5_dp

  la%o%epot = la%o%epot + e
  g%epot = g%epot + e

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


subroutine shoplane_set(g,k,r,j)
type(shoplane)      :: g
real(dp),intent(in) :: k,r
integer,intent(in)  :: j
g%k=k
g%r=r
g%j=j
end subroutine

subroutine shoplane_interact(g)
class(shoplane),intent(inout)   :: g 
type(atom_dclist),pointer         :: la
real(dp)                          :: x,e
integer                           :: i,j
                 
j=g%j
      
g%epot=0._dp

la => g%alist
do i = 1,g%nat
  la=>la%next
    
  ! Calcular el vector desplazamiento al eje
  x=la%o%pos(j)-g%r

  la%o%force(j) = la%o%force(j)-g%k*x
  e=g%k*x**2*0.5_dp

  la%o%epot = la%o%epot + e
  g%epot = g%epot + e

enddo

end subroutine
           
 
subroutine lucas1d_set(g,a1,a2,a3,b1,b2,b3,c1,c2,c3)
real(dp),intent(in)  :: a1,a2,a3,b1,b2,b3,c1,c2,c3
type(lucas1d)        :: g
g%prof1=a1
g%prof2=a2
g%prof3=a3
g%anch1=b1
g%anch2=b2
g%anch3=b3
g%cent1=c1
g%cent2=c2
g%cent3=c3
end subroutine
                  
subroutine lucas1d_interact(g)
class(lucas1d),intent(inout) :: g 
real(dp)                    :: x,e
real(dp)                    :: a1,a2,a3,b1,b2,b3,c1,c2,c3
type(atom_dclist),pointer   :: la
integer                     :: i
               
g%epot=0._dp

! parameters     
a1=g%prof1 
a2=g%prof2 
a3=g%prof3 
                 
b1=g%anch1 
b2=g%anch2 
b3=g%anch3 
                 
c1=g%cent1 
c2=g%cent2 
c3=g%cent3 

la => g%alist      ! sobre los atomos
do i = 1,g%nat
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
  g%epot = g%epot + e

enddo

end subroutine


end module gems_fields
