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

module gems_cvs

! Tengo la diea de definir un objeto CV que tenga informacion de el o los grupos
! asociados con la variable colectiva, los parametros (pueden ser arreglos!!!,
! revisar porque se usa vectores en los igr) el procedimiento que se debe
! efectuar para evaluar la CV y/o su jacobiano y los flags para determinar si es
! necesario actualizar su valor o no. La idea es que este objeto permitiria
! eliminar toda la informacion relativo a propiedades (o CVs) del objeto grupo.
! Sin embargo, el problema de que una CV puede ser simplemente una pequeña
! transformacion de otra CV no se bien como reflejarlo. Es decir, el
! inq_triaxialparam involucra un inq_inercia, y si el primero involucra un objeto
! CV el segundo deberia involucrar otro objeto CV diferente, pero no hay
! informacion de acceso a este otro CV. O al menos hacer 1 solo objeto CV por
! grupo que contenga todas las CVs asociadas con ese grupo? Pero en ese caso,
!       este objeto CV deberia tener un arreglo de punteros a los procedimientos
     
use gems_constants,only:dp,dm
! use gems_algebra,only:real_v,integer_v
use gems_groups, only:igroup,group_ap,atom,atom_dclist
use gems_inq_properties, only:group_inq_cmpos,inq_cm_vel,group_inq_rg
use gems_errors, only:werr

implicit none
private

public cv_eval_cmpos, cv_eval_cm, cv_calc_sho, cv_eval_rg
public cv_jaco_cmpos, cv_jaco_cm
  
type,extends(igroup),public :: cv
 
  ! Array of asociated groups (CV that use many groups)
  ! type(group_ap),allocatable  :: gs(:)
          
  ! Parametros reales de la variable colectiva
  ! type(real_v) :: pr
  real(dp),allocatable :: pr(:)
  
  ! Parametros enteros de la variable colectiva
  ! type(integer_v) :: ir
  integer              :: dm=0
  integer,allocatable  :: ir(:)
             
  ! real(dp),allocatable :: zt(:)     ! The CV
  real(dp),allocatable       :: t(:)     ! The CV
  real(dp),allocatable       :: tf(:)    ! The force in the CV given the restraint

  ! XXX: Should I use an atom instead of a CV here??
  ! after all, a CV is not mor than an atom with an arbitrary dimension
  real(dp),allocatable       :: z(:)     ! The restraint

  ! The jacobian, it has 3 dimensions, because the last 2 are atom number and
  ! coordinate, but note that is just a convention
  real(dp),allocatable       :: j(:,:,:) 
                       
  procedure(cv_eval),pointer :: eval=>null() 
  procedure(cv_eval),pointer :: jaco=>null() 
  procedure(cv_calc),pointer :: calc=>null() 

  ! Flag to avoid duplicated calculations
  logical                    :: b_flag=.false.

  contains

  procedure :: attach_atom => cv_attach
  ! Constructor
  ! procedure :: init => cv_init

  ! La interfaz del procedimiento para set no es conocido a priori ya que no
  ! se sabe el numero de parametros necesario. Por ello, hay que llamarlo con
  ! un call directo como se hace con los potenciales
  ! procedure :: cv_set

end type cv
             
abstract interface
 subroutine cv_eval(c)
   import cv
   class(cv)               :: c
 end subroutine
end interface
     
! #define _NODE cv_v
! #define _CLASS type(cv)
! #include "vector_header.inc"
            
! ! Collective variables that depends on the position
! type(cv_v),public     :: cvp
!         
! ! Collective variables that depends on the velocity
! type(cv_v),public     :: cvv
            
abstract interface
 subroutine cv_calc(c)
   import cv
   class(cv)               :: c
 end subroutine
end interface
                                  
! TODO: generalize the number of CVs
type(cv),public     :: cvs(10)
integer,public     :: ncvs=0

contains
       
! #define _NODE cv_v
! #define _CLASS type(cv)
! #include "vector_body.inc"
  
subroutine cv_attach(g,a,l_)
class(cv),target     :: g
class(atom),target   :: a
integer              :: n, m
integer,intent(out),optional :: l_
integer                      :: l

! Attempt to attach
call g%igroup_attach_atom(a,l)
if(present(l_)) l_=l
if(l==0) return
                                
! Reallocate if needed
if(allocated(g%j)) then
  if(g%nat<size(g%j,2)) return
  deallocate(g%j)
endif

n=g%nat+g%pad
m=g%dm

allocate(g%j(m,n,dm))
allocate(g%t(m))
allocate(g%tf(m))
allocate(g%z(m)) 
 
end subroutine cv_attach
            
!   
! subroutine cv_setg(c,ind,g)
! ! Set pointer to group
! type(group),target   :: g  
! class(cv)            :: c
! integer,intent(in)   :: ind
!
! call werr('group index out of cv range',size(c%gs)<ind)
! c%gs(ind)%o=>g
! end subroutine cv_setg

! subroutine cv_eval(c)
! class(cv)                  :: c
!
! if(c%b_flag) return
!
! c%t(:)=c%eval()
!
! c%b_flag=.true.
! end function
!       
! ! TODO: Aca hacer que tome CV y no grupos. Esto implica decir
! ! que nunca le voy a preguntar a un grupo de atomos una propiedad
! ! de forma directa. Que hacer con las dependencias?
! function inq_cmpos(g) result(cmpos)
! class(group)               :: g
! real(dp)                   :: cmpos(dm)
! type(atom_dclist),pointer  :: la
! integer                    :: i
!
! la => g%alist
! do i = 1,g%nat
!   la => la%next
!   cmpos(1:dm) = cmpos(1:dm) + la%o%mass * la%o%pos(1:dm)
! enddo
! cmpos(1:dm) = cmpos(1:dm)/inq_mass(g)
!
! end function
                           
subroutine cv_jaco_cmpos(c)
class(cv)                  :: c
type (atom_dclist),pointer :: la
integer                    :: i,j
logical,save               :: pragmaonce=.false.

! Solo se ejecuta una ves
if(pragmaonce) return

la => c%alist
do i=1,c%nat
  la => la%next 

  c%j(1:dm,i,1:dm) = 0._dp
  do j=1,dm
    c%j(j,i,j) = la%o%mass/c%mass
  enddo
enddo       

pragmaonce=.true.

end subroutine cv_jaco_cmpos
                                         
subroutine cv_jaco_cm(c)
class(cv)                  :: c
type (atom_dclist),pointer :: la
integer                    :: i,j
logical,save               :: pragmaonce=.false.

! Solo se ejecuta una ves
if(pragmaonce) return

j=c%ir(1)  
c%t(1)=c%cm_pos(j)

la => c%alist
do i=1,c%nat
  la => la%next 
  c%j(:,i,:) = 0._dp
  c%j(1,i,j) = la%o%mass/c%mass
enddo       
            

pragmaonce=.true.

end subroutine cv_jaco_cm
            

subroutine cv_eval_cm(c)
class(cv)                  :: c
integer                    :: j

! Add the COM pos/vel to the head
call group_inq_cmpos(c) 

j=c%ir(1)  
c%t(1)=c%cm_pos(j)

end subroutine cv_eval_cm
             
subroutine cv_eval_cmpos(c)
class(cv)                  :: c
type (atom_dclist),pointer :: la
integer                    :: i,j

! Add the COM pos/vel to the head
call group_inq_cmpos(c) 
c%t(1:dm)=c%cm_pos(1:dm)

! TODO: This part is only needed once
la => c%alist
do i=1,c%nat
  la => la%next 

  c%j(1:dm,i,1:dm) = 0._dp
  do j=1,dm
    c%j(j,i,j) = la%o%mass/c%mass
  enddo
enddo       

end subroutine cv_eval_cmpos

subroutine cv_eval_rg(c)
class(cv)                  :: c
type (atom_dclist),pointer :: la
integer                    :: i

! Add the COM pos/vel to the head
call group_inq_rg(c) 
c%t(1)=c%rg_pos

la => c%alist
do i=1,c%nat
  la => la%next 
  c%j(1,i,:)=(la%o%pos(:)-c%cm_pos(:))/(c%rg_pos*c%nat)
enddo       

end subroutine cv_eval_rg
              
subroutine cv_calc_sho(c)
class(cv)         :: c
real(dp)          :: k
type (atom_dclist),pointer :: la
integer           :: i,j

! Parameters
k=c%pr(1)
     
call c%eval()

! Hay ciertas CVs que esto no es necesario  
call c%jaco()

c%tf(:)=-2._dp*k*(c%t(:)-c%z(:))

la => c%alist
do i=1,c%nat
  la => la%next 
  do j=1,dm
    la%o%force(j) = la%o%force(j) + dot_product(c%j(:,i,j),c%tf(:))
  enddo
enddo       
    
end subroutine cv_calc_sho
          
end module gems_cvs
