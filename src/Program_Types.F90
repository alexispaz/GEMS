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


module gems_program_types
use gems_groups
use gems_constants,only:sp,dp,dm,find_io
use gems_elements,only:elements,ncsym,element,inq_z

implicit none
save
private

! The program and subprograms trace this variable to get the term order
! After the term_signal is read, is set to .false. too allow subprograms end
! without quit gems.
logical,public :: term_signal=.false., cycle_signal=.false.


! Program variables
real(dp),target,public   :: dm_steps=0._dp
real(dp),public,target   :: time=0.0_dp

integer,public            :: frame=0 ! Cuando escribe: el frame actual del outunit
integer,public            :: nstep=100 ! number of integration step

real(dp),target,public    :: dt=0.002

public :: box_setvars,box_expand,boxed
real(dp),public,dimension(dm,dm) :: tbox   =0.0_dp
real(dp),public,dimension(dm),target    :: box    =1.0e6_dp  
real(dp),public,dimension(dm)    :: box2   =5.0e5_dp  , &
                                    one_box=1.0e-6_dp , &
                                    one_box2=2.5e-5_dp, &
                                    box_old=0.0_dp      ! For pbcghost. Set this small to force initial ghost inclusion
real(dp),public                  :: box_vol=1.0e18_dp
logical                          :: boxed=.false.

real(dp),public             :: virial(dm,dm)=0._dp
logical,public              :: b_gvirial=.false.

! Objects
! =======

! Atoms
! -----

! Lista blanda con atomos fantasma
logical,public                      :: mic=.true.

public :: atom_distancetopoint
public :: atom_distancetoaxis
public :: vdistance
public :: do_pbc, set_pbc

! Groups
! ------
type(igroup),target,public  :: sys
    
contains
                      
! Box
! ===
           
subroutine box_setvars()
use gems_input_parsing

boxed=.true.

! For non-cubic boxes this should be box(:)=matmul(tbox(i,i),[1,1,1])
box(1)=tbox(1,1)
box(2)=tbox(2,2)
box(3)=tbox(3,3)

one_box(:) = 1.0_dp/box(:)

! For non-cubic boxes this should be the determinant of tbox
box_vol = box(1)*box(2)*box(3)

end subroutine box_setvars

subroutine box_expand(a,b,c)
real(dp)       :: a,b,c

tbox(1,1)=tbox(1,1)*a
tbox(2,2)=tbox(2,2)*b
tbox(3,3)=tbox(3,3)*c

call box_setvars()

end subroutine box_expand

! PBC
! ---

subroutine set_pbc(g,pbc)
! rota un grupo un angulo r sobre el punto v
class(group),intent(inout)  :: g
class(atom_dclist),pointer  :: la
logical,intent(in)         :: pbc(dm)
integer                    :: i

la => g%alist
do i = 1,g%nat
  la => la%next
  la%o%pbc = pbc
enddo

end subroutine set_pbc

subroutine do_pbc(g)
class(group),intent(inout)  :: g
class(atom_dclist),pointer :: la
integer                    :: i,k

la => g%alist
do i = 1,g%nat
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
     
end subroutine do_pbc
   
function vdistance(i,j,mic) result(vd)
!calculates the distance of two atoms with or without minimum image convention
real(dp),dimension(dm)  :: vd
type(atom),intent(in)   :: i,j
logical,intent(in)      :: mic
logical                 :: pbc(dm)
integer                 :: l
  
! Distancia
vd=i%pos-j%pos

! Convencion de imagen minima
if(.not.mic) return

! Mas rapido usar idnint que un if
pbc=i%pbc.or.j%pbc
do l = 1,dm
  if (pbc(l)) vd(l)=vd(l)-box(l)*idnint(vd(l)*one_box(l))
enddo

end function vdistance
           
! Distancias
! ----------

function atom_distancetopoint(a,r) result(vd)
!calculates the distance of two atoms also in pbc case
class(atom),intent(in)   :: a
real(dp),dimension(dm)  :: vd
real(dp),intent(in)     :: r(dm)
integer                 :: l

! Mas rapido usar idnint que un if
vd=r-a%pos
do l = 1,dm
  if (a%pbc(l)) vd(l)=vd(l)-box(l)*idnint(vd(l)*one_box(l))
enddo

end function atom_distancetopoint

function atom_distancetoaxis(a,p,r) result(vd)
! a es el atomo, p y r dan la recta(t)=p+t*r
!calculates the distance of two atoms also in pbc case
class(atom),intent(in)   :: a
real(dp),dimension(dm)  :: vd,aux
real(dp),intent(in)    :: r(dm),p(dm)

 aux=atom_distancetopoint(a,p)
 vd=aux-(dot_product(aux,r))*r

end function atom_distancetoaxis

end module gems_program_types
 
