! Copyright (c) 2020  Sergio Alexis Paz
!
!  This file is part of GEMS. GEMS is an Extensible Molecular Simulator.
!   .
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

module gems_interaction
use gems_output
use gems_neighbour

implicit none

! Point to the ofile element that is use to output this module
type(outfile),pointer         :: epotparts => null()

! Total energy
real(dp), public               :: tepot

contains

subroutine interact_new(ig)
use gems_constants,     only: linewidth
use gems_input_parsing, only: reada, readi
use gems_bias,          only: bias_new
use gems_fields,        only: fields_new
use gems_pairs,     only: pair_new
use gems_tb,            only: smatb_new
use gems_tersoff,       only: tsf_new
class(intergroup),pointer  :: ig
type(group)                :: g1,g2
logical                    :: under, feels
integer                    :: i,i1,i2
character(len=linewidth)   :: w1
   
! Read group for g1
call readi(i1)
call g1%init()
call g2%init()

! Read under/with/feels/set switch  
call reada(w1)
under=.false.
feels=.false. 
selectcase(w1)
case('@','under')
  under=.true.
  call g1%add(gr(i1))
case('<>','with')
  call readi(i2)
  if (i2==i1) then
    under=.true.
    call g1%add(gr(i1))
  else
    call g1%add(gr(i1))
    call g2%add(gr(i2))
  endif
case('<-','feels')
  feels = .true.
  call readi(i2)
  call g1%add(gr(i1))
  call g2%add(gr(i2))
case('->')
  feels = .true.
  call readi(i2)
  call g1%add(gr(i2))
  call g2%add(gr(i1)) 
case default
  call werr('Bad keyword. Use `<>`, `->`, `<-` or `@`')
end select

call werr('No se agrego ningun atomo',g1%nat==0)

call readl(w1)
selectcase(w1)
case('bias')
  call werr("Only self interaction, use `under` keyword.",.not.under) 
  call reada(w1)
  call bias_new(ig,g1,w1)
! case('ff')
!
!   call reada(w1)
!   call read_psf(w1)
!   call reada(w1)
!   call read_prm(w1)

case('wca', 'sm1','sho','shocm','shofix','slj','lj','cuw','rwca','plj') 
  if(under) then 
    call pair_new(ig,g1,w1)
  else
    call pair_new(ig,g1,w1,g2)
  endif 
case('halfsho_plane','sho2d','oscar2d','leiva1d','pozoa1d','voter2d','sho_plane','sho_line','lucas1d') 
  call werr("Only self interaction, use `under` keyword.",.not.under) 
  call fields_new(ig,g1,w1)
case('tb')
  if(under) then 
    call smatb_new(ig,g1)
  else
    call smatb_new(ig,g1,g2)
  endif
case('tersoff')
  call werr("Only self interaction, use `under` keyword.",.not.under)
  call tsf_new(ig,g1)
case default
  call wwan('I do not understand the last command')
endselect    

call g1%dest()
call g2%dest()
 
! Further setups  
if(feels) then
  igr_vop%o(igr_vop%size)%o%half=.false.
  igr_vop%o(igr_vop%size)%o%newton=.false.
endif
        
end subroutine interact_new

subroutine interact(b_out,pot,dontclean)
! calcula la fuerza y energia potencial de los atomos del systema
use gems_tb
use gems_forcefield
use gems_output, only: ptime,pnframe
use gems_cvs, only: cvs, ncvs
use gems_bias, only: bias, cepot

logical,optional,intent(in)    :: dontclean
logical,intent(in)             :: b_out
real(dp),optional,intent(out)  :: pot
integer                        :: i,j
real(dp)                       :: aux
type(boundgr_l),pointer        :: ln
type (atom_dclist), pointer    :: la
class(intergroup), pointer     :: ig
class(atom),pointer            :: ap
  
 ! hard constrain
 !do i = 1,natoms
 !  if (.not.a(i)%o%bconst) cycle
 !  rd=dot(a(i)%o%pos-a(i)%o%pconst,a(i)%o%vconst)
 !  if(rd<1.d-12) cycle
 !  if (a(i)%o%lconst) then
 !    a(i)%o%pos=a(i)%o%pconst+rd
 !  else
 !    a(i)%o%pos=a(i)%o%pos-rd
 !  endif 
 !enddo 
if(nomoreigs) then
  call test_update()
else
  ! First time here.
  ! No more new interaction groups are allowed
  nomoreigs=.true.
  ! Runing set up of ghost and interacting groups
  if(ghost) call pbcghost

  call update()
endif
    
        
bias = 0._dp
tepot = 0._dp
cepot = 0._dp
if(.not.present(dontclean)) then
  virial(:,:) = 0._dp
  do i = 1,natoms
    ap=>a(i)%o
    ap%force = 0._dp
    ap%epot = 0._dp
  enddo
  tepot=0._dp
endif

do i = 1, igr_vop%size
  ig => igr_vop%o(i)%o

  ! Check to skip
  if (ig%disable) cycle

  call ig%interact()

  tepot=tepot+ig%epot
  cepot=tepot
enddo 

! Force fields
ln => boundgrs
do while( associated(ln%next) )
  ln => ln%next
  call ln%o%interact()
  tepot=tepot+ln%o%epot
  cepot=tepot
enddo 

if (b_out) call write_out(2,dm_steps)

if(present(pot)) pot=tepot

! Soft constrain
do i = 1,natoms
  ap=>a(i)%o
  if (.not.ap%bconst) cycle
  if (ap%lconst) then
    aux=dot(ap%force,ap%vconst)
    ap%force = ap%vconst*aux
    aux=dot(ap%acel,ap%vconst)
    ap%acel  = ap%vconst*aux
    aux=dot(ap%vel,ap%vconst)
    ap%vel   = ap%vconst*aux
  else
    aux=dot(ap%force,ap%vconst)
    ap%force = ap%force-aux
    aux=dot(ap%acel,ap%vconst)
    ap%acel  = ap%acel-aux
    aux=dot(ap%vel,ap%vconst)
    ap%vel   = ap%vel-aux
  endif
enddo

! Agroup the forces in the head
do i=1,ncvs
  call cvs(i)%calc()
enddo
      
end subroutine interact
     
subroutine write_eparts(op)
class(outpropa)                :: op
integer                        :: i

call werr('Error. Use epart after interaction delcarations',size(op%f)/=igr_vop%size)

do i=1,igr_vop%size
  op%f(i) = op%f(i) + igr_vop%o(i)%o%epot*ui_ev
enddo  
   
end subroutine

function polvar_interact(var) result(g)
use gems_variables, only: polvar, polvar_find
use gems_errors, only: werr
character(*),intent(in)    :: var
type(polvar),pointer       :: pv
class(intergroup),pointer  :: g


call werr('Labels should start with colon `:` symbol',var(1:1)/=':')
pv=>polvar_find(var)

g=>null()
if(.not.associated(pv)) return

call werr('Variable not linked',pv%hard)
 
! Print
select type(v=>pv%val)
class is (intergroup)
  g=>v
class default
  call werr('I dont know how to return that')
end select

end function polvar_interact
    
end module gems_interaction
 
