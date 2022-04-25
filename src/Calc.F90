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
!
! ---
 
module gems_calc
use gems_groups, only: atom, group, group_vop, atom_dclist
use gems_cells, only: cgroup, map
use gems_constants, only:isp, dm, dp

implicit none

private
public    :: calc_cli
       
type, extends(cgroup) :: widom
  real(dp)            :: z1,z2
  real(dp)            :: mu
  integer             :: nadj
  contains
  procedure           :: calc => widom_calc
end type
                                      
type(group_vop),target,public   :: cindex

contains

subroutine calc_cli()
use gems_input_parsing, only:readl,getf,geti
use gems_neighbor, only:ngroup
use gems_constants, only:linewidth
use gems_errors, only:werr, wstd
use gems_strings, only: operator(.ich.)
use gems_variables, only: polvar_expand, polvar_hard
class(group),pointer      :: g
real(dp)                  :: rc
character(:),allocatable  :: cid,key,label
integer                   :: id
integer                   :: inpi(5)
real(dp)                  :: inpf(5)
        
! Read object label if given
call readl(label)

! Read object key and create label if was not given
if(label(1:1)==':') then
  call readl(key)
else  
  key=label
  label=':c'//.ich.(cindex%size+1)
  call wstd('New calc label is '//label)
endif

! Find object id if exist
cid=polvar_expand(label)

! Create object if does not exist
if(cid=='NODEF') then

  ! Allocate object
  select case(key)
  case('widom'); allocate(widom::g)
  case default
    call werr('Unkown calculation request.',.true.)
  end select  

  ! Index object
  call cindex%append()
  id=cindex%size
  cindex%o(id)%o=>g

  ! Save the new label
  call polvar_hard(label,id)

else

  read(cid,*) id
  g=>cindex%o(id)%o

endif

! Run CLI
select type(g)
type is(widom) 
  g%z1  =getf()
  g%z2  =getf()
  g%rcut = getf()
  g%nadj=geti()
class default
  call werr('Unkown requested calculation.',.true.)
endselect

end subroutine calc_cli
 
subroutine widom_calc(g)
! For rigid spheres and considering a particular area between two xy planes.
! TODO: Generalize with a boltzman energy
! TODO: Check overlap with a reference group (for mixtures)
use gems_groups, only: ghost_from_atom, useghost
use gems_neighbor, only: ngindex,ngroup
use gems_random, only: ranu,rang
use gems_constants, only: dm, kB_ui
class(widom)               :: g
real(dp)                   :: r(dm), dr
type(atom),pointer         :: o
integer                    :: n,m
integer                    :: i,ii,j,ic,k
integer                    :: nabor
real(dp)                   :: rd,vd(dm)
integer                    :: rc(dm),nc(dm)
type(atom_dclist),pointer  :: la
  
! ! Count particles in the control volume
! n=0
! la=>g%alist
! do j=1,g%nat
!   la=>la%next
!   if(la%o%pos(3)<g%z1.or.la%o%pos(3)>g%z2) cycle
!   n=n+1
! enddo
!     
! ! Point to an atom that will work as template
! ! in order to add new atoms into groups.
!        
! ! Attempted adjustments 
! adj: do i=1,g%nadj
!
!   ! Random coordinates
!   r(1)=ranu()*box(1)
!   r(2)=ranu()*box(2)
!   r(3)=ranu()*(g%z2-g%z1)+g%z1
!   
! ! Sort particles
! g%rcut=rad2
! call g%tessellate()
! call g%sort()
!   
! enddo adj
!       
end subroutine widom_calc
  
subroutine write_widom(of)
use gems_output, only: outfile
use gems_groups, only: atom_dclist, atom
use gems_elements, only: csym
class(outfile)       :: of
type(atom),pointer   :: o
type(atom_dclist),pointer :: la
integer              :: i,j

! select type(gr)
! type is(widom)
! write(of%un,*) gr%nat
! write(of%un,*) gr%nwidoms
!
! la => gr%alist
! do i = 1,gr%nat
!   la => la%next
!   o => la%o
!   write(of%un,'(a'//csym//',3(2x,e25.12),x,i0)') o%sym,(o%pos(j),j=1,dm),gr%label(i)
! enddo
!
! if(of%flush) call flush(of%un)
! end select
!
end subroutine

end module gems_calc
