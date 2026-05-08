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
public    :: calc_cli, write_calc
       
type, extends(cgroup) :: widom

  real(dp)            :: z1,z2
  real(dp)            :: mu
  integer             :: nadj

  contains

  procedure   :: calc => widom_calc
  procedure   :: write => widom_write

end type
                                      
type(group_vop),target,public   :: cindex

contains

subroutine calc_cli(gsel)
use gems_input_parsing, only:readl,getf,geti
use gems_neighbor, only:ngroup
use gems_constants, only:linewidth
use gems_errors, only:werr, wstd
use gems_strings, only: str
use gems_variables, only: polvars
class(group),intent(in)           :: gsel
class(group),pointer              :: g=>null()
character(:),allocatable          :: cid,key,label
integer                           :: id
        
! Read object label if given
call readl(label)

! Read object key and create label if was not given
if(label(1:1)==':') then
  call readl(key)
else  
  key=label
  label=':c'//str((cindex%size+1))
  call wstd('New calc label is '//label)
endif

! Find object id if exist
cid=polvars%expand(label)

! Create object if does not exist
if(cid=='NODEF') then

  ! Allocate object
  select case(key)
  case('widom')
    allocate(widom::g)
  case default
    call werr('Unkown calculation request.',.true.)
  end select  

  ! Initialize object
  call g%init()

  ! Index object
  call cindex%append()
  id=cindex%size
  cindex%o(id)%o=>g

  ! Save the new label
  call polvars%hard(label,id)

else

  read(cid,*) id
  g=>cindex%o(id)%o

endif

! Run CLI
select type(g)
type is(widom) 
  g%z1  =getf()
  g%z2  =getf()
  g%rcut=getf()
  g%nadj=geti()
  call g%attach(gsel)
class default
  call werr('Unkown requested calculation.',.true.)
endselect

end subroutine calc_cli
 
subroutine widom_calc(g)
! For rigid spheres and considering a particular area between two xy planes.
! TODO: Generalize with a boltzman energy
! TODO: Check overlap with a reference group (for mixtures)
use gems_program_types, only: box, distance
use gems_random, only: ranu
use gems_constants, only: dm, kB_ui, ui_ev
class(widom)               :: g
real(dp)                   :: r(dm), dr
type(atom),pointer         :: o
integer                    :: n,m
integer                    :: i,j
real(dp)                   :: vd(dm),v
type(atom_dclist),pointer  :: la
  
g%mu=3.14_dp
 
! Compute volume
v=box(1)*box(2)*(g%z2-g%z1)
 
! Count particles in the control volume
n=0
la=>g%alist
do j=1,g%nat
  la=>la%next
  if(la%o%pos(3)<g%z1.or.la%o%pos(3)>g%z2) cycle
  n=n+1
enddo
    
! Attempted insertions
m=0
adj: do i=1,g%nadj

  ! Random coordinates
  r(1)=ranu()*box(1)
  r(2)=ranu()*box(2)
  r(3)=ranu()*(g%z2-g%z1)+g%z1
  
  ! ! Sort particles
  ! call g%tessellate()
  ! call g%sort()
  
  ! Check overlap
  la=>g%alist
  do j=1,g%nat
    la=>la%next
    o => la%o

    ! ! Skip particles outside the control volume.
    ! FIXME: consider PBC
    ! if(o%pos(3)<z1-rc) cycle
    ! if(o%pos(3)>z2+rc) cycle

    ! Skip if overlapping
    vd(:) = distance(o%pos,r,o%pbc)
    dr = dot_product(vd,vd)
    if(dr<g%rcut*g%rcut) cycle adj

  enddo

  m=m+1     
enddo adj
  
! g%mu=log(n*g%rcut**3/v)-log(real(m,dp)/g%nadj)
! g%mu=g%mu/kB_ui*ui_ev

! Activdada
g%mu=n/v*g%nadj/m
    
end subroutine widom_calc
  
subroutine widom_write(g, unit, iotype, v_list, iostat, iomsg)
use gems_strings, only: str
class(widom),intent(in)    :: g
integer,intent(in)         :: unit,v_list(:)
integer,intent(out)        :: iostat
character(*),intent(in)    :: iotype
character(*),intent(inout) :: iomsg  

character(:),allocatable   :: wfmt

! Format descriptor
select case(size(v_list))
case(0) ! default
  wfmt = '(e25.12)'
case(2)
  wfmt = '(e'//str(v_list(1))//'.'//str(v_list(2))//')'
case default
  iostat = 1
  iomsg = 'wrong number of format descriptors'
  return
end select

! Output type 
select case(iotype)
case('DT','DTmu')
  call g%calc() ! FIXME
  write(unit,wfmt,iostat=iostat,iomsg=iomsg)  g%mu
case default
  iostat = 1
  iomsg = 'unexpected iotype'
end select

end subroutine
 
subroutine write_calc(of)
use gems_output
class(outfile)  :: of
integer         :: i

! FIXME
do i = 1,cindex%size
  ! write (*,'(dt"mu")') cindex%o(i)%o
  write (of%un,'(dt)') cindex%o(i)%o
enddo
if(of%flush) call flush(of%un)

end subroutine
                 
end module gems_calc
