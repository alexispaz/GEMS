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
! 
! This file incorporates work derived from the following code:
!
! GRAFPACK <https://people.sc.fsu.edu/~jburkardt/f_src/grafpack/grafpack.html>
! by John Burkard, covered by the following copyright and permission notice:  
!
! Copyright (c) 1999 John Burkard
!
!  The computer code and data files described and made available on this web 
!  page are distributed under the GNU LGPL license.
 
module gems_graphs 
use gems_groups, only:atom, atom_dclist
use gems_neighbor, only:ngroup
use gems_constants, only:isp, dm, dp

implicit none
private

! Many subroutines in this module are taken and adapted from:
!
! https://people.sc.fsu.edu/~jburkardt/f_src/grafpack/grafpack.html
!
! The computer code and data files described and made available on this web page are distributed under [the GNU LGPL
! license](https://www.gnu.org/licenses/lgpl-3.0.en.html).
    
type, extends(ngroup) :: graph_adj
  integer(isp),allocatable       :: adj(:,:)
  integer(isp)                   :: nblocks
  contains
  procedure   :: interact => graph_block
  procedure,nopass :: cli => graph_cli
end type
                                      
type, extends(ngroup) :: graph
  integer,allocatable       :: order(:)
  integer,allocatable       :: stack(:)
  integer,allocatable       :: label(:)
  integer                   :: ngraphs
  contains
  procedure   :: interact => graph_subgraphs
  procedure,nopass :: cli => graph_cli
end type
                                      
public  :: graph_new, write_graph

class(ngroup),pointer   :: gr=>null()

contains

subroutine graph_cli(g)
use gems_input_parsing, only:readf,reada
use gems_neighbor, only:ngroup
use gems_constants, only:linewidth
use gems_errors, only:werr
class(ngroup),intent(inout)  :: g
real(dp)                     :: rc
character(len=linewidth)     :: w1
integer                      :: i

i=g%amax
select type(g)
type is(graph_adj)  
  allocate(g%adj(i,i))
type is(graph) 
  allocate(g%order(i))
  allocate(g%stack(i))
  allocate(g%label(i))
class default
  call werr('Interaction type mismatch. Expected graph type.',.true.)
end select  


call readf(rc)
call g%setrc(rc)

end subroutine graph_cli
 
subroutine graph_new(g,w)
use gems_neighbor, only:ngroup
use gems_errors, only:wwan,werr
use gems_groups, only:group
class(ngroup),pointer :: g
character(*),intent(in)  :: w
                                            
! TODO: Necesito hacerlo para más de un
! graph type pero para eso tengo que acomodar
! las subrrutinas de escritura
call werr('Not possible yet...',associated(gr))
      
select case(w)
case('subgraphs') ; allocate(graph::g)
case('blocks')    ; allocate(graph_adj::g)
  call wwan("Output routine under construction.") 
  ! Pensar si se puede hacer una subrrutina dedicada de escritura por type
case default
  call werr('Graph function not found',.true.)
endselect
  
! FIXME Internal 
gr=>g
    
end subroutine graph_new
     
subroutine graph_block(g)
! Search neighbors and update adj matrix accordingly
! Then call to graph_adj_block to compute the biconected components (blocks).
use gems_program_types, only: mic
use gems_inq_properties, only: vdistance
class(graph_adj),intent(inout)  :: g
integer                         :: i,ii,j,l
real(dp)                        :: vd(dm),dr
type(atom),pointer              :: o1,o2
type(atom_dclist),pointer       :: la

! Set zeros
g%adj(:,:)=0
 
! Search edges and fill adj matrix
la => g%ref%alist
do ii = 1,g%ref%nat
  la => la%next
  o1 => la%o
  i = o1%gid(g)

  ! g%adj(i,i)=1
  do l = 1, g%nn(i)  ! sobre los vecinos

    j  = g%list(i,l)
    o2 =>g%a(j)%o

    ! Skip atoms in limbo
    if(g%b_limbo) then
      if(associated(o2,target=g%limbo)) cycle
    endif

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>g%rcut2) cycle

    g%adj(i,j)=5
 
  enddo     

enddo 

call graph_adj_block(g%adj, g%nblocks)
                            
end subroutine graph_block

subroutine graph_subgraphs ( g ) 
! Deep-first search algoritm to find subgraphs  
use gems_program_types, only: mic
use gems_inq_properties, only: vdistance
class(graph),intent(inout)    :: g
integer                       :: i,j,l,k,lstack
real(dp)                      :: vd(dm),dr
type(atom),pointer            :: o1,o2

g%order(:) = 0
g%stack(:) = 0
g%label(:) = 0

k = 0         ! Visiting order counter
i = 1         ! ID of visited atom
lstack = 0    ! Stack of visited atoms
g%ngraphs=1  ! Current subgraph number

! Update here since clean index is asumed below.
call g%update_index()

main: do
  k = k + 1
  o1=>g%a(i)%o

  ! If never visited, save order and increase stack
  if ( g%order(i) == 0) then
    ! The subgraph number
    g%label(i)=g%ngraphs
    
    ! Save order
    g%order(i) = k

    ! Increase stack and save the node ID
    lstack = lstack + 1
    g%stack(lstack) = i
  endif

  !  Check the next neighbor.
  do l = 1, g%nn(i)  ! sobre los vecinos

    j  = g%list(i,l)

    ! If already visited, skip
    if ( g%order(j) /= 0) cycle

    o2 =>g%a(j)%o

    ! Skip atoms in limbo
    if(g%b_limbo) then
      if(associated(o2,target=g%limbo)) cycle
    endif

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>g%rcut2) cycle
 
    ! Depth-first search. Follow the first (non visited) nieghbor found.
    ! Node labels will save the order of the root node. 
    i = j
    cycle main

  end do
  
  ! Check if there is nodes in the stack
  lstack = lstack - 1
  if ( lstack > 0 ) then

    ! Go back in the stack to the previous node and search neighboors from it (i.e. depth search from second neighbor).
    j = i
    i = g%stack(lstack)
    cycle main

  else

    ! Stack is exhausted, look for a node we haven't visited yet.
    do l = 1,g%nat

      if ( g%order(l) /= 0 ) cycle
      g%ngraphs = g%ngraphs+1
      i = l
      cycle main

    enddo 
 

    exit main
    
  end if

end do main

end subroutine graph_subgraphs
     
subroutine graph_adj_block ( adj, nblock ) 
!Blocks of an undirected graph from its adjacency list.
!
!Definition:
!        ngindex%o(i1)%o
!  A component of a graph is a connected subset of the graph.  If a node
!  is in the component, then all nodes to which it is connected are also
!  in the component.
!
!  An articulation point of a component of a graph is a node whose
!  removal causes the component to no longer be connected.
!
!  A component with no articulation points is called a block.
!
!Licensing:
!
!  This code is distributed under the GNU LGPL license. 
!
!Modified:
!
!  15 April 1999
!
!Reference:
!
!  Alan Gibbons,
!  Algorithmic Graph Theory,
!  Cambridge University Press, 1985,
!  ISBN 0-521-28881-9.
!
!Parameters:
!
!  Input/output, integer(isp) ADJ(LDA,NNODE).
!  On input, ADJ is the adjacency matrix.  ADJ(I,J) is
!  positive if there is an edge from node I to node J, and 0 otherwise.
!  On output, each positive entry of ADJ has been replaced
!  by the number of the block that the corresponding edge belongs to.
!
!  Input, integer(isp) LDA, the leading dimension of ADJ, which must
!  be at least NNODE.
!
!  Input, integer(isp) NNODE, the number of nodes.
!
!  Output, integer(isp) DAD(NNODE), DAD(I) is the node from which
!  node I is visited.  Node 1 is the first node in the search,
!  and has no predecessor, so DAD(1) is zero.  If there is
!  more than one connected component in the graph, then there
!  will be other nodes with DAD equal to zero.
!
!  Output, integer(isp) ORDER(NNODE).  ORDER(I) records the order
!  in which the node was visited during the depth-first search.
!  The first node, node 1, has ORDER(1) = 1.
!  Note, however, that any node which is an articulation point
!  will have the value of ORDER negated.
!
!  Workspace, integer STACK(NNODE).
!
!  Output, integer(isp) NBLOCK, the number of blocks.

integer(isp)   :: adj(:,:)
integer(isp)   :: lda
integer(isp)   :: nnode

integer(isp) dad(size(adj,2))
integer(isp) i
integer(isp) idir
integer(isp) ii
integer(isp) inode(size(adj,2)*(size(adj,2)-1)/2)
integer(isp) order(size(adj,2))
integer(isp) iroot
integer(isp) j
integer(isp) jedge
integer(isp) jj
integer(isp) jnode(size(inode))
integer(isp) k
integer(isp) l
integer(isp) label(size(adj,2))
integer(isp) lstack
integer(isp) nblock
integer(isp) stack(size(adj,2))


lda=size(adj,1)
nnode=size(adj,2)

dad(1:nnode) = 0
inode(1:nnode) = 0
order(1:nnode) = 0
stack(1:nnode) = 0
jnode(1:nnode) = 0
label(1:nnode) = 0

nblock = 0
k = 0
i = 1
lstack = 0
jedge = 0

!Find all descendants of the parent node in this connected component
!of the graph.

iroot = 1

main: do
  k = k + 1

  order(i) = k
  label(i) = k

  ! Increase stack and save the node ID
  lstack = lstack + 1
  stack(lstack) = i

  idir = + 1

   
  j = 0

  !  Check the next neighbor.
  neigh: do

    j = j + 1

    ! Depth-first search. Follow the first (non visited) nieghbor found.
    ! Node labels will save the order of the root node. 
    if ( nnode >= j ) then

      ! Cero matrix element. Not a neighboor, skip
      if ( adj(i,j) == 0 .or. adj(j,i) == 0 ) cycle

      ! Positive matrix element. Save the edge.
      if ( 0 < adj(i,j) .or. 0 < adj(j,i) ) then
        jedge = jedge + 1
        inode(jedge) = i
        jnode(jedge) = j
      end if

      if ( order(j) == 0 ) then

        ! If never visited before, continue neighbor search from it
        dad(j) = i
        i = j
        cycle main

      else

        if ( idir == +1 ) then
          label(i) = min ( label(i), abs ( order(j) ) )
        else
          label(i) = min ( label(i), label(j) )
        end if

        cycle

      end if

    endif

    ! Go back in the stack to the previous node. 
		! If its not the root, search neighboors from it (i.e. depth search from
    ! second neighbor).
		! If its the root...???
    ! If stack is exhausted, look for a node we haven't visited,
    ! which therefore belongs to a new connected component. 

    !Searched all directions from current node.  Back up one node,
    !or, 
     
    lstack = lstack - 1
    idir = -1

    if ( 0 < lstack ) then

      j = i
      i = stack(lstack)

      if ( abs ( order(i) ) <= label(j) ) then

        if ( 0 < order(i) ) then

          if ( i /= iroot ) then
            order(i) = - order(i)
          else
            iroot = 0
          end if

        end if

        nblock = nblock + 1

        do

          ii = inode(jedge)
          jj = jnode(jedge)
          jedge = jedge - 1
          adj(ii,jj) = - nblock
          adj(jj,ii) = - nblock

          if ( ii == i .and. jj == j )  exit

        end do

      end if

    else

      ! Search for new not visited node
      lstack = 0
      do l = 1, nnode
        if ( order(l) == 0 ) then
          i = l
          iroot = i
          cycle main
        end if
      end do

      exit main
      
    end if

  end do neigh
end do main



!Restore the positive sign of the adjacency matrix.

adj(1:nnode,1:nnode) = abs(adj(1:nnode,1:nnode))

end subroutine graph_adj_block

subroutine write_graph(of)
use gems_output, only: outfile
use gems_groups, only: atom_dclist, atom
use gems_elements, only: csym
class(outfile)       :: of
type(atom),pointer   :: o
type(atom_dclist),pointer :: la
integer              :: i,j

select type(gr)
type is(graph)
write(of%un,*) gr%nat
write(of%un,*) gr%ngraphs

la => gr%alist
do i = 1,gr%nat
  la => la%next
  o => la%o
  write(of%un,'(a'//csym//',3(2x,e25.12),x,i0)') o%sym,(o%pos(j),j=1,dm),gr%label(i)
enddo

if(of%flush) call flush(of%un)
end select

end subroutine

subroutine select_connected(g,gini,gout)
! Add to gout all atoms in graph conected to gini
use gems_groups, only: group
class(graph),intent(inout) :: g
type(group),intent(in)     :: gini
type(group),intent(inout)  :: gout
real(dp)                   :: rd,rad,ctr(dm)
type (atom),pointer        :: at
integer                    :: i
logical, allocatable       :: lmask(:)

! Update index... just in case (TODO: check)
call g%update_index()

! Compute subraphs  
call g%interact()
allocate(lmask(g%ngraphs))

! Get labels from atoms in gini and save them to lmask
lmask(:)=.false.
do i = 1,g%nat
  at=>g%a(i)%o
  if(at%gri(gini)/=-1) lmask(g%label(i))=.true.
enddo
 
! Search the igaph again and peek atoms by its labels.
do i = 1,g%nat
  at=>g%a(i)%o
  if(lmask(i)) call gout%attach(at)
enddo  

deallocate(lmask)

end subroutine select_connected

end module gems_graphs
