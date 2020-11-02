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

 
module gems_bondboost
!from 
!Y. Shibuta and S. Maruyama, "Bond-order potential for transition metal
!carbide cluster for the growth simulation of a single-walled carbon nanotube,"
!Computational Materials Science, vol. 39, 2007, pp. 842-848.

  use gems_atoms, only: atom,atom_dclist
  use gems_groups, only: group
  use gems_constants, only:dp,dm
  use gems_neighbour

  implicit none
  private

 
  ! type :: bondb
  !   real(dp)               :: re=0.d0,              & ! distancia de eqilibrio
  !                             rd,vd(dm),            & ! distancia actual
  !                             d=0.d0,               & ! distorcion abs(rd-re)/re
  !                             t=0.d0,               & ! tiempo caracteristico de cada intento
  !                             w=0.d0,               & ! weight
  !                             et=1.d0                 ! error relativo del tiempo caracteristico
  !   type(atom),pointer     :: i,j                     ! atomos del enlace
  ! end type bondb
  !


!  type(bondb),pointer  :: bb
!  real(dp)            :: sa(msubs,msubs)=0
!  real(dp)            :: boost=0.0_dp
!  public bb,sa,rq2,qlim,alpha,bondb
!
!
!  ! Esto estaba para el simplify bondb boost, pero ya no esta mas
!  real(dp)            :: eq=0.5_dp,eq2=0.025_dp
!
!
!  real(dp)            :: rq2(msubs,msubs)
  

  ! type(bondb),target,allocatable  :: bonds(:)
  ! integer                         :: nbonds

  !public generate_dispboost,longerdisp_lenght,bondb_disp_lenght

!  real(dp)             :: qlim=0
  real(dp)             :: alpha=0

  ! type :: linked_vector
  !   type(linked_vector),pointer  :: next  ! pointer to the next item
  !   real(dp),pointer             :: v(:)
  !   integer                      :: n=0
  !   real(dp)                     :: r=0
  !   real(dp)                     :: x2=0
  !   real(dp)                     :: x=0
  ! end type
  !
  !public linked_vector,list_vector_grow

!  contains
!
!! Eventos del tipo linked_vector
!
!  subroutine list_vector_grow(lv,p,e)
!    real(dp)                                  :: p(:)
!    real(dp),intent(in)                       :: e
!    type(linked_vector),pointer,intent(inout) :: lv
!    integer                                   :: n,i
!    logical                                   :: found
!    
!    ! A la entrada, la direccion donde comenzar el recorrido de la lista
!    ! A la salida, el objeto agregado
!
!    found=.true.
!    n=size(p)
!    do while (associated(lv))
!      if (n==size(lv%v)) then ! Si los tamaños coinciden 
!        !Me fijo si cada posicion coincide
!        do i = 1,n 
!          if (dabs(lv%v(i)-p(i))>e) then
!            found=.false.
!            exit
!          endif
!        enddo
!        if (found) then ! Si el do anterior se completo
!        !Posiciones iguales, suma a configuracion
!            lv%n=lv%n+1
!            return
!        endif
!      endif
!      lv => lv%next
!    enddo
!   
!    ! Si llegue hasta aca, configuración no encontrada o lista recien empieza
!    allocate(lv) 
!    nullify(lv%next) 
!    allocate(lv%v(n))
!    lv%v=p
!    lv%n=lv%n+1
! 
!  end subroutine list_vector_grow 
!
!! -------------Por desplazamiento  
!  
!  subroutine generate_dispboost(g) 
!    ! Genera el parametro re para cada par de enlaces
!    integer            :: i
!    type(group)        :: g
!    type(atom_dclist),pointer :: la
!    
!    if (allocated(bonds)) deallocate(bonds)
!
!    allocate(bonds(g%nat))
!    nbonds = g%nat
!
!    la => g%alist%next
!    do i = 1, g%nat
!      la%o%pos_eq=la%o%pos
!      bonds(i)%re = 0.0_dp
!      bonds(i)%d = 0.0_dp
!      !bonds(i)%k = sa(la%o%sid,la%o%sid)
!      bonds(i)%i => la%o
!      nullify(bonds(i)%j)
!      la => la%next
!    enddo
!
!    bb=>bonds(1)
!  end subroutine
!
!  subroutine longerdisp_lenght
!    integer             :: k
!    real(dp)            :: ebb
!
!    ! Busco el enlace mas deformado
!    ebb=0.0_dp
!    do k = 1, nbonds
!      bonds(k)%vd = bonds(k)%i%pos-bonds(k)%i%pos_eq
!      bonds(k)%rd = sqrt(dot_product(bonds(k)%vd,bonds(k)%vd))
!      bonds(k)%d = bonds(k)%rd
!      if (bonds(k)%d>ebb) then
!        ebb=bonds(k)%rd
!        bb=>bonds(k)
!      endif
!    enddo
!
!  end subroutine
! 
!  subroutine bondb_disp_lenght(b)
!    type(bondb),intent(inout)   :: b
!
!    b%vd = b%i%pos-b%i%pos_eq
!    b%rd = sqrt(dot_product(b%vd,b%vd))
!    b%d = b%rd
!
!  end subroutine bondb_disp_lenght
 
!! -------------Por enlace
!  
!  subroutine generate_bondsboost(g) 
!    ! Genera el parametro re para cada par de enlaces
!    integer            :: i,j,k
!    real(dp)           :: rd,e,vd(dm)
!    type(group)        :: g
!    type(atom_dclist),pointer :: la,lb
!    
!    if (allocated(bonds)) deallocate(bonds)
!
!    !Cuento cuantos enlaces voy a tener
!    nbonds = 0
!    la => g%alist%next
!    do i = 1, g%nat-1
!      lb => la%next
!      do j = i+1, g%nat
!        vd = vdistance( la%o,lb%o , mic)
!        rd = dot_product(vd,vd)
!        if (rd<rq2(la%o%sid,lb%o%sid)) nbonds = nbonds + 1
!        lb => lb%next
!      enddo
!      la => la%next
!    enddo
!    if (nbonds==0) print *, 'No se encontro ningun enlace busteable'
!
!    ! Allocateo la cantidad de enlaces
!    allocate(bonds(nbonds))
!
!    ! Vuelvo y guardo la info de los enlaces
!    k = 0
!    la => g%alist%next
!    do i = 1, g%nat-1
!      lb => la%next
!      do j = i+1, g%nat
!        vd = vdistance( la%o,lb%o , mic)
!        rd = dot_product(vd,vd) 
!
!        if (rd<rq2(la%o%sid,lb%o%sid)) then
!          k = k + 1 
!          bonds(k)%re = sqrt(rd)
!          bonds(k)%d = 0.0_dp
!          bonds(k)%k = sa(la%o%sid,lb%o%sid)
!          bonds(k)%i => la%o
!          bonds(k)%j => lb%o
!        endif
!
!        lb => lb%next
!      enddo
!      la => la%next
!    enddo
! 
!    bb=>bonds(1)
!
!  end subroutine
!
!  subroutine longerbond_lenght
!    integer             :: k
!    real(dp)            :: ebb
!
!    ! Busco el enlace mas deformado
!    ebb=0.0_dp
!    do k = 1, nbonds
!      bonds(k)%vd = vdistance(bonds(k)%j,bonds(k)%i, mic)
!      bonds(k)%rd = sqrt(dot_product(bonds(k)%vd,bonds(k)%vd))
!      bonds(k)%d = dabs(bonds(k)%rd-bonds(k)%re)/bonds(k)%re
!      if (bonds(k)%d>ebb) then
!        ebb=bonds(k)%d
!        bb=>bonds(k)
!      endif
!    enddo
!    
!
!
!  end subroutine
!

! ------ Otras cosas
!
!  subroutine shorterbond_energy
!    integer             :: k
!    real(dp)            :: boost
!
!    ! Busco el enlace menos boosteable
!    boost=0.0_dp
!    do k = 1, nbonds
!      bonds(k)%vd = vdistance(bonds(k)%j,bonds(k)%i, mic)
!      bonds(k)%rd = sqrt(dot_product(bonds(k)%vd,bonds(k)%vd))
!      bonds(k)%d = dabs(bonds(k)%rd-bonds(k)%re)/bonds(k)%re
!      bonds(k)%e = bonds(k)%k*(1-(bonds(k)%d*bonds(k)%d/eq2))
!      if (boost<bonds(k)%e) then
!        ! Calculo el boost del enlace
!        boost=bonds(k)%e
!        bb=>bonds(k)
!      endif
!    enddo
!    
!  end subroutine
!
!  subroutine simplified_bondboost
!    real(dp)            :: factor2(dm)   
!
!    call shorterbond_energy
!
!   ! if (mod(dm_steps,outeach)==0) then
!   !   write(9,*) 2
!   !   write(9,*) 
!   !   write(9,'(a2,2x,3(f20.15))') 'Au',b%i%pos
!   !   write(9,'(a2,2x,3(f20.15))') 'Au',b%j%pos
!   ! endif
!
!    ! force calculation (restricted)
!    factor2 = bb%k/eq2*bb%d *bb%vd(1:dm)/bb%rd
!    bb%i%force(1:dm) = bb%i%force(1:dm) + factor2
!    bb%j%force(1:dm) = bb%j%force(1:dm) - factor2
!     
!    ! energy calculation
!    bb%i%epot = bb%i%epot + bb%e*0.5
!    bb%j%epot = bb%j%epot + bb%e*0.5
!
!    boost=1.0e16_dp 
!  end subroutine
!
end module gems_bondboost
