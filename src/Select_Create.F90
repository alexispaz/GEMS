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

 
module gems_select_create
 use gems_errors
 use gems_groups, only: group,group_allatom_del,group_atom_add,dgroup_atom_add
 use gems_atoms, only: atom,atom_dclist,atom_setelmnt,atom_asign
 use gems_program_types, only: natoms
 use gems_constants,only:dp,grad_rad,find_io,dm
 use gems_set_properties
 use gems_inq_properties

  implicit none

  public
  private :: group,atom,atom_dclist ! program_types
  private :: group_allatom_del,group_atom_add ! program_types
  private :: atom_setelmnt     ! program_types
  private :: dm,natoms                                     ! program_types
  private :: dp,grad_rad                                  ! constants

 contains

!                                                    create
!----------------------------------------------------------

! Los comandos creacion, agregan atomos a un grupo dummy llamado gnew, con
! la particularidad de que los nodos de esta lista estan allocateados. La
! subrutina que realiza este tipo de crecimiento de la lista es la
! dgroup_atom_add.A su vez agrega los nuevos atomos al  grupo dummy gout.

subroutine create_atom(ctr,gnew,gout)
  type(group),intent(inout)         :: gnew
  type(group),intent(inout),optional  :: gout
  real(dp),intent(in)    :: ctr(dm)
  type(atom)             :: a
  type(atom),pointer     :: ptr_a

  call  a%init()
  a%pos = ctr
  call dgroup_atom_add( a , gnew )

  ptr_a => gnew%alist%prev%o
  if(present(gout)) call group_atom_add( ptr_a , gout )
  call  a%dest()

end subroutine create_atom
 
subroutine create_fill(sigma,num,gnew,opt_frame,opt_origin,gout,opt_pbc)
  ! Fill simulation box (or frame) with atoms without superpose them. The
  ! criteria for superposition is based in the sigma (`s`) property of the
  ! atoms. Each created atoms is instantanouesly set with a sigma equal to
  ! the input `sigma`.
  use gems_random,only: ranu
  use gems_program_types, only: box
  type(group),intent(inout)           :: gnew
  type(group),intent(inout),optional  :: gout
  real(dp),intent(in)                 :: sigma
  logical,intent(in),optional         :: opt_pbc
  real(dp),intent(in),optional        :: opt_frame(1:dm),opt_origin(1:dm)
  real(dp)                            :: frame(1:dm),origin(1:dm)
  real(dp)                            :: cut,rd
  integer,intent(in)                  :: num
  type(atom)                          :: a
  type(atom),pointer                  :: ptr_a
  type(atom_dclist),pointer           :: la
  integer                             :: i,j,k,n
  integer,parameter                   :: nmax=1000

           
  if(.not.present(opt_origin)) then
    origin(1:dm)=0.0_dp
  else
    origin(1:dm)=opt_origin(1:dm)
  endif
   
  if(.not.present(opt_frame)) then
    frame(1:dm)=box(1:dm)
  else
    frame(1:dm)=opt_frame(1:dm)
  endif

  call  a%init()
  a%s = sigma
 
  if(.not.present(opt_pbc)) then
    a%pbc(1:dm)=.false.
  else
    a%pbc(1:dm)=opt_pbc
  endif
            
  do i = 1,num

    try: do n = 1,nmax
      do j = 1,dm
        a%pos(j) = ranu()*frame(j)+origin(j)
      enddo
 
      la => gnew%alist 
      do k = 1,gnew%nat
       la=>la%next
       rd =  rdistance2( a, la%o )
       cut = a%s+la%o%s
       cut = cut*cut
       if (rd<cut) cycle try
      enddo     
              
      exit
    enddo try
    call werr('Max iteration reached, no space was found.',n>nmax)

    ! Add atom
    call dgroup_atom_add( a , gnew )
    ptr_a => gnew%alist%prev%o
    if(present(gout)) call group_atom_add( ptr_a , gout )
  enddo

  call  a% dest()

end subroutine create_fill
 
subroutine create_reply(gini,tras,times,gnew,gout)
  type(group),intent(inout)          :: gnew
  type(group),intent(inout),optional :: gout
  real(dp),intent(in)                :: tras(dm)
  integer,intent(in)                 :: times
  integer                            :: i,j
  type(atom)                         :: a
  type (atom_dclist),pointer         :: la
  type(group),intent(in)             :: gini
  type(atom),pointer                 :: ptr_a

  call  a%init()

  la => gini%alist%next
  do i=1,gini%nat
    ! a=la%o
    call atom_asign(a,la%o) 
    do j=1,times-1
      a%pos = la%o%pos + tras*j
      call dgroup_atom_add( a , gnew )
      ptr_a => gnew%alist%prev%o
      if(present(gout)) call group_atom_add( ptr_a , gout )
    enddo
    la=>la%next
  enddo

  call a%dest()

end subroutine create_reply

subroutine create_file(archivo,ext,gnew,gout,frame)
  use gems_elements,only:ncsym
  ! read atoms from file
  type(group),intent(inout)          :: gnew
  type(group),intent(inout),optional :: gout 
  integer                            :: i,j,u,io
  character(ncsym)                   :: sym
  character(3),intent(in)            :: ext
  character(*),intent(in)            :: archivo
  character(80)                      :: pdbword
  integer,intent(in)                 :: frame
  type(atom)                         :: a
  type(atom),pointer                 :: ptr_a

  call  a%init()

  u = find_io(30)
  open(u,action='read',file=archivo)

  select case(ext)
  case ('xyz')

    ! In case a specific frame wants to be loaded
    do i=1,frame-1
      read(u,*) j
      read(u,*)
      do j = 1, j
        read(u,*)
      enddo
    enddo

    ! Reading the frame
    read(u,*) j
    read(u,*)
    do i = 1, j
      read(u,*) sym, a%pos
      call atom_setelmnt(a,sym)
      
      call dgroup_atom_add( a , gnew )
      ptr_a => gnew%alist%prev%o
      if(present(gout)) call group_atom_add( ptr_a , gout )
    enddo

  case ('pdb')

    io=0
    do while(io==0)
      read(u,fmt='(a80)') pdbword

      select case(pdbword(1:6))
      case ('ATOM  ','HETATM')

        ! Set the atom element
        call atom_setelmnt(a,pdbword(13:16))

        read(pdbword(31:54),fmt='(3(f8.3))') a%pos
        !read(pdbword(55:60),fmt='(3(f8.3))') a%q
        !read(pdbword(23:26),fmt='(3(f8.3))') resid

        call dgroup_atom_add( a , gnew )
        ptr_a => gnew%alist%prev%o
        if(present(gout)) call group_atom_add( ptr_a , gout )
      case ('REMARK')
      case ('END   ')
        exit
      case default
        call wwan(); write(logunit,*) 'PDB keyword ',trim(pdbword),' is not supported'
      end select

    enddo 

  case default
    call werr(); write(logunit,*) 'File format ',trim(ext),' unknown'
  end select

  call a%dest()
  close(u)

end subroutine create_file
 
subroutine create_fillh(gini,gnew,gout)
  ! Agrega hidrogenos en atomos de Carbono
  use gems_constants,only:rad90,rad60,rad109
  type(group),intent(in)             :: gini
  type(group),intent(inout)          :: gnew
  type(group),intent(inout),optional :: gout 
  real(dp)                           :: rd
  real(dp),parameter                 :: chd=1.1_dp       ! Mas o menos la distancia C-H
  real(dp),parameter                 :: rc=1.2_dp*1.4_dp ! Mas o menos un radio de corte
  real(dp),parameter                 :: rc2=rc*rc       
  integer                            :: i,j,n,m
  type (atom_dclist),pointer         :: la,lb
  type(atom)                         :: b(4)
  real(dp),dimension(dm)             :: vd,vaux,vaux2
  type(atom),pointer                 :: ptr_a
     
  call b(1)%init()
  call b(2)%init()
  call b(3)%init()
  call b(4)%init()
       
  la => gini%alist
  do i=1,gini%nat
    la=>la%next

    if(la%o%z/=6) cycle

    ! Si no tiene hybrydizacion establecida aviso
    if (la%o%sp==0) then
      call wwan();write(logunit,'(a,i0,a)') 'El atomo ',la%o%id,' no tiene seteado el sp'
      cycle
    endif
    
    ! Busco los vecinos que ya existen
    n=0
    lb => gini%alist 
    do j = 1,gini%nat
      lb=>lb%next
      
      if(associated(lb%o,target=la%o)) cycle
      
      rd = rdistance2( lb%o , la%o )
      
      if(rd>rc2) cycle
      
      n=n+1
      
      if (n>4)  exit
      
      ! b(n)=lb%o
      call atom_asign(b(n),lb%o)
    enddo     
 

    ! Guardo el numero de los que faltan
    m=la%o%sp-n+1
 
    ! Si ya esta completo ciclo
    if (m==0) then

      cycle
    endif

    ! Si esta mas que completo, advierto
    if (m<0) then
      call wwan();write(logunit,'(a,i0,a)') 'El atomo ',la%o%id,' tiene mas de 4 vecinos'
      cycle
    endif

    ! Caso particular: Si faltan 3 agrego un atomo y paso al siguiente if
    if(m==3) then ! si faltan 3 

      ! Versor en la direccion al centro
      vd = vrs( la%o, b(1) )

      ! Obtengo un versor ortogonal
      vaux=ortogonal(vd)
      vaux=vaux/sqrt(dot_product(vaux,vaux))

      ! lo giro algun angulo 
      b(2)%pos(:) = chd*rodrigues_rotation(rad109,vaux,vd)+la%o%pos(:)

      ! Agrego un atomo
      call dgroup_atom_add( b(2) , gnew )
      ptr_a => gnew%alist%prev%o
      call ptr_a%setz(1)
      if(present(gout)) call group_atom_add( ptr_a , gout )

      n=n+1
      m=m-1

    endif
 

    if(m==2) then ! si faltan 2

      if(la%o%sp==3) then

        ! replico un atomo a su punto opuesto de i
        b(3)%pos = la%o%pos+vrs( la%o, b(1) )*chd

        ! replico el otro atomo a su punto opuesto de i
        b(4)%pos = la%o%pos+vrs( la%o, b(2) )*chd

        ! Consigo un eje que pasa por el punto medio de los dos nuevos atomos
        vaux2 = b(3)%pos+(b(4)%pos-b(3)%pos)*0.5 - la%o%pos

        ! Ahora roto estos dos 90 grados para formar el tetrahedro
        vaux = b(3)%pos-la%o%pos
        b(3)%pos = rodrigues_rotation(rad90,vaux2,vaux)+la%o%pos
        vaux = b(4)%pos-la%o%pos
        b(4)%pos = rodrigues_rotation(rad90,vaux2,vaux)+la%o%pos
        
      else     ! sp==2

        ! vector que apunta al atomo central
        vaux = vrs( la%o, b(1) )*chd


        ! Obtengo un versor ortogonal
        vaux2=ortogonal(vaux)
        vaux2=vaux2/sqrt(dot_product(vaux2,vaux2))
 
        ! lo giro algun angulo 
        b(2)%pos = rodrigues_rotation(rad60,vaux2,vaux)+la%o%pos
        b(3)%pos = rodrigues_rotation(-rad60,vaux2,vaux)+la%o%pos

      endif

    else  ! si falta 1

      if(la%o%sp==3) then
        b(4)%pos = vrs(la%o,b(1)) + vrs(la%o,b(2)) + vrs(la%o,b(3))
        b(4)%pos = chd*b(4)%pos/sqrt(dot_product(b(4)%pos,b(4)%pos))
        b(4)%pos = b(4)%pos+la%o%pos
      else if(la%o%sp==2) then 
        b(3)%pos = vrs(la%o,b(1)) + vrs(la%o,b(2))
        b(3)%pos = chd*b(3)%pos/sqrt(dot_product(b(3)%pos,b(3)%pos))
        b(3)%pos = b(3)%pos+la%o%pos
      else   ! sp==1
        b(2)%pos = chd*vrs(la%o,b(1))+la%o%pos
      endif
 
    endif
    
    do j=n+1,n+m
      call dgroup_atom_add( b(j) , gnew )
      ptr_a => gnew%alist%prev%o
      call ptr_a%setz(1)
      if(present(gout)) call group_atom_add( ptr_a , gout )
    enddo

  enddo
 
  call b(1)%dest()
  call b(2)%dest()
  call b(3)%dest()
  call b(4)%dest()
 
end subroutine create_fillh

 
!                                                    select
!----------------------------------------------------------

! Selecciones de grupos en grupos

! Los procesos de seleccion, seleccionan del grupo gini aquellos nodos que
! cumplan con una determinada condicion. Estos nodos constituyen la salida de la
! subrutina a travez del grupo gout

subroutine select_sphere(rad,ctr,gini,gout)
  ! recorta una esfera en el centro del groupo
  ! y lo suma a gini
  type(group),intent(in)     :: gini
  type(group),intent(inout)  :: gout
  real(dp)                   :: rd,rad,ctr(dm)
  type (atom_dclist),pointer :: la
  integer                    :: i

  call inq_pos_v(gini,ctr)
  
  rad = rad*rad

  la => gini%alist%next
  do i = 1,gini%nat
    rd = dot_product(la%o % pos_v,la%o % pos_v)
    if (rd < rad ) call group_atom_add(la%o,gout)
    la => la%next
  enddo

end subroutine select_sphere 

subroutine select_cylinder(rad,axis,ctr,gini,gout)
  ! recorta una esfera en el centro del groupo
  ! y lo suma a gini
  type(group),intent(in)     :: gini
  type(group),intent(inout)  :: gout
  real(dp)                   :: rd,rad,ctr(dm),axis(dm),aux(dm)
  type (atom_dclist),pointer :: la
  integer                    :: i

  call inq_pos_v(gini,ctr)
  
  rd=sqrt(dot_product(axis,axis))
  axis=axis/rd

  rad = rad*rad

  la => gini%alist%next
  do i = 1,gini%nat
    aux = la%o%pos_v-dot_product(la%o % pos_v,axis)*axis
    rd=dot_product(aux,aux)
    if (rd < rad ) call group_atom_add(la%o,gout)
    la => la%next
  enddo

end subroutine select_cylinder

subroutine select_near(ctr,gini,gout)
  ! selecciona el atomo mas cercano al punto
  type(group),intent(in)     :: gini
  type(group),intent(inout)  :: gout
  real(dp)                   :: rd,rad,ctr(dm)
  type (atom_dclist),pointer :: la
  integer                    :: i,j=0

  call inq_pos_v(gini,ctr)
  
  rad = 1.0e8_dp

  la => gini%alist%next
  do i = 1,gini%nat
    rd = dot_product(la%o % pos_v,la%o % pos_v)
    if (rd < rad ) then
      rad=rd
      j=i
    endif
    la => la%next
  enddo
  if (j == 0) return
  
  la => gini%alist%next
  do i = 1,j-1
    la => la%next
  enddo 

  call group_atom_add(la%o,gout)

end subroutine select_near

subroutine select_far(ctr,gini,gout)
  ! selecciona el atomo mas cercano al punto
  type(group),intent(in)     :: gini
  type(group),intent(inout)  :: gout
  real(dp)                   :: ctr(dm)

  call group_atom_add(inq_far(ctr,gini),gout)

end subroutine select_far

subroutine select_band(d,mn,mx,gini,gout)
  ! recorta una esfera en el centro del groupo
  ! y lo suma a gini
  type(group),intent(in)     :: gini
  type(group),intent(inout)  :: gout
  integer,intent(in)         :: d
  real(dp)                   :: mn,mx,aux
  type (atom_dclist),pointer :: la
  integer                    :: i

  if (mn>mx) then
    aux=mx
    mx=mn
    mn=aux
  endif
  
  la => gini%alist%next
  do i = 1,gini%nat
    if (mn>=la%o%pos(d)) goto 10
    if (la%o%pos(d)>=mx) goto 10
    call group_atom_add(la%o,gout)
10  continue    
    la => la%next
  enddo

end subroutine select_band
 
subroutine select_rectangle(vert1,vert2,gini,gout)
  ! recorta una esfera en el centro del groupo
  ! y lo suma a gini
  type(group),intent(in)     :: gini
  type(group),intent(inout)  :: gout
  real(dp)                   :: vert1(dm),vert2(dm)
  type (atom_dclist),pointer :: la
  integer                    :: i,j

  la => gini%alist
  do i = 1,gini%nat
    la => la%next
    do j = 1,dm
      if (la%o%pos(j)< vert1(j) .or. vert2(j)< la%o%pos(j)  ) goto 10
    enddo
    call group_atom_add(la%o,gout)
10  continue    
  enddo
  
end subroutine select_rectangle
 
subroutine select_element(z,gini,gout)
  ! selecciona el atomo mas cercano al punto
  type(group),intent(in)     :: gini
  type(group),intent(inout)  :: gout
  type (atom_dclist),pointer :: la
  integer,intent(in)         :: z
  integer                    :: i

  la => gini%alist%next
  do i = 1,gini%nat
    if (z == la%o%z ) then
      call group_atom_add(la%o,gout)
    endif
    la => la%next
  enddo

end subroutine select_element

  subroutine select_molid(id,gini,gout)
  type(group),intent(in)     :: gini
  type(group),intent(inout)  :: gout
  type (atom_dclist),pointer :: la
  integer                    :: i,id

  la => gini%alist%next
  do i = 1,natoms
    if (id == la%o%molid ) then
      call group_atom_add(la%o,gout)
    endif
    la => la%next
  enddo

  end subroutine select_molid

subroutine select_below(p1,p2,p3,gini,gout)
  ! Selecciona un semiplano despues de dividir el espacio con un plano.  
  ! Si el plano no corta el origen, Se queda con el semiespacio que contenga el
  ! origen. El plano se lo recibe como 3 puntos no colineales en el espacio y se
  ! calcula el vector normal N1 y con el la forma canonica del tipo X.N=p1.N. Si
  ! el producto de la posicion de un atomo r.N1<p1.N1, se acepta dicho atomo
  ! dentro del semiespacio. Esto es lo mismo que (r-p1).N1<0.  Dado que hay dos
  ! vectores N, elegir N1 o N2 determina el semiespacio elegido. Para elegir que
  ! sea el semiespacio que contenga el origen, se pide que (0-p1).N1<0, o lo que
  ! es lo mismo, p1.N>o, osea que N1 no apunta al origen.  
  ! Si el plano corta el origen.... no se que hacer
  type(group),intent(in)     :: gini
  type(group),intent(inout)  :: gout
  real(dp),dimension(dm)     :: p1,p2,p3,n,aux,aux2
  type (atom_dclist),pointer :: la
  integer                    :: i
#ifdef DIM3

  ! Hago las restas
  aux=p2-p1
  n=p3-p1

  ! Obtengo algun vector normal
  aux2=cross_product(aux,n)
  
  ! Elijo que sea el de la famila que no apunta al origen
  if(dot_product(p1,aux2)<0) then
    n=-aux2
  else
    n=aux2
  endif

  la => gini%alist%next
  do i = 1,gini%nat
    aux=la%o%pos-p1
    if (dot_product(aux,n)<=0) then
      call group_atom_add(la%o,gout)
    endif
    la => la%next
  enddo

#endif
end subroutine select_below

subroutine select_above(p1,p2,p3,gini,gout)
  ! Idem select_below pero usando el otro semiplano
  type(group),intent(in)     :: gini
  type(group),intent(inout)  :: gout
  real(dp),dimension(dm)     :: p1,p2,p3,n,aux,aux2
  type (atom_dclist),pointer :: la
  integer                    :: i
#ifdef DIM3

  ! Hago las restas
  aux=p2-p1
  n=p3-p1

  ! Obtengo algun vector normal
  aux2=cross_product(aux,n)
  
  ! Elijo que sea el de la famila que no apunta al origen
  if(dot_product(p1,aux2)<0) then
    n=-aux2
  else
    n=aux2
  endif

  la => gini%alist%next
  do i = 1,gini%nat
    aux=la%o%pos-p1
    if (dot_product(aux,n)>=0) then
      call group_atom_add(la%o,gout)
    endif
    la => la%next
  enddo

#endif
end subroutine select_above

subroutine select_bond_order_range(rini,rfin,r1,r2,gini,g2,gout)
  ! This will select the atoms of "gini" that have a bond order number respect
  ! to atoms of group "g2" that fall inside the range defined by "rini"--"rfin".
  ! All the neighbour inside a radious of "r1" contribute with 1 to the bond
  ! order, after that there is a shell of lenght "r2"-"r1" that smooth the
  ! neighbour contribution.
  type(group),intent(inout)    :: gini,g2
  type(group),intent(inout)    :: gout
  real(dp),intent(in)          :: rini,rfin,r1,r2
  type (atom_dclist),pointer   :: la
  integer                      :: i


  ! Compute the bond order of each atom
  call bond_order_groups(gini,g2,r1,r2)

  ! Select the atoms with de desired bond order range
  la => gini%alist
  do i = 1,gini%nat
    la => la%next
    if(la%o%border>=rini.and.la%o%border<=rfin) then
      call group_atom_add(la%o,gout)
    endif
  enddo

end subroutine select_bond_order_range

subroutine select_upper_bond_order(r1,r2,gini,g2,gout)
  ! This will select the atom of "gini" that have the lower bond order number
  ! respect to atoms of group "g2".  All the neighbour inside a radious of "r1"
  ! contribute with 1 to the bond order, after that there is a shell of lenght
  ! "r2"-"r1" that smooth the neighbour contribution.
  type(group),intent(inout)    :: gini,g2
  type(group),intent(inout)    :: gout
  real(dp),intent(in)          :: r1,r2
  real(dp)                     :: aux
  type (atom_dclist),pointer   :: la
  type (atom),pointer          :: pa
  integer                      :: i


  ! Compute the bond order of each atom
  call bond_order_groups(gini,g2,r1,r2)

  ! Serch for the lowest bond order atom "pa"
  aux=-1e8
  la => gini%alist
  do i = 1,gini%nat
    la => la%next
    if(la%o%border>aux) then
      pa=>la%o
      aux=la%o%border
    endif
  enddo

  ! Selection
  call group_atom_add(pa,gout)

end subroutine select_upper_bond_order
         
subroutine select_lower_bond_order(r1,r2,gini,g2,gout)
  ! This will select the atom of "gini" that have the lower bond order number
  ! respect to atoms of group "g2".  All the neighbour inside a radious of "r1"
  ! contribute with 1 to the bond order, after that there is a shell of lenght
  ! "r2"-"r1" that smooth the neighbour contribution.
  type(group),intent(inout)    :: gini,g2
  type(group),intent(inout)    :: gout
  real(dp),intent(in)          :: r1,r2
  real(dp)                     :: aux
  type (atom_dclist),pointer   :: la
  type (atom),pointer          :: pa
  integer                      :: i


  ! Compute the bond order of each atom
  call bond_order_groups(gini,g2,r1,r2)

  ! Serch for the lowest bond order atom "pa"
  aux=1e8
  la => gini%alist
  do i = 1,gini%nat
    la => la%next
    if(la%o%border<aux) then
      pa=>la%o
      aux=la%o%border
    endif
  enddo

  ! Selection
  call group_atom_add(pa,gout)

end subroutine select_lower_bond_order
      
subroutine select_atom(id,gini,gout)
! selecciona el atomo con indice id
  type(group),intent(in)     :: gini
  type(group),intent(inout)  :: gout
  type (atom_dclist),pointer :: la
  integer                    :: i,id

  la => gini%alist%next
  do i = 1,gini%nat
    if(id==la%o%id) then
      call group_atom_add(la%o,gout)
      return
    endif    
    la => la%next
  enddo


end subroutine select_atom

subroutine select_random(n,gini,gout)
  use gems_random, only: ranu
! selecciona n atomos random
  type(group),intent(in)     :: gini
  type(group),intent(inout)  :: gout
  type (atom_dclist),pointer :: la
  integer                    :: i,k,n
  logical,allocatable        :: picked(:)

  allocate(picked(gini%nat))
  picked(:)=.false.

  call werr('Not enough atoms', gini%nat<n)

  do while (gout%nat < n)

    k=int(gini%nat*ranu())+1

    if(picked(k)) cycle
    picked(k)=.true.
    
    la => gini%alist
    do i = 1,gini%nat
      la => la%next
      if(k==i) then
        call group_atom_add(la%o,gout)
        exit
      endif    
    enddo
  enddo

  deallocate(picked)

end subroutine select_random

! Selecciones de grupos en atomos

subroutine select_covbond(a,gini,gout)
  use gems_elements,only:elements
  ! busca en gini el grupo de atomos enlazados covalentemente al atomo a
  type(atom),intent(in)      :: a
  type(group),intent(in)     :: gini
  type(group),intent(inout)  :: gout
  real(dp)                   :: rd,cov
  type(atom_dclist),pointer  :: la
  integer                    :: i

  cov=elements%o(a%z)%rad

  la => gini%alist
  do i = 1,gini%nat
    la => la%next
    rd = rdistance(a,la%o)
    if(rd<=cov.and.a%id/=la%o%id) call group_atom_add(la%o,gout)
  enddo


end subroutine select_covbond

! subroutine vatom_addatom(v,a)
!   ! Cuidado: devuelve un vector 1 elemento mas grande
!   type(atom),intent(in)                :: a
!   type(atom),allocatable,intent(inout) :: v(:)
!   type(atom),allocatable               :: vaux(:)
!   integer                              :: n,i
!
!   if (allocated(v)) then
!     n = size(v)
!     allocate(vaux(n))
!     do i =1,n
!       call  vaux(i)% init()
!     enddo
!     deallocate(v)
!   else
!     n = 0
!     allocate(vaux(1))
!     call  vaux(1)% init()
!   endif
!
!   allocate(v(n+1))
!   do i =1,n+1
!     call  v(i)% init()
!   enddo
!   
!   v(1:n)=vaux
!   v(n+1)=a
!  
!   deallocate(vaux)
!
! end subroutine vatom_addatom
 
end module gems_select_create
