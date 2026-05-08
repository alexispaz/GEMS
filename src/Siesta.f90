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

 
module gems_siesta
  use gems_program_types
  use gems_input_parsing
  use gems_constants
  use gems_neighbor
  use gems_errors

  integer              :: igr_siesta=0

  integer              :: siesta_esp=0
  integer,parameter    :: siesta_mesp=10
 
  ! interaccion-----
  type(group)          :: gsiesta
  real(dp)             :: siesta_cut
  !-----------------
  
  ! Unidad que se conecta a los distintos archivos
  integer              :: u_siesta

  ! Aca alocateo con msubs, dado que el numero maximo de especies posibles esta
  ! dado por el numero de pseudopotenciales, el cual es uno por cada subsistema
  character(90)        :: pseudo(siesta_mesp)
  integer              :: esp_ss(siesta_mesp)  ! Por cada especia da el subsitema

  character(linewidth)   :: call1,call2,call3,call4
  character(90)          :: siesta_ini, siesta_log, siesta_tmp, siesta_fdf

  public        :: siesta_fdf

contains

  subroutine siesta_writehead
    integer              :: i,siesta_nat   
    character(linewidth) :: temp

    siesta_nat=ngindex%o(igr_siesta)%o%n(1)

    ! Los archivos involucrados en el llamado al siesta
    write(siesta_ini,*) trim(adjustl(ioprefix))//"_siesta.fdf"  ! El input... se genera automaticamente
    write(siesta_log,*) trim(adjustl(ioprefix))//"_siesta.log"  ! La salida... se genera automaticamente
    write(siesta_tmp,*) trim(adjustl(ioprefix))//"_siesta.tmp"  ! Para intercambiar info entre el siesta y el gems
    ! Es importante que estos archivos lleven el prefijo para que no se mezcle
    ! si corren programas simultaneos en la misma carpeta
 
    
    ! Los calls que se van a hacer para interaccionar con el siesta. TODO: Pasar
    ! a lecturas inrinsecas sin depender de callsistems
    write(call1,'(a)') "./siesta < " // trim(adjustl(siesta_ini)) //" > "// trim(adjustl(siesta_log))   ! Correr el siesta
    write(call2,'(a)') "grep 'Total =' "// trim(adjustl(siesta_log)) //" | tail -n 1 | awk '{print $4}'> " &
      // trim(adjustl(siesta_tmp)) ! Obtener la energia del log
    write(call3,'(a,i0,x,a,i0,a)') "grep 'Atomic forces' -A",siesta_nat, trim(adjustl(siesta_log)) &
      //" |   tail -n ",siesta_nat," | sed 's/siesta://g' >> "// trim(adjustl(siesta_tmp)) ! Obtener la fuerza del log
    write(call4,'(a)') "grep -A 3 'siesta: Stress tensor' " // trim(adjustl(siesta_log)) &
      //" | tail -n 3 >> "// trim(adjustl(siesta_tmp)) ! Obtener el estres

    ! TODO: en el call4 se podria poner un call que extraiga un rango de lineas
    ! de la salida del sisesta, asi seria independiente de que uno quiere
    ! extraer y podria ser general
 
    call wlog('SIESTA', 'Run data is set here, any other change will be ignored')
    call wlog('SIESTA','')
    call wlog('SIESTA','The command to run in each step will be:')
    call wlog('SIESTA',call1)
    call wlog('SIESTA','')
    call wlog('SIESTA', trim(adjustl(siesta_ini)) //  ' contain:')
    
    u_siesta=find_io(30)
    open( u_siesta , file=trim(adjustl(siesta_ini)) )
    ! Nombre del trabajo 
    write(temp,*) 'SystemName '//trim(adjustl(ioprefix))//'_siesta' ; write(u_siesta,*) trim(temp); call wlog('SIESTA',temp) 
    write(temp,*) 'SystemLabel '//trim(adjustl(ioprefix))//'_siesta'; write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    ! Cosas generales (bases y demas) en archivo fdf suministrado
    write(temp,*) '%include '//trim(adjustl(siesta_fdf)); write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    ! Especies quimicas y elementos
    write(temp,*) 'NumberOfAtoms ',siesta_nat   ; write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    write(temp,*) 'NumberOfSpecies ',siesta_esp ; write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    write(temp,*) '%block ChemicalSpeciesLabel' ; write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    do i = 1,siesta_esp
      write(temp,*) i, esp_ss(i),' ',trim(adjustl(pseudo(i)))
      write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    enddo
    write(temp,*) '%endblock ChemicalSpeciesLabel'; write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    ! La caja
    write(temp,*) '%block LatticeVectors'   ;write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    write(temp,*) box(1),' 0.00 0.00'       ;write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    write(temp,*) '0.00 ',box(2),' 0.00'    ;write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    write(temp,*) '0.00 0.00 ',box(3)       ;write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    write(temp,*) '%endblock LatticeVectors';write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    ! De donde lee las coordenadas atomicas
    write(temp,*) '%block AtomicCoordinatesAndAtomicSpecies < '//trim(adjustl(siesta_tmp))
    write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    ! Que me escriba las fuerzas
    write(temp,*) 'WriteForces T'; write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    ! Unidades
    write(temp,*) 'AtomCoorFormatOut Ang'      ;write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    write(temp,*) 'AtomicCoordinatesFormat Ang';write(u_siesta,*)trim(temp); call wlog('SIESTA',temp) 
    close(u_siesta)

    ! Uso de una densidad anteriormente convergida
    write(temp,*) 'DM.UseSaveDM T';call wlog('SIESTA',temp) 

    ! Corriendo sin la linea anterior para encontrar la primera densidad
    call wlog('SIESTA','')
    call wlog('SIESTA', 'Get .DM file with the actual configuration')
    call ngindex%o(igr_siesta)%o%interact()

    ! Escribiendo nuevamente con el uso de la densidad anterior
    u_siesta=find_io(30)
    open( u_siesta , file=trim(adjustl(siesta_ini)) )
    write(u_siesta,*) 'SystemName '//trim(adjustl(ioprefix))//'_siesta' 
    write(u_siesta,*) 'SystemLabel '//trim(adjustl(ioprefix))//'_siesta'
    write(u_siesta,*) '%include '//trim(adjustl(siesta_fdf))
    write(u_siesta,*) 'NumberOfAtoms ',siesta_nat   
    write(u_siesta,*) 'NumberOfSpecies ',siesta_esp 
    write(u_siesta,*) '%block ChemicalSpeciesLabel' 
    do i = 1,siesta_esp
      write(u_siesta,*) i, esp_ss(i),' ',trim(adjustl(pseudo(i)))
    enddo
    write(u_siesta,*) '%endblock ChemicalSpeciesLabel'
    write(u_siesta,*) '%block LatticeVectors'   
    write(u_siesta,*) box(1),' 0.00 0.00'       
    write(u_siesta,*) '0.00 ',box(2),' 0.00'    
    write(u_siesta,*) '0.00 0.00 ',box(3)       
    write(u_siesta,*) '%endblock LatticeVectors'
    write(u_siesta,*) '%block AtomicCoordinatesAndAtomicSpecies < '//trim(adjustl(siesta_tmp))
    write(u_siesta,*) 'WriteForces T'
    write(u_siesta,*) 'AtomCoorFormatOut Ang'      
    write(u_siesta,*) 'AtomicCoordinatesFormat Ang'
    write(u_siesta,*) 'DM.UseSaveDM T'
    close(u_siesta)
                               
  end subroutine siesta_writehead

  subroutine siesta_interaction(ig)
    class(ngroup),intent(inout)         :: ig
    integer                   :: i,k,l,j
    character(150)            :: w1
    type(atom_dclist),pointer :: la

    ig%epot=0.0_dp


    if(ig%nat==0) return

    ! Escribo coordenadas
    u_siesta=find_io(30)
    open( u_siesta , file=trim(adjustl(siesta_tmp)) )

    la => ig%alist      ! sobre los atomos
    do i = 1,ig%nat
      la=>la%next
 
      do k = 1,siesta_esp
        if(esp_ss(k)/=la%o%z) cycle
        write(u_siesta,'(3(f20.12,2x),i0)') (la%o%pos(l),l=1,dm),k
        exit
      enddo

    enddo
    close( u_siesta )

    ! Corro siesta
    call system(call1)

    ! Filtro la salida para obtener energía y fuerza
    call system(call2) ! Energia
    call system(call3) ! Fuerza
    call system(call4) ! Estress

    u_siesta=find_io(30)
    open( u_siesta , file=trim(adjustl(siesta_tmp)) )

    ! Leo la energía
    read(u_siesta,*) ig%epot
    ig%epot=ig%epot*ev_ui

    ! Leo las fuerzas (que estan en ev/a)
    l=0
    la => ig%alist      ! sobre los atomos
    do i = 1,ig%nat
      la=>la%next
      read(u_siesta,*) j,(la%o%force(l),l=1,dm)
      la%o%force = la%o%force*ev_ui 
      la%o%epot = ig%epot/ig%nat  ! POR COMPATIBLIDAD
    enddo

    ! Tensor de estress
    do i = 1,3
      read(u_siesta,'(a150)')  w1
      call wstd(); write(logunit,'(a150)')  w1
    enddo 

    close( u_siesta )

  end subroutine siesta_interaction

end module gems_siesta
