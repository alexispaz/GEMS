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

 
module gems_checkpoint
use gems_errors
use gems_program_types
use gems_constants, only: dm, linewidth
use gems_input_parsing
use gems_random, only: write_chpseed, read_chpseed
use gems_integration, only: write_chppiston, read_chppiston
use gems_interaction, only: interact
use gems_set_properties, only: posvel_changed
use gems_output, only: chpmode
#ifdef HAVE_MPI
 use mpi_f08
#endif
implicit none

integer               :: current=0

logical               :: b_ckp=.false.     ! Enabled chelpoint
integer               :: chpeach=1000000   ! Resolucion de escritura
character(linewidth)  :: chpfile=''        ! nombre del archivo

! Dado que se cierra y abre el archivo muchas veces, y para evitar buscar
! siempre una unidad libre, aca se asigna de forma manual una unidad que este
! siempre reservada para el checkpoint. Por eso debe ser inferior a 30, ya que
! como convencion, cuando se abre un archivo en el programa GEMS, siempre se usa
! la funcion find_io(30) que busca unidades libres por arriba de 30, cosa de
! permitir unidades reservadas por debajo de 30.
integer               :: chpunit=25      
contains

subroutine readwrite_chp(step,nsteps,control)
! Cada vez que esta subrrutina es llamada es porque hay que guardar, leer o
! acumular un checkpoint. Antes de cada algoritmo suseptible de tener un
! checkpoint se debe llamar.
! En modo lecutra nsteps es el numero total de pasos a
! ejecutar en el algoritmo desado, que puede devolverse disminuido dado que se lee el
! checkpoint. En modo escritura nsteps es el numero de pasos que ya se
! realizaron en el presente algoritmo y quedan salvados por el checkpoint.
integer,intent(inout)   :: step,nsteps
logical,intent(out),optional  ::  control

if(present(control)) control =.true.
current = current+1

! Me fijo si es hay un archivo de checkpoint para leer
if (chpmode) then

  ! Leo el archivo
  control=search_chp(step,nsteps)

#ifdef HAVE_MPI
!   ! Espero que todos los otros procesadores lean su propio archivo
! FIXME  call mpi_barrier()
#endif

  return
endif

! Escribo (o sobrescribo) el checkpoint
call write_chp(step,nsteps)
end subroutine readwrite_chp
      
function search_chp(step,nsteps,fname)
! Esta subrrutina es llamada cuando se esta por entrar al algoritmo en el cual
! se corto el calculo. Antes de entrar leo el paso del chepoint, tomo la
! configuracion del checkpoint, y disminuyo el numero de pasos acorde.
integer,intent(out)    :: step,nsteps
logical                :: search_chp
character(*),optional  :: fname      
integer                :: i,j,na
logical                :: ghosted

if(present(fname)) then
  open(chpunit, file=trim(adjustl(fname)),form='unformatted')
else
  open(chpunit, file=trim(adjustl(ioprefix))//".chp",form='unformatted')
endif
        
read(chpunit) i,step,nsteps,na

! Si hay un checkpoint guardado para leer en un algoritmo proximo pero no en el
! presente. Devuelvo false, como flag de control para poder evitar correr el
! presente siendo que mas adelante hay uno con checkpoint. Si hay un
! checkpoint guardado para leer en este algoritmo. Lo leo, y cambio el modo de
! checkpoint a escritura.

search_chp=(current==i)

if(search_chp) then
  
  call werr('Atoms in checkpoint are more or less tan the atoms in the system',na/=nlocal)

  nsteps=nsteps-step

  ! Read the box size
  read(chpunit) tbox(:,:)
  call box_setvars()

  ! Read random number generator state
  call read_chpseed(chpunit)
        
  ! Read piston information if is needed
  call read_chppiston(chpunit)
              
  read(chpunit) nframe,time,dm_steps,pnframe,ptime,frame
  do i =1,nlocal
    read(chpunit) (a(i)%o%pos(j),j=1,dm)
    read(chpunit) (a(i)%o%vel(j),j=1,dm)
    read(chpunit) (a(i)%o%acel(j),j=1,dm)
    read(chpunit) (a(i)%o%acel2(j),j=1,dm)
    read(chpunit) (a(i)%o%acel3(j),j=1,dm)
    read(chpunit) (a(i)%o%acel4(j),j=1,dm)
  enddo

  chpmode=.false.

  call posvel_changed()
   
  ! CHECK: I think this require full ghost updates
  if(ghost) then
    ghosted=fullghost
    fullghost=.true.
  endif
  call interact(.false.)
  if(ghost) then
    fullghost=ghosted
  endif
          

endif      

close(chpunit)

end function search_chp

subroutine write_chp(step,nsteps,fname)
integer,intent(in)    :: step,nsteps
character(*),optional :: fname      
integer               :: i,j

if(present(fname)) then
  open(chpunit, file=trim(adjustl(fname)),form='unformatted')
else
  open(chpunit, file=trim(adjustl(ioprefix))//".chp",form='unformatted')
endif

write(chpunit) current,step,nsteps,nlocal

! Read the box size
write(chpunit) tbox(:,:)

! Write random number generator state
call write_chpseed(chpunit)

! Read piston information if is needed
call write_chppiston(chpunit)

write(chpunit) nframe,time,dm_steps,pnframe,ptime,frame
do i =1,nlocal
  write(chpunit) (a(i)%o%pos(j),j=1,dm)
  write(chpunit) (a(i)%o%vel(j),j=1,dm)
  write(chpunit) (a(i)%o%acel(j),j=1,dm)
  write(chpunit) (a(i)%o%acel2(j),j=1,dm)
  write(chpunit) (a(i)%o%acel3(j),j=1,dm)
  write(chpunit) (a(i)%o%acel4(j),j=1,dm) 
enddo
close(chpunit)

end subroutine write_chp

subroutine read_chp(step,nsteps,fname)
! Esta subrrutina es llamada cuando se esta por entrar al algoritmo en el cual
! se corto el calculo. Antes de entrar leo el paso del chepoint, tomo la
! configuracion del checkpoint, y disminuyo el numero de pasos acorde.
integer,intent(out)    :: step,nsteps
character(*),optional  :: fname  
integer                :: i,j,na
logical                :: ghosted

if(present(fname)) then
  open(chpunit, file=trim(adjustl(fname)),form='unformatted')
else
  open(chpunit, file=trim(adjustl(ioprefix))//".chp",form='unformatted')
endif
        
read(chpunit) i,step,nsteps,na

call werr('Atoms in checkpoint are more than the atoms in the system',na>nlocal)
call wwan('Atoms in checkpoint are less than the atoms in the system',na<nlocal)

! Read the box size
read(chpunit) tbox(:,:)
call box_setvars()

! Read random number generator state
call read_chpseed(chpunit)
           
! Read piston information if is needed
call read_chppiston(chpunit)
                  
read(chpunit) nframe,time,dm_steps,pnframe,ptime,frame
do i =1,na
  read(chpunit) (a(i)%o%pos(j),j=1,dm)
  read(chpunit) (a(i)%o%vel(j),j=1,dm)
  read(chpunit) (a(i)%o%acel(j),j=1,dm)
  read(chpunit) (a(i)%o%acel2(j),j=1,dm)
  read(chpunit) (a(i)%o%acel3(j),j=1,dm)
  read(chpunit) (a(i)%o%acel4(j),j=1,dm)
enddo
close(chpunit)

call posvel_changed()


! CHECK: I think this require full ghost updates
if(ghost) then
  ghosted=fullghost
  fullghost=.true.
endif
call interact(.false.)
if(ghost) then
  fullghost=ghosted
endif
                    

end subroutine read_chp

end module gems_checkpoint


