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
!    ___________________   _____    _________
!   /  _____/\_   _____/  /     \  /   _____/
!  /   \  ___ |    __)_  /  \ /  \ \_____  \
!  \    \_\  \|        \/    Y    \/        \
!   \______  /_______  /\____|__  /_______  /
!          \/        \/         \/        \/ 
!       is an Extensible Molecular Simulator
!
program gems
 use gems_program_types
 use gems_clinterpreter
 use gems_random
 use gems_input_parsing
 use gems_errors
 use gems_programs
 use gems_hyperdynamics
 use gems_neighbor, only:polvar_neighbor
 use gems_fields
 use gems_elements
 use gems_output, only: chpmode
 use gems_mpi, only: mpi_pc
 use gems_variables, only:polvar_hard, polvars, polvar_link, polvar_expand, polvar_save, polvar_hard, polvar_readonly

#ifdef HAVE_MPI
 use gems_mpi, only: gems_mpi_init, mpi_tpc, ch_mpi_pc
 use mpi_f08
#endif

!$ use omp_lib 
 
 implicit none

 ! Command line argument
 integer                    :: narg,carg,larg !number of arguments, counter and length
 character(:),allocatable   :: arg            !argument

 integer      :: i

 ! Si se compila con soporte para mpi, inicializo los procesadores.
#ifdef HAVE_MPI
  call gems_mpi_init()
#endif

  ! Init wall time
  call system_clock(count_rate=sclock_rate)
  call system_clock(count_max=sclock_max)
  call system_clock(sclock_t1)


  !Check if any arguments are found
  narg=command_argument_count()

  !Check for optional arguments
  if(narg>0) then

    do carg=1,narg
      call get_command_argument(1,length=larg)
      allocate(character(larg) :: arg)
      call get_command_argument (1, value=arg)

      arg=adjustl(arg)

      select case(arg)
       case("--help","-h")
         write(*,*) "This is a test"
         exit
       ! case("-n")
       !  ! Quiero ser el procesador n!!
       case default
        stdin=arg
        stdinlen=larg
      end select

    end do

  endif
    
  ! Ask: Input from standar input or from a .gms file
  if (stdinlen>0) then
  ! if (stdin/='') then

    ! Set ioprefix as the input name
    if(stdin(stdinlen-3:stdinlen)==".gms") then
      ioprefix=stdin(1:stdinlen-4)
    else
      ioprefix=stdin(1:stdinlen)
    endif

    ! Inputfile
    inunit=250+mpi_pc*4
    open(inunit,file=trim(adjustl(stdin)))

#ifdef HAVE_MPI
    if (mpi_tpc>1) then
      ! Anexo el numero de procesador al ioprefix
      ioprefix = trim(ioprefix) //".mpi"// trim(ch_mpi_pc)
    endif
#endif

    ! File uits
    logunit=253+mpi_pc*4
    chpunit=254+mpi_pc*4

    ! Logfile
    logfile=trim(adjustl(ioprefix))//".log"
    open( logunit, file=trim(adjustl(logfile)))

    ! If stdin is a chp file, enable the checkpoint
    if (stdin(stdinlen-3:stdinlen)=='.chp') then

      ! Busco si hay un checkpoint para lectura.
      chpfile=trim(stdin)
      inquire(file=trim(adjustl(chpfile)),exist=chpmode)
      call werr('wrong checkpoint file',.not.chpmode)

      stdin=trim(stdin(1:stdinlen-4))//'.gms'

    else

      chpmode=.false.
      chpfile=trim(adjustl(ioprefix))//".chp"

    endif
 

  else

    ! Input from standar input
    
    ! Stndin is standar in
    stdin='("standar in")'
    stdinlen=3+4

    ! Log will be standar out
    logfile='("standar out")'

    ! ioprefix will be fixed
    ioprefix='gms'

    ! No checkpoint file
    chpfile='none'

  endif
  
  ! call get_environment_variable("GEMS_DOC",length=i1,status=stat)
  !
  ! if(stat>0) then
  !   !Not found, use default prefix 
  !   docdir=DOCDIR
  ! else
  !   allocate(character(i1) :: docdir)
  !   call get_environment_variable("GEMS_DOC",docdir)
  ! endif
  !
  ! call get_environment_variable("GEMS_PRM",prmdir,status=stat)


  call wlog('','gms file') ! VIM Syntax Highlight

  !Escribo el header del programa
  call wlog('',''                                          )
  call wlog('','  ___________________   _____    _________') 
  call wlog('',' /  _____/\_   _____/  /     \  /   _____/') 
  call wlog('','/   \  ___ |    __)_  /  \ /  \ \_____  \ ') 
  call wlog('','\    \_\  \|        \/    Y    \/        \') 
  call wlog('',' \______  /_______  /\____|__  /_______  /') 
  call wlog('','        \/        \/         \/         \/') 
  call wlog('','      is an Extensible Molecular Simulator')
  call wlog('','')
  call wlog('','Copyright (C) 2020 Sergio Alexis Paz                                    ')
  call wlog('','This program comes with ABSOLUTELY NO WARRANTY; for details type `warranty`. ')
  call wlog('','This is free software, and you are welcome to redistribute it                ')
  call wlog('','under certain conditions; type `license` for details.                        ')
  call wlog('','')
  call wlog('','Version:'//PACKAGE_VERSION                   )
  call wlog('','I/O Units:'                                  )
  call wlog('','  -input: '//trim(adjustl(stdin))            )
  call wlog('','  -log:   '//logfile                         )
  call wlog('','  -chp:   '//chpfile                         )
  call wlog('','    ---> trying to recover from chpfile ',chpmode)

#ifdef HAVE_MPI
  call wlog(''); write(logunit,'(a,i0,a,i0)') "MPI: Im thread ", mpi_pc, " of ", mpi_tpc
#else  
  call wlog(''); write(logunit,'(a)') "MPI: Not compiled for MPI"
#endif

  !$OMP PARALLEL
  !$OMP SINGLE
  !$ call wlog(''); write(logunit,'(a,i0,a)') 'OpenMP: Using ',omp_get_num_threads(),' threads'
  !$OMP END SINGLE
  !$OMP END PARALLEL

! INCIALIZACION DE OBJETOS Y VARIABLES

  ! Init default random seed
  call std_init()
  call init_ran()

  ! Init default elements
  call elements_init()

  ! Init default output precision
  prf='20.7'
  pri='0'

  ! Variable caracter con la dimension de compilacion
  write(cdm,'(i1)') dm

  ! Variables interpretadas (Metavariables)
  call polvars%init()

  ! Variables de modulos
  call polvar_neighbor()

  ! Guardo algunas variables utiles

  ! $jobname$
  ioprefix=trim(adjustl(ioprefix))
  call polvar_hard('jobname',ioprefix)

  ! $i$
  call polvar_hard('i',0)
  ! write(var_value%o(var_value%size),fmt='(i0)') 0
     
  ! $ci$ TESTME Puede que cambie el indice de las demas al meterlo aca
  call polvar_hard('ci','0')
  ! write(var_value%o(var_value%size),fmt='(i0)') 0

  call polvar_link('boxx',box(1))
  call polvar_link('boxy',box(2))
  call polvar_link('boxz',box(3))
  call polvar_link('time',time)
  call polvar_link('step',dm_steps)
  call polvar_link('dt',dt)
  
  call polvar_link('pi',pi);call polvar_readonly('pi');

#ifdef HAVE_MPI

    ! $pc$
    call polvar_link('pc',mpi_pc)

    ! $tpc$
    call polvar_link('tpc',mpi_tpc)
    ! write(arg,fmt='(i0)') mpi_tpc
    ! call var_name%append('tpc',trim(arg))

#endif

!ALLOW EXECUTION BLOQUES INSIDE PROGRAMS
exec => execute_command
              
!ALLOW VARIABLES AND LABELS INSIDE INPUT PARSING
! Variables that are set on execution will be hard and character. I only care about characters, since for parsing expansion that is
! enough.
var_save => polvar_save
var_set => polvar_hard
var_expand => polvar_expand

! Groups
! ------

! gindex
call gindex%init()

! System
call sys%init()
   
! Ghost
call ghost%init()
   
! Selection
call gsel%init()
  
! FIXME:
call gmeta%init()

! CLI memory groups
do i = 1,mgr
  call gr(i)%init()
enddo

! Initialize Integration, Interaction and OutputFiles vectors
call ngindex%init()
call of_vop%init()
call its%init()

! This is the command interpreter driver. In general execute_command subroutine is
! call to parse the command. However, there are a few commands that are directly
! here. The rason is that this commands call a subprogram. I call subprogram to
! that portion of code that can call internal to the command interpreter too.
! (e. g. when a dinamic is run, a block of code can be interpreted each a
! certain amounght of time.). In general this subprograms also have a
! checkpoint capability.  Esto es un problema porque dentro de  un bloque no se
! puede llamar a estos programas, pero es necesario para no crear lio con el
! bloque save. Entonces el bloque save es incompatible con el repeat.

  truelogunit=logunit
  printunit=>logunit
  
  opts => gems_iopts
  call opts%init(echo_lines=.true.,error_flag=1,in=inunit,out=logunit)

  do
    call read_line(eof)
    call flush(logunit)
    if (eof.or.term_signal) exit

    call readl(w1)
    call execute_command(trim(adjustl(w1)))

  enddo


!FIN
  call cpu_time(time1)
  if (time1<60) then
    call wstd(); write(logunit,'("cpu time: ",f10.3," s")') time1
  elseif(time1<3600) then
    call wstd(); write(logunit,'("cpu time: ",i0," m ",f10.3," s")') int(time1/60.0_dp),mod(time1,60.d0)
  else
    call wstd(); write(logunit,'("cpu time: ",i0," h ",i0," m")') int(time1/3600.0_dp),int(mod(time1,3600.d0)/60.0_dp)
  endif 

  call system_clock(sclock_t2)
  time1=(sclock_t2-sclock_t1)/real(sclock_rate)
  if (time1<60) then
    call wstd(); write(logunit,'("wall time: ",f10.3," s")') time1
  elseif(time1<3600) then
    call wstd(); write(logunit,'("wall time: ",i0," m ",f10.3," s")') int(time1/60.0_dp),mod(time1,60.d0)
  else
    call wstd(); write(logunit,'("wall time: ",i0," h ",i0," m")') int(time1/3600.0_dp),int(mod(time1,3600.d0)/60.0_dp)
  endif 
 
!  call wstd(); write(logunit,*) 'with ',rupdate,' neighbor list acualizations'
  call wstd(); write(logunit,'("vecinos actualizados: ",i0," veces")') nupd_vlist

  if(b_ckp) call system('rm -rf '//chpfile)

  close(chpunit)


#ifdef HAVE_MPI
  call mpi_finalize()
#endif

! Free default random seed
call free_ran()

end program gems


