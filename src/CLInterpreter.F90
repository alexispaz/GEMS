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

 
module gems_clinterpreter

! This module content the regular command that can be interpreted by GEMS. 
! The main subroutine is execute_command that interpret the main commands. In
! general, all the others subroutine interpret the childs command of this.

use gems_program_types
use gems_constants, only: dp,linewidth

use gems_input_parsing
use gems_errors

use gems_inq_properties 
use gems_set_properties

use gems_select_create
use gems_checkpoint
use gems_neighbour
use gems_integration
use gems_interaction
use gems_output

use gems_neb

use gems_quasi_newton

use gems_random
use gems_algebra
use gems_variables, only:polvar_expand

implicit none

integer                   :: ip(2)
real(dp)                  :: fv(dm),fv2(dm),fv3(dm)
logical                   :: b1,bv1(dm)!,b2,b3

! Variables auxiliares para lectura, parsing y demas. 
character(len=linewidth)  :: w1,w2,w3
integer                   :: i1,i2,i3,i4,i5,i6,i7,i8
real(dp)                  :: f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16 

! character(:)              :: help_file='DOCDIR/help.md'


contains

! EXECUTE COMMANDS

recursive subroutine execute_command(com)
use gems_programs
use gems_hyperdynamics
#ifdef HAVE_MPI
use gems_replicaexchange, only: parallel_tempering,parallel_tempering2
#endif
use gems_metadynamics, only:  metadynamics, wtmd2D_set,dm_cv_set,&
                               wall2D_set,wallauxCore_set,wallauxShell_set,&
                               wall1D_set,wtmetad_set,posicion1d_set
                              
character(*)  :: com
type(group)   :: gsel_aux

select case(com)
case('license')
  call wlog('','GEMS is free software: you can redistribute it and/or modify         ')
  call wlog('','it under the terms of the GNU General Public License as published by ')
  call wlog('',' the Free Software Foundation, either version 3 of the License, or   ')
  call wlog('','(at your option) any later version.                                  ')
  call wlog('','')
  call wlog('','You should have received a copy of the GNU General Public License ')
  call wlog('','along with GEMS.  If not, see <https://www.gnu.org/licenses/>.    ')
case('warranty')
  call wlog('','GEMS is distributed in the hope that it will be useful,        ')
  call wlog('','but WITHOUT ANY WARRANTY; without even the implied warranty of ')
  call wlog('','MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  ')
  call wlog('','GNU General Public License for more details.                   ')
  call wlog('','')
  call wlog('','You should have received a copy of the GNU General Public License ')
  call wlog('','along with GEMS.  If not, see <https://www.gnu.org/licenses/>.    ')
case('if')
  call readi(i1)
  if(i1==1) then
    call readl(w1)
    call execute_command(w1)
  endif

case('exit','stop')

  ! This allow to terminate the procedure that are interpreting commands.
  term_signal=.true.
   
case('cycle')

  ! This allow to cycle in a repeat block
  cycle_signal=.true.
   
case('<') ! Crea
  call create_commands(gnew)

case('+<') ! Crea y agrega a la seleccion previa

  ! WARNING FIXME Aquellas selecciones en donde gini=gsel, a medida que
  ! agreguen a gs (gout) incrementan el gini.... esto puede traer conflicto
  ! call create_commands(gnew,gsel)

  call gsel_aux%init('gsax')
  call create_commands(gnew,gsel_aux)
  call gsel%add(gsel_aux)
  call wstd(); write(logunit,*) gsel%nat, 'particles selected'
  call wwan('empty selection',gsel%nat==0)  
  call gsel_aux%dest()

case('><') ! Crea y selecciona lo creado

  call gsel_aux%init('gsax')
  call create_commands(gnew,gsel_aux)
  call group_allatom_del(gsel)
  call gsel%add(gsel_aux)
  call wstd(); write(logunit,*) gsel%nat, 'particles selected'
  call wwan('empty selection',gsel%nat==0) 
  call gsel_aux%dest()

! Seleccion en la creacion

case('.>') ! Seleccioname esto de la creacion

  call group_allatom_del(gsel)
  call select_commands(gnew,gsel)
  call wstd(); write(logunit,*) gsel%nat, 'particles selected'
  call wwan('empty selection',gsel%nat==0)

case('.+') ! Agrega esto a mi seleccion (Union de conjuntos)

  call select_commands(gnew,gsel)
  call wstd(); write(logunit,*) gsel%nat, 'particles selected'
  call wwan('empty selection',gsel%nat==0) 

! Seleccion en el sistema

case('>') ! Seleccioname esto del sistema

  call group_allatom_del(gsel)
  call select_commands(sys,gsel)
  call wstd(); write(logunit,*) gsel%nat, 'particles selected'
  call wwan('empty selection',gsel%nat==0)

case('+') ! Agrega esto a mi seleccion (Union de conjuntos)

  call select_commands(sys,gsel)
  call wstd(); write(logunit,*) gsel%nat, 'particles selected'
  call wwan('empty selection',gsel%nat==0) 

! Seleccion tanto en creacion como en el sistema

case('>>') ! Seleccioname esto de mi seleccion  (Interseccion de conjuntos)

  call gsel_aux%init('gsax')
  call select_commands(gsel,gsel_aux)
  call group_allatom_del(gsel)
  call gsel%add(gsel_aux) 
  call gsel_aux%dest()

  call wstd(); write(logunit,*) gsel%nat, 'particles selected'
  call wwan('empty selection',gsel%nat==0)

case('-')  ! Borrame esto de mi seleccion (resta de conjuntos)

  call unselect_commands(gsel)
  call wstd(); write(logunit,*) gsel%nat, 'particles selected'
  call wwan('empty selection',gsel%nat==0)

! sets globales
case('dimension')
  call readi(i1)
  write(cdm,'(i1)') dd
case('temp') ! Solo para usar unidades de kT
  call readf(temp)
  kt=dsqrt(kB_ui*temp)
case('prng')
  call prng_commands
! otros comandos
case('prueba')
  call write_state
case('constrain')
  call readl(w1)
  call constrain_commands(w1)
case('set')
  call set_commands
case('interact')
  call interacciones
case('getin') ! 
  call get_commands
case('sys')
  call sys_commands
case('group')
  call group_commands
case('mpi')
  call mpi_commands
case('log')
  call log_commands
case('outfile')
  call outfile_commands
case('print')
  i2=print_commands()
  do i1 =1, i2
    write(w1,'(a,i0,a)') '_ans[',i1,']'
    w2=polvar_expand(w1)
    call wprt(w2)
  enddo
case('time')
  call time_commands
case('walltime')
  call cpu_time(time1)
  if (time1<60) then
    call wstd(); write(logunit,'(f10.3,"s from start")') time1
  elseif(time1<3600) then
    call wstd(); write(logunit,'(i0,"m ",f10.3,"s from start")') int(time1/60.0_dp),mod(time1,60.d0)
  else                    
    call wstd(); write(logunit,'(i0,"h ",i0,"m from start")') int(time1/3600.0_dp),int(mod(time1,3600.d0)/60.0_dp)
  endif 
case('box')
  call box_commands
case('element')
  call element_commands
case('out')
  call out_commands
case('checkpoint')
  call checkpoint_commands
case('evolve')
  call evolve_commands
case('cv')
  call cv_commands
case('clean')
  call clean_commands
case('table')
  call table_commands
case('mtable')
  call mtable_commands



case('neb_clean')
  deallocate(pneb) 
  deallocate(fiximg) 
case('neb_nimg')
  call readi(nimg)   ! Imagenes
  nebdim = gsel%nat*dm
  allocate(pneb(nimg,nebdim)) 
  allocate(fiximg(nimg)) 
  nfiximg = 0
case('neb_img')
  call readi(i1)
  call neb_fiximg(gsel,i1)
case('neb_run')
  call readi(i1)   ! Iteracciones
  call readf(f1)   ! Constante resorte
  call neb(gsel,f1,i1,.true.)
case('lbfgs')
  b1=.true.
  if (nitems>item) call readb(b1)
  call lbfgs_minimizator(gsel,b1)
case('fix_lbfgs')
  call fix_lbfgs_minimizator(gsel,.true.)
! case('escape')
!  call readf(f1)
!  call escape_dinamic(gsel,f1,.true.) 
! case('hyperd_escape')
!  call readi(i1)  ! escapes
!  call readf(f1,ev_ui)  ! E
!  call readf(f2) ! factor para el poso
!  call readi(i2)  ! numero de medidas de eprom 
!  call readi(i3)  ! pasos en cada medida de eprom
!  call hyperd_escape(i1,f1,f2,i2,i3,.true.)  
!case('hyperd_eprom')
!  ! Ojo con este promedio que es sobre la funcion a la cual se aplica bias
!  call readi(i1)  ! numero de pasos
!  hd_eprom=hyperd_eprom(i1,.true.,sigma=f2)  
!  call wstd(); write(logunit,*) 'Energia promedio:'
!  call wstd(); write(logunit,*) '  eprom ', hd_eprom*ui_ev,'+-',f2*ui_ev
case('hybrid_wb','hybrid_ea')
  if(com=='hybrid_wb') then
    call readi(i4)  ! numero de grados de libertad
  endif
  call readi(i1)  ! numero de bloques
  call readi(i2)  ! pasos por bloque
  call readi(i3)  ! pasos de equilibracion
  call readf(f1)  ! eprime o w
  call readf(f4)  ! aprime o B
  if(com=='hybrid_wb') then
    call wb2ea(f1,f4,f2,f3,i4)
    f1=f2
    f4=f3
  endif
  call readf(f2)  ! Factor para decidir cuando una cosa es un plato
  call readf(f3)  ! Probabilidad para el maximo bias permitido
  call readf(f5)  ! Temperatura asociada con el termino de potencial a boostear

  call readwrite_chp(i3,i1,b1) 
  !if(b1) call lp_hyperdinamics(i1,i2,f1,f2,f4,f3,.true.)  
  !              nsteps,msteps,lsteps,hd_f,eprime,aprime,z,b_out
  if(b1) call hybrid_ea(i1,i2,i3,f2,f1,f4,f3,f5,.true.)  
  
case('dinamica_from_xyz')
  call readi(i1)
  call reada(w1)
  call dinamic_from_xyz(i1,w1,.true.)
case('dinamica')
  call readi(i1)
  i2=1
  call readwrite_chp(i2,i1,b1)
  if(chpmode)stop
  if(b1) call dinamic(i1,.true.,.true.) 
           
case('setwtmd1d')
  !Lectura de los parametros necesarios para Well tempered Metadynamics
  call readf(f1)   ! altura inicial de las gaussianas
  call readf(f2)   ! ancho de las gaussianas
  call readf(f3)   ! Parametro de la WTMD
  call readf(f4)   ! tau (frecuencia)
  call readi(i5)   ! cada cuanto grabo la dCM
  call readi(i6)   ! cantidad de particulas en el core
  !Pasar los parametros al modudo Metadynamics         
  i3=1
    call werr('No integration group',its%size<1)
    if(its%size>1) then
      call werr('There is more than one integration group, specify which to exchange',item==nitems)
      call readi(i3)
    endif  
  call wtmetad_set(f1,f2,f3,f4,i5,i6,i3)

case('setpos1d')
  !Lectura de los parametros necesarios para Well tempered Metadynamics
  call readf(f1)   ! altura inicial de las gaussianas
  call readf(f2)   ! ancho de las gaussianas
  call readf(f3)   ! Parametro de la WTMD
  call readf(f4)   ! tau (frecuencia)              /wt
  call readi(i5)   ! cada cuanto grabo la dCM
  !Pasar los parametros al modudo Metadynamics  
  i3=1
    call werr('No integration group',its%size<1)
    if(its%size>1) then
      call werr('There is more than one integration group, specify which to exchange',item==nitems)
      call readi(i3)
    endif  
  call posicion1d_set(f1,f2,f3,f4,i5,i3)

 case('setwall1d')
  call readf(f5)   ! potencial inicial
  call readf(f6)   ! potencial final
  call readi(i4)   ! potencial bines
  call readf(f9)   ! radio de la burbuja
  call wall1D_set(f5,f6,i4,f9)

case('wtdcm')
  call readi(i1)
  i2=1
  call readwrite_chp(i2,i1,b1)
  if(chpmode)stop
  !well-tempered Metadynamics 1D en dCM
  if(b1) call metadynamics(i1,.true.,.false.,.false.,.false.,.false.,.true.,.false.,.false.,.false.)


case('wtx')
  call readi(i1)
  i2=1
  call readwrite_chp(i2,i1,b1)
  if(chpmode)stop
  !well-tempered Metadynamics en x
  if(b1) call metadynamics(i1,.true.,.false.,.false.,.false.,.false.,.false.,.false.,.true.,.false.)

case('wtpos1d')
  call readi(i1)
  i2=1
  call readwrite_chp(i2,i1,b1)
  if(chpmode)stop
  !well-tempered Metadynamics 1D en dCM
  if(b1) call metadynamics(i1,.true.,.false.,.false.,.false.,.false.,.false.,.true.,.false.,.false.)

case('setwtmd2d')
  !Lectura de los parametros necesarios para Well tempered Metadynamics
  call readf(f1)   ! altura inicial de las gaussianas
  call readf(f2)   ! ancho de dcm
  call readf(f3)   ! ancho de rg
  call readf(f4)   ! Parametro de la WTMD
  call readf(f5)   ! tau (frecuencia)
  call readi(i6)   ! cada cuanto grabo la dCM
  call readf(f15)  ! tolerancia en las gaussianas
  call readi(i8)   ! cantidad de particulas en el core
  !Pasar los parametros al modudo Metadynamics 
  i3=1
    call werr('No integration group',its%size<1)
    if(its%size>1) then
      call werr('There is more than one integration group, specify which to exchange',item==nitems)
      call readi(i3)
    endif  
  call wtmd2D_set(f1,f2,f3,f4,f5,i6,f15,i8,i3)

case('setwall2d')
  call readf(f6)   ! potencial inicial dcm
  call readf(f7)   ! potencial final dcm
  call readi(i4)   ! potencial bines dcm
  call readf(f10)  ! potencial inicial rg
  call readf(f11)  ! potencial final rg
  call readi(i5)   ! potencial bines rg
  call readf(f14)  ! radio de la burbuja
  call wall2D_set(f6,f7,i4,f10,f11,i5,f14)


case('setwallauxcore')
  call readf(f6)   ! potencial inicial dcm
  call readf(f7)   ! potencial final dcm
  call wallauxCore_set(f6,f7)
case('setwallauxshell')
  call readf(f6)   ! potencial inicial dcm
  call readf(f7)   ! potencial final dcm
  call wallauxShell_set(f6,f7)

case('wtdcmrgco')
      
  call readi(i1)
  i2=1
  call readwrite_chp(i2,i1,b1)
  if(chpmode)stop
  !well-tempered Metadynamics 2D
  if(b1) call metadynamics(i1,.true.,.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false.)
!!!!!
case('wtxy')
      
  call readi(i1)
  i2=1
  call readwrite_chp(i2,i1,b1)
  if(chpmode)stop
  !well-tempered Metadynamics 2D
  if(b1) call metadynamics(i1,.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.true.)
!!!!!

case('wtdcmrgau')
      
  call readi(i1)
  i2=1
  call readwrite_chp(i2,i1,b1)
  if(chpmode)stop
  !well-tempered Metadynamics 2D
  if(b1) call metadynamics(i1,.true.,.false.,.false.,.false.,.true.,.false.,.false.,.false.,.false.)
!!!!!

case('wtdcmrgtotal')
      
  call readi(i1)
  i2=1
  call readwrite_chp(i2,i1,b1)
  if(chpmode)stop
  
  !well-tempered Metadynamics 2D en dCm y Rg total
  if(b1) call metadynamics(i1,.true.,.false.,.false.,.true.,.false.,.false.,.false.,.false.,.false.)
!!!!!


case('dmcv')
  call readi(i1)
  i2=1
  call readwrite_chp(i2,i1,b1)
  if(chpmode)stop
  !Lectura de los parametros necesarios para Parabola (Umb_Samp)
  call readi(i5)   ! cada cuanto grabo la dCM
  call readf(f6)   ! burbuja
  call readi(i6)   ! cantidad de particulas en el core
  !Pasar los parametros al modudo Metadynamics
  call DM_CV_set(i5,f6,i6)
  !Umbrella sampling
  if(b1) call metadynamics(i1,.true.,.false.,.true.,.false.,.false.,.false.,.false.,.false.,.false.)
!!!!!

case('partemp2')
#ifdef HAVE_MPI
  call readi(i1)
  call readi(i2)
  i3=1
  call werr('No integration group',its%size<1)
  if(its%size>1) then
    call werr('There is more than one integration group, specify which to exchange',item==nitems)
    call readi(i3)
  endif
  call parallel_tempering2(i1,i2,i3,.true.,.true.) 
#else
  call werr('Compile GEMS with MPI (./configure --with-mpi)')
#endif  
    
case('partemp')
#ifdef HAVE_MPI
  call readi(i1)
  call readi(i2)
  call readl(w1)
  select case (w1)
   case('dmcv')
    call readi(i5)   ! cada cuanto grabo la dCM
    call readf(f6)   ! burbuja
    call readi(i6)   ! cantidad de particulas en el core
    !Pasar los parametros al modudo Metadynamics
    call DM_CV_set(i5,f6,i6)
   !PTWTMD
   case('wtdcm')
    !Lectura de los parametros necesarios para Well tempered Metadynamics
    call readf(f1)   ! altura inicial de las gaussianas
    call readf(f2)   ! ancho de las gaussianas
    call readf(f3)   ! Parametro de la WTMD
    call readf(f4)   ! tau (frecuencia)
    call readi(i5)   ! cada cuanto grabo la dCM
    call readi(i6)   ! cantidad de particulas en el core
    !Pasar los parametros al modudo Metadynamics
    i3=1
    call werr('No integration group',its%size<1)
    if(its%size>1) then
      call werr('There is more than one integration group, specify which to exchange',item==nitems)
      call readi(i3)
    endif
    call wtmetad_set(f1,f2,f3,f4,i5,i6,i3)
    case default
    call werr('I do not understand the last command')  
  endselect
  
  call parallel_tempering(i1,i2,i3,.true.,.true.,w1) 
#else
  call werr('Compile GEMS with MPI (./configure --with-mpi)')
#endif  

case('help')

  call readl(w1)
  call help(w1)

case default

  ! Nothing knew... 
  call wwan('Unknown command')

end select 

endsubroutine execute_command

subroutine interacciones
use gems_elements, only: inq_z
use gems_pair_pot0
use gems_tb
use gems_tersoff
use gems_analitic_pot
use gems_forcefield
use gems_bias, only: bias_new, bias_cli, igb
   
class(intergroup), pointer     :: ig
type(group)    :: g1,g2
logical        :: under,feels
type(atom_dclist),pointer :: la
integer        :: i

call g1%init()
call g2%init()

call readi(i1)
call reada(w1)

under=.false.
feels=.false.

selectcase(w1)
case('under')
  under=.true.
  call g1%add(gr(i1))
case('with')
  call readi(i2)
  if (i2==i1) then
    under=.true.
    call g1%add(gr(i1))
  else
    call g1%add(gr(i1))
    call g2%add(gr(i2))
  endif
case('feels')
  feels = .true.
  call readi(i2)
  call g1%add(gr(i1))
  call g2%add(gr(i2))
case('set')

  ! FIXME: insted of reading the type, GEMS should take it from the interaction id.
  call reada(w1)
  
  selectcase(w1)
  case('halfsho_plane')
    call werr('Interaction not found',i1>igr_vop%size)
    call analitical_cli(igr_vop%o(i1)%o,w1)
  case default
    call werr('Unknown command')
  end select

  call g1%dest()
  call g2%dest()
   
  return
case default
  call werr('Bad keyword. Use `under`, `with` or `feels`')
end select

call werr('No se agrego ningun atomo',g1%nat==0)



call readl(w1)
selectcase(w1)
case('shofix') 
  
  call werr("Only self interaction, use `under` keyword.",.not.under) 
  call readf(f1,ev_ui)
  call shofix_set(g1,f1)

case('bias')

  call werr("Only self interaction, use `under` keyword.",.not.under) 
  call werr("Can't call this interaction twice.",associated(igb))
  call reada(w1)
  ig=>bias_new(w1,g1)
  call bias_cli(ig,w1)
            
case('ff')

  call reada(w1)
  call read_psf(w1)
  call reada(w1)
  call read_prm(w1)
  
case('cuw') 

  call readf(f1) !epsilon
  call readf(f2) !sigma
  call readf(f3) !start
  call readf(f4) !radious
  call readi(i1) !coordinate
  call readf(f5) !rcut
  
  if(under) then 
    call cuw_set(g1,f1,f2,f3,f4,i1,f5)
  else                                
    call cuw_set(g1,f1,f2,f3,f4,i1,f5,g2)
  endif

case('rwca') 

  call readf(f1) !epsilon
  call readf(f2) !sigma
  call readf(f3) !start
  
  if(under) then 
    call rwca_set(g1,f1,f2,f3)
  else
    call rwca_set(g1,f1,f2,f3,g2)
  endif
        
case('wca') 

  call readf(f1)
  call readf(f2)
  
  if(under) then 
    call wca_set(g1,f1,f2)
  else
    call wca_set(g1,f1,f2,g2)
  endif
                 
case('sm1') 

  call readf(f1)
  call readf(f2)
  call readf(f3,ev_ui)
  
  if(under) then 
    call sm1_set(g1,f1,f2,f3)
  else                      
    call sm1_set(g1,f1,f2,f3,g2)
  endif
                  
case('slj') 

  call readf(f1)
  call readf(f2)
  call readf(f3)
  
  if(under) then 
    call slj_set(g1,f1,f2,f3)
  else
    call slj_set(g1,f1,f2,f3,g2)
  endif
                   
case('lj') 

  call readf(f1)
  call readf(f2)
  call readf(f3)
  
  if(under) then 
    call lj_set(g1,f1,f2,f3)
  else
    call lj_set(g1,f1,f2,f3,g2)
  endif
                          
case('plj') 

  call readf(f1)
  call readf(f2)
  
  if(under) then 
    call plj_set(g1,f1,f2)
  else
    call plj_set(g1,f1,f2,g2)
  endif

case('sho','shocm') 
    
  if(under) then 
    ig=>sho_new(w1,g1)
  else
    ig=>sho_new(w1,g1,g2)
  endif
      
  call sho_cli(ig)
                    
case('halfsho_plane','sho2d','oscar2d','leiva1d','pozoa1d','voter2d') 
  
  call werr("Only self interaction, use `under` keyword.",.not.under) 
  ig=>analitical_new(w1,g1)
  call analitical_cli(ig,w1)

case('sho_plane') 
  
  call werr("Only self interaction, use `under` keyword.",.not.under) 

  call readf(f1,ev_ui)
  call readf(f2)
  call readi(i1)

  call werr("Dimension out of bound.",i1>dm)
  call sho_plane_init(g1,f1,f2,i1)

case('sho_line') 
  
  call werr("Only self interaction, use `under` keyword.",.not.under) 

  call readf(f1,ev_ui)
  call readf(fv) ! axis_p
  call readf(fv2) ! axis_v
  
  call sho_line_set(g1,f1,fv,fv2)

case('lucas1d')

  ! Solo autointeraccion 
  call werr("Only self interaction, use `under` keyword.",.not.under)
  call readf(f1,kjm_ui)
  call readf(f2,kjm_ui)
  call readf(f3,kjm_ui)
  call readf(f4)
  call readf(f5)
  call readf(f6)
  call readf(f7)
  call readf(f8)
  call readf(f9)
  call lucas1d_set(g1,f1,f2,f3,f4,f5,f6,f7,f8,f9)

case('tb')

  call reada(w1)
  if(under) then 
    call tb_set(w1,g1)
  else
    call tb_set(w1,g1,g2)
  endif

case('tersoff')
  ! TODO ERR: SOLO PUEDE LLAMARSE UNA VEZ

  ! Solo autointeraccion (esto podria cambiar)
  call werr("Only self interaction, use `under` keyword.",.not.under)

  call tsf_set(g1)

case default
  call wwan('I do not understand the last command')
endselect

call g1%dest()
call g2%dest()


if(feels) then
  igr_vop%o(igr_vop%size)%o%half=.false.
  igr_vop%o(igr_vop%size)%o%newton=.false.
endif


! Actualizar lista de vecinos
! Si bien podría actualizar solamente la de la interaccion declarada
! necesito hacer para los atomos involucrados pos_old=pos y al hacer
! esto pyedo estar desactivando la actualizacion de alguna otra interaccion.
! Por ejemplo, si:
!   interact 1 ....
!   out state
!   dinamica 4
!   interact 2 ....
! Entonces si el comando `interact 2 ...` ejecuta pos_old=pos necesariamente
! hay que actualizar la lista de `interact 1 ...`.

 
! First time selection of neighboor list. Unless some particular list
! is requested in the set procedure of the force field, all the neighboor list
! are sablished for first time here.
do i=1,igr_vop%size
  ig => igr_vop%o(i)%o
  if(.not.associated(ig%lista)) call ig%setcells()
enddo


end subroutine interacciones
 

subroutine help(w)
  character(*),intent(in)   :: w
  ! character(:)              :: line
  integer                   :: i,j
  integer                   :: stat 
  character(50)             :: buffer  
  integer :: nch
     
  j = find_io(30)
  ! open(j,file=help_file)
  !
  !   line = ''
  !   do
  !
  !     read (j, "(a)", advance='no', iostat=stat, size=nch) buffer
  !     if (stat > 0) stop 'error reading file.'
  !     line = line // buffer(:nch)
  !
  !     ! end of record or file
  !     if (stat < 0) return
  !
  !   end do
  !
  !   print *, line
  ! close(j)

endsubroutine help

subroutine sys_commands

call readl(w1)
selectcase(w1)
case('add')    !FIXME Makeme remember the previously group selecction...

  call werr('No se agrego ningun atomo',gsel%nat==0)
             
  ! Por si se les quiere cambiar el tipo
  ! FIXME: Esto antes estaba debajo para evitar cambiar el tipo a la
  ! seleccion y solo a lo que se agrega al systema
  if(item<nitems) then
    call readelement(i1)
    call set_element(gsel,i1)
  endif
          
  ! ! Agrego los atomos
  ! call system_group_add(gsel)

  ! Agrego los atomos XGHOST
  call atoms_group_add(gsel)
          
  call inq_mass(sys) 
  call wstd(); write(logunit,*) natoms, 'particles in the system '

case default  
  call wwan('I do not understand the last command')  
endselect      
endsubroutine sys_commands

subroutine evolve_commands
use gems_integration,only: polvar_integrate
use gems_variables,only: polvar_link
use gems_integration,only: integrate_cli, integrate
type(integrate),pointer   :: it
integer                   :: i

! Get label name
call readl(w1)

! Find previous label
it=>null()
if(w1(1:1)==':') then
  it=>polvar_integrate(w1)
  if (associated(it)) call readl(w1)
endif

! If label is not found create new one
if (.not.associated(it)) then

  ! Create new interaction
  call its%append()
  it=>its%o(its%size)

  ! Create new label
  if(w1(1:1)==':') then
    ! Using label defined by user
    call polvar_link(trim(w1),it)
    call readl(w1)
  else
    ! Using label defined by default
    write(w2,'(a,i0)') ':',its%size
    call polvar_link(trim(w2),it)
    call wstd('The interaction label is '//trim(w2))
  endif

  ! Initialize the integrate
  call it%init_ext(gsel)

endif

call integrate_cli(it,trim(adjustl(w1)))

endsubroutine evolve_commands

subroutine cv_commands
use gems_cvs

! ! Interacciones
! ! Esto es beta, pero la idea es ir agregando atomos al grupo de interaccion
! case('igr')
!   call readi(i1)
!   call igr(i1)%adda(gsel)

! CVS
call readl(w2)
selectcase(w2)
case('cm')
  ncvs=ncvs+1

  call werr('Too many cvs',ncvs>10)
      
  call readl(w1)
  select case(w1)
  case('x') 
    call cvs(ncvs)%init(1,gsel)
    allocate(cvs(ncvs)%ir(1))
    cvs(ncvs)%ir(1)=1
    cvs(ncvs)%eval => cv_eval_cm
    cvs(ncvs)%jaco => cv_jaco_cm
  case('y') 
    call cvs(ncvs)%init(1,gsel)
    allocate(cvs(ncvs)%ir(1))
    cvs(ncvs)%ir(1)=2
    cvs(ncvs)%eval => cv_eval_cm
    cvs(ncvs)%jaco => cv_jaco_cm
  case('z') 
    call cvs(ncvs)%init(1,gsel)
    allocate(cvs(ncvs)%ir(1))
    cvs(ncvs)%ir(1)=3
    cvs(ncvs)%eval => cv_eval_cm
    cvs(ncvs)%jaco => cv_jaco_cm
  case('xyz','zxy','yzx','xzy','yxz','zyx') 
    call cvs(ncvs)%init(3,gsel)
    cvs(ncvs)%eval => cv_eval_cmpos
    cvs(ncvs)%jaco => cv_jaco_cmpos
  case default  
    call wwan('I do not understand the CV selected')  
  endselect  
         
  ! Save 1 real parameter
  call readf(f1)
  allocate(cvs(ncvs)%pr(1))
  cvs(ncvs)%pr(1)=f1*eV_ui
  
  call cvs(ncvs)%jaco()
  call cvs(ncvs)%eval()
  cvs(ncvs)%z(:) = cvs(ncvs)%t(:) 

  cvs(ncvs)%calc => cv_calc_sho
                      
case default  
  call wwan('I do not understand the CV selected')  
endselect  

endsubroutine cv_commands

subroutine clean_commands
call readl(w1)
selectcase(w1)
case('creation')
  call group_allatom_del(gsel)
  call gsel%add(sys)
  call group_allatom_del(gnew)
  call wstd(); write(logunit,*) 'sys is selected'
case default  
  call wwan('I do not understand the last command')  
endselect      
endsubroutine clean_commands
 
subroutine group_commands
call readi(i1)
call readl(w1)
selectcase(w1)
case('add')
  call gr(i1)%add(gsel) 
  call wstd(); write(logunit,*) gr(i1)%nat, 'particles in group', i1
case('clean')
  call group_allatom_del(gr(i1)) 
  call wstd(); write(logunit,*) gr(i1)%nat, 'particles in group', i1 
case default  
  call wwan('I do not understand the last command')  
endselect      
endsubroutine group_commands

subroutine create_commands(gn,gs)
! Agrega lo nuevo a gn y agrega en seleccion a gs si este esta presente
! Hay que tener en cuenta que gsel puede ser usado, asique no es cuestion 
! de borrarlo asi como asi, lo mismos sys
type(group),intent(inout)            :: gn
type(group),intent(inout),optional   :: gs
integer     :: i

call readl(w1)
selectcase(w1)
case('reply')
  call readf(fv)
  call readi(i1)
  call create_reply(gsel,fv,i1,gn,gs)
case('atom')
  call readf(fv)
  call create_atom(fv,gn,gs)
case('fillpbc')     ! Llena con n atomos agregados sequencialmente al pbc
  call wlog(''); write(logunit,*) "box: ",box
  call readi(i1)
  call readf(f1)
  if (item < nitems) then
    call readf(fv)
    call readf(fv2)
    call create_fill(f1,i1,gn,fv,fv2,gs,opt_pbc=.true.)
  else
    call create_fill(f1,i1,gn,gout=gs,opt_pbc=.true.)
  endif 
case('fill')     ! Llena con n atomos separados por un cierto radio
  call wlog(''); write(logunit,*) "box: ",box
  call readi(i1)
  call readf(f1)
  if (item < nitems) then
    call readf(fv)
    call readf(fv2)
    call create_fill(f1,i1,gn,fv,fv2,gs)
  else
    call create_fill(f1,i1,gn,gout=gs)
  endif
case('fillh')    ! Llena con hidrogenos los atomos de carbono.
  call create_fillh(gsel,gn,gs) 
case('read')

  ! Allow to especify the extension?
  ! call reada(ext)
  ! select case(ext)
  ! case ('xyz')
  ! case default
  !   call reread(-1)
  ! end select

  call reada(w1)
  i=len_trim(adjustl(w1))

  i1=1
  if (item < nitems) call readi(i1) ! The frame
   
  call create_file(w1,w1(i-2:i+1),gn,gs,i1)

case default  
call wwan('I do not understand the last command')  
endselect      
call wstd(); write(logunit,*) gnew%nat, 'particles in creation'

endsubroutine create_commands

  subroutine constrain_commands(w)
    integer                    :: i,j,k
    character(*)               :: w
    character(1)               :: c
    type (atom_dclist),pointer :: la
    logical                    :: lconst=.false.
    real(dp)                   :: aux

 
    selectcase(w)
    case('axis') 
      lconst=.true.
    case('plane') 
      lconst=.false.
    endselect

    call readf(fv)
    aux=dot_product(fv,fv) 
    fv=fv/sqrt(aux)

    la => gsel%alist
    do i = 1,gsel%nat
      la=>la%next
      if(.not.allocated(la%o%pconst)) allocate(la%o%pconst(dm))
      la%o%pconst=la%o%pos
      la%o%vconst=fv
      la%o%bconst=.true.
      la%o%lconst=lconst
    enddo 
   
    ! Lo de abajo creo es obsoleto....

    return

    selectcase(w)
    ! case('fix') 
    !   call readl(w2)
    !   selectcase(c) 
    !   case('x')
    !     k=0
    !   case('y')
    !     k=1
    !   case('z')
    !     k=2
    !   endselect 
    !
    !   call gr_fix%add(gsel)
    !   la => gsel%alist%next
    !   do i = 1,gsel%nat
    !     j = la%o%idv
    !     fix(j+k) = .true.
    !     fix_pos(j+k) = la%o%pos(k+1)
    !     la=>la%next
    !   enddo 
    !
    case('join')
      call readl(w2)
      selectcase(c) 
      case('x')
        k=0
      case('y')
        k=1
      case('z')
        k=2
      endselect 
    
      njoin=0
      la => gsel%alist%next
      do i = 1,gsel%nat
        j = la%o%idv
        join(j+k) = .true.
        njoin=njoin+1
        la=>la%next
      enddo 
    endselect

  endsubroutine constrain_commands
 
  subroutine unselect_commands(g)
    integer                    :: i
    type(group)                :: g
    type(group)                :: aux
    type (atom_dclist),pointer :: la
  
    call aux%init()
    call select_commands(g,aux)

    la => aux%alist%next
    do i = 1,aux%nat
      call group_atom_del(la%o, g)
      la=>la%next
    enddo

    call aux%dest()  ! Sin esto, habia un segmentation full

  endsubroutine unselect_commands
   
subroutine select_commands(gini,gout)
type(group),intent(inout)          :: gout
type(group),intent(inout),optional :: gini
integer                            :: i

call readl(w1)

! Seleccion donde gini entra por la ventana

selectcase(w1)
case('group')
  call readi(i1)
  call gout%add(gr(i1))
  return
case('creation')
  call gout%add(gnew)
  return
case('sys')
  call gout%add(sys)
  return
endselect

! Seleccion donde gini entra por la puerta

call werr('No existe seleccion previa',.not.present(gini))
call werr('No existe seleccion previa',gini%nat<1)

selectcase(w1)
case('cylinder')
  call readf(f1)
  call readf(fv2)
  if (item < nitems) then
    call readf(fv)
  else
    call inq_cg(gini)
    fv=gini%cg_pos
  endif
  call select_cylinder(f1,fv2,fv,gini,gout) 
case('sphere')
  call readf(f1)
  if (item < nitems) then
    call readf(fv)
  else
    call inq_cg(gini)
    fv=gini%cg_pos
  endif
  call select_sphere(f1,fv,gini,gout)
case('bo_range')
  
  call readi(i1) ! The second group related with bo calculation
  call readf(f1) ! Lower bound of the range
  call readf(f2) ! Upper bound of the range

  ! Default values for cut radious of bond order function
  f3=4.0_dp  ! minor cut radious
  f4=5.0_dp  ! mayor cut radious
  if (item < nitems) call readf(f3)  ! User value for minor cut radious
  if (item < nitems) call readf(f4)  ! User value for mayor cut radious

  ! Run selecction according to the second group selected
  if (i1==-1) then
    call select_bond_order_range(f1,f2,f3,f4,gini,gini,gout) 
  elseif (i1==0) then
    call select_bond_order_range(f1,f2,f3,f4,gini,sys,gout) 
  else
    call select_bond_order_range(f1,f2,f3,f4,gini,gr(i1),gout) 
  endif

case('bo_min')
   
  call readi(i1) ! The second group related with bo calculation

  ! Default values for cut radious of bond order function
  f3=4.0_dp  ! minor cut radious
  f4=5.0_dp  ! mayor cut radious
  if (item < nitems) call readf(f3)  ! User value for minor cut radious
  if (item < nitems) call readf(f4)  ! User value for mayor cut radious

  ! Run selecction according to the second group selected
  if (i1==-1) then
    call select_lower_bond_order(f3,f4,gini,gini,gout) 
  elseif (i1==0) then
    call select_lower_bond_order(f3,f4,gini,sys,gout) 
  else               
    call select_lower_bond_order(f3,f4,gini,gr(i1),gout) 
  endif

 case('bo_max')
   
  call readi(i1) ! The second group related with bo calculation

  ! Default values for cut radious of bond order function
  f3=4.0_dp  ! minor cut radious
  f4=5.0_dp  ! mayor cut radious
  if (item < nitems) call readf(f3)  ! User value for minor cut radious
  if (item < nitems) call readf(f4)  ! User value for mayor cut radious

  ! Run selecction according to the second group selected
  if (i1==-1) then
    call select_upper_bond_order(f3,f4,gini,gini,gout) 
  elseif (i1==0) then
    call select_upper_bond_order(f3,f4,gini,sys,gout) 
  else               
    call select_upper_bond_order(f3,f4,gini,gr(i1),gout) 
  endif
  
case('zrange')
  call readf(f1)
  call readf(f2)
  call select_band(3,f1,f2,gini,gout) 
case('yrange')
  call readf(f1)
  call readf(f2)
  call select_band(2,f1,f2,gini,gout) 
case('xrange')
  call readf(f1)
  call readf(f2)
  call select_band(1,f1,f2,gini,gout) 
case('above')
  call readf(fv)
  call readf(fv2)
  call readf(fv3)
  call select_above(fv,fv2,fv3,gini,gout) 
case('below')
  call readf(fv)
  call readf(fv2)
  call readf(fv3)
  call select_below(fv,fv2,fv3,gini,gout) 
case('rectangle')
  call readf(fv)
  call readf(fv2)
  call select_rectangle(fv,fv2,gini,gout) 
case('near')
  call readf(fv)
  call select_near(fv,gini,gout) 
case('far')
  call readf(fv)
  call select_far(fv,gini,gout) 
case('element')
  call readelement(i1)
  call select_element(i1,gini,gout)
case('atom')
  call readi(i1)
  i2=i1
  if (nitems>item) call readi(i2)
  do i = i1,i2
    call select_atom(i,gini,gout)
  enddo

!   case('from')
!       call readi(i1)
!       call readl(w2)
!       call readi(i2)
!       call readi(i3)
!       gini%nat = gini%nat + i3-i2
!       do j = i2,i3
!         call group_atom_add(ss(i1)%a(j),gini)
!       enddo 
case('random')
  call readi(i1)
  call select_random(i1,gini,gout)
case default  
call wwan('I do not understand the last command')  
endselect

endsubroutine select_commands

subroutine out_commands
use gems_output, only:prf,pri

call readl(w1)
  
selectcase(w1)
  case('posxyz')
    write(w1,*) "positions.xyz"
    if (nitems>item) call reada(w1)
    call write_screenshot(w1,gsel)   

  case('state')
    call interact(.true.)
    call write_out_force(1)
    
  case('float')
    call readl(prf)

  case('integer')
    call readl(pri)

end select
       
end subroutine out_commands

subroutine checkpoint_commands
call readl(w1)
selectcase(w1)
case('write')
  call reada(w1)
  call interact(.true.)
  call write_chp(1,1,w1)
case('read')
  call reada(w1)
  call interact(.true.)
  call read_chp(i1,i2,w1)
case('on')
  call werr('Not checkpoint aviable when standar input is used',stdin=='standar in') 
  b_ckp=.true.
case('off')
  b_ckp=.false.
case('each')
  call readi(chpeach)
case default  
  call wwan('I do not understand the last command')  
endselect
endsubroutine checkpoint_commands

  subroutine set_commands
    use gems_forcefield
    use gems_elements, only: inq_z

    call readl(w1)
    select case(w1)
             
    case('mass')
      call readf(f1)
      call set_mass(gsel,f1)
              
    case('cm_vel')
      call readf(fv)
      call set_cm_vel(gsel,fv)

    case('pbc') ! Establezco condiciones periodicas para ese grupo
      call readb(bv1)
      call set_pbc(gsel,bv1)
    case('tempgdist') ! distribución gaussiana 
      call readf(f1)
      call set_gdist(gsel,f1)
     ! call write_velocities_files(ss,'gdist') may be for reproducibilidad
    case('cm_pos') 
      call readf(fv)
      call set_cm_pos(gsel,fv)
    case('sigma') 
      call readf(f1)
      call set_sigma(gsel,f1)
    case('cg_pos') 
      call readf(fv)
      call set_cg_pos(gsel,fv)
    case('minpos') 
      call readl(w2)
      selectcase(w2)
      case('x')
        i1=1
      case('y')
        i1=2
      case('z')
        i1=3
      endselect 
      call readf(f1) 
      call set_minpos(gsel,f1,i1) 
    case("maxpos") 
      call readl(w2)
      selectcase(w2)
      case('x')
        i1=1
      case('y')
        i1=2
      case('z')
        i1=3
      endselect 
      call readf(f1) 
      call set_maxpos(gsel,f1,i1)  
    case("posdist") ! read position from file
      call reada(w1)
      i1=1
      if (item < nitems) call readi(i1) ! The frame
      call set_pos_from_file(gsel,trim(adjustl(w1)),i1) 
    case("veldist") ! read velocity from file
      call reada(w1)
      i1=1
      if (item < nitems) call readi(i1) ! The frame
      call set_vel_from_file(gsel,trim(adjustl(w1)),i1)        
    case("element")
      call readl(w1)
      call set_z(gsel,inq_z(w1))
    case("min_interdist")
      call readf(f1)
      call minterdist(gsel,f1)
    case('align')
      call align(gsel)
    case('align_mass')
      call align_mass(gsel)
    case('alignxy_mass')
      call alignxy_mass(gsel)
    case('align_by')
      call readi(i1)
      call bialign(gsel,gr(i1))
    case('align_massxy_by')
      call readi(i1)
      call bialignxy_mass(gsel,gr(i1))
    case('align_tutor')
      call readf(fv)
      call readf(fv2)
      call align_tutor(gsel,fv,fv2)
    case('rotateaxis')
      call readf(f1)
      call readf(fv)
      call readf(fv2)
      f1=f1*pi/180.0_dp
      call rotate_rodrigues(gsel,f1,fv,fv2)
   case('rotate')
     call readl(w2)
     selectcase(w2)
     case('x')
       ip=[3,2]
     case('y')
       ip=[1,3]
     case('z')
       ip=[2,1]
     endselect 
     call readf(f1)
     if (item == nitems-dm) then
       call readf(fv)
       call rotate_axis_ref(gsel,f1,ip,fv)
     else
       call givens_rotation(gsel,f1,ip)
     endif
    case('expand')
      fv2=0.0_dp   ! Default en el origen?
      call readf(fv)
      call readf(fv2)
      call expand(gsel,fv,fv2)
    case('move')
      call readf(fv)
      call move(gsel,fv)   
    case('do_pbc')
      call do_pbc(gsel)   
     case('sp')
      call readi(i1)
      call set_sp(gsel,i1)   
 
    case default
      call wwan('I do not understand the last command')
    endselect 
  endsubroutine set_commands

subroutine get_commands
use gems_random
use gems_variables,only:polvar_hard, polvar_expand
integer                :: i,j
character(linewidth)   :: var, vari

! Reding the variable name
call readl(var) 

! Las variables que empiezan con underscore son solo internas
! no para que definan o accedan los usuarios
call werr('Variable name should not start with underscore',var(1:1)=='_')

i2=print_commands()

if (i2>1) then
  do i =1, i2

    write(vari,'(a,i0,a)') trim(var)//'[',i,']'
    write(w1,'(a,i0,a)') '_ans[',i,']'

    w2=polvar_expand(w1)
    call polvar_hard(trim(vari),trim(w2))
    call wstd(trim(vari)//' = '//trim(w2))
   
  enddo 
  return
endif

w2=polvar_expand('_ans[1]')
call polvar_hard(trim(var),trim(w2))
call wstd(trim(var)//' = '//trim(w2))

endsubroutine get_commands

  subroutine mpi_commands
    use gems_random
    use gems_mpi, only: mpi_pc, mpi_tpc, ch_mpi_pc
    call readl(w1)
#ifdef HAVE_MPI
#else
    select case(w1)
    case('like') ! Para establecer la semilla distinta en cada procesador
      call readi(mpi_pc)
      call readi(mpi_tpc)
      call wstd(); write(logunit,*) 'Taken like MPI_PC=',mpi_pc,' MPI_TPC=',mpi_tpc 
      return
    endselect  
#endif

    select case(w1)
    case ('only')   ! Para correr una sintaxis en un procesador determinado
      call readi(i1)  
      call readl(w1)  
      if(mpi_pc==i1) call execute_command(w1)
    case default  
      call wwan('I do not understand the last command')  
    endselect    

  end subroutine  mpi_commands

  subroutine log_commands
  
    if(logunit/=truelogunit) close(logunit) 

    call readl(w1)
    select case(w1)
    case('>>')
      call reada(w1)
      logunit=find_io(30)
      open(logunit,file=trim(w1), position='append') 
    case('>')
      call reada(w1)
      logunit=find_io(30)
      open(logunit,file=trim(w1)) 
    case('std')
      logunit=6
    case('-')
      logunit=truelogunit
    endselect   
    call opts%init(out=logunit)

  end subroutine log_commands

function print_commands() result(ans)
use gems_random
use gems_tables, only: etable
use gems_variables,only:polvar_hard
integer                :: ans
type(atom_dclist),pointer :: la
type(etable) :: t
integer             :: i

!if (gsel%nat==int(gsel%mass)) call wwan('Atomo/s sin elemento definidos')

call pos_changed()

! Reading print order
call readl(w1)

ans=1

select case(w1)
case('log')
  if(.not.associated(printunit,target=logunit)) then
     close(printunit) 
     deallocate(printunit)
  endif
  printunit=>logunit
case('file')
  if(.not.associated(printunit,target=logunit)) then
     close(printunit) 
     deallocate(printunit)
  endif
  call reada(w1)
  allocate(printunit)
  printunit=find_io(30)
  open(printunit,file=trim(w1)) 
case('std')
  if(.not.associated(printunit,target=logunit)) then
     close(printunit) 
     deallocate(printunit)
  endif
  allocate(printunit)
  printunit=6

! Estos case es para cuando el resultado se puede guardar en ans

case('dm_steps')
  call polvar_hard('_ans[1]',dm_steps)
case('selnat') 
  call polvar_hard('_ans[1]',gsel%nat)
case('temp') 
  call inq_temperature(gsel)
  call polvar_hard('_ans[1]',gsel%temp)
case('cm_vel')
  call inq_cm_vel(gsel)
  call polvar_hard('_ans[1]',gsel%cm_vel(1))
  call polvar_hard('_ans[2]',gsel%cm_vel(2))
  call polvar_hard('_ans[3]',gsel%cm_vel(3))
  ans=3
case('cm_pos')
  call group_inq_cmpos(gsel)
  call polvar_hard('_ans[1]',gsel%cm_pos(1))
  call polvar_hard('_ans[2]',gsel%cm_pos(2))
  call polvar_hard('_ans[3]',gsel%cm_pos(3))
  ans=3
case('cm_diff2','cm_diff')
  call readi(i1) ! The second group
  ! Run selecction according to the second group selected
  call group_inq_cmpos(gsel)
  if (i1==-1) then
    fv = gsel%cm_pos-gsel%cm_pos
  else if (i1==0) then
    call pos_changed()
    call group_inq_cmpos(sys)
    fv = gsel%cm_pos-sys%cm_pos
  else
    call pos_changed()
    call group_inq_cmpos(gr(i1))
    fv = gsel%cm_pos-gr(i1)%cm_pos
  endif
  if(w1=='cm_diff2') then
    call polvar_hard('_ans[1]',dot_product(fv,fv))
  else
    call polvar_hard('_ans[1]',sqrt(dot_product(fv,fv)))
  endif
case('rg_pos')
  call group_inq_rg(gsel)
  call polvar_hard('_ans[1]',gsel%rg_pos)
case('minpos') 
  call inq_boundingbox(gsel)
  call polvar_hard('_ans[1]',gsel%minpos(1))
  call polvar_hard('_ans[2]',gsel%minpos(2))
  call polvar_hard('_ans[3]',gsel%minpos(3))
  ans=3
case('maxpos') 
  call inq_boundingbox(gsel) 
  call polvar_hard('_ans[1]',gsel%maxpos(1))
  call polvar_hard('_ans[2]',gsel%maxpos(2))
  ans=3
  call polvar_hard('_ans[3]',gsel%maxpos(3))
case('norm') 
  call readi(i1)
  call readi(i2)
  call readi(i3)
  ! BUG: Por alguna extraña razon no es lo mismo si se intercambia i2 con i1
  fv=cross_product(vdistance(a(i1)%o,a(i2)%o,ghost),vdistance(a(i3)%o,a(i2)%o,mic))
  f1=dot_product(fv,fv)
  fv=fv/sqrt(f1)
  call polvar_hard('_ans[1]',fv(1))
  call polvar_hard('_ans[2]',fv(2))
  call polvar_hard('_ans[3]',fv(3))
  ans=3
case('axis') 
  call readi(i1)
  call readi(i2)
  fv=vdistance(a(i1)%o,a(i2)%o, mic)
  f1=dot_product(fv,fv)
  fv=fv/sqrt(f1)
  call polvar_hard('_ans[1]',fv(1))
  call polvar_hard('_ans[2]',fv(2))
  call polvar_hard('_ans[3]',fv(3))
  ans=3
case('ptriaxial') 
  write(ans,fmt='(e15.8)') inq_triaxial_param(gsel)
case('mayordist') 
  ! Obtiene la distancia mas grande entre los atomos del grupo seleccionado y
  ! otro grupo
  if (item < nitems) then
    ! Si hay mas items, esta especificando otro grup

    call readi(i1) ! The second group

    ! Run selecction according to the second group selected
    if (i1==-1) then
      fv(:)=inq_mayordistance(gsel,gsel)
    else if (i1==0) then
      fv(:)=inq_mayordistance(gsel,sys)
    else
      fv(:)=inq_mayordistance(gsel,gr(i1))
    endif
    call polvar_hard('_ans[1]',fv(1))
    call polvar_hard('_ans[2]',fv(2))
    call polvar_hard('_ans[3]',fv(3))
    ans=3

  end if
case('border') 
  ! Obtiene el bond order value para el atomo seleccionado. Sale un error si
  ! se selecciono mas de 1 atomo
  call werr('Only one atom must be in the selection',gsel%nat/=1)
  if (item < nitems) then
    ! Si hay mas items, esta especificando un calculo de bond order nuevo, de
    ! lo contrario solo se imprime lo que ya esta guardado en la variable

    call readi(i1) ! The second group related with bo calculation

    ! Default values for cut radious of bond order function
    f1=4.0_dp  ! minor cut radious
    f2=5.0_dp  ! mayor cut radious
    if (item < nitems) call readf(f1)  ! User value for minor cut radious
    if (item < nitems) call readf(f2)  ! User value for mayor cut radious

    ! Run selecction according to the second group selected
    if (i1==0) then
      call bond_order_groups(gsel,sys,f1,f2)
    else
      call bond_order_groups(gsel,gr(i1),f1,f2)
    endif

  end if
  call polvar_hard('_ans[1]',gsel%alist%next%o%border)
case('index') 
  call werr('Only one atom must be in the selection',gsel%nat/=1)
  call readi(i1) ! The group 
  la => gr(i1)%alist
  do i = 1,gr(i1)%nat
    la => la%next
    if (.not.associated(la%o,target=gsel%alist%next%o)) cycle
    call polvar_hard('_ans[1]',i)
    exit
  enddo 

case('rang')
  call rang(f1)
  call polvar_hard('_ans[1]',f1)
case('ranu')
  call polvar_hard('_ans[1]',ranu())
case('max')
  call readf(f1)
  call readf(f2)
  write(ans,fmt='(e15.8)') max(f1,f2)
  call polvar_hard('_ans[1]',max(f1,f2))
case('min')
  call readf(f1)
  call readf(f2)
  call polvar_hard('_ans[1]',min(f1,f2))
case('int')
  call readf(f1)
  call polvar_hard('_ans[1]',int(f1))
case('floor')
  call readf(f1)
  call polvar_hard('_ans[1]',floor(f1)) 
case('ceil','ceiling')
  call readf(f1)
  call polvar_hard('_ans[1]',ceiling(f1))   
case("box") 
  call polvar_hard('_ans[1]',box(1))
  call polvar_hard('_ans[2]',box(2))
  call polvar_hard('_ans[3]',box(3))
  ans=3


! Estos case es para cuando el resultado no se puede guardar en ans
!             
! case("pdistrad") 
!   call readf(fv)
!   call readf(f1)
!
!   call t%init_histo(inq_pdisrad(fv,gsel),f1)
!   call t%wprt()
!
case default
  call polvar_hard('_ans[1]',trim(w1))
endselect 


end function print_commands

  subroutine time_commands
    call readl(w1)
    selectcase(w1)
    case('cero')
      time=0.0_dp
      ptime=0.0_dp
      dm_steps=0.0_dp
      nframe=0.0_dp
      call wstd('time set to cero')
    case('step')
      call readf(dt)
      dtaux=dt
    case default  
      call wwan('I do not understand the last command')  
    endselect 
  endsubroutine time_commands
             
subroutine element_commands
use gems_elements, only: add_z 

call readl(w1)
selectcase(w1)
case('add')
  call reada(w1)
  call readf(f2)
  call add_z(trim(adjustl(w1)),f2,0._dp,1._dp)
case default  
  call wwan('I do not understand the last command')  
endselect 

endsubroutine element_commands

subroutine prng_commands
call readl(w1)
selectcase(w1)
case('lcg')
  call lcg_init()
case('std')
  call std_init()
case('seed')
  call readi(i1)
  call init_ran(lseed=i1)
case('seed_first')
#ifdef SPRNG
  call readi(i1)
  call readi(i2)
  call init_ran(lstreamnum=i1,lnstreams=i2)
#else
  call werr('Feature not enable, configure with --with-sprng')
#endif

case('stream')
#ifdef SPRNG
  call readi(i1)
  if (nitems>item) then
    call readi(i2)
    call spawn_ran(i1,i2)
  else
    call spawn_ran(i1)
  endif
#else
  call werr('Feature not enable, configure with --with-sprng')
#endif

case default  
  call wwan('I do not understand the last command')  
endselect 

endsubroutine prng_commands
      
  subroutine box_commands
    integer   :: i,j

    call readl(w1)
    selectcase(w1)
    case('move')
      call inq_cm_vel(sys)
      fv=-sys%cm_vel
      call set_add_cmvel(sys,fv)
    case('tsize')
      cubic=.false.
      do i = 1,dm
        call read_line(eof)
        if (eof) exit
        do j = 1,dm
          call readf(tbox(j,i))
        enddo
        tbox(:,i)=tbox(:,i)/box(i)
        tbox(i,i)=0.0_dp
      enddo
      call box_setvars()
    case('size')
      call werr('Only cubic boxes for now',.not.cubic)
      tbox(:,:)=0._dp
      call readf(tbox(1,1))
      call readf(tbox(2,2))
      call readf(tbox(3,3))
      call box_setvars()
    
    case('expand')
      call readf(f1)
      call readf(f2)
      call readf(f3)
      call box_expand(f1,f2,f3)
      w1='(A,'//cdm//'(f10.5))'
      call wstd(); write(logunit,trim(adjustl(w1))) '  box size:', (box(i),i=1,dm)
    case('mic') ! Establezco condiciones periodicas para ese grupo
      call readb(b1)
      mic=b1
      ghost=.not.b1
    case default  
      call wwan('I do not understand the last command')  
    endselect 

  endsubroutine box_commands

  subroutine table_commands
    type(etable) t1

    call reada(w1) ! Archivo
    call readi(i2) ! Columna x
    call readi(i3) ! Columna y
    call etable_read(t1,w1,i2,i3)

    call readl(w2) ! Comando 
    select case(w2)
    !case("cspline")
    !  call readi(i1)
    !  t2=cspline_t(i1*(t1%pnt-1)+t1%pnt,t1)
    !case("deriv3")
    !  call table_deriv_tres(t1,t2)
    !case("deriv5")
    !  call table_deriv_cinco(t1,t2)
    !case("fpplot")
    !  call etable_fpplot(t1,t2,t3)
    !  write(w2,*) trim(adjustl(w1))//'.var'
    !  call table_write(w2,t2) 
    !  write(w2,*) trim(adjustl(w1))//'.err'
    !  call table_write(w2,t3) 
    case default  
      call wwan('I do not understand the last command')  
    endselect 

    !if (nitems>item) then
    !  call readl(w2)
    !  call reada(w2)
    !  call write_simple_table(w2,t2) 
    !  deallocate(t1%y,t2%y)
    !endif

  endsubroutine table_commands

  subroutine mtable_commands
    !type(etable) t1,t2,t3

    call readl(w2) ! Comando 
    !select case(w2)
    !case("linreg")
    !  call readl(w1) ! Archivos
    !  call readi(i2) ! Columna x
    !  call readi(i3) ! Columna y    case("linreg")
    !  do i=1,size(chkey)
    !    call etable_read(t1,chkey(i),i2,i3)
    !    call etable_lreg(t1,f1,f2,f3)
    !    call wstd(); write(logunit,*) 'Lin.Reg for ',trim(adjustl(chkey(i))
    !    call wstd(); write(logunit,*) ' -Slope (m)',f1
    !    call wstd(); write(logunit,*) ' -Ycut (b)',f2
    !    call wstd(); write(logunit,*) ' -Correlation (r)',f3
    !    call wstd(); write(logunit,*) 
    !  enddo 
    !case("yprom")
    !  call readl(w1) ! Archivos
    !  call readi(i2) ! Columna x
    !  call readi(i3) ! Columna y
    !  call etable_read(t1,chkey(1),i2,i3)
    !  do i=2,size(chkey)
    !    call etable_read(t3,chkey(i),i2,i3)
    !    t2=ysum_t(t1,t3)
    !    if(i<size(chkey)) t1%y=t2%y
    !  enddo
    !  t2%y(:)=t2%y(:)/size(chkey)
    !case default  
    !  call wwan('I do not understand the last command')  
    !endselect 
    !
    !if (nitems>item) then
    !  call readl(w2)
    !  call reada(w2)
    !  call table_write(w2,t2) 
    !  deallocate(t1%y,t2%y)
    !endif

  endsubroutine mtable_commands
 
subroutine outfile_commands
  use gems_output
  use gems_hyperdynamics
  use gems_bias, only: write_bias
  use gems_forcefield
  use gems_tersoff
  use gems_interaction
  use gems_neb
  use gems_metadynamics, only: write_cvs,write_E_1D
  use gems_analitic_pot, only: write_halfsho
  class(outfile),pointer  :: of=>null()
  class(outpropa),pointer  :: op
  type(outpropa_l),pointer :: ln
  integer                 :: j

  ! Selecciono el archivo de salida
  call readi(i1)
  if (i1>ofiles_max) stop 'demasiados archivos de salida'
  of => ofiles(i1)

  call readl(w1)
    
  selectcase(w1)
  case('prom')
    of%prom=.true.
    call readi(of%outeach)
    call werr('Integer must be grater than 0',of%outeach<=0)
    of%w => outfile_prom

  case('ddda')
    of%ddda=.true.
    call readi(of%outeach)
    call werr('Integer must be grater than 0',of%outeach<=0)
    call of%reopen()
    of%w => outfile_ddda
     
    ! werr declare columns first

    ! Force fields
    ln => of%p
    do while( associated(ln%next) )
      ln => ln%next

      if(.not.associated(ln%o%f)) cycle
      if(associated(ln%o%d)) cycle

      allocate(ln%o%d(size(ln%o%f)))
      do j = 1,size(ln%o%d)
        call ln%o%d(j)%init()
      enddo
      
    enddo 
 
  case('each')
    call readi(of%sumeach)

  case('flush')
    call reada(w2)
    selectcase(w2)
    case('on')
      of%flush=.true.
    case('off')
      of%flush=.false.
    case default
      call wwan('I do not understand the last command')  
    endselect

  case('name')
    call reada(w2)
    of%name = trim(adjustl(w2))
    call of%reopen()

  case('off')
    of%enable=.false.

  case('on')
    call encode_logicalvector(of%enable,of%ecode)

  case('at')
    of%enable=.false.

    ! Para indicar en que momento (algoritmo) debe tratar de escribirse el archivos
    call readl(w1)
    select case(w1)
    case ('common' )
      of%enable(1)=.true. 
    case ('eparts' )
      of%enable(2)=.true.
    case ('forcefield' )
      of%enable(3)=.true.
    case ('tersoff' )
      of%enable(4)=.true.
    case ('hd' )
      of%enable(5)=.true.
    end select
  
    of%ecode=decode_logicalvector(of%enable)
    
  case('cols')
    

    of%g=>null()

    ! Inicializo
    call of%reopen()
    !TODO:! Destruir of%p
    ! call of%p%init()

    !Adding a node, 
    !FIXME: using add_hardcpy(op) didnt work with a very strange error

    ! ANTES
      !case ('energy'    ); call outprop_init(op_energy   ,dm    ,'Energy.'   ,'.dat' ,grout(i1))
      !case ('border'    ); call outprop_init(op_border   ,dm    ,'Border.'   ,'.dat' ,grout(i1))
      !case ('nebi'      ); call outprop_init(neb_opi     ,dm    ,'NebI.'     ,'.xyz' ,grout(i1))
      !case ('nebf'      ); call outprop_init(neb_opf     ,dm    ,'NebF.'     ,'.xyz' ,grout(i1))
      !case ('dihedrals' ); get_dihedral  =.true.
      !! Tsf
      !case ('tsfangs'           ); get_tsfangs   =.true.
 
    ! Leo las propiedades y las asigno al archivo
    do while (item<nitems)
      call reada(w1)

      ln=>of%p%append_hard()
      allocate(ln%o)
      op=>ln%o
          
      ! Read (or not) the group ID
      select case(w1)
      case('eparts','time','step','bias','box')
        ! These does not require a selection group
      case('e_tors','e_bend','e_stretch','e_neb','neb_info','e_halfsho')
      case  default
        call werr('not group selected',item==nitems)
        call readi(i2)
      end select
     
      select case(w1)
      case('eparts'   );  call op%init(igr_vop%size,write_eparts)
      !case('dihedrals'); call op%init(,write_dihedrals)
      !case('angulos'  ); call op%init(,write_tsfangs)
      case('time'     ); call op%init(1,write_time)
      case('step'     ); call op%init(1,write_step)
      case('bias'     ); call op%init(3,write_bias)
      case('box'      ); call op%init(3,write_box)
      case('e_tors'   ); call op%init(1,write_etors)
      case('e_bend'   ); call op%init(1,write_ebend)
      case('e_halfsho'); call op%init(2,write_halfsho)
      case('e_stretch'); call op%init(1,write_estretch)
      case('e_neb'    ); call op%init(1,write_eneb)
      case('neb_info' ); call op%init(1,write_eneb)
      case('tmoment'   ); call op%init(3,write_tmoment       ,gr(i2))
      case('amoment'   ); call op%init(3,write_amoment       ,gr(i2))
      case('angvel'    ); call op%init(3,write_angvel        ,gr(i2))
      case('pressure'  ); call op%init(dm*dm,write_pressure  ,gr(i2))
      case('virial'    ); call op%init(dm*dm,write_virial    ,gr(i2))
      case('girrad'    ); call op%init(1,write_girrad        ,gr(i2))
      case('cvs'       ); call op%init(6,write_cvs           ,gr(i2))
      case('inercia'   ); call op%init(dm*dm,write_inercia   ,gr(i2))
      case('covariance'); call op%init(dm*dm,write_covariance,gr(i2))
      case('mainaxis'  ); call op%init(3,write_mainaxis      ,gr(i2))
      case('ptriaxial' ); call op%init(1,write_ptriaxial     ,gr(i2))
      case('globalerr' ); call op%init(2,write_globalerror   ,gr(i2))
      case('energy'    ); call op%init(3,write_energy        ,gr(i2))
      case('absenergy' ); call op%init(3,write_absenergy     ,gr(i2))
      case('epot'      ); call op%init(1,write_epot          ,gr(i2))
      case('ecin'      ); call op%init(1,write_ecin          ,gr(i2))
      case('absecin'   ); call op%init(1,write_absecin       ,gr(i2))
      case('energypa'  ); call op%init(3,write_energypa      ,gr(i2))
      case('aenergy'   ); call op%init(3,write_aenergy       ,gr(i2))
      case('temp'      ); call op%init(1,write_temp          ,gr(i2))
      case('tempall'   ); call op%init(3,write_tempall       ,gr(i2))
      !case('ctime'     ); op%w => write_ctime      ; n =1))   ; op%adv = .true. ; b1=.false.  
      ! call wwan('Just remember that ctime must be load each 1 steps to have meaning')
      !case('hd_f'      ); op%w => write_f          ; n =1))   ; op%adv = .false.; b1=.false.  
!      case ('nebi'      ); ne%w => writeb_opi
!      case ('nebf'      ); ne%w => writeb_opf       
      case  default
        call werr('incorrect imput')
      end select

      ! call of%p%add_hardcpy(op)
      ! call op%destroy()
      of%w => outfile_write

    enddo

  case default

    b1=.true.

    !TODO of%p%destroy()
    of%ddda = .false.
    of%prom = .false.

    !xyz
    select case(w1)
    case('pos'        ); of%w => write_pos     
    case('free_en_1d' ); of%w => write_E_1D     
    case('poscr'      ); of%w => write_poscr
    case('hd_fpp'     ); of%w => write_fpp     ; b1=.false.  
    case('vel'        ); of%w => write_vel     
    case('vel_rot'    ); of%w => write_vel_rot 
    case('vel_vib'    ); of%w => write_vel_vib 
    case('fce'        ); of%w => write_fce     
    case('pes'        ); of%w => write_pes     
    case('pose'       ); of%w => write_pose
    case('dist'       ); of%w => write_dist    
    case('charge'     ); of%w => write_charge
   !%case ('border'   ); of%w => write_border    
    case  default
      call werr('incorrect imput')
    end select
    
    call of%reopen()

    if(b1) then
      call werr('not group selected',item==nitems)
      call readi(i2)
      of%g => gr(i2)
    endif
        
    call werr('many inputs. this file not admit mor columns',item<nitems)

  endselect

  endsubroutine outfile_commands

end module gems_clinterpreter


