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
! References
!
! Paz & Leiva (2015). Time Recovery for a Complex Process Using Accelerated
! Dynamics. JCTC 11 1725. http://doi.org/10.1021/ct5009729
!
! Reigada et. al. (1999). One-dimensional arrays of oscillators: Energy
! localization in thermal equilibrium. JCP 111(4) 1373.
! http://doi.org/10.1063/1.479396  
 
module gems_hyperdynamics
use gems_program_types
use gems_constants
use gems_inq_properties 
use gems_set_properties
use gems_quasi_newton, only: lbfgs_minimizator
use gems_integration, only: integration_stepa, integration_stepb 
use gems_output
use gems_errors
use gems_programs
use gems_checkpoint
use gems_neighbor
use gems_bias, only: compress_below

implicit none
private

! Intento para el corte
integer              :: bsteps=0,nbsteps=0
real(dp)             :: maxbias,bias_highc,bias_lowc,btime,nbtime
logical              :: over_lowc=.false.,over_highc=.false.
real(dp)             :: sigma,vprom

type(group),pointer   :: ghd=>null()

! lp hyperdinamics
type(decorrelation)           :: dvp,dbp,dvbp
real(dp)                      :: f
integer                       :: ctime

type(compress_below),pointer  :: igb=>null()

real(dp)    :: hd_p1, hd_p2, desp

!  public :: hyperd_escape
public :: hyperd_eprom
public :: hybrid_ea
public :: wb2ea
public :: write_fpp
public :: write_f
public :: write_ctime
                        
contains

! subroutine hyperd_cli(it,w)
!
! end subroutine hyperd_cli(it,w)

function boostfactor(x) result(f)
real(dp)            :: f
real(dp),intent(in) :: x
!integer             :: i

f = hd_p1*exp(x*(hd_p2-x/2._dp))/cdf_snorm(hd_p2-x)
end function 
            
function boostfactor2(eprime,aprime) result(f)
! Compute boost factor according to equation 15 of (Paz2015)
real(dp)            :: f
real(dp),intent(in) :: eprime,aprime

f=1._dp-cdf_snorm(eprime)
f=f+exp(-aprime*(eprime-aprime/2._dp))*cdf_snorm(eprime-aprime)
f=1._dp/f
end function 

! subroutine hyperdinamics_cli()
! use gems_bias, only: compress_below
! use gems_neighbor, only: ngroup_dl
! use gems_interaction, only: polvar_interact(w1)
! class(ngroup), pointer      :: ig
! class(compress_below), pointer  :: p
!
! selectcase(w1)
! case('bias') 
!   call reada(w1)
!   igb=>polvar_interact(w1)
!   select type(igb)
!   type is(compress_below)
!   case default
!     call werr('Bias shoulw be of type lpe to allow hybrid HD-DM algorithm')
!   end select
! case default
!   call werr('Bad keyword. Use `under`, `with` or `feels`')
! end select
!
! end subroutine
!
     
! METODO HIBRIDO HD_DM

subroutine hybrid_ea(nsteps,msteps,lsteps,hd_f,eprime,aprime,z,t,b_out)
! eprom es la energia de referencia de E (i.e. enería promedio del sistema)
! lp1 y lp2 deben estar establecidos
use gems_fields
use gems_input_parsing
use gems_tb
use gems_bias, only: compress_below
use gems_interaction, only: polvar_interact, ngroup
logical,intent(in)        :: b_out
integer,intent(in)        :: nsteps,msteps,lsteps ! numero de bloques, pasos por bloque, equilibracion
real(dp),intent(in)       :: hd_f,eprime,aprime,z,t
integer                   :: i,n
integer,save              :: uhd=0,udm=0
real(dp)                  :: bprom,vbprom,prob,dumy,zfact!,dvbmin,aux
logical                   :: acelerar=.false.!,dummy
class(ngroup),pointer :: ig


call dvp%init()
call dbp%init()
call dvbp%init()
    
kt=(kB_ui*t) ! The temperature asociated with the mcpot part of the energy
beta = 1._dp/kt
acelerar=.false.
maxbias=0._dp
btime=0._dp
nbtime=0._dp

! Get igb (bias) association
ig=>polvar_interact(':hybrid_ea')
call werr('An interaction with `:hybrid_ea` label is required',.not.associated(ig))
! call werr('',.not.same_type_as(igb,igb))
select type(ig)
type is (compress_below)
  igb=>ig
class default
  call werr(':hybrid_ea should be of type compress_below')
end select

! Ensure the parameters to be inside the range
igb%alpha=boostfactor2(eprime,aprime)
igb%e=1._dp-igb%alpha+igb%alpha*cdf_snorm(eprime)

call wlog('HYD'); write(logunit,*) 'El boost factor será:', igb%alpha
call wlog('HYD'); write(logunit,*) 'y la integral w será:', igb%e
call flush(logunit)

zfact=inv_cdf_snorm(igb%e)-inv_cdf_snorm(eprime*z)

! Pongo parametros para que el bias sea cero
igb%alpha = 1._dp
igb%e = -1e16
bias_highc = 1._dp


if(uhd==0) then
  uhd=find_io(30)
  open(uhd,file=trim(ioprefix)//'.hd')
  udm=find_io(30)
  open(udm,file=trim(ioprefix)//'.dm')
endif

do i=1,nsteps
                  
  ! Corro dinamica o hyperdinamica
  call hyperd_eprom(msteps,b_out)

  !Calculo de cantidades interesantes

  vprom=dvp%med()                      !<V>b

  call dvp%var(n,dumy)

  !Tiene que ser -1... nos e porque anda mejor el -2
  f=float(n)/(dvp%size-2)              ! fraccion plato

  bprom=dbp%med()                      ! <DVb>b   o     <0>
  
  vbprom=dvbp%med()                    ! <Vb>b    o     <V>

  sigma=dvbp%disp()                    ! s(<Vb>b) o    s(<V>)

  ! vmin=lbfgs_getenergy(gr_erm,.false.) ! vmin

  !call bias_func(dvbmin,fdummy,vmin,dummy)    ! bias en vmin dvbmin

  prob=float(bsteps)/float(nbsteps+bsteps)  ! fraccion de pasos con bias

  ! Escribo el output
  write(udm,fmt='(i0,9(1x,e25.12))') i,dm_steps,time,f,vprom*ui_ev,sigma*ui_ev

  ! Lisent to term signal
  if (term_signal) exit

  if(acelerar) then ! Vengo de DA

    !if(.not.over_lowc) then ! Estoy debajo del bfact deseado
    !  if(f>hd_f) then ! Me canse, aumento HD
    !    !e = inv_cdf_norm(eprime,sigma=sigma) ! E-<Vb>b
    !    !e = (e-(vmin+dvbmin)*(1._dp-a)+vbprom)/a
    !    !aux = vmin + (e-vmin)*(a/(a-daumento))
    !    !e=aux 
    !    !a = a-daumento
    !    !bias_highc = inv_cdf_norm(z,sigma=sigma*afact,mu=(E-vprom)*afact)
    !    bias_highc = 1e10
    !    initial = .true.
    !  endif
    !endif

    write(uhd,fmt='(i0,9(1x,e25.12))') &
      i,dm_steps,igb%e*ui_ev,igb%alpha,prob,vbprom*ui_ev, &
      bprom*ui_ev,sigma*ui_ev,maxbias*ui_ev,bias_highc*ui_ev

    if(over_highc) then ! Demasiado bias apago la HD
      acelerar=.false.
      igb%e=-1.e16_dp
      igb%alpha=1._dp
      over_lowc=.false.
      over_highc=.false.
      maxbias=0._dp
      write(uhd,*)
    endif    
     
  else ! Vengo de DM

    ! TODO si es que vengo de DM hago el test de cansancio. Si vengo de DA no lo
    ! hago.
    ! if( e = -1e16) acelerar=(f>hd_f)
    acelerar=(f>hd_f)

    if(acelerar) then ! Me canse, pasa a HD

      ! Establezco el a y e inicial
      igb%e=sigma*eprime+vbprom       
      igb%alpha=1._dp-aprime/(beta*sigma) 

      !bias_highc = (inv_cdf_snorm(hd_phi)*sigma+(hd_e-vprom))*(1._dp-hd_a)/hd_a
      bias_highc = (1._dp-igb%alpha)*sigma*zfact
      !bias_highc = 1e10

      ! Put a blank line to plot with discontinuities
      write(udm,*)

      ! Hago la HD de inicializacion
      write(uhd,fmt='(i0,9(1x,e25.12))') i,dm_steps,igb%e*ui_ev,igb%alpha,prob,&
            vbprom*ui_ev,bprom*ui_ev,sigma*ui_ev,maxbias*ui_ev,bias_highc*ui_ev
      call hyperd_init(lsteps,.false.)
      write(uhd,fmt='(i0,9(1x,e25.12))') i,dm_steps,igb%e*ui_ev,igb%alpha,prob,&
            vbprom*ui_ev,bprom*ui_ev,sigma*ui_ev,maxbias*ui_ev,bias_highc*ui_ev

      over_lowc=.false.
      over_highc=.false.

    else  ! soy paciente

      igb%alpha = 1._dp
      igb%e = -1e16

    endif

  endif 
            
  call dvp%empty()
  call dbp%empty()
  call dvbp%empty()

enddo

! FIXME Esto es porque las variables no se calculan cada vez en el write y queda
! el numero trabado
! bias = 0._dp
! biased = 0._dp
! bfact = 1._dp 

call dvp%dest()
call dbp%dest()
call dvbp%dest()
          
end subroutine hybrid_ea

subroutine hyperd_init(steps,b_out)
use gems_programs, only: dinamic
! Realiza `steps` pasos de HD para termalizar la configuracion inicial en el
! potencial con bias
integer                       :: steps
logical,intent(in)            :: b_out
real(dp)                      :: auxtime

auxtime = time
maxbias=0._dp

call dinamic(steps,b_out,.false.) 

time = auxtime

end subroutine
                
subroutine hyperd_eprom(steps,b_out)
! Realiza steps pasos de DM o HD acumulando la energia media y demas variables
! estadisticas
use gems_bias, only: biased
use gems_integration, only: integration_reversea
use gems_input_parsing, only: load_blk, execute_block, bloques
use gems_ddda
integer,intent(in)  :: steps
logical,intent(in)  :: b_out
integer             :: ns!,ct
real(dp)            :: aux

call interact(.true.) 

! Reseteo los contadores y flags
nbsteps = 0
bsteps = 0

do ns = 1,steps

  dm_steps=dm_steps+1._dp

  call dvp%add(biased)
  call dbp%add(igb%epot)
  call dvbp%add(biased+igb%epot)
      
  call integration_stepa()
  call interact(b_out) 
  
  ! Me fijo si es muy zarpado el boost antes de seguir
  maxbias=max(maxbias,igb%epot)

  if(.not.over_lowc) over_lowc=(igb%epot>=bias_lowc)

  if(igb%epot>=bias_highc) then
    ! Agrego al alarma y corto esta corrida
    over_highc=.true.
    ! Vuelvo atras este paso.
    call integration_reversea()
    return
  endif

  aux=time
  call integration_stepb()

  if(igb%epot>0) then
    bsteps = bsteps+1
    btime = btime+time-aux
  else
    nbsteps = nbsteps+1
    nbtime = nbtime+dt
  endif
   
  !! Cuento el intervalo de tiempo hasta que el bias toca cero de nuevo
  !if(b_bias)  then
  !  ct = ct + 1
  !else
  !  ctime = ct
  !  ct = 0
  !endif
          
  !Escribo aca para que coincida interacción y configuracion
  if (b_out) call write_out(1,dm_steps)

  ! Checkpoint
  if (b_ckp) then
    if (mod(dm_steps,real(chpeach,dp))==0._dp)  call write_chp(ns,steps)
  endif
       
  ! Command interpreter
  if (b_load) then
    if (mod(dm_steps,real(load_each,dp))==0._dp)  call execute_block(bloques(load_blk),1,1,1,.false.)

    ! Lisent to term signals
    if (term_signal) exit

  endif

enddo

! Escribo el grafico de fp
if (b_out) call write_out_force(5)
   
end subroutine hyperd_eprom

subroutine wb2ea(w,b,e,alpha,n)
! convierte w y b en aprime y eprime
! n es el numero de grados de libertad del termino a boostear
! esto debe ser un input porque los numeros de atomos involucrados en un
! potencial no necesariamente reflejan el grado de libertad. Por ejemplo, un
! dihedro necesita 4 atomos pero es solo un dihedro. XXX: pensar bien esto
use gems_fields
real(dp), intent(in)        :: w, b
integer, intent(in)         :: n
real(dp), intent(out)       :: e, alpha

! Ensure the parameters to be inside the range
call werr('w must be between 0 and 1',w>1._dp.or.w<0._dp)
call werr('Really? boost factor lower than 1??',b<1._dp)

! Saving w
hd_p1=w

! Saving eprime, obtained from equation 15 for w in (Paz2015)
hd_p2 = inv_cdf_snorm((b+w-1._dp)/b)

! The of the energy variance of an harmonic oscilator is (kT)**2 (Reigada1999).
! From the equipartition of the energy, this should arise from the variance of
! Ecin and Epot, so the variance of Epot should be (kT)**2/2 and considering n
! oscilators this will be n(kT)**2/2 so the maximum value of alpha prime will be 
e = sqrt(0.5*n)

alpha = boostfactor(e)

call wlog('HYD');write(logunit,*) 'Maximum boostfactor for this w is near to:', alpha
call werr('The boost factor is too near to the maximum allowed',alpha<=b)

! Establezco el alpha inicial
alpha = inverte_bifuncrr(b,0._dp,e,boostfactor)
e = hd_p2

call wlog('HYD');write(logunit,*) 'eprime is:', e
call wlog('HYD');write(logunit,*) 'aprime is:', alpha

end subroutine wb2ea

! OTROS

! subroutine hyperd_escape(nsteps,e,m,nmed,lmed,b_out)
! use gems_programs 
!
! logical,intent(in)          :: b_out
! real(dp),intent(in)         :: e,m
! integer,intent(in)          :: lmed,nmed,nsteps
! real(dp)                    :: desv,emin
! real(dp),target             :: pmin(ghd%nat*dm)
! integer                     :: i,flag
! logical                     :: switched
!
! !FIXME
!
! ! ghd es el grupo que se va a mirar para ver cuando escapa
!
! !Hypervector de fuerza
! allocate(hd_fce(natoms*dm))
!  
! ! Guardo el minimo en pmin
! call lbfgs_minimizator(ghd,.false.,pmin,emin)
! kt=(kB_ui*ghd%fixtemp)
! hd_e=e
! hd_a=(hd_e-emin)*m*sqrt(kt)/((hd_e-emin)-m*sqrt(kt))
! beta = 1._dp/kt
! desp=1.5_dp
!
! do i = 1,nsteps
!   call dinamic(nmed,.false.) ! Descorrelacion
!   flag=hyperd_eprom_safe(nmed,lmed,hd_eprom,desv,.false.)  
!   if(flag==0) then
!     call wstd(); write(logunit,*) 'TRAPED! BOOST ON:-----------------------'
!     call wstd(); write(logunit,*) '  frame:',frame,' time:',time
!     call wstd(); write(logunit,*) '  eprom:', hd_eprom*ui_ev,'+-',desv*ui_ev
!     hd_e = hd_eprom+hd_p1
!     hd_a = hd_e*hd_p2
!     do while(hyperd_dinamic_safe(nmed,lmed,hd_a,b_out))
!       call wstd(); write(logunit,*) 'MORE BOOST:'
!       !call wstd(); write(logunit,*) '  frame:',frame,' time:',time
!       !call wstd(); write(logunit,*) '  eprom:',media*ui_ev,'+-',desviacion*ui_ev
!       hd_e=e+hd_a
!       hd_a=(hd_e-emin)*m*sqrt(kt)/((hd_e-emin)-m*sqrt(kt))
!     enddo
!   endif
!
!   call wstd(); write(logunit,*) 'MOVING! BOOST OFF:----------------------'
!   bias._dp
!   call write_out_force(1)
!
!   ! Guardo el minimo en pmin
!   call lbfgs_minimizator(ghd,.false.,pmin,emin)
!   call wstd(); write(logunit,*) 'Jump!:', sqrt(dot_product(pmin-ghd%pp,pmin-ghd%pp))
! enddo
!
! ! Retomo el modo de almacenamiento anterior
! if(switched) call group_switch_objeto(ghd)
! deallocate(hd_fce)
!        
! end subroutine hyperd_escape

function hyperd_dinamic_safe(m,n,media,b_out)
! ghd es el grupo que se utiliza en esta dinamica
! eprom es la energia de referencia de E (i.e. enería promedio del sistema)
use gems_bias, only: biased
use gems_fields
use gems_tb

logical                     :: hyperd_dinamic_safe
logical,intent(in)          :: b_out
integer,intent(in)          :: m,n
real(dp)                    :: medida,suma,suma2
real(dp)                    :: timeaux,desviacion
integer                     :: i,j
real(dp),target             :: pmin(ghd%nat*dm)
!real(dp),target             :: fce(ghd%nat*dm),pos(ghd%nat*dm)
!real(dp)                    :: vmin
real(dp),intent(out)        :: media

hyperd_dinamic_safe=.true.

! Imprimiendo datos utiles
call wstd(); write(logunit,*) '  E', igb%e*ui_ev
call wstd(); write(logunit,*) '  alpha', igb%alpha*ui_ev

!Corro la dinamica
timeaux=0

suma=0
suma2=0

call interact(b_out) 

do j=1,n

  medida=0._dp

  do i=1,m

    dm_steps = dm_steps + 1._dp

    !Acumulo para promedio y desviacion
    medida=medida+biased+igb%epot

    ! Avanzo
    call integration_stepa
    call interact(b_out) 
    call integration_stepb

    ! Escribo aca para que coincida interacción y configuracion
    if (b_out) call write_out(1,dm_steps)
                       
  enddo

  medida=medida/m
  suma=suma+medida
  suma2=suma2+medida*medida
  print *, medida*ui_ev

  call lbfgs_minimizator(ghd,.false.,pmin)
  if (diff_vect(ghd%pp,pmin,desp))  then
    hyperd_dinamic_safe=.false.
    return
  endif

enddo

media=suma/n
desviacion=sqrt(suma2/n-media*media)

end function hyperd_dinamic_safe
 
! function hyperd_eprom_safe(m,n,media,desviacion,b_out) result(flag)
! ! Determina la energía potencial promedio del termino a boostear y la
! ! desviacion estandar. Utiliza el grupo de ermak. 
! ! m Numero de pasos para determinar una medida de energía potencial 
! ! n Numero de medidas para determinar la media
! use gems_programs
! use gems_bias, only: noboost, biased
! real(dp),intent(out)        :: media,desviacion
! integer                     :: i,j
! integer                     :: flag
! real(dp),target             :: pmin(ghd%nat*dm)
! integer,intent(in)          :: m,n
! real(dp)                    :: medida,suma,suma2
! logical,intent(in)    :: b_out
!
! suma=0
! suma2=0
! flag=1
!
! ! No quiero boostear las fuerzas aca
! noboost=.true.  
!
! ! Guardo el minimo en pmin
! call lbfgs_minimizator(ghd,.false.,pmin)
! call interact(.false.)
!   
! do j = 1,n
!
!   medida=0._dp
!   do i = 1,m
!     dm_steps=dm_steps+1._dp
!     
!     !Acumulo para promedio y desviacion
!     medida=medida+biased
!
!     !Dinamica Step
!     call integration_stepa
!     call interact(.false.)
!     call integration_stepb
!
!     !Avance del tiempo
!     time = time + dt
!
!     ! Escribo aca para que coincida interacción y configuracion
!     if (b_out) call write_out(1,dm_steps)
!   enddo
!
!   ! Salgo porque no estoy en una basija
!   call lbfgs_minimizator(ghd,.false.,pmin)
!   if (diff_vect(ghd%pp,pmin,desp)) return
!
!   medida=medida/m
!   suma=suma+medida
!   suma2=suma2+medida*medida
! enddo
!
! media=suma/n
! desviacion=sqrt(suma2/n-media*media)
!  
! ! Sin errores
! flag=0
!
! end function
!
! SALIDA

subroutine write_fpp(of)
! Realiza el grafico de Flyvbjerg and Petersen: desviacion estandar en funcion
! de tamaño de bloque estadistico. El plato en el grafico indica que el tamaño
! de bloque consigue perder la correlacion
use gems_ddda
use gems_output
class(outfile)                  :: of
integer                         :: j,n,m
class(statistic_dclist),pointer :: ls
real(dp)                       :: err,errerr,med,plato!,bigerr


!call dvp%var(n,var)
!write(un,fmt='(i0,a,i0,f5.2)') n,' puntos de ',dvp%size-2,float(n)/(dvp%size-2)
!write(un,fmt='(4(a,e25.12))')  'media: ',dvp%med(),'+-',var,'dispersion: ',sqrt(dvp%disp())

! El siguiente do hace log_2(N)-1 lineas, donde N es el numero de puntos
! agregados en el d. El menos uno es porque para calcula var se divide en
! N-1 y esto provocaria un NAN 

n=dvp%blocks%o%nsamples
med=dvp%med()

call dvp%var(m,plato)

! Separador de bloque (Enable `plot "file" i 3` in gnuplot)
write(of%un,*)
write(of%un,*)

!Tiene que ser -1... nos e porque anda mejor el -2
!El cero tambien cuenta
ls =>dvp%blocks%prev
do j =1,dvp%size-2
  ls => ls%next
  ! var = ls%o%vars()

  !Ver DDDA decorrelation_variance
  err = sqrt(ls%o%vars()/(ls%o%nsamples-1))*ui_ev
  errerr=err/(sqrt(2._dp*(ls%o%nsamples-1)))

  !write(of%un,fmt='(i0,x,2(e25.12,2x),i0)')  dm_steps+j*(n/(dvp%size-2)),sqrt(var)*ui_ev+med*ui_ev,sqrt(err)*ui_ev,ls%o%nsamples
  write(of%un,fmt='(i10,1x,3(e25.12,2x))')  ls%o%nsamples,err,errerr,plato*ui_ev
enddo

!ATENCION: esto se lee asi
!Ejemplo 2. Supongamos que tengo 18 datos. 
! nsize=5 pero tengo 6 nodos(contando la cabeza)
!  Estas columnas no estan  |  Estas si
!    id    tama bloque      !   nsamples     error  errorerror
!    1 (head)    1          !    18        
!    2           2          !     9        
!    3           4          !     4        
!    4           8          !     2        
!    5          16          !     1        
!    6 (prev)    -          !     0        

!Es util la siguiente representación
!nodos            -(head)   -   -   -   x  x
!nsamples que
!aumentan size              1   2   4   8  16
!size                       1   2   3   4  5 
!bloque              1      2   4   8  16  32


end subroutine write_fpp
  
subroutine write_f(op,un)!(op,name)
use gems_output
class(outpropa)                :: op
integer,intent(in)             :: un

write(un,fmt='(e25.12)',advance='no')  f

end subroutine write_f
   
subroutine write_ctime(op,un)!(op,name)
use gems_output
class(outpropa)                :: op
integer,intent(in)             :: un

! Cuento el intervalo de tiempo hasta que el bias toca cero de nuevo
if(ctime>0) write(un,fmt='(e25.12)')  ctime*dt
ctime=0

end subroutine write_ctime
  
end module gems_hyperdynamics

