! Copyright (c) 2020  Lucas Farigliano
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


module gems_metadynamics
use gems_constants,          only: dp,pi,ev_ui,dm,kB_ui
use gems_groups,      only: group,atom,atom_dclist
use gems_strings, only: str
use gems_integration,     only: integration_stepa, integration_stepb, its, integrate
use gems_interaction,     only: interact
use gems_output
use gems_errors
use gems_neighbor, only: ngroup
use gems_checkpoint


implicit none

integer,parameter    :: kjau_pro= 100
real(dp)             :: kBuCM = 1.0E6
real(dp)             :: kBuRg = 1.0E6
real(dp)             :: kBu   = 1.0E6
real(dp)             :: dT_WT,rburb,tauWT
real(dp)             :: rcmA(3),rcmB(3),rcm(3)
integer              :: printD,cant_parta,tauwt_pasos
!variables colectivas
real(dp)             :: dCM,Rg,RgAu,Rgtotal,pos_x
!Wt 2D
real(dp)             :: wWTini,sigCV1,sigCV2,dpot_CV1,dpot_CV2
real(dp)             :: potiniCV1,potfinCV1,potiniCV2,potfinCV2,parA,parB
real(dp)             :: potiniaux1,potfinaux1,potiniaux2,potfinaux2
real(dp)             :: tol_Bias,deltapared,bias_ant,biasf 
integer              :: potbinCV2,potbinCV1
!WT 1D
real(dp),allocatable :: interpolation_x(:)
real(dp),allocatable :: pot_int_y(:)
real(dp),allocatable :: for_int_y(:)
real(dP),allocatable :: dfor_int_y(:)
!Apertura de archivo para las gaussianas 2d
integer              :: unidad

real(dp)       :: temp_md
type(integrate),pointer :: g   ! grupo del integrador para calcular temp 
real(sp), allocatable    :: B_WTMD_2D(:,:)
real(sp), allocatable    :: F_WTMD_2D_CV1(:,:)
real(sp), allocatable    :: F_WTMD_2D_CV2(:,:)
private
public         :: metadynamics
public         :: wtmd2D_set,dm_cv_set,wall2D_set,wallauxCore_set,wallauxShell_set
public         :: wtmetad_set,wall1D_set, bias_point_1D,Collective_Variable     
public         :: write_cvs,write_E_1D,dCM,posicion1d_set,wtxy_2D
 
! FIXME:
type(group),target,public   :: gmeta
 
contains


subroutine metadynamics(steps,b_out,b_wtmd_2D,b_DMCV,b_dCM_RgTotal,b_dCM_RgAu,b_dCM,b_pos1d,b_x,b_xy)
use gems_neighbor, only:nupd_vlist 
use gems_input_parsing, only:execute_block,load_blk, bloques
! use gems_errors, only: timer_start, timer_dump

! esta subrutina funciona con subsystemas
integer,intent(in)    :: steps
integer               :: ns
logical,intent(in)    :: b_out,b_wtmd_2D,b_DMCV,b_dCM_RgTotal,b_dCM_RgAu,b_dCM,b_pos1d,b_x,b_xy
 
! Timing
! call timer_start(time)
! Define the temperature globally

! Escribo aca para que coincida interacción y configuracion
if(b_wtmd_2D) call wtmetad_2D(0,steps) ! MetaD 2D con rgCore y cm
if(b_dCM_RgAu) call wtdcmrgau(0,steps) ! MetaD 2D con rgShell y cm
if(b_dCM_RgTotal) call dCM_RgTotal(0,steps) ! MetaD 2D con rgTotal y cm
if(b_dCM) call wtmetad(0,steps) ! MetaD 1D con cm
if(b_pos1d) call posicion1d(0,steps) ! MetaD 1D con cm
if(b_x) call wtx(0,steps) ! MetaD sobre coordenada x
if(b_xy) call wtxy_2D(0,steps) ! MetaD sobre coordenada x e y 
if(b_DMCV) call DM_CV(0)        ! MD con barrera guardando CV 
! if (b_out) call write_out(1,0)

do ns = 1,steps
  dm_steps=dm_steps+1._dp

  ! Integración
  call integration_stepa
  call interact(b_out)
  if(b_wtmd_2D) call wtmetad_2D(ns,steps)
  if(b_dCM_RgAu) call wtdcmrgau(ns,steps)
  if(b_dCM_RgTotal) call dCM_RgTotal(ns,steps)
  if(b_DMCV) call DM_CV(ns)
  if(b_dCM) call wtmetad(ns,steps)
  if(b_x) call wtx(ns,steps) ! MetaD sobre coordenada x
  if(b_xy) call wtxy_2D(ns,steps) ! MetaD sobre coordenada x e y 
  if(b_pos1d) call posicion1d(ns,steps) ! MetaD 1D con cm
  call integration_stepb

  ! Escribo aca para que coincida interacción y configuracion
  if (b_out) call write_out(1,dm_steps)

 ! Checkpoint
  if (b_ckp) then
    if (mod(dm_steps,real(chpeach,dp))==0._dp)  call write_chp(ns,steps)
  endif
   
  ! Command interpreter
  if (b_load) then
    if (mod(dm_steps,real(load_each,dp))==0._dp)  call execute_block(bloques(load_blk),1,1,1,.false.)

    ! Lisent to term signal
    if (term_signal) then
      term_signal=.false.
      exit
    endif
  endif
   
  ! ! Timming information
  ! if (b_time) call timer_dump(ns,time,nupd_vlist)

enddo

end subroutine
         
subroutine posicion1d_set(wt1,wt2,wt3,wt4,wt5,wt6)
use gems_strings, only: str
real(dp),intent(in)    :: wt1    ! Altura inicial de las gaussianas
real(dp),intent(in)    :: wt2    ! Ancho de las gaussianas
real(dp),intent(in)    :: wt3    ! Parametro de WTMD
real(dp),intent(in)    :: wt4    ! tau (frecuencia)
integer ,intent(in)    :: wt5    ! cada cuanto imprimir la CV
integer ,intent(in)    :: wt6    ! id group

!seleccion del grupo de integracion (cuales se mueven)
g => its%o(wt6)
call werr('Integration group should be in a constant T ensamble',.not.g%b_fixt)
temp_md=g%fixt

wWTini      = wt1
sigCV1      = wt2
dT_WT       = wt3
tauWT       = wt4
printD      = wt5
cant_parta  = 1

call wlog('WTMD1D_Pos', 'Altura inicial de la Gauss (KJ/mol)='//str(wWTini))
call wlog('WTMD1D_Pos', 'Ancho de medio pico de la Gauss (A)='//str(sigCV1))
call wlog('WTMD1D_Pos', 'Parametrop de la WTMetaD           ='//str(dT_WT))
call wlog('WTMD1D_Pos', 'Frecuencia de incorporacion (ps)   ='//str(tauWT))
call wlog('WTMD1D_Pos', 'Cantidad Particulas en el core     ='//str(cant_parta))
call wlog('WTMD1D_Pos', 'Cantidad Particulas en el core     ='//str(temp_md))
call wlog('WTMD1D_Pos', 'printCV                            ='//str(printD))

wWTini  = wWTini * real(kjau_pro,dp)
tauWT_pasos= int(tauWT/dt)
call ini_1D
end subroutine posicion1d_set
 subroutine pared_1d_max(pared)
  implicit none
   real(dp)                     :: ff1d,pared
   type (atom_dclist),pointer   :: la
   integer                      :: i
   ff1d = -kBuCM*(pos_x-pared)
   la => gmeta%alist
   do i = 1,cant_parta
    la => la%next
    la%o%force(1) = la%o%force(1) + ff1d
   enddo
  end subroutine pared_1d_max
  
 subroutine pared_1d_min(pared)
  implicit none
   real(dp)                     :: ff1d,pared
   type (atom_dclist),pointer   :: la
   integer                      :: i
   ff1d = -kBuCM*(pos_x-pared)

   la => gmeta%alist
   do i = 1,cant_parta
    la => la%next
    la%o%force(1) = la%o%force(1) + ff1d
   enddo
end subroutine pared_1d_min 

subroutine posicion1d(n,nstep)
type (atom_dclist),pointer   :: la
integer                      :: i
integer                      :: binCV1_1,binCV1_2,n,nstep
real(dp)                     :: ffcte

pos_x=0

! Apertura  de archivos y poner a 0 las variables del potencial
la => gmeta%alist
do i = 1,cant_parta
  la => la%next
  pos_x=la%o%pos(1) 
enddo

if(pos_x >= potiniCV1 .AND. pos_x <= potfinCV1) then 
  binCV1_1= int ((pos_x-potiniCV1) / real  (dpot_CV1,dp)) + 1
  binCV1_2= binCV1_1 + 1
  
  ! Cada tiempo tau incorporo una gaussiana a la matriz del bias y de las fuerzas
  if(mod (n,tauWT_pasos)==0)  call sumagauss1D(pos_x,binCV1_1,binCV1_2,n)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                      
  !Coeficientes para el calculo de las fuerzas para todos los atomos
  ffcte= sline_interpotation(pos_x,interpolation_x(binCV1_1),for_int_y(binCV1_1),&
      interpolation_x(binCV1_2),for_int_y(binCV1_2),dfor_int_y(binCV1_1),dfor_int_y(binCV1_2))
  !calculo de la fuerza del bias en la cordenada colectiva dCM
  la => gmeta%alist
  do i = 1,cant_parta
    la => la%next
    la%o%force(1) = la%o%force(1) + real (ffcte,dp)
  enddo
else 
  !Potencial externo sobre max dCM
  if (pos_x > (potfinCV1)) call pared_1d_max(potfinCV1)
  !Potencial externo sobre min dcm
  if (pos_x < (potiniCV1)) call pared_1d_min(potiniCV1)
endif

!!!!!!!!!!!!!!!!IMPRIMIR!!!!!!!!!!!!!!!!!!!!!!!
if(mod (n,printD)==0)  write(236,'(f12.3,f12.6)') real (n*dt,dp), pos_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if (n == nstep) then
!   open(UNIT=235,FILE= 'E_libre.dat')
!   fact=(-(temp_md+dT_WT)/(dT_WT))*0.01
!   do i=1, potbinCV1
!    write(235,*) interpolation_x(i),fact*pot_int_y(i)
!  enddo
!  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 end subroutine posicion1d
!!!!!!!!!!!!! Set de parametros para correr wtmd 1D 2D con o sin barreras!!!!!!!
!Parametros WTMD 2D
subroutine wtmd2D_set(wt1,wt2,wt3,wt4,wt5,wt6,wt7,wt8,wt9)
real(dp),intent(in)    :: wt1    ! Altura inicial de las gaussianas
real(dp),intent(in)    :: wt2    ! Ancho CV1
real(dp),intent(in)    :: wt3    ! Ancho CV2
real(dp),intent(in)    :: wt4    ! Parametro de WTMD
real(dp),intent(in)    :: wt5    ! tau (frecuencia)
integer ,intent(in)    :: wt6    ! cada cuanto imprimir las variables Colectivas
real(dp),intent(in)    :: wt7    ! tolerancia para las gaussianas         
integer ,intent(in)    :: wt8    ! cantidad de particulas en el core
integer ,intent(in)    :: wt9    ! id de grupo de integracion donde se monta la Meta
         
!seleccion del grupo de integracion (cuales se mueven)
g => its%o(wt9)
call werr('Integration group should be in a constant T ensamble',.not.g%b_fixt)
temp_md=g%fixt
wWTini      = wt1
sigCV1      = wt2
sigCV2      = wt3
dT_WT       = wt4
tauWT       = wt5
printD      = wt6
tol_Bias    = wt7
cant_parta  = wt8

call wlog ('WTMD2D'); write(logunit,*) 'w->',wWTini,'sig1->',sigCV1,'sig2->',sigCV2
call wlog ('WTMD2D'); write(logunit,*) 'dT->',dT_WT,'tau->',tauWT
call wlog ('WTMD2D'); write(logunit,*) 'printCV->',printD, 'temp', temp_md
call wlog ('WTMD2D'); write(logunit,*) 'tol->',tol_Bias,'N_partcore->',cant_parta

biasf= real (((temp_md+dT_WT)/(dT_WT)),dp)
parA= 1.0_dp/real (cant_parta,dp)
parB= 1.0_dp/real ((gmeta%nat-cant_parta),dp)
wWTini  = wWTini * real(kjau_pro,dp)
tol_Bias = tol_Bias * real(kjau_pro,dp)  
tauWT_pasos= int(tauWT/dt)
end subroutine wtmd2D_set
!Parametros Limites CV1
subroutine wall2D_set(wt1,wt2,wt3,wt4,wt5,wt6,wt7)
real(dp),intent(in)    :: wt1    ! potencial inicial CV1
real(dp),intent(in)    :: wt2    ! potencial final CV1
integer ,intent(in)    :: wt3    ! potencial bines CV1
real(dp),intent(in)    :: wt4   ! potencial inicial CV2
real(dp),intent(in)    :: wt5   ! potencial final CV2
integer ,intent(in)    :: wt6   ! potencial bines CV2
real(dp),intent(in)    :: wt7   ! rburb

potiniCV1   = wt1
potfinCV1   = wt2
potbinCV1   = wt3
potiniCV2   = wt4
potfinCV2   = wt5
potbinCV2   = wt6
rburb       = wt7

dpot_CV1= (potfinCV1 - potiniCV1) / real (potbinCV1,dp)
dpot_CV2= (potfinCV2 - potiniCV2) / real (potbinCV2,dp)

deltapared=0.2
call wlog('W2D', 'Ventana Well Tempered Metadynamics (1D)')
call wlog('W2D'); write(logunit,'(2(a,f12.4,2x),a,i6)') 'Pot_ini_CV1 (A) =',potiniCV1, & 
         'Pot_fin_CV1 (A) =',potfinCV1, 'Bins CV1 =',potbinCV1
call wlog('W2D'); write(logunit,'(2(a,f12.4,2x),a,i6)') 'Pot_ini_CV2 (A) =',potiniCV2, &
         'Pot_fin_CV1 (A) =',potfinCV2, 'Bins CV1 =',potbinCV2
call wlog('W2D'); write(logunit,'(2(a,f12.4),2x)') 'burb ext (A) =',rburb, &
         'delta pared (A)   =',deltapared

potbinCV1=potbinCV1+1
potbinCV2=potbinCV2+1

allocate(B_WTMD_2D(potbinCV1,potbinCV2))
allocate(F_WTMD_2D_CV1(potbinCV1,potbinCV2))
allocate(F_WTMD_2D_CV2(potbinCV1,potbinCV2))
  !Control, escribimos en el .log
call wlog ('W2D'); write(logunit,*) 'WCV1min->',potiniCV1,'WCV1max->',potfinCV1,'WCV1bin->',potbinCV1
call wlog ('W2D'); write(logunit,*) 'WCV2min->',potiniCV2,'WCV2max->',potfinCV2,'WCV1bin->',potbinCV2
call wlog ('W2D'); write(logunit,*) 'rburb->',rburb
end subroutine wall2D_set

!Parametros de barreras auxiliares para Rg Core 
subroutine wallauxCore_set(wt1,wt2)
real(dp),intent(in)    :: wt1    ! potencial inicial CV1
real(dp),intent(in)    :: wt2    ! potencial final CV1

potiniaux1   = wt1
potfinaux1   = wt2

call wlog('WAC', 'Ventana Aux Core Well Tempered Metadynamics')
call wlog('WAC'); write(logunit,'(2(A19,f12.4,2x))') 'Paux_ini Core (A) =',potiniaux1, &
    'Paux_fin Core (A) =',potfinaux1
end subroutine wallauxCore_set

subroutine wallauxShell_set(wt1,wt2)
real(dp),intent(in)    :: wt1    ! potencial inicial CV1
real(dp),intent(in)    :: wt2    ! potencial final CV1

potiniaux2   = wt1
potfinaux2   = wt2

call wlog('WAS','Ventana Aux Core Well Tempered Metadynamics')
call wlog('WAS'); write(logunit,'(2(A19,f12.4,2x))') 'Paux_ini Core (A) =',potiniaux2 ,'Paux_fin Core (A) =',potfinaux2
end subroutine wallauxShell_set

!!!Well Tempered Metadynamics 2D !!!!!
!dCM-RgCore
subroutine wtmetad_2D(n,nstep)
  type (atom_dclist),pointer   :: la
   integer                     :: n,nstep,i,binCV1_1,binCV1_2,binCV2_1,binCV2_2   
   real(dp)                    :: CV1_ffcte,CV2_ffcte,Cv1_ff(3),CV2_ff(3)
!
!Calculo mis variables colectivas en el paso n de la simulacion  
  call Collective_Variable

!Apertura de archivos de datos de la simulación e inicializacion!!
  if(n == 0) call ini_2D(n)

  if(dCM >= potiniCV1 .AND. dCM <= potfinCV1) then
    if(Rg >= potiniCV2 .AND. Rg <= potfinCV2) then
    
!Calculo en que lugar de la grilla del potencial y las fuerzas estoy parado
  binCV1_1= int ((dCM-potiniCV1)/ real (dpot_CV1,dp)) + 1
  binCV1_2= binCV1_1 + 1
  binCV2_1= int ((Rg-potiniCV2)/ real (dpot_CV2,dp)) + 1 
  binCV2_2= binCV2_1 + 1

 
  if(mod (n,tauWT_pasos)==0)  call sumagauss2D(dCM,Rg,binCV1_1,binCV1_2,binCV2_1,binCV2_2,n)
  
  CV1_ffcte= interpolacion_bilineal(dCM,Rg,binCV1_1,binCV1_2,binCV2_1,binCV2_2,&
      F_WTMD_2D_CV1(binCV1_1,binCV2_1),F_WTMD_2D_CV1(binCV1_2,binCV2_1),&
      F_WTMD_2D_CV1(binCV1_1,binCV2_2),F_WTMD_2D_CV1(binCV1_2,binCV2_2))
  CV2_ffcte = interpolacion_bilineal(dCM,Rg,binCV1_1,binCV1_2,binCV2_1,binCV2_2,&
      F_WTMD_2D_CV2(binCV1_1,binCV2_1),F_WTMD_2D_CV2 (binCV1_2,binCV2_1),&
      F_WTMD_2D_CV2 (binCV1_1,binCV2_2),F_WTMD_2D_CV2 (binCV1_2,binCV2_2))
 
 
  CV1_ff= (CV1_ffcte*(rcmA-rcmB))/ real (dCM,dp)
  CV2_ff = CV2_ffcte /real ((cant_parta*Rg),dp)

!calculo de las fuerzas de los bias 
  la => gmeta%alist
  do i = 1,cant_parta
   la => la%next
   la%o%force = la%o%force + (CV1_ff*parA)
   la%o%force = la%o%force + (CV2_ff*(la%o%pos-rcmA))   
  enddo
  do i=cant_parta+1, gmeta%nat
   la => la%next
   la%o%force= la%o%force - (CV1_ff*parB)
  enddo
  endif
  endif
!Potencial externo sobre max dCM
  if (dCM > (potfinCV1-deltapared)) call pared_dCM_max((potfinCV1-deltapared))
!Potencial externo sobre min dcm
  if ((dCM < (potiniCV1+deltapared)) .AND.(potiniCV1 /= 0.0)) call &
  pared_dCM_min((potiniCV1+deltapared))
  !Potencial externo sobre eg min  
  if (Rg < potiniCV2+deltapared) call pared_RgCo_min(potiniCV2+deltapared)
!Potencial externo sobre rg maxima
  if (Rg > potfinCV2-deltapared) call pared_RgCo_max(potfinCV2-deltapared)
!Burbuja sobre la posicion de las particulas
call pared_burb(n)

!!!!!!!!!!!!!!IMPRIMIR!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if(mod (n,printD)==0)  then
!  write(236,'(f12.3,4(f12.6))') real (n*dt,dp), dCM, Rg, RgAu, Rgtotal
! endif
 end subroutine wtmetad_2D
!dCM-RgShell
 subroutine wtdcmrgau(n,nstep)
  type (atom_dclist),pointer   :: la
   integer                     :: n,nstep,i,binCV1_1,binCV1_2,binCV2_1,binCV2_2   
   real(dp)                    :: CV1_ffcte,CV2_ffcte,Cv1_ff(3),CV2_ff(3)
!
!Calculo mis variables colectivas en el paso n de la simulacion  
  call Collective_Variable
!Calculo en que lugar de la grilla del potencial y las fuerzas estoy parado
  binCV1_1= int ((dCM-potiniCV1)/ real (dpot_CV1,dp)) + 1
  binCV1_2= binCV1_1 + 1
  binCV2_1=  int ((RgAu-potiniCV2)/ real (dpot_CV2,dp)) + 1 
  binCV2_2= binCV2_1 + 1


!Apertura de archivos de datos de la simulacion e inicializacion!!!!!!!!
  if(n == 0) call ini_2D(n)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(floor(float(binCV1_1)/potbinCV1)==0) then
    if(floor(float(binCV2_1)/potbinCV2)==0) then

                   
!Sumo las gauss dada una frecuencia
  if(mod (n,tauWT_pasos)==0)  call sumagauss2D(dCM,RgAu,binCV1_1,binCV1_2,binCV2_1,binCV2_2,n)
  
  CV1_ffcte= interpolacion_bilineal(dCM,RgAu,binCV1_1,binCV1_2,binCV2_1,binCV2_2,&
      F_WTMD_2D_CV1(binCV1_1,binCV2_1),F_WTMD_2D_CV1(binCV1_2,binCV2_1),&
      F_WTMD_2D_CV1(binCV1_1,binCV2_2),F_WTMD_2D_CV1(binCV1_2,binCV2_2))
  CV2_ffcte = interpolacion_bilineal(dCM,RgAu,binCV1_1,binCV1_2,binCV2_1,binCV2_2,&
      F_WTMD_2D_CV2(binCV1_1,binCV2_1),F_WTMD_2D_CV2 (binCV1_2,binCV2_1),&
      F_WTMD_2D_CV2 (binCV1_1,binCV2_2),F_WTMD_2D_CV2 (binCV1_2,binCV2_2))
 
          
  CV1_ff= (CV1_ffcte*(rcmA-rcmB))/ real (dCM,dp)
  CV2_ff = CV2_ffcte /real (((gmeta%nat-cant_parta)*RgAu),dp)

!calculo de las fuerzas de los bias 
  la => gmeta%alist
  do i = 1,cant_parta
   la => la%next
   la%o%force = la%o%force + (CV1_ff*parA)
  enddo
  do i=cant_parta+1, gmeta%nat
   la => la%next
   la%o%force= la%o%force - (CV1_ff*parB)
   la%o%force = la%o%force + (CV2_ff*(la%o%pos-rcmB))   
  enddo
  endif
  endif

  
  !Potencial externo sobre max dCM
  if (dCM > potfinCV1) call pared_dCM_max(potfinCV1)
  !Potencial externo sobre min dcm
  if (dCM < potiniCV1) call pared_dCM_min(potiniCV1)
!Potencial externo sobre eg min  
  if (RgAu < potiniCV2) call pared_RgShell_min(potiniCV2)
!Potencial externo sobre rg maxima
  if (RgAu > potfinCV2)  call pared_RgShell_max(potfinCV2)
!Burbuja sobre la posicion de las particulas
  call pared_burb(n)
!!!!!!!!!!!!!!IMPRIMIR!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if(mod (n,printD)==0)  then
!  write(236,'(f12.3,4(f12.6))') real (n*dt,dp), dCM, Rg, RgAu, Rgtotal
! endif
 end subroutine wtdcmrgau
!dCM-RgTotal
subroutine dCM_RgTotal(n,nstep)
  type (atom_dclist),pointer   :: la
   integer                     :: n,nstep,i,binCV1_1,binCV1_2,binCV2_1,binCV2_2   
   real(dp)                    :: CV1_ffcte,CV2_ffcte,CV1_ff(3),CV2_ff(3)
!Calculo mis variables colectivas en el paso n de la simulacion  
  call Collective_Variable


!Apertura de archivos de datos de la simulacion e inicializacion!!!!!!!!
  if(n == 0) call ini_2D(n)
  
  if(dCM >= potiniCV1 .AND. dCM <= potfinCV1) then
    if(Rgtotal >= potiniCV2 .AND. Rgtotal <= potfinCV2) then
    
!Calculo en que lugar de la grilla del potencial y las fuerzas estoy parado
  binCV1_1= int ((dCM-potiniCV1)/ real (dpot_CV1,dp)) + 1
  binCV1_2= binCV1_1 + 1
  binCV2_1= int ((Rgtotal-potiniCV2)/ real (dpot_CV2,dp)) + 1 
  binCV2_2= binCV2_1 + 1

  !Sumo las gauss dada una frecuencia
  if(mod (n,tauWT_pasos)==0)  call sumagauss2D(dCM,Rgtotal,binCV1_1,binCV1_2,binCV2_1,binCV2_2,n)
  
  CV1_ffcte= interpolacion_bilineal(dCM,Rgtotal,binCV1_1,binCV1_2,binCV2_1,binCV2_2,&
      F_WTMD_2D_CV1(binCV1_1,binCV2_1),F_WTMD_2D_CV1(binCV1_2,binCV2_1),&
      F_WTMD_2D_CV1(binCV1_1,binCV2_2),F_WTMD_2D_CV1(binCV1_2,binCV2_2))
  CV1_ff= (CV1_ffcte*(rcmA-rcmB))/ real (dCM,dp)
      
  CV2_ffcte = interpolacion_bilineal(dCM,Rgtotal,binCV1_1,binCV1_2,binCV2_1,binCV2_2,&
      F_WTMD_2D_CV2(binCV1_1,binCV2_1),F_WTMD_2D_CV2 (binCV1_2,binCV2_1),&
      F_WTMD_2D_CV2 (binCV1_1,binCV2_2),F_WTMD_2D_CV2 (binCV1_2,binCV2_2))
  CV2_ff = CV2_ffcte /real ((gmeta%nat*Rgtotal),dp)
!calculo de las fuerzas de los bias 
  la => gmeta%alist
  do i = 1,cant_parta
   la => la%next
   la%o%force = la%o%force + (CV1_ff*parA)
   la%o%force = la%o%force + (CV2_ff*(la%o%pos-rcm))   
  enddo
  do i=cant_parta+1, gmeta%nat
   la => la%next
   la%o%force= la%o%force - (CV1_ff*parB)
   la%o%force = la%o%force + (CV2_ff*(la%o%pos-rcm))   
  enddo                         
  endif
  endif
!Potencial externo sobre max dCM
  if (dCM > (potfinCV1-deltapared)) call pared_dCM_max((potfinCV1-deltapared))
!Potencial externo sobre min dcm
  if ((dCM < (potiniCV1+deltapared)) .AND. (potiniCV1 /= 0.0)) then
  call pared_dCM_min((potiniCV1+deltapared))
  endif  
  !Potencial externo sobre eg min  
  if (Rgtotal < potiniCV2+deltapared) call pared_RgT_min(potiniCV2+deltapared)
!Potencial externo sobre rg maxima
  if (Rgtotal > potfinCV2-deltapared) call pared_RgT_max(potfinCV2-deltapared)
!Burbuja sobre la posicion de las particulas
call pared_burb(n)
!Potencial externo sobre eg min  
  if (Rg < potiniaux1) call pared_RgCo_min(potiniaux1)
!Potencial externo sobre rg maxima
  if (Rg > potfinaux1) call pared_RgCo_max(potfinaux1)
!Potencial externo sobre eg min  
  if (RgAu < potiniaux2) call pared_RgShell_min(potiniaux2)
!Potencial externo sobre rg maxima
  if (RgAu > potfinaux2) call pared_RgShell_max(potfinaux2)


!!!!!!!!!!!!IMPRIMIR!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! if(mod (n,printD)==0)  then
!    write(236,'(f12.3,4(f12.6))') real (n*dt,dp), dCM, Rg, RgAu, Rgtotal
!  endif

 end subroutine dCM_RgTotal



  
  
!! !!!!!!!!!!! Set de parametros para correr wtmd 1D con o sin barreras!!!!!!!
subroutine wtmetad_set(wt1,wt2,wt3,wt4,wt5,wt6,wt7)
real(dp),intent(in)    :: wt1    ! Altura inicial de las gaussianas
real(dp),intent(in)    :: wt2    ! Ancho de las gaussianas
real(dp),intent(in)    :: wt3    ! Parametro de WTMD
real(dp),intent(in)    :: wt4    ! tau (frecuencia)
integer ,intent(in)    :: wt5   ! cada cuanto imprimir la CV
integer ,intent(in)    :: wt6   ! cantidad de particulas en el core
integer ,intent(in)    :: wt7   ! id de grupo de integracion donde se monta la Meta
         

!seleccion del grupo de integracion (cuales se mueven)
g => its%o(wt7)
call werr('Integration group should be in a constant T ensamble',.not.g%b_fixt)
temp_md=g%fixt
                 
wWTini      = wt1
sigCV1      = wt2
dT_WT       = wt3
tauWT       = wt4
printD      = wt5
cant_parta  = wt6


call wlog('WTMD1D_Pos', 'Altura inicial de la Gauss (KJ/mol)='//str(wWTini))
call wlog('WTMD1D_Pos', 'Ancho de medio pico de la Gauss (A)='//str(sigCV1))
call wlog('WTMD1D_Pos', 'Parametrop de la WTMetaD           ='//str(dT_WT))
call wlog('WTMD1D_Pos', 'Frecuencia de incorporacion (ps)   ='//str(tauWT))
call wlog('WTMD1D_Pos', 'Cantidad Particulas en el core     ='//str(cant_parta))
call wlog('WTMD1D_Pos', 'Cantidad Particulas en el core     ='//str(temp_md))
call wlog('WTMD1D_Pos', 'printCV                            ='//str(printD))
 
parA= 1.0_dp/real (cant_parta,dp)
parB= 1.0_dp/real ((gmeta%nat-cant_parta),dp)
                                        
wWTini  = wWTini * real(kjau_pro,dp)
tauWT_pasos= int(tauWT/dt)
call ini_1D
end subroutine wtmetad_set
      
  
subroutine wall1D_set(wt1,wt2,wt3,wt4)
  real(dp),intent(in)    :: wt1    ! potencial inicial CV1
  real(dp),intent(in)    :: wt2    ! potencial final CV1
  integer ,intent(in)    :: wt3    ! potencial bines CV1
  real(dp),intent(in)    :: wt4   ! rburb
  potiniCV1   = wt1
  potfinCV1   = wt2
  potbinCV1   = wt3
  rburb       = wt4


  dpot_CV1=(potfinCV1 - potiniCV1) / real(potbinCV1,dp)

  call wlog('W1D', 'potencial inicial (A) ='//str(potiniCV1))
  call wlog('W1D', 'potencial final (A)   ='//str(potfinCV1))
  call wlog('W1D', 'Cantidad de Bins      ='//str(potbinCV1))
  call wlog('W1D', 'Radio burb (A)        ='//str(rburb))
  
  potbinCV1=potbinCV1+1
  allocate(interpolation_x(potbinCV1))
  allocate(pot_int_y(potbinCV1))
  allocate(for_int_y(potbinCV1))
  allocate(dfor_int_y(potbinCV1))
end subroutine wall1D_set

!!!Well Tempered Metadynamics 1D !!!!!
!dCM
subroutine wtmetad(n,nstep)
  type (atom_dclist),pointer   :: la
  integer          :: i
  integer          :: binCV1_1,binCV1_2,n,nstep
  real(dp)         :: ff(3),ffcte
  
  call Collective_Variable
! Apertura  de archivos y poner a 0 las variables del potencial
  
  
  if(dCM >= potiniCV1 .AND. dCM <= potfinCV1) then
  
  binCV1_1= int ((dCM-potiniCV1) / real  (dpot_CV1,dp)) + 1
  binCV1_2= binCV1_1 + 1
  
! Cada tiempo tau incorporo una gaussiana a la matriz del bias y de las fuerzas
  if(mod (n,tauWT_pasos)==0)  call sumagauss1D(dCM,binCV1_1,binCV1_2,n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                      
! Coeficientes para el calculo de las fuerzas para todos los atomos
  ffcte= sline_interpotation(dCM,interpolation_x(binCV1_1),for_int_y(binCV1_1),&
      interpolation_x(binCV1_2),for_int_y(binCV1_2),dfor_int_y(binCV1_1),dfor_int_y(binCV1_2))
  ff= (ffcte*(rcmA-rcmB))/real (dCM,dp)
!calculo de la fuerza del bias en la cordenada colectiva dCM
  la => gmeta%alist
  do i = 1,cant_parta
   la => la%next
   la%o%force = la%o%force + real (ff*parA,dp)
  enddo
  do i=cant_parta+1, gmeta%nat
   la => la%next
   la%o%force= la%o%force - real (ff*parB,dp)
  enddo
  else
!Potencial externo sobre max dCM
  if (dCM > potfinCV1) call pared_dCM_max(potfinCV1)
!Potencial externo sobre min dcm
  if (dCM < potiniCV1) call pared_dCM_min(potiniCV1)
!Burbuja sobre la posicion de las particulas
endif
call pared_burb(n)
  if (Rg < potiniaux1) call pared_RgCo_min(potiniaux1)
!Potencial externo sobre rg maxima
  if (Rg > potfinaux1) call pared_RgCo_max(potfinaux1)
!Potencial externo sobre eg min  
  if (RgAu < potiniaux2) call pared_RgShell_min(potiniaux2)
!Potencial externo sobre rg maxima
  if (RgAu > potfinaux2) call pared_RgShell_max(potfinaux2)


!!!!!!!!!!!!!!!!IMPRIMIR!!!!!!!!!!!!!!!!!!!!!!!
!  if(mod (n,printD)==0)  then
!  write(236,'(f12.3,4(f12.6))') real (n*dt,dp), dCM, Rg, RgAu, Rgtotal
!  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if (n == nstep) then
!   open(UNIT=235,FILE= 'E_libre.dat')
!   fact=(-(temp_md+dT_WT)/(dT_WT))*0.01
!   do i=1, potbinCV1
!    write(235,*) interpolation_x(i),fact*pot_int_y(i)
!  enddo
!  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 end subroutine wtmetad


!Parametros de la DM siguiendo las CV
subroutine DM_CV_set(wt1,wt2,wt3)
integer ,intent(in)    :: wt1    ! cada cuanto imprimir 
real(dp),intent(in)    :: wt2    ! burbuja 
integer ,intent(in)    :: wt3    ! cantidad de atomos en el core 

printD     = wt1
rburb      = wt2
cant_parta = wt3

call wlog('WTMD1D_Pos','Parametros para Din. Mol. siguiendo CV')
call wlog('WTMD1D_Pos','Burbujas Potencial externo     ='//str(rburb))
call wlog('WTMD1D_Pos','Cantidad Particulas en el core ='//str(cant_parta))

end subroutine DM_CV_set

  
  !DM siguiendo CV
 subroutine DM_CV(n)
 integer n
  call Collective_Variable

  call pared_burb(n)
  
  !!!!!!!!!!!!!!IMPRIMIR!!!!!!!!!!!!!!!!!!!!!!!
 ! if(n == 0) then
 !  open(UNIT=222,FILE= 'CV.dat')
 !  write(222,*) '# tiempo     dCM      Rg_Co     Rg_Au     Rg_total'
 ! endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! if(mod (n,printD)==0)  then
 !   write(222,'(f12.3,4(f12.6))') n*dt, dCM, Rg, RgAu, Rgtotal
 ! endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine DM_CV
 

  
 function gauss_2D_Bias(x_CV1,x_CV2,CVs1,CVs2,wWT) result(resultado)
  implicit none
  real(dp) CVs1,CVs2,x_CV1,x_CV2,wWT
  real(dp) a_CV1s,b_CV2s
  real(dp)resultado
  a_CV1s = ((x_CV1-CVs1)*(x_CV1-CVs1))/ real ((2.0_dp*sigCV1*sigCV1),dp)
  b_CV2s = ((x_CV2-CVs2)*(x_CV2-CVs2))/ real ((2.0_dp*sigCV2*sigCV2),dp)
  resultado = wWT*exp(-(a_CV1s+b_CV2s)) 
 end function gauss_2D_Bias
 
  function gauss_2D_Force_CV1 (x_CV1,x_CV2,CVs1,CVs2,wWT) result(resultado)
  implicit none
  real(dp) CVs1,CVs2,x_CV1,x_CV2,wWT
  real(dp) a_CV1s,b_CV2s
  real(dp) resultado
  a_CV1s = ((x_CV1-CVs1)*(x_CV1-CVs1))/ real ((2.0_dp*sigCV1*sigCV1),dp)
  b_CV2s = ((x_CV2-CVs2)*(x_CV2-CVs2))/ real ((2.0_dp*sigCV2*sigCV2),dp)
  resultado = wWT*((x_CV1-CVs1)/real ((sigCV1*sigCV1),dp))*exp(-(a_CV1s+b_CV2s)) 
 end function gauss_2D_Force_CV1

 function gauss_2D_Force_CV2 (x_CV1,x_CV2,CVs1,CVs2,wWT) result(resultado)
  implicit none
  real(dp) CVs1,CVs2,x_CV1,x_CV2,wWT
  real(dp) a_CV1s,b_CV2s
  real(dp) resultado
  a_CV1s = ((x_CV1-CVs1)*(x_CV1-CVs1))/ real ((2.0_dp*sigCV1*sigCV1),dp)
  b_CV2s = ((x_CV2-CVs2)*(x_CV2-CVs2))/ real ((2.0_dp*sigCV2*sigCV2),dp)
  resultado = wWT*((x_CV2-CVs2)/real ((sigCV2*sigCV2),dp))*exp(-(a_CV1s+b_CV2s)) 
 end function gauss_2D_Force_CV2

 function interpolacion_bilineal(CVs1,CVs2,binCV1_1,binCV1_2,&
     binCV2_1,binCV2_2,Q_11,Q_21,Q_12,Q_22) result(resultado)
  implicit none
  integer binCV1_1,binCV1_2,binCV2_1,binCV2_2
  real(dp) CVs1,CVs2
  real Q_11,Q_21,Q_12,Q_22
  real(dp) x1_IB,x2_IB,y1_IB,y2_IB,fq11,fq21,fq12,fq22,term
  real (dp) resultado
  x1_IB= real (((binCV1_1*dpot_CV1)+potiniCV1),dp)
  x2_IB= real (((binCV1_2*dpot_CV1)+potiniCV1),dp)
  y1_IB= real (((binCV2_1*dpot_CV2)+potiniCV2),dp)
  y2_IB= real (((binCV2_2*dpot_CV2)+potiniCV2),dp)
  term=1.0_dp/ real (((x2_IB-x1_IB)*(y2_IB-y1_IB)),dp)
  fq11= Q_11*(x2_IB-CVs1)*(y2_IB-CVs2)
  fq21= Q_21*(CVs1-x1_IB)*(y2_IB-CVs2)
  fq12= Q_12*(x2_IB-CVs1)*(CVs2-y1_IB)
  fq22= Q_22*(CVs1-x1_IB)*(CVs2-y1_IB)
  resultado = term*(fq11+fq21+fq12+fq22)
 end function interpolacion_bilineal



 function sline_interpotation(xWT,x1WT,y1WT,x2WT,y2WT,dx1WT,dx2WT) result(resultado)
  implicit none
  real(dp) aWT,bWT,tWT,xWT,x1WT,x2WT,y1WT,y2WT,dx1WT,dx2WT,resultado
  aWT= dx1WT*(x2WT-x1WT)-(y2WT-y1WT)
  bWT= -dx2WT*(x2WT-x1WT)+(y2WT-y1WT)
  tWT= (xWT-x1WT)/(x2WT-x1WT) 
  resultado = ((1-tWT)*y1WT+(tWT*y2WT))+((tWT*(1-tWT))*((aWT*(1-tWT))+(bWT*tWT)))
 end function sline_interpotation

 function interpolacion_lineal(x_IL,x1_IL,x2_IL,y1_IL,y2_IL) result(resultado)
  real(dp) x_IL,x1_IL,x2_IL,y1_IL,y2_IL
  real(dp) a_IL,resultado
  a_IL= (y2_IL-y1_IL)/(x2_IL-x1_IL)
  resultado= (a_IL*(x_IL-x1_IL))+ y1_IL
 end function interpolacion_lineal

 function gauss1(i , CVs, wWT,sigCVs) result(resultado) 
  implicit none
  integer i
  real(dp) CVs,wWT,resultado,sigCVs
  resultado= wWT*exp((-(interpolation_x(i)-CVs)*(interpolation_x(i)-CVs))/(2.0*sigCVs*sigCVs))
 end function gauss1

 function gauss_forza1(i,CVs, wWT,sigCVs) result(resultado)
  implicit none
  integer i
  real(dp) CVs,wWT,resultado,sigCVs
  resultado=((wWT*(interpolation_x(i)-CVs))/(sigCVs*sigCVs))*&
     exp((- (interpolation_x(i)-CVs)*(interpolation_x(i)-CVs))/(2.0*sigCVs*sigCVs))
 end function gauss_forza1

 function gauss_dforza1 (i ,CVs,wWT,sigCVs) result(resultado)
  implicit none
  integer i
  real(dp) CVs, wWT, resultado,sigCVs
  resultado= (-wWT/(sigCVs*sigCVs))*(exp((-(interpolation_x(i)-CVs)*(interpolation_x(i)-CVs))/&
        (2.0*sigCVs*sigCVs)))*((((interpolation_x(i)-CVs)*(interpolation_x(i)-CVs))/(sigCVs*sigCVs))-1.0)
 end function gauss_dforza1


subroutine Collective_Variable()
  type (atom_dclist),pointer   :: la
  integer                      :: i
  real(dp)                     :: mas1,mas2,sumaRgAu,sumaRg,sumaRgtotal

  !inicializaciones
  dCM=0.0_dp
  rcmA=0.0_dp
  rcmB=0.0_dp
  RgAu=0.0_dp
  sumaRgAu=0.0_dp   
  Rg=0.0_dp
  sumaRg=0.0_dp   
  Rgtotal=0.0_dp
  sumaRgtotal=0.0_dp   

 !   Calculo la posicion del CM  de cada especie
 la => gmeta%alist
  do i=1, cant_parta
   la => la%next
   rcmA = rcmA + (la%o%pos*la%o%mass)
   rcm = rcm + (la%o%pos*la%o%mass)
   if(i/=1) cycle
   mas1= la%o%mass
  enddo
  do i=cant_parta+1, gmeta%nat
   la => la%next
   rcmB = rcmB + (la%o%pos*la%o%mass)
   rcm = rcm + (la%o%pos*la%o%mass)
   if(i/=(cant_parta+1)) cycle
   mas2= la%o%mass
  enddo
  rcmA = rcmA/ real ((cant_parta*mas1),dp)
  rcmB = rcmB/ real (((gmeta%nat-cant_parta)*mas2),dp)
  rcm = rcm/ real (((gmeta%nat-cant_parta)*mas2)+(cant_parta*mas1),dp)
  
  la => gmeta%alist
  do i=1, cant_parta
   la => la%next
   sumaRg= sumaRg + (dot_product(la%o%pos-rcmA, la%o%pos-rcmA))
   sumaRgtotal= sumaRgtotal + (dot_product(la%o%pos-rcm, la%o%pos-rcm))
  enddo
  do i=cant_parta + 1, gmeta%nat 
   la => la%next
   sumaRgAu= sumaRgAu + (dot_product(la%o%pos-rcmB, la%o%pos-rcmB))
   sumaRgtotal= sumaRgtotal + (dot_product(la%o%pos-rcm, la%o%pos-rcm))
  enddo

  Rgtotal=sqrt(sumaRgtotal / real (gmeta%nat,dp))
  RgAu=sqrt(sumaRgAu / real (gmeta%nat - cant_parta,dp))
  Rg=sqrt(sumaRg / real (cant_parta,dp))
  dCM=sqrt(((rcmB(1)-(rcmA(1)))*(rcmB(1)-rcmA(1)))+((rcmB(2)-rcmA(2))*(rcmB(2)-rcmA(2)))+((rcmB(3)-rcmA(3))*(rcmB(3)-rcmA(3))))

end subroutine Collective_Variable

 function bias_point_1D(CV) result(resultado) 
  implicit none
  real(dp) CV,resultado
  integer binCV1_1,binCV1_2
  binCV1_1= int ((CV-potiniCV1) / real  (dpot_CV1,dp)) + 1
  binCV1_2= binCV1_1 + 1
   resultado= sline_interpotation(CV,interpolation_x(binCV1_1),pot_int_y(binCV1_1),interpolation_x(binCV1_2),&
       pot_int_y(binCV1_2),-for_int_y(binCV1_1),-for_int_y(binCV1_2))
 end function bias_point_1D



 subroutine pared_dCM_max(pared)
  implicit none
   real(dp)                     :: ffCM(3),pared
   type (atom_dclist),pointer   :: la
   integer                      :: i
   ffCM = -kBuCM*(dCM-pared)*((rcmA-rcmB)/dCM)
   la => gmeta%alist
   do i = 1,cant_parta
    la => la%next
    la%o%force = la%o%force + (ffCM*parA)
   enddo
   do i=cant_parta+1, gmeta%nat
    la => la%next
    la%o%force = la%o%force - (ffCM*parB)
   enddo
  end subroutine pared_dCM_max
   
 subroutine pared_dCM_min(pared)
  implicit none
   real(dp)                     :: ffCM(3),pared
   type (atom_dclist),pointer   :: la
   integer                      :: i
   ffCM = -kBuCM*(dCM-pared)*((rcmA-rcmB)/dCM)

   la => gmeta%alist
   do i = 1,cant_parta
    la => la%next
    la%o%force = la%o%force + (ffCM*parA)
   enddo
   do i=cant_parta+1, gmeta%nat
    la => la%next
    la%o%force = la%o%force - (ffCM*parB)
   enddo
 end subroutine pared_dCM_min 
 
 subroutine pared_RgCo_min(pared)
  implicit none
  real(dp)                     :: ffRgT,pared
  integer                      :: i
  type (atom_dclist),pointer   :: la
   ffRgT = (-kBuRg*(Rg-pared))/real ((Rg*cant_parta),dp)
   la => gmeta%alist
   do i = 1,cant_parta
    la => la%next 
    la%o%force = la%o%force + (ffRgT*(la%o%pos-rcmA))   
   enddo
   end subroutine pared_RgCo_min
   
 subroutine pared_RgShell_min(pared)
  implicit none
  real(dp)                     :: ffRgT,pared
  integer                      :: i
  type (atom_dclist),pointer   :: la
   ffRgT = (-kBuRg*(RgAu-pared))/real ((RgAu*(gmeta%nat-cant_parta)),dp)
   la => gmeta%alist
   do i = cant_parta+1,gmeta%nat
    la => la%next 
    la%o%force = la%o%force + (ffRgT*(la%o%pos-rcmB))   
   enddo
   end subroutine pared_RgShell_min
 
   subroutine pared_RgCo_max(pared)
  implicit none
  real(dp)                     ::ffRgT,pared
  integer                      :: i
  type (atom_dclist),pointer   :: la
   ffRgT = (-kBuRg*(Rg-pared))/real ((Rg*cant_parta),dp)
   la => gmeta%alist
   do i = 1,cant_parta
    la => la%next 
    la%o%force = la%o%force + (ffRgT*(la%o%pos-rcmA))   
   enddo  
  end subroutine pared_RgCo_max
 
  subroutine pared_RgShell_max(pared)
  implicit none
  real(dp)                     ::ffRgT,pared
  integer                      :: i
  type (atom_dclist),pointer   :: la
   ffRgT = (-kBuRg*(RgAu-pared))/real ((RgAu*(gmeta%nat-cant_parta)),dp)
   la => gmeta%alist
   do i = cant_parta+1,gmeta%nat
    la => la%next 
    la%o%force = la%o%force + (ffRgT*(la%o%pos-rcmB))   
   enddo  
  end subroutine pared_RgShell_max
   
 subroutine pared_RgT_min(pared)
  implicit none
  real(dp)                     :: ffRgT,pared
  integer                      :: i
  type (atom_dclist),pointer   :: la
   ffRgT = (-kBuRg*(Rgtotal-pared))/real ((Rgtotal*gmeta%nat),dp)
   la => gmeta%alist
   do i = 1,cant_parta
    la => la%next 
    la%o%force = la%o%force + (ffRgT*(la%o%pos-rcm))   
   enddo
   do i = cant_parta +1,gmeta%nat
    la => la%next 
    la%o%force = la%o%force + (ffRgT*(la%o%pos-rcm))   
   enddo  
   end subroutine pared_RgT_min
   
 subroutine pared_RgT_max(pared)
  implicit none
  real(dp)                     ::ffRgT,pared
  integer                      :: i
  type (atom_dclist),pointer   :: la
   ffRgT = (-kBuRg*(Rgtotal-pared))/real ((Rgtotal*gmeta%nat),dp)
   la => gmeta%alist
   do i = 1,cant_parta
    la => la%next 
    la%o%force = la%o%force + (ffRgT*(la%o%pos-rcm))   
   enddo  
   do i = cant_parta +1,gmeta%nat
    la => la%next 
    la%o%force = la%o%force + (ffRgT*(la%o%pos-rcm))   
   enddo  
  end subroutine pared_RgT_max

  subroutine pared_burb(n)
  implicit none
  integer                      :: i,n
  real(dp)                     :: dist_part2,dist_part,ffburb
  type (atom_dclist),pointer   :: la
  la => gmeta%alist

  do i = 1,gmeta%nat
   la => la%next
   dist_part2 = dot_product(la%o%pos,la%o%pos)
   if (dist_part2 < (rburb*rburb)) cycle
    dist_part = sqrt(dist_part2)
    ffburb = (-kBu*(dist_part-rburb))/dist_part
    la%o%force= la%o%force + (ffburb * la%o%pos)
  enddo
  end subroutine pared_burb




  subroutine sumagauss2D(CV1,CV2,binCV1_1,binCV1_2,binCV2_1,binCV2_2,n)
   integer                      :: i,binCV1_1,binCV1_2,binCV2_1,binCV2_2,n,j
   integer                      :: bin_ini_CV1,bin_fin_CV1
   integer                      :: bin_ini_CV2,bin_fin_CV2
   real(dp)                     :: CV1,CV2,bias_ant,wWT,potential_point
   real(dp)                     :: force_point_cv1,force_point_cv2
   real(dp)                     :: x_ini,x_fin,y_ini,y_fin,bincte_CV1,bincte_CV2 

   bias_ant = interpolacion_bilineal(CV1,CV2,binCV1_1,binCV1_2,binCV2_1,binCV2_2,&
       B_WTMD_2D(binCV1_1,binCV2_1),B_WTMD_2D(binCV1_2,binCV2_1),&
       B_WTMD_2D(binCV1_1,binCV2_2),B_WTMD_2D(binCV1_2,binCV2_2))
   wWT=wWTini*exp(-bias_ant/(kB_ui*dT_WT))
   !Grabo cada una de las gaussianas que deposito para luego armar el perfil de energia libre
   if(n /= 0) write(unidad,'(f12.3,3(f12.6))')  real(n*dt,dp), CV1, CV2, wWT*biasf*0.01
   !limitar desde que punto y hasta que punto sumo las gaussianas en la direccion de la distancia al centro de masa
   bincte_CV1=sqrt (-2.0_dp*sigCV1*sigCV1*log(tol_Bias/ real (wWT)))
    x_ini= CV1 - bincte_CV1 
    if(x_ini<potiniCV1) x_ini=potiniCV1 
    x_fin= CV1 + bincte_CV1
    if(x_fin>potfinCV1) x_fin= potfinCV1
    bin_ini_CV1= int ((x_ini-potiniCV1)/ real (dpot_CV1,dp)) + 1
    bin_fin_CV1= int ((x_fin-potiniCV1)/ real (dpot_CV1,dp)) + 1
   !limitar desde que punto y hasta que punto sumo las gaussianas en la direccion del radio de giro
    bincte_CV2=sqrt (-2.0_dp*sigCV2*sigCV2*log(tol_Bias/ real (wWT)))
    y_ini= CV2 - bincte_CV2
    if(y_ini<potiniCV2) y_ini= potiniCV2  
    y_fin= CV2 + bincte_CV2
    if(y_fin>potfinCV2) y_fin=potfinCV2
    bin_ini_CV2= int ((y_ini-potiniCV2)/ real (dpot_CV2,dp)) + 1
    bin_fin_CV2= int ((y_fin-potiniCV2)/ real (dpot_CV2,dp)) + 1
   do i=bin_ini_CV1,bin_fin_CV1
    do j=bin_ini_CV2,bin_fin_CV2
     potential_point= gauss_2D_Bias (((real (i,dp)*dpot_CV1)+potiniCV1),&
         ((real (j,dp)*dpot_CV2)+potiniCV2),CV1,CV2,wWT) 
     B_WTMD_2D(i,j) = B_WTMD_2D(i,j) + potential_point 
     force_point_cv1 = gauss_2D_Force_CV1 (((real (i,dp)*dpot_CV1)+potiniCV1),&
         ((real (j,dp)*dpot_CV2)+potiniCV2),CV1,CV2,wWT)
     F_WTMD_2D_CV1(i,j) = F_WTMD_2D_CV1(i,j) + force_point_cv1
     force_point_cv2 = gauss_2D_Force_CV2 (((real (i,dp)*dpot_CV1)+potiniCV1),&
         ((real (j,dp)*dpot_CV2)+potiniCV2),CV1,CV2,wWT)
     F_WTMD_2D_CV2(i,j)  = F_WTMD_2D_CV2(i,j) + force_point_cv2
    enddo
   enddo
 end subroutine sumagauss2D



 subroutine sumagauss1D(CV1,binCV1_1,binCV1_2,n)
 integer                       :: i,binCV1_1,binCV1_2,n
 real(dp)                      :: bias_ant,CV1,wWT,biasf
   bias_ant= sline_interpotation(CV1,interpolation_x(binCV1_1),pot_int_y(binCV1_1),interpolation_x(binCV1_2),&
       pot_int_y(binCV1_2),-for_int_y(binCV1_1),-for_int_y(binCV1_2))
   wWT=wWTini*exp(-bias_ant/(kB_ui*dT_WT))
   biasf= real (((temp_md+dT_WT)/(dT_WT)),dp)
!   if(n /= 0) write(235,'(f12.3,3(f12.6))')  real(n*dt,dp), CV1, sigCV1,wWT*biasf*0.01

   do i=1,potbinCV1
    pot_int_y(i) = pot_int_y(i)  + gauss1(i,CV1,wWT,sigCV1)
    for_int_y(i) = for_int_y(i)  + gauss_forza1(i,CV1,wWT,sigCV1)
    dfor_int_y(i)= dfor_int_y(i) + gauss_dforza1 (i,CV1,wWT,sigCV1)
   enddo
  end subroutine sumagauss1D


 subroutine ini_2D(n)
   use gems_input_parsing, only: ioprefix
   use gems_constants, only:find_io
   integer         :: n,j,i!,ierr,bincv1_1,bincv1_2,bincv2_1,bincv2_2
   ! real(dp)        :: c_1,c_2,ww,w_old,t,s_1,s_2,fac 
   
   unidad=find_io(10)
   open(unidad,FILE= 'Par_E_libre.'//trim(ioprefix)//'.dat')
   write(unidad,*) '#!Sigma_CV1',sigCV1
   write(unidad,*) '#!Sigma_CV2',sigCV2
   write(unidad,*) '#!Tau', tauWT 
   write(unidad,*) '#!Temp', temp_md
   write(unidad,*) '#!Tol', tol_Bias
   if(potiniCV1==0.0) then
   write(unidad,*) '#!CV_1_min', potiniCV1
   else
   write(unidad,*) '#!CV_1_min', potiniCV1+deltapared
   endif
   write(unidad,*) '#!CV_1_max', potfinCV1-deltapared
   write(unidad,*) '#!CV_2_min', potiniCV2+deltapared
   write(unidad,*) '#!CV_2_max', potfinCV2-deltapared
!   write(236,* ) '# tiempo    dCM    RgCo    RgAu   Rgtotal'
   bias_ant=0.0_dp
   do i=1,potbinCV1
     do j=1,potbinCV2
       B_WTMD_2D(i,j)=0.0_dp
       F_WTMD_2D_CV1(i,j)=0.0_dp
       F_WTMD_2D_CV2(i,j)=0.0_dp
     enddo
   enddo   
   !cargando si tenemos un archivo de gaussianas
!   k=0
!   open(UNIT=658,FILE= 'Gauss_old.dat',status='old',err=999)
!     read(658,*)      
!     read(658,*)      
!   do
!     read(658,*,IOSTAT=iErr) t, C_1, C_2, s_1,s_2,w_old, fac             
!     IF (iErr.NE.0) EXIT 
!     if(k==0) ww=w_old
!     k=1
     !Calculo en que lugar de la grilla del potencial y las fuerzas estoy parado
!       binCV1_1= int ((C_1-potiniCV1)/ real (dpot_CV1,dp)) + 1
!       binCV1_2= binCV1_1 + 1
!       binCV2_1=  int ((C_2-potiniCV2)/ real (dpot_CV2,dp)) + 1 
!       binCV2_2= binCV2_1 + 1
!     !Control si me fui afuera de la grilla
!       call sumagauss2D(C_1,C_2,binCV1_1,binCV1_2,binCV2_1,binCV2_2,n)
!   end do
!999 if(k==0) then
!      write(*,*) 'Metadynamics new'
!      else
!      write(*,*) 'Metadynamics restart'
!      if(s_1==sigCV1) write(*,*) 'sigma CV1 iguales'
!      if(s_2==sigCV2) write(*,*) 'sigma CV2 iguales'
!      if(ww==(wWTini*0.01*((temp_md+dT_WT)/dT_WT))) write(*,*) 'w-dT-T iguales'
!      endif
      end subroutine ini_2D

      subroutine ini_1D()
        integer         :: i!,ierr,bincv1_1,bincv1_2
       ! real(dp)        :: c_1,ww,w_old,t,s_1
!       open(UNIT=235,FILE= 'Gauss.dat')
!       write(235,*) '#! FIELDS time dCM Rg sigma_dCM  sigma_Rg height biasf'
!       write(235,*) '#! SET multivariate false' 
!       open(UNIT=236,FILE= 'VC.dat')
!       write(236,* ) '# tiempo    dCM     Rg_Co    Rg_Au   Rg_total'

       do i=1,potbinCV1
         interpolation_x(i)=(dpot_CV1*(real (i,dp)))+potiniCV1
         pot_int_y(i) =0.0_dp
         for_int_y(i) =0.0_dp
         dfor_int_y(i)=0.0_dp
       enddo
!       k=0
!       open(UNIT=658,FILE= 'Gauss_old.dat',status='old',err=999)
!       read(658,*)      
!       read(658,*)      
!   do
!       read(658,*,IOSTAT=iErr) t, C_1, s_1,w_old
!       IF (iErr.NE.0) EXIT 
!       if(k==0) ww=w_old
!       k=1
!       !Calculo en que lugar de la grilla del potencial y las fuerzas estoy parado
!       binCV1_1= int ((C_1 - potiniCV1) / real  (dpot_CV1,dp)) + 1
!       binCV1_2= binCV1_1 + 1
!       call sumagauss1D(C_1,binCV1_1,binCV1_2,n)
!   end do
!999 if(k==0) then
!      write(*,*) 'Metadynamics new'
!      else
!      write(*,*) 'Metadynamics restart'
!      if(s_1==sigCV1) write(*,*) 'sigma CV1 iguales'
!      if(ww==(wWTini*0.01*((temp_md+dT_WT)/dT_WT))) write(*,*) 'w-dT-T iguales'
!      endif



      end subroutine ini_1D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Lucas !!!!!!!!!!
subroutine write_cvs(op)
  class(outpropa)     :: op
 ! call Collective_variable                    
  op%f(1) = op%f(1) + dCM
  op%f(2) = op%f(2) + Rg
  op%f(3) = op%f(3) + RgAu
  op%f(4) = op%f(4) + Rgtotal
  op%f(5) = op%f(5) + rcmA(1)
  op%f(6) = op%f(6) + rcmA(2)
end subroutine
!!!!!!!!!!!!!!!!!!!!!! 

subroutine write_E_1D(of)
  class(outfile)     :: of
  integer            :: i
  real(dp)           :: fact 

  fact=(-(temp_md+dT_WT)/(dT_WT))*0.01
   do i=1, potbinCV1
   write(of%un,'(e25.12,x,e25.12)') interpolation_x(i),fact*pot_int_y(i)
   enddo
   
   write(of%un,'(a)')

  if(of%flush) call flush(of%un)

end subroutine
  
subroutine wtx(n,nstep)
  type (atom_dclist),pointer   :: la
  integer                      :: i
  integer binCV1_1,binCV1_2,n,nstep
  real(dp)  :: ffcte
  call Collective_Variable
! Apertura  de archivos y poner a 0 las variables del potencial
  
  
  if(rcmA(1) >= potiniCV1 .AND. rcmA(1) <= potfinCV1) then
    binCV1_1= int ((rcmA(1)-potiniCV1) / real  (dpot_CV1,dp)) + 1
    binCV1_2= binCV1_1 + 1
  
! Cada tiempo tau incorporo una gaussiana a la matriz del bias y de las fuerzas
    if(mod (n,tauWT_pasos)==0)  call sumagauss1D(rcmA(1),binCV1_1,binCV1_2,n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                      
! Coeficientes para el calculo de las fuerzas para todos los atomos
    ffcte= sline_interpotation(rcmA(1),interpolation_x(binCV1_1),for_int_y(binCV1_1),&
      interpolation_x(binCV1_2),for_int_y(binCV1_2),dfor_int_y(binCV1_1),dfor_int_y(binCV1_2))
!calculo de la fuerza del bias en la cordenada colectiva dCM
    la => gmeta%alist
    do i = 1,cant_parta
      la => la%next
      la%o%force(1) = la%o%force(1) + real (ffcte*parA,dp)
    enddo
  else
!Potencial externo sobre max dCM
    if (rcmA(1) > potfinCV1) call pared_dCM_max(potfinCV1)
!Potencial externo sobre min dcm
    if (rcmA(1) < potiniCV1) call pared_dCM_min(potiniCV1)
!Burbuja sobre la posicion de las particulas
  endif
  call pared_burb(n)
  if (rcmA(2) < potiniaux1) call pared_RgCo_min(potiniaux1)
!Potencial externo sobre rg maxima
  if (rcmA(2) > potfinaux1) call pared_RgCo_max(potfinaux1)
 end subroutine wtx





subroutine wtxy_2D(n,nstep)
  type (atom_dclist),pointer   :: la
  integer                     :: n,nstep,i,binCV1_1,binCV1_2,binCV2_1,binCV2_2   
  real(dp)                    :: CV1_ffcte,CV2_ffcte

!Calculo mis variables colectivas en el paso n de la simulacion  
  call Collective_Variable
!Apertura de archivos de datos de la simulación e inicializacion!!
  if(n == 0) call ini_2D(n)
! write(*,*) rcmA(1), rcmA(2)
  if(rcmA(1) >= potiniCV1 .AND. rcmA(1) <= potfinCV1) then
    if(rcmA(2) >= potiniCV2 .AND. rcmA(2) <= potfinCV2) then
      !Calculo en que lugar de la grilla del potencial y las fuerzas estoy parado
      binCV1_1= int ((rcmA(1)-potiniCV1)/ real (dpot_CV1,dp)) + 1
      binCV1_2= binCV1_1 + 1
      binCV2_1= int ((rcmA(2)-potiniCV2)/ real (dpot_CV2,dp)) + 1 
      binCV2_2= binCV2_1 + 1
      if(mod (n,tauWT_pasos)==0)  call sumagauss2D(rcmA(1),rcmA(2),binCV1_1,binCV1_2,binCV2_1,binCV2_2,n)
      CV1_ffcte= interpolacion_bilineal(rcmA(1),rcmA(2),binCV1_1,binCV1_2,binCV2_1,binCV2_2,&
        F_WTMD_2D_CV1(binCV1_1,binCV2_1),F_WTMD_2D_CV1(binCV1_2,binCV2_1),&
        F_WTMD_2D_CV1(binCV1_1,binCV2_2),F_WTMD_2D_CV1(binCV1_2,binCV2_2))
      CV2_ffcte = interpolacion_bilineal(rcmA(1),rcmA(2),binCV1_1,binCV1_2,binCV2_1,binCV2_2,&
        F_WTMD_2D_CV2(binCV1_1,binCV2_1),F_WTMD_2D_CV2 (binCV1_2,binCV2_1),&
        F_WTMD_2D_CV2 (binCV1_1,binCV2_2),F_WTMD_2D_CV2 (binCV1_2,binCV2_2))
   
      la => gmeta%alist
      do i = 1,cant_parta
        la => la%next
        la%o%force(1) = la%o%force(1) + real (CV1_ffcte*parA,dp)
        la%o%force(2) = la%o%force(2) + real (CV2_ffcte*parA,dp)
      enddo
    endif
  endif
  !Potencial externo sobre max dCM
  if (rcmA(1) > (potfinCV1-deltapared)) call pared_dCM_max((potfinCV1-deltapared))
!Potencial externo sobre min dcm
  if ((rcmA(1) < (potiniCV1+deltapared)) .AND.(potiniCV1 /= 0.0)) call &
  pared_dCM_min((potiniCV1+deltapared))
  !Potencial externo sobre eg min  
  if (rcmA(2) < potiniCV2+deltapared) call pared_RgCo_min(potiniCV2+deltapared)
!Potencial externo sobre rg maxima
  if (rcmA(2) > potfinCV2-deltapared) call pared_RgCo_max(potfinCV2-deltapared)
!Burbuja sobre la posicion de las particulas
  call pared_burb(n)
 end subroutine wtxy_2D

end module gems_metadynamics

