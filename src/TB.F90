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



!SETS DE PARAMETROS HOMOATOMICOS: 
!   (  Co-Co; Cu-Cu;  Pd-Pd; Ag-Ag; Pt-Pt; Au-Au   ) 

!Para Interacciones HOMO, se tiene que :

!r0=Parametro de red / dsqrt(2)  Primeros Vecinos
!rcut_1=r0*dsqrt(2)              Segundos Vecinos (Donde se da el Empalme con el polinomio soft-cutt)
!rcut_2=r0*dsqrt(3)              Terceros Vecinos (Donde el polinomio soft-cutt se hace cero)

!SETS DE PARAMETROS HETEROATOMICOS:
!   (  Co-Ag; Ni-Ag; Cu-Ag; Pd-Ag;  )
!   (  Co-Au; Pt-Au                 )
 
!Para interacciones HETERO, se tiene que :

! r0=[r0(A-A)+r0(B-B)] / 2
! rcut_1= Parametro de red (A)   
! rcut_2= Parametro de red (B)
! Si se cumple que : Parametro de red (A) < Parametro de red (B)


module gems_tb
 use gems_program_types
 use gems_constants
 use gems_constants
 use gems_algebra 
 use gems_tables
 use gems_neighbour
!$ use omp_lib 
 implicit none

 private
 public tb_preinteraction,tb_interaction,tb_set, tb_parameter_set

 !Lo hago publico, para que lo vea el modulo donde esta el Intergroup puesto que lo necesita.
 public rcut_ext

 ! NOTA: Si se desea declarar un tipo derivado del integroup que contenga los
 ! parametros de la interaccion, se tiene el problema de que la subrrutina
 ! interact tiene que tener como passed argument la clase intergroup, forzando a
 ! usar select type dentro de ella. Se podría declarar al intergroup como tipo
 ! abstract, permitiendo hacer que interact sea un deferred procedure, pero en
 ! este caso todas las interacciones deben tener un tipo derivados, ya que los
 ! abstract types no pueden ser instanciados.
 ! ver http://stackoverflow.com/a/25413107/1342186
 ! El problema de esto es que cada procedimiento de interact fuerza a definir un
 ! nuevo tipo derivado ya que no se puede reapuntar el interact puesto que no es
 ! puntero


 ! Lista privada de interacciones tb.
 type(intergroup_dl),target  :: igrtb_dl

 type :: tb_parameter

    ! Parametros de TB
    real(dp)    ::  a   = 0._dp,&
                    eps = 0._dp,&
                    p   = 0._dp,&
                    q   = 0._dp,&
                    r0  = 0._dp

    ! Regla de suavizado               
    real(dp)    ::  rci = 0._dp,&
                    rce = 0._dp,&
                    b1  = 0._dp,&
                    b2  = 0._dp,&
                    b3  = 0._dp,&
                    r1  = 0._dp,&
                    r2  = 0._dp,&
                    r3  = 0._dp

    ! logical     ::  set_b = .false.

 end type tb_parameter

 type(tb_parameter)   :: tbp(118,118)

 real(dp),allocatable :: band(:),eband(:)

 !Este es el radio de corte mayor utilizado para el TB.
 !TODO: Que si se seleciona otro juego de parametros que reduce el radio de
 !corte esto se tenga en cuenta
 real(dp)   :: rcut_ext=0.0_dp    
       

contains

subroutine tb_parameter_set(z1,z2,a,eps,p,q,r0,rci,rce)
integer,intent(in) :: z1,z2
real(dp),intent(in) :: a,eps,p,q,r0,rce,rci
real(dp)            :: dcut,ab,bb,cb,aux

tbp(z1,z2)%a   = a  
tbp(z1,z2)%eps = eps
tbp(z1,z2)%p   = p  
tbp(z1,z2)%q   = q  
tbp(z1,z2)%r0  = r0 

!Polinomio de suavizado. 
!Logsdail etal. RSC Advances 2(2012)5863. 10.1039/c2ra20309j (en el supporting esta el polinomio)
!Logsdail dice que este polinomio machea el suavisado de Baletto etal. (2002). JCP, 116(9), 3856–3863. 10.1063/1.1448484. En este
!paper Balleto menciona que usa entre segundos y terceros vecinos.
!En Rapallo etal. JCP 122(2005)194308. 10.1063/1.1898223 (en este ultimo dice que rext es el 3er vecino)
tbp(z1,z2)%rci = rci
tbp(z1,z2)%rce = rce

! Aleaciones vieja tirada: Como estaba antes
! tbp(j,j)%rci = 4.0729350596D0 !rci del Au
! tbp(j,j)%rce = 4.4271218640D0

!TODO: Que si se seleciona otro juego de parametros que reduce el radio de
!corte esto se tenga en cuenta
rcut_ext=max(rcut_ext,tbp(z1,z2)%rce)

! Polinomial coefficient for band energy (ver supporting info of Logsdail et al.).
! Notese aqui que aux es la raiz de la energía de banda (ver supporting info of Logsdail et al.).
dcut=tbp(z1,z2)%rce-tbp(z1,z2)%rci
ab=-   1.0_dp/dcut**3
bb=-(q/r0)   /dcut**2
cb=-(q/r0)**2/dcut
aux=eps*dexp(-q*(tbp(z1,z2)%rci/r0-1.0_dp))
tbp(z1,z2)%b1=aux*(12.0_dp*ab-6.0_dp*bb+cb)/(2.0_dp*dcut**2)
tbp(z1,z2)%b2=aux*(15.0_dp*ab-7.0_dp*bb+cb)/dcut
tbp(z1,z2)%b3=aux*(20.0_dp*ab-8.0_dp*bb+cb)/2.0_dp

! Polynomia coeffciens for repulstion energy (ver supporting info of Logsdail et al.).
aux=a*dexp(-p*(tbp(z1,z2)%rci/r0-1.0_dp))
bb=-(p/r0)   /dcut**2
cb=-(p/r0)**2/dcut
tbp(z1,z2)%r1=aux*(12.0_dp*ab-6.0_dp*bb+cb)/(2.0_dp*dcut**2)
tbp(z1,z2)%r2=aux*(15.0_dp*ab-7.0_dp*bb+cb)/dcut
tbp(z1,z2)%r3=aux*(20.0_dp*ab-8.0_dp*bb+cb)/2.0_dp
  
! In case of an aleation
tbp(z2,z1)=tbp(z1,z2)

end subroutine

subroutine tb_set(prmfile,g1,g2)
use gems_input_parsing
use gems_inq_properties, only: inq_pure
character(*),intent(in)     :: prmfile
type(group),intent(in)    :: g1
type(group),intent(in),optional    :: g2
type(intergroup),pointer    :: igr

! Initialize the integroup
allocate(igr)
if(present(g2)) then
  call igr%init(g1=g1,g2=g2)
else
  call igr%init(g1=g1)
endif

! Add this interaction to the interal TB interaction list
call igrtb_dl%add_after()
call igrtb_dl%next%point(igr)

! Set interaction  
igr%interact => tb_interaction

! Alloc the band, eband
allocate(eband(natoms))
allocate(band(natoms))

! if(allocated(band)) then
!   call wwan('Note that TB energy of a system can not be expresed as the sum of two subsystems TB energies.'  )
!   return
! endif

! Parameters default
call tb_readprm(prmfile) 

! Radio de corte
call igr%setrc(rcut_ext)
       
end subroutine tb_set
          
subroutine tb_readprm(prmfile)
! Subrroutine to read the parameter file
use gems_input_parsing
use gems_elements, only:inq_z
use gems_errors, only:wlog
use gems_algebra, only:sort_int
character(*),intent(in)     :: prmfile
type(input_options), target :: iopts
character(90)               :: clase,w1,w2
real(dp)                    :: a,eps,p,q,r0,rci,rce
integer                     :: ix(2),u

! Redirigiendo el input parsing a las opciones locales
u = find_io(30)
open(u,action='read',file=prmfile)
call iopts%init(skip_blank_lines=.true.,in=u,echo_lines=.false.,comments='!')
opts=>iopts

do
  call read_line(eof)
  if (eof) exit

  call reada(w1)
  if (w1(1:1)=='!') cycle

  ! Selecciono la clase
  selectcase(w1)
  case('GUPTA','SMATB')
    clase=w1
    cycle
  endselect

  selectcase(clase)
  case('GUPTA','SMATB')

    !a: 
    !eps:
    !p:
    !q:
    !r0: A
    !rci: A
    !rce: A

    call reada(w2) !nomb2
    call readf(a)
    call readf(eps)
    call readf(p)
    call readf(q)
    call readf(r0)
    call readf(rci)
    call readf(rce)

    ! Solving symmetry in bonds
    ix(1)=inq_z(w1)
    ix(2)=inq_z(w2)
    call sort_int(ix(1:2))

    ! Add parameters to the corresponding bondgr
    call tb_parameter_set(ix(1),ix(2),a,eps,p,q,r0,rci,rce)
  endselect
enddo
 
! Devuelvo el input parsing al gems
opts=>gems_iopts
close(u) 
 
end subroutine tb_readprm
 
subroutine tb_preinteraction
  real(dp)                  :: dr,vd(dm),banda,r0
  real(dp)                  :: drm,drm2,drm3,drm4,drm5
  integer                   :: z1,z2
  integer                   :: i,j,l
  type(intergroup),pointer    :: ig
  type(atom),pointer          :: o1,o2
  type(intergroup_dl),pointer :: node


  if(.not.allocated(band)) return

  ! Redimensiono las variables si es necesario
  if(natoms>size(band)) then
    deallocate(band )
    deallocate(eband)
    allocate(band(natoms))
    allocate(eband(natoms))
  endif
          
  ! Set zeros
  band(:)=0.0_dp
  eband(:)=0.0_dp

  ! Interaccion
  node => igrtb_dl
  do while( associated(node%next) )
    node => node%next
    ig => node%o

    !$OMP PARALLEL DO DEFAULT(NONE) &
    !$OMP& FIRSTPRIVATE(ig)    &
    !$OMP& PRIVATE(l,z1,j,z2,dr,vd,banda,drm,drm2,drm3,drm4,drm5,r0,o1,o2) &
    !$OMP& SHARED(tbp,a,band,mic,nlocal) 

    ! Sobre los atomos
    do i = 1,ig%n(1)

      o1 => ig%at(i)%o
      z1 = o1%z
      
      ! sobre los vecinos
      do l = 1, ig%nn(i)

        j  = ig%list(i,l)
        o2 => ig%at(j)%o
        z2 = o2%z

        vd = vdistance( o2, o1, mic)
        dr =  dsqrt(dot_product(vd,vd)) 
 
        if(dr>tbp(z1,z2)%rce) cycle

        if(dr<tbp(z1,z2)%rci)then

          r0=tbp(z1,z2)%r0 
        
          banda=tbp(z1,z2)%eps*tbp(z1,z2)%eps*exp(-2.0_dp*tbp(z1,z2)%q*(dr/r0-1.0_dp))
      
        else
        
          drm=dr-tbp(z1,z2)%rce
          drm2=drm *drm
          drm3=drm2*drm
          drm4=drm3*drm
          drm5=drm4*drm
        
          banda=tbp(z1,z2)%b1*drm5+tbp(z1,z2)%b2*drm4+tbp(z1,z2)%b3*drm3
          banda=banda*banda

        end if
 
        !OMP: Un thread puede estar accediendo al mismo sector del indice i a
        !partir del indice j mas abajo. Me pregunto si es necesario ponerlo aca
        !o no 
        !!$OMP ATOMIC
        band(i)=band(i)+banda

        if(.not.ig%newton) cycle
        if(j>nlocal) cycle

        !!$OMP ATOMIC
        band(j)=band(j)+banda
                  
      enddo 

    enddo    

    !$OMP END PARALLEL DO

  enddo 

  ! band for the ghost
  do i=1+ig%n(1),ig%n(2)
    j = ig%at(i)%o%tag
    band(i) = band(j)
  enddo

  ! Aca tuvo que terminar todas las interacciones con todos
  ! eband=band**(-0.5_dp)
  ! !!$OMP WORKSHARE
  eband(:)=1.0_dp/sqrt(band(:))
  ! !!$OMP END WORKSHARE
                         
end subroutine tb_preinteraction

subroutine tb_interaction(ig)
  class(intergroup),intent(inout)   :: ig
  real(dp)                          :: dr,vd(dm),r0,aux1,aux2,lev_ui,repul
  real(dp)                          :: drm,drm2,drm3,drm4,drm5,factor,factor2(dm),epot
  integer                           :: z1,z2
  integer                           :: i,j,l,n
  type(atom),pointer                :: o1,o2
 
  ! Make ev_ui inside the lexical extent of OMP
  lev_ui=ev_ui
                                
  epot=0._dp

  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP& PRIVATE(i,o1,z1,l,j,o2,z2,vd,dr,n,aux1,aux2,repul)  &
  !$OMP& PRIVATE(factor,factor2,r0,drm,drm2,drm3,drm4,drm5)  &
  !$OMP& SHARED(ig,lev_ui,b_gvirial,virial,eband)       &
  !$OMP& SHARED(tbp,a,band,mic,nlocal,one_box,box)      &
  !$OMP& REDUCTION(+:epot)
                    
  do i = 1,ig%n(1)
     
    o1 => ig%at(i)%o
    z1 = o1%z
  
    do l = 1, ig%nn(i)  ! sobre los vecinos

      j  = ig%list(i,l)
      o2 => ig%at(j)%o
      z2 = o2%z 

      vd = vdistance( o2, o1 , mic) ! respetar el orden
      dr =  dsqrt(dot_product(vd,vd)) 
 
      factor=0.0_dp
  
      if(dr>tbp(z1,z2)%rce) cycle

      if(dr<=tbp(z1,z2)%rci)then       
        
        r0=tbp(z1,z2)%r0  
  
        ! Energia de repulsion
        repul=tbp(z1,z2)%a*exp(-tbp(z1,z2)%p*(dr/r0-1.0_dp))
                                   
        ! Fuerza de repulsion y de banda
        factor=2.0_dp*repul*(tbp(z1,z2)%p/r0)                 &
              -tbp(z1,z2)%eps*tbp(z1,z2)%eps*1.0_dp/2.0_dp    &
               *exp(-2.0_dp*tbp(z1,z2)%q*(dr/r0-1.0_dp))      &
               *(2.0_dp*tbp(z1,z2)%q/r0)*(eband(i)+eband(j))
  
      else if(dr<=tbp(z1,z2)%rce)then 
       
        drm=dr-tbp(z1,z2)%rce
        drm2=drm *drm
        drm3=drm2*drm
        drm4=drm3*drm
        drm5=drm4*drm
      
        ! Fuerza de banda
        factor= (  tbp(z1,z2)%b1*drm5 &
                  +tbp(z1,z2)%b2*drm4 &
                  +tbp(z1,z2)%b3*drm3)&
               *(  5.0_dp*tbp(z1,z2)%b1*drm4 &
                  +4.0_dp*tbp(z1,z2)%b2*drm3 &
                  +3.0_dp*tbp(z1,z2)%b3*drm2)&
               *(eband(i)+eband(j))
  
        ! Fuerza de repulsion
        factor=factor-2.0_dp*(5.0_dp*tbp(z1,z2)%r1*drm4+4.0_dp*tbp(z1,z2)%r2*drm3+3.0_dp*tbp(z1,z2)%r3*drm2)
  
        ! Energia de repulsion
        repul=tbp(z1,z2)%r1*drm5+tbp(z1,z2)%r2*drm4+tbp(z1,z2)%r3*drm3
      end if
      
      factor2(1:dm) = factor *lev_ui * vd(1:dm)/dr
      
      o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 
      o1%epot = o1%epot + repul*lev_ui 
      epot = epot + repul*lev_ui 
           
      if(j<=nlocal) then

        if(.not.ig%newton) cycle

        o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
        o2%epot = o2%epot + repul*lev_ui 
        epot = epot + repul*lev_ui 

        if(b_gvirial) then
          do n = 1,3
            virial(n,:) = virial(n,:) + o2%pos(n)*factor2(:) 
          enddo
        endif

      else if (o2%tag>o1%tag) then  
                
        o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
        o2%epot = o2%epot + repul*lev_ui 
        epot = epot + repul*lev_ui 
                   
        if(b_gvirial) then
          do n = 1,3
            virial(n,:) = virial(n,:) + box(n)*floor(o2%pos(n)*one_box(n))*factor2(:) 
          enddo
        endif
              
      endif

      !!$OMP END CRITICAL
      
    enddo 
    
    aux2=dsqrt(band(i))*lev_ui 

    ! Energia de coohesion por atomo
    o1%epot = o1%epot - aux2

    ! Energia total
    epot = epot - aux2

  enddo

  !$OMP END PARALLEL DO

  ig%epot=epot

end subroutine  

end module gems_tb
