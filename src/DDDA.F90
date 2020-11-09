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

module gems_DDDA
 
! Dynamic Distributable Decorrelation Algorithm (DDDA). 
! Algorithm Ref: Journal of Computaciontal Chemistry Kent et, al, Vol. 28 No. 14.

! Este algoritmo calcula on the fly el error en la media de un dato estadistico.
! Requiere O(N+mlog_2(N)) operaciones para calcular la varianza m veces durante
! una simulacion de N paso
! Si se programa parallelo requiere O(log_2N) de datos a comunicar y el bond
! procedure.

! A medida que el dato es generado durante el calculo, es sumado a un objeto
! de la clase decorrleation. Hay que tener en cuenta que el dato se pierde, en
! el sentido de que no se puede retirar despues (no se puede hacer un buffer
! que vaya avanzando)

! Calculos en paralelo:
! Si la idea es que funcione de manera serial, solo
! un objeto decorrelation basta. Si uno quiere que funcione de manera paralelo,
! cada procesador, que genera sus datos, debe tener su propio objeto de
! decorrelacion y luego se los suma a travez del proceso decorrelation_addition
! para obtener la desviacion estandar de propiedad calculada en paralelo.

! La implementacion que desarrolle , esta basada en el pseudocodigo de
! ese paper. Con algunas modificaciones:

! - En el pseudocodigo hay un error. Un if confuso que no se sab deonde
! termina (o si, pero no deberia temrinar ahi. El if Nsamples>=2 Ŝize termina
! con el tercer append , antes del add. 

! - El waiting_sample y waiting_sample_exist no los hago propiedad de la clase
! Decorrelacion, sino de la Statistic. Esto es para no tener que manejar tantas
! listas dinamicas, sino una sola.
!
! 
!                                                               Alexis Paz

! Como para probar lo escribi en fortran 2003 con el dummy class
use gems_errors
use gems_constants,only: sqrt2,dp,idp
!integer,parameter             :: dp=kind(0.0d0)
!integer,parameter             :: logunit=6


implicit none

private 
public   :: decorrelation
public   :: statistic_dclist

type        :: statistic
  integer   :: nsamples = 0
  real(dp)  :: suma=0.0_dp, sumsq=0.0_dp
  real(dp)  :: var=0.0_dp,errvar=0.0_dp
  real(dp)  :: waiting=0.0_dp
  logical   :: waiting_exists=.false.
  contains
  procedure :: vars => statistic_var
  procedure :: add => statistic_adddata
  procedure :: empty => statistic_empty
  procedure :: statistic_assign, statistic_addition
  generic   :: assignment(=) => statistic_assign
  generic   :: operator(+) => statistic_addition
endtype

! Esta lista va a ser usada dentro de los objetos decorrlation La lista es
! circular para poder rapidamente obtener la desviacion estandar del ultimo
! item   
#define _NODE statistic_dclist
#define _CLASS type(statistic)
#define _ALLOCATABLE
#include "cdlist_header.inc"
    
type :: decorrelation
  integer                         :: size=0, nsamples=0
  integer                         :: fixsize=0
  logical                         :: fixsize_exists=.false.
  type(statistic_dclist),pointer  :: blocks=>null()
  contains
  procedure :: destroy => decorrelation_destroy
  procedure :: empty => decorrelation_empty
  procedure :: init => decorrelation_init
  procedure :: add => decorrelation_adddata
  procedure :: fix => decorrelation_setfixsize
  procedure :: var => decorrelation_variance
  procedure :: med => decorrelation_media
  procedure :: disp => decorrelation_dispersion
  procedure :: write => decorrelation_write
  procedure :: plato
  !generic   :: write(formatted) => decorrelation_write
end type decorrelation


contains


#define _NODE statistic_dclist
#define _CLASS type(statistic)
#include "cdlist_body.inc"
                    
subroutine statistic_assign (a,b )
class(statistic),intent(out):: a
class(statistic),intent(in) :: b

a%nsamples       = b%nsamples       
a%suma           = b%suma           
a%var            = b%var            
a%waiting        = b%waiting        
a%waiting_exists = b%waiting_exists 

end subroutine
 
subroutine statistic_adddata(sb,sample)
class(statistic)       :: sb
real(dp),intent(in)   :: sample
sb%nsamples = sb%nsamples + 1
sb%suma = sb%suma + sample
sb%sumsq = sb%sumsq + sample*sample
end subroutine statistic_adddata

subroutine statistic_empty ( sb )
! Esta es una lista doble circular
class(statistic)       :: sb

sb%nsamples = 0
sb%suma=0.0_dp
sb%sumsq=0.0_dp
sb%var=0.0_dp
sb%errvar=0.0_dp
sb%waiting=0.0_dp
sb%waiting_exists=.false. 

end subroutine statistic_empty
 
function statistic_addition(a,b) result(c)
class(statistic),intent(in)   :: a,b
type(statistic)               :: c
c%nsamples = a%nsamples + b%nsamples
c%suma = a%suma + b%suma
c%sumsq = a%sumsq + b%sumsq
end function
 
function statistic_var(sb) result(var)
class(statistic)               :: sb
integer                        :: nsamples_1
real(dp)                       :: var,media

media=sb%suma/sb%nsamples
nsamples_1=(sb%nsamples-1)

var=sb%sumsq/sb%nsamples-media*media
 
! Hasta aca la varianza. La desviacion estandar (sd) de la propiedad (media)
! sera sqrt(var) El error estandar (SE) de la propiedad es sd sobre raiz del numero
! de muestras (n) por el factor de confianza (supongamos 1.92 para el 95% de
! confianza)
! SE=sd/sqrt(n)=sqrt(var/n)
!  Upper 95% Limit = \av{x} + (SE*1.96)
!  Lower 95% Limit = \av{x} - (SE*1.96) 
             
!  Flyvbjerg, H., & Petersen, H. G. (1989). Error estimates on
!  averages of correlated data. The Journal of Chemical Physics,
!  91(1), 461. doi:10.1063/1.457480

! SE tambien tiene su propio error estandar. La  desviacion estandar de SE 
! es la misma SE por raiz de (dos divido n-1).
! SESE=SE*dsqrt(2.0_dp/nsamples-1)

end function
         


subroutine decorrelation_init ( d )
class(decorrelation)    :: d
if(associated(d%blocks)) return

! Create the statistics list
allocate(d%blocks)
call d%blocks%init()
allocate(d%blocks%o)

end subroutine decorrelation_init
 
subroutine decorrelation_destroy ( d )
class(decorrelation)    :: d

call d%empty()

! Borro la lista
deallocate(d%blocks%o)
deallocate(d%blocks)
d%blocks => null()
 
end subroutine decorrelation_destroy 


subroutine decorrelation_empty ( d )
class(decorrelation)    :: d

d%size=0
d%nsamples=0
d%fixsize=0
d%fixsize_exists=.false.

! Borro la lista
call d%blocks%destroy_all()

end subroutine decorrelation_empty

subroutine decorrelation_setfixsize(d,fixsize)
! Si ya se cual es el tamaño de bloque que devuelve datos descorrelacionados,
! me basta con calcular propagar hasta ese tamaño.
class(decorrelation)  :: d
integer,intent(in)    :: fixsize

call dddawarning(2001,d%fixsize_exists) 
d%fixsize = fixsize
d%fixsize_exists = .true.

! Aca podria preguntarme si tengo una lista mas larga y cortarla. 

end subroutine decorrelation_setfixsize
 
subroutine decorrelation_adddata(d,sample)
class(decorrelation)            :: d
real(dp)                        :: sample
real(dp)                        :: carry
class(statistic_dclist),pointer  :: ls
logical                         :: lengthen
integer                         :: i

!Es util la siguiente representación
!nodos            -(head)   -   -   -   x  x
!nsamples que
!aumentan size              1   2   4   8  16
!size                       1   2   3   4  5 
!bloque              1      2   4   8  16  32
 
! Nuevo numero de muestras
d%nsamples = d%nsamples + 1

! Test por si tengo que agrandar el rancho
lengthen = (d%nsamples >= 2**d%size)

! Si no quiero propagar mas alla de un dado tamaño
if(d%fixsize_exists) then
  if(d%fixsize<d%size) then
    lengthen = .false.
    d%size = d%fixsize
  endif
endif

! Lengthen the vectors, when necessary to accommodate all entered data
if (lengthen) then
  d%size = d%size + 1
  call d%blocks%add_before()
  call d%blocks%prev%alloc()
endif

! Agrego el dato al nivel cero d%blocks. El nivel cero siempre acepta.
call d%blocks%o%add(sample)
carry = sample

! Propagate the new sample up trought the data structure
ls=>d%blocks
do i=1,d%size
  ls=>ls%next

  ! Pregunto si este objeto statistic es de capa completa (nsamples par)
  if ( ls%o % waiting_exists ) then

    ! Computo el promedio entre el que faltaba con el nuevo y lo guardo
    sample = (ls%o % waiting + carry)*0.5_dp
    carry = sample

    ! Agrego el promedio al objeto
    call ls%o%add(sample)
    ls%o % waiting_exists = .false.

  else

    ! Guardo el valor para la espera
    ls%o % waiting_exists = .true.
    ls%o % waiting = carry
    exit
  endif

enddo

end subroutine decorrelation_adddata
 
subroutine decorrelation_variance(d,n,err)
!  Flyvbjerg, H., & Petersen, H. G. (1989). Error estimates on
!  averages of correlated data. The Journal of Chemical Physics,
!  91(1), 461. doi:10.1063/1.457480

! Cuando uno hace un grafico de Flyvbjerg-Petersen, del estimador de la
! desviacion estandar con su error en funcion del log_2(tamaño de bloque)
! encuentra esto:

! 0.0013 ++----------------------------------------------------------------+
!        |                                                    ***          |
! 0.0012 ++                                               ***  *      ***  |
! 0.0011 ++                                                *   *       *   |
!        |                                       ***       *   *  ***  *   |
!  0.001 ++   err=                       *** ***  A  ***   A   *   *   *   |
!        |   sqrt(var/N-1)           ***  A   A   *   A    *   A   *   *   |
! 0.0009 ++                          *A* *** *** ***  *    *   *   *   *   |
!        |                       *A*                 ***  ***  *   *   *   |
! 0.0008 ++                                                    *   A   A   |
! 0.0007 ++                  *A*                              ***  *   *   |
!        |                                                         *   *   |
! 0.0006 ++                          errerr=err/(sqrt(2(n-1))      *   *   |
!        |               *A*                                       *   *   |
! 0.0005 ++                                                       ***  *   |
!        |          *A*                                                *   |
! 0.0004 ++                                                           ***  |
! 0.0003 ++     *A*                                                        |
!        +  *A*  +        +       +       +       +        +       +       +
! 0.0002*A*------+--------+-------+-------+-------+--------+-------+-------+
!        0       2        4       6       8       10       12      14      16

! La idea de esta funcion es un criterio para definir on the fly el plato del
! grafico. Primero parto del bloque mas grande y con mayor error y voy
! intersectando el intervalo de error hacia los bloques mas pequeños. Cuando
! la interseccion es vacia, devuelvo el numero de tamaños de bloque que se
! intersectaron. 
! Ademas la subrrutina devuelve la varianza obtenida con este metodo.
! No se si existe otro metodo mejor, con algun fiteo, pero algo es algo.

! Imprime la penultima. La ultima da Nan al dividir por N-1.
class(decorrelation),intent(in) :: d
real(dp),intent(out)            :: err
integer,intent(out)             :: n
class(statistic_dclist),pointer :: ls
real(dp)                        :: errerr,y1,y2,y1n,y2n,y1o,y2o
integer                         :: i

! Pocos datos para calcular nada
if(d%size<3) return 

y1o=-1e16_dp
y2o=1e16_dp
n=0

! El ultimo objeto statistics siempre tiene nsamples=0, asique lo descarto
! El anteultimo tiene nsamples=1 y para calcular la varianza uno hace /N-1
! asique tampoco sirve. El anterior a ese tambien tiene aveces nsample=1,
! asique lo descarto tambien
ls => d%blocks%prev%prev

 !Ejemplo 1. Supongamos que tengo 15 datos. 
 !nsize=4 pero tengo 5 nodos (contando la cabeza)
 !  id      n bloque    nsamples    descartar
 !  1 (head)   1          15           v
 !  2          2          7            v
 !  3          4          3            v
 !  4          8          1            x
 !  5(prev)   16          0            x
          
 !Ejemplo 2. Supongamos que tengo 18 datos. 
 ! nsize=5 pero tengo 6 nodos(contando la cabeza)
 !  id        n bloque    nsamples    descartar 
 !  1 (head)     1         18            v
 !  2            2          9            v      
 !  3            4          4            v      
 !  4            8          2            v     
 !  5           16          1            x      
 !  6 (prev)     -          0            x      

 !Es util la siguiente representación
 !nodos            -(head)   -   -   -   x  x
 !nsamples que
 !aumentan size              1   2   4   8  16
 !size                       1   2   3   4  5 
 !bloque              1      2   4   8  16  32


!  Comienzo las intersecciones con los conjuntos anteriores 
!Este loop tiene que ser -1 porque hay d%size+1 nodos
!El cero tambien cuenta
!Tiene que ser -1... nos e porque anda mejor el -2
do i =1,d%size-2
  ls => ls%prev

  !Segun la definicion de Chi como la varianza sobre n-1 o el errerror estandar
  !Kent, I., David, R., Muller, R. P., Anderson, A. G., Goddard III, W. A., &
  !Feldmann, M. T. (2007). Efficient Algorithm for “‘On-the-Fly’” errerror
  !Analysis of Local or Distributed Serially Correlated Data. Journal of
  !Computational Chemistry, 28(14), 2309–2316. doi:10.1002/jcc

  !!chi es el error estandar al cuadrado
  !chi = ls%o%vars()/(ls%o%nsamples-1)  !Chi

  !! Desviacion estandar de chi
  !chi_err=chi*sqrt(2.0_dp/(ls%o%nsamples-1)) 

  !! Error en chi
  !chi_err=chi_err/sqrt(ls%o%nsamples-1)

  !! Aca uso la propagacion de errores
  !sqrtchi = sqrt(chi)                       
  !sqrtchi_err=sqrtchi*0.5_dp*chi_err/chi    

  ! RESUMEN
  !chi = ls%o%vars()/(ls%o%nsamples-1)
  !chi_err=chi*sqrt(2/(ls%o%nsamples-1))
  !sqrtchi = sqrt(chi)                        
  !sqrtchi_err=0.5_dp*chi_err/sqrtchi         
  !sqrtchi = sqrt(ls%o%vars()/(ls%o%nsamples-1))
  !sqrtchi_err= 0.5_dp*sqrtchi*sqrt(2/(ls%o%nsamples-1))
 
  !y como la sqrtchi es el error de la medicion entonces
  err = sqrt(ls%o%vars()/(ls%o%nsamples-1))
  !errerr= 0.5_dp*err*sqrt2/(ls%o%nsamples-1)

  !Esta es lo que usan en el paper de FP
  errerr=err/sqrt(2.0_dp*(ls%o%nsamples-1))


  ! Nuevo intervalo
  y1n=err-errerr
  y2n=err+errerr

  ! Me quedo con el rango que mantenga la interseccion
  y1=max(y1o,y1n)
  y2=min(y2o,y2n)

  ! Me fijo si la interseccion es vacia
  if((y2-y1)<=0.0_dp) exit

  y1o=y1
  y2o=y2

  ! Avanzo el intervalo
  n=n+1

enddo

! Encuentro el plato en el valor que esta al medio del intervalo final
err = (y2o+y1o)*0.5

end subroutine decorrelation_variance

!subroutine decorrelation_variance2(d,f,a,chi2,sigma)
!  !Lo mismo que lo anterior pero tratando de mejorar el criterio.  Se toman los
!  !ultimos datos del grafico y se ajusta una recta. Se informa el valor de la
!  !pendiente, lo que ayuda a definir un criterio para plato o no plato. La recta
!  !se ajusta teniendo en cuenta el error. Para saber cuantos datos del final
!  !tomar se usa f. f es un parametro suministrado como la fraccion del grafico a
!  !utilizar para defeinir si es o no un plato.
!  use nr,only: fit
!  class(decorrelation),intent(in) :: d
!  real(dp),intent(in)             :: f
!  real(dp),intent(out)            :: sigma
!  class(statistic_dclist),pointer :: ls
!  real(dp)                        :: err,y1,y2
!  real(dp)                        :: a,b,sigb,siga,chi2,q
!  real(dp),allocatable            :: x(:),y(:),sig(:)
!  integer                         :: i,n
!
!  ! Pocos datos para calcular nada
!  if(d%size<3) return 
!
!  y1=-1e16_dp
!  y2=1e16_dp
!  n=0
!
!  ! El ultimo objeto statistics siempre tiene nsamples=0, asique lo descarto
!  ! El anteultimo tiene nsamples=1 y para calcular la varianza uno hace /N-1
!  ! asique tampoco sirve. El anterior a ese tiene aveces nsample=1, asique lo
!  ! descarto tambien
!  ls => d%blocks%prev%prev%prev
!
!  n=f*d%size
!  allocate(x(n))
!  allocate(y(n))
!  allocate(sig(n))
!
!  !  Comienzo las intersecciones con los conjuntos anteriores 
!  do i =1,n
!    ls => ls%prev
!
!    ! Me quedo con el rango que mantenga la interseccion
!    x(i) = i
!    y(i) = sqrt(ls%o%vars)
!    sig(i) = sqrt(err)
!  
!    call fit(x,y,a,b,siga,sigb,chi2,q,sig)
!
!  enddo
!
!  ! Encuentro el plato como el valor de y de la recta en la mitad de su dominio
!  sigma = a*x(int(n*0.5))+b
!
!
!end subroutine decorrelation_variance2
!
function plato(d,f) result(res)
! Devuelve true si se encontro un plato en el objeto
class(decorrelation)             :: d
logical                          :: res
real(dp)                         :: var
integer                          :: n
real(dp),intent(in)              :: f

call d%var(n,var)
res=( float(n)/(d%size-2) > f )

end function plato
 
function decorrelation_media(d) result(res)
class(decorrelation)             :: d
real(dp)                         :: res
class(statistic_dclist),pointer  :: ls

ls => d%blocks
res = ls%o%suma/ls%o%nsamples

end function decorrelation_media
 
function decorrelation_dispersion(d) result(res)
class(decorrelation)             :: d
real(dp)                         :: res
class(statistic_dclist),pointer  :: ls

ls => d%blocks
! FIXME PORQUE ESTABA ESTO??? 
!res = sqrt(ls%o%vars()*(ls%o%nsamples-1))
res = sqrt(ls%o%vars())

end function decorrelation_dispersion
        
subroutine dddawarning(code,logica)
use gems_errors,only: logunit
integer                :: code
logical,intent(in)     :: logica

if(.not.logica) return

call wstd(); write(logunit,advance='no',fmt='(a)') '# GMD(DDDA): '

selectcase(code)
case(2001) 
  call wstd(); write(logunit,*) 'IGNORED!!!!! you cant fix the box twice in the same run (DDDA is on the fly!) '
endselect 

endsubroutine dddawarning


! Esta subrrutina sirve para cuando se paraleliza.
! Lo dejo para despues.
! 
!function decorrelation_addition(a,b) result(c)
!  type(decorrelation),intent(in)   :: a,b
!  type(decorrelation)              :: c
!  type(statistic)                  :: statA
!  logical                          :: waiting_existA
!  real(dp)                         :: waitingA
!  type(statistic)                  :: statB
!  logical                          :: waiting_existB 
!  real(dp)                         :: waitingB
!  logical                          :: carry_exist
!  real(dp)                         :: carry
!  integer                          :: i
!  
!  c%nsamples = a%nsamples + b%nsamples
!
!  ! Make c big enough to hold all data from a and b
!  do while(c%nsamples >= 2.0_dp*c%size)
!    c%size=c%size+1
!    call statistic_list_append(c%blocks)
!    call real_list_append(c%waiting)
!    call logical_list_append(c%waiting_exists)
!    carry_exist = .false.
!    carry=0.0_dp
!    
!    lsa=>a%blocks
!    lra=>a%waiting
!    lla=>a%waiting_exists
!    lsb=>b%blocks
!    lrb=>b%waiting
!    llb=>b%waiting_exists
!    lsc=>c%blocks
!    lrc=>c%waiting
!    llc=>c%waiting_exists
!    do i = 0,c%size-1
!
!      ! Estos alias los genero para mejorar la legibilidad
!      Ablocks             =>lsa%s
!      Awaiting      =>lra%r
!      Awaiting_exists=>lla%l
!      Bblocks             =>lsb%s
!      Bwaiting      =>lrb%r
!      Bwaiting_exists=>llb%l
!      Cblocks             =>lsc%s
!      Cwaiting      =>lrc%r
!      Cwaiting_exists=>llc%l
!
!      if (i<=a%size) then
!        statA = Ablocks
!        waitingA = Awaiting
!        waiting_existsA = Awaiting_exists
!      else
!        statA = !cero!
!        waitingA = 0.0_dp
!        waiting_existsA = .false.
!      endif
!        
!      if (i<=b%size) then
!        statB = Bblocks
!        waitingB = Bwaiting
!        waiting_existsB = Bwaiting_exists
!      else
!        statB = !cero!
!        waitingB = 0.0_dp
!        waiting_existsB = .false.
!      endif
!        
!      Cblocks=statistic_addition(statA,statB)
!
!      if (carry_exist.and.waiting_existA.and.waiting_existB) then
!        ! Three samples to handle
!        call statistica_adddata(Cblocks,(waitingA+waitingB)*0.5_dp)
!        Cwaiting = 0.0_dp
!        Cwaiting_exists = .false.
!        carry_exist=.true.
!        carry=(waitingA+waitingB)*0.5_dp
!      elseif (carry_exist.and..not.waiting_existA.and.waiting_existB) then
!        ! Two samples to handle
!        call statistica_adddata(Cblocks,(carry+waitingB)*0.5_dp)
!        Cwaiting = 0.0_dp
!        Cwaiting_exists = .false.
!        carry_exist=.true.
!        carry=(carry+waitingB)*0.5_dp
!      elseif (carry_exist.and.waiting_existA.and..not.waiting_existB) then
!        ! Two samples to handle
!        call statistica_adddata(Cblocks,(carry+waitingA)*0.5_dp)
!        Cwaiting = 0.0_dp
!        Cwaiting_exists = .false.
!        carry_exist=.true.
!        carry=(carry+waitingA)*0.5_dp
!      elseif (carry_exist.or.waiting_existA.or.waiting_existB) then
!        ! One samples to handle
!        Cwaiting = carry + waitingA + waitingB
!        Cwaiting_exists = .true.
!        carry_exist=.false.
!        carry = 0.0_dp
!      else
!        ! No samples to handle
!        Cwaiting = 0.0_dp
!        Cwaiting_exists = .false.
!        carry_exist=.false.
!        carry = 0.0_dp
!      endif
!
!      lsa=>lsa%next
!      lra=>lra%next
!      lla=>lla%next
!      lsb=>lsb%next
!      lrb=>lrb%next
!      llb=>llb%next
!      lsc=>lsc%next
!      lrc=>lrc%next
!      llc=>llc%next
!    enddo
!      
!      
!  enddo
!
!end function
!

!subroutine decorrelation_write(d, iotype, v_list, iostat, iomsg )
subroutine decorrelation_write(d, un)
use gems_constants, only: ui_ev
! Realiza el grafico de Flyvbjerg and Petersen: desviacion estandar en funcion
! de tamaño de bloque estadistico. El plato en el grafico indica que el tamaño
! de bloque consigue perder la correlacion
class(decorrelation),intent(in)      :: d
integer,               intent(in)    :: un
!character(len=*),      intent(in)    :: iotype
!integer, dimension(:), intent(in)    :: v_list
!integer,               intent(out)   :: iostat
!character(len=*),      intent(inout) :: iomsg
integer                        :: j,n
class(statistic_dclist),pointer :: ls
real(dp)                       :: var,err


!call d%var(n,var)
!write(un,fmt='(i0,a,i0,f5.2)') n,' puntos de ',d%size-2,float(n)/(d%size-2)
!write(un,fmt='(4(a,e25.12))')  'media: ',d%med(),'+-',var,'dispersion: ',sqrt(d%disp())

! El siguiente do hace log_2(N)-1 lineas, donde N es el numero de puntos
! agregados en el d. El menos uno es porque para calcula var se divide en
! N-1 y esto provocaria un NAN 
 
n=d%blocks%o%nsamples

ls =>d%blocks%prev
do j =1,d%size-2
  ls => ls%next
  var = ls%o%vars()
  err = sqrt(var)/(ls%o%nsamples-1)*ui_ev
  write(un,fmt='(i0,1x,2(e25.12,2x))') ls%o%nsamples,err,err*sqrt2/(ls%o%nsamples-1)
enddo

end subroutine decorrelation_write
          
end module gems_DDDA
