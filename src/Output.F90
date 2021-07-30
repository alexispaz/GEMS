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

! Object outpropa collect a set of properties, the way to compute them, the
! group of atoms needed to compute them, the DDDA associated object

! Outfile collect a vector of properties and the information of the file to
! ouptut them as columns OR collect a particular property that needs a formated
! output.

! The pointer assigment for writing procedures for both outfile and outprop is
! in Main.f90. It is not possible to include is here.  Since different modules
! might have write_??? subrroutines, it is not possible to include the pointer
! assigment to this subrroutines inside output.f90. The thing is that this
! modules needs elements in output module to construct write_ subrroutines and
! the output module needs the write_ subrroutines from those modules.

! TODO: Me parece que se podrían sacar los outfiles y dejar solo las outprop
! que tengan una unidad asociada y que se haga estilo:

module gems_output
 use gems_inq_properties
 use gems_program_types
 use gems_groups
 use gems_constants
 use gems_strings
 use gems_set_properties
 use gems_tables
 use gems_input_parsing
 use gems_ddda
 use gems_elements, only: ncsym,csym

implicit none

! TODO: generalizar outpropa para apuntar a grupos, ngroups, interactions
!       quizas se puede declarar distintos tipos de outpropa y que el outfile
!       apunte a todas ellas. No se si hacer un deferred porque hay muchas propiedades
!       parecidas (epot, ecin, temp, etc) y no vale la pena apuntar a todas ellas, pero si 2 o 3 objetos heredaros, uno para
!       ngroups, otro para integrations, etc. otro para cosas como pos, vel, etc que se escriben sin buffer, otras con buffer.
!       todas que esten dentro de un outfile, y que sea el outfile el que distinga cuando le ponen de un tipo incompatible con otro
!       tipo (e.g. epot con pos)
type :: outpropa
  type(group),pointer                :: g=>null()     ! Uno de los grout(:)
  !TODO integer                      :: units=0       ! Unidades. 0-A ps, eV 1-A ps kcm, etc
  integer                            :: n=0           ! Numero de columnas
  real(dp),pointer                   :: f(:)=>null()  ! Las columnas
  type(decorrelation),pointer        :: d(:)=>null()  ! El objeto decorrelation del DDDA- 
  procedure(write_outpropa),pointer  :: w=>null()     ! Regla de calculo

  contains
    procedure        :: outprop_assign
    procedure        :: init => outprop_init
    procedure        :: destroy => outprop_destroy 
    ! procedure        :: outprop_write
    generic          :: assignment(=) => outprop_assign
    ! generic          :: write(formatted) => outprop_write
    ! final            :: outprop_destroy
endtype
                
interface 
  subroutine write_outpropa(op)
    import
    class(outpropa)     :: op
  end subroutine
end interface
 
#define _NODE outpropa_l
#define _CLASS type(outpropa)
#include "list_header.inc"

type :: outfile

  ! Esto es comun a los arhcivos
  logical                  :: open=.false.       ! El archivo esta abierto
  logical                  :: enable(10)=.false. ! El archivo esta activo, discriminado por lugar de escritura
  integer                  :: ecode=0            ! Numero binario que indica el estado de enable (que depende de 'at' command)
  integer                  :: outeach=1          ! Cada cuantos pasos de dm_step se imprime
  character(linewidth)     :: name=''            ! El nombre del archivo
  integer                  :: un=0               ! La unidad del archivo
  procedure(ofw),pointer   :: w=>null()          ! Regla de escritura
  logical                  :: flush=.false.      ! Regla de flushing
                                            
  type(outpropa_l),pointer :: p=>null()          ! properties
  integer                  :: sumeach=9999999    ! Cada cuantos pasos de dm_step se calcula (e.g. promedios, ddda)

  !Eso es para los archivos con decorrelacion
  logical                  :: ddda=.false.  ! Indica si se esta usando ddda en este archivo
  integer                  :: und=0         ! La unidad del archivo de decorrelacion

  ! Esto es para los archivos promedio
  logical                  :: prom=.false.  ! promedio
  integer                  :: n=0           ! los puntos sumados

  !Esto es para los archivos en bloque
  type(group),pointer      :: g=>null()     ! Uno de los grout(:)

  contains

    procedure :: reopen => outfile_reopen

endtype
             
interface
  subroutine ofw(of)
    import outfile
    class(outfile)          :: of
  end subroutine
end interface

#define _NODE outfile_aop
#define _CLASS class(outfile)
#include "arrayofptrs_header.inc"

#define _NODE outfile_vop
#define _TYPE type(outfile_aop)
#include "vector_header.inc"

! VOP to collect the ofiles declared inside modules
type(outfile_vop),public :: of_vop
             
! Precision de salida
character(:),allocatable,public    :: prf!='20.7'
character(:),allocatable,public    :: pri!='0'
  
public :: outfile, polvar_outfile
public :: outfile_write,outfile_prom,outfile_ddda

! Checkpoint mode, para abrir los archivos pero no sobreescribir
logical,public        :: chpmode=.false.   

integer              :: fstart=1,fend=1e8,fstep=1, npoints
! character(1)         :: tips
character(ncsym)      :: symbol
logical              :: read_error
real(dp)             :: long
integer              :: contador


! Unidades de salida
real(dp)                :: time_out=1.0_dp

! Dinamica acelerada
integer                 :: bbunit

! Sistema de buffer
logical                 :: vbuffer=.false.
integer                 :: kvel,vresto,vstage
real(dp)                :: dtframe=0.0_dp
                                
! Sistema de promedio
integer                 :: fprom=1
logical                 :: geterr=.false.

! Auxiliares
integer                 :: disrad_i
real(dp)                :: disrad_f
complex(dp),allocatable :: caux(:)
real(dp)                :: eini

!Archivos de salida

! Para arreglar esto habria que armar una pila con estas salidas, y asi recien
! poder distribuir cada una a su modulo correpondiente-


logical        :: get_dihedral=.false.,op_dihedral=.true.

interface write_out
   module procedure :: write_out_int,write_out_real
end interface

contains 

#define _NODE outpropa_l
#define _CLASS type(outpropa)
#include "list_body.inc"
                  
#define _NODE outfile_vop
#define _TYPE type(outfile_aop)
#include "vector_body.inc"
                 
! El vector necesita de la asignacion
subroutine outprop_assign (a,b )
class(outpropa),intent(out):: a
class(outpropa),intent(in) :: b

a%g   =>b%g
if(associated(b%f)) allocate(a%f(size(b%f)))
a%f   = b%f
!a%d  = b%d
a%w   =>b%w
a%n   = b%n

end subroutine

! subroutine outprop_basics ( )
! type(outprop),pointer :: op
!   
! call op_kv%ensure(50)
! call op_v%ensure(50)
!
! call  op_kv%put(1,'time'); op=>op_v%o(1);   op%n=1,  op%w=>write_time
! call opg_kv%put(1,'temp'); opg=>opg_v%o(1); opg%n=1, opg%w=>write_temp
! call op3_kv%put(1,'pos');  op3=>op3_v%o(1); op3%n=1, op3%w=>write_temp
! ... etc
! 
! To use in CLI:
! read(w1)
! i=op_kv%find(w1)
! if (i/=0) then
!   op%n=op_v%o(i)%n
!   op%n=op_v%o(i)%n
!   call op%init(op_v%o(i)%n,op_v%o(i)%w)
!   return  
! endif
!
! i=opg_kv%find(w1)
! if (i/=0) then
!   opg%n=op_v%o(i)%n
!   opg%n=op_v%o(i)%n
!   call opg%init(op_v%o(i)%n,op_v%o(i)%w,op_v%g)
!   return  
! endif
!
! end subroutine

subroutine outprop_init(op,n,w,g)
class(outpropa)             :: op
integer,intent(in)          :: n
procedure(write_outpropa)   :: w
type(group),target,optional :: g

op%n=n
allocate(op%f(n))
op%f=0.0_dp
op%w => w
if(present(g)) op%g => g
          
end subroutine

subroutine outprop_destroy (op )
class(outpropa)   :: op
! type(outpropa)   :: op
if(associated(op%d)) deallocate(op%d)
if(associated(op%f)) deallocate(op%f)
op%g=>null()   
op%n=0
op%w => null()
end subroutine

subroutine outfile_write (of)
class(outfile)              :: of
type(outpropa_l),pointer    :: ln
character(100)              :: fm
     
! Calculo y escribo todas las propiedades 
ln => of%p
do while( associated(ln%next) )
  ln => ln%next

  call ln%o%w

  fm='('//trim(.ich.ln%o%n)//'(x,e'//trim(prf)//'))'

  write(of%un,fmt=fm,advance='no') ln%o%f(:)
  ln%o%f(:)=0.0_dp
enddo

write(of%un,*)

if(of%flush) call flush(of%un)

end subroutine outfile_write

subroutine outfile_prom (of)
class(outfile)              :: of
type(outpropa_l),pointer    :: ln
character(100)              :: fm
     
! Calculo todas las propiedades 
ln => of%p
do while( associated(ln%next) )
  ln => ln%next
  call ln%o%w
enddo
     
of%n=of%n+1

! Me fijo si se cumple el modulo
if(mod(of%n,of%outeach)/=0) return

! Divido por n
ln => of%p
do while( associated(ln%next) )
  ln => ln%next
  ln%o%f(:)=ln%o%f(:)/of%n
enddo 
of%n=0

! Escribo todas las columnas
ln => of%p
do while( associated(ln%next) )
  ln => ln%next
  fm='('//trim(.ich.ln%o%n)//'(x,e'//trim(prf)//'))'

  write(of%un,fmt=fm,advance='no') ln%o%f
  ln%o%f(:)=0.0_dp
enddo

write(of%un,*)

if(of%flush) call flush(of%un)

end subroutine outfile_prom

subroutine outfile_ddda (of)
class(outfile)              :: of
type(outpropa_l),pointer    :: ln, lk
class(statistic_dclist),pointer :: ls
type(outpropa),pointer      :: op
real(dp)                    :: err,errerr
integer                     :: n,l,j
     
! Calculo todas las propiedades 
ln => of%p
do while( associated(ln%next) )
  ln => ln%next
  op => ln%o
  call op%w
  do l = 1,size(op%d)
    call op%d(l)%add( op%f(l) )
    op%f(l)=0._dp
  enddo 
enddo
     
of%n=of%n+1

! Me fijo si se cumple el modulo
if(mod(of%n,of%outeach)/=0) return

ln => of%p
do while( associated(ln%next) )
  ln => ln%next
  op => ln%o
  do l=1,size(op%d)
    write(of%un,fmt='(e25.12,1x)',advance='no') op%d(l)%med()
  enddo 
enddo
write(of%un,*)  
     

! loop between window blocks
ln => of%p%next
do j =1,ln%o%d(1)%size-2
       
  !select to the current window block
  ls => ln%o%d(1)%blocks%prev
  do n=1,j
    ls => ls%next
  enddo
  write(of%und,fmt='(i10,1x)',advance='no')  ls%o%nsamples
 
  !loop trough properties
  lk => of%p
  do while( associated(lk%next) )
    lk => lk%next
    op => lk%o
    
    !loop trough each ddda column
    do l =1,size(op%d)

      !select to the current window block
      ls => op%d(l)%blocks%prev
      do n=1,j
        ls => ls%next
      enddo
  
      err = sqrt(ls%o%vars()/(ls%o%nsamples-1))
      errerr=err/(sqrt(2.0_dp*(ls%o%nsamples-1)))

      write(of%und,fmt='(2(e25.12,2x))',advance='no')  err,errerr
    enddo

    write(of%und,*)
 
  enddo 
enddo

! Separador de bloque (Enable `plot "file" i 3` in gnuplot)
write(of%und,*)
write(of%und,*)

if(of%flush) then
  call flush(of%un)
  call flush(of%und)
endif

end subroutine outfile_ddda

subroutine outfile_reopen ( of )
use gems_errors
use gems_algebra, only:decode_logicalvector
use gems_constants, only: find_io
class(outfile)    :: of
type(outpropa_l),pointer :: ln
integer           :: j

if(of%open) close(of%un)

if(.not.associated(of%p)) allocate(of%p)

! Si no tiene nombre, error
call werr('El archivo necesita un nombre',trim(of%name)=='') 

! Si no tiene posicion le doy una ahora
if(.not.any(of%enable))  of%enable(1)=.true.
                                         
of%un = find_io(30)
of%open=.true.
 
if(chpmode) then
  open( of%un , file=trim(adjustl(of%name)),action="write",status="replace", position='append' )
else
  open( of%un , file=trim(adjustl(of%name)),action="write",status="replace" )
endif

of%ecode=decode_logicalvector(of%enable)

! If dda
if(of%ddda) then
                
  if(of%und/=0) close(of%und)
  of%und = find_io(30)
  of%open=.true.
   
  if(chpmode) then
    open( of%und , file='FP_'//trim(adjustl(of%name)), position='append' )
  else
    open( of%und , file='FP_'//trim(adjustl(of%name)) )
  endif
            
  ln => of%p
  do while( associated(ln%next) )
    ln => ln%next
    if(associated(ln%o%d)) cycle
    allocate(ln%o%d(size(ln%o%f)))
    do j = 1,size(ln%o%d)
       call ln%o%d(j)%init()
    enddo
  enddo 

endif


end subroutine outfile_reopen
 
subroutine del_flag ( flag )
character(*),intent(in)   :: flag
integer                   :: i,j,k
type(outfile),pointer     :: of

! Abro el out de nuevo
!close( outunit )
!open( outunit, file=trim(ioprefix)//".out" ,form='unformatted' )
!call wstd(); write(logunit,*) "-redirecting out to ",trim(ioprefix)//".out" 
do i = 1,of_vop%size
  of => of_vop%o(i)%o

  if(.not.any(of%enable)) cycle

  j=index(of%name,flag)
  k=j+len(flag)
  of%name = of%name(:j-2)//of%name(k:)

  call of%reopen()
enddo
!op_dihedral=.true.

end subroutine 

subroutine set_flag ( flag )
character(*),intent(in)   :: flag
integer                   :: i
type(outfile),pointer     :: of

! Abro el out de nuevo
!close( outunit )
!open( outunit, file=trim(ioprefix)//".out" ,form='unformatted' )
!call wstd(); write(logunit,*) "-redirecting out to ",trim(ioprefix)//".out" 
do i = 1,of_vop%size
  of => of_vop%o(i)%o

  if(.not.any(of%enable)) cycle

  ! Le agrego el flag
  of%name=trim(adjustl(of%name))//"."//trim(adjustl(flag))

  call of%reopen
enddo
!op_dihedral=.true.

end subroutine 

!ESCRITURA SALIDA

subroutine write_out_real(j,n)
real(dp),intent(in)      :: n
integer,intent(in)       :: j  ! Para discriminar si es llamada de distintoslugares
integer                  :: i
type(outfile),pointer     :: of

! Mecanismo de checkpoint. No escribo hasta que no se haya leido.
if(chpmode) return

! Promedios y demas
do i = 1,of_vop%size
  of => of_vop%o(i)%o

  ! Ignorar archivos cerrados
  if(.not.of%open) cycle

  ! Me fijo si se cumple el modulo de sampling
  if(int(mod(n,real(of%sumeach,dp)))/=0) cycle

  ! Computar si se encuentra abierto
  if(of%enable(j)) call of%w()

enddo

endsubroutine

subroutine write_out_int(j,n)
integer,intent(in),optional :: n
integer,intent(in)      :: j  ! Para discriminar si es llamada de distintoslugares
integer                 :: i,l
type(outpropa_l),pointer :: ln
type(outfile),pointer     :: of

! Mecanismo de checkpoint. No escribo hasta que no se haya leido.
if(chpmode) return

! Promedios y demas
do i = 1,of_vop%size
  of => of_vop%o(i)%o

  if(.not.of%open) cycle

  ! Calculo los ddda
  if(of%ddda) then
    ln => of%p
    do while( associated(ln%next) )
      ln => ln%next
      call ln%o%w()
      do l = 1,size(ln%o%d)
        call ln%o%d(l)%add( ln%o%f(l) )
      enddo 
    enddo 
  endif

  ! Si promedio en cada paso
  if(of%prom) then 
    ! Calculo todas las propiedades 
    ln => of%p
    do while( associated(ln%next) )
      ln => ln%next
      call ln%o%w()
    enddo
    of%n=of%n+1
  endif

enddo       

do i = 1,of_vop%size
  of => of_vop%o(i)%o

  if(.not.of%open) cycle

  ! Me fijo si se cumple el modulo
  if(present(n)) then 
    if(mod(n,of%sumeach)/=0) cycle
  endif

  if(of%enable(j)) call of%w()

  ! Dejo correr
  call flush(of%un)

enddo

endsubroutine

subroutine write_out_force(j)
integer,intent(in)       :: j  ! Para discriminar si es llamada de distintoslugares 
integer      :: i
type(outfile),pointer     :: of

frame=frame+1
pnframe=frame
ptime=time*time_out

! Mecanismo de checkpoint. No escribo hasta que no se haya leido.
if(chpmode) return

do i = 1,of_vop%size
  of => of_vop%o(i)%o
  if(of%enable(j))  call of%w()
enddo

endsubroutine

subroutine write_screenshot(archivo,g)
type(group),intent(in)      :: g
character(*)                :: archivo
type (atom_dclist),pointer  :: la
integer                     :: i,u
real(dp)                    :: f(3)

f(:)=0._dp 

u = find_io(30)

open(u, file=trim(adjustl(archivo)),position='rewind')
write(u,*) g%nat
write(u,*)
la => g%alist
do i = 1,g%nat
  la => la%next
  f(1:dm)=la%o%pos(1:dm)
  write(u,'(a'//csym//',3(e13.5),2x,i2)') la%o%sym,f(:),la%o%molid
enddo
close(u)

endsubroutine write_screenshot

subroutine write_state
integer                    :: bkpunit

! Elimino los subsystemas

bkpunit = find_io(30)
open( bkpunit, file=trim(adjustl(ioprefix))//'.bkp' ,form='unformatted')
write(bkpunit) time
close(bkpunit)

end subroutine write_state

!subroutine write_tcorrvel(dump)
!  !  logical,intent(in)          :: dump
!  real(dp)                    :: f=0.0_dp
!
!  if(abrir) then
!    if (u==0) then
!      u=find_io(30)
!    else
!      close(u)
!    endif        
!    open( u , file='Ptriaxial.' //trim(ioprefix)//'.dat' )
!    abrir=.false.
!  endif
!
!  !Calculo la correlacion de la primera del buffer con las demas
!  do j=1,kvel
!    f(1,j) = f(1,j) + (j-1)
!    f(2,j) = f(2,j) + (j-1)*dtframe
!    f(3,j) = f(3,j) + dot_product(vel(1,:),vel(j,:))/(3*g%nat)
!  enddo
! 
!  if(dump) then
!
!    ! cspline de la TCF
!    f(1,1:kvel) = f(1,:)
!    raux = (f(1,kvel)-f(1,1))/(kvel*100)
!    do i = 1,kvel*100
!      f(1,i) = f(1,1)+raux*(i-1)
!    enddo
!    raux = (f(2,kvel)-f(2,1))/(kvel*100)
!    do i = 1,kvel*100
!      f(2,i) = f(2,1)+raux*(i-1)
!    enddo
!    f(3,:) = cspline_v(kvel*100,io_tcvl%f(2,:),io_tcvl%f(3,:))
!
!     ! FFT de la TCF smootalista
!     do i=1,100*kvel
!       !if(units.EQ.'ev')omega=2.0_dp*Pi*(i-1)*hbarra/(dim1*dt)
!       io_irsp%f(2,i)=2.0_dp*Pi*(i-1)/(kvel*100*dt)
!     enddo
!     io_irsp%f(2,:) = 1239.84_dp/io_irsp%f(2,:) ! Salida en nm
!     caux = io_ctcv%f(3,:)
!     caux = fft(caux)
!     io_irsp%f(1,:) = real(caux)

!     !write(u,"(i0)") ubound(iof2(i)%f,2)
!     write(u,*) !pnframe,fprom,ptime
!     do j=1,ubound(f,2) ! N Filas
!       write(u,fmt=fm) f(:,j)/n
!     enddo !    

!     f=0.0_dp
!  endif
!  
!end subroutine     


!subroutine write_disrad(dump)
!  logical,save  :: abrir=.true.
!  logical,intent(in)          :: dump
!  integer,save  :: u
!  f = f + inq_disrad(disrad_i,disrad_f,g)
!end subroutine   


!    ! Escribo 
!    ptime = ptime+(time-ptime)/2.0_dp
!    pnframe = nframe-float(fprom-1)/2.0_dp
 
!ARCHIVOS PROPIEDADES

subroutine write_tmoment(op)
class(outpropa)     :: op
call inq_cm_vel(op%g) 
op%f(1:dm) = op%f(1:dm) + op%g%cm_vel(1:dm)*op%g%mass
op%f(dm+1) = op%f(dm+1) + sqrt(dot_product(op%g%cm_vel(1:dm),op%g%cm_vel(1:dm))*op%g%mass)
end subroutine

subroutine write_amoment(op)
class(outpropa)     :: op
call inq_angular_mom(op%g) 
op%f(1:3) = op%f(1:3) + op%g%ang_mom
op%f(4) = op%f(4) + sqrt(dot_product(op%g%ang_mom,op%g%ang_mom))
end subroutine

subroutine write_angvel(op)
class(outpropa)     :: op
call inq_angular_vel(op%g)                    
op%f(1:3) = op%f(1:3) + op%g%ang_vel
op%f(4) = op%f(4) + sqrt(dot_product(op%g%ang_vel,op%g%ang_vel))
end subroutine

subroutine write_virial(op)
use gems_program_types, only:b_gvirial,virial
class(outpropa)     :: op
if(.not.b_gvirial) then
  call inq_virial(op%g)
  virial(:,:)=op%g%virial(:,:)
endif
op%f(:) = op%f(:) + reshape(virial,(/dm*dm/))*ui_ev
end subroutine
                
subroutine write_girrad(op)
class(outpropa)     :: op
call group_inq_rg(op%g)
op%f(1) = op%f(1) + op%g%rg_pos
end subroutine
               
subroutine write_pressure(op)
class(outpropa)     :: op
call inq_pressure(op%g)
op%f(:) = op%f(:) + reshape(op%g%pressure,(/dm*dm/))*ui_pa*pa_bar*bar_atm
end subroutine

subroutine write_inercia(op)
class(outpropa)     :: op
call inq_inercia(op%g)
op%f = op%f + reshape(op%g%inercia,(/dm*dm/))
end subroutine

subroutine write_covariance(op)
class(outpropa)     :: op
call inq_covariance(op%g) 
op%f = op%f + reshape(op%g%covar,(/dm*dm/))
end subroutine

subroutine write_mainaxis(op)
class(outpropa)     :: op
call inq_principal_geometric_axis(op%g) 
op%f = op%f + op%g%mainaxis
end subroutine

subroutine write_ptriaxial(op)
class(outpropa)     :: op
op%f(1) = op%f(1) + inq_triaxial_param(op%g)
end subroutine

subroutine write_time(op)
class(outpropa)     :: op
op%f(1) = op%f(1) + time
end subroutine write_time

subroutine write_step(op)
class(outpropa)     :: op
op%f(1) = op%f(1) + dm_steps
end subroutine write_step

subroutine write_epot(op)
class(outpropa)     :: op
call inq_pot_energy(op%g)
op%f(1) = op%f(1) + op%g%epot*ui_ev
end subroutine write_epot

subroutine write_ecin(op)
class(outpropa)     :: op
call inq_kin_energy(op%g)
op%f(1) = op%f(1) + op%g%ekin*ui_ev
end subroutine write_ecin
 
subroutine write_absecin(op)
class(outpropa)     :: op
call inq_abskin_energy(op%g)
op%f(1) = op%f(1) + op%g%ekin*ui_ev
end subroutine write_absecin
  
subroutine write_energy(op)
class(outpropa)     :: op
call inq_kin_energy(op%g)
call inq_pot_energy(op%g)
op%f(1) = op%f(1) + op%g%epot*ui_ev
op%f(2) = op%f(2) + op%g%ekin*ui_ev
op%f(3) = op%f(1) + op%f(2) 
end subroutine write_energy

subroutine write_box(op)
class(outpropa)     :: op
op%f(:) = op%f(:) + box(:)
end subroutine write_box

subroutine write_absenergy(op)
class(outpropa)     :: op
call inq_abskin_energy(op%g)
call inq_pot_energy(op%g)
op%f(1) = op%f(1) + op%g%epot*ui_ev
op%f(2) = op%f(2) + op%g%ekin*ui_ev
op%f(3) = op%f(1) + op%f(2) 
end subroutine write_absenergy
                   
subroutine write_globalerror(op)
class(outpropa)     :: op
logical,save        :: first=.true.
real(dp)            :: eini=0.0_dp,ge=0.0_dp
real(dp)            :: aux

!Primera energia
if(first) then
  call inq_kin_energy(op%g)
  call inq_pot_energy(op%g)
  eini=op%g%epot+op%g%ekin
  first=.false.
endif

call inq_kin_energy(op%g)
call inq_pot_energy(op%g)
ge = ge+(eini-(op%g%epot+op%g%ekin))**2
aux = sqrt(ge)
op%f(1)=aux
op%f(2)=aux/(dm_steps+1._dp)**2
end subroutine write_globalerror

subroutine write_energypa(op)
class(outpropa)     :: op
call inq_kin_energy(op%g)
call inq_pot_energy(op%g)
op%f(1) = op%f(1) + op%g%epot*ui_ev/op%g%nat
op%f(2) = op%f(2) + op%g%ekin*ui_ev/op%g%nat
op%f(3) = op%f(3) + (op%g%epot+op%g%ekin)*ui_ev/op%g%nat 
end subroutine write_energypa

subroutine write_aenergy(op)
class(outpropa)     :: op
call inq_angular_energy(op%g) 
op%f(1) = op%f(1) + op%g%erot*ui_ev
op%f(2) = op%f(2) + op%g%evib*ui_ev
op%f(3) = op%f(3) + (op%g%erot+op%g%evib)*ui_ev
end subroutine

subroutine write_temp(op)
class(outpropa)     :: op
call inq_temperature(op%g)
op%f(1) = op%f(1) + op%g%temp
end subroutine

subroutine write_tempall(op)
class(outpropa)     :: op
call inq_temperature(op%g)
op%f(1) = op%f(1) + op%g%temprot
op%f(2) = op%f(2) + op%g%tempvib
op%f(3) = op%f(3) + op%g%temp
end subroutine

subroutine write_pos(of)
use gems_program_types
class(outfile)     :: of
type(atom_dclist),pointer   :: la
integer                     :: i
real(dp)                    :: f(3)

f(:)=0._dp 
 
write(of%un,*) of%g%nat
write(of%un,*) nframe,time
la => of%g%alist
do i=1,of%g%nat
  la => la%next
  !call group_inq_cmpos(of%g)
  !write(un,'(a'//csym//',3(2x,e25.12))') la%o%sym,(la%o%pos(j)-of%g%cm_pos(j),j=1,dm),(0._dp,j=dm,2)
  f(1:dm)=la%o%pos(1:dm)
  write(of%un,'(a'//csym//',3(2x,e25.12))') la%o%sym,f(:)
enddo

if(of%flush) call flush(of%un)

end subroutine

subroutine write_poscr(of)
class(outfile)            :: of
type(atom_dclist),pointer :: la
integer                   :: i,g(3)
real(dp)                  :: f(3)

f(:)=0._dp
g(:)=0._dp

write(of%un,*) of%g%nat
write(of%un,*) nframe,time
la => of%g%alist
do i=1,of%g%nat
  la => la%next
  f(1:dm)=la%o%pos(1:dm)
  g(1:dm)=la%o%boxcr(1:dm)
  write(of%un,'(a'//csym//',3(2x,e25.12),(3(2x,i0)))') la%o%sym,f(:),g(:)
enddo
 
if(of%flush) call flush(of%un)

end subroutine 
                               
subroutine write_dist(of)
use gems_program_types
class(outfile)     :: of
type(atom_dclist),pointer   :: la,lb
real(dp)                    :: rd
integer                     :: i,j

write(of%un,*) of%g%nat
write(of%un,*) nframe,time

la => of%g%alist
do i=1,of%g%nat
  la => la%next

  lb => la
  do j=i,of%g%nat-1
    lb => lb%next

    rd=rdistance2(la%o,lb%o)
    if(rd<100.0_dp) write(of%un,'(2(a'//csym//',1x),e25.12)')  la%o%sym,lb%o%sym,sqrt(rd)

  enddo
enddo

if(of%flush) call flush(of%un)

end subroutine

subroutine write_fce(of)
use gems_program_types
class(outfile)     :: of
type(atom_dclist),pointer   :: la
integer                     :: i
real(dp)                    :: f(3)

f(:)=0._dp 
 
write(of%un,*) of%g%nat
write(of%un,*) nframe,time

la => of%g%alist
do i=1,of%g%nat
  la => la%next
  f(1:dm)=la%o%force(1:dm)*ui_kcm
  write(of%un,'(a'//csym//',3(2x,e25.12))') la%o%sym,f(:)
enddo

if(of%flush) call flush(of%un)

end subroutine

subroutine write_vel(of)
use gems_program_types
class(outfile)     :: of
type(atom_dclist),pointer   :: la
integer                     :: i
real(dp)                    :: f(3)

f(:)=0._dp 
 
write(of%un,*) of%g%nat
write(of%un,*) nframe,time

la => of%g%alist
do i=1,of%g%nat
  la => la%next
  f(1:dm)=la%o%vel(1:dm)
  write(of%un,'(a'//csym//',3(2x,e25.12))') la%o%sym,f(:)
enddo

if(of%flush) call flush(of%un)

end subroutine

subroutine write_pes(of)
use gems_program_types
class(outfile)     :: of
type(atom_dclist),pointer   :: la
integer                     :: i,j
 
la => of%g%alist
do i=1,of%g%nat
  la => la%next
  write(of%un,fmt='(e25.12,'//cdm//'(x,e25.12))') (la%o%pos(j),j=1,dm),la%o%epot*ui_ev
enddo

if(of%flush) call flush(of%un)

end subroutine

subroutine write_pose(of)
use gems_program_types
class(outfile)     :: of
type(atom_dclist),pointer   :: la
integer                     :: i
real(dp)                    :: f(3)

f(:)=0._dp 
 
write(of%un,*) of%g%nat
write(of%un,*) nframe,time

la => of%g%alist
do i=1,of%g%nat
  la => la%next
  f(1:dm)=la%o%pos(1:dm)
  write(of%un,'(a'//csym//',3(2x,e25.12))') la%o%sym,f(:),la%o%epot*ui_ev
enddo

if(of%flush) call flush(of%un)

end subroutine

subroutine write_charge(of)
use gems_program_types
class(outfile)     :: of
type(atom_dclist),pointer   :: la
integer                     :: i

la => of%g%alist
do i=1,of%g%nat
  la => la%next
  write(of%un,'(a'//csym//',2x,e25.12)') la%o%sym,la%o%q
enddo                   
write(of%un,*)

if(of%flush) call flush(of%un)

end subroutine                           

subroutine write_border(of)
use gems_program_types
class(outfile)     :: of
type(atom_dclist),pointer   :: la
integer                     :: i
real(dp)                    :: r1=2.9_dp,r2=3.8_dp

 
call inq_bondorder(of%g,r1,r2)

la => of%g%alist
do i=1,of%g%nat
  la => la%next
  write(of%un,'(a'//csym//',2x,e25.12)') la%o%sym,la%o%border
enddo                   
write(of%un,*)

if(of%flush) call flush(of%un)

end subroutine                           

subroutine write_vel_rot(of)
use gems_program_types
class(outfile)     :: of
type(atom_dclist),pointer   :: la
integer                     :: i
real(dp)                    :: f(3)

f(:)=0._dp 
 
call inq_angular_energy(of%g) 

write(of%un,*) of%g%nat
write(of%un,*) nframe,time

la => of%g%alist
do i=1,of%g%nat
  la => la%next
  f(1:dm)=la%o%vel_rot(1:dm)
  write(of%un,'(a'//csym//',3(2x,e25.12))') la%o%sym,f(:)
enddo

if(of%flush) call flush(of%un)

end subroutine

subroutine write_vel_vib(of)
use gems_program_types
class(outfile)     :: of
type(atom_dclist),pointer   :: la
integer                     :: i
real(dp)                    :: f(3)

f(:)=0._dp 

call inq_angular_energy(of%g) 

write(of%un,*) of%g%nat
write(of%un,*) nframe,time

la => of%g%alist
do i=1,of%g%nat
  la => la%next
  f(1:dm)=la%o%vel_vib(1:dm)
  write(of%un,'(a'//csym//',3(2x,e25.12))') la%o%sym,f(:)
enddo

if(of%flush) call flush(of%un)

end subroutine

function polvar_outfile(var) result(g)
use gems_variables, only: polvar, polvar_find
use gems_errors, only: werr
character(*),intent(in) :: var
type(polvar),pointer    :: pv
class(outfile),pointer  :: g


call werr('Labels should start with colon `:` symbol',var(1:1)/=':')
pv=>polvar_find(var)

g=>null()
if(.not.associated(pv)) return

call werr('Variable not linked',pv%hard)
 
! Print
select type(v=>pv%val)
class is (outfile)
  g=>v
class default
  call werr('I dont know how to return that')
end select

end function polvar_outfile
     
end module gems_output

