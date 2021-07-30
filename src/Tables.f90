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

 
module gems_tables
 ! a table with a fix step, is used for potentials for example. in this module
 ! there are subroutines to operate (i.e linear_inter) with tables.

 ! Se definio dos tipos de tablas, la quiespaciada y la que no lo es. Esto es
 ! porque muchos procedimientos son muy diferentes dependiendo que condición se
 ! cumple. Sin embargo, hay procedimientos (i.e. etable_autocorr,stable_autcorr)
 ! que necesitarian de polimorfismo, dado que no hace falta distinguir entre
 ! estas tablas. Existe un paper que puede indicar como implementar polymorfismo
 ! en fortran 90:
 !   Decyk, V. 1998. How to support inheritance and run-time polymorphism in
 !   Fortran 90. Computer Physics Communications 115, no. 1 (December): 9-17.
 !   doi:10.1016/S0010-4655(98)00101-5.
 !   http://linkinghub.elsevier.com/retrieve/pii/S0010465598001015.
 ! sin embargo, esto es mucho mas facil en fortran 2003. Si decido escribir en
 ! 2003, esto se hace muy facil con class(*... no lo hago pq todavia no me
 ! animo.)

 use gems_errors
 use gems_constants,only:sp,dp,find_io
 use gems_input_parsing

 type :: stable
   !Tabla noequiespaciada
   real(dp),pointer  :: y(:),x(:) 
 endtype 
 
 type :: etable
   !Tabla equiespaciada
   real(dp)          :: ini,fin,stp
   real(dp),pointer  :: y(:)
   contains
     procedure :: init1=>etable_build_ifh
     procedure :: init2=>etable_build_fhn
     procedure :: init3=>etable_build_hni
     procedure :: init4=>etable_build_nif
     generic   :: init => init1,init2,init3,init4
     procedure :: init_histo => etable_build_histo
     procedure :: wprt => etable_wprt
 endtype 

 interface table_cspline
   module procedure etable_cspline
 end interface
   
 interface table_write
   module procedure stable_write,etable_write
 end interface
 
 interface table_read
   module procedure stable_read,etable_read
 end interface
 
 interface table_deriv_tres
   module procedure etable_deriv_tres
 end interface
 
 interface table_lreg
   module procedure etable_lreg,stable_lreg
 end interface
 
 interface table_invert
   module procedure stable_invert
 end interface
 
 interface  table_build
   module procedure etable_build, stable_build
 end interface
   
 interface  table_destroy
   module procedure etable_destroy, stable_destroy
 end interface

 interface  table_histo_add
   module procedure etable_histo_add  
 end interface

 interface  table_ceroclean
   module procedure etable_stable_ceroclean
 end interface

 interface  assignment (=)
   module procedure etable_assign, stable_assign,etable_stable_assign,stable_etable_asign
 end interface

 contains

 ! Eventos de la clase etable

  subroutine etable_build(t,npt,ini,fin,stp)
  ! Da vuelta una tabla
    integer,intent(in)           :: npt
    real(dp),intent(in),optional :: ini,fin,stp
    type(etable),intent(inout) :: t

    if(associated(t%y)) call table_destroy(t)

    if(.not.present(ini)) then
      t%ini=0.0_dp
    else
      t%ini=ini
    endif

    if(.not.present(stp)) then
      t%stp=0.0_dp
      if(present(fin)) t%stp=(fin-t%ini)/(npt-1)
    else
      t%stp=stp
    endif

    if(.not.present(fin)) then
      t%fin=0.0_dp
      if(present(stp)) t%fin=t%ini+stp*(npt-1)
    else
      t%fin=fin
    endif

    allocate(t%y(npt))

  end subroutine
 
  subroutine etable_build_ifh(t,i,f,h)
    class(etable),intent(out)    :: t
    real(dp),intent(in)          :: i,f,h

    t%ini=i
    t%fin=f
    t%stp=h
    allocate(t%y(int((f-i)/h)+1))

  end subroutine
    
  subroutine etable_build_fhn(t,f,h,n)
    class(etable),intent(out)    :: t
    integer,intent(in)           :: n
    real(dp),intent(in)          :: f,h

    t%ini=f-(n-1)*h
    t%fin=f
    t%stp=h
    allocate(t%y(n))

  end subroutine
     
  subroutine etable_build_hni(t,h,n,i)
    class(etable),intent(out)    :: t
    integer,intent(in)           :: n
    real(dp),intent(in)          :: i,h

    t%ini=i
    t%fin=i+(n-1)*h
    t%stp=h
    allocate(t%y(n))

  end subroutine
     
  subroutine etable_build_nif(t,n,i,f)
    class(etable),intent(out)    :: t
    integer,intent(in)           :: n
    real(dp),intent(in)          :: i,f

    t%ini=i
    t%fin=f
    t%stp=(f-i)/(n-1)
    allocate(t%y(n))

  end subroutine
     
  subroutine etable_destroy(t)
    type(etable),intent(out) :: t
    if(associated(t%y)) deallocate(t%y)
    t%stp=0
    t%ini=0.0_dp
    t%fin=0.0_dp
  end subroutine

  subroutine etable_assign(t2,t)
    type(etable),intent(in)    :: t  !original
    type(etable),intent(inout) :: t2  !copia

    if(associated(t2%y)) call table_destroy(t2)
    if(.not.associated(t%y)) return

    allocate(t2%y(size(t%y)))
    t2%y=t%y
    t2%stp=t%stp
    t2%ini=t%ini
    t2%fin=t%fin
  end subroutine
 
   subroutine etable_wprt(t)
   use gems_errors
    class(etable)                :: t
    real(dp)                     :: xt
    integer                      :: i

    xt=t%ini+t%stp*0.5_dp
    do i=1,size(t%y)
      call wprt(); write(printunit,*) xt,t%y(i)
      xt = xt + t%stp
    enddo

   end subroutine etable_wprt
 
   subroutine etable_write(tfile,t)
   real(dp)                     :: xt
   character(*)                 :: tfile
   type(etable)                 :: t
   integer                      :: i,u

    u = find_io(30)
    open(u,file=trim(adjustl(tfile)))

    xt=t%ini
    do i=1,size(t%y)
      write(u,*) xt,t%y(i)
      xt = xt + t%stp
    enddo

    close(u)  
   end subroutine etable_write
  
  subroutine etable_read(t,ifile,coly,colx)
  ! w must be lower case
  integer,optional,intent(inout) :: colx
  integer,intent(inout)          :: coly
  character(*)                   :: ifile
  type(etable)                   :: t
  real(dp)                       :: p1,p2
  integer                        :: i,j,n,io,ios
 
    if (present(colx)) then
      i=max(colx,coly)
      j=min(colx,coly)
      colx=j
      coly=i
    endif

    io = find_io(30)
    open(io,file=ifile)
 
    ! Cuento las lineas
    n=-1 
    ios=0
    do while(ios==0)
      read(io,iostat=ios,fmt=*)
      n=n+1
    enddo
 
    rewind(io)

    call table_build(t,n)
 
    if(present(colx)) then
      read(io,*) (p1,j=1,colx),(t%y(1),j=colx+1,coly)
      do i = 2,n-1
        read(io,*) (t%y(i),j=1,coly)
      enddo
      read(io,*) (p2,j=1,colx),(t%y(1),j=colx+1,coly)
    else
      p1=1
      do i = 1,n
        read(io,*) (t%y(i),j=1,coly)
      enddo
      p2=n
    endif

    t%ini=p1
    t%stp=(p2-p1)/(n-1)
    t%fin=p2  

    close(io)

  end subroutine
 
  subroutine charge_from_function_sp(f,t)  ! INTERFACE ASOCIATED
    ! Pasa una funcion a tabla
    real(sp),dimension(1)     :: x,dx,aux
    type(etable),intent(inout) :: t
    integer                   :: i
    interface
      function f(x)
        integer, parameter    :: sp = kind(1.0) 
        real(sp), dimension(:), intent(in) :: x
        real(sp), dimension(size(x)) :: f
      end function f
    end interface

    x=real(t%ini,sp)
    dx=real(t%stp,sp)
    do i=1,size(t%y)
      aux=f(x)
      t%y(i)=aux(1)
      x=x+dx
    enddo
  end subroutine  
 
  subroutine etable_cspline(t,t2,npt)
    ! Esta funcion devuelve un vector mas grande o no uniforme en x con un
    ! interpolamiento cubic spline. No hace falta que la tabla de entrada sea
    ! equi 
    ! Ref.: From Numath Library By Tuan Dang Trong in     
    !       Fortran 77 [BIBLI 18].                        
    !       F90 Release By J-P Moreau, Paris. 
    real(dp)                  :: rx,dx
    type(etable),intent(in)   :: t
    type(etable)              :: t2
    integer,optional          :: npt
    integer                   :: n,m
    integer                   :: l
    real(dp),allocatable      :: b(:),c(:),d(:),x(:),y(:)


    m=size(t%y)

    if(present(npt)) then
      n=npt
    else
      n=m
    endif

    allocate(b(m),c(m),d(m),x(m),y(m))

    dx = t%stp
    do l = 1,m
      x(l) = t%ini + (l-1)*dx
      y(l) = t%y(l)
    enddo
    call table_build(t2,n,ini=t%ini,fin=t%fin)

    call cn_spline (x,y,b,c,d)
    
    dx = t2%stp
    rx = t2%ini
    do l = 1,n
      t2%y(l) = seval(m,rx,x,y,b,c,d)
      rx = rx + dx
    enddo

  end subroutine etable_cspline
 
  subroutine etable_lreg(t,m,b,r)
  !Source:   Dr. David G. Simpson                                           
  !          Department of Physical Science                                 
  !          Prince George's Community College                              
  !          Largo, Maryland  20774                                         
  !Date:         January 21, 2002                                               
  !This program performs a linear regression analysis for a set
  !of data given as (x,y) pairs.  The output from  the program is the
  !slope and y-intercept of the least-squares best fit straight line through
  !the data points.      
  !Modificado por alexis paz

  ! Las rectas de regresión son las rectas que mejor se ajustan a la nube de
  ! puntos (o también llamado diagrama de dispersión) generada por una
  ! distribución binomial. Matemáticamente, son posibles dos rectas de máximo
  ! ajuste La recta de regresión de Y sobre X y la de X sobre Y.  La correlación
  ! ("r") de las rectas determinará la calidad del ajuste. Si r es cercano o
  ! igual a 1, el ajuste será bueno; si r es cercano o igual a 0, se tratará de
  ! un ajuste malo. Ambas rectas de regresión se intersecan en un punto llamado
  ! centro de gravedad de la distribución.

  real(dp)                :: sumx,sumx2,sumxy,sumy,sumy2,x,y 
  real(dp),intent(out)    :: m,b,r
  integer                 :: l
  type(etable),intent(in) :: t          ! input table data

    sumx  = 0.0_dp
    sumx2 = 0.0_dp
    sumxy = 0.0_dp
    sumy  = 0.0_dp
    sumy2 = 0.0_dp
    do l = 1,size(t%y)
      x = t%ini + (l-1)*t%stp
      y = t%y(l)
      sumx  = sumx  + x    ! compute sum of x
      sumx2 = sumx2 + x*x  ! compute sum of x**2
      sumxy = sumxy + x*y  ! compute sum of x * y
      sumy  = sumy  + y    ! compute sum of y
      sumy2 = sumy2 + y*y  ! compute sum of y**2
    enddo 

    ! compute slope
    m = (size(t%y)*sumxy  -  sumx*sumy) / (size(t%y)*sumx2 - sumx**2)
    ! compute y-intercept
    b = (sumy*sumx2 - sumx*sumxy) / (size(t%y)*sumx2 - sumx**2)
    ! compute correlation coefficient
    r = (sumxy - sumx*sumy/size(t%y)) / sqrt((sumx2 - sumx**2/size(t%y))*(sumy2 - sumy**2/size(t%y)))

  end subroutine

  subroutine etable_deriv_tres(t,dt)
   type(etable)      :: t,dt
   integer           :: i

    if(.not.associated(t%y)) return
    if(associated(dt%y)) call table_destroy(dt)
    call table_build(dt,size(t%y),ini=t%ini,stp=t%stp,fin=t%fin)

    dt%y(1) = ( t%y(2)-t%y(1) )/t%stp !euler hacia adelante

    do i = 2, size(t%y)-1
      dt%y(i) = t%y(i+1)-t%y(i-1) ! tres puntos
    enddo
    dt%y(2:size(t%y)-1) = dt%y(2:size(t%y)-1)/(2.0_dp*t%stp)

    dt%y(size(t%y)) = ( t%y(size(t%y))-t%y(size(t%y)-1) )/t%stp  !euler hacia atras

  end subroutine etable_deriv_tres

  subroutine etable_deriv_cinco(t,dt)
   type(etable)       :: t,dt
   integer           :: i

    if(.not.associated(t%y)) return
    if(associated(dt%y)) call table_destroy(dt)
    call table_build(dt,size(t%y),ini=t%ini,stp=t%stp,fin=t%fin)

    dt%y(1) = ( t%y(2)-t%y(1) )/t%stp !euler hacia adelante
    dt%y(2) = ( t%y(3)-t%y(1) )/(2.0_dp*t%stp) ! tres puntos

    do i = 3, size(t%y)-2
      dt%y(i) = (t%y(i-2)-8.0_dp*t%y(i-1)+8.0_dp*t%y(i+1)-t%y(i+2)) ! Cinco puntos
    enddo

    dt%y(3:size(t%y)-2) = dt%y(3:size(t%y)-2)/(12.0_dp*t%stp)

    dt%y(size(t%y)-1) = ( t%y(size(t%y))-t%y(size(t%y)-2) )/(2.0_dp*t%stp)! tres puntos
    dt%y(size(t%y)) = ( t%y(size(t%y))-t%y(size(t%y)-1) )/t%stp  !euler hacia atras

  end subroutine etable_deriv_cinco
 
  subroutine etable_histo_add(t,rx)
   ! Si el valor esta fuera del rango del histrograma, simplemente no se cuenta
   ! La decision de que no se cuente es porque hay distribuciones con dominio
   ! todos los reales... y por ende puede ocurrir un evento en infinito, pero
   ! hay que tener cuidado.
   type(etable)      :: t
   real(dp)          :: rx,dx
   integer           :: ix
    dx = (rx-t%ini)/t%stp        !dx is the real index for a x continius grid (ideal).
    ix = idnint(dx)+1            !ix is the integer index (knowed, and lower close) of the rho discrete grid.
    if (ix>size(t%y).or.ix<0) return ! Si el valor esta fuera del rango, no se cuenta
    t%y(ix) = t%y(ix)+1          !sumo al histograma

  end subroutine etable_histo_add

  subroutine etable_build_histo(et,r,h)
   class(etable)       :: et
   real(dp),intent(in) :: r(:),h
   integer             :: ix,i

   call werr('bin histogram lower than cero',h<0)
   call et%init(minval(r),maxval(r),h)
   et%y=0.0_dp

   !make the histogram
   do i = 1, size(r)
     ix = int((r(i)-et%ini)/et%stp)+1
     et%y(ix) = et%y(ix) + 1
   enddo
   
   !volfac = 4.0_dp*pi/dm
   
   !!normalize the histogram
   !do i = 1,npoints
   !  inq_disrad(1,i) = dr*(i-0.5_dp)
   !  inq_disrad(2,i) = et%y(i)/(volfac*(i*dr)**dm)
   !enddo
   !
   !inq_disrad(2,:) = inq_disrad(2,:) / maxval(inq_disrad)
 
  end subroutine etable_build_histo
               
  function linear_inter(t,rx)
   ! No se puede interpolar tablas no equipaciadas
   ! Hay que convertirla antes (ej, cspline).
   ! Si esta fuera del rango de la tabla satara un segmentation fault
   type(etable)      :: t
   real(dp)          :: rx,dx,linear_inter
   integer           :: ix

    dx = (rx-t%ini)/t%stp        !dx is the real index for a x continius grid (ideal).
    ix = int(dx)                 !ix is the integer index (knowed, and lower close) of the rho discrete grid.
    dx = dx-ix                   !now dx is the delta index (0 to 1) and is used to interpolate.
    ix = ix+1                    !the index of y begin whit 1

    if(ix>size(t%y)) then
      print *, ' warning: interpolation out of range'
      linear_inter = t%y(size(t%y))
    else
      linear_inter = t%y(ix) + (t%y(ix+1)-t%y(ix))*dx 
    endif
  end function linear_inter 

subroutine etable_autocorr(a,t)
  ! Suministrar valores a a la subrutina. Estos se guardan en el vector buffer.
  ! cuando se han suministrado tantos valores como tamaño tenga la tabla t, se
  ! calcula la tcf y se la guarda en la tabla t. Logicamente el valor
  ! estadistico de la tcf a ventanas de tiempo grande es poco. Asique se puede
  ! seguir suministrando valores de a y la funcion de autocrrelacion mejora.
  real(dp),intent(in)          :: a
  type(etable),intent(inout)   :: t     
  integer,save                 :: j=0,n=0
  integer                      :: i
  real(dp),allocatable,save    :: buffer(:)


  !Si esta vacio el buffer
  if(j<size(t%y)) then

    !lo allocateo y relleno
    if(n==0) then
      n=size(t%y)
      allocate(buffer(n))
    endif
    buffer(j)=a
    j=j+1

  else

    !calculo y lo corro un paso
    do i=1,size(t%y)
      t%y(i) = t%y(i) + buffer(1)*buffer(i)
    enddo

    !Adelanto el buffer
    buffer = cshift(buffer,shift=1) 
    buffer(size(t%y)) = a

  endif

end subroutine     
 
!subroutine etable_fpplot(signal,var,err)
!  ! Realiza el grafico de Flyvbjerg and Petersen: desviacion estandar en funcion
!  ! de tamaño de bloque estadistico. El plato en el grafico indica que el tamaño
!  ! de bloque consigue perder la correlacion. Usa el algoritmo DDDA
!  use ddda
!  type(etable),intent(in)        :: signal
!  type(etable),intent(out)       :: var,err
!  real(dp)                       :: v,e
!  integer                        :: i
!  class(statistic_dclist),pointer :: ls
!  type(decorrelation)            :: d  ! El objeto decorrelation del DDDA
!  
!  ! No quiero nada de aca
!  call table_destroy(var)
!  call table_destroy(err)
!  call d%init()
!
!  do i = 1,size(signal%y)
!    call d%attach(signal%y(i))
!  enddo
!
!  ! El siguiente do hace log_2(N)-1 lineas, donde N es el numero de puntos
!  ! agregados en el d. El menos uno es porque para calcula var se divide en
!  ! N-1 y esto provocaria un NAN 
!
!  call table_build(var,d%size-1)
!  call table_build(err,d%size-1)
!
!  ls =>d%blocks
!  call wstd(); write(logunit,fmt='(a,e25.12)')  '# la media es ',ls%o%suma/ls%o%nsamples
!  do i=1,d%size-1 
!    v = ls%o%vars(e)
!    var%y(i)=sqrt(v)
!    err%y(i)=sqrt(e)
!    ls => ls%next
!  enddo
!
!end subroutine     
 
 ! Eventos de la clase stable

  subroutine stable_build(t,npt)
  ! Da vuelta una tabla
    integer,intent(in)         :: npt
    type(stable),intent(inout) :: t

    if(associated(t%x)) call table_destroy(t)

    allocate(t%y(npt))
    allocate(t%x(npt))
  end subroutine
 
  subroutine stable_destroy(t)
    type(stable),intent(out) :: t

    deallocate(t%y)
    deallocate(t%x)
  end subroutine
 
  subroutine stable_assign(t2,t)
    type(stable),intent(in)    :: t  !original
    type(stable),intent(inout) :: t2  !copia

    if(associated(t2%y)) call table_destroy(t2)
    if(.not.associated(t%y)) return

    allocate(t2%x(size(t%x)))
    allocate(t2%y(size(t%y)))
    t2%y=t%y
    t2%x=t%x
  end subroutine
 
  subroutine stable_monovalued(t,e)
    !Si la tabla presenta valores de x muy cercanos (en el sentido que no
    !difieren mas de un valor e) y consecutivos , combina los dos en uno,
    !promediando los valores de y
    type(stable),intent(inout) :: t 
    real(dp),intent(in)        :: e
    integer                    :: n,m,i,j
    real(dp)                   :: x(size(t%y)),y(size(t%y)),aux,auy

    if(.not.associated(t%x)) return

    !n sera el numero total de valores x
    !j sera un indice que vaya por los valores x
    !i sera un indice que vaya por los valores x distintos
    n=size(t%y)
    i=1
    j=1
    do while (j<n)
      m=1
      auy=t%y(j)
      aux=t%x(j)
      do while (dabs(t%x(j)-t%x(j+1))<e) 
         m=m+1
         j=j+1
         auy=auy+t%y(j) ! Sumo para promediar
         aux=aux+t%x(j) ! Sumo para promediar
         if (j==n) exit
      enddo
      y(i)=auy/m  ! promedio 
      x(i)=aux/m  ! promedio 
      i=i+1
      j=j+1
    enddo

    ! En el caso de que el valor nesimo ya haya sido promediado, j=n+1, de lo
    ! contrario j=n
    if (j==n) then
      y(i)=t%y(j)
      x(i)=t%x(j)
      n=i
    else
      n=i-1
    endif

    if (n/=size(t%y)) then

      call table_destroy(t)
      call table_build(t,n)
      
      do i =1,n
        t%y(i)=y(i)
        t%x(i)=x(i)
      enddo
       
    endif

  end subroutine stable_monovalued
 
  subroutine stable_invert(t)
    real(dp)                  :: aux
    type(stable),intent(inout) :: t
    integer                   :: i

    do i=1,size(t%y)
      aux=t%x(i)
      t%x(i)=t%y(i)
      t%y(i)=aux
    enddo
  end subroutine

   subroutine stable_write(tfile,t)
   character(*)                 :: tfile
   type(stable)                 :: t
   integer                      :: i,u
 
    u = find_io(30)
    open(u,file=trim(adjustl(tfile)))

    do i=1,size(t%y)
      write(u,*) t%x(i),t%y(i)
    enddo 

    close(u)  

  end subroutine stable_write

  subroutine stable_read(t,ifile,coly,colx)
  ! w must be lower case
  integer,intent(inout) :: coly
  integer,optional,intent(inout) :: colx
  character(*)          :: ifile
  type(stable)          :: t
  integer               :: i,j,n,io,ios
 
    if (present(colx)) then
      i=max(colx,coly)
      j=min(colx,coly)
      colx=j
      coly=i
    endif
 
    io = find_io(30)
    open(io,file=ifile)
 
    ! Cuento las lineas
    n=-1
    ios=0
    do while(ios==0)
      read(io,iostat=ios,fmt=*)
      n=n+1
    enddo
 
    rewind(io)

    call table_build(t,n)
 
    if (present(colx)) then
    do i = 1,n
      read(io,*) (t%x(i),j=1,colx),(t%y(i),j=colx+1,coly)
    enddo
    else
      do i = 1,n
        read(io,*) (t%y(i),j=1,coly)
        t%x(i)=i
      enddo
    endif
    close(io)

  end subroutine
 
  subroutine stable_lreg(t,m,b,r)
  ! ver etable_lreg para los comentarios
  real(dp)                :: sumx,sumx2,sumxy,sumy,sumy2,x,y 
  real(dp),intent(out)    :: m,b,r
  integer                 :: l
  type(stable),intent(in) :: t

    sumx  = 0.0_dp
    sumx2 = 0.0_dp
    sumxy = 0.0_dp
    sumy  = 0.0_dp
    sumy2 = 0.0_dp
    do l = 1,size(t%y)
      x = t%x(l)
      y = t%y(l)
      sumx  = sumx  + x    ! compute sum of x
      sumx2 = sumx2 + x*x  ! compute sum of x**2
      sumxy = sumxy + x*y  ! compute sum of x * y
      sumy  = sumy  + y    ! compute sum of y
      sumy2 = sumy2 + y*y  ! compute sum of y**2
    enddo 

    ! compute slope
    m = (size(t%y)*sumxy  -  sumx*sumy) / (size(t%y)*sumx2 - sumx**2)
    ! compute y-intercept
    b = (sumy*sumx2 - sumx*sumxy) / (size(t%y)*sumx2 - sumx**2)
    ! compute correlation coefficient
    r = (sumxy - sumx*sumy/size(t%y)) / sqrt((sumx2 - sumx**2/size(t%y))*(sumy2 - sumy**2/size(t%y)))

  end subroutine
 
 ! Eventos de comunicacion entre las clases stable y etable

  subroutine etable_stable_assign(st,et)
    type(etable),intent(in)    :: et !original
    type(stable),intent(inout) :: st  !copia
    real(dp)                   :: dx,x
    integer                    :: i

    if(.not.associated(et%y)) return

    if(associated(st%y)) call table_destroy(st)
    call table_build(st,size(et%y))

    x=et%ini
    dx=et%stp
    do i =1,size(et%y)
      st%x(i)=x
      st%y(i)=et%y(i)
      x=x+dx
    enddo
  end subroutine
 
  subroutine stable_etable_asign(et,st)
    ! Esta funcion devuelve un vector mas grande o no uniforme en x con un
    ! interpolamiento cubic spline. No hace falta que la tabla de entrada sea
    ! equi 
    ! Ref.: From Numath Library By Tuan Dang Trong in     
    !       Fortran 77 [BIBLI 18].                        
    !       F90 Release By J-P Moreau, Paris. 
    real(dp)                  :: rx,dx
    type(stable),intent(in)   :: st
    type(etable),intent(inout):: et
    integer                   :: m,l
    real(dp),allocatable      :: b(:),c(:),d(:),x(:),y(:)

    if(.not.associated(st%y)) return
    if(associated(et%y)) call table_destroy(et)
    m=size(st%y)
    call table_build(et,m,ini=st%x(1),fin=st%x(m))

    allocate(b(m),c(m),d(m),x(m),y(m))

    x = st%x
    y = st%y

    ! Ojo, se necsita que los valores de x difieran un cierto valor para que el
    ! cspline no devuelva nan
    call cn_spline (x,y,b,c,d)
    
    dx = et%stp
    rx = et%ini
    do l = 1,m
      et%y(l) = seval(m,rx,x,y,b,c,d)
      rx = rx + dx
    enddo

  end subroutine stable_etable_asign

  subroutine etable_stable_ceroclean(et,st)
   ! Borra los puntos que sean cero de la tabla. Util para tomar log de
   ! histogramas
   type(etable),intent(in)   :: et
   type(stable),intent(out)  :: st
   real(dp)                  :: x(size(et%y)),y(size(et%y))
   integer                   :: i,j
   
   j=0
   do i = 1, size(et%y)
     if (et%y(i)/=0) then
       j=j+1
       x(j) = et%ini + (i-1)*et%stp
       y(j) = et%y(i)
     endif
   enddo
   
   call table_build(st,j)
   st%x=x
   st%y=y
    
  end subroutine etable_stable_ceroclean

! Auxiliares para cspline
  
  function cspline_v(nx,x,y) result(v)
    integer,intent(in)        :: nx
    integer                   :: l,n
    real(dp)                  :: rx,dx
    real(dp)                  :: x(:),y(:),ini
    real(dp),allocatable      :: b(:),c(:),d(:),v(:)


    n=ubound(x,1)
    if(n/=ubound(y,1)) then
      call wstd(); write(logunit,*) 'Cubic spline error'
      stop
    endif

    allocate(b(n))
    allocate(c(n))
    allocate(d(n))

    ini = x(1)
    dx = (x(n)-x(1))/(nx-1)

    call cn_spline (x,y,b,c,d)

    allocate(v(nx))
    v(1)  = y(1)
    v(nx) = y(n)
    do l = 1,nx
      rx = ini + (l-1)*dx
      v(l) = seval (n,rx,x,y,b,c,d)
    enddo

  end function cspline_v
 
  function ysum_t(t1,t2) result(tr)
    type(etable),intent(in)    :: t1,t2
    type(etable)               :: tr
    integer                    :: i

    if(size(t1%y)/=size(t2%y).or.t1%stp/=t2%stp.or.t1%ini/=t2%ini) then
      call wstd(); write(logunit,*) 'Las tablas no son ysumables'
      stop
    endif

    tr%ini=t1%ini
    tr%stp=t1%stp
    tr%fin=t1%fin

    allocate(tr%y(size(t1%y)))
    do i = 1,size(tr%y)
      tr%y(i) = t1%y(i) + t2%y(i)
    enddo

  end function
 
  function seval (n,u,x,y,b,c,d)
!------------------------------------------------------------------------
!     evaluate a cubic spline interpolation of a discreet function f(x),
!     given in n points x(i), y(i). the b, c and d coefficients defining
!     the best cubic spline for the given points, are calculated before
!     by the spline subroutine.
!
!     inputs:
!     n       number of points of curve y = f(x)
!     u       abscissa of point to be interpolated
!     x,y     tables of dimension n, storing the coordinates
!             of curve f(x)
!     b,c,d   tables storing the coefficients defining the
!             cubic spline
!
!     outputs:
!     seval   interpolated value
!             = y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))
!             with dx = u-x(i), u between x(i) and x(i+1)
!
!     reference :
!     forsythe,g.e. (1977) computer methods for mathematical
!     computations. prentice-hall,inc.
!------------------------------------------------------------------------
      integer,intent(in)         :: n
      real(dp) b(n),c(n),d(n),x(n),y(n),u,dx,seval
      integer           :: i=1,j,k

!     binary search

      if (i>=n) i = 1
      if (u<x(i)) go to 10
      if (u<=x(i+1)) go to 30
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if (u<x(k)) j = k
      if (u>=x(k)) i = k
      if (j>i+1) go to 20

!     spline evaluation

   30 dx = u-x(i)
      seval = y(i)+dx*(b(i)+dx*(c(i)+dx*d(i)))

      end function 

      subroutine cn_spline (x,y,b,c,d)
!     Calculates the coefficients b,c,d of a cubic spline to best approximate a
!     discreet fonction given by n points
!
!     inputs:
!     x,y     vectors of dimension n, storing the coordinates
!             of function f(x)
!
!     outputs:
!     a,b,c   vectors of dimension n, storing the coefficients
!             of the cubic spline
!
!     Ref:
!         forsythe,g.e. (1977) computer methods for mathematical
!         computations. prentice-hall,inc.

      implicit real(dp)(a-h,o-z)
      integer              :: n
      integer              :: i,nm1,l
      real(dp)             :: b(:),c(:),d(:),x(:),y(:),t

      n=size(x)
      nm1 = n-1

      if (n.lt.2) return
      if (n.lt.3) go to 50

!     build the tridiagonal system
!     b (diagonal), d (upperdiagonal) , c (second member)
      d(1) = x(2)-x(1)
      c(2) = (y(2)-y(1))/d(1)
      do i = 2,nm1
        d(i) = x(i+1)-x(i)
        b(i) = 2.0_dp*(d(i-1)+d(i))
        c(i+1) = (y(i+1)-y(i))/d(i)
        c(i) = c(i+1)-c(i)
      enddo

!     conditions at limits
!     third derivatives obtained by divided differences

      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.0_dp
      c(n) = 0.0_dp
      if (n.eq.3) go to 15
      c(1) = c(3)/(x(4)-x(2))-c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2))-c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)*d(1)/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))

!     forward elimination

   15 do 20 i = 2,n
      t = d(i-1)/b(i-1)
      b(i) = b(i)-t*d(i-1)
      c(i) = c(i)-t*c(i-1)
   20 continue

!     back substitution

      c(n) = c(n)/b(n)
      do 30 l = 1,nm1
      i = n-l
      c(i) = (c(i)-d(i)*c(i+1))/b(i)
   30 continue

!     coefficients of 3rd degree polynomial

      b(n) = (y(n)-y(nm1))/d(nm1)+d(nm1)*(c(nm1)+2.0_dp*c(n))
      do 40 i = 1,nm1
      b(i) = (y(i+1)-y(i))/d(i)-d(i)*(c(i+1)+2.0_dp*c(i))
      d(i) = (c(i+1)-c(i))/d(i)
      c(i) = 3.0_dp*c(i)
   40 continue
      c(n) = 3.0_dp*c(n)
      d(n) = d(nm1)
      return

!     cas n = 2

   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.0_dp
      d(1) = 0.0_dp
      b(2) = b(1)
      c(2) = 0.0_dp
      d(2) = 0.0_dp

      return
      end subroutine cn_spline

! Autocorrelacion





!  subroutine denprob_table(f,t,m)
!  ! Devuelve una tabla con la inversa de la integral acumulada de f
!  ! Esto es usado para un random sampling
!    type(table),intent(inout) :: t
!    integer                :: i
!    integer,intent(in)     :: m
!    real(dp)               :: x
!    interface
!      function f(x)
!        real(dp), intent(in) :: x
!        real(dp)             :: f 
!      end function f
!    end interface 
! 
!    call table_equi_build(t,m,ini=0.0_dp,fin=1.0_dp)
!    t%y(1)=1/(f(0.0_dp)*m)
!      print *, 0,f(0.0_dp)
!      stop
!    do i = 1,m-1
!      x=t%y(i)
!      print *, x,f(x)
!      t%y(i+1) = 1/(m*f(x))+t%y(i)
!    enddo
!
!  end subroutine 
 


! Regresiones

!  subroutine read_table(tabfile,t)
!   ! the table dont cant to have more than 1000 rows
!   namelist/tab/taby,tabini,tabfin,tabpnt,tabpas
!   real(dp)                     :: tabini,tabfin,tabpas
!   real(dp),dimension(1000)     :: taby
!   integer                      :: tabpnt
!   character(*)                 :: tabfile
!   type(table)                  :: t
!
!    open(910,file=tabfile)
!    read(910,nml=tab)
!
!    tabpas=0.0_dp
!
!    if (tabpas==0.0_dp) tabpas=(tabfin-tabini)/tabpnt
!
!    t%tabfile=tabfile
!    t%x(1)=tabini
!    t%x(2)=tabpas
!    t%x(3)=tabfin
!    allocate(t%y(tabpnt))
!    t%y=taby(1:tabpnt)
!
!    close(910)
!
!  end subroutine read_table
 

end module gems_tables
