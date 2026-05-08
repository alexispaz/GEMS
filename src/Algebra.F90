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

 
module gems_algebra
  ! a parameter dp is used to set te precision of the module
  use gems_constants, only: pi,dp,dm,sqrt2
  use gems_random,only: ranu
  use gems_errors

  ! convention of matrix fortran to matrix algebraic expresion???? 
  !             | a(1,1) a(2,1) a(3,1) |
  !a(3,2)=a3x2= | a(1,2) a(2,2) a(3,2) |
  ! me parece que en fortran los indices estan transpuestros???- las matrices rotacion estan correctas segun esto

  implicit none

  public
  private :: pi,dp ! constants

  private :: partition
        
! Vectors of reals and integers are always needed

#define _NODE real_v
#define _TYPE real(dp)
#include "vector_header.inc"
              
#define _NODE integer_v
#define _TYPE integer
#include "vector_header.inc"
       
  interface dabs
    module procedure norm_vect
  end interface

 ! interface cross_product
 !   module procedure cross_product_1,cross_product_2,cross_product_3
 ! end interface

 interface rotmatrix
   module procedure rotmatrix_2x2,rotmatrix_3x3
 end interface

 ! This is a list for the parameters of the tapper function
 type(real_v)                  :: swc
       
 contains
 
#define _NODE real_v
#define _TYPE real(dp)
#include "vector_body.inc"
               
#define _NODE integer_v
#define _TYPE integer
#include "vector_body.inc"
                    
 ! Utilidades programativas

  function decode_logicalvector(x) result(y)
    ! Considerando el vector logico `x` como un numero binario, devuelve el
    ! correspondiente un numero entero.
    integer                      :: y
    logical,intent(in)           :: x(:)
    integer                      :: i,base
    
    y=0
    base=1
    do i = 1,size(x)
     if (x(i)) y=y+base
     base=2*base
    enddo
    
  end function

  subroutine encode_logicalvector(x,yin)
    ! Convierte `yin` en un numero binario que luego asigna al vector logico `x`
    integer,intent(in)           :: yin
    logical                      :: x(:)
    integer                      :: i,n

    n=size(x)
    call werr('No space in logical vector to store this number',n<bit_size(yin))

    do i=1,n
        x(i) = btest(yin,i-1) ! bit intrinsics index from 0
    end do
    
  end subroutine

recursive subroutine sort_int(a)
  ! Ordena 2 enteros de menor a mayor
  integer,intent(inout)  :: a(:)
  integer :: iq

  if (size(a) == 1) return

  call partition_integer(a, iq)
  call sort_int(a(:iq-1))
  call sort_int(a(iq:))

end subroutine sort_int
                              
          

      
 ! Analisis
  
  function inverte_bifuncrr(y,xx1,xx2,func) result(x)
    ! Return x for a given value of a monotonic function f(x) 
    real(dp)                     :: x
    real(dp),intent(in)          :: y,xx1,xx2
    real(dp),parameter           :: tol=1.0e-8_dp
    real(dp)                     :: y1,y2,x1,x2,ym
    integer                      :: i
    interface 
      function func1(x) 
        import              :: dp
        real(dp)            :: func1
        real(dp),intent(in) :: x
      end function
    end interface
    procedure(func1)             :: func
 
    x=xx1
    y1=func(x)
    if(abs(y-y1)<tol) return

    x=xx2
    y2=func(x)
    if(abs(y-y2)<tol) return

    call werr('x value is beyond the limit',y>max(y2,y1))
    call werr('x value is beyond the limit',y<min(y2,y1))

    if(y1<y2) then
      x1=min(xx1,xx2)
      x2=max(xx1,xx2)
    else
      x1=max(xx1,xx2)
      x2=min(xx1,xx2)
    endif

    do i=1,1000000
      x=x1+(x2-x1)/2.0_dp
      ym=func(x)
      if(abs(y-ym)<tol) return

      if(y>ym) then
        x1=x
      else
        x2=x
      endif

    enddo

    call werr('bisection didnt converge',.true.)
    
  end function
  
  
 ! Statistical functions
  
  function cdf_norm(x,sigma,mu) result(f)
    use gems_constants, only:oneos2
    ! Cumulative distribution function of the normal distribution
    ! f=0.5+0.5*erf(x-mu/sqrt(2)/sigma)
    real(dp)                     :: f
    real(dp),intent(in)          :: x
    real(dp),optional,intent(in) :: sigma,mu
    real(dp)                     :: s,m

    s=1.0_dp
    m=0.0_dp
    if(present(sigma)) s=sigma
    if(present(mu)) m=mu
    !f=cdf_snorm((x-m)/s)
    !f=0.5_dp*(1.0_dp+erf((x-m)/(s*sqrt2)))

    !Use erfc to avoid the sum involve in the form with erf
    f=0.5_dp*erfc((m-x)/(s*sqrt2))
  end function
   
  function inv_cdf_norm(x,sigma,mu) result(f)
    ! Inverse of the cumulative distribution function of the normal distribution
    ! f=sqrt(2)*inv_erf(2*x-1)*sigma+mu
    real(dp) :: f
    real(dp),intent(in) :: x
    real(dp),optional,intent(in) :: sigma,mu
    real(dp)               :: s,m

    s=1.0_dp
    m=0.0_dp
    if(present(sigma)) s=sigma
    if(present(mu)) m=mu
    f=inv_cdf_snorm(x)*s+m
  end function
                          
  function cdf_snorm(x) result(f)
    use gems_constants, only:oneos2
    ! Cumulative distribution function of the standar normal distribution

    !   1 ++--+------+-------+------+-------+------+------+---******************
    !     |   +      +       +      +       +      +   ******** norm(x) ****** |
    !     |                                         ****                       |
    !     |                                        **                          |
    ! 0.8 ++                                     **                           ++
    !     |                                     **                             |
    !     |                                    **                              |
    !     |                                   **                               |
    ! 0.6 ++                                **                                ++
    !     |                                **                                  |
    !     |                               **                                   |
    ! 0.4 ++                             **                                   ++
    !     |                             **                                     |
    !     |                            **                                      |
    !     |                           **                                       |
    ! 0.2 ++                        ***                                       ++
    !     |                       ***                                          |
    !     |                    ****                                            |
    !   0 **********************    +       +      +      +       +      +    ++
    !     +---+------+-------+------+-------+------+------+-------+------+-----+
    !        -4     -3      -2     -1       0      1      2       3      4

    real(dp) :: f
    real(dp),intent(in) :: x
    !f=0.5_dp*(1.0_dp+erf(x*oneos2))
    !Use erfc to avoid the sum involve in the form with erf
    f=0.5_dp*erfc(-x*oneos2)
  end function
   
  function inv_cdf_snorm(p) result(x)
    ! Inverse of the cdf_snorm function

    !      +---------+---------------+----------------*----------------+--------+                                                                                                                   
    !    2 ++        +               +               **       invnorm(x) ******++                                                                                                                   
    !      |                                         *                          |                                                                                                                   
    !  1.5 ++                                       **                         ++                                                                                                                   
    !      |                                       **                           |                                                                                                                   
    !      |                                     ***                            |                                                                                                                   
    !    1 ++                                  ***                             ++                                                                                                                   
    !      |                                ****                                |                                                                                                                   
    !  0.5 ++                            ****                                  ++                                                                                                                   
    !      |                          ****                                      |                                                                                                                   
    !    0 ++                      ****                                        ++                                                                                                                   
    !      |                    ****                                            |                                                                                                                   
    ! -0.5 ++                ****                                              ++                                                                                                                   
    !      |               ***                                                  |
    !   -1 ++            ***                                                   ++
    !      |            **                                                      |
    ! -1.5 ++         **                                                       ++
    !      |          *                                                         |
    !   -2 ++        **              +                +                +       ++
    !      +---------*---------------+----------------+----------------+--------+
    !               0              0.5               1               1.5

    !The relative error of the approximation has absolute value less than 1.15e-9
    !One iteration of Halley’s rational method (third order) gives full machine
    !precision. The idea with a pseudocode by Peter John Acklam:
    !e <- 0.5 * erfc(-x/sqrt(2)) - p    
    !u <- e * sqrt(2*pi) * exp(x^2/2)
    !x <- x - u/(1 + x*u/2)
    
    !Here Ren-Raw Chen (associate professor at Rutgers University in New Brunswick,
    !New Jersey, USA), wrote a Fortran implementation based on John Herrero’s VB
    !implementation. And now I (alexis) adapt this code to fortran 90.
       
    real(dp),parameter :: a1=-39.6968302866538_dp, a2=220.946098424521_dp, a3=-275.928510446969_dp   
    real(dp),parameter :: a4=138.357751867269_dp, a5=-30.6647980661472_dp, a6=2.50662827745924_dp    
    real(dp),parameter :: b1=-54.4760987982241_dp, b2=161.585836858041_dp, b3=-155.698979859887_dp   
    real(dp),parameter :: b4=66.8013118877197_dp, b5=-13.2806815528857_dp
    real(dp),parameter :: c1=-0.00778489400243029_dp, c2=-0.322396458041136_dp, c3=-2.40075827716184_dp
    real(dp),parameter :: c4=-2.54973253934373_dp, c5=4.37466414146497_dp, c6=2.93816398269878_dp    
    real(dp),parameter :: d1=0.00778469570904146_dp, d2=0.32246712907004_dp, d3=2.445134137143_dp
    real(dp),parameter :: d4=3.75440866190742_dp
    real(dp),parameter :: p_low=0.02425_dp, p_high=1.0_dp-p_low
    real(dp) :: x
    real(dp),intent(in) :: p
    real(dp) :: q,r
    
    call werr('NAN on erfi',p>1.0_dp.or.p<0.0_dp)
    
    if(p<p_low) then
      q=sqrt(-2.0_dp*log(p))
      x=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    elseif(p<p_high) then
      q=p-0.5_dp
      r=q*q
      x=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
    else
      q=sqrt(-2.0_dp*log(1.0_dp-p))
      x=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    endif
    
  end function
 
function inv_erf(y) result(x)
    ! FIXME! LA inv_cdf_snorm DA MEJOR QUE inv_erf(2*y-1)
    ! Inverse of the erf function
    !      +---------+---------------+----------------*----------------+--------+                                                                                                                   
    !    2 ++        +               +               **        inverf(x) ******++                                                                                                                   
    !      |                                         *                          |                                                                                                                   
    !  1.5 ++                                       **                         ++                                                                                                                   
    !      |                                       **                           |                                                                                                                   
    !      |                                     ***                            |                                                                                                                   
    !    1 ++                                  ***                             ++                                                                                                                   
    !      |                                ****                                |                                                                                                                   
    !  0.5 ++                            ****                                  ++                                                                                                                   
    !      |                          ****                                      |                                                                                                                   
    !    0 ++                      ****                                        ++                                                                                                                   
    !      |                    ****                                            |                                                                                                                   
    ! -0.5 ++                ****                                              ++                                                                                                                   
    !      |               ***                                                  |
    !   -1 ++            ***                                                   ++
    !      |            **                                                      |
    ! -1.5 ++         **                                                       ++
    !      |          *                                                         |
    !   -2 ++        **              +                +                +       ++
    !      +---------*---------------+----------------+----------------+--------+
    !               -1               0                1                2

      
    !rational approx coefficients
    real(dp),parameter  :: a(4) = [ 0.886226899_dp, -1.645349621_dp,  0.914624893_dp, -0.140543331_dp ]
    real(dp),parameter  :: b(4) = [-2.118377725_dp,  1.442710462_dp, -0.329097515_dp,  0.012229801_dp ]
    real(dp),parameter  :: c(4) = [-1.970840454_dp, -1.624906493_dp,  3.429567803_dp,  1.641345311_dp ]
    real(dp),parameter  :: d(2) = [ 3.543889200_dp,  1.637067800_dp ]

    real(dp),parameter  :: y0=0.7_dp
    real(dp),intent(in) :: y
    real(dp) :: x, z

      call werr('NAN on erfi',abs(y)>=1.0_dp) 
    
      if(y<-y0) then
        z = sqrt(-log((1.0_dp+y)/2.0_dp))
        x = -(((c(4)*z+c(3))*z+c(2))*z+c(1))/((d(2)*z+d(1))*z+1.0_dp)
      else
        if (y<y0) then
          z = y*y
          x = y*(((a(4)*z+a(3))*z+a(2))*z+a(1))/((((b(4)*z+b(4))*z+b(2))*z+b(1))*z+1.0_dp)
        else
          z = sqrt(-log((1.0_dp-y)/2.0_dp))
          x = (((c(4)*z+c(3))*z+c(2))*z+c(1))/((d(2)*z+d(1))*z+1.0_dp)
        endif
        !polish x to full accuracy
        x = x - (erf(x)-y) / (2.0_dp/sqrt(pi) * exp(-x*x))
        x = x - (erf(x)-y) / (2.0_dp/sqrt(pi) * exp(-x*x))
      endif
    
  end function 

 ! function outerprod_r(a,b)
 !   real(sp), dimension(:), intent(in) :: a,b
 !   real(sp), dimension(size(a),size(b)) :: outerprod_r
 !   outerprod_r = spread(a,dim=2,ncopies=size(b)) * spread(b,dim=1,ncopies=size(a))
 ! end function outerprod_r
          
! Smooth step functions

function fcut(r,r1,r2)
use gems_constants, only: pi
real(dp)       :: r,r1,r2
real(dp)       :: fcut,aux

if(r<r1) then
  fcut=1.0_dp
elseif(r>r2) then
  fcut=0.0_dp
else
  aux=pi*(r-(r1+r2)*0.5_dp)/(r2-r1)
  fcut =0.5_dp-0.5_dp*sin(aux)
endif      

end function 

function fcut_dfcut(r,r1,r2,dfcut) result(fcut)
use gems_constants, only: pi,pio2
real(dp)       :: r,r1,r2
real(dp)       :: fcut,dfcut,aux

if(r<r1) then
  fcut=1.0_dp
  dfcut=0.0_dp
elseif(r>r2) then
  fcut=0.0_dp
  dfcut=0.0_dp
else
  aux=pi*(r-(r1+r2)*0.5_dp)/(r2-r1)
  fcut =0.5_dp-0.5_dp*sin(aux)
  !dfcut=-0.5_dp*cos(aux)*pio2*r/(r2-r1)
  dfcut=-cos(aux)*pio2/(r2-r1)
endif      

end function 

function fcut2(r,r1,r2) result(fcut)
! The flipped fcut
use gems_constants, only: pi
real(dp)       :: r,r1,r2
real(dp)       :: fcut,aux

if(r<r1) then
  fcut=0.0_dp
elseif(r>r2) then
  fcut=1.0_dp
else
  aux=pi*((r1+r2)*0.5_dp-r)/(r2-r1)
  fcut =0.5_dp-0.5_dp*sin(aux)
endif      

end function 

function fcut2_dfcut(r,r1,r2,dfcut) result(fcut)
! The flipped fcut and its derivative (dfcut)
use gems_constants, only: pi,pio2
real(dp),intent(out)   :: dfcut
real(dp),intent(in)    :: r,r1,r2
real(dp)               :: fcut,aux

if(r<r1) then
  fcut=0.0_dp
  dfcut=0.0_dp
elseif(r>r2) then
  fcut=1.0_dp
  dfcut=0.0_dp
else
  aux=pi*((r1+r2)*0.5_dp-r)/(r2-r1)
  fcut =0.5_dp-0.5_dp*sin(aux)
  dfcut=cos(aux)*pio2/(r2-r1)
endif      

end function 

function bond_order(r,r1,r2)
use gems_constants, only: pi
real(dp),intent(in)  :: r,r1,r2
real(dp)             :: aux
real(dp)             :: bond_order

if(r<r1) then
  bond_order=1.0_dp
elseif(r>=r2) then
  bond_order=0.0_dp
else
  aux=0.5_dp*(r1+r2)
  aux=pi*(r-aux)/(r2-r1) 
  bond_order = 0.5_dp-sin(aux)*0.5_dp
endif      
end function 

function smooth(r,r1,r2,dfs)
use gems_constants, only: pi
real(dp)                       :: r,r1,r2
real(dp)                       :: smooth
real(dp),optional,intent(out)  :: dfs
                      
 if(r<r1)then
   smooth=1.0_dp
   if(present(dfs)) dfs=0.0_dp
 else
   smooth=0.5_dp+0.5_dp*dcos(pi*(r-r1)/(r2-r1))
   if(present(dfs)) dfs=-0.5_dp*pi/(r2-r1)*dsin(pi*(r-r1)/(r2-r1))
 end if

end function 

! Vectors

!No tiene mucho sentido definir un producto cruz entre dos vecoters que no sea en 3D. Cualquier
!definicion seria incosistente. 

! Sin embargo esto plante un problema cuando se pide esta operacion en ciertas
! subrutinas pero se quiere compliar el codigo en otra dimension. 

! Necesitamos truquear el significado de este vector. Se puede
! generalizar de distintas formas. 

! Una de las generalizaciones implica perder el ser operador binario.
! Para un n arbitrario se necesitan n-1 vectores. Otra generalizacion puede ser
! pensando en que sentido fisico pueda tener. Abajo pongo algunas

  !function cross_product_1(a,b) result(c)
  !  ! Para compatibilidad
  !  real(dp)                :: c
  !  real(dp),intent(in)     :: a
  !  real(dp),intent(in)     :: b

  !  c = 0.0_dp

  !end function
 
  !function cross_product_2(a,b) result(c)
  !  ! Para compatibilidad, devuelvo un escalar
  !  real(dp)                :: c
  !  real(dp),intent(in)     :: a(2)
  !  real(dp),intent(in)     :: b(2)

  !  c = a(1)*b(2) - a(2)*b(1)

  !end function

  !function cross_product_2(a,b) result(c)
  !  real(dp)                :: c(3)
  !  real(dp),intent(in)     :: a(3)
  !  real(dp),intent(in)     :: b(3)
  !
  !  c(1) = 0.0_dp
  !  c(3) = a(1)*b(2) - a(2)*b(1)
  !end function 
 

  !function cross_product_2(a) result(c)
  !  ! Para compatibilidad con el de obtener vector ortogonal, devuelvo un
  !  vector
  !  real(dp)                :: c
  !  real(dp),intent(in)     :: a(2)

  !  c(1) = a(2)
  !  c(2) = -a(1)

  !end function
    
  function norm_vect(a)
   ! the dsqrt(dot_product(a,a)) is the intrinsic way, but i dont know what is better
   real(dp)                :: norm_vect
   real(dp),intent(in)     :: a(:)
   norm_vect = dsqrt(dot_product(a,a))
  end function

  function cross_product(a,b) result(c)
    real(dp)                :: c(3),a1(3),b1(3)
    real(dp),intent(in)     :: a(:)
    real(dp),intent(in)     :: b(:)
    integer                 :: i

    ! make 3d vectors
    a1(1:size(a))=[ (a(i),i=1,size(a)) ]
    b1(1:size(b))=[ (b(i),i=1,size(b)) ]
    if(size(a)<3) a1(size(a)+1:3)=[ (0.0_dp,i=size(a)+1,3) ]
    if(size(b)<3) b1(size(b)+1:3)=[ (0.0_dp,i=size(b)+1,3) ]

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function

  function diff_vect(a,b,d) 
    logical                 :: diff_vect
    real(dp),intent(in)     :: a(:)
    real(dp),intent(in)     :: b(:)
    real(dp),intent(in)     :: d
    integer                 :: i

 
    diff_vect = .false.
    do i = 1,size(a)
      if (abs(a(i)-b(i))>=d) then
        diff_vect = .true.
        exit
      endif
    enddo       

  end function

  function dot(a,b) result(c)
    ! Muy raro: dot_product tiene un problema de precision si un vector es
    ! unos 4 ordenes de magnitud mas chico que el otro. Asi lo arregle yo.
    real(dp)                :: c
    real(dp),intent(in)     :: a(:)
    real(dp),intent(in)     :: b(:)

    !c=dot_product(a,a)
    !c=dot_product(c*a,b)/c
    c=sum(a*b)

  end function

! triangular matrix
!------------------

   function st_wd_index(i,j,n)
   ! Indice de una matriz simetrica triangular con diagonal
   ! Dado un par i,j devuelve el valor M_ij de la matriz:
   !  1 2 3
   !  2 4 5
   !  3 5 6
   ! Si la entrada es incoherente,devuelve 0
     integer,intent(in)         :: i,j,n
     integer                    :: i2,j2
     integer                    :: st_wd_index

     j2=max(i,j)
     if (n<j2) then
       st_wd_index=0
       return
     endif
     i2=min(i,j)
     if (i2<=0) then
       st_wd_index=0
       return
     endif
     st_wd_index=j2+(i2-1)*(n-i2/2)

   end function

! deteminant
!-----------

  function det_3x3(a)
    real(dp)                :: det_3x3
    real(dp),intent(in)     :: a(3,3)
    det_3x3=          a(1,1) *( a(2,2)*a(3,3) - a(3,2)*a(2,3) )
    det_3x3=det_3x3 - a(1,2) *( a(2,1)*a(3,3) - a(2,3)*a(3,1) )
    det_3x3=det_3x3 + a(1,3) *( a(2,1)*a(3,2) - a(3,1)*a(2,2) )
  end function

! ax=b
!-----

  function cramers_3x3(a,b)
    real(dp)                :: cramers_3x3(3),deta
    real(dp),intent(in)     :: a(3,3), b(3)
    deta=det_3x3(a)
    cramers_3x3(1) = det_3x3( reshape([  b   ,a(:,2),a(:,3)],[3,3]) )/deta
    cramers_3x3(2) = det_3x3( reshape([a(:,1),  b   ,a(:,3)],[3,3]) )/deta
    cramers_3x3(3) = det_3x3( reshape([a(:,1),a(:,2),  b   ],[3,3]) )/deta
  end function 

! Se puede usar gaussj del numerical recipes

! tridiagonal matrix
!-------------------


! lineal transformation
!----------------------

  function ortogonal(v) result(vrot)
   ! Me da un vector ortogonal a v. Estilo Gram-Schmidt
   real(dp),intent(in)  :: v(:)
   real(dp)             :: vaux(size(v)),vrot(size(v))

   do 
     vaux(1)=ranu()
     vaux(2)=ranu()
     vaux(3)=ranu()
     vrot=vaux-dot_product(vaux,v)*v/dot_product(v,v)
     if(dabs(dot_product(vrot,v))<1.0e-3_dp) return
   enddo

  end function ortogonal

  function rodrigues_rotation(ang,eje,v) result(vrot)
   ! Utiliza la formula de Rodrigues para rotar un v alrededor de un eje un
   ! ang determinado. 

   real(dp),intent(in)  :: ang,v(:),eje(3)
   real(dp)             :: vrot(size(v)),cosa,aux(3)


   cosa=cos(ang)
   aux=cross_product(eje,v)
   vrot=v*cosa+aux(1:dm)*sin(ang)+eje*dot_product(eje(1:dm),v)*(1-cosa)

  end function rodrigues_rotation

  function rotmatrix_3x3(angulo)
   ! calcula las matrices de rotacion con el angulo dado. si la rotacion es una
   ! rotacion con componentes en cada eje esta subrutina realiza primero la
   ! rotacion en el orden z, y, x. (tener en cuenta que las rotaciones no son
   ! conmutativas.)

   real(dp),intent(in)  :: angulo(3)
   real(dp)             :: rotx(3,3),roty(3,3),rotz(3,3),rotmatrix_3x3(3,3),rad(3)

    rad = angulo*pi/180.0_dp

    rotx=0.0_dp
    rotx(1,1)= 1.0_dp
    rotx(2,2)= 1.0_dp
    rotx(3,3)= 1.0_dp

    roty=0.0_dp
    roty(1,1)= 1.0_dp
    roty(2,2)= 1.0_dp
    roty(3,3)= 1.0_dp

    rotz=0.0_dp
    rotz(1,1)= 1.0_dp
    rotz(2,2)= 1.0_dp
    rotz(3,3)= 1.0_dp

    rotx(2,2)=dcos(rad(1))
    rotx(2,3)=-dsin(rad(1))
    rotx(3,2)=dsin(rad(1))
    rotx(3,3)=dcos(rad(1))

    roty(1,1)=dcos(rad(2)) 
    roty(1,3)=dsin(rad(2)) 
    roty(3,1)=-dsin(rad(2))
    roty(3,3)=dcos(rad(2)) 

    rotz(1,1)=dcos(rad(3)) 
    rotz(1,2)=-dsin(rad(3)) 
    rotz(2,1)=dsin(rad(3))
    rotz(2,2)=dcos(rad(3)) 

    rotmatrix_3x3=matmul(matmul(rotx,roty),rotz)

  end function rotmatrix_3x3

  function rotmatrix_2x2(angulo)
   ! calcula las matriz de rotacion en el plano xy con el angulo dado.

   real(dp),intent(in)  :: angulo
   real(dp)             :: rotmatrix_2x2(2,2),rad

    rad = angulo*pi/180.0_dp

    rotmatrix_2x2(2,2)=dcos(rad)
    rotmatrix_2x2(1,1)=dcos(rad) 
    rotmatrix_2x2(1,2)=-dsin(rad) 
    rotmatrix_2x2(2,1)=dsin(rad)

  end function rotmatrix_2x2

! Coordenadas polares
!--------------------

  subroutine xyz_polares(r)
   ! devuelve las coordenadas cartesianas dado las polares

   real(dp),intent(inout)  :: r(:)
   real(dp)                :: v(size(r))
   real(dp)                :: aux
    
    !if (dm==3)
    aux=r(3)*sin(r(2))
    v(1)=aux*cos(r(1))
    v(2)=aux*sin(r(1))
    v(3)=r(3)*cos(r(2))
    r=v

  end subroutine xyz_polares



! found eigenvalue ande eigenvectors
!-----------------------------------

recursive subroutine qsortc(a)
  ! recursive fortran 95 quicksort routine
  ! sorts real numbers into ascending numerical order
  ! author: juli rew, scd consulting (juliana@ucar.edu), 9/03
  ! based on algorithm from cormen et al., introduction to algorithms,
  ! 1997 printing

  ! made f conformant by walt brainerd 

  ! use: call qsortc(vector)

  real, intent(in out), dimension(:) :: a
  integer :: iq

  if(size(a) > 1) then
     call partition(a, iq)
     call qsortc(a(:iq-1))
     call qsortc(a(iq:))
  endif
end subroutine qsortc

subroutine partition(a, marker)
  real, intent(in out), dimension(:) :: a
  integer, intent(out) :: marker
  integer :: i, j
  real :: temp
  real :: x      ! pivot point
  x = a(1)
  i= 0
  j= size(a) + 1

  do
     j = j-1
     do
        if (a(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (a(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange a(i) and a(j)
        temp = a(i)
        a(i) = a(j)
        a(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine partition

subroutine partition_integer(a, marker)
  integer, intent(in out), dimension(:) :: a
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
  integer :: x      ! pivot point
  x = a(1)
  i= 0
  j= size(a) + 1

  do
     j = j-1
     do
        if (a(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (a(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange a(i) and a(j)
        temp = a(i)
        a(i) = a(j)
        a(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine partition_integer


! Arrays
!=======

! Binary search
! -------------

subroutine binleft(a,val,m,found)
! Find the index m where a(m)<=val<a(m+1) using a binary search.
! The function return true if a(m)=val and false otherwise.
! a(:) must be sorted.
integer,intent(in)  :: a(:),val
integer             :: r,l,j
integer,intent(out) :: m
logical,intent(out) :: found

! Left and right bounds
l=lbound(a,1)-1
r=ubound(a,1)+1

! Shrink bounds
do while (r>l+1)
  m=(r+l)/2
  j=a(m)
  
  if(val==j) then
    found=.true.
    return
  endif

  if(val>j) then
    l=m
  else
    r=m
  endif
enddo  

! Worst case
if(r<=size(a)) then
  if(a(r)==val) then
    m=r
    found=.true.
    return
  endif
endif

! Not found
found=.false.
m=l 

end subroutine binleft
  
! Linear search
! -------------
  
function linsearch(a,val) result(m)
! Find the index m where a(m)==val using a linear search.
! Return 0 if val is not found in a(:).
! a(:) must be sorted.
integer,intent(in)  :: a(:),val
integer             :: m,j

! Search
do m=1,size(a)
  j=a(m)
  if(val>=j) exit
enddo  

! Check success
if(val==j) return

! Not found
m=0
end function linsearch
 

 
end module gems_algebra
