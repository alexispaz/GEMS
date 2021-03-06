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
! This file incorporates work derived from the following codes:
!
! L-BFGS code <http://users.iems.northwestern.edu/~nocedal/lbfgs.html>
! by Jorge Nocedal, covered by the following copyright and permission notice:  
!
! Copyright (c)  Jorge Nocedal
! 
!	 L-BFGS is a limited-memory quasi-Newton code for unconstrained optimization. The code has been developed at the Optimization
!	 Center, a joint venture of Argonne National Laboratory and Northwestern University.
!  .
!	 You are welcome to grab the full Unix distribution, containing source code, makefile, and user guide.
!  .
!	 Condition for Use: This software is freely available for educational or commercial purposes. We expect that all publications
!	 describing work using this software quote at least one of the references given below. This software is released under the BSD
!	 License. 
!  .
!	 Authors
!	 Jorge Nocedal
!  .
!	 References
!	 - J. Nocedal. Updating Quasi-Newton Matrices with Limited Storage (1980), Mathematics of Computation 35, pp. 773-782.
!	 - D.C. Liu and J. Nocedal. On the Limited Memory Method for Large Scale Optimization (1989), Mathematical Programming B, 45, 3, pp. 503-528.
!
! Modification of sdrive.f <https://sourceforge.net/p/maxima/code/ci/master/tree/share/lbfgs/sdrive.f>
! by Robert Dodier, covered by the following copyright and permission notice:  
!
! Copyright (c) 2006  Robert Dodier
! 
!  Modification of sdrive.f as retrieved 1997/03/29 from 
!  ftp://ftp.netlib.org/opt/lbfgs_um.shar
!  .
!  This version copyright 2006 by Robert Dodier and released
!  under the terms of the GNU General Public License.


module gems_quasi_newton
use gems_program_types
use gems_set_properties
use gems_inq_properties
use gems_interaction
use gems_output
use gems_errors

save
private 
   
!The only variables that are machine-dependent are XTOL, STPMIN and STPMAX.

!real(dp),parameter :: ftol= 1.0d-4 ! used in decrece wolf condition 

!STPMIN and STPMAX are non-negative which specify lower and uper bounds for the
!step in the line search.  Their default values are 1.e-20 and 1.e+20,
!respectively. These values need not be modified unless the exponents are too
!large for the machine being used, or unless the problem is extremely badly
!scaled (in which case the exponents should be increased).  
real(dp) :: stpmin =1.e-20_dp, stpmax =1.e20_dp

!LP is used as the unit number for the printing of error messages. This
!printing may be suppressed by setting LP to a non-positive value.
integer  :: lp=6

!GTOL controls the accuracy of the line search routine MCSRCH. If the function
!and gradient evaluations are inexpensive with respect to the cost of the
!iteration (which is sometimes the case when solving very large problems) it
!may be advantageous to set GTOL to a small value. A typical small value is
!0.1.  Restriction: GTOL should be greater than 1.D-04.
real(dp) :: gtol= 9.0d-01 
!real(dp) :: gtol= 1.e-3_dp


! XXX: No estoy muy seguro que hace esto
integer,parameter                    :: msave=5

! Since the callback systems seems to me very obscure to understand, I have
! created an internal subroutine to compute f and g. To avoid burry deeper the
! input variables to call this routine, I give the needed variables (like the
! group used to compute f and g) globally. This are those variables.
type(group),pointer            :: gro
logical                        :: b_fgout

public :: lbfgs,lp,lbfgs_minimizator,fix_lbfgs_minimizator,lbfgs_getenergy

 contains

! DRIVERS

function lbfgs_getenergy(g,b_out) result(emin)
  ! Minimiza el grupo g y devuelve solamente la energia minimio. Si g no es
  ! suministrada minimiza la posicion actual
  !integer,parameter       :: msave=7 real(dp),parameter      ::
  type(group),target                   :: g
  real(dp)                             :: emin
  logical,intent(in)                   :: b_out
  real(dp),parameter                   :: eps=1e-8_dp,xtol=1.e-16_dp
  real(dp),dimension(dm*g%nat)         :: gr,diag,scache
  real(dp)                             :: w(dm*g%nat*(2*msave +1)+2*msave),f
  real(dp),dimension(dm*g%nat),target  :: auxp,auxv,auxf,auxa
  integer                              :: iprint,iflag,icall,m,n
  logical,parameter                    :: diagco=.false.
  logical                              :: switched

  ! Needed by get_fg
  gro=>g
  b_fgout=b_out

  ! Cambio al modo de almacenamiento vectorial
  call group_switch_vectorial(g,switched)

  ! Hago checkpoint del estado actual
  auxp=g%pp
  auxv=g%pv
  auxf=g%pf
  auxa=g%pa
  
  ! Set the output unit and print mode
  if (b_out) then
    lp = logunit
    iprint=1
  else
    lp = -1
    iprint=-1
  endif

  n = dm*g%nat
  iflag=0
  m=5

  call get_fg(f,gr)
  call lbfgs(n,m,g%pp,f,gr,iprint,eps,xtol,w,iflag,scache,get_fg)

  ! Salida y retomo el checkpoint
  emin=f*ev_ui
  g%pp=auxp
  g%pv=auxv
  g%pf=auxf
  g%pa=auxa
  
  ! Retomo el modo de almacenamiento anterior
  if(switched) call group_switch_objeto(g)

end function lbfgs_getenergy
 
subroutine lbfgs_minimizator(g,b_out,pmin,emin)
! Minimiza la posicion p.
! Si no es suministrada minimiza la posicion actual
!integer,parameter       :: msave=7
!real(dp),parameter      :: eps=1d-6,xtol=1.0d-16
type(group),target                   :: g
real(dp),intent(out),optional        :: pmin(dm*g%nat),emin
integer,parameter                    :: msave=5
real(dp),parameter                   :: eps=1e-8_dp,xtol=1.e-15_dp
real(dp),dimension(dm*g%nat)         :: gr,diag,scache
real(dp)                             :: w(dm*g%nat*(2*msave +1)+2*msave),f
real(dp),dimension(dm*g%nat),target  :: auxp,auxv,auxf,auxa
integer                              :: iprint,iflag,icall,m,n
logical,parameter                    :: diagco=.false.
logical,intent(in)                   :: b_out
logical                              :: switched, ghosted

! Needed by get_fg
gro=>g
b_fgout=b_out

! Sudden atom movements require fullghost
if(ghost) then
  ghosted=fullghost
  fullghost=.true.
endif

! Cambio al modo de almacenamiento vectorial
call group_switch_vectorial(g,switched)

! Hago checkpoint del estado actual
if (present(pmin)) then
   auxp=g%pp
   auxv=g%pv
   auxf=g%pf
   auxa=g%pa
endif

! Set the output unit
if (b_out) then
  lp = logunit
  iprint=1
else
  lp = -1
  iprint=-1
endif

n = dm*g%nat
iflag=0
m=5

call get_fg(f,gr)
call lbfgs(n,m,g%pp,f,gr,iprint,eps,xtol,w,iflag,scache,get_fg)

! Establezco la modalidad de salida
if (present(emin)) emin=f*ev_ui
if (present(pmin)) then
  pmin=g%pp
  g%pp=auxp
  g%pv=auxv
  g%pf=auxf
  g%pa=auxa
endif

! Retomo el modo de almacenamiento anterior
if(switched) call group_switch_objeto(g)
        
! Retomo el modo ghost anterior
if(ghost) then
  fullghost=ghosted
endif

end subroutine lbfgs_minimizator
 
subroutine fix_lbfgs_minimizator(g,b_out,pmin,emin)
! Minimiza la posicion p.
! Si no es suministrada minimiza la posicion actual
!integer,parameter       :: msave=7
!real(dp),parameter      :: eps=1d-6,xtol=1.0d-16
type(group),target                   :: g
real(dp),intent(out),optional        :: pmin(dm*g%nat),emin
integer,parameter                    :: msave=5
real(dp),parameter                   :: eps=1e-8_dp,xtol=1.e-16_dp
real(dp),dimension(dm*g%nat)         :: gr,diag,scache
real(dp)                             :: w(dm*g%nat*(2*msave +1)+2*msave),f
real(dp),dimension(dm*g%nat),target  :: fixpos
real(dp),dimension(dm*g%nat),target  :: auxp,auxv,auxf,auxa
logical,dimension(dm*g%nat),target   :: fixes
integer                              :: iprint,iflag,icall,m,n
logical,parameter                    :: diagco=.false.
logical,intent(in)                   :: b_out
type (atom_dclist),pointer           :: la
integer                              :: i,j
logical                              :: switched

! ! Needed by get_fg
! gro=>g
! b_fgout=b_out
!
! ! Cambio al modo de almacenamiento vectorial
! call group_switch_vectorial(g,switched)
!
! ! Hago checkpoint del estado actual
! if (present(pmin)) then
!    auxp=g%pp
!    auxv=g%pv
!    auxf=g%pf
!    auxa=g%pa
! endif
!
! n = dm*g%nat
!
! ! Set the output unit
! if (b_out) then
!   lp = logunit
! else
!   lp = -1
! endif
!
! ! Guardo la posicion como la posicion a cuidar
! fixpos=g%pp
!
! ! Genero vector del tamaño del grupo con fijacion
! la => g%alist%next  
! do i = 1,g%nat
!   j=la%o%idv
!   fixes((i-1)*dm+1:i*dm) = fix(j:j+dm) 
!   la => la%next
! enddo 
!
! if (b_out) then
!   iprint=1
! else
!   iprint=-1
! endif
!
! iflag=0
! icall=0
! m=5
!
! !we allow at most nfevalmax evaluations of f and g
! do icall=1,nfevalmax
!   ! Enegia y demas yerba
!   call pos_changed() ! Fundamental!
!   where(fixes)
!      g%pp = fixpos
!   end where
!
!   call interact(.false.)
!   where(fixes)
!      g%pf = 0.0_dp
!   end where
!
!   call inq_pot_energy(g)
!   f = g%epot!*ui_ev
!   gr = -g%pf
!
!   call lbfgs(n,m,g%pp,f,gr,iprint,eps,xtol,w,iflag,scache,get_fg)
!
!   if (b_out) call write_out(1,icall)
!
!   if(iflag<=0) exit
!
! enddo
!
! if (b_out) then
!   if (ical==nfevalmax+1) then
!     call wwan('Limit of function evaluation reached. Try increasing nfevalmax.');
!   else
!     call wlog('LBFGS');
!     write(logunit,'(a,i4,a)') 'Terminated after ',icall-1,' evaluations.'
!   endif
! endif
!    
! ! Establezco la salida
! if (present(emin)) emin=f
! if (present(pmin)) then
!   pmin=g%pp
!   g%pp=auxp
!   g%pv=auxv
!   g%pf=auxf
!   g%pa=auxa
! endif
!
! ! Retomo el modo de almacenamiento anterior
! if(switched) call group_switch_objeto(g)
!
end subroutine fix_lbfgs_minimizator
    
subroutine lbfgs(n,m,x,f,g,iprint,eps,xtol,w,iflag,scache,get_fg,get_diag)
 
! Limited memory BFGS method for large scale optimization by Jorge Nocedal (July
! 1990). 

! Alexis Paz (October 2018) have replaced the callback system of the original
! routine to clarify the execution flow of the code. Now, the code expect the
! existance of two wrap subroutines to compute Hk or f and g as needed by the
! method. The user must specify this routines with the specific interfaces
! declared (any required variables must be supplied globally). This system
! allows to remove a lot of goto in the code and introduce do loops for the
! iterations. 
!
! Here the original routine header by Jorge Nocedal:

! This subroutine solves the unconstrained minimization problem

!                 min F(x),    x= (x1,x2,...,xN),

! using the limited memory BFGS method. The routine is especially effective on
! problems involving a large number of variables. In a typical iteration of this
! method an approximation Hk to the inverse of the Hessian is obtained by
! applying M BFGS updates to a diagonal matrix Hk0, using information from the
! previous M steps.  The user specifies the number M, which determines the
! amount of storage required by the routine. The user may also provide the
! diagonal matrices Hk0 if not satisfied with the default choice.  The algorithm
! is described in "On the limited memory BFGS method for large scale
! optimization", by D. Liu and J. Nocedal, Mathematical Programming B 45 (1989)
! 503-528.
!
! The user is required to calculate the function value F and its gradient G. The
! steplength is determined at each iteration by means of the line search routine
! MCVSRCH, which is a slight modification of the routine CSRCH written by More'
! and Thuente.

integer,intent(in)     :: n,&  ! Number of variables
                          m    ! Number of corrections in the BFGS update. 
                               ! 3<=M<=7 is recommended; large values result in excessive computing time.
real(dp),intent(inout) :: x(n) ! In: initial estimate of the solution. Out: best point found (usually a solution).               

real(dp)   :: w(n*(2*m+1)+2*m) ! Work array
integer, parameter  :: itermax=1000 ! We allow itermax LBFGS iterations
real(dp)   :: xtol             ! Tolerance for the relative width of the interval of uncertainty
real(dp)   :: eps              ! Solution accuracy. Subroutine terminates when |G|<eps*max(1,|X|) 
real(dp)   :: f                ! Function F at X
real(dp)   :: g(n)             ! Gradient G at X
real(dp)   :: scache(n)
real(dp)   :: diag(n) ! Diagonal matrix Hk0 provided by the user at each iteration (only when diagco=true)
logical    :: diagco  ! If true, routine will return back (IFLAG=2) to get diag(:) matrix.

integer    :: iprint  !  Control output. <0: off, 0: first and last iteration, > 0: every iprint iterations

integer    :: iflag   ! It may have the following values:
                      !  0  No errors (must be set to cero on initial entry).
                      ! -1  The line search routine MCSRCH failed. The
                      !     parameter INFO provides more detailed information
                      !     (see also the documentation of MCSRCH)
                      ! -2  The i-th diagonal element of the diagonal inverse
                      !     Hessian approximation, given in DIAG, is not
                      !     positive.
                      ! -3  Improper input parameters for LBFGS (N or M are
                      !     not positive).

real(dp)            :: gnorm,stp1,ftol,stp,ys,yy,sq,yr,beta,xnorm
integer             :: iter,nfun,point,ispt,iypt,maxfev,info, &
                       bound,npt,cp,i,nfev,inmc,iycn,iscn

logical             :: finish
  
! Interface for the subroutine to compute f and g
interface
  subroutine get_fg(f,g)
  import dp
  real(dp),intent(out)   :: g(:),f
  end subroutine       
end interface
   
! Interface for the subroutine to compute the diagonal Hk0 matrix
interface
  subroutine get_diag_interface(diag)
  import dp
  real(dp),intent(out)   :: diag(:)
  end subroutine 
end interface
procedure(get_diag_interface), optional :: get_diag

save
 
! Set diagco if get_diag is supply
if(present(get_diag)) diagco=.true.

! Check for input errors
if(n<=0.or.m<=0) then
  if(lp>0) then
    call wwan('Wrong input (n or m are not positive)')
  endif
  iflag= -3
  return
endif

! Restrict gtol above 1e-4
if(gtol<=1.e-4_dp) then
  if(lp>0) then
    call wwan('Since gtol<=1e-4 it will be reset to 0.9.')
  endif
  gtol=9.d-01
endif

! Check for diag as input
if(diagco) then
  do i=1,n
    if (diag(i)<=0._dp) then
      if(lp>0) then
        call wwan(); write(lp,'(i5,a)') i,'-th diag element of the inverse hessian approx is < 0'
      endif
      iflag=-2
      return
    endif
  enddo
else
   diag(1:n)= 1._dp
endif

! Initialize variables
nfun= 1
point= 0
finish= .false.
ispt= n+2*m
iypt= ispt+n*m
              
! Parameters for line search routine
ftol=1.e-4_dp
maxfev= 100

!The work vector w is divided as follows:
! w(1:n)             gradient and other temporary information.
! w(n+1:n+m)         scalars rho.
! w(n+m+1:ispt)      numbers alpha used in the formula that computes h*g.
! w(ispt+1:iypt)     last m search steps.
! w(iypt+1:iypt+nm)  last m gradient differences.
!
! the search steps and gradient differences are stored in a
! circular order controlled by the parameter point.

w(ispt+1:ispt+n)= -g(1:n)*diag(1:n)

gnorm= dsqrt(dot_product(g,g))
stp1= 1._dp/gnorm


call lb1(iprint,iter,nfun,gnorm,n,m,x,f,g,1._dp,finish)

!Main iteration loop.
do iter=1,itermax

  ! First iteration goes to 165
  info=0
  bound=iter-1
  if(iter==1) go to 165

  if (iter>m) bound=m

  ! "Circular order controlled by the parameter point" (npt=point*n)
  ! Dot product between gradient difference and the search step
  ys= dot_product(w(iypt+npt+1:iypt+npt+n),w(ispt+npt+1:ispt+npt+n))

  if(.not.diagco) then
     ! The norm of the gradient difference
     yy= dot_product(w(iypt+npt+1:iypt+npt+n),w(iypt+npt+1:iypt+npt+n))
     diag(:)= ys/yy
  else
     call get_diag(diag)
  endif

  ! Check for error
  if(diagco) then
    do i=1,n
      if (diag(i)<=0.0_dp) then
        iflag=-2
        if(lp>0) then
          call wwan(); write(lp,'(a,i5,a)') 'iflag= -2.',i,'-th diag element of the inverse hessian approx is < 0'
        endif
        return
      endif
    enddo
  endif

  !COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
  !"Updating quasi-Newton matrices with limited storage",
  !mathematics of computation, vol.24, no.151, pp. 773-782.
  cp= point
  if (point==0) cp=m
  w(n+cp)= 1.0_dp/ys

  w(1:n)= -g(1:n)

  cp= point
  do i= 1,bound
      cp=cp-1
      if (cp== -1)cp=m-1
      ! Dot product between the search step and the gradient 
      sq= dot_product(w(ispt+cp*n+1:ispt+cp*n+n),w(1:n))
      inmc=n+m+cp+1
      iycn=iypt+cp*n
      w(inmc)= w(n+cp+1)*sq
      call daxpy(n,-w(inmc),w(iycn+1:iycn+n),w(1:n))
  enddo

  do i=1,n
    w(i)=diag(i)*w(i)
  enddo

  do i=1,bound
     ! Dot product between the gradient difference and the gradient 
     yr= dot_product(w(iypt+cp*n+1:iypt+cp*n+n),w(1:n)) 
     beta= w(n+cp+1)*yr
     inmc=n+m+cp+1
     beta= w(inmc)-beta
     iscn=ispt+cp*n
     call daxpy(n,beta,w(iscn+1:iscn+n),w(1:n))
     cp=cp+1
     if (cp==m)cp=0
  enddo

  !store the new search direction
  do i=1,n
    w(ispt+point*n+i)= w(i)
  enddo

  !Obtain the one-dimensional minimizer of the function by using the line
  !search routine mcsrch
165    stp=1.0_dp
  if (iter==1) stp=stp1

  w(1:n)=g(1:n)

  call mcsrch(n,x,f,g,w(ispt+point*n+1),stp,ftol,xtol,maxfev,info,nfev,diag)

  ! End with error
  if (info /= 1) then
    iflag=-1
    if(lp>0) then
      call wwan(); write(lp,'(a,i2,a)') 'iflag= -1. line search failed: info= ',info,'see doc of routine mcsrch +++'
      call wwan('possible causes: function or gradient are incorrect or incorrect tolerances')
    endif
    return
  endif

  nfun= nfun + nfev

  !compute the new step and gradient change
  npt=point*n
  do i=1,n
    w(ispt+npt+i)= stp*w(ispt+npt+i)
    w(iypt+npt+i)= g(i)-w(i)
  enddo 
  point=point+1
  if (point==m) point=0

  !termination test
  gnorm= sqrt(dot_product(g,g))
  xnorm= sqrt(dot_product(x,x))
  xnorm= int(max(1._dp,xnorm))
  if (gnorm/xnorm <= eps) finish=.true.

  call lb1(iprint,iter,nfun,gnorm,n,m,x,f,g,stp,finish)

  !We've completed a line search. Cache the current solution vector.
  !      
  !At the end of searching, the solution cache is different from
  !the current solution if the search is terminated in the
  !middle of a line search. That is the case if, for example,
  !the number of function evaluations is the criterion to stop
  !the search instead of the number of line searches or convergence
  !to a prescribed tolerance.

  scache(1:n) = x(1:n)

  if (finish) then
     iflag=0
     call wlog('LBFGS');
     write(logunit,'(a,i4,a)') 'Terminated after ',nfun,' evaluations.'
     return
  endif

enddo

! Maximum iteration reached
if(lp>0) then
  call wwan('Limit of iteration evaluation reached. Try increasing itermax.');
endif
                   

end subroutine

subroutine lb1(iprint,iter,nfun,gnorm,n,m,x,f,g,stp,finish)
! this routine prints monitoring information. the frequency and
! amount of output are controlled by iprint.
integer,intent(in)  :: n,m,iter,nfun,iprint
real(dp),intent(in) :: x(n),g(n),f,gnorm,stp
logical,intent(in)  :: finish
integer             :: i

if(iprint<0) return 

if (iter==0)then
  call wlog('LBFGS'); write(logunit,'(3x,a,5x,a,20x,a,19x,a)') 'I  NFN','FUNC','GNORM','STEPLENGTH'
endif

if ((iprint==0).and.(iter/=1.and..not.finish)) return

if (iprint/=0)then
  if(mod(iter-1,iprint)==0.or.finish)then
    call wlog('LBFGS'); write(logunit,'(2(i4,1x),3x,3(e22.15,2X))') iter,nfun,f,gnorm,stp
  else
    return
  endif

else

  call wlog('LBFGS'); write(logunit,'(2(i4,1x),3x,3(e22.15,2X))') iter,nfun,f,gnorm,stp

endif

if (finish) call wlog('LBFGS', 'Terminated without errors. IFLAG = 0')

end  subroutine

    subroutine daxpy(n,da,dx,dy)
    !constant times a vector plus a vector.
    !uses unrolled loops for increments equal to one.
    !jack dongarra, linpack, 3/11/78.
    integer,intent(in)  :: n
    real(dp)            :: dx(:),dy(:),da
    integer             :: i,ix,iy,m,mp1

    if(n<=0)return
    if (da == 0.0d0) return

    m = mod(n,4)
    if( m == 0 ) go to 40
    do i = 1,m
      dy(i) = dy(i) + da*dx(i)
    enddo
    if( n < 4 ) return
40  mp1 = m + 1
    do i = mp1,n,4
      dy(i) = dy(i) + da*dx(i)
      dy(i + 1) = dy(i + 1) + da*dx(i + 1)
      dy(i + 2) = dy(i + 2) + da*dx(i + 2)
      dy(i + 3) = dy(i + 3) + da*dx(i + 3)
    enddo
    return
    end  subroutine

subroutine mcsrch(n,x,f,g,s,stp,ftol,xtol,maxfev,info,nfev,wa)
!A slight modification of the subroutine CSRCH of More' and Thuente. The
!changes are to allow reverse communication, and do not affect the performance
!of the routine.

!The purpose of mcsrch is to find a step which satisfies a sufficient decrease
!condition and a curvature condition.

!At each stage the subroutine updates an interval of uncertainty with endpoints
!stx and sty. The interval of uncertainty is initially chosen so that it
!contains a minimizer of the modified function
!
!     f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).
!
!If a step is obtained for which the modified function has a nonpositive
!function value and nonnegative derivative, then the interval of uncertainty is
!chosen so that it contains a minimizer of f(x+stp*s).

!The algorithm is designed to find a step which satisfies the sufficient
!decrease condition
!
!      f(x+stp*s) <= f(x) + ftol*stp*(gradf(x)'s),
!  
!and the curvature condition
!
!      abs(gradf(x+stp*s)'s)) <= gtol*abs(gradf(x)'s).
!
!If ftol is less than gtol and if, for example, the function is bounded below,
!then there is always a step which satisfies both conditions. If no step can be
!found which satisfies both conditions, then the algorithm usually stops when
!rounding errors prevent further progress. In this case stp only satisfies the
!sufficient decrease condition.
          
! Global variables expected:
! - gtol: the sufficient directional derivative condition.
! - stpmin and stpmax: lower and upper bounds for the step. 

!ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
!JORGE J. MORE', DAVID J. THUENTE
 
integer,intent(in)    :: n,      & ! Number of variables
                         maxfev    ! Max number of calls to fcn allowed by the end of an iteration.

integer,intent(out)   :: nfev  , & ! The number of calls to fcn.   
                         info  ! Information on the termination reason. The labels are:
                               ! 0: Improper input 
                               ! 1: gtol and ftol are fulfilled
                               ! 2: xtol fulfilled
                               ! 3: fcn==maxfev
                               ! 4: step is at the lower bound stpmin.
                               ! 5: step is at the upper bound stpmax.
                               ! 6: Rounding errors prevent further progress. There may not be a step which
                               ! satisfies the sufficient decrease and curvature conditions. Tolerances may
                               ! be too small.

real(dp),intent(inout):: x(n), & ! In: the base point for the line search. Out: x+stp*s
                         g(n), & ! The gradient of f. In: g(x). Out: g(x+stp*s)
                         f,    & ! In: f(x). Out: f(x+stp*s)
                         stp,  & ! In: initial estimate of a satisfactory step. Out: final estimate.
                         wa(n)   ! Work array.

real(dp),intent(in)   :: s(n), & ! The search direction
                         ftol, & ! Tolerance for the sufficient decrease condition.
                         xtol    ! Tolerance for the relative width of the interval of uncertainty.

save

integer            :: infoc,j
logical            :: brackt,stage1
real(dp)           :: dg,dgm,dginit,dgtest,dgx,dgxm,dgy,dgym, &
                      finit,ftest1,fm,fx,fxm,fy,fym,stx,sty,  &
                      stmin,stmax,width,width1
real(dp),parameter :: p66=0.66_dp,xtrapf=4.0_dp,p5=0.5_dp

infoc = 1

!Check the input parameters for errors.
if (n <= 0 .or. stp <= 0.0_dp .or. ftol < 0.0_dp .or. &
     gtol < 0.0_dp .or. xtol < 0.0_dp .or. stpmin < 0.0_dp &
     .or. stpmax < stpmin .or. maxfev <=0) return

!Compute the initial gradient in the search direction.
dginit=dot_product(g,s)

!Check that s(:) is a descent direction.                                   
if (dginit>=0._dp) then
   call wstd(); write(lp,fmt='(/"  the search direction is not a descent direction")')
   return
endif

!Initialize local variables.
brackt = .false.
stage1 = .true.
nfev = 0
finit = f
dgtest = ftol*dginit
width = stpmax - stpmin
width1 = width/p5
wa(:) = x(:)

!The variables stx, fx, dgx contain the values of the step, function, and
!directional derivative at the best step. The variables sty, fy, dgy contain
!the value of the step, function, and derivative at the other endpoint of the
!interval of uncertainty. The variables stp, f, dg contain the values of the
!step, function, and derivative at the current step.
stx = 0._dp
fx  = finit
dgx = dginit
sty = 0._dp
fy  = finit
dgy = dginit


do nfev=1,maxfev

  !Set the minimum and maximum steps to the present interval of uncertainty.
  if (brackt) then
    stmin = min(stx,sty)
    stmax = max(stx,sty)
  else
    stmin = stx
    stmax = stp + xtrapf*(stp - stx)
  end if

  !Force the step to be within the bounds stpmax and stpmin.
  stp = max(stp,stpmin)
  stp = min(stp,stpmax)

  !If an unusual termination is to occur then let stp be the lowest point
  !obtained so far.
  if ((brackt .and. (stp <= stmin .or. stp >= stmax)) &
     .or. nfev == maxfev .or. infoc == 0            &
     .or. (brackt .and. stmax-stmin <= xtol*stmax)) stp = stx

  !Evaluate the function and gradient at stp and compute the directional
  !derivative.
  x(:) = wa(:) + stp*s(:)
  call get_fg(f,g)
  info=0

  dg=dot_product(g,s)
  ftest1 = finit + stp*dgtest

  !test for convergence.
  if ((brackt .and. (stp<=stmin .or. stp>=stmax)) .or. infoc==0) info = 6
  if (stp==stpmax .and. f<=ftest1 .and. dg<=dgtest)              info = 5
  if (stp==stpmin .and. (f>ftest1 .or. dg>=dgtest))              info = 4
  if (nfev==maxfev)                                              exit
  if (brackt .and. stmax-stmin<=xtol*stmax)                      info = 2
  if (f<=ftest1 .and. abs(dg)<=gtol*(-dginit))                  info = 1

  ! Check for termination.
  if (info/=0) return

  !In the first stage we seek a step for which the modified function has a
  !nonpositive value and nonnegative derivative.
  if (stage1 .and. f<=ftest1 .and. dg>=min(ftol,gtol)*dginit) stage1 = .false.

  !A modified function is used to predict the step only if we have not obtained
  !a step for which the modified function has a nonpositive function value and
  !nonnegative derivative, and if a lower function value has been obtained but
  !the decrease is not sufficient.
  if (stage1 .and. f<=fx .and. f>ftest1) then

    !Define the modified function and derivative values.
    fm = f - stp*dgtest
    fxm = fx - stx*dgtest
    fym = fy - sty*dgtest
    dgm = dg - dgtest
    dgxm = dgx - dgtest
    dgym = dgy - dgtest

    !Update the interval of uncertainty and compute the new step.
    call mcstep(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm,brackt,stmin,stmax,infoc)

    !Reset the function and gradient values for f.
    fx = fxm + stx*dgtest
    fy = fym + sty*dgtest
    dgx = dgxm + dgtest
    dgy = dgym + dgtest
  else

    !Update the interval of uncertainty and compute the new step.
    call mcstep(stx,fx,dgx,sty,fy,dgy,stp,f,dg,brackt,stmin,stmax,infoc)
  end if

  !Force a sufficient decrease in the size of the interval of uncertainty.
  if (brackt) then
    if (dabs(sty-stx) >= p66*width1) stp = stx + p5*(sty - stx)
    width1 = width
    width = dabs(sty-stx)
  end if

enddo

! Maxfev reached
info = 3

end  subroutine mcsrch

subroutine mcstep(stx,fx,dx,sty,fy,dy,stp,fp,derp,brackt,stpmin,stpmax,info)
integer   :: info
real(dp)  :: stx,fx,dx,sty,fy,dy,stp,fp,derp,stpmin,stpmax
logical   :: brackt,bound

!The purpose of mcstep is to compute a safeguarded step for
!a linesearch and to update an interval of uncertainty for
!a minimizer of the function.

!The parameter stx contains the step with the least function
!value. the parameter stp contains the current step. it is
!assumed that the derivative at stx is negative in the
!direction of the step. if brackt is set true then a
!minimizer has been bracketed in an interval of uncertainty
!with endpoints stx and sty.

!  stx, fx, and dx are variables which specify the step,
!    the function, and the derivative at the best step obtained
!    so far. the derivative must be negative in the direction
!    of the step, that is, dx and stp-stx must have opposite
!    signs. on output these parameters are updated appropriately.

!  sty, fy, and dy are variables which specify the step,
!    the function, and the derivative at the other endpoint of
!    the interval of uncertainty. on output these parameters are
!    updated appropriately.

!  stp, fp, and derp are variables which specify the step,
!    the function, and the derivative at the current step.
!    if brackt is set true then on input stp must be
!    between stx and sty. on output stp is set to the new step.

!  brackt is a logical variable which specifies if a minimizer
!    has been bracketed. if the minimizer has not been bracketed
!    then on input brackt must be set false. if the minimizer
!    is bracketed then on output brackt is set true.

!  stpmin and stpmax are input variables which specify lower
!    and upper bounds for the step.

!  info is an integer output variable set as follows:
!    if info = 1,2,3,4,5, then the step has been computed
!    according to one of the five cases below. otherwise
!    info = 0, and this indicates improper input parameters.

!argonne national laboratory. minpack project. june 1983
!jorge j. more', david j. thuente

real(dp)      :: gamm,p,q,r,s,sgnd,stpc,stpf,stpq,theta

info = 0

!Check the input parameters for errors.
if ((brackt .and. (stp <= min(stx,sty) .or. &
    stp >= max(stx,sty))) .or. &
    dx*(stp-stx) >= 0._dp .or. stpmax < stpmin) return

!Determine if the derivatives have opposite sign.
sgnd = derp*(dx/abs(dx))

!First case. A higher function value.
!the minimum is bracketed. if the cubic step is closer
!to stx than the quadratic step, the cubic step is taken,
!else the average of the cubic and quadratic steps is taken.
if (fp > fx) then
  info = 1
  bound = .true.
  theta = 3*(fx - fp)/(stp - stx) + dx + derp
  s = max(abs(theta),abs(dx),abs(derp))
  gamm = s*sqrt((theta/s)**2 - (dx/s)*(derp/s))
  if (stp < stx) gamm = -gamm
  p = (gamm - dx) + theta
  q = ((gamm - dx) + gamm) + derp
  r = p/q
  stpc = stx + r*(stp - stx)
  stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx)
  if (abs(stpc-stx) < abs(stpq-stx)) then
    stpf = stpc
  else
    stpf = stpc + (stpq - stpc)/2
  end if
  brackt = .true.

!second case. a lower function value and derivatives of
!opposite sign. the minimum is bracketed. if the cubic
!step is closer to stx than the quadratic (secant) step,
!the cubic step is taken, else the quadratic step is taken.
else if (sgnd < 0._dp) then
    info = 2
    bound = .false.
    theta = 3*(fx - fp)/(stp - stx) + dx + derp
    s = max(abs(theta),abs(dx),abs(derp))
    gamm = s*sqrt((theta/s)**2 - (dx/s)*(derp/s))
    if (stp > stx) gamm = -gamm
    p = (gamm - derp) + theta
    q = ((gamm - derp) + gamm) + dx
    r = p/q
    stpc = stp + r*(stx - stp)
    stpq = stp + (derp/(derp-dx))*(stx - stp)
    if (abs(stpc-stp) > abs(stpq-stp)) then
       stpf = stpc
    else
       stpf = stpq
       end if
    brackt = .true.

!third case. a lower function value, derivatives of the
!same sign, and the magnitude of the derivative decreases.
!the cubic step is only used if the cubic tends to infinity
!in the direction of the step or if the minimum of the cubic
!is beyond stp. otherwise the cubic step is defined to be
!either stpmin or stpmax. the quadratic (secant) step is also
!computed and if the minimum is bracketed then the the step
!closest to stx is taken, else the step farthest away is taken.
else if (abs(derp) < abs(dx)) then
  info = 3
  bound = .true.
  theta = 3*(fx - fp)/(stp - stx) + dx + derp
  s = max(abs(theta),abs(dx),abs(derp))

  ! the case gamm = 0 only arises if the cubic does not tend
  ! to infinity in the direction of the step.

  gamm = s*sqrt(max(0._dp,(theta/s)**2 - (dx/s)*(derp/s)))
  if (stp > stx) gamm = -gamm
  p = (gamm - derp) + theta
  q = (gamm + (dx - derp)) + gamm
  r = p/q
  if (r < 0._dp .and. gamm /= 0._dp) then
    stpc = stp + r*(stx - stp)
  else if (stp > stx) then
    stpc = stpmax
  else
    stpc = stpmin
    end if
  stpq = stp + (derp/(derp-dx))*(stx - stp)
  if (brackt) then
    if (abs(stp-stpc) < abs(stp-stpq)) then
       stpf = stpc
    else
       stpf = stpq
    end if
  else
    if (abs(stp-stpc) > abs(stp-stpq)) then
       stpf = stpc
    else
       stpf = stpq
    end if
  end if

!fourth case. a lower function value, derivatives of the
!same sign, and the magnitude of the derivative does
!not decrease. if the minimum is not bracketed, the step
!is either stpmin or stpmax, else the cubic step is taken.
else
  info = 4
  bound = .false.
  if (brackt) then
    theta = 3*(fp - fy)/(sty - stp) + dy + derp
    s = max(abs(theta),abs(dy),abs(derp))
    gamm = s*sqrt((theta/s)**2 - (dy/s)*(derp/s))
    if (stp > sty) gamm = -gamm
    p = (gamm - derp) + theta
    q = ((gamm - derp) + gamm) + dy
    r = p/q
    stpc = stp + r*(sty - stp)
    stpf = stpc
  else if (stp > stx) then
    stpf = stpmax
  else
    stpf = stpmin
  end if
end if

!Update the interval of uncertainty. this update does not
!depend on the new step or the case analysis above.
if (fp > fx) then
  sty = stp
  fy = fp
  dy = derp
else
  if (sgnd < 0._dp) then
     sty = stx
     fy = fx
     dy = dx
  end if
  stx = stp
  fx = fp
  dx = derp
end if

!compute the new step and safeguard it.
stpf = min(stpmax,stpf)
stpf = max(stpmin,stpf)
stp = stpf
if (brackt .and. bound) then
  if (sty > stx) then
     stp = min(stx+0.66_dp*(sty-stx),stp)
  else
     stp = max(stx+0.66_dp*(sty-stx),stp)
  end if
end if

end  subroutine
 

subroutine get_fg(f,g)
! Compute gradient and write output  
real(dp),intent(out)   :: g(:),f
integer,save           :: icall=0
icall=icall+1
call pos_changed() ! Fundamental!
call interact(.false.,f)
f = f*ui_ev
g(:) = -gro%pf(:)*ui_ev
if (b_fgout) call write_out(1,icall)
end subroutine 

                  

end module gems_quasi_newton
