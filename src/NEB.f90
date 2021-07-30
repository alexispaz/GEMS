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

 
module gems_neb

 ! G. Henkelman and H. Jonsson, Improved tangent estimate in the nudged
  ! elastic band method for finding minimum energy paths and saddle points,
  ! J. Chem. Phys., 113, 9978 (2000). 

  use gems_program_types
  use gems_integration, only: integrate, integrate_cli
  use gems_interaction
  use gems_inq_properties
  use gems_set_properties
  use gems_errors
  use gems_output
  use gems_algebra
  use gems_constants,only:find_io
  
  implicit none
  private

  real(dp),allocatable     :: ks(:)
  real(dp),allocatable     :: eimg(:),eimg_old(:)  ! Energia de cada imagen
  real(dp)                 :: sum_eimg=0.0_dp        ! Energia de todas las imagenes
  real(dp)                 :: etrans=0.0_dp          ! Energia de transicion
  integer                  :: nimg,nebdim          ! Esto sera igual a gsel%nat*3
 
  real(dp)                    :: mass, one_mass, facc
  
  real(dp),allocatable        :: tplus(:), tmin(:), qneb(:)
  real(dp),target,allocatable :: pneb(:,:),fneb(:,:),vneb(:,:)
  integer,allocatable         :: fiximg(:)         ! Guarda en orden los indices de las imagenes fijas (entre las que se va a interpolar)
  integer                     :: nfiximg            ! Numero de imagenes fijas

  
  public  :: ks,fiximg,nfiximg,nimg,nebdim,pneb,neb,neb_fiximg,write_eneb 
    
  type(integrate)           :: nebit

  contains

subroutine neb(g,k,iter,b_out)
logical,intent(in)        :: b_out
real(dp),intent(in)       :: k ! Constante del resorte en eV
integer,intent(in)        :: iter ! Constante del resorte en eV y numero de interacciones
type(group),intent(inout) :: g
integer                   :: ns,i
real(dp)                  :: delr,fprom
real(dp)                  :: aux(g%nat*dm)
logical                   :: switched, ghosted

! Sudden atom movements require fullghost
ghosted=fullghost
fullghost=.true.

! Cambio al modo de almacenamiento vectorial
call group_switch_vectorial(g,switched)

! Algorithmo de integracion
call nebit%init_ext(g)
call integrate_cli(nebit,'v_verlet')

! Allocateo los vectores importantes
! pneb ya deberia estar allocateado, con el comando neb_fin
allocate(vneb(nimg,nebdim),fneb(nimg,nebdim)) !Estado
allocate(eimg(nimg))                          !Energia
allocate(eimg_old(nimg))                      !Energia auxiliar
eimg(:) = 0.0_dp

! Construyo imagenes intermedias por extrapolacion lineal y calculo la energía
call neb_linin
call set_flag('NebI')
if(b_out) call write_neb(g)
call del_flag('NebI')

!Spring force parameters, first and last images are fixed (no force)
allocate(ks(nimg))        !Constante del resorte 'ev/A**2'
ks(1) = 0.0_dp
ks(nimg) = 0.0_dp
ks(:) = k*ev_ui  
mass = 1.0_dp             !Masa 
one_mass = 1.0_dp/mass
call change_atom_mas(mass)

!Otras
vneb = 0.0_dp
allocate(tplus(nebdim))
allocate(tmin(nebdim))
allocate(qneb(nebdim))         !Fuerza tangente?

do i =1 , nimg-1
  g%pp = pneb(i,:)  ! Cargo la posocion
  call pos_changed()
  call interact(.false.)
  call inq_pot_energy(g)
  eimg(i) = g%epot
  fneb(i,:) = g%pf  ! Cargo la fuerza
enddo
    
call nebfce         !Agrega la fuerza del neb

eimg_old=eimg

!estimacion de dt en la primera iteracion       
! fmax = maxval(dabs(fneb))
! dt = dsqrt(2.0_dp*mass*0.005_dp/fmax)

!estimacion de dt en la primera iteracion (segun exneb)      
fprom=0.0_dp
do i =2,nimg-1
  fprom=fprom+dot_product(fneb(i,:),fneb(i,:))
enddo
fprom=dsqrt(fprom/(real(nimg-2)*real(nebdim)))
dt = dsqrt(2.d0*mass*5.d-5/fprom) ! 5d-5 debe ser el delta x

do ns = 1,iter
  !'desplazamiento del paso      =',delr
  !'energia del paso (promedio)  =',sum_eimg/real(nimg)

  ! Velocidades necesita las fuerzas
  call nebvel 

  delr = 0.0_dp
  do i =2 , nimg-1

    facc = facc + dot_product(fneb(i,:),fneb(i,:))
    aux(:) = pneb(i,:) 

    g%pp = pneb(i,:) ! Cargo la posocion
    g%pf = fneb(i,:) ! Cargo la fuerza
    g%pv = vneb(i,:) ! Cargo la velocidad
    call nebit%stepa()
    pneb(i,:) = g%pp  ! Cargo la posocion

    aux(:) = pneb(i,:)-aux(:)
    delr = delr + dot_product(aux,aux)

    call pos_changed()
    call interact(.false.)
    fneb(i,:) = g%pf  ! Cargo la fuerza
    vneb(i,:) = g%pv  ! Cargo la velocidad

    call inq_pot_energy(g)
    eimg(i) = g%epot

  enddo

  !Convergencia del neb e info util
  sum_eimg = sum(eimg) 
  etrans=maxval(eimg)-eimg(1)
  facc = dsqrt(facc/dfloat((nimg-2)*nebdim))
  eimg_old=eimg
  delr = dsqrt(delr/(dfloat(nimg-2)*dfloat(nebdim)))

  if(b_out) then
    call wlog('NEB'); write(logunit,'(i0,4(1x,e25.12))') ns, delr, sum_eimg*ui_ev , facc, etrans*ui_ev
  endif

  !if(delr<0.01) exit

  ! Fuerzas y velocidades
  call nebfce
  call nebvel

  do i =2 , nimg-1
    g%pf = fneb(i,:) ! Cargo la fuerza
    g%pv = vneb(i,:) ! Cargo la velocidad
    call nebit%stepb()
    vneb(i,:) = g%pv  ! Cargo la velocidad
  enddo

enddo

call set_flag('NebF')
if(b_out) call write_neb(g)
call del_flag('NebF')


! la => gsel%alist%next
! do i = 1,gsel%nat
!   la%o%one_mass = 1.0_dp/z_mass(la%o%z)
!   la => la%next
! enddo     
nebdim = 0
nimg = 0

if(switched) call group_switch_objeto(nebit)

fullghost=ghosted

deallocate(vneb,fneb) !Estado
deallocate(eimg)      !Energia
deallocate(eimg_old)  !Energia auxiliar
deallocate(tplus)
deallocate(tmin)
deallocate(qneb)      !Fuerza tangente?
deallocate(ks)        !Constante del resorte 'ev/A**2'

end subroutine neb 
  
  subroutine change_atom_mas(m)
    type(atom_dclist),pointer :: la
    integer                   :: i
    real(dp)                  :: m

    la => gsel%alist%next
    do i = 1,gsel%nat
      !señalo al hypervector desde cada atomo
      la%o%one_mass = 1.0_dp/m
      la => la%next
    enddo
 
  end subroutine change_atom_mas

  subroutine nebfce
    ! Estos vectores son solo sobre el numero de atomos moviles
 
    ! G. Henkelman and H. Jonsson, Improved tangent estimate in the nudged
    ! elastic band method for finding minimum energy paths and saddle points,
    ! J. Chem. Phys., 113, 9978 (2000). 

    implicit none       
    real(dp)         :: tplusmod, tminmod
    real(dp)         :: delvplus, delvmin, fspring
    real(dp)         :: norm, dotvq, delvmax, delvinf
    integer          :: i
 

    facc = 0.0_dp
    do i = 2, nimg-1
 
      ! Tao plus an tao minus (eq 9, J. Chem. Phys., 113, 9978 (2000))
      tplus = pneb(i+1,:) - pneb(i  ,:)
      tmin  = pneb(i  ,:) - pneb(i-1,:)
      delvmin =dabs(eimg(i)-eimg(i-1))

      ! Energies deltas (eq 11, J. Chem. Phys., 113, 9978 (2000))
      delvplus=dabs(eimg(i+1)-eimg(i))
      delvmax = max(delvmin,delvplus)
      delvinf = min(delvmin,delvplus)

      ! Tanget vector (eq 10 and 8, J. Chem. Phys., 113, 9978 (2000))
      if(eimg(i+1)>eimg(i))then
        if(eimg(i)>eimg(i-1))then
          qneb = tplus
        else
          if(eimg(i+1)>eimg(i-1))then
            qneb =  delvmax*tplus+delvinf*tmin
          elseif(eimg(i+1)<eimg(i-1))then
            qneb = delvinf*tplus+delvmax*tmin
          endif
        endif
      else
        if(eimg(i)>eimg(i-1))then
          if(eimg(i+1)>eimg(i-1))then
            qneb =  delvmax*tplus+delvinf*tmin
          elseif(eimg(i+1)<eimg(i-1))then
            qneb = delvinf*tplus+delvmax*tmin
          endif
        else
          qneb = tmin
        endif 
      endif 
      
      ! "Finally, the tanget vector needs to be normalized. With this modified
      ! tangent, the elastic band is well behaved and con-verges rigorously to
      ! the MEP if sufficient number of images are included in the band." (p.
      ! 9981, J.  Chem. Phys., 113, 9978 (2000))
      norm = dot_product(qneb,qneb)
      qneb = qneb/dsqrt(norm)
 
      ! force calculation
      dotvq = dot_product(fneb(i,:),qneb)

      ! Spring force module (eq 12, J. Chem. Phys., 113, 9978 (2000)) 
      ! "This ensure equal spacing of the image (when the same spring constant
      ! is used for the spring) even in regions of high curvature where the
      ! angle between tplusmod and tminmod is large " (J. Chem. Phys., 113, 9978
      ! (2000))
      tplusmod= dsqrt(dot_product(tplus,tplus))
      tminmod = dsqrt(dot_product(tmin,tmin))
      fspring = ks(i)*(tplusmod-tminmod)

      ! Image force (eq 3, J. Chem. Phys., 113, 9978 (2000)) 
      fneb(i,:) = fneb(i,:)-dotvq*qneb(:)+qneb(:)*fspring
      ! Esta otra expresion tenia problemas de convergencia numerica...
      !fneb(i,:) = fneb(i,:)-qneb(:)*(dotvq+fspring)

    enddo

  end subroutine nebfce

  subroutine nebvel
    !composition of the velocities, depending if they point in the direction of
    !the force or not
    real(dp)         :: dotvf
    integer          :: i
 
    do i = 2, nimg-1
      dotvf = dot_product(vneb(i,:),fneb(i,:))
      if(dotvf>0.0_dp)then
        vneb(i,:) = dotvf*fneb(i,:)/dot_product(fneb(i,:),fneb(i,:))
      else
        vneb(i,:)=0.0_dp 
      endif 
    enddo
 
  end subroutine nebvel

  subroutine neb_linin
    !Construyo imagenes intermedias por extrapolacion lineal y les calculo la
    !energía
    integer                   :: i,j
    real(dp),allocatable      :: aux(:)

    call werr('Insufficient fixed images',nfiximg<2)
    allocate(aux(nebdim))

    !Extrapolacion lineal entre las imagenes fijas
    call sort_int(fiximg(1:nfiximg))
    do j = 1,nfiximg-1

      aux=(pneb(fiximg(j+1),:)-pneb(fiximg(j),:))
      call werr('Images to close',sum(abs(aux))<1e-8)
      aux=aux/float(fiximg(j+1)-fiximg(j))

      do i = fiximg(j)+1,fiximg(j+1)-1 
        pneb(i,:) = pneb(i-1,:) + aux
      enddo
    enddo
    
    deallocate(aux)
    
  end subroutine neb_linin

  subroutine neb_fiximg(g,i)
    type(group),intent(in)     :: g
    type (atom_dclist),pointer :: la
    integer,intent(in)         :: i
    integer                    :: j

    call werr('Images must have the same number of atoms',g%nat*dm/=nebdim)

    nfiximg = nfiximg + 1 
    fiximg(nfiximg) = i

    la => g%alist
    do j = 1,g%nat
      la => la%next
      pneb(i,(j-1)*dm+1:j*dm) =  la%o%pos
    enddo 

  end subroutine neb_fiximg

subroutine write_neb(g)  
! Output the energy of the NEB images and its reaction coordinates
! TODO: they will always be X-spaced? So this is not needed
type(group),intent(inout) :: g
real(dp),allocatable         :: gcr(:)
real(dp)                     :: aux,sumrad,bkptime
integer                      :: i
real(dp),allocatable         :: aux2(:) ! Auxiliares

! Calculo de la coordenada de reaccion
allocate(gcr(nimg))
allocate(aux2(nebdim))
gcr(1)=0.0_dp   
sumrad=0.0_dp
do i = 2, nimg  
  ! Calculo la distancia
  aux2(:) = pneb(i,:)-pneb(i-1,:)
  aux = dsqrt(dot_product(aux2,aux2))
  gcr(i) = gcr(i-1) + aux
  sumrad=sumrad+aux
enddo       

!Salida general
bkptime=time
do i =1, nimg
  time=gcr(i)/sumrad/time_out
  g%pp = pneb(i,:) ! Cargo la posocion
  g%pf = fneb(i,:) ! Cargo la fuerza
  call pos_changed()
  call interact(.true.)
  call write_out_force(1)
enddo
time=bkptime

deallocate(aux2)
deallocate(gcr)


end subroutine write_neb

subroutine write_eneb(op) 
! Output the energy of the NEB images
class(outpropa)           :: op
op%f(1:nimg)  = op%f(1:nimg) + eimg(1:nimg)
end subroutine write_eneb
 
end module gems_neb
