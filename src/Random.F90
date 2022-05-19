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

 
module gems_random
use gems_constants, only:dp
use gems_errors
implicit none


private 


! F2003 standar
integer             :: seedsz=1
integer,allocatable :: seedv(:)

! LCG
integer,save          :: idum=-12345789
integer, parameter    :: k4b=selected_int_kind(9)
integer(k4b), save    :: ix=-1,iy=-1

public :: std_init, lcg_init

public :: ranu, rang
public :: init_ran, free_ran, write_chpseed, read_chpseed

procedure(std_init_ran),pointer       :: init_ran
procedure(std_free_ran),pointer       :: free_ran
procedure(std_write_chpseed),pointer  :: write_chpseed
procedure(std_ranu),pointer           :: ranu
    
contains

! Using F2003 standar especifications. The PRNG algorithm used depend on the
! chosen compiler.

subroutine std_init()
  init_ran      => std_init_ran     
  free_ran      => std_free_ran     
  write_chpseed => std_write_chpseed
  ranu          => std_ranu
  call wwan('F2003 standar do not specify PRNG implementation, so this is a compiler dependent choice.')
end subroutine
                           
function std_ranu()
  real(dp)                    :: std_ranu
  call random_number(std_ranu)
end function

subroutine std_init_ran(lseed)
  integer,optional,intent(in)  :: lseed
  integer                      :: i
 
  if(allocated(seedv)) then
    call free_ran()
  endif
  
  ! Find the size of the seed
  call random_seed(size=seedsz)

  ! Allocate the seed
  allocate(seedv(seedsz))

  ! Establishes the user seed... not sure about this
  if(present(lseed)) then
    seedv(:) = [(i,i=1,seedsz)]
    seedv(1) = lseed
    call random_seed(put=seedv)
  else
    call random_seed()
  endif

  ! Print random seed information. TODO: cut the integer in the used bites
  call random_seed(get=seedv)
  call wlog(''); write(logunit,'(a)') "Random Seed: "
  do i=1,seedsz
    call wlog(''); write(logunit,'(2x,i0)') seedv(i)
  enddo
  ! call wlog(''); write(logunit,*) "First random number: ",ranu(stream)
end subroutine
    
subroutine std_free_ran()
  if(.not.allocated(seedv)) return
  deallocate(seedv)
end subroutine

subroutine std_write_chpseed(chpunit)
  integer,intent(in)      :: chpunit
  integer,allocatable     :: seednow(:)

  allocate(seednow(seedsz))
  call random_seed(get=seednow)

  write(chpunit) 'std'
  write(chpunit) seedsz
  write(chpunit) seednow(:)

  deallocate(seednow)

end subroutine std_write_chpseed
                
subroutine std_read_chpseed(chpunit)
  integer,intent(in)      :: chpunit
  integer                 :: sz

  read(chpunit) sz
  if (seedsz/=sz) then
    call free_ran()
    call init_ran()
    call werr('Seed sizes does not match',seedsz/=sz)
  endif

  read(chpunit) seedv(:)
  call random_seed(put=seedv)

end subroutine std_read_chpseed
            

! Linear congruential generator (LCG).
! Not recomended. Use only to enforce reproducibility of test simulations.

subroutine lcg_init()
  init_ran      => lcg_init_ran     
  free_ran      => lcg_free_ran     
  write_chpseed => lcg_write_chpseed
  ranu          => lcg_ranu         
  call wwan('LCG not recommended. Use only to enforce reproducibility in test simulations.')
end subroutine
                        
function lcg_ranu()
  real(dp)               :: lcg_ranu
  integer, parameter     :: ia=16807,im=2147483647,iq=127773,ir=2836
  real(dp), parameter    :: am=nearest(1.0,-1.0)/real(im,dp)
  integer(k4b)           :: k
  if (idum <= 0 .or. iy < 0) then
    iy=ior(ieor(888889999,abs(idum)),1)
    ix=ieor(777755555,abs(idum))
    idum=abs(idum)+1
  end if
  ix=ieor(ix,ishft(ix,13))
  ix=ieor(ix,ishft(ix,-17))
  ix=ieor(ix,ishft(ix,5))
  k=iy/iq
  iy=ia*(iy-k*iq)-ir*k
  if (iy < 0) iy=iy+im
  lcg_ranu=am*ior(iand(im,ieor(ix,iy)),1)  
end function

subroutine lcg_init_ran(lseed)
  integer,optional,intent(in)  :: lseed
  idum=-123456789
  if(present(lseed)) idum=-abs(lseed)

  ! Print random seed information. TODO: cut the integer in the used bites
  call wlog(''); write(logunit,'(a,i0)') "Random Seed: ",idum
  ! call wlog(''); write(logunit,*) "First random number: ",ranu() 
end subroutine
    
subroutine lcg_free_ran()
  idum=123456789
end subroutine

subroutine lcg_write_chpseed(chpunit)
  integer,intent(in)      :: chpunit

  write(chpunit) 'lcg'
  write(chpunit) idum, iy, ix

end subroutine lcg_write_chpseed

subroutine lcg_read_chpseed(chpunit)
  integer,intent(in)      :: chpunit
  read(chpunit) idum, iy, ix

end subroutine lcg_read_chpseed


! Common rutines

subroutine read_chpseed(chpunit)
integer,intent(in)      :: chpunit
character(3)            :: prng

read(chpunit) prng
select case(prng)
case('lcg')
  call lcg_init()
  call lcg_read_chpseed(chpunit)
case('std')
  call std_init()
  call std_read_chpseed(chpunit)
case default
  call werr('checkpoint file corrupted',.true.)
end select

end subroutine read_chpseed
                                       
subroutine  rang(r1,r2)
! Como la gasdev obtiene numeros con distribucion gaussiana de a pares,
! esta subrrutina se asegura que siempre pida numeros de a pares. De no ser
! asi, uno de los numeros se tira a la basura. Esto es para poder recuperar
! exactamente la secuencia de numeros random, dado que la variable de fase
! (gaus_stored) del numerical es interna
  real(dp),intent(out)            :: r1
  real(dp),intent(out),optional   :: r2
  real(dp)                        :: aux

  r1=gasdev()
  aux=gasdev()
  if(present(r2)) r2=aux
end subroutine  rang
  
function gasdev()
  real(dp)                  :: rsq,v1,v2
  real(dp), save            :: g
  real(dp)                  :: gasdev
  logical, save             :: gaus_stored=.false.

  if (gaus_stored) then
    gasdev=g
    gaus_stored=.false.
  else
    do
      v1=2.0_dp*ranu()-1.0_dp
      v2=2.0_dp*ranu()-1.0_dp
      rsq=v1**2+v2**2
      if (rsq > 0._dp .and. rsq < 1._dp) exit
    end do
    rsq=sqrt(-2.0_dp*log(rsq)/rsq)
    gasdev=v1*rsq
    g=v2*rsq
    gaus_stored=.true.
  end if

end function gasdev

end module gems_random
