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


!SPRNG permits the use of several random number streams on each processor. 
#include "sprng_f.h"

 private 

 ! A single process may have several different random number streams available.
 ! These different streams are distinguished by unique ID's. In C++, these IDs are
 ! the SPRNG objects that we use to call SPRNG functions. In FORTRAN, these are
 ! implemented as pointers to memory locations where states of the respective
 ! streams are stored. Since standard FORTRAN 77 does not have a pointer type, we
 ! can store a pointer as an integer of the same size as a C++ pointer. We have
 ! defined a macro called SPRNG_POINTER in the file sprng_f.h that automatically
 ! defines an integer of the correct size on platforms on which SPRNG is
 ! supported. A FORTRAN programmer can then use the type SPRNG_POINTER just as if
 ! it were a FORTRAN data type. However, if the flag -DPOINTERSIZE was used in
 ! building SPRNG FORTRAN executables, then this flag should also be used in
 ! compiling the user's program for this feature to work correctly.
 SPRNG_POINTER,target             :: stream
 SPRNG_POINTER,target,allocatable :: streams(:)

 ! Since SPRNG_POINTER is just an integer pointer and not a fortran pointer I
 ! should have a way to ckeck if it is asociated.
 SPRNG_POINTER,pointer :: sprng_stream=>null()

 !The seed to the random number generator. It is not the starting state of the
 !sequence; rather, it is an encoding of the starting state. It is acceptable
 !(and recommended) to use the same seed for all the streams. Distinct streams
 !are returned based on their seed and the stream number. Only the 31 least
 !significant bits of seed are used in determining the initial starting state of
 !the stream. Higher order bits that are set will be ignored. No warning message
 !is printed if the higher order bits are set.
 integer               :: seed = 985456376

 ! nstreams is the number of distinct streams that will be initialized across
 ! all the processes and must be greater than 0. Otherwise it is reset to 1 and
 ! a warning message is sent to stderr stating that the number of streams has
 ! been reset.
 integer               :: nstreams=1

 ! streamnum is the stream number, typically the processor number, and must be
 ! in [0,nstreams-1]. If it is not in the acceptable range, then an error
 ! message is sent to stderr and the function returns a 0 in C++ or a NULL
 ! pointer in FORTRAN. Note that the number of independent streams for each
 ! type of generator is limited (but at least 10^5). If streamnum is larger
 ! than this number, then a warning message is sent to stderr stating that the
 ! independence of streams cannot be guaranteed. 
 integer               :: streamnum=0
           
 ! gtype indicates the available generators:
 ! 0 - Lagged Fibonacci Generator
 ! 1 - Linear Congruential Generator
 ! 2 - lcg64 
 ! 3 - cmrg  
 ! 4 - mlfg  
 ! 5 - pmlcg 
 integer               :: sprng_type=0

 ! The argument param selects the appropriate parameters (for example, the
 ! multiplier for a Linear Congruential Generator or the lag for a Lagged
 ! Fibonacci Generator). The macroSPRNG_DEFAULT, defined in the SPRNG
 ! header files, can be used to choose the default parameters. If an invalid
 ! parameter is passed to this function, then a warning message is sent to
 ! stderr and the default parameter is used.     
 integer               :: sprng_param=SPRNG_DEFAULT
 
 ! public :: sprng_param, sprng_type
 ! public :: seed, streamnum, sprng_stream, nstreams
 public :: ranu, rang, seed
 public :: init_ran, free_ran, spawn_ran, write_chpseed, read_chpseed

contains

function ranu()
  ! Wraper to avoid to include sprng_f.h
  real(dp)                    :: ranu
  ranu=sprng(sprng_stream)
end function

subroutine init_ran(ltype,lstreamnum,lnstreams,lseed,lparam)
  ! FIXME: Sometimes I get segfault when I try to init_sprng twice (even after
  ! free_sprng) with different parameters. Sometimes I run aclocal again
  ! in the library and the error disappear. Some times I get the segfault in the
  ! case of spawn and not in the init_sprng.
  integer,optional,intent(in)  :: ltype
  integer,optional,intent(in)  :: lstreamnum
  integer,optional,intent(in)  :: lnstreams
  integer,optional,intent(in)  :: lseed
  integer,optional,intent(in)  :: lparam
 
  ! call werr("sprng already initializated",associated(sprng_stream))

  if(associated(sprng_stream)) then
    call free_ran()
  endif

  if(present(lseed))      seed        = lseed
  if(present(lstreamnum)) streamnum   = lstreamnum
  if(present(ltype))      sprng_type  = ltype
  if(present(lparam))     sprng_param = lparam
  if(present(lnstreams))  nstreams    = lnstreams
    
  ! This function retuns the ID of the stream when it completely successfully.
  ! It returns a NULL pointer and sends an error mesage to stderr in case of
  ! failure.
  stream = init_sprng(sprng_type,streamnum,nstreams,seed,sprng_param)
  sprng_stream=>stream
      
  ! Print random seed information. TODO: cut the integer in the used bites
  call wlog(''); write(logunit,'(a,i0)') "Random Seed: ", seed
  ! call wlog(''); write(logunit,*) "First random number: ",ranu(stream)
end subroutine
    
subroutine spawn_ran(lstreamnum,lnstreams)
  ! FIXME: Sometimes I get segfault when I try to init_sprng twice (even after
  ! free_sprng) with different parameters. Sometimes I run aclocal again
  ! in the library and the error disappear. Some times I get the segfault in the
  ! case of spawn and not in the init_sprng.
  integer,intent(in)          :: lstreamnum
  integer,optional,intent(in) :: lnstreams
  integer                     :: n,m

  ! Spawn (only once!)
  if(present(lnstreams)) then
    n=lnstreams-1

    call werr("sprng already spawned",allocated(streams))
    call werr("You should spawn more than 1 stream",n<1)
    
    allocate(streams(n))
    m=spawn_sprng(stream,n,streams)
    call werr("spawn failed",m/=n)

  endif
          
  n=lstreamnum-1

  call werr("Never spawned",.not.allocated(streams))
  call werr("Not such a stream",n>size(streams))

  ! Select the thread
  if(n==0) then
    sprng_stream=>stream
  else        
    sprng_stream=>streams(n)
  endif
end subroutine
  
subroutine free_ran()
  integer    :: junk,i
  junk=free_sprng(stream)
  if(.not.allocated(streams)) return
  do i=1,size(streams)
    junk=free_sprng(streams(i))
  enddo
end subroutine


subroutine  write_chpseed(chpunit)
  ! The programmer packs the state of the stream into an array and then saves it
  ! to a file. This state can later be retrieved by calling unpack_sprng, which is
  ! explained below. pack_sprng can also be used to pass a stream to another
  ! process. That process will unpack the packed array to obtain the stream.
  !
  ! Note: SPRNG does not free the memory associated with a stream when it packs
  ! it. If users do not plan to use the stream that has been packed, then they can
  ! explicitly call free_sprng in order to free the memory.
  integer,intent(in)            :: chpunit
  integer                       :: bsize
  character                     :: buffer(MAX_PACKED_LENGTH)
  character(MAX_PACKED_LENGTH)  :: actual_seed

  ! pack_sprng packs the state of the stream with IDstream into an array
  ! and returns the number of bytes actually required for the storage.
  ! Calls to this function involve memory allocation within this function.
  bsize = pack_sprng(sprng_stream,buffer)

  write(chpunit) buffer(1:bsize)

end subroutine write_chpseed


subroutine read_chpseed(chpunit)
  integer,intent(in)     :: chpunit
  integer                :: sz

  character              :: packed(MAX_PACKED_LENGTH)
  integer                :: i

  read(chpunit)   packed
  if(associated(sprng_stream)) i = free_sprng()
  sprng_stream = unpack_sprng(packed)
                                    
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
      v1=2.0_dp*sprng(sprng_stream)-1.0_dp
      v2=2.0_dp*sprng(sprng_stream)-1.0_dp
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
