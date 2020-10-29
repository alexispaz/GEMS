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

module gems_mpi
   
#ifdef HAVE_MPI  
use mpi_f08

type(mpi_comm)        :: mpi_world
type(mpi_datatype)    :: mpi_dp, mpi_ie1

! Passing MPI_STATUS_IGNORE for the status argument causes IBM® PE MPI to skip
! filling in the status fields. By passing this value for status, you can avoid
! having to allocate a status object in programs that do not need to examine the
! status fields.
type(mpi_status)      :: mpi_st

 
#endif

logical,public         :: mpi_on=.false.
integer,target,public  :: mpi_pc=0, mpi_tpc=1
integer,public         :: mpi_err
character(6),public    :: ch_mpi_pc='',ch_mpi_tpc=''

#ifdef HAVE_MPI  
contains

subroutine gems_mpi_init
use gems_constants, only:dp
use gems_strings, only:int2char0

! Initialize the MPI execution environment
call mpi_init()

! Determines the size of the group associated with a communicator
call mpi_comm_rank(MPI_COMM_WORLD, mpi_pc)

! Determines the rank of the calling process in the communicator
call mpi_comm_size(MPI_COMM_WORLD, mpi_tpc)

! set mpi mode
mpi_on=.true.

! set ch_mpi_pc for future string manipulation
write(ch_mpi_pc,'(a)') trim(int2char0(mpi_pc,mpi_tpc))
     
! Create types
call mpi_type_create_f90_real(dp,MPI_UNDEFINED,mpi_dp)
call mpi_type_create_f90_integer(1,mpi_ie1)

! Set the world
mpi_world = MPI_COMM_WORLD

end subroutine gems_mpi_init
#endif

end module gems_mpi

