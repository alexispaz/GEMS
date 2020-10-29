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

 
module gems_strings
  use gems_constants, only:dp

  implicit none

  private

  interface operator(.ich.)
    module procedure int2char,float2char
  end interface

  public    :: upcase,locase,operator(.ich.),int2char0,int2char
  public    :: contar_cifras

  character(len=*),parameter,public ::    &
    chset_u="ABCDEFGHIJKLMNOPQRSTUVWXYZ", & ! upper
    chset_l="abcdefghijklmnopqrstuvwxyz", & ! lower
    chset_n="0123456789",                 & ! number
    chset_ln=chset_l//chset_n,            & ! lower and numbers
    chset_lun=chset_l//chset_u//chset_n,  & ! upper lower and numbers
    chset_var=chset_ln//"_[]:"             ! allowed in a variable name

 contains
 
  function contar_cifras(num) result(c)
   integer,intent(in)          :: num
   integer                     :: a,c
   a=num
   c=0
   if(a==0) then
     c=1
     return
   endif
   do while (a>0)
     a=int(a/10)
     c=c+1
   enddo
  end function  
 
  function int2char(num)
   ! Convert integer to string
   integer,intent(in)           :: num
   character(20)                :: int2char
   write(int2char,'(I0)') num 
  end function  

  function float2char(num)
   ! Convert float to string
   real(dp),intent(in)          :: num
   character(20)                :: float2char
   write(float2char,'(f10.6)') num 
  end function   

  subroutine upcase(word)
  ! Change the word to upercase
  character(len=*), intent(inout) :: word
  integer :: i,k

  do i=1,len(word)
    k=index(chset_l,word(i:i))
    if (k .ne. 0) word(i:i)=chset_u(k:k)
  end do

  end subroutine upcase

  subroutine locase(word)
  ! Change the word to lowercase
  character(len=*), intent(inout) :: word
  integer :: i,k

  do i=1,len(word)
    k=index(chset_u,word(i:i))
    if (k .ne. 0) word(i:i)=chset_l(k:k)
  end do
  end subroutine locase
 

  function int2char0(num1,num2) result(c)
   ! Like integer_char but add ceros on the left until 
   ! fill the same digit number of num2.
   integer,intent(in)          :: num1,num2
   integer                     :: a,i
   character(20)               :: c
   a=contar_cifras(num2)-contar_cifras(num1)
   c = ''
   do i=1,a
     c(i:i) = '0' 
   enddo
   write(c,'(a)') trim(c) // trim(.ich.num1)
  end function  
 


end module gems_strings 


