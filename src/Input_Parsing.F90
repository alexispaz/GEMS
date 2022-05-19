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
! ---
! 
! This file incorporates work derived from the following codes:
!
! Fortran input module <http://www-stone.ch.cam.ac.uk/programs.html>
! by Anthony J. Stone, covered by the following copyright and permission notice:  
!
! Copyright (c) 2005  Anthony J. Stone
! 
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!  .
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  .
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to
!  the Free Software Foundation, Inc., 51 Franklin Street,
!  Fifth Floor, Boston, MA 02110-1301, USA, or see
!  http://www.gnu.org/copyleft/gpl.html
!
! Code posted on stackoverflow <https://stackoverflow.com/a/34752340/1342186>
! by IanH <https://stackoverflow.com/users/1234550/ianh> 
! used under CC BY-SA 3.0.
 
module gems_input_parsing    
use gems_constants
use gems_algebra, only: xyz_polares
use gems_strings, only: locase, upcase
use gems_variables, only: polvars
use gems_errors

implicit none
private
 
 
! Variables auxiliares para lectura, parsing y demas. 
character(:),allocatable :: w1,w2
integer                  :: i1,i2

! The width
integer      :: lrecl = linewidth

! This is to allow all the modules to execute commands declared by
! CLinterpreter module, regardless if CLinterpreter also use those modules.
procedure(iface_exe),pointer   :: exec
interface
  recursive subroutine iface_exe(com)
  character(*)  :: com 
  end subroutine
end interface
  
public :: execute_block, exec

type input_options

  integer :: or=6
  integer :: ir=5


  logical :: clear=.true.,    &
             skipbl=.true.,  &
             echo=.false.
  
  integer :: nerror=0

  character(1) :: comment='#'
  character(1) :: concat="\"

  contains
    procedure   :: init=>input_options_init
end type

! The input options for gems input
type(input_options), public, target     :: gems_iopts
                                     
! Point to the input options to be used in the module
! TODO: Use instead a structure?
! The deffault should be gems_iopts
type(input_options), public, pointer   :: opts=>null()

! El instructivo original
! 
!  free-format input routines
!     call read_line(eof[,inunit,or])
!  read next input record from unit ir into buffer, and analyse for data
!  items, null items (enclosed by commas), and comment (enclosed in
!  parentheses, which may be nested, and occurring where a new item may
!  occur). data items are terminated by space or comma, unless enclosed
!  in single or double quotes.
!  if the optional argument inunit is provided, the record is read from
!  there instead.
!  if a line ends with the concatenation sequence (default "+++")
!  possibly followed by spaces, the sequence and any following spaces
!  are removed and the next line concatenated in its place.
!  this is repeated if necessary until a line is read that does not end
!  with the concatenation sequence.
!  the logical variable eof is set true if end of file is encountered,
!  otherwise false.
!  the public module variable item is set to zero and nitems to the number
!  of items in the record.

!     call parse_items(string)
!  parse the string provided in the same way as the read_line routine
!  (which uses this routine itself) and leave the details in the buffer
!  as for a line read from the input. input directives are not
!  interpreted.

!     call readx(v)
!  read an item of type x from the buffer into variable v:
!     call readf   single or real(dp), depending on the type of v
!     call readi   integer
!     call reada   character string
!     call readu   character string, uppercased
!     call readl   character string, lowercased
!     call readb   booleano (add by alexis paz)
!  all of these return null or zero if there are no more items on the
!  current line. otherwise item is incremented so that it gives the
!  number of the item that has just been read.
!     call readf(v,factor)
!  if the optional second argument is supplied to readf, the variable v
!  is divided (multipli!!! alexis.) by it after it has been read. this is convenient for converting
!  from external to internal units.
!     call reread(k)
!  k>0  prepare to read item k
!  k<0  go back |k| items
!  k=0  same as k=-1, i.e. reread last item.

!     call getx
!  same as the corresponding readx, but a new record is read if there are
!  no more items in the current one.

!     call read_colour(fmt,col,clamp)
!  read a colour definition, in a form specified by fmt:
!  fmt="grey": read a single number between 0 and 1
!  fmt="rgb": read 3 numbers, which are rgb colour values between 0 and 1
!  fmt="rgb255": read 3 numbers, which are rgb colour values between 0 and 255
!  fmt="rgbx": read a single 6-character string giving the rgb colour
!        values in hexadecimal.
!  in the last two cases the colour values are scaled to the range from 0 to 1.
!  clamp is an optional logical argument. if present and true, the colour
!        values are clamped (after scaling, if appropriate) to the range 0 to 1.

!  item    is the number of the last item read from the buffer

!  nitems  is the number of items in the buffer

!  char(i:i) is the ith character in the buffer

!  loc(i)  is the starting position of the ith item
!  end(i)  is the position of the last character of the ith item

!  line    is the number of lines (records) read

!  if skipbl is set to true, lines containing no items other than
!  comment are skipped automatically, so that input will always return
!  a line with at least one item (possibly null).  if skipbl is false
!  (default) no skipping occurs and the next data line is returned
!  regardless of what it contains.

!  if clear is set to true (default) then null items will be returned 
!  as zero or blank. if an attempt is made to read more than
!  nitems items from a line, the items are treated as null. if
!  clear is false, a variable into which a null item is read is
!  left unaltered.

!  nerror specifies the treatment of errors found while reading numbers:
!    0  hard error - print message and stop (default).
!    1  soft error - print message and return zero.
!    2  soft error - no message, return zero, set nerror to -1.
!  if nerror is set to 2 it is important to test and reset it after reading
!  a number.

!  ir is the input stream from which the data are to be read (default 5).

!  if echo is true, the input line will be reflected to standard output.

!  last gives the position of the last non-blank character in the line.
            


!#ifdef nagf95
!use f90_unix_env, only: getarg
!#endif


! Variables relacionadas con los nombres del input output log etc.
character(:),allocatable   :: ioprefix
character(:),allocatable   :: stdin,logfile
integer                :: stdinlen
public                 :: ioprefix,stdin,logfile,eof

logical                :: eof=.false.

character(:), allocatable, save :: char

integer, save :: item=0, nitems=0, loc(0:80)=0, end(80)=0,               &
    line(0:10)=0, level=0, unit(0:10)

character, parameter :: space = " ", squote = "'",             &
    comma = ",",  dquote = '"', tab=achar(9),                  &
    plus="+", minus="-", dot="."

character(linewidth) :: file(10)=''

integer, save :: lc=3
 
 
interface readb
  module procedure readb, readb_vector
end interface

public :: item, nitems, stdinlen,read_line, stream, reada,read_magnitud,read_energy, readu, readl,        &
    readf, get, readi, getf, geti, reread, input_options, readblock,     &
    read_colour, try_get,                 &
    parse_items, readb, readelement, readia
    public :: line,level,char,end,loc


interface get
  module procedure readf, readi
endinterface
 
interface try_get
  module procedure try_readf, try_readi
endinterface
           

!                                                                 block control
!------------------------------------------------------------------------------

#define _NODE cmdline_v
#define _TYPE character(len=linewidth)
!   real, parameter :: growth_rate = 1.1
#include "vector_header.inc"

type,extends(cmdline_v) :: cmd_v
  contains
  procedure :: find => cmd_v_find_first
end type
 
! vector de bloques de comando para interpretar despues
type(cmd_v)  :: bloques(10),bloque

! This variable are related to control the execute commands block inside a
! subprogram.
logical,public   :: b_load=.false.,load_silent=.false.
integer,public   :: load_each, load_blk

public  :: bloques

integer          :: in

contains
 
#define _NODE cmdline_v
#define _TYPE character(len=linewidth)
#define _DUMMY character(*)
#include "vector_body.inc"
               
! TODO: send this routine to vectr_body
function cmd_v_find_first( vec, data ) result(i)
! Return the first index of "vec" that match with the "indx" input.
! If there is not matching return 0
class(cmd_v),intent(in)  :: vec
character(*),intent(in)  :: data
integer            :: i

do i = 1,vec%size
  if (trim(vec%o(i)) == trim(data)) return
enddo
i=0

end function cmd_v_find_first
  
subroutine read_line(eof,inunit)

!  read next input record from unit ir into buffer, and analyse for data
!  items, null items (enclosed by commas), and comment (enclosed in
!  parentheses, which may be nested, and occurring where spaces may
!  occur).
!  if the optional argument inunit is specified, read a line from that
!  unit instead of ir.

!  stream-switching commands may occur in the data:
!      #include file-name
!         switches input to be from the file specified;
!      #revert
!         (or end-of-file) reverts to the previous file.
!  however it does not make sense for stream-switching commands to
!  be used when the input unit is specified explicitly. if they are given
!  in this case, they apply to the default input stream.

!  '\' is the line concatenation string: if it is found as the last
!  non-blank character string on a line, the following line is read and
!  appended to the current line, overwriting the concatenation string.
!  this procedure is repeated if necessary until a line is found that does
!  not end with the concatenation string.

logical, intent(out) :: eof
integer, intent(in), optional :: inunit

integer :: l, m, flag, iostat 

character(:), allocatable :: aux, iomsg, aux2

iomsg=''
eof=.false.
if (present(inunit)) then
  in=inunit
else
  in=opts%ir
end if
 

lines: do
   
  ! Get line
  line(level)=line(level)+1
  call get_line(in, aux, iostat, iomsg)
  aux=trim(aux)
  char=aux
  
  ! Handle EOF
  if(is_iostat_end(iostat)) then
    if (level > 0) then
      !  revert to previous input
      close(in)
      level=level-1
      opts%ir=unit(level)
      in=opts%ir
      cycle lines
    else
      !  end of input
      char=''
      item=0
      nitems=-1
      eof=.true.
      return
    endif
  endif
   
  ! Echo
  if (opts%echo) write(opts%or,"(a)") aux
  
  ! Empty line
  m=len(char)
  if(m==0) cycle
         
  ! Concatenate with next line/s
  do while (char(m:m)==opts%concat)
    call get_line(in, aux, iostat, iomsg)
    aux=trim(aux)
    if (opts%echo) write(opts%or,"(a)") aux
    char=char(:m-1)//adjustl(aux)
    m=len(char)
  enddo
 
  ! Replace tab by single space
  do while (index(char,tab) > 0)
    l=index(char,tab)
    char(l:l)=space
  end do                 

  ! Remove leading spaces
  char=adjustl(char)
                
  ! Search for special directives
  flag=parse_special(char)
  if(flag==0) cycle lines
         
  ! Skip commentaries
  if(char(1:1)==opts%comment) cycle lines
                      
  ! Parsing asignaciones de variables
  l=index(char,':=')
  if(l/=0) then
 
    ! Las variables que empiezan con underscore son solo internas
    ! no para que definan o accedan los usuarios
    call werr('Variable name should not start with underscore',adjustl(char(:1))=='_')
 
    aux=trim(adjustl(char(:l-1)))
    aux2=parse_expand(trim(adjustl(char(l+2:))))
    call polvars%save(aux,aux2)
    cycle lines
  else 
    l=index(char,'=')
    if(l/=0) then
 
      ! Las variables que empiezan con underscore son solo internas
      ! no para que definan o accedan los usuarios
      call werr('Variable name should not start with underscore',adjustl(char(:1))=='_')
            
      aux=trim(adjustl(char(:l-1)))
      aux2=trim(adjustl(char(l+1:)))
      call polvars%save(aux,aux2)
      cycle lines
    endif
  endif
       
  ! Expando la linea
  char=parse_expand(char)
    
  ! Parsing words
  call parse_items(char)

  !  blank except for comment?
  if (nitems == 0 .and. opts%skipbl) then
    cycle lines   !  read another line
  end if
                     
  ! Parsing blocks
  flag=parse_blocks()
  if(flag==0) cycle lines

  exit lines    !  finished

end do lines

end subroutine read_line

!-----------------------------------------------------------------------

subroutine parse_items(string)
! Split a line into words (items)
character(len=*),intent(inout) :: string
integer :: l,i,j

l=0
item=0
nitems=0
char=string
do
  
  ! Busqueda del comienzo de un item
  i=verify(string(l+1:),space//tab)
  if(i==0) exit
  i=i+l

  select case(string(i:i))
  case(',') 

    ! Establezco siguiente posicion
    l=i+1

    ! Null item
    i=0
    j=0

  case("'")
    
    ! Busco el cierre
    j=scan(string(i+1:),"'")
    call werr('missing quote',j==0)
    j=j+i

    ! Establezco siguiente posicion
    l=j+1

    ! Le resto las comillas
    i=i-1
    j=j-1

  case default  

    ! Busqueda del final de un item
    j=scan(string(i+1:),space//tab//comma)
    j=j+i-1
       
    ! Establezco siguiente posicion
    l=j+1
       
  end select

  ! Guardo el item
  nitems=nitems+1
  loc(nitems)=i
  end(nitems)=j
 
  ! print *, nitems,i,j,string(i:j)
end do

end subroutine parse_items

function parse_blocks() result(flag)
use gems_strings
integer              :: k
integer              :: flag
integer              :: repini,repfin,repstep
logical              :: reverse
character(:), allocatable :: aux,w


flag=1

! Get first item (without change item counter)
w=char(loc(1):end(1))

! Parsing bloques
if(w/='bloque') return

! If it is a block, it doesn't matters to change item counter
call readl(w) ! Already known from above w assignment line
flag=0

! Get block kind (without change item counter)
call readl(w) 
select case(w)
case('load')
  call readi(k)
  call execute_block(bloques(k),1,1,1,.false.)
  return
case('load_silent')
  call readb(load_silent)
  return 
case('load_each')
  b_load=.true.
  call readi(load_each)
  call readi(load_blk)
  return 
case('load_off')
  b_load=.false.
  return  

case('repeat')

  ! Get iteration numbers
  call readi(i1)
  if(item<nitems) then
    call readi(i2)
  else
    i2=i1
    i1=1
  endif

   ! Taking the arguments 
  repstep=1
  repini=min(i1,i2)
  repfin=max(i1,i2)
  if(i1>i2) then
    reverse=.true.
    aux=trim(int2char(repfin))
    call polvars%hard('ci',aux)
    call polvars%hard('i',repfin)
  else
    reverse=.false.
    aux=trim(int2char(repfin))
    call polvars%hard('ci',aux)
    aux=trim(int2char(repini))
    call polvars%hard('i',aux)
  endif
  if(item<nitems) call readi(repstep)
 
  ! Read the block
  call readblock(bloque)
   
  ! Execute the block
  call execute_block(bloque,repini,repfin,repstep,reverse)
 
case('save')

  call readi(i1)  ! 1 of 10
  call readblock(bloques(i1))

case default  
  call wwan('I do not understand the last command')  

endselect      
 

end function parse_blocks
     
function parse_math(string) result(ans)
  use fparser, only: initf, destf, parsef, evalf!, EvalErrMsg, EvalErrType
  ! Reemplaza la strings por el valor numerico
  character(*),intent(inout) :: string
  character(*), dimension(0),  parameter :: var  = "a"
  real(dp), dimension(0),  parameter :: val  = 0._dp
  character(linewidth)       :: ans
  real(dp)                   :: res

  ! Run the calculator plugin.

  ! Use bc (benchmark: 13.1 sec with 5000 times of 's(0.5)*10^5')
  ! bc don't hand exponential notation
  !call system('echo '''//trim(w1)//''' | bc -l >  gems'//trim(ch_mpi_pc)//'.ans')

  ! Initialize function parser for 1 simultaneous evaluation
  call initf (1)

  ! Parse and bytecompile ith function string 
  call parsef(1, string, var)

  ! Interprete bytecode representation of ith function
  res = evalf(1, val)

  write(ans,*) res
 
  ! Initialize function parser for 1 simultaneous evaluation
  call destf()

end function parse_math

function parse_vars(string) result(ans)
  use gems_random, only:ranu
  use gems_strings, only: chset_var
  ! Toma un string que puede contener variables y reemplaza las mismas por su
  ! valor
  character(*) :: string       ! string with variables to expand
  character(linewidth) :: med  ! value of an expanded variable 
  character(linewidth) :: ans  ! expanded string
  integer :: l                 ! position in string after an expansion take place
  integer :: i, j              ! Start and end position of the variable to expand
  integer :: k, m              ! Auxiliary integer
                                
  ans=string

  l=1
  do

    ! El nombre de la variable no incluye signos $
    ! e.g. $pepe tiene como nombre pepe
    ! e.g. $pepe$ tiene como nombre pepe
    ! e.g. $pe[3] tiene como nombre pe[3]

    ! Busco el primer caracter del nombre
    i=scan(ans(l:),'$')
    if (i==0) return
    i=l+i
         
    ! Busco el ultimo caracter del nombre
    j=verify(ans(i:),chset_var)
    if (j==0) j=len(ans(i:))
    j=i+j-2
                
    ! Check for array variable syntax (like a[2] or a[pepe])
    ! TODO allow multiple array (like a[2,pepe])
    k=scan(ans(i:j),'[') ! Inicio array
    if(k/=0) then
      m=i+k
      call werr('unexpected [ symbol',scan(ans(m:j),'[')/=0)
      k=scan(ans(m:j),']')
      call werr('Missing ]',k==0)
      call werr('Unexpected character after ]',m+k-1/=j)
    else
      call werr('unexpected ] symbol',scan(ans(i:j),']')/=0)
    endif
  
    ! Reemplazo por el valor si se encuentra definida la variable
    if(ans(i:j)=="rnd") then
      write(med,fmt='(e13.6)') ranu()
    else
      med=polvars%expand(ans(i:j))
    endif
    med=adjustl(med)
 
    ! Compute the next position after expansion
    j=j+1
    if(ans(j:j)=='$') j=j+1

    ! Escribo la cadena
    ans=ans(:i-2)//trim(med)//ans(j:)

    ! Next position by expansion size
    l=i-1+len(trim(med))
                         
  enddo
 
end function parse_vars

function parse_special(w) result(flag)
  ! It start reading from the begining, and
  ! after find a special symbol it search its end.
  character(*),intent(in)    :: w
  character(linewidth)       :: ans
  integer :: flag, fail

  flag=1

  if(len(w)==0) then
    flag=0
    return
  endif

  ! To let w intent(in)
  ans=adjustl(w)

  if (w(1:1) == "#") then

    flag=0

    ! Just in case opts%comment is #
    if(w(1:1)==opts%comment) then
      if(len(w)==1)  return
      if(w(2:2) == " ") return
    endif

    ! Le saco los espacios
    ans=adjustl(ans(2:))

    ! Search for words
    call parse_items(ans)

    call readu(w1)
    select case(w1)
    case("INCLUDE")   ! take input from specified file

      call reada(w2)
      
      if (level == 0) unit(0)=opts%ir
      level=level+1
      line(level)=0
      opts%ir=find_io(91)
      unit(level)=opts%ir

      open(unit=opts%ir,file=trim(w2),status="old",iostat=fail)
      call werr(trim(w2)//" could not be opened",fail/=0)
      in=opts%ir

      file(level)=trim(w2)

    case("CONCAT")

      call reada(w2)
      opts%concat=trim(w2)
      lc=len(trim(opts%concat))

    case("REVERT")
      close(in)
      file(level)=''
      level=level-1
      opts%ir=unit(level)
      in=opts%ir

    case("WIDTH")

      call reada(w2)
      read (unit=w2,fmt="(i6)") lrecl

    case("ECHO")

      opts%echo=.true.
      call readl(w2)
      if (w2 == "off") opts%echo=.false.

    case default
      call werr("unrecognized directive "//trim(w)//"in input",opts%comment/='#')
    end select

  endif
   
end function parse_special

function parse_expand(string) result(ans)
  ! It start reading from the begining, and
  ! after find a special symbol it search its end.
  character(*),intent(in)    :: string
  character(linewidth)       :: med,ans
  integer :: l, f, i, j, k, m 

  ans=adjustl(string)
  l=1

  ! No needed, skiped by read_line
  ! if(ans(1:1)==opts%comment) return

  do

    !Busco primer caracter del contenido (incluyendo el parentesis)
    i=scan(ans(l:),"'{(")
    if (i==0) exit
    i=i+l-1

    select case(ans(i:i))
    case("'")

      ! Busco el cierre
      j = scan(ans(i+1:),"'")
      call werr('missing quote',j==0)
      j=j+i

      ! Establezco siguiente posicion
      l = j+1

    case("{")

      ! Busco el cierre
      j=scan(ans(i+1:),"}")
      call werr('missing }',j==0) 
      j=j+i

      ! Proceso en busqueda de variables
      med=parse_vars(ans(i+1:j-1))
      med=parse_math(med)

      ! Parsing del contenido
      ans=ans(:i-1)//trim(med)//ans(j+1:)
  
      ! Establezco siguiente posicion
      l=i+len(trim(med))
 
    case("(")

      ! Busco el cierre
      j=scan(ans(i+1:),")")
      call werr('missing )',j==0) 
      j=j+i

      ! Busco aperturas internas por si estan anidados
      ! Esto permite usar parentesis en los comentarios:
      ! (bla bla bla (e.g. asd) bla bla )
      k=0
      m=scan(ans(i+1:j-1),"(")
      f=i+m
      do while(m/=0)
        k=k+1
        m=scan(ans(f+1:j-1),"(")
        f=f+m
      end do

      ! Busco cierre considerando anidados.
      do f=1,k 
        m=scan(ans(j+1:),")")
        j=j+m
        ! Si no hay tantos cierres ignoro el inconveniente (total son comentarios).
        ! call werr('mismatch parenthesis',m==0)
      end do

      ! Parsing del contenido
      ans=ans(:i-1)//ans(j+1:)
 
      ! Establezco siguiente posicion
      l=i
                             
    end select

  end do


  ans=trim(parse_vars(ans))

  call werr('desbordamiento',len(ans)>linewidth)  ! FIX

end function parse_expand

subroutine input_options_init(opts,clear_if_null,skip_blank_lines, echo_lines, error_flag, concat_string, comments, out, in)
class(input_options)          :: opts
integer, intent(in), optional :: out,in
logical, optional             ::  clear_if_null, skip_blank_lines, echo_lines
integer, optional             :: error_flag
character(len=*), optional    :: concat_string
character(1), optional        :: comments

if (present(comments))         opts%comment = comments
if (present(clear_if_null))    opts%clear =clear_if_null
if (present(skip_blank_lines)) opts%skipbl=skip_blank_lines
if (present(echo_lines))       opts%echo  =echo_lines
if (present(error_flag))       opts%nerror=error_flag
if (present(out))              opts%or=out
if (present(in))               opts%ir=in

if (present(concat_string)) then
  call werr("concatenation string must be 8 characters or fewer",len(trim(concat_string)) > 8)
  opts%concat=concat_string
  lc=len(trim(concat_string))
endif
 
end subroutine input_options_init

subroutine stream(n)
!  set the input stream for subsequent data to be n.
integer, intent(in) :: n
opts%ir=n
end subroutine stream

! Read
!-----

subroutine reada(m)
! Copy characters from the next item into `m`.
! if the first character is a single or double quote, the string is
! terminated by a matching quote and the quotes are removed.
character(:),allocatable,intent(out) :: m

if (opts%clear) m=''
!  if there are no more items on the line, m is unchanged
! if (item >= nitems) return
call werr('Expected string',item>=nitems)

item=item+1

!  null item?
if (loc(item) == 0) return

m=char(loc(item):end(item))

end subroutine reada

subroutine readf_sbrace(f,m)
! lee un escrito de la forma flotante[string] y devuelve flotante en a y string
! en m. Util para las unidades.
character(:),allocatable, intent(out)  :: m
integer :: l,j,i
real(dp) :: f

if (opts%clear) m=''
!  if there are no more it/sems on the line, m is unchanged
! if (item >= nitems) return
call werr('Expected array',item>=nitems)

item=item+1
!  null item?
if (loc(item) == 0) return

l=loc(item)
if (char(l:l) == squote .or. char(l:l) == dquote) then
  m=char(l+1:end(item)-1)
else
  m=char(l:end(item))
endif

j=end(item)
end(item)=l+scan(m,'[')-2

if(end(item)/=0) then

  item=item-1
  call readf(f)

  i=loc(item)
  loc(item)=end(item)+2
  end(item)=j-1
  item=item-1
  call readl(m)

  loc(item)=i
  end(item)=j
else
  end(item)=j
  call readf(f)
  m='ui'
endif 


end subroutine readf_sbrace
 
subroutine read_energy(a,flag,factor)
! lee una magnitud del tipo 5[eV] y convierte a su equivalente en ui. Si [..] no se
! especifica asume que es el valor en ui.
integer,optional                :: flag
real(dp), intent(inout) :: a
real(dp), intent(in), optional :: factor
character(:),allocatable  :: m

if (present(flag)) flag=0
call readf_sbrace(a,m)
          
flag=int(factor) ! FIXME this variable is here only for compatibility, is not used

select case(trim(m))
case('kt')
  a=a*kt
case('ev')
  a=a*ev_ui
case('j')
  a=a*joule_ui
case default
  call reread(-1)
  if (present(flag)) flag=1  ! No es una energia valida
end select

end subroutine read_energy

subroutine read_magnitud(a,factor,flag)
! lee una magnitud del tipo 5[eV] y convierte a su equivalente en ui. Si [..] no se
! especifica asume que es el valor en ui.
integer,optional                :: flag
real(dp), intent(inout) :: a
real(dp), intent(in), optional :: factor

flag=int(factor) ! FIXME this variable is here only for compatibility, is not used

flag=0
call read_energy(a,flag)
call read_energy(a,flag)
call read_energy(a,flag)

!  call reread(-1)
!  flag=1  ! No es una magnitud valida

end subroutine read_magnitud
    
subroutine readb(b)
! subrutina agregada por alexis paz
!  read an integer from the current record

logical, intent(inout) :: b
character(:),allocatable  :: string

if (opts%clear) b=.false.

!  if there are no more items on the line, i is unchanged
! if (item >= nitems) return
call werr('Expected loigcal',item>=nitems)

string=''
call reada(string)
!  if the item is null, i is unchanged
if (string == "") return
read (unit=string,fmt=*,err=99) b
return

99 b=.false.
select case(opts%nerror)
case(-1,0)
  call werr("error while reading bollean number",.true.)
case(1)
  write(opts%or,"(2a)") "error while reading bollean. input is ", trim(string)
case(2)
  opts%nerror=-1
end select

end subroutine readb

subroutine readb_vector(b)
! subrutina agregada por alexis paz
!  read an integer from the current record

logical, intent(inout) :: b(:)
integer           :: i
character(:),allocatable  :: string

if (opts%clear) b=.false.

do i = 1,size(b)

  !  if there are no more items on the line, i is unchanged
  ! if (item >= nitems) return
  call werr('Expected vector',item>=nitems)

  string=''
  call reada(string)
  !  if the item is null, i is unchanged
  if (string == "") return
  read (unit=string,fmt=*,err=99) b(i)

enddo

return

99 b=.false.
select case(opts%nerror)
case(-1,0)
  call werr("error while reading bollean number",.true.)
case(1)
  write(opts%or,"(2a)") "error while reading bollean. input is ", trim(string)
case(2)
  opts%nerror=-1
end select

end subroutine readb_vector
 
subroutine readelement(i)
use gems_elements, only: inq_z
! Lee un elemento atomico desde su nombre o numero atomico (i.e. 'he' or 2 ) y
! devuelve su numero atomico.  basicamente es la readi pero no se enoja si no es
! un entero.
integer, intent(inout) :: i
character(:),allocatable  :: string

if (opts%clear) i=0

!  if there are no more items on the line, i is unchanged
! if (item >= nitems) return
call werr('Expected element',item>=nitems)

string=''
call reada(string)
!  if the item is null, i is unchanged
if (string == "") return
read (unit=string,fmt=*,err=99) i
return

99 i=0
select case(opts%nerror)
case(-1,0)
  call werr("error while reading integer number",.true.)
case(1)
  i = inq_z(string)
  if(i==-1) write(opts%or,"(2a)") "error while reading atomic element. input is ", trim(string)
case(2)
  opts%nerror=-1
end select

end subroutine readelement
 
subroutine readia(i,m)
! Lee un input que a priori puede ser entero o string. 
! Si es entero devuelve en m un string nulo
! Si es string devuelve m y no cambia i
! Basicamente es el readi, pero no chilla si no es entero, sino que devuelve lo
! que lee
integer, intent(inout) :: i
character(:),allocatable, intent(out)  :: m

if (opts%clear) i=0

!  if there are no more items on the line, i is unchanged
! if (item >= nitems) return
call werr('Expected integer or string',item>=nitems)

m=''
call reada(m)
!  if the item is null, i is unchanged
if (m == "") return
read (unit=m,fmt=*,err=99) i
return

99 i=0
select case(opts%nerror)
case(-1,0)
  call werr("error while reading integer number",.true.)
case(1)
  return
case(2)
  opts%nerror=-1
end select

end subroutine readia
 
impure elemental function getf() result(a)
! Read the next item from the buffer as `real(dp)`.
! if the optional argument factor is present, the value read should be
! divided by it. (external value = factor*internal value)
use gems_strings, only: operator(.ich.)
real(dp)                  :: a
character(:),allocatable  :: string

if (opts%clear) a=0._dp

! Check end of items
call werr('Expected float',item>=nitems)

! Read number as string
string=''
call reada(string)
if (string == "") return

read (unit=string,fmt=*,err=99) a
return

99 a=0._dp
select case(opts%nerror)
case(-1,0)
  call werr("error while reading real number",.true.)
case(1)
  call werr("Expected float but got: "//string,.true.)
case(2)
  opts%nerror=-1
end select

end function getf

impure elemental function geti() result(i)
! Read the next item from the buffer as `integer`.
use gems_strings, only: operator(.ich.)
integer                   :: i
real(dp)                  :: f
character(:),allocatable  :: string

if (opts%clear) i=0

! Check end of items
call werr('Expected integer',item>=nitems)

! Read number as string
string=''
call reada(string)
if (string == "") return

! This way to read an integer allow for exponential notation like "1e10"
read (unit=string,fmt=*,err=99) f
call werr('Integer beyond kind boundaries: +-'//.ich.huge(i),huge(i)<abs(f))
i=int(f)
call werr('Expected integer but got float.',abs(f-i)>epsilon(f))

return

99 i=0
select case(opts%nerror)
case(-1,0)
  call werr("error while reading integer number",.true.)
case(1)
  call werr("Expected integer but got: "//string,.true.)
case(2)
  opts%nerror=-1
end select

end function geti
     
subroutine readu(m)
character(:),allocatable  :: m

call reada(m)
call upcase(m)

end subroutine readu

subroutine readl(m)
character(:),allocatable  :: m

call reada(m)
call locase(m)

end subroutine readl

! Read if possible
!-----------------
 
impure elemental subroutine try_readf(a,factor)
! Read the next item from the buffer as a real (real(dp)) number.
! if the optional argument factor is present, the value read should be
! divided by it. (external value = factor*internal value)
use gems_errors, only: silent, errf
real(dp), intent(inout)        :: a
real(dp)                       :: b
real(dp), intent(in), optional :: factor

! Save default value
b=a
silent=.true.
if (present(factor)) then
  call readf(a,factor) 
else
  call readf(a) 
endif       
if(errf) a=b
silent=.false.
 
end subroutine try_readf
  
impure elemental subroutine try_readi(i)
! Read the next item from the buffer as a real (real(dp)) number.
! if the optional argument factor is present, the value read should be
! divided by it. (external value = factor*internal value)
use gems_errors, only: silent, errf
integer, intent(inout)         :: i
integer                        :: j

! Save default value
j=i
silent=.true.
call readi(i) 
if(errf) i=j
silent=.false.
 
end subroutine try_readi
 

!-----------------------------------------------------------------------

subroutine reread(k)

integer, intent(in) :: k
!  k>0  reread from item k
!  k<0  go back |k| items
!  k=0  same as k=-1, i.e. reread last item.

if (k < 0) then
  item=item+k
else if (k == 0) then
  item=item-1
else
  item=k-1
endif
if (item < 0) item=0

end subroutine reread

!----------------------------------------------------------------

subroutine read_colour(fmt, colour, clamp)
character(len=*), intent(in)  :: fmt
real(dp), intent(out)         :: colour(3)
logical, intent(in), optional :: clamp
character(:),allocatable  :: x
integer      :: i, r, g, b
real(dp)     :: c

select case(fmt)
case("grey","gray")
  call readf(c)
  colour(:)=c
case("rgb")
  call readf(colour(1))
  call readf(colour(2))
  call readf(colour(3))
case("rgb255")
  call readf(colour(1),255._dp)
  call readf(colour(2),255._dp)
  call readf(colour(3),255._dp)
case("rgbx")
  call readu(x)
  read (x(1:2),"(z2)") r
  colour(1)=r/255_dp
  read (x(3:4),"(z2)") g
  colour(2)=g/255_dp
  read (x(5:6),"(z2)") b
  colour(3)=b/255_dp
case default
  call werr('colour keyword not recognised',.true.)
end select

if (present(clamp)) then
  if (clamp) then
    do i=1,3
      if (colour(i)>1d0) colour(i)=1d0
      if (colour(i)<0._dp) colour(i)=0._dp
    end do
  end if
end if

end subroutine read_colour


!----------------------------------------------------------------

subroutine readblock(bloque)
    ! Asumiendo que se leyo alguna palabra que detecto un comienzo de bloque, se
    ! puede leer tod el bloque llamando a esta subrutina. Esta va a guardar el
    ! encabezado
    type(cmd_v)         :: bloque
    character(len=linewidth) :: w

    ! Inicializo el bloque si no estaba inicializado
    call bloque%init()

    ! Guardo el encabezado en el elemento 1
    call bloque%append(char) 

    ! Guardo las lineas del bloque mediante el guardado de la variable noparsed.
    ! TODO: Esto por ahora no permite includes dentro del bloque. 
    do
      read (opts%ir,"(a)") w
      if (trim(w)=='fin') exit 
      if (opts%echo) write(opts%or,"(a)") trim(w)
      call bloque%append(w)
    enddo

    if (opts%echo) write(opts%or,"(a)") trim(w)
 
end subroutine readblock
     
subroutine execute_block(bloque,repini,repfin,repstep,reverse)
use gems_strings
  class(cmd_v),intent(in)   :: bloque
  integer,intent(in)        :: repini,repfin,repstep
  logical,intent(in)        :: reverse
  integer                   :: i,j,l
  character(:),allocatable  :: w,aux

  ! Proforming the repeat loop
  outer: do j = repini,repfin,repstep
  
    if(reverse) then
      i=repfin-(j-repini)
      aux=trim(int2char0(i,repfin))
      call polvars%hard('ci',aux)
      call polvars%hard('i',i)
    else
      aux=trim(int2char0(j,repfin))
      call polvars%hard('ci',aux)
      call polvars%hard('i', j)
    endif 
  
    if (repfin-repini>0) then
      call wstd(); write(logunit,*)  '>>> iterando: $i$=' // trim(adjustl(polvars%expand('ci')))
    endif
 
    do i = 2,bloque%size

      w=trim(adjustl(bloque%o(i)))

      ! Parsing comentarios
      if(w(1:1)==opts%comment) cycle

      ! Parsing asignaciones de variables
      l=index(w,":=")
      if(l/=0) then
        call polvars%save(trim(adjustl(w(:l-1))),parse_expand(trim(adjustl(w(l+2:)))))
        cycle
      else 
        l=index(w,"=")
        if(l/=0) then
          call polvars%save(trim(adjustl(w(:l-1))),trim(adjustl(w(l+1:))))
          cycle
        endif
      endif
           
      ! Expando la linea
      w=parse_expand(w)
        
      ! Parsing words
      call parse_items(w)
                  
      !  blank except for comment?
      if (nitems == 0 .and. opts%skipbl) then
        cycle 
      end if
              
      call readl(w1)

      select case(w1)
      case('cycle')
        call wstd('# > cycled signal captured')
        exit
      case('exit')
        call wstd('# > exit signal captured')
        exit outer
      ! case('term')
      !   call wstd('repeat: termt signal captured')
      !   term_signal=.true.
      !   exit outer
      endselect

      call exec(w1)

    enddo        

  enddo outer
  
endsubroutine execute_block

subroutine get_line(lun, line, iostat, iomsg)
! Derived from https://stackoverflow.com/a/34752340/1342186 
! code by IanH <https://stackoverflow.com/users/1234550/ianh> 
! used under CC BY-SA 3.0.
integer, intent(in)           :: lun
character(:), intent(out), allocatable :: line
integer, intent(out)          :: iostat
character(*), intent(inout)   :: iomsg

integer, parameter            :: buffer_len = 80
character(len=buffer_len)     :: buffer
integer                       :: size_read

line = ''
do
  read ( lun, '(A)',    &
      iostat = iostat,  &
      iomsg = iomsg,    &
      advance = 'no',   &
      size = size_read ) buffer
  if (is_iostat_eor(iostat)) then
    line = line // buffer(:size_read)
    iostat = 0
    exit
  else if (iostat == 0) then
    line = line // buffer
  else
    exit
  end if
end do

end subroutine get_line

! Wrappers
! -------
 
impure elemental subroutine readf(a,factor)
! Wrap getf into subrroutine to create a generic interface
real(dp), intent(inout)        :: a
real(dp), intent(in), optional :: factor
a=getf()
if (present(factor)) a=a*factor
end subroutine readf
     
impure elemental subroutine readi(i)
! Wrap geti into subrroutine to create a generic interface
use gems_strings, only: operator(.ich.)
integer, intent(inout) :: i
i=geti()
end subroutine readi
 
end module gems_input_parsing

