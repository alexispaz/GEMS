program main
use intrinsic_class
use lion_class
integer,parameter         :: dp=8

real(dp),target           :: a=3.14_dp

type(real_dlist), target   :: realdl
class(real_dlist), pointer  :: node

type(lion_list), target   :: lions
class(lion_list), pointer  :: thislion, prev

type(lion_cdlist), target   :: clions
class(lion_cdlist), pointer  :: thisclion
 

type(lion), target        :: lion1, lion2, lion3
type(lion), pointer       :: lionptr
type(lionking)            :: simba
         

print *, 'Real double list:'
print *, 'First 3 nodes have real numbers'
print *, '4th node have a pointer to a real'

! Real list
call realdl%add_hardcpy( 1.0d0 )
call realdl%add_hardcpy( 2.6d0 )
call realdl%add_hardcpy( 3.5d0 )
call realdl%add_soft( a )

! Simple loop  
node => realdl
do while( associated(node%next) )
    node => node%next
    print *, node%o
enddo
print *, ''

print *, 'Deleting second element'
node=>realdl%next
node=>node%next
call node%deattach()
call node%destroy_node()
deallocate(node)

node => realdl
do while( associated(node%next) )
    node => node%next
    print *, node%o
enddo
print *, ''
print *, ''
       


print *, 'Lion list:'
print *, 'First 2 nodes are pointers to the lion objects'
print *, '3rd node has allocated a lion object'
print *, '4th node is a pointer to a lionking object'

! Add lion 1 after the head
lion1%color = "orange"
call lions%add_soft( lion1 )

! Add lion 2 after the head
lion2%color = "red"
call lions%add_soft( lion2 )

! Add lion 3 after the head
lion3%color = "brown"
call lions%add_hardcpy( lion3 )

! Add lionking after the head
simba%color = "gold"
simba%name  = "simba"
call lions%add_soft( simba )

print *, 'Access common properties in the class'
thislion => lions
do while( associated(thislion%next) )
    thislion => thislion%next

    print *, 'a lion ', thislion%o%color

enddo
print *, ''
              
print *, 'Access all properties using select type.'
thislion => lions
do while( associated(thislion%next) )
    thislion => thislion%next

    ! Print the lion name if is a lionking
    select type (a=>thislion%o)
    type is (lionking)
      print *, a%name, a%color
    class default  
      ! Print the lion color
      print *, 'simple lion ', thislion%o%color
    end select

enddo
print *, ''


print *, 'Access all properties using internal procedures (cleaner in my opinion)'
thislion => lions
do while( associated(thislion%next) )
    thislion => thislion%next
    call thislion%o%printme
enddo
print *, ''

print *, 'Removing pointer to the red lion'
thislion => lions
do while( associated(thislion%next) )
    prev => thislion
    thislion => thislion%next

    if (associated(thislion%o,target=lion2)) then
      call thislion%deattach(prev)
      call thislion%destroy_node()
      deallocate(thislion)
      thislion=>prev
      cycle
    endif

    call thislion%o%printme
enddo
! print *, ''




print *, 'Lion circular double list:'
print *, 'First 2 nodes are pointers to the lion objects'
print *, '3rd node has allocated a lion object'
print *, '4th node is a pointer to a lionking object'

! Initialize the list
call clions%init()

! Add lions
call clions%add_soft( lion1 )
call clions%add_soft( lion2 )
call clions%add_hardcpy( lion3 )
call clions%add_soft( simba )

print *, 'Access all properties using internal procedures (cleaner in my opinion)'
thisclion => clions
do while( associated(thisclion%next%o) )
    thisclion => thisclion%next
    call thisclion%o%printme
enddo
print *, ''

print *, 'Deallocating the brown lion' 
thisclion => clions
do while( associated(thisclion%next%o) )
    thisclion => thisclion%next
    if (thisclion%o%color=='brown') then
      call thisclion%deattach()
      call thisclion%destroy_node()
      exit
    endif
enddo

thisclion => clions
do while( associated(thisclion%next%o) )
    thisclion => thisclion%next
    call thisclion%o%printme
enddo
print *, '' 
          

end program main

