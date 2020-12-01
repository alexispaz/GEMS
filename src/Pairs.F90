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

 
module gems_pairs
use gems_program_types
use gems_algebra 
use gems_constants,     only: dp,ev_ui,ui_ev,dm
use gems_inq_properties
use gems_neighbour

implicit none

private
public  :: pair_new
       
type,extends(intergroup) :: lj
  real(dp)  :: e,s
  contains
  procedure :: interact => lj_interact
  procedure,nopass :: cli => pair_cli
end type

type,extends(intergroup) :: plj
  real(dp)  :: e,s,prc,frc
  contains
  procedure :: interact => plj_interact
  procedure,nopass :: cli => pair_cli
end type
 
type,extends(intergroup) :: slj
  real(dp)  :: e,s,sa4,sa8,rs
  contains
  procedure :: interact => slj_interact   
  procedure,nopass :: cli => pair_cli
end type
 
type,extends(intergroup) :: wca
  real(dp)  :: e,s
  contains
  procedure :: interact => wca_interact
  procedure,nopass :: cli => pair_cli
end type
 
type,extends(intergroup) :: rwca
  real(dp)  :: e,s,r
  contains
  procedure :: interact => rwca_interact
  procedure,nopass :: cli => pair_cli
end type
 
type,extends(intergroup) :: cuw
  real(dp)    :: ks,rs,radc,radi,rc
  integer     :: eje
  contains
  procedure :: interact => cuw_interact
  procedure,nopass :: cli => pair_cli
end type
  
type,extends(intergroup) :: sho
  real(dp)  :: k,r
  contains
  procedure :: interact => sho_interact
  procedure,nopass :: cli => pair_cli
end type
 
type,extends(intergroup) :: shocm
  real(dp)  :: k,r
  contains
  procedure :: interact => shocm_interact
  procedure,nopass :: cli => pair_cli
end type
        
type,extends(intergroup) :: shofix
  real(dp)              :: k
  real(dp),allocatable  :: r(:,:)
  contains
  procedure :: interact => shofix_interact
  procedure,nopass :: cli => pair_cli
end type
                                      
type,extends(intergroup) :: sm1
  real(dp)              :: r1,e
  contains
  procedure :: interact => sm1_interact   
  procedure,nopass :: cli => pair_cli
end type
                            
contains


subroutine pair_new(ig,g1,w,g2)
character(*),intent(in)   :: w
type(group)               :: g1
type(group),intent(in),optional :: g2
class(intergroup),pointer :: ig

select case(w)
case('lj')     ; allocate(lj::ig)
case('plj')    ; allocate(plj::ig)
case('slj')    ; allocate(slj::ig)
case('wca')    ; allocate(wca::ig)
case('rwca')   ; allocate(rwca::ig)
case('cuw')    ; allocate(cuw::ig)
case('sho')    ; allocate(sho::ig)
case('shocm')  ; allocate(shocm::ig)
case('shofix') ; allocate(shofix::ig)
  call werr("Only self interaction, use `under` keyword.",present(g2)) 
case('sm1')    ; allocate(sm1::ig)
case default
  call werr('Pair potential not found')
endselect
   
if(present(g2)) then
  call ig%init(g1=g1,g2=g2)
else
  call ig%init(g1=g1)
endif

end subroutine pair_new

subroutine pair_cli(ig)
use gems_input_parsing
class(intergroup),intent(inout) :: ig
integer                         :: i1
real(dp)                        :: f1,f2,f3,f4,f5
           
select type(ig)
type is(lj) 
  call readf(f1)
  call readf(f2)
  call readf(f3)
  call lj_set(ig,f1,f2,f3)
type is(plj)
  call readf(f1)
  call readf(f2)
  call plj_set(ig,f1,f2,10._dp)
type is(slj)
  call readf(f1)
  call readf(f2)
  call readf(f3)
  call slj_set(ig,f1,f2,f3)
type is(wca)
  call readf(f1)
  call readf(f2)
  call wca_set(ig,f1,f2)
type is(rwca)
  call readf(f1) !epsilon
  call readf(f2) !sigma
  call readf(f3) !start
  call rwca_set(ig,f1,f2,f3)
type is(cuw)
  call readf(f1) !epsilon
  call readf(f2) !sigma
  call readf(f3) !start
  call readf(f4) !radious
  call readi(i1) !coordinate
  call readf(f5) !rcut
  call cuw_set(ig,f1,f2,f3,f4,i1,f5)
type is(sho)
  call readf(f1)
  call readf(f2)
  call sho_set(ig,f1,f2)
type is(shocm)
  call readf(f1)
  call readf(f2)
  call shocm_set(ig,f1,f2)
type is(shofix)
  call readf(f1,ev_ui)
  call shofix_set(ig,f1)
  ig%lista => intergroup0_empty
type is(sm1)
  call readf(f1)
  call readf(f2)
  call readf(f3,ev_ui)
  call sm1_set(ig,f1,f2,f3)
class default
  call werr('Pair potential not found')
endselect

                           
end subroutine
          
subroutine lj_set(ig,e,s,rc)
! Lennard Jones with parabolic smooth function
real(dp),intent(in) :: e,s,rc
type(lj)            :: ig
ig%e=e
ig%s=s
call ig%setrc(rc)
end subroutine lj_set
  
subroutine lj_interact(ig)
class(lj),intent(inout) :: ig
integer                :: i,j,n,l
real(dp)               :: vd(dm),dr,factor2(dm)
real(dp)               :: p,f,s_dr
type(atom),pointer     :: o1,o2
 

ig%epot = 0._dp

do i = 1,ig%n(1)
  o1=>ig%at(i)%o

  do l = 1, ig%nn(i)  ! sobre los vecinos

    j  = ig%list(i,l)
    o2 =>ig%at(j)%o

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>ig%rcut2) cycle

    dr=sqrt(dr)

    s_dr=(ig%s/dr)**6

    p=4._dp*ig%e*(s_dr*s_dr-s_dr)
    f=ig%e*(48._dp*s_dr-24._dp)*s_dr/dr   ! -dp/dr

    p=p*ev_ui
    f=f*ev_ui

    o1%epot = o1%epot + p*0.5_dp
    ig%epot = ig%epot + p*0.5_dp

    factor2(1:dm) = f*vd(1:dm)/dr      ! (-dp/dr)*(-grad_1dr)
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 

    if(j<=nlocal) then

      if(.not.ig%newton) cycle
  
      o2%epot = o2%epot + p*0.5_dp
      ig%epot = ig%epot + p*0.5_dp

      o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
      if(b_gvirial) then
        do n = 1,3
          virial(n,:) = virial(n,:) + o2%pos(n)*factor2(:) 
        enddo
      endif

    else if (o2%tag>o1%tag) then
      ! Una sola vez por par para que se pueda calcular el virial

      o2%epot = o2%epot + p*0.5_dp
      ig%epot = ig%epot + p*0.5_dp
       
      o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
      if(b_gvirial) then
        do n = 1,3
          virial(n,:) = virial(n,:) + box(n)*floor(o2%pos(n)*one_box(n))*factor2(:) 
        enddo
      endif
    
    endif
    
  enddo 

enddo     

end subroutine lj_interact 
                             

subroutine plj_set(ig,e,s,ljcut)
! Lennard Jones with parabolic smooth function
real(dp),intent(in)      :: e,s,ljcut
type(plj)                :: ig
ig%e=e
ig%s=s
call ig%setrc(ljcut)
ig%prc=4._dp*e*((s/ljcut)**12-(s/ljcut)**6)
ig%frc=24._dp*e*(2._dp*(s**12/ljcut**13)-(s**6/ljcut**7))
end subroutine plj_set
  
subroutine plj_interact(ig)
class(plj),intent(inout) :: ig
integer                  :: i,j,n,l
real(dp)                 :: vd(dm),dr,factor2(dm)
real(dp)                 :: p,f,s_dr
type(atom),pointer       :: o1,o2
 
ig%epot = 0._dp

do i = 1,ig%n(1)
  o1=>ig%at(i)%o

  do l = 1, ig%nn(i)  ! sobre los vecinos

    j  = ig%list(i,l)
    o2 =>ig%at(j)%o

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>ig%rcut2) cycle

    dr=sqrt(dr)

    s_dr=(ig%s/dr)**6

    p=4._dp*ig%e*(s_dr*s_dr-s_dr)-ig%prc-ig%frc*(dr-ig%rcut)
    f=ig%e*(48._dp*s_dr-24._dp)*s_dr/dr+ig%frc        ! (-dp/dr)

    p=p*ev_ui
    f=f*ev_ui

    o1%epot = o1%epot + p*0.5_dp
    ig%epot = ig%epot + p*0.5_dp

    factor2(1:dm) = f*vd(1:dm)/dr               ! (-dp/dr)*(-grad_1dr)
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 
                 
    if(b_gvirial) then
      do n = 1,3
        virial(n,:) = virial(n,:) - o1%pos(n)*factor2(:) 
      enddo
    endif
       
    if(j<=nlocal) then

      if(.not.ig%newton) cycle
  
      o2%epot = o2%epot + p*0.5_dp
      ig%epot = ig%epot + p*0.5_dp

      o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
      if(b_gvirial) then
        do n = 1,3
          virial(n,:) = virial(n,:) + o2%pos(n)*factor2(:) 
        enddo
      endif

    else if (o2%tag>o1%tag) then
      ! Una sola vez por par para que se pueda calcular el virial

      o2%epot = o2%epot + p*0.5_dp
      ig%epot = ig%epot + p*0.5_dp
       
      o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
      if(b_gvirial) then
        do n = 1,3
          virial(n,:) = virial(n,:) + box(n)*floor(o2%pos(n)*one_box(n))*factor2(:) 
        enddo
      endif
    
    endif
    
  enddo 

enddo     

end subroutine plj_interact 

                               

subroutine slj_set(ig,e,s,rs)
! Lennard Jones with splined smooth function
! The equations for this potential can be found in
! Grønbech-Jensen, N., & Farago, O. (2014). Constant pressure and temperature
! discrete-time Langevin molecular dynamics. The Journal of Chemical Physics,
! 141(19), 194108. http://doi.org/10.1063/1.4901303
real(dp),intent(in)  :: e,s,rs
real(dp)             :: s_dr,urs,dudrs,rc
type(slj)            :: ig

ig%e=e
ig%s=s
ig%rs=rs

! cutoff
s_dr=(s/rs)**6
urs=4._dp*e*(s_dr*s_dr-s_dr)
dudrs=-e*(48._dp*s_dr-24._dp)*s_dr/rs

rc=rs-32._dp*urs/(11._dp*dudrs)

call ig%setrc(rc)
          
! other parameters
ig%sa4=(8._dp*urs+(ig%rcut-rs)*dudrs)/(4*(ig%rcut-rs)**4)
ig%sa8=-(4._dp*urs+(ig%rcut-rs)*dudrs)/(4*(ig%rcut-rs)**8)

end subroutine slj_set
  
subroutine slj_interact(ig)
class(slj),intent(inout) :: ig
integer                  :: i,j,l,n
real(dp)                 :: vd(dm),dr,factor2(dm)
real(dp)                 :: p,f,s_dr
type(atom),pointer       :: o1,o2
 
ig%epot = 0._dp

do i = 1,ig%n(1)
  o1=>ig%at(i)%o

  do l = 1, ig%nn(i)  ! sobre los vecinos

    j  = ig%list(i,l)
    o2 =>ig%at(j)%o

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>ig%rcut2) cycle

    dr=sqrt(dr)

    if(dr<ig%rs) then

      s_dr=(ig%s/dr)**6
      p=4._dp*ig%e*(s_dr-1._dp)*s_dr
      f=ig%e*(48._dp*s_dr-24._dp)*s_dr/dr      ! -dp/dr

    else
     
      s_dr=(dr-ig%rcut)**4
      p=ig%sa4*s_dr+ig%sa8*s_dr*s_dr
      f=-(4._dp*ig%sa4*s_dr+8._dp*ig%sa8*s_dr*s_dr)/(dr-ig%rcut)        ! -dp/dr
      
    endif

    p=p*ev_ui
    f=f*ev_ui

    o1%epot = o1%epot + p*0.5_dp
    ig%epot = ig%epot + p*0.5_dp

    factor2(1:dm) = f*vd(1:dm)/dr           ! (-dp/dr)*(-grad_1dr)
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 

    if(b_gvirial) then
      do n = 1,3
        virial(n,:) = virial(n,:) - o1%pos(n)*factor2(:) 
      enddo
    endif
    
    if(j<=nlocal) then

      if(.not.ig%newton) cycle

      o2%epot = o2%epot + p*0.5_dp
      ig%epot = ig%epot + p*0.5_dp

      o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
      if(b_gvirial) then
        do n = 1,3
          virial(n,:) = virial(n,:) + o2%pos(n)*factor2(:) 
        enddo
      endif

    else if (o2%tag>o1%tag) then
      ! Una sola vez por par para que se pueda calcular el virial

      o2%epot = o2%epot + p*0.5_dp
      ig%epot = ig%epot + p*0.5_dp
       
      o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
      if(b_gvirial) then
        do n = 1,3
          virial(n,:) = virial(n,:) + box(n)*floor(o2%pos(n)*one_box(n))*factor2(:) 
        enddo
      endif
    
    endif
         
  enddo     

enddo 
  

              
end subroutine  


               
subroutine wca_set(ig,e,s)
! The Weeks-Chandler-Andersen potential: a Lennard-Jones potential shifted and
! cut at its minimum to retain a purely repulsive potential.
real(dp),intent(in)  :: e,s
type(wca)            :: ig
                                 
! parametros
ig%e=e
ig%s=s
call ig%setrc(2._dp**(1./6.)*s)
             
end subroutine wca_set
  
subroutine wca_interact(ig)
class(wca),intent(inout) :: ig
integer                  :: i,j,n,l
real(dp)                 :: vd(dm),dr,factor2(dm)
real(dp)                 :: p,f,s_dr
type(atom),pointer       :: o1,o2
 
! Por info ver LJ0.f90
ig%epot = 0._dp

do i = 1,ig%n(1)
  o1=>ig%at(i)%o

  do l = 1, ig%nn(i)  ! sobre los vecinos

    j  = ig%list(i,l)
    o2 =>ig%at(j)%o

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 


    if(dr>ig%rcut2) cycle

    dr=sqrt(dr)

    s_dr=(ig%s/dr)**6

    p=4._dp*ig%e*(s_dr*s_dr-s_dr)+ig%e
    f=ig%e*(48._dp*s_dr-24._dp)*s_dr/dr         ! -dp/dr

    p=p*ev_ui
    f=f*ev_ui
         
    o1%epot = o1%epot + p*0.5_dp
    ig%epot = ig%epot + p*0.5_dp
      
    factor2(1:dm) = f*vd(1:dm)/dr            ! (-dp/dr)*(-grad_1dr)
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 

    if(b_gvirial) then
      do n = 1,3
        virial(n,:) = virial(n,:) - o1%pos(n)*factor2(:) 
      enddo
    endif
    
    if(j<=nlocal) then

      if(.not.ig%newton) cycle

      o2%epot = o2%epot + p*0.5_dp
      ig%epot = ig%epot + p*0.5_dp

      o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
      if(b_gvirial) then
        do n = 1,3
          virial(n,:) = virial(n,:) + o2%pos(n)*factor2(:) 
        enddo
      endif

    else if (o2%tag>o1%tag) then
      ! Una sola vez por par para que se pueda calcular el virial

      o2%epot = o2%epot + p*0.5_dp
      ig%epot = ig%epot + p*0.5_dp
       
      o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
      if(b_gvirial) then
        do n = 1,3
          virial(n,:) = virial(n,:) + box(n)*floor(o2%pos(n)*one_box(n))*factor2(:) 
        enddo
      endif
    
    endif
 
  enddo 

enddo     

end subroutine  
        


subroutine rwca_set(ig,e,s,r)
! Reversed Weeks-Chandler-Andersen potential
real(dp),intent(in)  :: e,s,r
type(rwca)           :: ig
ig%e=e
ig%s=s 
ig%r=r*r

ig%rcut=2._dp**(1./6.)*s+r
ig%rcut2=ig%rcut*ig%rcut
              
end subroutine rwca_set
  
subroutine rwca_interact(ig)
! Reversed Weeks-Chandler-Andersen potential
class(rwca),intent(inout) :: ig
integer                   :: i,j,l
real(dp)                  :: vd(dm),dr,factor2(dm)
real(dp)                  :: p,f,s_dr
type(atom),pointer        :: o1,o2
 
! Por info ver LJ0.f90
ig%epot = 0._dp

! Sobre los atomos
do i = 1,ig%n(1)
  o1 => ig%at(i)%o

  ! Sobre los vecinos
  do l = 1, ig%nn(i) 

    j  = ig%list(i,l)
    o2 => ig%at(j)%o

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>ig%rcut2) cycle
    if(dr<ig%r) cycle

    dr=sqrt(dr)
    dr=ig%rcut-dr
    
    s_dr=(ig%s/dr)**6

    p=4._dp*ig%e*(s_dr*s_dr-s_dr)+ig%e
    p=p*ev_ui

    o1%epot = o1%epot + p*0.5_dp
    o2%epot = o2%epot + p*0.5_dp
    ig%epot = ig%epot + p

    f=ig%e*(48._dp*s_dr-24._dp)*s_dr/dr   ! dp/dr (TODO: chequear signos)
    f=-f*ev_ui

    factor2(1:dm) = f*vd(1:dm)/dr              ! (-dp/dr)*(-grad_1dr)  (TODO: chequear signos)
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 
    o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
    
  enddo 

enddo     

end subroutine  
                                
 
subroutine cuw_set(ig,ks,rs,radc,radi,eje,rc)
! Reversed cylindrically confined Weeks-Chandler-Andersen potential
real(dp),intent(in) :: ks,rs,radc,radi,rc
integer,intent(in)  :: eje
type(cuw)           :: ig
ig%ks=ks   
ig%rs=rs*rs
ig%radc=radc 
ig%radi=radi 
ig%eje=eje  
call ig%setrc(rc)
end subroutine cuw_set
  
subroutine cuw_interact(ig)
! Reversed cylindrically confined Weeks-Chandler-Andersen potential
use gems_algebra, only:fcut2_dfcut
class(cuw),intent(inout) :: ig
integer                  :: i,j,l
real(dp)                 :: vd(dm),dr,factor2(dm)
real(dp)                 :: p,proj,w,r1,r2,norm
type(atom),pointer       :: oi,oj
real(dp)                 :: aux,df1,df2
 
ig%epot = 0._dp

do i = 1,ig%n(1)
  oi => ig%at(i)%o

  ! Compute the normalization factor
  norm=0._dp
  p=0._dp
  do l = 1, ig%nn(i)  ! sobre los vecinos

    j  = ig%list(i,l)
    oj => ig%at(j)%o

    vd = vdistance(oj,oi, mic)

    ! Skip this if the z distance is below the boundary
    if(vd(3)<ig%rs) cycle

    dr = dot_product(vd,vd) 

    ! Skip using neighboor lsit
    if(dr>ig%rcut2) cycle

    dr=sqrt(dr)
    proj=dr/vd(ig%eje)

    ! Skip if the projection is outside the circle
    if(proj>ig%radc+ig%radi) cycle

    ! Compute the weight
    r1=proj-ig%radi
    r2=proj+ig%radi
    w=fcut2_dfcut(ig%radc,r1,r2,df1)!-fcut2(-radc,r1,r2,df2)
    
    ! Accumulate and save the unnormalized weight
    norm=norm+w
    
  enddo 
    
  do l = 1, ig%nn(i)  ! sobre los vecinos

    j  = ig%list(i,l)
    oj => ig%at(j)%o
        
    vd = vdistance(oj,oi, mic)

    ! Skip this if the z distance is below the boundary
    if(vd(3)<ig%rs) cycle

    dr = dot_product(vd,vd) 

    ! Skip using neighboor lsit
    if(dr>ig%rcut2) cycle

    dr=sqrt(dr)
    proj=dr/vd(ig%eje)

    ! Skip if the projection is outside the circle
    if(proj>ig%radc+ig%radi) cycle
              
    ! Compute the weight
    r1=proj-ig%radi
    r2=proj+ig%radi

    !TODO: If radi<radc then fcut2(-radc,r1,r2,df2)=0 always
    w=fcut2_dfcut(ig%radc,r1,r2,df1)-fcut2_dfcut(-ig%radc,r1,r2,df2)
                 
    aux=(vd(3)-ig%rs)
    p=w*ig%ks*aux*aux*0.5_dp/norm

    oi%epot = oi%epot + p*0.5_dp
    oj%epot = oj%epot + p*0.5_dp
    ig%epot = ig%epot + p

    factor2(1:2) = (df1-df2)*vd(1:2)/dr
    factor2(1:2) = p*(1._dp/factor2(1:2)+factor2(1:2)/w)
    factor2(3) = -w*ig%ks*aux*vd(3)/dr

    oi%force(1:3) = oi%force(1:3) + factor2(1:3) 
    oj%force(1:3) = oj%force(1:3) - factor2(1:3)
    
  enddo 

enddo     

end subroutine  
                          
      
subroutine sho_set(ig,k,r)
type(sho)   :: ig
real(dp)    :: k,r
ig%k=k
ig%r=r
call ig%setrc(100._dp)
end subroutine        
      
subroutine sho_interact(ig)
class(sho),intent(inout)    :: ig
integer                     :: i,j,l
real(dp)                    :: vd(dm),dr,factor2(dm)
real(dp)                    :: p,f
type(atom),pointer          :: o1,o2
 
ig%epot = 0._dp

do i = 1,ig%n(1)
  o1 => ig%at(i)%o

  do l = 1, ig%nn(i)  ! sobre los vecinos
    j  = ig%list(i,l)
    o2 => ig%at(j)%o

    vd = vdistance( o1, o2 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>ig%rcut2) cycle

    dr=sqrt(dr)-ig%r

    p=-0.5_dp*ig%k*dr*dr

    o1%epot = o1%epot + p*0.5_dp
    o2%epot = o2%epot + p*0.5_dp
    ig%epot = ig%epot + p

    f=-ig%k*dr

    factor2(1:dm) = -f*vd(1:dm)   ! el dividido dr esta en f
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 
    o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
    
  enddo 

enddo     

end subroutine  
       
subroutine shocm_set(ig,k,r)
type(shocm) :: ig
real(dp)    :: k,r
ig%k=k
ig%r=r
call ig%setrc(100._dp)
end subroutine        
 
subroutine shocm_interact(ig)
class(shocm),intent(inout) :: ig
real(dp)                   :: vd(dm),dr,paux,faux,m,m2
type(atom_dclist),pointer  :: la,alist
integer                    :: i
 
ig%epot = 0._dp

vd = atom_dclist_inq_cmpos(ig%a,m)  &
    -atom_dclist_inq_cmpos(ig%b,m2)  ! Ojo el orden
dr = sqrt(dot_product(vd,vd))

paux = dr-ig%r
faux = ig%k*paux
paux = paux*faux*0.5_dp
vd=faux*vd/dr

! Sobre los dos grupos
alist => ig%a
do i = 1,2

  ! sobre los atomos
  la => alist%next
  do while(.not.associated(la,target=alist))
    la%o%epot = la%o%epot + paux*0.5_dp
    la%o%force = la%o%force - vd*la%o%mass/m
    la=>la%next
  enddo     

  vd=-vd
  m=m2
  alist => ig%b

enddo

ig%epot = ig%epot + paux*ig%n(1) + paux*(ig%n(3)-ig%n(2))

end subroutine  



subroutine shofix_set(ig,k)
real(dp),intent(in)       :: k
integer                   :: i
type(atom_dclist),pointer :: la
type(shofix)              :: ig
   
ig%k=k

allocate(ig%r(ig%n(1),3))

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next
  ig%r(i,:)=la%o%pos(:)
enddo     
             
end subroutine shofix_set
  
subroutine shofix_interact(ig)
use gems_program_types, only: atom_distancetopoint
class(shofix),intent(inout)   :: ig
integer                       :: i
real(dp)                      :: vd(dm),dr,factor2(dm)
real(dp)                      :: p,f
type(atom_dclist),pointer :: la
 
ig%epot = 0._dp

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next

  vd = atom_distancetopoint(la%o,ig%r(i,:))  ! vector a-->p
  dr = dot_product(vd,vd) 

  p=-0.5_dp*ig%k*dr*dr
     
  la%o%epot = la%o%epot + p*0.5_dp
  ig%epot = ig%epot + p

  f=-ig%k*dr

  factor2(1:dm)=f*vd(1:dm)   ! el dividido dr esta en f

  la%o%force(1:dm) = la%o%force(1:dm) - factor2(1:dm) 
  
enddo     

end subroutine  

subroutine sm1_set(ig,r1,r2,e)
! A smooth function. It has to reach cero (so the rcut hav sense)
type(sm1)           :: ig
real(dp),intent(in) :: r1,r2,e
ig%r1=r1
ig%e=e
call ig%setrc(r2)
end subroutine sm1_set
  
subroutine sm1_interact(ig)
use gems_constants, only:pi,pio2
class(sm1),intent(inout) :: ig
integer                  :: i,j,l,n
real(dp)                 :: vd(dm),dr,factor2(dm)
real(dp)                 :: p,f
real(dp)                 :: r1,r2,aux
type(atom),pointer       :: o1,o2
 
r1=ig%r1
r2=ig%rcut

ig%epot = 0._dp

do i = 1,ig%n(1)
  o1=>ig%at(i)%o

  do l = 1, ig%nn(i)  ! sobre los vecinos

    j  = ig%list(i,l)
    o2 =>ig%at(j)%o

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>ig%rcut2) cycle

    dr=sqrt(dr)

    if(dr>r1) then
      aux=pi*(dr-(r1+r2)*0.5_dp)/(r2-r1)
      p=0.5_dp+0.5_dp*sin(aux)
      f=cos(aux)*pio2*dr/(r2-r1)      
      p=p*ig%e-ig%e
      f=f*ig%e
    else
      p=-ig%e
      f=0
    endif

    o1%epot = o1%epot + p*0.5_dp
    ig%epot = ig%epot + p*0.5_dp

    factor2(1:dm) = -f*vd(1:dm)/dr   
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 

    if(b_gvirial) then
      do n = 1,3
        virial(n,:) = virial(n,:) - o1%pos(n)*factor2(:) 
      enddo
    endif
          
    if(j<=nlocal) then

      if(.not.ig%newton) cycle

      o2%epot = o2%epot + p*0.5_dp
      ig%epot = ig%epot + p*0.5_dp

      o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
      if(b_gvirial) then
        do n = 1,3
          virial(n,:) = virial(n,:) + o2%pos(n)*factor2(:) 
        enddo
      endif

    else if (o2%tag>o1%tag) then
      ! Una sola vez por par para que se pueda calcular el virial

      o2%epot = o2%epot + p*0.5_dp
      ig%epot = ig%epot + p*0.5_dp
       
      o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
      if(b_gvirial) then
        do n = 1,3
          virial(n,:) = virial(n,:) + box(n)*floor(o2%pos(n)*one_box(n))*factor2(:) 
        enddo
      endif
    
    endif
         
  enddo     

enddo 
  
              
end subroutine  

    
end module gems_pairs
