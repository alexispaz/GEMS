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
use gems_neighbor

implicit none

private
public  :: pair_new
       
type,extends(ngroup) :: lj
  real(dp)  :: e,s
  contains
  procedure :: interact => lj_interact
  procedure,nopass :: cli => pair_cli
end type

type,extends(ngroup) :: plj
  real(dp)  :: e,s,prc,frc
  contains
  procedure :: interact => plj_interact
  procedure,nopass :: cli => pair_cli
end type
 
type,extends(ngroup) :: slj
  real(dp)  :: e,s,sa4,sa8,rs
  contains
  procedure :: interact => slj_interact   
  procedure,nopass :: cli => pair_cli
end type
 
type,extends(ngroup) :: wca
  real(dp)  :: e,s
  contains
  procedure :: interact => wca_interact
  procedure,nopass :: cli => pair_cli
end type
 
type,extends(ngroup) :: rwca
  real(dp)  :: e,s,r
  contains
  procedure :: interact => rwca_interact
  procedure,nopass :: cli => pair_cli
end type
 
type,extends(ngroup) :: cuw
  real(dp)    :: ks,rs,radc,radi,rc
  integer     :: eje
  contains
  procedure :: interact => cuw_interact
  procedure,nopass :: cli => pair_cli
end type
  
type,extends(ngroup) :: sho
  real(dp)  :: k,r
  contains
  procedure :: interact => sho_interact
  procedure,nopass :: cli => pair_cli
end type
 
type,extends(ngroup) :: shocm
  real(dp)  :: k,r
  contains
  procedure :: interact => shocm_interact
  procedure,nopass :: cli => pair_cli
end type
                     
type,extends(ngroup) :: sm1
  real(dp)              :: r1,e
  contains
  procedure :: interact => sm1_interact   
  procedure,nopass :: cli => pair_cli
end type
                            
contains


subroutine pair_new(g,w)
character(*),intent(in)  :: w
class(ngroup),pointer    :: g

select case(w)
case('lj')     ; allocate(lj::g)
case('plj')    ; allocate(plj::g)
case('slj')    ; allocate(slj::g)
case('wca')    ; allocate(wca::g)
case('rwca')   ; allocate(rwca::g)
case('cuw')    ; allocate(cuw::g)
case('sho')    ; allocate(sho::g)
case('shocm')  ; allocate(shocm::g)
case('sm1')    ; allocate(sm1::g)
case default
  call werr('Pair potential not found')
endselect

end subroutine pair_new

subroutine pair_cli(g)
use gems_input_parsing
class(ngroup),intent(inout) :: g
integer                         :: i1
real(dp)                        :: f1,f2,f3,f4,f5
           
select type(g)
type is(lj) 
  call readf(f1)
  call readf(f2)
  call readf(f3)
  call lj_set(g,f1,f2,f3)
type is(plj)
  call readf(f1)
  call readf(f2)
  call plj_set(g,f1,f2,10._dp)
type is(slj)
  call readf(f1)
  call readf(f2)
  call readf(f3)
  call slj_set(g,f1,f2,f3)
type is(wca)
  call readf(f1)
  call readf(f2)
  call wca_set(g,f1,f2)
type is(rwca)
  call readf(f1) !epsilon
  call readf(f2) !sigma
  call readf(f3) !start
  call rwca_set(g,f1,f2,f3)
type is(cuw)
  call readf(f1) !epsilon
  call readf(f2) !sigma
  call readf(f3) !start
  call readf(f4) !radious
  call readi(i1) !coordinate
  call readf(f5) !rcut
  call cuw_set(g,f1,f2,f3,f4,i1,f5)
type is(sho)
  call readf(f1)
  call readf(f2)
  call sho_set(g,f1,f2)
type is(shocm)
  call readf(f1)
  call readf(f2)
  call shocm_set(g,f1,f2)
type is(sm1)
  call readf(f1)
  call readf(f2)
  call readf(f3,ev_ui)
  call sm1_set(g,f1,f2,f3)
class default
  call werr('Pair potential not found')
endselect

                           
end subroutine
          
subroutine lj_set(g,e,s,rc)
! Lennard Jones with parabolic smooth function
real(dp),intent(in) :: e,s,rc
type(lj)            :: g
g%e=e
g%s=s
call g%setrc(rc)
end subroutine lj_set
  
subroutine lj_interact(g)
class(lj),intent(inout)   :: g
integer                   :: i,ii,j,n,l
real(dp)                  :: vd(dm),dr,factor2(dm)
real(dp)                  :: p,f,s_dr
type(atom),pointer        :: o1,o2
type(atom_dclist),pointer :: la

g%epot = 0._dp

la => g%ref%alist
do ii = 1,g%ref%nat
  la => la%next
  o1 => la%o
  i = o1%gid(g%id)

  do l = 1, g%nn(i)  ! sobre los vecinos

    j  = g%list(i,l)
    o2 =>g%a(j)%o

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>g%rcut2) cycle

    dr=sqrt(dr)

    s_dr=(g%s/dr)**6

    p=4._dp*g%e*(s_dr*s_dr-s_dr)
    f=g%e*(48._dp*s_dr-24._dp)*s_dr/dr   ! -dp/dr

    p=p*ev_ui
    f=f*ev_ui

    o1%epot = o1%epot + p*0.5_dp
    g%epot = g%epot + p*0.5_dp

    factor2(1:dm) = f*vd(1:dm)/dr      ! (-dp/dr)*(-grad_1dr)
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 

    if (associated(o2%ghost)) then
      if (o2%gid(g%id)>o1%gid(g%id)) then

        o2%epot = o2%epot + p*0.5_dp
        g%epot = g%epot + p*0.5_dp
         
        o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
        if(b_gvirial) then
          do n = 1,3
            virial(n,:) = virial(n,:) + box(n)*floor(o2%pos(n)*one_box(n))*factor2(:) 
          enddo
        endif
      
      endif
    endif
    
  enddo 

enddo     

end subroutine lj_interact 
                             

subroutine plj_set(g,e,s,ljcut)
! Lennard Jones with parabolic smooth function
real(dp),intent(in)      :: e,s,ljcut
type(plj)                :: g
g%e=e
g%s=s
call g%setrc(ljcut)
g%prc=4._dp*e*((s/ljcut)**12-(s/ljcut)**6)
g%frc=24._dp*e*(2._dp*(s**12/ljcut**13)-(s**6/ljcut**7))
end subroutine plj_set
  
subroutine plj_interact(g)
class(plj),intent(inout) :: g
integer                  :: i,ii,j,n,l
real(dp)                 :: vd(dm),dr,factor2(dm)
real(dp)                 :: p,f,s_dr
type(atom),pointer       :: o1,o2
type(atom_dclist),pointer   :: la

g%epot = 0._dp

la => g%ref%alist
do ii = 1,g%ref%nat
  la => la%next
  o1 => la%o
  i = o1%gid(g%id)

  do l = 1, g%nn(i)  ! sobre los vecinos

    j  = g%list(i,l)
    o2 =>g%a(j)%o

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>g%rcut2) cycle

    dr=sqrt(dr)

    s_dr=(g%s/dr)**6

    p=4._dp*g%e*(s_dr*s_dr-s_dr)-g%prc-g%frc*(dr-g%rcut)
    f=g%e*(48._dp*s_dr-24._dp)*s_dr/dr+g%frc        ! (-dp/dr)

    p=p*ev_ui
    f=f*ev_ui

    o1%epot = o1%epot + p*0.5_dp
    g%epot = g%epot + p*0.5_dp

    factor2(1:dm) = f*vd(1:dm)/dr               ! (-dp/dr)*(-grad_1dr)
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 
                 
    if(b_gvirial) then
      do n = 1,3
        virial(n,:) = virial(n,:) - o1%pos(n)*factor2(:) 
      enddo
    endif
       
    if (associated(o2%ghost)) then
      if (o2%gid(g%id)>o1%gid(g%id)) then
        ! Una sola vez por par para que se pueda calcular el virial

        o2%epot = o2%epot + p*0.5_dp
        g%epot = g%epot + p*0.5_dp
         
        o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
        if(b_gvirial) then
          do n = 1,3
            virial(n,:) = virial(n,:) + box(n)*floor(o2%pos(n)*one_box(n))*factor2(:) 
          enddo
        endif
      
      endif
    endif
    
  enddo 

enddo     

end subroutine plj_interact 

                               

subroutine slj_set(g,e,s,rs)
! Lennard Jones with splined smooth function
! The equations for this potential can be found in
! Grønbech-Jensen, N., & Farago, O. (2014). Constant pressure and temperature
! discrete-time Langevin molecular dynamics. The Journal of Chemical Physics,
! 141(19), 194108. http://doi.org/10.1063/1.4901303
real(dp),intent(in)  :: e,s,rs
real(dp)             :: s_dr,urs,dudrs,rc
type(slj)            :: g

g%e=e
g%s=s
g%rs=rs

! cutoff
s_dr=(s/rs)**6
urs=4._dp*e*(s_dr*s_dr-s_dr)
dudrs=-e*(48._dp*s_dr-24._dp)*s_dr/rs

rc=rs-32._dp*urs/(11._dp*dudrs)

call g%setrc(rc)
          
! other parameters
g%sa4=(8._dp*urs+(g%rcut-rs)*dudrs)/(4*(g%rcut-rs)**4)
g%sa8=-(4._dp*urs+(g%rcut-rs)*dudrs)/(4*(g%rcut-rs)**8)

end subroutine slj_set
  
subroutine slj_interact(g)
class(slj),intent(inout) :: g
integer                  :: i,ii,j,l,n
real(dp)                 :: vd(dm),dr,factor2(dm)
real(dp)                 :: p,f,s_dr
type(atom),pointer       :: o1,o2
type(atom_dclist),pointer   :: la
 
g%epot = 0._dp

la => g%ref%alist
do ii = 1,g%ref%nat
  la => la%next
  o1 => la%o
  i = o1%gid(g%id)

  do l = 1, g%nn(i)  ! sobre los vecinos

    j = g%list(i,l)
    o2 =>g%a(j)%o

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>g%rcut2) cycle

    dr=sqrt(dr)

    if(dr<g%rs) then

      s_dr=(g%s/dr)**6
      p=4._dp*g%e*(s_dr-1._dp)*s_dr
      f=g%e*(48._dp*s_dr-24._dp)*s_dr/dr      ! -dp/dr

    else
     
      s_dr=(dr-g%rcut)**4
      p=g%sa4*s_dr+g%sa8*s_dr*s_dr
      f=-(4._dp*g%sa4*s_dr+8._dp*g%sa8*s_dr*s_dr)/(dr-g%rcut)        ! -dp/dr
      
    endif

    p=p*ev_ui
    f=f*ev_ui

    o1%epot = o1%epot + p*0.5_dp
    g%epot = g%epot + p*0.5_dp

    factor2(1:dm) = f*vd(1:dm)/dr           ! (-dp/dr)*(-grad_1dr)
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 

    if(b_gvirial) then
      do n = 1,3
        virial(n,:) = virial(n,:) - o1%pos(n)*factor2(:) 
      enddo
    endif
    
    if (associated(o2%ghost)) then
      if (o2%ghost%gid(g%id)>o1%gid(g%id)) then
        ! Una sola vez por par para que se pueda calcular el virial

        ! print *, i,j,l
        o2%epot = o2%epot + p*0.5_dp
        g%epot = g%epot + p*0.5_dp
         
        o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
        if(b_gvirial) then
          do n = 1,3
            virial(n,:) = virial(n,:) + box(n)*floor(o2%pos(n)*one_box(n))*factor2(:) 
          enddo
        endif
      
      endif
    endif
         
  enddo     

enddo 

              
end subroutine  


               
subroutine wca_set(g,e,s)
! The Weeks-Chandler-Andersen potential: a Lennard-Jones potential shifted and
! cut at its minimum to retain a purely repulsive potential.
real(dp),intent(in)  :: e,s
type(wca)            :: g
                                 
! parametros
g%e=e
g%s=s
call g%setrc(2._dp**(1./6.)*s)
             
end subroutine wca_set
  
subroutine wca_interact(g)
class(wca),intent(inout) :: g
integer                  :: i,ii,j,n,l
real(dp)                 :: vd(dm),dr,factor2(dm)
real(dp)                 :: p,f,s_dr
type(atom),pointer       :: o1,o2
type(atom_dclist),pointer   :: la
 
! Por info ver LJ0.f90
g%epot = 0._dp

la => g%ref%alist
do ii = 1,g%ref%nat
  la => la%next
  o1 => la%o
  i = o1%gid(g%id)

  do l = 1, g%nn(i)  ! sobre los vecinos

    j  = g%list(i,l)
    o2 =>g%a(j)%o

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 


    if(dr>g%rcut2) cycle

    dr=sqrt(dr)

    s_dr=(g%s/dr)**6

    p=4._dp*g%e*(s_dr*s_dr-s_dr)+g%e
    f=g%e*(48._dp*s_dr-24._dp)*s_dr/dr         ! -dp/dr

    p=p*ev_ui
    f=f*ev_ui
         
    o1%epot = o1%epot + p*0.5_dp
    g%epot = g%epot + p*0.5_dp
      
    factor2(1:dm) = f*vd(1:dm)/dr            ! (-dp/dr)*(-grad_1dr)
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 

    if(b_gvirial) then
      do n = 1,3
        virial(n,:) = virial(n,:) - o1%pos(n)*factor2(:) 
      enddo
    endif
    
    if (associated(o2%ghost)) then
      if (o2%gid(g%id)>o1%gid(g%id)) then

        ! Una sola vez por par para que se pueda calcular el virial

        o2%epot = o2%epot + p*0.5_dp
        g%epot = g%epot + p*0.5_dp
         
        o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
        if(b_gvirial) then
          do n = 1,3
            virial(n,:) = virial(n,:) + box(n)*floor(o2%pos(n)*one_box(n))*factor2(:) 
          enddo
        endif
      
      endif
    endif
 
  enddo 

enddo     

end subroutine  
        


subroutine rwca_set(g,e,s,r)
! Reversed Weeks-Chandler-Andersen potential
real(dp),intent(in)  :: e,s,r
type(rwca)           :: g
g%e=e
g%s=s 
g%r=r*r

g%rcut=2._dp**(1./6.)*s+r
g%rcut2=g%rcut*g%rcut
              
end subroutine rwca_set
  
subroutine rwca_interact(g)
! Reversed Weeks-Chandler-Andersen potential
class(rwca),intent(inout) :: g
integer                   :: i,ii,j,l
real(dp)                  :: vd(dm),dr,factor2(dm)
real(dp)                  :: p,f,s_dr
type(atom),pointer        :: o1,o2
type(atom_dclist),pointer   :: la
 
! Por info ver LJ0.f90
g%epot = 0._dp

la => g%ref%alist
do ii = 1,g%ref%nat
  la => la%next
  o1 => la%o
  i = o1%gid(g%id)

  ! Sobre los vecinos
  do l = 1, g%nn(i) 

    j  = g%list(i,l)
    o2 => g%a(j)%o

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>g%rcut2) cycle
    if(dr<g%r) cycle

    dr=sqrt(dr)
    dr=g%rcut-dr
    
    s_dr=(g%s/dr)**6

    p=4._dp*g%e*(s_dr*s_dr-s_dr)+g%e
    p=p*ev_ui

    o1%epot = o1%epot + p*0.5_dp
    o2%epot = o2%epot + p*0.5_dp
    g%epot = g%epot + p

    f=g%e*(48._dp*s_dr-24._dp)*s_dr/dr   ! dp/dr (TODO: chequear signos)
    f=-f*ev_ui

    factor2(1:dm) = f*vd(1:dm)/dr              ! (-dp/dr)*(-grad_1dr)  (TODO: chequear signos)
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 
    o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
    
  enddo 

enddo     

end subroutine  
                                
 
subroutine cuw_set(g,ks,rs,radc,radi,eje,rc)
! Reversed cylindrically confined Weeks-Chandler-Andersen potential
real(dp),intent(in) :: ks,rs,radc,radi,rc
integer,intent(in)  :: eje
type(cuw)           :: g
g%ks=ks   
g%rs=rs*rs
g%radc=radc 
g%radi=radi 
g%eje=eje  
call g%setrc(rc)
end subroutine cuw_set
  
subroutine cuw_interact(g)
! Reversed cylindrically confined Weeks-Chandler-Andersen potential
use gems_algebra, only:fcut2_dfcut
class(cuw),intent(inout) :: g
integer                  :: i,ii,j,l
real(dp)                 :: vd(dm),dr,factor2(dm)
real(dp)                 :: p,proj,w,r1,r2,norm
type(atom),pointer       :: oi,oj
real(dp)                 :: aux,df1,df2
type(atom_dclist),pointer   :: la
 
g%epot = 0._dp

la => g%ref%alist
do ii = 1,g%ref%nat
  la => la%next
  oi => la%o
  i = oi%gid(g%id)

  ! Compute the normalization factor
  norm=0._dp
  p=0._dp
  do l = 1, g%nn(i)  ! sobre los vecinos

    j  = g%list(i,l)
    oj => g%a(j)%o

    vd = vdistance(oj,oi, mic)

    ! Skip this if the z distance is below the boundary
    if(vd(3)<g%rs) cycle

    dr = dot_product(vd,vd) 

    ! Skip using neighboor lsit
    if(dr>g%rcut2) cycle

    dr=sqrt(dr)
    proj=dr/vd(g%eje)

    ! Skip if the projection is outside the circle
    if(proj>g%radc+g%radi) cycle

    ! Compute the weight
    r1=proj-g%radi
    r2=proj+g%radi
    w=fcut2_dfcut(g%radc,r1,r2,df1)!-fcut2(-radc,r1,r2,df2)
    
    ! Accumulate and save the unnormalized weight
    norm=norm+w
    
  enddo 
    
  do l = 1, g%nn(i)  ! sobre los vecinos

    j  = g%list(i,l)
    oj => g%a(j)%o
        
    vd = vdistance(oj,oi, mic)

    ! Skip this if the z distance is below the boundary
    if(vd(3)<g%rs) cycle

    dr = dot_product(vd,vd) 

    ! Skip using neighboor lsit
    if(dr>g%rcut2) cycle

    dr=sqrt(dr)
    proj=dr/vd(g%eje)

    ! Skip if the projection is outside the circle
    if(proj>g%radc+g%radi) cycle
              
    ! Compute the weight
    r1=proj-g%radi
    r2=proj+g%radi

    !TODO: If radi<radc then fcut2(-radc,r1,r2,df2)=0 always
    w=fcut2_dfcut(g%radc,r1,r2,df1)-fcut2_dfcut(-g%radc,r1,r2,df2)
                 
    aux=(vd(3)-g%rs)
    p=w*g%ks*aux*aux*0.5_dp/norm

    oi%epot = oi%epot + p*0.5_dp
    oj%epot = oj%epot + p*0.5_dp
    g%epot = g%epot + p

    factor2(1:2) = (df1-df2)*vd(1:2)/dr
    factor2(1:2) = p*(1._dp/factor2(1:2)+factor2(1:2)/w)
    factor2(3) = -w*g%ks*aux*vd(3)/dr

    oi%force(1:3) = oi%force(1:3) + factor2(1:3) 
    oj%force(1:3) = oj%force(1:3) - factor2(1:3)
    
  enddo 

enddo     

end subroutine  
                          
      
subroutine sho_set(g,k,r)
type(sho)   :: g
real(dp)    :: k,r
g%k=k
g%r=r
call g%setrc(100._dp)
end subroutine        
      
subroutine sho_interact(g)
class(sho),intent(inout)    :: g
integer                     :: i,ii,j,l
real(dp)                    :: vd(dm),dr,factor2(dm)
real(dp)                    :: p,f
type(atom),pointer          :: o1,o2
type(atom_dclist),pointer   :: la
 
g%epot = 0._dp

la => g%ref%alist
do ii = 1,g%ref%nat
  la => la%next
  o1 => la%o
  i = o1%gid(g%id)

  do l = 1, g%nn(i)  ! sobre los vecinos
    j  = g%list(i,l)
    o2 => g%a(j)%o

    vd = vdistance( o1, o2 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>g%rcut2) cycle

    dr=sqrt(dr)-g%r

    p=-0.5_dp*g%k*dr*dr

    o1%epot = o1%epot + p*0.5_dp
    o2%epot = o2%epot + p*0.5_dp
    g%epot = g%epot + p

    f=-g%k*dr

    factor2(1:dm) = -f*vd(1:dm)   ! el dividido dr esta en f
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 
    o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
    
  enddo 

enddo     

end subroutine  
       
subroutine shocm_set(g,k,r)
type(shocm) :: g
real(dp)    :: k,r
g%k=k
g%r=r
call g%setrc(100._dp)
end subroutine        
 
subroutine shocm_interact(g)
class(shocm),intent(inout) :: g
real(dp)                   :: vd(dm),dr,paux,faux,m,m2
type(atom_dclist),pointer  :: la,alist
integer                    :: i
 
g%epot = 0._dp

vd = atom_dclist_inq_cmpos(g%alist,m)  &
    -atom_dclist_inq_cmpos(g%b%alist,m2)  ! Ojo el orden
dr = sqrt(dot_product(vd,vd))

paux = dr-g%r
faux = g%k*paux
paux = paux*faux*0.5_dp
vd=faux*vd/dr

! Sobre los dos grupos
alist => g%ref%alist
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
  alist => g%b%alist

enddo

! FIXME: g%b%nat does not count ghost?
g%epot = g%epot + paux*g%nat + paux*(g%b%nat)

end subroutine  



subroutine sm1_set(g,r1,r2,e)
! A smooth function. It has to reach cero (so the rcut hav sense)
type(sm1)           :: g
real(dp),intent(in) :: r1,r2,e
g%r1=r1
g%e=e
call g%setrc(r2)
end subroutine sm1_set
  
subroutine sm1_interact(g)
use gems_constants, only:pi,pio2
class(sm1),intent(inout) :: g
integer                  :: i,ii,j,l,n
real(dp)                 :: vd(dm),dr,factor2(dm)
real(dp)                 :: p,f
real(dp)                 :: r1,r2,aux
type(atom),pointer       :: o1,o2
type(atom_dclist),pointer   :: la
 
r1=g%r1
r2=g%rcut

g%epot = 0._dp

la => g%ref%alist
do ii = 1,g%ref%nat
  la => la%next
  o1 => la%o
  i = o1%gid(g%id)

  do l = 1, g%nn(i)  ! sobre los vecinos

    j  = g%list(i,l)
    o2 => g%a(j)%o

    vd = vdistance( o2, o1 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>g%rcut2) cycle

    dr=sqrt(dr)

    if(dr>r1) then
      aux=pi*(dr-(r1+r2)*0.5_dp)/(r2-r1)
      p=0.5_dp+0.5_dp*sin(aux)
      f=cos(aux)*pio2*dr/(r2-r1)      
      p=p*g%e-g%e
      f=f*g%e
    else
      p=-g%e
      f=0
    endif

    o1%epot = o1%epot + p*0.5_dp
    g%epot = g%epot + p*0.5_dp

    factor2(1:dm) = -f*vd(1:dm)/dr   
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 

    if(b_gvirial) then
      do n = 1,3
        virial(n,:) = virial(n,:) - o1%pos(n)*factor2(:) 
      enddo
    endif
          
    if (associated(o2%ghost)) then
      if (o2%gid(g%id)>o1%gid(g%id)) then

        o2%epot = o2%epot + p*0.5_dp
        g%epot = g%epot + p*0.5_dp
         
        o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
        if(b_gvirial) then
          do n = 1,3
            virial(n,:) = virial(n,:) + box(n)*floor(o2%pos(n)*one_box(n))*factor2(:) 
          enddo
        endif
      
      endif
    endif
         
  enddo     

enddo 
  
              
end subroutine  

    
end module gems_pairs
