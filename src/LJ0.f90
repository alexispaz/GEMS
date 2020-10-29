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

 
module gems_pair_pot0
use gems_program_types
use gems_algebra 
use gems_constants,     only: dp,ev_ui,ui_ev,dm
use gems_inq_properties
use gems_neighbour

implicit none

private
public sig_arule,eps_arule
public rclj,sho_new,sho_cli

real(dp)  :: rclj=20._dp
 
logical           :: sig_arule,eps_arule
   
integer,parameter         :: mm=10
       
! Lennard Jones with parabolic smooth function
public plj_set
real(dp),parameter        :: ljcut=10._dp
       
! Lennard Jones
public lj_set

! Weeks-Chandler-Andersen potential
public wca_set
        
! Reversed Weeks-Chandler-Andersen potential
public rwca_set
integer                   :: rwca_n=0
type(intergroup),target   :: rwcagr(mm)
real(dp),dimension(mm)    :: rtsig,rteps,rtrmi
         
! Reversed cylindrically confined Weeks-Chandler-Andersen potential
public cuw_set
       

! Lennard Jones with splined smooth function
! Grønbech-Jensen, N., & Farago, O. (2014). Constant pressure and temperature
! discrete-time Langevin molecular dynamics. The Journal of Chemical Physics,
! 141(19), 194108. http://doi.org/10.1063/1.4901303
public slj_set

public sm1_set
        

! SHO potential
real(dp),allocatable        :: shofixr(:,:)
public shofix_interaction,shofix_set


! ! SHO centro de masas
! public shocm_interaction

contains


subroutine lj_set(g1,e,s,rc,g2)
! Lennard Jones with parabolic smooth function
real(dp),intent(in)             :: e,s,rc
type(group),intent(in)          :: g1
type(group),intent(in),optional :: g2
type(intergroup),pointer        :: igr
     
! Initialize the integroup
allocate(igr)
if(present(g2)) then
  call igr%init(g1=g1,g2=g2)
else
  call igr%init(g1=g1)
endif
igr%interact => lj_interaction
                
! if(item<nitems) then
!   call readl(w1)
!   if(w1=='intermolecular') then
!     if(under) then
!       igr(nigr)%lista => verlet_self_crossmol
!     else
!       igr(nigr)%lista => verlet_cross_crossmol
!     endif
!   else
!     call wwan('I do not understand the last command')
!   endif
! endif
                                     
! parametros
call igr%setrc(rc)
call igr%p%put(1,e)
call igr%p%put(2,s)
                           
end subroutine lj_set
  
subroutine lj_interaction(ig)
class(intergroup),intent(inout) :: ig
integer                         :: i,j,n,l
real(dp)                        :: vd(dm),dr,factor2(dm)
real(dp)                        :: p,f,s_dr,e,s
type(atom),pointer              :: o1,o2
 
! Parameters
e  =ig%p%o(1)
s  =ig%p%o(2)

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

    s_dr=(s/dr)**6

    p=4._dp*e*(s_dr*s_dr-s_dr)
    f=e*(48._dp*s_dr-24._dp)*s_dr/dr   ! -dp/dr

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

end subroutine lj_interaction 
                             

subroutine plj_set(g1,e,s,g2)
! Lennard Jones with parabolic smooth function
real(dp),intent(in)         :: e,s
type(group),intent(in)          :: g1
type(group),intent(in),optional :: g2
real(dp)                        :: prc,frc
type(intergroup),pointer        :: igr
     
! Initialize the integroup
allocate(igr)
if(present(g2)) then
  call igr%init(ljcut,g1=g1,g2=g2)
else
  call igr%init(ljcut,g1=g1)
endif
igr%interact => plj_interaction
 
! if(item<nitems) then
!   call readl(w1)
!   if(w1=='intermolecular') then
!     if(under) then
!       igr(nigr)%lista => verlet_self_crossmol
!     else
!       igr(nigr)%lista => verlet_cross_crossmol
!     endif
!   else
!     call wwan('I do not understand the last command')
!   endif
! endif
                                     
! parametros
call igr%p%put(1,e)
call igr%p%put(2,s)

prc=4._dp*e*((s/ljcut)**12-(s/ljcut)**6)
frc=24._dp*e*(2._dp*(s**12/ljcut**13)-(s**6/ljcut**7))
call igr%p%put(3,prc)
call igr%p%put(4,frc)
                  
end subroutine plj_set
  
subroutine plj_interaction(ig)
class(intergroup),intent(inout) :: ig
integer                         :: i,j,n,l
real(dp)                        :: vd(dm),dr,factor2(dm)
real(dp)                        :: p,f,s_dr,e,s,prc,frc
type(atom),pointer              :: o1,o2
 
! Parameters
e  =ig%p%o(1)
s  =ig%p%o(2)
prc=ig%p%o(3)
frc=ig%p%o(4)

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

    s_dr=(s/dr)**6

    p=4._dp*e*(s_dr*s_dr-s_dr)-prc-frc*(dr-ljcut)
    f=e*(48._dp*s_dr-24._dp)*s_dr/dr+frc        ! (-dp/dr)

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

end subroutine plj_interaction 

                               

subroutine slj_set(g1,e,s,rs,g2)
! Lennard Jones with splined smooth function
! The equations for this potential can be found in
! Grønbech-Jensen, N., & Farago, O. (2014). Constant pressure and temperature
! discrete-time Langevin molecular dynamics. The Journal of Chemical Physics,
! 141(19), 194108. http://doi.org/10.1063/1.4901303

real(dp),intent(in)         :: e,s,rs
real(dp)                    :: s_dr,urs,dudrs,sa4,sa8,rc
type(group),intent(in)          :: g1
type(group),intent(in),optional :: g2
type(intergroup),pointer    :: igr

! Initialize the integroup
allocate(igr)
if(present(g2)) then
  call igr%init(g1=g1,g2=g2)
else
  call igr%init(g1=g1)
endif
igr%interact => slj_interaction   
       
! if(item<nitems) then
!   call readl(w1)
!   if(w1=='intermolecular') then
!     if(under) then
!       igr(nigr)%lista => verlet_self_crossmol
!     else
!       igr(nigr)%lista => verlet_cross_crossmol
!     endif
!   else
!     call wwan('I do not understand the last command')
!   endif
! endif
                                      
! parametros
call igr%p%put(1,e)
call igr%p%put(2,s)
call igr%p%put(3,rs)

! cutoff
s_dr=(s/rs)**6
urs=4._dp*e*(s_dr*s_dr-s_dr)
dudrs=-e*(48._dp*s_dr-24._dp)*s_dr/rs

rc=rs-32._dp*urs/(11._dp*dudrs)

call igr%setrc(rc)
          
! other parameters
sa4=(8._dp*urs+(igr%rcut-rs)*dudrs)/(4*(igr%rcut-rs)**4)
sa8=-(4._dp*urs+(igr%rcut-rs)*dudrs)/(4*(igr%rcut-rs)**8)
call igr%p%put(4,sa4)
call igr%p%put(5,sa8)

end subroutine slj_set
  
subroutine slj_interaction(ig)
class(intergroup),intent(inout) :: ig
integer                         :: i,j,l,n
real(dp)                        :: vd(dm),dr,factor2(dm)
real(dp)                        :: p,f,s_dr
real(dp)                        :: e,s,rs,sa4,sa8
type(atom),pointer              :: o1,o2
 
! Parameters
e  =ig%p%o(1)
s  =ig%p%o(2)
rs =ig%p%o(3)
sa4=ig%p%o(4)
sa8=ig%p%o(5)

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

    if(dr<rs) then

      s_dr=(s/dr)**6
      p=4._dp*e*(s_dr-1._dp)*s_dr
      f=e*(48._dp*s_dr-24._dp)*s_dr/dr      ! -dp/dr

    else
     
      s_dr=(dr-ig%rcut)**4
      p=sa4*s_dr+sa8*s_dr*s_dr
      f=-(4._dp*sa4*s_dr+8._dp*sa8*s_dr*s_dr)/(dr-ig%rcut)        ! -dp/dr
      
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


               
subroutine wca_set(g1,e,s,g2)
! The Weeks-Chandler-Andersen potential: a Lennard-Jones potential shifted and
! cut at its minimum to retain a purely repulsive potential.
real(dp),intent(in)         :: e,s
type(group),intent(in)          :: g1
type(group),intent(in),optional :: g2
type(intergroup),pointer    :: igr

! Initialize the integroup
allocate(igr)
if(present(g2)) then
  call igr%init(g1=g1,g2=g2)
else
  call igr%init(g1=g1)
endif
igr%interact => wca_interaction   
       
! if(item<nitems) then
!   call readl(w1)
!   if(w1=='intermolecular') then
!     if(under) then
!       igr(nigr)%lista => verlet_self_crossmol
!     else
!       igr(nigr)%lista => verlet_cross_crossmol
!     endif
!   else
!     call wwan('I do not understand the last command')
!   endif
! endif
                                      
! parametros
                                  
! parametros
call igr%p%put(1,e)
call igr%p%put(2,s)
call igr%setrc(2._dp**(1./6.)*s)
             
end subroutine wca_set
  
subroutine wca_interaction(ig)
class(intergroup),intent(inout) :: ig
integer                         :: i,j,n,l
real(dp)                        :: vd(dm),dr,factor2(dm)
real(dp)                        :: p,f,e,s,s_dr
type(atom),pointer              :: o1,o2
 
! Parameters
e  =ig%p%o(1)
s  =ig%p%o(2)
     
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

    s_dr=(s/dr)**6

    p=4._dp*e*(s_dr*s_dr-s_dr)+e
    f=e*(48._dp*s_dr-24._dp)*s_dr/dr         ! -dp/dr

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
        


subroutine rwca_set(g1,e,s,r,g2)
! The Weeks-Chandler-Andersen potential: a Lennard-Jones potential shifted and
! cut at its minimum to retain a purely repulsive potential.
real(dp),intent(in)             :: e,s,r
type(group),intent(in)          :: g1
type(group),intent(in),optional :: g2
type(intergroup),pointer        :: igr

rwca_n=rwca_n+1
igr => rwcagr(rwca_n)
        
! Init the igr
if(present(g2)) then
  call igr%init(g1=g1,g2=g2)
else
  call igr%init(g1=g1)
endif
igr%interact => rwca_interaction
                
! if(item<nitems) then
!   call readl(w1)
!   if(w1=='intermolecular') then
!     if(under) then
!       igr(nigr)%lista => verlet_self_crossmol
!     else
!       igr(nigr)%lista => verlet_cross_crossmol
!     endif
!   else
!     call wwan('I do not understand the last command')
!   endif
! endif
                                      
! parametros
rteps(rwca_n)=e
rtsig(rwca_n)=s 
rtrmi(rwca_n)=r*r

igr%rcut=2._dp**(1./6.)*s+r
igr%rcut2=igr%rcut*igr%rcut
              
end subroutine rwca_set
  
subroutine rwca_interaction(ig)
class(intergroup),intent(inout) :: ig
integer                         :: i,j,l,m
real(dp)                        :: vd(dm),dr,factor2(dm)
real(dp)                        :: p,f,s_dr
type(atom),pointer              :: o1,o2
 
m=ig%id

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
    if(dr<rtrmi(m)) cycle

    dr=sqrt(dr)
    dr=ig%rcut-dr
    
    s_dr=(rtsig(m)/dr)**6

    p=4._dp*rteps(m)*(s_dr*s_dr-s_dr)+rteps(m)
    p=p*ev_ui

    o1%epot = o1%epot + p*0.5_dp
    o2%epot = o2%epot + p*0.5_dp
    ig%epot = ig%epot + p

    f=rteps(m)*(48._dp*s_dr-24._dp)*s_dr/dr   ! dp/dr (TODO: chequear signos)
    f=-f*ev_ui

    factor2(1:dm) = f*vd(1:dm)/dr              ! (-dp/dr)*(-grad_1dr)  (TODO: chequear signos)
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 
    o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
    
  enddo 

enddo     

end subroutine  
                                
 
subroutine cuw_set(g1,ks,rs,radc,radi,eje,rc,g2)
!Cilindrical sectioned upper wall
real(dp),intent(in)             :: ks,rs,radc,radi,rc
integer,intent(in)              :: eje
type(group),intent(in)          :: g1
type(group),intent(in),optional :: g2
type(intergroup),pointer        :: igr

! Initialize the integroup
allocate(igr)
if(present(g2)) then
  call igr%init(g1=g1,g2=g2)
else
  call igr%init(g1=g1)
endif
igr%interact => cuw_interaction
                               
! parametros
call igr%p%put(1,ks    )
call igr%p%put(2,rs*rs )
call igr%p%put(3,radc  )
call igr%p%put(4,radi  )
call igr%i%put(1,eje   )
    
call igr%setrc(rc)
          
end subroutine cuw_set
  
subroutine cuw_interaction(ig)
use gems_algebra, only:fcut2_dfcut
class(intergroup),intent(inout) :: ig
integer                         :: i,j,l,m
real(dp)                        :: vd(dm),dr,factor2(dm)
real(dp)                        :: p,proj,w,r1,r2,norm
type(atom),pointer              :: oi,oj
real(dp)                        :: radc,radi,ks,rs,aux,df1,df2
integer                         :: eje
 
m=ig%id
   
! Parameters
ks =ig%p%o(1)
rs =ig%p%o(2)
radc=ig%p%o(3) ! The radious of the cilinder
radi=ig%p%o(4) ! The radious of each particle
eje=ig%i%o(1)

      
! Por info ver LJ0.f90
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
    if(vd(3)<rs) cycle

    dr = dot_product(vd,vd) 

    ! Skip using neighboor lsit
    if(dr>ig%rcut2) cycle

    dr=sqrt(dr)
    proj=dr/vd(eje)

    ! Skip if the projection is outside the circle
    if(proj>radc+radi) cycle

    ! Compute the weight
    r1=proj-radi
    r2=proj+radi
    w=fcut2_dfcut(radc,r1,r2,df1)!-fcut2(-radc,r1,r2,df2)
    
    ! Accumulate and save the unnormalized weight
    norm=norm+w
    
  enddo 
    
  do l = 1, ig%nn(i)  ! sobre los vecinos

    j  = ig%list(i,l)
    oj => ig%at(j)%o
        
    vd = vdistance(oj,oi, mic)

    ! Skip this if the z distance is below the boundary
    if(vd(3)<rs) cycle

    dr = dot_product(vd,vd) 

    ! Skip using neighboor lsit
    if(dr>ig%rcut2) cycle

    dr=sqrt(dr)
    proj=dr/vd(eje)

    ! Skip if the projection is outside the circle
    if(proj>radc+radi) cycle
              
    ! Compute the weight
    r1=proj-radi
    r2=proj+radi

    !TODO: If radi<radc then fcut2(-radc,r1,r2,df2)=0 always
    w=fcut2_dfcut(radc,r1,r2,df1)-fcut2_dfcut(-radc,r1,r2,df2)
                 
    aux=(vd(3)-rs)
    p=w*ks*aux*aux*0.5_dp/norm

    oi%epot = oi%epot + p*0.5_dp
    oj%epot = oj%epot + p*0.5_dp
    ig%epot = ig%epot + p

    factor2(1:2) = (df1-df2)*vd(1:2)/dr
    factor2(1:2) = p*(1._dp/factor2(1:2)+factor2(1:2)/w)
    factor2(3) = -w*ks*aux*vd(3)/dr

    oi%force(1:3) = oi%force(1:3) + factor2(1:3) 
    oj%force(1:3) = oj%force(1:3) - factor2(1:3)
    
  enddo 

enddo     

end subroutine  
                          
                     
function sho_new(w,g1,g2) result(igr)
! TODO: Generalize for any plane
type(group),intent(in)     :: g1
type(group),intent(in),optional    :: g2
type(intergroup),pointer   :: igr
character(*),intent(in)    :: w

! Bulid the igr
allocate(igr)
if(present(g2)) then
  call igr%init(g1=g1,g2=g2)
else
  call igr%init(g1=g1)
endif
 
select case(w)
case('sho')
  igr%interact => sho_interaction
case('shocm')
  igr%interact => shocm_interaction
case default
  call werr('SHO type not found')
endselect

end function

subroutine sho_cli(igr)
use gems_input_parsing
type(intergroup)           :: igr
real(dp)                   :: f1,f2

call readf(f1)
call readf(f2)

call igr%p%put(1,f1)
call igr%p%put(2,f2)

call igr%setrc(100._dp)

end subroutine        

subroutine sho_interaction(ig)
class(intergroup),intent(inout)         :: ig
integer     :: i,j,l,m
real(dp)    :: vd(dm),dr,factor2(dm)
real(dp)    :: p,f,k,r
type(atom),pointer          :: o1,o2
 
k=ig%p%o(1)
r=ig%p%o(2)

ig%epot = 0._dp

do i = 1,ig%n(1)
  o1 => ig%at(i)%o

  do l = 1, ig%nn(i)  ! sobre los vecinos
    j  = ig%list(i,l)
    o2 => ig%at(j)%o

    vd = vdistance( o1, o2 , mic) ! respetar el orden
    dr = dot_product(vd,vd) 

    if(dr>ig%rcut2) cycle

    dr=sqrt(dr)-r

    p=-0.5_dp*k*dr*dr

    o1%epot = o1%epot + p*0.5_dp
    o2%epot = o2%epot + p*0.5_dp
    ig%epot = ig%epot + p

    f=-k*dr

    factor2(1:dm) = -f*vd(1:dm)   ! el dividido dr esta en f
    o1%force(1:dm) = o1%force(1:dm) - factor2(1:dm) 
    o2%force(1:dm) = o2%force(1:dm) + factor2(1:dm)
    
  enddo 

enddo     

end subroutine  

subroutine shocm_interaction(ig)
class(intergroup),intent(inout)         :: ig
real(dp)                  :: vd(dm),dr,paux,faux,m,m2,k,r
type(atom_dclist),pointer :: la,alist
integer                   :: i
 
k=ig%p%o(1)
r=ig%p%o(2)

ig%epot = 0._dp

vd = atom_dclist_inq_cmpos(ig%a,m)  &
    -atom_dclist_inq_cmpos(ig%b,m2)  ! Ojo el orden
dr = sqrt(dot_product(vd,vd))

paux = dr-r
faux = k*paux
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



subroutine shofix_set(g1,s,g2)
type(group),intent(in)          :: g1
type(group),intent(in),optional :: g2
real(dp),intent(in)         :: s
integer                     :: k
type(atom_dclist),pointer   :: la
type(intergroup),pointer    :: igr

! Initialize the integroup
allocate(igr)
if(present(g2)) then
  call igr%init(g1=g1,g2=g2)
else
  call igr%init(g1=g1)
  igr%lista => intergroup0_empty
endif
igr%interact => shofix_interaction
    
! parametros
call igr%p%put(1,s)
! FIXME: call igr%p%put(1,r)

allocate(shofixr(igr%n(1),3))

la => igr%a      ! sobre los atomos
do k = 1,igr%n(1)
  la=>la%next
  shofixr(k,:)=la%o%pos(:)
enddo     
             
end subroutine shofix_set
  
subroutine shofix_interaction(ig)
class(intergroup),intent(inout)         :: ig
integer     :: i,m
real(dp)    :: vd(dm),dr,factor2(dm)
real(dp)    :: p,f,k,r
type(atom_dclist),pointer :: la
 
k=ig%p%o(1)

! FIXME: out of bound  
r=ig%p%o(2)

ig%epot = 0._dp

la => ig%a      ! sobre los atomos
do i = 1,ig%n(1)
  la=>la%next

  vd = la%o%distance(shofixr(i,:))  ! vector a-->p
  dr = dot_product(vd,vd) 

  dr=sqrt(dr)-r

  p=-0.5_dp*k*dr*dr
     
  la%o%epot = la%o%epot + p*0.5_dp
  ig%epot = ig%epot + p

  f=-k*dr

  factor2(1:dm)=f*vd(1:dm)   ! el dividido dr esta en f

  la%o%force(1:dm) = la%o%force(1:dm) - factor2(1:dm) 
  
enddo     

end subroutine  


      

subroutine sm1_set(g1,r1,r2,e,g2)
  ! A smooth function. It has to reach cero (so the rcut hav sense)
  real(dp),intent(in)             :: r1,r2,e
  type(group),intent(in)          :: g1
  type(group),intent(in),optional :: g2
  type(intergroup),pointer        :: igr


  ! Initialize the integroup
  allocate(igr)
  if(present(g2)) then
    call igr%init(g1=g1,g2=g2)
  else
    call igr%init(g1=g1)
  endif
  igr%interact => sm1_interaction   
                             
  ! parametros
  call igr%p%put(1,r1)
  call igr%p%put(2,e)

  ! cutoff
  call igr%setrc(r2)

end subroutine sm1_set
  
subroutine sm1_interaction(ig)
use gems_constants, only:pi,pio2
class(intergroup),intent(inout) :: ig
integer                         :: i,j,l,n
real(dp)                        :: vd(dm),dr,factor2(dm)
real(dp)                        :: p,f
real(dp)                        :: e,r1,r2,aux
type(atom),pointer              :: o1,o2
 
! Parameters
r1=ig%p%o(1)
r2=ig%rcut
e=ig%p%o(2)

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
      p=p*e-e
      f=f*e
    else
      p=-e
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

    
end module gems_pair_pot0
