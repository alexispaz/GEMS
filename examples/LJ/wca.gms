(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))
(                    INPUT SYSTEM FOR GeMS PROGRAM                 )
(                   comentaries are between brackets               )
(                      caps can be lower or upper                  )
(               caps must be upper for potential names             )
(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))




dimension 3
box size 5.78 5.78 5.78

>+ atom 1 1 1
^+ atom 2.2 1 1
set element H

> all
set pbc F F F 
group 1 add
 
> atom 2  
group 2 add

> atom 1  
group 3 add


box move (make the total velocity zero)

# Epsilon=1eV=1E0
# Rmin=1A=1r0
interact 2 < 3 pair wca 1. 0.890898718140339

# About time unit
#   1eV=9648.61 uma*A**2/ps**2
#   1E0=9648.61 m0*r0**2/ps**2
#   1ps**2=  9648.61 m0*r0**2/E0
#   1ps=sqrt(9648.61) tau 
#   1ps=98.22 tau
tau_ps:=0.01018
       
time step {0.002*$tau_ps$}   (integration timestep [ps])
                
                        
outfile :f1 name Energy.$jobname$.dat
outfile :f1 cols energy 1
outfile :f1 each 10
# This should be in 1
# plot [][:] 'Energy.gms.dat' u ($0*100*0.002):($4*kB_eV) w l 

outfile :f3 name Pos.$jobname$.xyz
outfile :f3 pos 1
outfile :f3 each 1
         

> all

out state

> group 2
set move -1.2 0 0


# This should give 1./2. (half epsilon, because of the feel command)
set move 0.890898718140339 0 0
out state
set move -0.890898718140339 0 0

set move 0.9 0 0

bloque repeat 250
> group 2
set move 0.01 0 0
> all
out state
# lbfgs
fin

> group 2  
set move -2.5 0 0
set move -0.9 0 0
set move 0.95 0 0

> all
evolve v_verlet
dinamica 1000

