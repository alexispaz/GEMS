(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))
(                    INPUT SYSTEM FOR GeMS PROGRAM                 )
(                   comentaries are between brackets               )
(                      caps can be lower or upper                  )
(               caps must be upper for potential names             )
(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))


# This LJ parameters was used in:  
# Grønbech-Jensen, N., & Farago, O. (2014). Constant pressure and
# temperature discrete-time Langevin molecular dynamics. The Journal of
# Chemical Physics, 141(19), 194108. http://doi.org/10.1063/1.4901303


(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))
(                    Estableciendo la  configuracion               )
(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))


dimension 3

# WARNING, LCG not recommended. 
# Use only to enforce reproducibility in test simulations.
prng lcg
 
# WARNING: Comment this line when doing your simulations. This is used here
# only to faciliate the code test
prng seed {123456*$pc$}

box size 5.78 5.78 5.78

>< atom 1 1 1
+< atom 2.2 1 1
# m0=1uma
sys add H

> sys
set move 1. 1. 1.
set pbc F F F 
set mass 3.
group 1 add
 
> atom 2  
group 2 add

((((((((((( caja )))))))))))


box move (make the total velocity zero)

# Epsilon=1eV=1E0
# Rmin=1A=1r0
interact 1 lj 1. 0.890898718140339 2

# About time unit
#   1eV=9648.61 uma*A**2/ps**2
#   1E0=9648.61 m0*r0**2/ps**2
#   1ps**2=  9648.61 m0*r0**2/E0
#   1ps=sqrt(9648.61) tau 
#   1ps=98.22 tau
tau_ps:=0.01018
       
time step {0.002*$tau_ps$}   (integration timestep [ps])

> sys

# MPI example 1: Note the variable substitution using rank
temp:={$pc$*5.+30}
set tempgdist $temp$

# MPI example 2: Note the code executed by a certain MPI rank
mpi only 0 evolve ermak 30 {1/$tau_ps$} 
mpi only 1 evolve ermak 35 {1/$tau_ps$} 
mpi only 2 evolve ermak 40 {1/$tau_ps$}  
mpi only 3 evolve ermak 45 {1/$tau_ps$} 

dinamica 1000                   
                        
outfile :f1 name Energy.$jobname$.dat
outfile :f1 cols energy 1 
outfile :f1 each 10
# This should be in 1
# plot [][:] 'Energy.gms.dat' u ($0*100*0.002):($4*kB_eV) w l 
                         
# outfile :f2 name Temp.$jobname$.dat
# outfile :f2 cols temp 1
# outfile :f2 each 10
         
# outfile :f3 name Pos.$jobname$.xyz
# outfile :f3 pos 1
# outfile :f3 each 1
# outfile :f3 flush on
         

> sys
out state
partemp2 100 100

