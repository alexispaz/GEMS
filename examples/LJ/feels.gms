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
box size 5.78 5.78 5.78

>< atom 1 1 1
+< atom 2.2 1 1
+< atom 1 1 5
+< atom 2.2 1 5
+< atom 3.4 1 1
# m0=1uma
sys add H

> sys
set pbc F F F 
group 1 add
 
> sys
- atom 5
group 2 add
 
> atom 2 
+ atom 4
group 3 add

> atom 5
group 4 add

((((((((((( caja )))))))))))


box move (make the total velocity zero)

# Epsilon=1eV=1E0
# Rmin=1A=1r0
interact 2 pair slj 1. 0.890898718140339 1.108683
interact 4 < 3 pair slj 1. 0.890898718140339 1.108683

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

# outfile :f3 name Pos.$jobname$.xyz
# outfile :f3 pos 1
# outfile :f3 each 1
         

> sys

out state
evolve v_verlet
dinamica 1000

