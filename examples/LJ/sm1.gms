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

>< atom 1 1 1
+< atom 1 1 1
# m0=1uma
set element H

> sys
set pbc F F F 
group 1 add
 
> atom 2  
group 2 add

((((((((((( caja )))))))))))


box move (make the total velocity zero)

# Epsilon=1eV=1E0
# Rmin=1A=1r0
interact 1 pair sm1 1. 2. 0.01

time step 0.001
                
                        
outfile :f1 name Energy.$jobname$.dat
outfile :f1 cols energy 1
outfile :f1 each 1

# outfile :f3 name Pos.$jobname$.xyz
# outfile :f3 pos 1
# outfile :f3 each 1
         

> sys

out state

bloque repeat 300
> group 2
set move 0.01 0 0
> sys
out state
# lbfgs
fin

> group 2  
set move -2 0 0
set cm_vel 10 0 0

> sys
evolve v_verlet 

dinamica 1000

