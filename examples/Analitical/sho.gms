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
box size 20 20 20

# WARNING, LCG not recommended. 
# Use only to enforce reproducibility in test simulations.
prng lcg
prng seed 123456
 
>+ fill 10 0.25 $boxx$ $boxy$ 0.1   0. 0. 3.3  
set element H

>+ fill 10 0.25 0.1 0.1 $boxz$      5. 5. 0.
set element C

> all
set pbc T T T 
group 1 add
out posxyz 

> element H
group 2 add

> element C
group 3 add

((((((((((( caja )))))))))))


box move (make the total velocity zero)

interact 2 field sho_plane 1. 3.3 3
interact 3 field sho_line  1.  5. 5. 0.   0. 0. 1

time step 0.001
                
                        
outfile :f1 name Energy.$jobname$.dat
outfile :f1 cols energy 1
outfile :f1 each 100

# outfile :f3 name Pos.$jobname$.xyz
# outfile :f3 pos 1
# outfile :f3 each 1

> all
evolve v_verlet
set tempgdist 300
dinamica 1000

