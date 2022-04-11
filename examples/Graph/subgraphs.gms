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
 
>< fillpbc 10 1.0 10. 10. 10. 0. 0. 0. 
+< fillpbc 5 1.0 {$boxx-12} {$boxy-12} {$boxz-12} 12. 12. 12. 
set element H

> sys
group 1 add
set pbc T T T
out posxyz asd.xyz

box move (make the total velocity zero)

interact :gr 1 graph subgraphs 3.

outfile :f1 name Graph.$jobname$.dat
outfile :f1 graph 1

out state

