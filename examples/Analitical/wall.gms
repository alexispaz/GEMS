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
 
>< fill 10 0.25 $boxx $boxy 0.1  0. 0. 3.3  
set element H

>< fill 10 0.25 0.1 0.1 $boxz    5. 5. 0.
set element C

> sys
group 1 add
out posxyz 

> element H
set pbc T T F
group 2 add

> element C
set pbc T T T
group 3 add

((((((((((( caja )))))))))))


box move (make the total velocity zero)


interact :pared 2 field halfsho_plane 1. 3.0 -3
interact 2 field halfsho_plane 1. $boxz 3

bloque save 1
  r0 := {3.0+(1.*$time)}
interact :pared 1. $r0$  -3
fin


bloque load_silent T
bloque load_each 10 1
 
time step 0.001
                
                        
outfile :f1 name Energy.$jobname$.dat
outfile :f1 cols energy 1 e_halfsho
outfile :f1 each 100

outfile :f3 name Pos.$jobname$.xyz
outfile :f3 pos 1
outfile :f3 each 1

> sys
evolve v_verlet
set tempgdist 300
dinamica 5000

