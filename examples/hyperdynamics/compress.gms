(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))
(                    INPUT SYSTEM FOR GeMS PROGRAM                 )
(                   comentaries are between brackets               )
(                      caps can be lower or upper                  )
(               caps must be upper for potential names             )
(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))


dimension 3
 
# WARNING, LCG not recommended. 
# Use only to enforce reproducibility in test simulations.
prng lcg
prng seed 123456
 
>< read coords_min.xyz
sys add Ag

> element Ag
group 1 add

time step 0.001d0  (integration timestep [ps])

interact 1 under tb read ../../parameters/baletto_2003.prm
interact 1 under bias compress 0.8
 
# outfile 1 name E.$jobname$.dat
# outfile 1 cols time e_reax
# outfile 1 each 10

outfile 2 name Energy.$jobname$.dat
outfile 2 cols time energy 1 bias
outfile 2 each 10
             
> sys
evolve :termo ermak 300 1
evolve :termo voter .true.
dinamica 300
 
