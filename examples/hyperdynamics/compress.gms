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
 
>+ read coords_min.xyz
set element Ag

> element Ag
group 1 add

time step 0.001d0  (integration timestep [ps])

interact 1 tb read ../../parameters/baletto_2003.prm
interact 1 bias compress 0.8
 
# outfile :f1 name E.$jobname$.dat
# outfile :f1 cols time e_reax
# outfile :f1 each 10

outfile :f2 name Energy.$jobname$.dat
outfile :f2 cols time energy 1 bias
outfile :f2 each 10
             
> all
evolve :termo ermak 300 1
evolve :termo voter .true.
dinamica 300
 
