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

>< read coords.xyz
sys add Ag

> element Ag
group 1 add

time step 0.001d0  (integration timestep [ps])

interact 1 under tb read ../../parameters/baletto_2003.prm  

                
outfile 1 name Energy.$jobname$.dat
outfile 1 cols energy 1
outfile 1 each 1
                 
outfile 2 name Eprom.$jobname$.dat
outfile 2 cols energy 1
outfile 2 each 1
outfile 2 prom 10
                 
outfile 3 name Eddda.$jobname$.dat
outfile 3 cols energy 1
outfile 3 each 1
outfile 3 ddda 10
 
> sys
evolve ermak 300 1
                
dinamica 100
              
