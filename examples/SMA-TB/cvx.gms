(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))
(                    INPUT SYSTEM FOR GeMS PROGRAM                 )
(                   comentaries are between brackets               )
(                      caps can be lower or upper                  )
(               caps must be upper for potential names             )
(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))


dimension 3

>< read coords.xyz
sys add Ag

> element Ag
group 1 add

time step 0.001d0  (integration timestep [ps])

interact 1 under tb ../../parameters/baletto_2003.prm

outfile 1 name Energy.$jobname$.dat
outfile 1 cols energy 1
outfile 1 each 100

# outfile 2 name Pos.$jobname$.xyz
# outfile 2 Pos 1
# outfile 2 each 1
 
lbfgs
out posxyz asd.xyz

> sys
set cm_vel 1. 1. 0.1
evolve v_verlet 
cv cm x 50.
cv cm y 50.
                
dinamica 1000

