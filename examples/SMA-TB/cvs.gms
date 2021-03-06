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

interact 1 under tb read ../../parameters/baletto_2003.prm

outfile :f1 name Energy.$jobname$.dat
outfile :f1 cols energy 1
outfile :f1 each 100
                
# outfile :f2 name Pos.$jobname$.xyz
# outfile :f2 Pos 1
# outfile :f2 each 10
 
lbfgs
out posxyz asd.xyz

> sys
set cm_vel 1. 0. 0.
evolve v_verlet 
cv cm xyz 5.
                
dinamica 1000

