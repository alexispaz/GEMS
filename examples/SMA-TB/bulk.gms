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

# Primitive cell for the generation of a face-centered diamond-cubic
# (fcdc) crystal structure.
>< atom  0      0      0    
+< atom  0.5    0.5    0
+< atom  0    0.5    0.5  
+< atom  0.5    0    0.5  
                          
#    y
#    ^
#    |              
# 0.5+ 0.5--------0
#    |  |         |
#    |  |         |
#    |  |         |
#    |  |         |
#   0+  0--------0.5   
#    |         
#    +--+---------+----> x
#       0         0.5         
                          
# The distance between first neighbors.
getin d0 {sqrt(0.5)}

# All the surfaces exposed are (100). 
 
+< reply   1 0 0     25
+< reply   0 1 0     25
+< reply   0 0 1     25
sys add Ag
box size 30 30 30

> sys
set move 2 2 2 
getin fact {2.89/$d0$}
set expand $fact$ $fact$ $fact$
box expand $fact$ $fact$ $fact$
# set pbc T T T
group 1 add

interact 1 under tb read ../../parameters/baletto_2003.prm

outfile 1 name Energy.$jobname$.dat
outfile 1 cols energy 1
outfile 1 each 10
 
outfile 2 name Pos.$jobname$.xyz
outfile 2 pos 1
outfile 2 each 50

lbfgs

