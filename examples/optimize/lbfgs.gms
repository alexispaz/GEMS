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
>+ atom  0      0      0    
^+ atom  0.5    0.5    0
^+ atom  0    0.5    0.5  
^+ atom  0.5    0    0.5  
                          
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
d0:={sqrt(0.5)}

# All the surfaces exposed are (100). 
 
^+ reply   1 0 0     5
^+ reply   0 1 0     5
^+ reply   0 0 1     5
set element Ag
box size 12 12 12

> all
set move 1 1 1 
# fact:={2.89/$d0$}
fact:={3.89/$d0$}
set expand $fact$ $fact$ $fact$
box expand $fact$ $fact$ $fact$
# set pbc T T T
group 1 add

interact 1 tb read ../../parameters/baletto_2003.prm

outfile :f1 name Energy.$jobname$.dat
outfile :f1 cols energy 1
outfile :f1 each 20
 
# outfile :f2 name Pos.$jobname$.xyz
# outfile :f2 pos 1
# outfile :f2 each 20

lbfgs

