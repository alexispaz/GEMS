# Input file for GeMS

getin pi   3.14159265358979
getin s60  {sin(60*$pi$/180)}
getin c60  0.5

# Primitive cell for the generation of a graphene/graphite
>< atom   0.00      $s60$  0.00
+< atom {2*$c60$+1} $s60$  0.00
+< atom   $c60$      0.00  0.00
+< atom {$c60$+1}    0.00  0.00

#    y
#    ^
#    |              
# s60+  0                0
#    |   \              /
#    |    \            / 
#    |     \          /  
#   0+      0--------0   
#    |
#    +--+---+--------+---+----> x
#       0  0.5      1.5  2

# The distance between first neighbors.
getin d0 1
                
# The primitive cell size
getin a {2*$c60$+2}
getin b {2*$s60$}
getin c 1
            
# To cap the graphene borders with H
# set sp 2
# getin r1 1.2*$r0$
# getin r2 1.3*$r0$
# - bo_range -1 0.9 1.0 $r1$ $r2$
# +< fillh
                                         
        
# Along x-axis, the border exposed is the arm-chair. 
# Along y-axis, the border exposed is the zig-zag. 

