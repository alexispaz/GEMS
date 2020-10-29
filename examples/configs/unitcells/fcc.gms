# Input file for GeMS

# Primitive cell for the generation of a face-centered diamond-cubic
# (fcdc) crystal structure.

# Parameter of the primitive cell:
a0:=1

# Distance between first neighbors:
d0:={sqrt(0.5)}

# Cell size
a:=1
b:=1
c:=1
        
# All the surfaces exposed are (100). 
 

aux:={0.5*$a0}
>< atom  0       0      0    
+< atom  $aux    $aux   0
+< atom  0       $aux   $aux  
+< atom  $aux    0      $aux  
                          
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

