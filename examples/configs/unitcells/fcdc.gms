# Input file for GeMS

# Define a primitive cell for the generation of a face-centered diamond-cubic
# (fcdc) crystal structure with faces (100).

# fcdc can be viewed as a pair of intersecting face-centered cubic (fcc)
# lattices, with each separated by 1/4 of the width of the unit cell in each
# dimension.

# First primitive fcc cell
>< atom  0      0      0    
+< atom  0.5    0.5    0
+< atom  0    0.5    0.5  
+< atom  0.5    0    0.5  

# 1/4 offset
set move 0.25 0.25 0.25

# Second primitive fcc cell
+< atom  0      0      0    
+< atom  0.5    0.5    0
+< atom  0    0.5    0.5
+< atom  0.5    0    0.5

#     y
#     ^
#     |              
# 0.75+     0.75--------0.25   
#     |       |           |    
#     |       |           |    
#  0.5+ 0.5-- +---0.5     |    
#     |  |    |     |     |    
#     |  |    |     |     |    
# 0.25+  |  0.25----+---0.75   
#     |  |          |         
#     |  |          |         
#    0+  0---------0.5          
#     |                     
#     +--+----+-----+-----+---> x
#        0  0.25   0.5  0.75     

# The distance between first neighbors.
d0:={sqrt(0.25)}
              
# The primitive cell size
a:=1
b:=1
c:=1
          
# All the surface are (100). 
                                      
