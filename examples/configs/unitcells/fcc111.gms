# Input file for GeMS

# Periodic cell for the generation of a
# face-centered diamond-cubic (fcdc) crystal
# structure.

# Parameter of the primitive cell:
a0:=1

# Distance between first neighbors:
d0:={sqrt(0.5)}
 
# z-surface exposed is (111). 
       
# Cell size
a:=$d0
b:={2*$d0*sin($pi/3)}
c:={9*$d0*tan($pi/6)/2}
     
# Creacion de celda fcc con bordes 111
x:={$c/9}
# y={s60/sqrt(2)}
>+ atom  0.0   0.0  0.0
^+ atom  {$d0/2}  {$b/2}  0.
^+ atom  {$d0/2}     $x  {$b/2}
^+ atom  0.0      {4*$x} {$b/2}
^+ atom  0.0      {2*$x}   $b
^+ atom  {$d0/2}  {5*$x}   $b

