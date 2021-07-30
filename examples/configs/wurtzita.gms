# Input file for GeMS

gr_rad:={$pi$/180}
              
# Celda unidad
# http://www.pveducation.org/pvcdrom/materials/CdSe
# 
# Considering a cell that has 1x1x1 dimensions and an angle of 60 between x and
# y but 90 between x or y and z:
# 
#          (0,1)          (1,1)
#           x--------------x  
#          /              /  A en coordenadas fraccionales 
#         /              /            (s,t)
#        /    A         /      
#       t--- x         /     A en coordenadas cartesianas:
#      /    /         /         (t+s*cos(60),s*sin(60))
#     x----s---------x 
#     (0,0)         (1,0)
# 
 
>< atom    0       0        0
s:=1./3. 
t:=1./3. 
+< atom  {$t$+$s$*cos(60*$gr_rad$)}  {$s$*sin(60*$gr_rad$)}  {1./2.}
set element Cd

>< atom    0       0     {3./8.}
+< atom  {$t$+$s$*cos(60*$gr_rad$)}  {$s$*sin(60*$gr_rad$)}  {7./8.}
set element Se


> creation
out posxyz wurzita_cell.xyz

# Set the number of repetitions of this cells in each direction
# Esto es mejor que sea impar, para que quede una celda centrada en el centro
# que despues podemos usar para elegir donde cortar
nx:= 11
ny:= 11
nz:= 11

# Repeat
> creation
#                 lattice vectors (x y z)            (nr times)       
+< reply           1              0               0      $nx$
+< reply  {cos(60*$gr_rad$)}  {sin(60*$gr_rad$)}  0      $ny$
+< reply           0              0               1      $nz$


# Lattice parameters from:
# Xu, Y.-N., & Ching, W. Y. (1993). Electronic, optical, and structural
# properties of some wurtzite crystals. Physical Review B, 48(7).
# http://doi.org/10.1103/PhysRevB.48.4335
a:=4.2985
c:=7.0152

# Expand the system using the lattice parameters
set expand $a$ $a$ $c$

# Save bulk xyz
out posxyz wurzita_bulk.xyz


# # Perform the cutting operation using planes
# >> below   3.5  3    3     1  2.5  0.5   1  0.5 2.5
# >> below   3.5  3    3     1  2.5  0.5   3  0.5 0.5 
# >> below   3.5  3    3     1  0.5  2.5   3  0.5 0.5
# >> above   1   2.5  0.5    1  0.5  2.5   3  0.5 0.5

# Get a position in the fraction coordinate system of the CdSe cell 1x1x1, for
# instance, the center:
sc:={2./3.}
tc:={2./3.}
zc:={(1./2.+7./8.)/2.}

# Selecciono el numero entero de celdas a añadir
snc:={int(($nx$-1)/2)}
tnc:={int(($ny$-1)/2)}
znc:={int(($nz$-1)/2)}

# Get that position in cartesian coordinates and in the cell that is in the
# center of the bulk
sc:={$a$*($snc$ +$sc$)}
tc:={$a$*($tnc$ +$tc$)}
zc:={$c$*($znc$ +$zc$)}

# Convert it to the cartesian coordinates
xc:={$tc$+$sc$*cos(60*$gr_rad$)} 
yc:={$sc$*sin(60*$gr_rad$)}

# >< atom $xc$ $yc$ $zc$
# set element Xe

# Perform the cutting operation using a sphere of radius 5
> creation
>> sphere 4.0 $xc$ $yc$ $zc$
print cm_pos

# Output the sphere
out posxyz Cd6Se6.xyz  
     


# Perform the cutting operation using a sphere of radius 5
> creation
>> sphere 6.0 $xc$ $yc$ $zc$
print cm_pos

# Output the sphere
out posxyz Cd15Se15.xyz  
     
