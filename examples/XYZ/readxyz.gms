(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))
(                    INPUT SYSTEM FOR GeMS PROGRAM                 )
(                   comentaries are between brackets               )
(                      caps can be lower or upper                  )
(               caps must be upper for potential names             )
(((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))


dimension 3

box size    13.1399  12.64396  20.0000 

>< read Pos.md.xyz
sys add
box mic f

> sys
group 1 add

> sys
set pbc t t t

interact 1 under reax

outfile 1 name E.$jobname$.dat
outfile 1 cols e_reax
outfile 1 each 1
#                 
# outfile 2 name Pos.$jobname$.xyz
# outfile 2 Pos  1
# outfile 2 each 1

outfile 3 name Energy.$jobname$.dat
outfile 3 cols epot 1
outfile 3 each 1
         
# outfile 4 name Ch.$jobname$.xyz
# outfile 4 charge 1
# outfile 4 each 1
#        
# outfile 5 name Caja.$jobname$.dat
# outfile 5 cols box
# outfile 5 each 1
     
> sys

dinamica_from_xyz 6 Pos.md.xyz
