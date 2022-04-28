#!/bin/bash 
set -o nounset
set -o pipefail

# Set timing information
TIMEFORMAT="%2Rsec %2Uusr %2Ssys (%P%% cpu)."
start=`date +%s.%N`
SECONDS=0

# Redirect stdout ( > ) into a named pipe ( >() ) running "tee"
exec > >(tee -i test.log)

# And error messages too  
exec 2>&1
 
usage()
{
cat <<'WXYZ'

    Use:
    
        test.sh [directory]
    
    Description:

      Execute `directory` test. If directory is not given, execute all tests. 

    Flags:
      -h this help
      -i interactive mode
      -l list aviable tests

WXYZ
exit
}
 
rm -rf test.log
function exe(){
  # valgrind --leak-check=yes --track-origins=yes ../../src/gems $1 &>> test.log
  ../../src/gems $1
}

mpiexe="mpirun -n 4 ../../src/gems"
mpiexe=$(echo exit | ../src/gems | grep "^#  MPI: Not compiled for MPI$" > /dev/null && echo no || echo $mpiexe)

pass='PASS'
fail='FAIL'

# Option parsing
while getopts ":hic" option; do
  case $option in
   h) usage;; 
   c) pass='\033[0;32mPASS\033[0m'
      fail='\033[0;31mFAIL\033[0m'
      ;; 
   i) prompt(){
        read -p "Test $1? [Y/N] (Y)" answer
        case $answer in
         [nN]* ) return 1;;
         [yY]* | * ) return 0;;
        esac
      } ;;
   :)  echo "Error: -$OPTARG requires an argument";usage;exit 1;;
   ?)  echo "Error: unknown option -$OPTARG"      ;usage;exit 1;;
  esac
done
shift $(($OPTIND - 1))
     
# Check for argument  
if [[ -z ${1:+x} ]]; then
  prompt(){
    echo ""
    echo "Testing $1"
  }
else
  input=${1%/}
  prompt(){
    [[ $1 == $input ]]
  }
fi 
           
if prompt configs; then
  cd configs
  rm -f *xyz

  time exe orthorhombic_crystal.gms
  diff  gold_cube.xyz         ref/gold_cube.xyz         && echo -e $pass || echo -e $fail
 
  time exe graphene.gms
  diff  graphene_ribon.xyz    ref/graphene_ribon.xyz    && \
  diff  graphene_triangle.xyz ref/graphene_triangle.xyz && echo -e $pass || echo -e $fail

  time exe graphito.gms
  diff  graphito.xyz          ref/graphito.xyz          && echo -e $pass || echo -e $fail
  
  time exe tetrahedron_fcc.gms
  diff  tetrahedron_fcc.xyz   ref/tetrahedron_fcc.xyz   && echo -e $pass || echo -e $fail

  cd ..
fi
            
if prompt syntax; then
  cd syntax
  rm -f *log

  time exe bloques.gms
  diff <( grep -v time bloques.log | grep -v '^# ' ) \
       <( grep -v time ref/bloques.log | grep -v '^# ' )  && echo -e $pass || echo -e $fail
 
  cd ..
fi
             
if prompt Analitical; then
  cd Analitical
  rm -f Energy.sho.dat
  time exe sho.gms
  paste ref/Energy.sho.dat Energy.sho.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'

  time exe wall.gms
  diff -rq ref/Energy.wall.dat Energy.wall.dat  && echo -e $pass || echo -e $fail

  cd ..
fi         
              
if prompt Graph; then
  cd Graph
  rm -f Graph.subgraphs.dat
  time exe subgraphs.gms
  diff -rq ref/Graph.subgraphs.dat Graph.subgraphs.dat  && echo -e $pass || echo -e $fail
  cd ..
fi         
     

if prompt WTMD; then
  cd WTMD
  rm -f E_libre.dat
  time exe lj_WTMD_1D.gms
  paste ref/E_libre.dat E_libre.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$3)^2;a+=($2-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
  cd ..
fi         
        
if prompt SMA-TB; then
  cd SMA-TB

  echo "... potential "
  rm -f Energy.md.dat
  time exe md.gms
  paste ref/Energy.md.dat Energy.md.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
           
  echo "... potential with linked cells "
  # Me pasa que dependiendo si compilo con OpenMP, aveces el lbfgs me hace un
  # pasito mas y me invalida el test. Si uno mira los numeros que reporta el
  # lbfgs, son ligeramente distintos en las ultimas cifras si uno corre con en
  # openMP con 1 thread que si uno corre serial. 
  # En  https://stackoverflow.com/a/27680902/1342186 explica de forma interesante
  # como la suma de flotantes no es asociativa en programación y que se puede
  # hacer la `Kahan summation` para evitar que esto introduzca diferencias
  # cuando se cambia el numero de threads.
  # Pero aun usando 1 solo thread veo estas diferencias. Entonces debe ser otra
  # cosa.
  # En este thread https://stackoverflow.com/a/27680902/1342186 dice que "The
  # compiler may optimize the code without the OpenMP statements differently." 
  # asique debe de ser eso. Por lo pronto le saco en la comparacion la utlima
  # linea, pero etiqueto con FIXME hasta confirmar que es eso.
  rm -f Energy.bulk.dat
  time exe bulk.gms
  paste ref/Energy.bulk.dat Energy.bulk.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
           
  echo "... CVS "
  rm -f Energy.cvs.dat
  time exe cvs.gms
  paste ref/Energy.cvs.dat Energy.cvs.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}' 

  rm -f Energy.cvx.dat
  time exe cvx.gms
  paste ref/Energy.cvx.dat Energy.cvx.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'  

  cd ..
fi
   
if prompt optimize; then
  cd optimize
  rm -f Energy.lbfgs.dat 
  time exe lbfgs.gms
  paste ref/Energy.lbfgs.dat Energy.lbfgs.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
 
  rm -f Energy.minvol.dat 
  time exe minvol.gms
  paste ref/Energy.minvol.dat Energy.minvol.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$5)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'

  cd ..
fi         
 
if prompt RepExchange; then

  if [[ $mpiexe == no ]]; then
    echo "To run this test compile GeMS for MPI"
  else
    cd RepExchange
         
    rm -f Energy.lj.mpi*.dat
    time $mpiexe lj.gms 2>&1 || echo 'Did you complie GeMS for MPI?'
    cat <(paste ref/Energy.lj.mpi0.dat Energy.lj.mpi0.dat) \
        <(paste ref/Energy.lj.mpi1.dat Energy.lj.mpi1.dat) \
        <(paste ref/Energy.lj.mpi2.dat Energy.lj.mpi2.dat) \
        <(paste ref/Energy.lj.mpi3.dat Energy.lj.mpi3.dat) \
        | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
               {a+=($1-$4)^2;b+=$1}
               END{if(abs(a/b)<1e-6){printf "%-20s %s","partemp...","'$pass'\n"}else{printf "morse... '$fail'\n"}}'
    cd .. 
  fi
fi
             
if prompt LJ; then
  cd LJ
         
  rm -f Energy.lj.dat
  time exe lj.gms
  paste ref/Energy.lj.dat Energy.lj.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "%-20s %s","lj...","'$pass'\n"}else{printf "lj... '$fail'\n"}}'
                
  rm -f Energy.wca.dat
  time exe wca.gms
  paste ref/Energy.wca.dat Energy.wca.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "%-20s %s","wca...","'$pass'\n"}else{printf "wca... '$fail'\n"}}'
                
  rm -f Energy.slj.dat
  time exe slj.gms
  paste ref/Energy.slj.dat Energy.slj.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "%-20s %s","slj...","'$pass'\n"}else{printf "slj... '$fail'\n"}}'
  
  rm -f Energy.sm1.dat
  time exe sm1.gms
  paste ref/Energy.sm1.dat Energy.sm1.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "%-20s %s","sm1...","'$pass'\n"}else{printf "sm1... '$fail'\n"}}' 
 
  rm -f Energy.feels.dat
  time exe feels.gms
  paste ref/Energy.feels.dat Energy.feels.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "%-20s %s","feels...","'$pass'\n"}else{printf "feels... '$fail'\n"}}' 

  cd ..                      
fi
       
if prompt NVE; then
  cd NVE

  echo "... velocity verlet "
  rm -f Energy.v_verlet.dat
  time exe v_verlet.gms
  paste ref/Energy.v_verlet.dat Energy.v_verlet.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
       
  echo "... predictor corrector 5to orden "
  rm -f Energy.pred_corr_5.dat
  time exe pred_corr_5.gms
  paste ref/Energy.pred_corr_5.dat Energy.pred_corr_5.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
 
  cd ..
fi

if prompt NVT; then
  cd NVT

  echo "... ermak "
  rm -f Energy.ermak.dat
  time exe ermak.gms
  paste ref/Energy.ermak.dat Energy.ermak.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$5)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
             
  echo "... checkpoint "
  rm -f Energy.ermak_a.dat Energy.ermak_b.dat
  exe ermak_a.gms
  exe ermak_b.gms
  cat Energy.ermak_a.dat Energy.ermak_b.dat \
    | diff -q - ref/Energy.ermak.dat \
    && echo -e $pass || echo -e $fail
             
  cd ..
fi
   
if prompt NPT; then
  cd NPT
                
  echo "... isotropic G-JF "
  rm -f Energy.lgf.dat Press.lgf.dat  Viri.lgf.dat 
  time exe lgf.gms
  paste ref/Energy.lgf.dat Energy.lgf.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$5)^2;a+=($4-$8)^2;b+=$1}
           END{if(abs(a/b)<1e-6){print "Energy... '$pass'"}else{print "Energy... '$fail'"}}'
  paste ref/Press.lgf.dat Press.lgf.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$10)^2;b+=$1}
           END{if(abs(a/b)<1e-6){print "Pressure... '$pass'"}else{print "Pressure... '$fail'"}}'
  paste ref/Viri.lgf.dat Viri.lgf.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$10)^2;a+=($4-$13)^2;a+=($2-$11)^2;b+=$1}
           END{if(abs(a/b)<1e-6){print "Virial... '$pass'"}else{print "Virial... '$fail'"}}'
                  
  echo "... fully flexible G-JF "
  rm -f Energy.lgf_flex.dat Press.lgf_flex.dat Viri.lgf_flex.dat 
  time exe lgf_flex.gms
  paste ref/Energy.lgf_flex.dat Energy.lgf_flex.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$5)^2;a+=($4-$8)^2;b+=$1}
           END{if(abs(a/b)<1e-6){print "Energy... '$pass'"}else{print "Energy... '$fail'"}}'
  paste ref/Press.lgf_flex.dat Press.lgf_flex.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$10)^2;b+=$1}
           END{if(abs(a/b)<1e-6){print "Pressure... '$pass'"}else{print "Pressure... '$fail'"}}'
  paste ref/Viri.lgf_flex.dat Viri.lgf_flex.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$10)^2;a+=($4-$13)^2;a+=($2-$11)^2;b+=$1}
           END{if(abs(a/b)<1e-6){print "Virial... '$pass'"}else{print "Virial... '$fail'"}}'
  
  echo "... flexible G-JF on x / Ermak on y and z"
  rm -f Energy.lgf_x.dat Press.lgf_x.dat Viri.lgf_x.dat 
  time exe lgf_x.gms
  paste ref/Energy.lgf_x.dat Energy.lgf_x.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$5)^2;a+=($4-$8)^2;b+=$1}
           END{if(abs(a/b)<1e-6){print "Energy... '$pass'"}else{print "Energy... '$fail'"}}'
  paste ref/Press.lgf_x.dat Press.lgf_x.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$10)^2;b+=$1}
           END{if(abs(a/b)<1e-6){print "Pressure... '$pass'"}else{print "Pressure... '$fail'"}}'
  paste ref/Caja.lgf_x.dat Caja.lgf_x.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;a+=($2-$5)^2;a+=($3-$6)^2;b+=$1}
           END{if(abs(a/b)<1e-6){print "Box... '$pass'"}else{print "Box... '$fail'"}}'

          
  echo "... checkpoint "
  rm -f Energy.lgf_flex_a.dat Energy.lgf_flex_b.dat
  exe lgf_flex_a.gms
  exe lgf_flex_b.gms
  cat Energy.lgf_flex_a.dat Energy.lgf_flex_b.dat \
    | diff -q - ref/Energy.lgf_flex.dat \
    && echo -e $pass || echo -e $fail
     
  cd ..
fi
 
if prompt mVT; then
  cd mVT

  rm -f Energy.main.dat
  time exe main.gms

  echo "... gcmc "
  paste ref/Energy.main.dat Energy.main.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$5)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
                 
  echo "... widom "
  paste ref/Calc.main.dat Calc.main.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$2)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
                 
  cd ..
fi
           
if prompt DDDA; then
  cd DDDA

  rm -f Eddda.md.dat  Eprom.md.dat  FP_Eddda.md.dat
  time exe md.gms

  echo "... Promedio "
  paste ref/Eprom.md.dat Eprom.md.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
                             
  echo "... Promedio DDDA"
  paste ref/Eddda.md.dat Eddda.md.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$4)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
                              
  echo "... FP DDDA"
  paste ref/FP_Eddda.md.dat FP_Eddda.md.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$8)^2;a+=($2-$9)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
 
  cd ..
fi
          
if prompt XYZ; then
  cd XYZ

  rm -f Energy.readxyz.dat
  time exe readxyz.gms
  paste ref/Energy.readxyz.dat Energy.readxyz.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$2)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
 
  cd ..
fi
               
if prompt hyperdynamics; then
  cd hyperdynamics
  rm -f Energy.lpe.dat
           
  echo "... compress "
  time exe compress.gms
  paste ref/Energy.compress.dat Energy.compress.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$8)^2;a+=($5-$12)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
            
  echo "... compress_below "
  time exe compress_below.gms
  paste ref/Energy.compress_below.dat Energy.compress_below.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$8)^2;a+=($5-$12)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
             
  echo "... hybrid HDDM "
  time exe lpe.gms
  paste ref/Energy.lpe.dat Energy.lpe.dat \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($1-$7)^2;a+=($5-$11)^2;b+=$1}
           END{if(abs(a/b)<1e-6){printf "'$pass'\n"}else{printf "'$fail'\n"}}'
           
  cd ..
fi
                 
if prompt NEB; then
  cd NEB
  rm -f Energy.neb.dat.NebF
           
  time exe neb.gms
  paste ref/Energy.neb.dat.NebF Energy.neb.dat.NebF \
    | awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
           {a+=($2-$6)^2;a+=($1-$5)^2;b+=$1}
           END{if(abs(a/b)<1e-4){printf "'$pass'\n"}else{printf "'$fail'\n"}}'

  cd ..
fi
     
# Print total time spent in the test. Despite the Nanosecond format, the
# precision will be around the millisedon or probably less.
echo ""
end=`date +%s.%N`
# echo "Total time of the test script $SECONDS"
echo "Total time of the test script $(echo "$end - $start" | bc -l)"

