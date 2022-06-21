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
      -c higlight with colors

WXYZ
exit
}
 
rm -rf test.log
function exe(){
  # valgrind --leak-check=yes --track-origins=yes ../../src/gems $1 &>> test.log
	line="TEST... $1 "
	echo -ne "$line \r"

	timing=$( { time for i in "${@:2}"; do ../../src/gems $i 2>&1; done; } 2>&1)	

	line="XXXX : ${timing%.} : $1 \r"
	echo -ne "$line \r"
}
 
pass='PASS'
fail='FAIL'
          
function odif(){
# Compare last value within a certain tolerance

	ec=""
	for i in ${@:2}; do
		awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
	    (NR==FNR){for(i=1;i<=NF;i++){a[i,NR]=$i};next} 
			END{for(i=1;i<=NF;i++){if(abs($i-a[i,FNR])>'$1'){exit 1}}}' $i ref/$i
		(( ec = ec || $? ))
  done

	[ $ec -eq 0 ]	&& echo -e $line$pass || echo -e $line$fail
}
         

		
     
function tdif(){
# Compare within a certain tolerance

	ec=""
	for i in ${@:2}; do
    # https://stackoverflow.com/a/61406365
		awk 'function trunc(x){return sprintf("'$1'", x)}
	    (NR==FNR){for(i=1;i<=NF;i++){a[i,NR]=$i};next} 
			{for(i=1;i<=NF;i++){if(trunc($i)!=trunc(a[i,FNR])){exit 1}}}' $i ref/$i
		(( ec = ec || $? ))
  done

	[ $ec -eq 0 ]	&& echo -e $line$pass || echo -e $line$fail
}
          
     
function sdif(){
# Compare within a certain tolerance

	ec=""
	for i in ${@:2}; do
		awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
	    (NR==FNR){for(i=1;i<=NF;i++){a[i,NR]=$i};next} 
			{for(i=1;i<=NF;i++){if(abs($i-a[i,FNR])>'$1'){exit 1}}}' $i ref/$i
		(( ec = ec || $? ))
  done

	[ $ec -eq 0 ]	&& echo -e $line$pass || echo -e $line$fail
}
          

function dif(){
  # valgrind --leak-check=yes --track-origins=yes ../../src/gems $1 &>> test.log

	ec=""
	for i in $@; do
  		diff -q $i ref/$i; (( ec = ec || $? ))
  done

# awk 'function abs(x){return (((x < 0.0) ? -x : x) + 0.0)}
#      {a+=($1-$4)^2;b+=$1}
#      END{if(abs(a/b)<1e-6){printf "%-20s %s","partemp...","'$pass'\n"}else{printf "morse... '$fail'\n"}}'
	[ $ec -eq 0 ]	&& echo -e $line$pass || echo -e $line$fail
}
 
mpiexe="mpirun -n 4 ../../src/gems"
mpiexe=$(echo exit | ../src/gems | grep "^#  MPI: Not compiled for MPI$" > /dev/null && echo no || echo $mpiexe)

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

  exe orthorhombic_crystal orthorhombic_crystal.gms
  dif gold_cube.xyz
 
  exe graphene graphene.gms
  dif  graphene_ribon.xyz    
  dif  graphene_triangle.xyz 

  exe graphito graphito.gms
  dif  graphito.xyz       
  
  exe tetrahedron_fcc tetrahedron_fcc.gms
  sdif '1e-7' tetrahedron_fcc.xyz

  cd ..
fi
            
if prompt syntax; then
  cd syntax
  rm -f *log

  exe bloques bloques.gms
  diff <( grep -v time bloques.log | grep -v '^# ' ) \
       <( grep -v time ref/bloques.log | grep -v '^# ' )  && echo -e $pass || echo -e $fail
 
  cd ..
fi
             
if prompt Analitical; then
  cd Analitical

  rm -f Energy.sho.dat
  exe sho sho.gms
  dif Energy.sho.dat

  rm -f Energy.wall.dat
  exe wall wall.gms
  dif Energy.wall.dat

  cd ..
fi         
              
if prompt Graph; then
  cd Graph
  rm -f Graph.subgraphs.dat
  exe subgraphs subgraphs.gms
  dif Graph.subgraphs.dat
  cd ..
fi         
     

if prompt WTMD; then
  cd WTMD
  rm -f E_libre.dat
  exe lj_WTMD_1D lj_WTMD_1D.gms
  sdif '1e-7' E_libre.dat
  cd ..
fi         
        
if prompt SMA-TB; then
  cd SMA-TB

  rm -f Energy.md.dat
  exe md md.gms
  sdif '1e-7' Energy.md.dat 
           
  # echo "... potential with linked cells "
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
  exe bulk bulk.gms
  dif Energy.bulk.dat
           
  rm -f Energy.cvs.dat
  exe cvs cvs.gms
  sdif '1e-7' Energy.cvs.dat

  rm -f Energy.cvx.dat
  exe cvx cvx.gms
  sdif '1e-7' Energy.cvx.dat

  cd ..
fi
   
if prompt optimize; then
  cd optimize
  rm -f Energy.lbfgs.dat 
  exe lbfgs lbfgs.gms
  odif '1d-7' Energy.lbfgs.dat
 
  rm -f Energy.minvol.dat 
  exe minvol minvol.gms
  dif Energy.minvol.dat

  cd ..
fi         
 
if prompt RepExchange; then

  if [[ $mpiexe == no ]]; then
    echo "To run this test compile GeMS for MPI"
  else
    cd RepExchange
         
    rm -f Energy.lj.mpi*.dat
    time $mpiexe lj.gms 2>&1 || echo 'Did you complie GeMS for MPI?'
		diff Energy.lj.mpi*.dat
    cd .. 
  fi
fi
             
if prompt LJ; then
  cd LJ
    
	list="lj wca slj sm1 feels"
	for i in $list; do
		rm -f Energy.$i.dat
		exe $i $i.gms
		dif Energy.$i.dat
	done

  cd ..                      
fi
       
if prompt NVE; then
  cd NVE
    
	list="v_verlet pred_corr_5"
	for i in $list; do
		rm -f Energy.$i.dat
		exe $i $i.gms
		sdif '1e-5' Energy.$i.dat
	done

  cd ..
fi

if prompt NVT; then
  cd NVT

  rm -f Energy.ermak.dat
  exe ermak ermak.gms
  dif Energy.ermak.dat
             
  rm -f Energy.ermak_a.dat Energy.ermak_b.dat Energy.ermak.dat
  exe "checkpoint" ermak_a ermak_a.gms ermak_b ermak_b.gms
  cat Energy.ermak_a.dat Energy.ermak_b.dat > Energy.ermak.dat
  dif Energy.ermak.dat 
             
  cd ..
fi
   
if prompt NPT; then
  cd NPT
                
	list="lgf lgf_flex"
	for i in $list; do
		rm -f Energy.$i.dat
		exe $i $i.gms
		dif Energy.$i.dat Viri.$i.dat Press.$i.dat
	done
     
	rm -f Energy.lgf_x.dat
	exe lgf_x lgf_x.gms
	dif Energy.lgf_x.dat Press.lgf_x.dat Caja.lgf_x.dat

  rm -f Energy.lgf_flex_a.dat Energy.lgf_flex_b.dat Energy.lgf_flex.dat
  exe "checkpoint" lgf_flex_a.gms lgf_flex_b.gms
  cat Energy.lgf_flex_a.dat Energy.lgf_flex_b.dat > Energy.lgf_flex.dat
  dif Energy.lgf_flex.dat
     
  cd ..
fi
 
if prompt mVT; then
  cd mVT

  rm -f Energy.main.dat
  exe "gcmc" main.gms
	sdif '1e-7' Energy.main.dat Calc.main.dat 

  cd ..
fi
           
if prompt DDDA; then
  cd DDDA

  rm -f Eddda.md.dat Eprom.md.dat FP_Eddda.md.dat
  exe md md.gms
	sdif '1e-7' Eddda.md.dat Eprom.md.dat FP_Eddda.md.dat
 
  cd ..
fi
          
if prompt XYZ; then
  cd XYZ

  rm -f Energy.readxyz.dat
  exe readxyz readxyz.gms
  dif Energy.readxyz.dat
 
  cd ..
fi
               
if prompt hyperdynamics; then
  cd hyperdynamics
           
	list="compress compress_below"
	for i in $list; do
		rm -f Energy.$i.dat
		exe $i $i.gms
		sdif '1e-7' Energy.$i.dat
	done
             
	rm -f Energy.lpe.dat
	exe lpe lpe.gms
	sdif '1e-2' Energy.lpe.dat

  cd ..
fi
                 
if prompt NEB; then
  cd NEB
           
	rm -f Energy.neb.dat
	exe neb neb.gms
	tdif '%.3e' Energy.neb.dat.NebF

  cd ..
fi
     
# Print total time spent in the test. Despite the Nanosecond format, the
# precision will be around the millisedon or probably less.
echo ""
end=`date +%s.%N`
# echo "Total time of the test script $SECONDS"
command -V bc &> /dev/null && { echo "Total time of the test script $(echo "$end - $start" | bc -l)"; exit; }
command -V perl &> /dev/null && { echo "Total time of the test script $(perl -E "say $end-$start")"; exit; }

