#!/bin/bash 

# This file is part of GEMS: Extensible Molecular Simulator
# Copyright (C) 2020 Sergio Alexis Paz
# 
# This code is distributed under the GNU General Public License (GPL) v3 or later.
# See <https://www.gnu.org/licenses/> for the full license text.

set -o nounset
set -o pipefail

# Environment
GEMS="../../usr/bin/gems"

# Log into named pipe
rm -rf test.log
exec > >(tee -i test.log)
exec 2>&1
 
# Timing
TIMEFORMAT="%2Rsec %2Uusr %2Ssys (%P%% cpu)."
start=$(date +%s.%N)
SECONDS=0
 
# Labels
pass='PASS'
fail='FAIL'
  
usage() {
cat <<'WXYZ'
Use:
    test.sh [directory]

Description:
  Execute `directory` test. If directory is not given, execute all tests.

Flags:
  -h this help
  -i interactive mode
  -c highlight with colors
WXYZ
exit
}

exe() {
  # valgrind --leak-check=yes --track-origins=yes ../../src/gems $1 &>> test.log
  local name="$1"; shift

  line="TEST... $name "
  echo -ne "$line \r"

  timing=$( time $GEMS "$name" 2>&1 )
	if [ $? -ne 0 ]; then
    line="ERROR ON EXECUTION"
    echo -ne "$line \r"
    return 1
  fi
	
  line="XXXX : ${timing%.} : $name \r"
  echo -ne "$line \r"
}

check() {
  local mode="$1"; shift
  local tol="$1"; shift
  local files=("$@")
  local ec=0

  for f in "${files[@]}"; do
    case $mode in

      exact)
        diff -q "$f" "ref/$f" || ec=1 ;;

      tol_abs)
        awk -v tol="$tol" '
          function abs(x){return x<0?-x:x}
          (NR==FNR){for(i=1;i<=NF;i++)a[i,NR]=$i;next}
          {for(i=1;i<=NF;i++) if(abs($i-a[i,FNR])>tol) exit 1}
        ' "$f" "ref/$f" || ec=1 ;;
				
      tol_trunc)
        awk -v fmt="$tol" '
          function trunc(x){return sprintf(fmt, x)}
          (NR==FNR){for(i=1;i<=NF;i++)a[i,NR]=$i;next}
          {for(i=1;i<=NF;i++) if(trunc($i)!=trunc(a[i,FNR])) exit 1}
        ' "$f" "ref/$f" || ec=1 ;;
    esac
  done

  [ $ec -eq 0 ] && echo -e "$line$pass" || echo -e "$line$fail"
}

run_test() {
  local name=( $(echo "$1" | tr ',' ' ') )
  local mode="$2"
  local tol="$3"
  local files=( $(echo "$4" | tr ',' ' ') ) 

  # Ejecuto
  for n in "${name[@]}"; do
		 [[ "$n" == "none" ]] && continue
     rm -f "${files[@]}"
     exe "$n"
  done 

  if [[ "$mode" == "lambda" ]]; then
	  $catb -f job0_ch0.bsp -ol 3 2>/dev/null > "${files[0]}"
		check "tol_abs" "$tol" "${files[0]}"
  else
		check "$mode" "$tol" "${files[@]}"
  fi
}
 
mpiexe="mpirun -n 4 $GEMS"
mpiexe=$(echo exit | LJ/$GEMS | grep "^#  MPI: Not compiled for MPI$" > /dev/null && echo no || echo $mpiexe)

# Flags
while getopts ":hic" option; do
case $option in
   h) usage;;
   c) pass='\033[0;32mPASS\033[0m'; fail='\033[0;31mFAIL\033[0m';;
   i) prompt(){ read -p "Test $1? [Y/N] (Y) " ans; [[ $ans =~ ^[nN] ]] && return 1 || return 0; };;
   :) echo "Error: -$OPTARG requires argument"; usage;;
   ?) echo "Error: unknown option -$OPTARG"; usage;;
  esac
done
shift $((OPTIND-1))

# Prompt
if [[ -z ${1:+x} ]]; then
  prompt(){ echo ""; echo "Testing $1"; }
else
  input=${1%/}
  prompt(){ [[ $1 == $input ]]; }
fi

# Test declarations
# Format: [folder]="input:check mode:check parameter:check file"
declare -A tests=(
  [configs]="orthorhombic_crystal.gms:exact:0:gold_cube.xyz \
             graphene.gms:exact:0:graphene_ribon.xyz,graphene_triangle.xyz \
             graphito.gms:exact:0:graphito.xyz \
             tetrahedron_fcc.gms:sdif:1e-7:tetrahedron_fcc.xyz"
  [syntax]="bloques.gms:diff:0:bloques.log"
  [Analitical]="sho.gms:exact:0:Energy.sho.dat \
                wall.gms:exact:0:Energy.wall.dat"
  [Graph]="subgraphs.gms:exact:0:Graph.subgraphs.dat"
  [WTMD]="lj_WTMD_1D.gms:sdif:1e-7:E_libre.dat"
  [SMA-TB]="md.gms:sdif:1e-7:Energy.md.dat \
            bulk.gms:exact:0:Energy.bulk.dat \
            cvs.gms:sdif:1e-7:Energy.cvs.dat \
            cvx.gms:sdif:1e-7:Energy.cvx.dat"
  [optimize]="lbfgs.gms:odif:1d-7:Energy.lbfgs.dat \
              minvol.gms:exact:0:Energy.minvol.dat"
  [LJ]="lj.gms:exact:0:Energy.lj.dat \
        wca.gms:exact:0:Energy.wca.dat \
        slj.gms:exact:0:Energy.slj.dat \
        sm1.gms:exact:0:Energy.sm1.dat \
        feels.gms:exact:0:Energy.feels.dat"
  [NVE]="v_verlet.gms:sdif:1e-5:Energy.v_verlet.dat \
         pred_corr_5.gms:sdif:1e-5:Energy.pred_corr_5.dat"
  [NVT]="ermak.gms:exact:0:Energy.ermak.dat
         ermak_chp.gms:exact:0:Energy.ermak_chp.dat"
  [NPT]="lgf.gms:exact:0:Energy.lgf.dat,Viri.lgf.dat,Press.lgf.dat \
         lgf_flex.gms:exact:0:Energy.lgf_flex.dat,Viri.lgf_flex.dat,Press.lgf_flex.dat \
         lgf_flex_chp.gms:exact:0:Energy.lgf_flex_chp.dat,Viri.lgf_flex_chp.dat,Press.lgf_flex_chp.dat \
         lgf_x.gms:exact:0:Energy.lgf_x.dat,Press.lgf_x.dat,Caja.lgf_x.dat"
  [mVT]="gcmc.gms:sdif:1e-7:Energy.main.dat,Calc.main.dat"
  [DDDA]="md.gms:sdif:1e-7:Eddda.md.dat,Eprom.md.dat,FP_Eddda.md.dat"
  [XYZ]="readxyz.gms:exact:0:Energy.readxyz.dat"
  [hyperdynamics]="compress.gms:sdif:1e-7:Energy.compress.dat \
                   compress_below.gms:sdif:1e-7:Energy.compress_below.dat \
                   lpe.gms:sdif:1e-2:Energy.lpe.dat"
  [NEB]="neb.gms:tdif:%.3e:Energy.neb.dat.NebF"
)


# Test execution
for section in "${!tests[@]}"; do
  if prompt "$section"; then
    cd "$section"
    # separar las entradas por espacios
    for entry in ${tests[$section]}; do
      IFS=":" read -r name mode tol files <<<"$entry"
      # expandir comas en espacios para pasar a run_test
      run_test "$name" "$mode" "$tol" "$files"
    done
    cd ..
  fi
done


# Print total time spent in the test. 
# Despite the Nanosecond format, the precision will be
# around the millisecond or probably less.
end=$(date +%s.%N)
if command -v bc &>/dev/null; then
  echo "Total time: $(echo "$end - $start" | bc -l)"
elif command -v perl &>/dev/null; then
  echo "Total time: $(perl -E "say $end-$start")"
fi

