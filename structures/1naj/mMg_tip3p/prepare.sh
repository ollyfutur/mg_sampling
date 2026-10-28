#!/bin/bash
from_rep=0
to_rep=2
#####

echo "PDBFixer ..."
pdbfixer --pdbid=1naj --add-atoms=heavy --keep-heterogens=none --output=1naj_fixed.pdb

echo "Making topology ..."
$GMX pdb2gmx -f 1naj_fixed.pdb -ff amber14sb_mMg_tip3p -water tip3p -noh &>> prepare.log

echo "Solvating ..."
echo "System" | $GMX editconf -f conf.gro -princ -bt cubic -d 2 -o box.gro &>> prepare.log
$GMX solvate -cp box.gro -o box_solv.gro -p topol.top &>> prepare.log

n_water=$(awk '/^SOL/ { print $2 }' topol.top)
echo "Total water: $n_water"

# Negative charge of RNA = number of phosphates
n_rna=$(awk '/^\s*[0-9]+[AUCG]\s+P\s*[0-9]/ { count++ } END { print count }' conf.gro)
echo "N Phosphates (RNA charge): $n_rna"

n_mg_added=11

# 3) To be added
echo "Genion will add: $n_mg_added MG"

### Prepare N independent replicas
for i in $(seq -w $from_rep $to_rep); do
   printf -v i "%02d" $i
   echo "Replica $i ..."
   mkdir rep$i
   cp topol.top topol$i.top
   (
      cd rep$i
      $GMX grompp -f ${PROJECT_DIR}/common/ions.mdp -c ../box_solv.gro -p ../topol$i.top -o ions2.tpr &>> prepare.log
      echo "SOL" | $GMX genion -s ions2.tpr -o ions2.gro -p ../topol$i.top -pname MG -pq 2 -np $n_mg_added &>> prepare.log
   )
done

printf -v i0 "%02d" $from_rep
[ -f topol${i0}.top ] && mv topol${i0}.top topol.top

for i in $(seq -w $((from_rep+1)) $to_rep); do
   printf -v i "%02d" $i
   diff topol.top topol$i.top && rm topol$i.top || echo "Error: Topologies should be identical…"
done

rm \#*
