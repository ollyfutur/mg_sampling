#!/bin/bash
from_rep=0
to_rep=7
#####

echo "PDBFixer ..."
# pdbfixer --pdbid=7ez0 --add-atoms=heavy --keep-heterogens=none --output=7ez0_fixed.pdb

echo "Making topology ..."
$GMX pdb2gmx -f 7ez0_fixed.pdb -ff amber14sb -water tip3p -noh &>> prepare.log

echo "Solvating ..."
echo "System" | $GMX editconf -f conf.gro -princ -d 2 -resnr 1 -o box.gro &>> prepare.log
$GMX solvate -cp box.gro -o box_solv.gro -p topol.top &>> prepare.log

n_water=$(awk '/^SOL/ { print $2 }' topol.top)
echo "Total water: $n_water"

### Prepare N independent replicas
for i in $(seq -w $from_rep $to_rep); do
   i=$(printf "%02d" $i)
   echo "Replica $i ..."
   mkdir rep$i
   cp box_solv.gro rep$i/ions2.gro
done
