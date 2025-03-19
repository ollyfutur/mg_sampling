#!/bin/bash

n_mg_added=8

echo "Extract PDB motif ..."
pdb_fetch -biounit 2r8s | pdb_selchain -R | pdb_selres -182:188,1001,1002 | pdb_tidy | pdb_chainbows > motif.pdb

# Remove 5' phosphate
sed -i '/MODEL/,+3d' motif.pdb

echo "Making topology ..."
$GMX pdb2gmx -f motif.pdb -ff amber14sb_mMg_tip3p -water tip3p -noh -o conf.gro -p topol.top &>> prepare.log

echo "Solvating ..."
echo "System" | $GMX editconf -f conf.gro -bt dodecahedron -d 1 -o box.gro &>> prepare.log
$GMX solvate -cp box.gro -o box_solv.gro -p topol.top &>> prepare.log

n_water=$(awk '/^SOL/ { print $2 }' topol.top)
echo "Total water: $n_water"

# Add ions
$GMX grompp -f ${PROJECT_DIR}/common/ions.mdp -c box_solv.gro -p topol.top -o ions.tpr &>> prepare.log
echo "SOL" | $GMX genion -s ions.tpr -p topol.top -pname MG -pq 2 -np ${n_mg_added} -nname CL -nq -1 -neutral -o ions.gro &>> prepare.log

rm \#*
