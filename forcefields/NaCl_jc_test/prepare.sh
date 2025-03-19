$GMX solvate -box 4 -o box_solv.gro
$GMX pdb2gmx -f box_solv.gro -ff amber14sb_OL15_O8scaled_v2 -water spce -o box_solv.gro

# Add ions
$GMX grompp -f ${PROJECT_DIR}/common/ions.mdp -c box_solv.gro -p topol.top -o ions.tpr &>> prepare.log
echo "SOL" | $GMX genion -s ions.tpr -p topol.top -pname NA -pq 1 -nname CL -nq -1 -conc 0.5 -o ions.gro &>> prepare.log

sed -i "s/^NA   /NA_jc/g" topol.top
sed -i "s/^CL   /CL_jc/g" topol.top


