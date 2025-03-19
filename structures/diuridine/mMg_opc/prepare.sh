$GMX pdb2gmx -f ../diuridine.pdb -ff amber14sb_mMg_opc -water tip4p

# We don’t need this
rm topol_Ion2.itp
rm posre_Ion2.itp

# Update water topology
sed -i 's/tip4p.itp/opc.itp/g' topol.top

# Rename this
sed -i 's/^Ion2/MG/g' topol.top

# Discard this
sed -i '/Ion2/d' topol.top

# Renumber residues
$GMX editconf -f conf.gro -resnr 1 -o conf.gro
rm \#conf.gro.1#

# Rename HO4 to SOL for consistency
sed -i 's/HO4/SOL/g' conf.gro

# Minimize and process single topology file
# warning: non-zero charge
$GMX grompp -f ../min.mdp -c conf.gro -p topol.top -pp processed.top -o min.tpr -maxwarn 1
$GMX mdrun -deffnm min

$GMX grompp -f ../equil.mdp -c conf.gro -p processed.top -o equil.tpr -maxwarn 1
$GMX mdrun -deffnm equil
