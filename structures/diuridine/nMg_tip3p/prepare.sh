$GMX pdb2gmx -f ../diuridine.pdb -ff amber14sb_nMg_tip3p -water tip3p

# We don’t need this
rm topol_Ion2.itp
rm posre_Ion2.itp

# Rename this
sed -i 's/^Ion2/MG/g' topol.top

# Discard this
sed -i '/Ion2/d' topol.top

# Renumber residues
$GMX editconf -f conf.gro -resnr 1 -o conf.gro
rm \#conf.gro.1#

# Rename HOH to SOL for consistency
sed -i 's/HOH/SOL/g' conf.gro

# Minimize and process single topology file
# warning: non-zero charge
$GMX grompp -f ../min.mdp -c conf.gro -p topol.top -pp processed.top -o min.tpr -maxwarn 1
$GMX mdrun -deffnm min

$GMX grompp -f ../equil.mdp -c conf.gro -p processed.top -o equil.tpr -maxwarn 1
$GMX mdrun -deffnm equil
