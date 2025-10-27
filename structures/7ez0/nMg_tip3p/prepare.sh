#!/bin/bash
from_rep=0
to_rep=15
#####

echo "PDBFixer ..."
pdbfixer --pdbid=7ez0 --add-atoms=heavy --keep-heterogens=all --output=7ez0_fixed.pdb

echo "Making topology ..."
$GMX pdb2gmx -f 7ez0_fixed.pdb -ff amber14sb_nMg_tip3p -water tip3p -noh &>> prepare.log

### We need to ensure PDB ions have the correct atom type
sed -i 's/\(^\s*[0-9]*\s*\)MG /\1nMg/g' topol_Ion_chain_N2.itp

echo "Solvating ..."
echo "System" | $GMX editconf -f conf.gro -princ -d 2 -resnr 1 -o box.gro &>> prepare.log
$GMX solvate -cp box.gro -o box_solv.gro -p topol.top &>> prepare.log

n_water=$(awk '/^SOL/ { print $2 }' topol.top)
echo "Total water: $n_water"

# Negative charge of RNA = number of phosphates
n_rna=$(awk '/^\s*[0-9]+[AUCG]\s+P\s*[0-9]/ { count++ } END { print count }' conf.gro)
echo "N Phosphates (RNA charge): $n_rna"
n_mg_init=$(awk '/^\s*[0-9]+MG/ { count++ } END { print count }' conf.gro)
echo "Initial MG: $n_mg_init"

# Compute number of ions required
# 1) Neutralize RNA with Mg²⁺ and Na⁺ with 3 to 1 ratio
# (paper ion competition: https://doi.org/10.1016%2Fj.bpj.2019.08.007)

# 2*n_mg + n_na = n_rna
# 7*n_na = n_rna
n_mg=$((3*n_rna/7))
echo "Total MG on RNA: $n_mg"
n_na=$((n_rna-2*n_mg))
echo "Total NA on RNA: $n_na"

# 2) Add ions in bulk to target concentrations 10 mM MgCl₂ and 50 mM NaCl
# n_water_bulk ~ 0.86*n_water (estimated from [Cl] in previous iterations)
n_mg_bulk=$(python -c "print(round(0.01/55.5*0.86*$n_water))")
echo "Total MG in bulk: $n_mg_bulk"
n_na_bulk=$((5*n_mg_bulk))
echo "Total NA in bulk: $n_na_bulk"
n_cl=$((2*n_mg_bulk+n_na_bulk))
echo "Total CL in bulk: $n_cl"

# 3) To be added
n_mg_added=$((n_mg+n_mg_bulk-n_mg_init))
n_na_added=$((n_na+n_na_bulk))
n_cl_added=$n_cl
echo "Genion will add: $n_mg_added MG, $n_na_added NA, $n_cl_added CL"

### Prepare N independent replicas
for i in $(seq -w $from_rep $to_rep); do
   echo "Replica $i ..."
   mkdir rep$i
   cp topol.top topol$i.top
   (
      cd rep$i
      $GMX grompp -f ${PROJECT_DIR}/common/ions.mdp -c ../box_solv.gro -p ../topol$i.top -o ions1.tpr &>> prepare.log
      echo "SOL" | $GMX genion -s ions1.tpr -o ions1.gro -p ../topol$i.top -pname MG -pq 2 -np $n_mg_added -nname CL -nq -1 -nn $n_cl_added &>> prepare.log

      $GMX grompp -f ${PROJECT_DIR}/common/ions.mdp -c ions1.gro -p ../topol$i.top -o ions2.tpr &>> prepare.log
      echo "SOL" | $GMX genion -s ions2.tpr -o ions2.gro -p ../topol$i.top -pname NA -pq 1 -np $n_na_added &>> prepare.log
   )
done

printf -v i0 "%02d" $from_rep
[ -f topol${i0}.top ] && mv topol${i0}.top topol.top

for i in $(seq -w $((from_rep+1)) $to_rep); do
   diff topol.top topol$i.top && rm topol$i.top || echo "Error: Topologies should be identical…"
done

rm \#*
