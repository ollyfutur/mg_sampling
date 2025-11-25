# Extract trajectory frames

systems="7ez0_noion"
forcefields="nMg_tip3p"

for system in $systems; do
for forcefield in $forcefields; do
   workdir=plain_md_${system}_${forcefield}
   mkdir -p $workdir
   path=$PROJECT_DIR/plain_md/$system/$forcefield
   cp $path/rep00/input.gro $workdir/input.gro
   i=0
   for rep in $path/rep??; do
      $GMX trjconv -f $rep/run.xtc \
                   -s $rep/run.tpr \
                   -b 150000 \
                   -dt 2000 \
                   -center \
                   -pbc mol \
                   -o $workdir/extracted_${i}.xtc << EOM
RNA
System
EOM
   i=$((i+1))
   done
   $GMX trjcat -f $workdir/extracted_*.xtc -cat -o $workdir/extracted.xtc
done
done
