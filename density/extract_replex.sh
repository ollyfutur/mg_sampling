# Extract trajectory frames

systems="7ez0 7ez0_noion"
forcefields="mMg_tip3p"
replicates="0"

for system in $systems; do
for forcefield in $forcefields; do
for replicate in $replicates; do
   workdir=replex_${system}_${forcefield}_${replicate}
   mkdir -p $workdir
   path=$PROJECT_DIR/replex/$system/$forcefield/replicate$replicate
   cp $path/rep00/input.gro $workdir/input.gro
   $GMX trjconv -f $path/rep00/run.xtc \
                -s $path/rep00/run.tpr \
                -b 50000 \
                -dt 100 \
                -center \
                -pbc mol \
                -o $workdir/extracted.xtc << EOM
RNA
System
EOM
done
done
done
