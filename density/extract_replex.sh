# Extract trajectory frames

systems="7ez0 7ez0_noion"
forcefields="mMg_tip3p"
#replicates="0"
replicates="1"

for system in $systems; do
for forcefield in $forcefields; do
for replicate in $replicates; do
   workdir=replex_${system}_${forcefield}_${replicate}
   mkdir -p $workdir
   path=$PROJECT_DIR/replex/$system/$forcefield/replicate$replicate
   if [ $replicate -eq 0 ]; then rep="00"; else rep="08"; fi
   cp $path/rep$rep/input.gro $workdir/input.gro
   $GMX trjconv -f $path/rep$rep/run.xtc \
                -s $path/rep$rep/run.tpr \
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
