systems="7ez0"
forcefields="nMg_tip3p nMg_tip3p_hard mMg_tip3p mMg_tip3p_hard"
replicates="rep00 rep01 rep02 rep03"

for system in $systems; do
   for ff in $forcefields; do
      dir=$PROJECT_DIR/density/plainmd_${system}_${ff}
      mkdir -p $dir
      cp $system/$ff/rep00/input.gro $dir/input.gro
      for rep in $replicates; do
         $GMX trjconv \
            -f $system/$ff/$rep/run.xtc \
            -s $system/$ff/$rep/run.tpr \
            -b 150000 \
            -e 400000 \
            -dt 1000 \
            -center \
            -pbc mol \
            -o $dir/${system}_${ff}_${rep}.xtc << EOM
RNA
System
EOM
      done
      $GMX trjcat -f $dir/${system}_${ff}_rep*.xtc -cat -o $dir/extracted.xtc
      rm $dir/${system}_${ff}_rep*.xtc
   done
done
