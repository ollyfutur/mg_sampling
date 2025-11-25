systems="7ez0_noion"
forcefields="nMg_tip3p"
replicates="rep00 rep01 rep02 rep03"

for system in $systems; do
   for ff in $forcefields; do
      for rep in $replicates; do
         mkdir -p mddb/$system/$ff/$rep
         echo "non-Water" | $GMX convert-tpr \
            -s ../$system/$ff/$rep/run.tpr \
            -n ../$system/$ff/index.ndx \
            -o $system/$ff/$rep/solute.tpr
         echo "System" | $GMX editconf \
            -f $system/$ff/$rep/solute.tpr \
            -o $system/$ff/$rep/solute.pdb
         echo "non-Water" | $GMX trjconv \
            -f ../$system/$ff/$rep/run.xtc \
            -s ../$system/$ff/$rep/run.tpr \
            -n ../$system/$ff/index.ndx \
            -o $system/$ff/$rep/solute.xtc
      done
   done
done
