systems="7ez0_noion"
forcefields="nMg_tip3p"
replicates="rep00 rep01 rep02 rep03"

for system in $systems; do
   for ff in $forcefields; do
      for rep in $replicates; do
         $GMX trjconv \
            -f $system/$ff/$rep/run.xtc \
            -s $system/$ff/$rep/run.tpr \
            -b 250000 \
            -e 500000 \
            -dt 1000 \
            -center \
            -pbc mol \
            -o extracted/${system}_${ff}_${rep}.xtc << EOM
RNA
System
EOM
      done
      $GMX trjcat -f extracted/${system}_${ff}_*.xtc -cat -o extracted/${system}_${ff}.xtc
   done
done
