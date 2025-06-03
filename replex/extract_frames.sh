systems="7ez0 7ez0_noion 7ez0_2 7ez0_noion_2"
forcefields="mMg_opc mMg_tip3p nMg_tip3p"

for system in $systems; do
   for ff in $forcefields; do
      $GMX trjconv \
         -f $system/$ff/rep00/run.xtc \
         -s $system/$ff/rep00/run.tpr \
         -b 100000 \
         -dt 100 \
         -center \
         -pbc mol \
         -o extracted/${system}_${ff}.xtc << EOM
RNA
System
EOM
   done
done
