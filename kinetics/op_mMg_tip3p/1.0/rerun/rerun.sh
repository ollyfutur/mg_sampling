#!/bin/bash

$GMX grompp -f $PROJECT_DIR/kinetics/run.mdp -p processed.top -c equil.gro -o rerun_nogroups.tpr -maxwarn 1
cp $PROJECT_DIR/kinetics/run.mdp run_groups.mdp
echo "energygrps = RNA Water" >> run_groups.mdp
$GMX grompp -f run_groups.mdp -p processed.top -c equil.gro -o rerun_groups.tpr -maxwarn 1

$GMX mdrun -s rerun_nogroups.tpr -e rerun_nogroups.edr -notunepme -rerun run.xtc
$GMX mdrun -s rerun_groups.tpr -e rerun_groups.edr -notunepme -rerun run.xtc
