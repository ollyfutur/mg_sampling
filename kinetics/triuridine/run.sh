if [ ! -f run.tpr ]; then
  $GMX grompp -f $PROJECT_DIR/kinetics/triuridine/run.mdp -p ../processed.top -c ../equil.gro -o run.tpr -maxwarn 1
fi

$GMX mdrun -deffnm run -cpi run.cpt -plumed plumed.dat -maxh 23.9
