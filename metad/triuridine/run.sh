export GMX_ENABLE_DIRECT_GPU_COMM=true

if [ ! -f run.tpr ]; then
  $GMX grompp -f $PROJECT_DIR/common/run_metad.mdp -p processed.top -c equil.gro -o run.tpr -maxwarn 1
fi

$GMX mdrun -deffnm run -cpi run.cpt -plumed plumed.dat
