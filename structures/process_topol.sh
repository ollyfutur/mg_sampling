#!/bin/bash

$GMX grompp -f $PROJECT_DIR/common/run.mdp -p topol.top -pp processed.top -c rep00/equil.gro
