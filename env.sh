#!/bin/bash

export PROJECT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
export GMX=${GMX:-$(command -v gmx || command -v gmx_mpi)}
[[ -n $GMX ]] && echo "Gromacs executable: $GMX" || echo "Gromacs executable not found."
export FF_DIR=$PROJECT_DIR/forcefields/
export GMXLIB=$FF_DIR

export PATH=$PATH:$PROJECT_DIR/scripts/

[[ -d $PROJECT_DIR/.venv ]] && source $PROJECT_DIR/.venv/bin/activate
