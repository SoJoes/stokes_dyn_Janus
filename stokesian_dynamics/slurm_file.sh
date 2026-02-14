#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1

#SBATCH -p cpu
#SBATCH --qos=debug
#SBATCH --job-name=dbxl46_pytential

#SBATCH -e stderr-file
#SBATCH -o stdout-file

VENV=/home3/dbxl46/pytential_stokes/pytential_stokes/myenv
source $VENV/bin/activate

export PATH="$VENV/bin:$PATH"
export PYOPENCL_CTX='0'

echo "Using interpreter:"
which python
python --version

echo "Installed packages:"
python -m pip list

# Run your script
python run_simulation.py 11 10 1 2 fte