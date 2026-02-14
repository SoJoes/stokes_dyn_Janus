#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1

#SBATCH -p cpu
#SBATCH --qos=debug
#SBATCH --job-name=dbxl46_pytential

#SBATCH -e stderr-file
#SBATCH -o stdout-file

export PYOPENCL_CTX='0'
source $VENV/bin/activate

echo "Using interpreter:"
which python
python --version

# Run your script
python run_simulation.py 11 10 1 2 fte