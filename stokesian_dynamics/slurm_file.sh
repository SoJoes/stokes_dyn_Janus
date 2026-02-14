#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1

#SBATCH -p cpu
#SBATCH --qos=debug
#SBATCH --job-name=dbxl46_pytential

#SBATCH -e stderr-file
#SBATCH -o stdout-file

export PYOPENCL_CTX='0'

VENV=/home3/dbxl46/pytential_stokes/pytential_stokes/myenv

echo "Using interpreter:"
$VENV/bin/python --version

$VENV/bin/python run_simulation.py 11 10 1 2 fte