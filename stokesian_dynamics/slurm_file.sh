#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1

#SBATCH -p cpu
#SBATCH --qos=debug
#SBATCH --job-name=dbxl46_pytential

#SBATCH -e stderr-file
#SBATCH -o stdout-file

export PYOPENCL_CTX='0'

# Source the bash profile (required to use the module command)
source /etc/profile
source home3/dbxl46/pytential_stokes/pytential_stokes/myenv/bin/activate

python --version
python run_simulation.py 11 10 1 2 fte