#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1

#SBATCH -p ug-gpu-small
#SBATCH --qos=debug
#SBATCH --job-name=dbxl46_pytential
#SBATCH --gres=gpu

export PYOPENCL_CTX='0'

# Source the bash profile (required to use the module command)
source /etc/profile

# Run your program (replace this with your program)
source ~/pytential_stokes/pytential_stokes/venv/bin/activate

python run_simulation.py 11 10 1 2 fte