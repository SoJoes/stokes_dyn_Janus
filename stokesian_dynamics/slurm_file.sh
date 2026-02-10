#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1

#SBATCH -p ug-gpu-small
#SBATCH --qos="debug"
#SBATCH -t 00-00:30:00
#SBATCH --job-name=dbxl46_pytential

# Source the bash profile (required to use the module command)
source /etc/profile

# Run your program (replace this with your program)
source ~/pytential_stokes/pytential_stokes/venv/bin/activate

python run_simulation.py 11 10 1 2 fte