#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH --gres=gpu

#SBATCH -p ug-gpu-small
#SBATCH --qos=debug
#SBATCH --job-name=dbxl46_pytential

#SBATCH -e stderr-file
#SBATCH -o stdout-file

module avail

source /etc/profile
module load intel-oneapi/vtune
module load intel-oneapi/mpi
module load intel-oneapi/compiler
module list

VENV=/home3/dbxl46/pytential_stokes/pytential_stokes/myenv
source $VENV/bin/activate

export PATH="$VENV/bin:$PATH"

python -c "import pyopencl as cl; print(cl.get_platforms())"

export PYOPENCL_CTX='0'
ls $VENV/lib/python3.8/site-packages | head -20

# Run your script
python run_simulation.py 11 10 1 2 fte