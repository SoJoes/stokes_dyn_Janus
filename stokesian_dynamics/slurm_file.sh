#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1

#SBATCH -p cpu
#SBATCH --qos=debug
#SBATCH --job-name=dbxl46_pytential

#SBATCH -e stderr-file
#SBATCH -o stdout-file

export PYOPENCL_CTX='0'

source /etc/profile

echo "Before activation:"
which python
python --version

source /home3/dbxl46/pytential_stokes/pytential_stokes/myenv/bin/activate

echo "After activation:"
which python
python --version
echo "VIRTUAL_ENV=$VIRTUAL_ENV"

/home3/dbxl46/pytential_stokes/pytential_stokes/myenv/bin/python --version

/home3/dbxl46/pytential_stokes/pytential_stokes/myenv/bin/python \run_simulation.py 11 10 1 2 fte