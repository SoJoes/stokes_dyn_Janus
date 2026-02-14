#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1

#SBATCH -p cpu
#SBATCH --qos=debug
#SBATCH --job-name=dbxl46_pytential

#SBATCH -e stderr-file
#SBATCH -o stdout-file

source /etc/profile
module load intel-oneapi

VENV=/home3/dbxl46/pytential_stokes/pytential_stokes/myenv
source $VENV/bin/activate

export PATH="$VENV/bin:$PATH"
export PYOPENCL_CTX='0'

echo "Using interpreter:"
which python
python --version

echo "Installed packages:"
python -m pip list

echo "Python location:"
python -c "import sys; print(sys.executable)"
echo "sys.path:"
python -c "import sys; print(sys.path)"
echo "site-packages contents:"
ls $VENV/lib/python3.8/site-packages | head -20

# Run your script
python run_simulation.py 11 10 1 2 fte