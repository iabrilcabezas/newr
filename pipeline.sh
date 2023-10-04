#!/bin/bash

module load python
conda activate /global/cfs/cdirs/act/software/iabril/condaenvs/newr
cd /global/common/software/act/iabril/python/newr

# # first: run theory CMB if seed has changed:
# python cl_theory.py

# run foreground:
echo foreground
python ./QEfgs/cl_fgs.py

# run ingredients:
# (comment out cmb lines if already computed)
echo ingredients
python ./QEfgs/ingredients.py

# run plots:
echo plots
python ./QEfgs/plot.py