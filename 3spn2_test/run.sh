#!/bin/bash

#SBATCH --job-name=3spn2_test
#SBATCH --output=3spn2_test.out
#SBATCH --error=3spn2_test.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pi-depablo
#SBATCH --partition=depablo-gpu
#SBATCH --gres=gpu:1

module load boost/1.55+python-2.7-2014q1
module load gcc/4.7
module load cuda/8.0

python 3spn2_test_charges.py
