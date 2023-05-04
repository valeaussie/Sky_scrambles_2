#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH -o /fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles/%x.out
#SBATCH --mem=1gb
#SBATCH --tmp=1gb

source ~/.bashrc
export OMP_NUM_THREADS=1
ml conda

conda activate gw

cd /fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles

touch "Running IPTA super scrambles; linear frequencies"

echo "Running IPTA super scrambles; linear frequencies"

mpirun python /fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles/superscramble_val_mod_IPTA.py

cd /fred/oz002/users/mmiles/PTA_misspec/Sky_scrambles

rm -f "Running IPTA super scrambles; linear frequencies"


echo done
