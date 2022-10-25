#!/bin/bash
#SBATCH --job-name=compile_Enzo   # Job name
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=g.kerex@gmail.com     # Where to send mail	
#SBATCH --time=48:00:00               # Time limit hrs:min:sec
#SBATCH -p cmbas -c 10

pwd; hostname; date

./make-rusty.sh

date

