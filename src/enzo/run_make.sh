#!/bin/bash
#SBATCH --job-name=compile_Enzo   # Job name
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=g.kerex@gmail.com     # Where to send mail	
#SBATCH --time=00:30:00               # Time limit hrs:min:sec
#SBATCH -p genx -c 2
###SBATCH -p gpu --gpus=1 -c 1  -C a100
###SBATCH -p gpu --gpus=1 -c 1  -C v100

pwd; hostname; date

module add cuda
moudle add cudnn

./make-rusty.sh

date

