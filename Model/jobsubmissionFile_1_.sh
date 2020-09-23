#!/bin/sh
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=victor.onink@climate.unibe.ch
#SBATCH --job-name=Stochastic_ST=10_RT=69_y2010_ENSEMBLE=1_run=1_restart=0
#SBATCH --output=runOutput/Stochastic_ST=10_RT=69_y2010_ENSEMBLE=1_run=1_restart=0.o%j
#SBATCH --mem-per-cpu=6G
#SBATCH --time=00:10:00
#SBATCH --partition=debug
source /home/ubelix/climate/vo18e689/.bash_profile
source /home/ubelix/climate/vo18e689/anaconda3/bin/activate py3_parcels_v2_2
cd "/home/ubelix/climate/vo18e689/codes/Next-Stage-Plastic-Beaching/Model/"
python generalBeachingScenarios.py -p 10 -v
