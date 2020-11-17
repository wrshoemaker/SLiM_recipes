#!/bin/bash
#$ -N run_analyses
#$ -e /u/home/w/wrshoema/project-ngarud/slim_recipes/slim_error
#$ -o /u/home/w/wrshoema/project-ngarud/slim_recipes/slim_output
#$ -l h_data=8G
#$ -l time=24:00:00
#$ -l highp
#$ -m bea
#$ -t 1:30
. /u/local/Modules/default/init/modules.sh

module unload python
module load python/3.6.1

module load anaconda

source activate slim
#conda create -n slim python=3.6
python3 /u/home/w/wrshoema/project-ngarud/slim_recipes/run_slim.py --generations 200000 --genome_size 1000000 --population_size 100000
