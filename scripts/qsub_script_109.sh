#!/bin/bash
#$ -N run_analyses_109
#$ -e /Users/wrshoemaker/GitHub/slim_recipes/data/negative_selection/qsub_script_error_109
#$ -o /Users/wrshoemaker/GitHub/slim_recipes/data/negative_selection/qsub_script_output_109
#$ -l h_data=32G
#$ -l time=84:00:00
#$ -l highp
#$ -m bea
. /u/local/Modules/default/init/modules.sh

module unload python
module load python/3.6.1
module load anaconda
source activate slim

python3 /Users/wrshoemaker/GitHub/slim_recipes/scripts/run_nonWF_gene_conversion.py --iterations 2 --simulation_number 108 --slim_data_directory /Users/wrshoemaker/GitHub/slim_recipes/data/negative_selection/
