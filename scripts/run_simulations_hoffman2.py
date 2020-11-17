import os
import argparse
import config
import sys
import math

import slim_utils

import itertools

data_directory = config.data_directory
scripts_directory = config.scripts_directory

# make the output directory

slim_output_name = "negative_selection"

simulation_data_directory = '%s%s/' % (data_directory, slim_output_name)

os.system('mkdir -p %s' % simulation_data_directory)


slim_path = "%snonWF_ConstantSize_Selection.slim" % scripts_directory

rescale_factor = 10
# rescale total mutation input (N/Q) * (Q*mu) * genome size
# rescale total recombination events (N/Q) * (per-base prob * Q) * (genome size / Q) * (Q * recombination fragment/genome size)


iterations = 20

population_size = int(1e6) / rescale_factor

genome_size = 1e6


s_beneficial = 0.01

#s_deleterious_list = [-1e-6, -1e-5, -1e-4]
s_deleterious_list = [-1e-6]

proportion_beneficial = 0
proportion_deleterious = 0.5
proportion_neutral = 1 - (proportion_beneficial + proportion_deleterious)

per_base_mutation_rate = 1e-8 * rescale_factor

per_base_recombination_rate_list = [1.53e-8]
#per_base_recombination_rate_list = [1e-8, 1e-9, 1e-10]

#tract_lengths = [122000.0]
tract_length_list = [100000]
#tract_length_list = [100000.0, 10000.0, 1000.0]

number_generations = 2*population_size


# just a count of number of scripts for reference
simulation_counts = 100

for s_deleterious in s_deleterious_list:

    s_deleterious = s_deleterious * rescale_factor

    for per_base_recombination_rate in per_base_recombination_rate_list:

        per_base_recombination_rate = per_base_recombination_rate * rescale_factor

        for tract_length in tract_length_list:

            # put parameters in another file
            simulation_count_filename = '%ssimulation_parameters_%d.csv' % (simulation_data_directory, simulation_counts)
            simulation_count_file = open(simulation_count_filename, 'w')

            simulation_count_file.write(",".join(['population_size', str(int(population_size))]) + '\n' )
            simulation_count_file.write(",".join(['genome_size', str(int(genome_size))]) + '\n' )
            simulation_count_file.write(",".join(['number_generations', str(int(number_generations))]) + '\n' )

            simulation_count_file.write(",".join(['s_beneficial', str(s_beneficial)]) + '\n' )
            simulation_count_file.write(",".join(['s_deleterious', str(s_deleterious)]) + '\n' )

            simulation_count_file.write(",".join(['proportion_neutral', str(proportion_neutral)]) + '\n' )
            simulation_count_file.write(",".join(['proportion_beneficial', str(proportion_beneficial)]) + '\n' )
            simulation_count_file.write(",".join(['proportion_deleterious', str(proportion_deleterious)]) + '\n' )

            simulation_count_file.write(",".join(['per_base_mutation_rate', str(per_base_mutation_rate)]) + '\n' )
            simulation_count_file.write(",".join(['per_base_recombination_rate', str(per_base_recombination_rate)]) + '\n' )
            simulation_count_file.write(",".join(['tract_length', str(int(tract_length))]) + '\n' )
            simulation_count_file.write(",".join(['rescale_factor', str(rescale_factor)]) + '\n' )

            simulation_count_file.close()


            #outputFull_path = "'%stest_out_nonWF_%d.txt'" % (simulation_data_directory, simulation_counts)


            #slim_cmd_list = ['slim',
            #            "-t", # print SLiM's total execution time (in user clock time)
            #            "-m", # print SLiM's peak memory usage
            #            #"-d", """"runId='{}'" """.format(os.path.join(out_dir, run_id)),
            #            #"-d", """"runIdShort='{}'" """.format(run_id),
            #            "-d", """"N_generations={}" """.format(int(number_generations)),
            #            "-d", """"genomeSize={}" """.format(int(genome_size)),
            #            "-d", """"Ne={}" """.format(int(population_size)),

            #            "-d", """"Mu={}" """.format(per_base_mutation_rate),
            #            "-d", """"Rho={}" """.format(per_base_recombination_rate),
            #            "-d", """"tractlen={}" """.format(int(tract_length)),

            #            "-d", """"s_beneficial={}" """.format(float(s_beneficial)),
            #            "-d", """"s_deleterious={}" """.format(float(s_deleterious)),

            #            "-d", """"proportion_neutral={}" """.format(float(proportion_neutral)),
            #            "-d", """"proportion_beneficial={}" """.format(float(proportion_beneficial)),
            #            "-d", """"proportion_deleterious={}" """.format(float(proportion_deleterious)),

            #            "-d", """"outputFull_path={}" """.format(outputFull_path),

            #            #"-d", """"sampleSize={}" """.format(int(sample_size)),
            #            #"-d", """"gcBurnin={}" """.format(gcBurnin),
            #            slim_path]


            #slim_cmd = " ".join(slim_cmd_list)

            python_cmd_list = ['python3',
                            "%srun_nonWF_gene_conversion.py" % (scripts_directory),
                            "--iterations %d" % iterations,
                            "--simulation_number %d" % simulation_counts,
                            "--slim_data_directory %s" % simulation_data_directory ]


            python_cmd = " ".join(python_cmd_list)



            #if os.geteuid() == 501:

            #    os.system(slim_cmd)


            # creates shell script

            run_name = 'run_analyses_%s_%d' % (slim_output_name, simulation_counts)

            error_path = '%sqsub_script_%s_%d_error' % (simulation_data_directory, slim_output_name, simulation_counts)
            output_path = '%sqsub_script_%s_%d_output' % (simulation_data_directory, slim_output_name, simulation_counts)
            shell_script_filename = "%sqsub_script_%s_%d.sh" % (scripts_directory, slim_output_name, simulation_counts)

            for fileee in [error_path, output_path, shell_script_filename]:

                try:
                    os.remove(fileee)
                except OSError:
                    pass


            shell_script_file = open(shell_script_filename, 'w')

            shell_script_file.write('#!/bin/bash' + '\n')
            shell_script_file.write(('#$ -N %s' % run_name) + '\n')
            shell_script_file.write(('#$ -e %s' % error_path) + '\n')
            shell_script_file.write(('#$ -o %s' % output_path) + '\n')
            shell_script_file.write('#$ -l h_data=32G' + '\n')
            shell_script_file.write('#$ -l time=36:00:00' + '\n')
            #shell_script_file.write('#$ -l highp ngarud' + '\n')
            shell_script_file.write('#$ -m bea' + '\n')
            #shell_script_file.write('#$ -t 1:30' + '\n')
            shell_script_file.write('. /u/local/Modules/default/init/modules.sh' + '\n')

            shell_script_file.write('\n')

            shell_script_file.write('module unload python' + '\n')
            shell_script_file.write('module load python/3.6.1' + '\n')
            shell_script_file.write('module load anaconda' + '\n')
            shell_script_file.write('source activate slim' + '\n')

            shell_script_file.write('\n')
            #shell_script_file.write(slim_cmd + '\n')
            shell_script_file.write(python_cmd + '\n')

            shell_script_file.close()

            os.system('qsub %s' % shell_script_filename)


            # cnsi cnsi_msa.q (highp)
            #mem/node=32G slots/node=8  slots=80


            simulation_counts += 1
