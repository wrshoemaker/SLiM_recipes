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

simulation_data_directory = '%stest_selection_HGT/' % data_directory

os.system('mkdir -p %s' % simulation_data_directory)

#parser = argparse.ArgumentParser()

#parser.add_argument('--N_generations', required=False, type=int, default=int(1e3))
#parser.add_argument('--genomeSize', required=False, type=int, default=int(1e5))
#parser.add_argument('--Ne', required=False, type=int, default=int(1e4))

#parser.add_argument('--Mu', required=False, type=float, default=1e-8)
#parser.add_argument('--Rho', required=False, type=float, default=1e-7)
#parser.add_argument('--tractlen', required=False, type=int, default=500)

#parser.add_argument('--s_beneficial', required=False, type=float, default=0.001)
#parser.add_argument('--s_deleterious', required=False, type=float, default=-0.001)

#parser.add_argument('--proportion_neutral', required=False, type=float, default=1)
#parser.add_argument('--proportion_beneficial', required=False, type=float, default=0)
#parser.add_argument('--proportion_deleterious', required=False, type=float, default=0)


#args = parser.parse_args()

#N_generations = args.N_generations
#genomeSize = args.genomeSize
#Ne = args.Ne

#Mu = args.Mu
#Rho = args.Rho
#tractlen = args.tractlen

#s_beneficial = args.s_beneficial
#s_deleterious = args.s_deleterious

#proportion_neutral = args.proportion_neutral
#proportion_beneficial = args.proportion_beneficial
#proportion_deleterious = args.proportion_deleterious





slim_path = "%snonWF_ConstantSize_Selection.slim" % scripts_directory


rescale_factor = 10
# rescale total mutation input (N/Q) * (Q*mu) * genome size
# rescale total recombination events (N/Q) * (per-base prob * Q) * (genome size / Q) * (Q * recombination fragment/genome size)


iterations = 10

population_size = int(1e6) / rescale_factor

genome_size = 1e6

s_beneficial = 0.01

s_deleterious_list = [-1e-6, -1e-5, -1e-4]


proportion_beneficial = 0.5
proportion_deleterious = 0.5
proportion_neutral = 1 - (proportion_beneficial + proportion_deleterious)

per_base_mutation_rate = 1e-9 * rescale_factor

#per_base_recombination_rates = [1.53e-8]
per_base_recombination_rate_list = [1e-8, 1e-9, 1e-10]

#tract_lengths = [122000.0]
tract_length_list = [100000]
#tract_length_list = [100000.0, 10000.0, 1000.0]

number_generations = 2*population_size


# just a count of number of scripts for reference
simulation_counts = 0


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

            simulation_counts += 1


            slim_cmd_list = ['slim',
                        "-t", # print SLiM's total execution time (in user clock time)
                        "-m", # print SLiM's peak memory usage
                        #"-d", """"runId='{}'" """.format(os.path.join(out_dir, run_id)),
                        #"-d", """"runIdShort='{}'" """.format(run_id),
                        "-d", """"N_generations={}" """.format(int(number_generations)),
                        "-d", """"genomeSize={}" """.format(int(genome_size)),
                        "-d", """"Ne={}" """.format(int(population_size)),

                        "-d", """"Mu={}" """.format(per_base_mutation_rate),
                        "-d", """"Rho={}" """.format(per_base_recombination_rate),
                        "-d", """"tractlen={}" """.format(int(tract_length)),

                        "-d", """"s_beneficial={}" """.format(float(s_beneficial)),
                        "-d", """"s_deleterious={}" """.format(float(s_deleterious)),

                        "-d", """"proportion_neutral={}" """.format(float(proportion_neutral)),
                        "-d", """"proportion_beneficial={}" """.format(float(proportion_beneficial)),
                        "-d", """"proportion_deleterious={}" """.format(float(proportion_deleterious)),

                        "-d", """"outputFull_path={}" """.format(outputFull_path),

                        #"-d", """"sampleSize={}" """.format(int(sample_size)),
                        #"-d", """"gcBurnin={}" """.format(gcBurnin),
                        slim_path]


            slim_cmd = " ".join(slim_cmd_list)


            for iteration in range(iterations):

                outputFull_path = "'%stest_out_nonWF_%d_%d.txt'" % (simulation_data_directory, simulation_counts, iteration)


                if os.geteuid() == 501:

                    os.system(slim_cmd)


                    mutation_frequency_dict, individual_frequency_dict, genomes_dict, population_size, = slim_utils.read_slim_outputFull(data_filename)

                    mutation_frequency_dict_negative = {key: value for key, value in mutation_frequency_dict.items() if value[1] < 0}
                    mutation_frequency_dict_neutral = {key: value for key, value in mutation_frequency_dict.items() if value[1] == 0}


                else:

                    # creates shell script

                    run_name = 'run_analyses_%d' % simulation_counts

                    error_path = '%sqsub_script_error_%d' % (simulation_data_directory, simulation_counts)
                    output_path = '%sqsub_script_output_%d' % (simulation_data_directory, simulation_counts)
                    shell_script_filename = "%sqsub_script_%d.sh" % (scripts_directory, simulation_counts)

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
                    shell_script_file.write('#$ -l time=84:00:00' + '\n')
                    shell_script_file.write('#$ -l highp' + '\n')
                    shell_script_file.write('#$ -m bea' + '\n')
                    #shell_script_file.write('#$ -t 1:30' + '\n')
                    shell_script_file.write('. /u/local/Modules/default/init/modules.sh' + '\n')

                    shell_script_file.write('\n')

                    shell_script_file.write('module unload python' + '\n')
                    shell_script_file.write('module load python/3.6.1' + '\n')
                    shell_script_file.write('module load anaconda' + '\n')
                    shell_script_file.write('source activate slim' + '\n')

                    shell_script_file.write('\n')
                    shell_script_file.write(slim_cmd + '\n')

                    shell_script_file.close()

                    #os.system('qsub %s' % shell_script_filename)







#slim_line = 'slim -d N_generations={} -d genomeSize={} -d Ne={} -d  Mu={} -d Rho={} -d tractlen={} -d s_beneficial={} -d s_deleterious={} -d proportion_neutral={} -d proportion_beneficial={} -d proportion_deleterious={} -d outputFull_path="{}" {}'.format(N_generations, genomeSize, Ne, Mu, Rho, tractlen, s_beneficial, s_deleterious, proportion_neutral, proportion_beneficial, proportion_deleterious, outputFull_path, slim_path)


#sys.stderr.write("Starting simulation...\n")


# write the shell script for hoffman

#tshell_lines = [, ,  '#$ -o /u/home/w/wrshoema/project-ngarud/slim_recipes/slim_output']
