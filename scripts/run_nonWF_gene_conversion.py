import os
import argparse
import config
import sys
import math
import gzip

import slim_utils

import itertools

data_directory = config.data_directory
scripts_directory = config.scripts_directory

# make the output directory
intermediate_filename_template = '%s%s.txt.gz'


parser = argparse.ArgumentParser()
#parser.add_argument("--slim_cmd", type=str, help="full slim script")
parser.add_argument("--iterations", type=int, help="number of iterations of the simulation", default=10)
parser.add_argument("--simulation_number", type=int, help="reference number for SLiM output")
parser.add_argument("--slim_data_directory", type=str, help="directory for slim run")

args = parser.parse_args()

#slim_cmd = args.slim_cmd
iterations = args.iterations
simulation_number = args.simulation_number
slim_data_directory = args.slim_data_directory



#simulation_data_directory = '%stest_selection_HGT/' % data_directory


metadata_filename = "%s/simulation_parameters_%d.csv" % ( slim_data_directory, simulation_number)
metadata_dict = slim_utils.get_slim_metadata(metadata_filename)

outputFull_path =  "'%stest_out_nonWF_%d.txt'" % (slim_data_directory, simulation_number)



slim_path = "%snonWF_ConstantSize_Selection.slim" % scripts_directory

slim_cmd_list = ['slim',
            "-t", # print SLiM's total execution time (in user clock time)
            "-m", # print SLiM's peak memory usage
            #"-d", """"runId='{}'" """.format(os.path.join(out_dir, run_id)),
            #"-d", """"runIdShort='{}'" """.format(run_id),
            "-d", """"N_generations={}" """.format(int(metadata_dict['number_generations'])),
            "-d", """"genomeSize={}" """.format(int(metadata_dict['genome_size'])),
            "-d", """"Ne={}" """.format(int(metadata_dict['population_size'])),

            "-d", """"Mu={}" """.format(metadata_dict['per_base_mutation_rate']),
            "-d", """"Rho={}" """.format(metadata_dict['per_base_recombination_rate']),
            "-d", """"tractlen={}" """.format(int(metadata_dict['tract_length'])),

            "-d", """"s_beneficial={}" """.format(float(metadata_dict['s_beneficial'])),
            "-d", """"s_deleterious={}" """.format(float(metadata_dict['s_deleterious'])),

            "-d", """"proportion_neutral={}" """.format(float(metadata_dict['proportion_neutral'])),
            "-d", """"proportion_beneficial={}" """.format(float(metadata_dict['proportion_beneficial'])),
            "-d", """"proportion_deleterious={}" """.format(float(metadata_dict['proportion_deleterious'])),

            "-d", """"outputFull_path={}" """.format(outputFull_path),

            #"-d", """"sampleSize={}" """.format(int(sample_size)),
            #"-d", """"gcBurnin={}" """.format(gcBurnin),
            slim_path]



slim_cmd = " ".join(slim_cmd_list)


distances_neutral_all = []
distances_negative_all = []

unbiased_sigmasquared_sums_neutral_all = []
unbiased_sigmasquared_sums_negative_all = []

record_strs = []

for iteration in range(iterations):

    os.system(slim_cmd)


    mutation_frequency_dict, individual_frequency_dict, genomes_dict, population_size, = slim_utils.read_slim_outputFull(outputFull_path.strip("\'"))

    mutation_frequency_dict_negative = {key: value for key, value in mutation_frequency_dict.items() if value[1] < 0}
    mutation_frequency_dict_neutral = {key: value for key, value in mutation_frequency_dict.items() if value[1] == 0}

    distances_neutral, unbiased_sigmasquared_sums_neutral = slim_utils.calculate_unbiased_sigmasquared(mutation_frequency_dict_neutral, individual_frequency_dict, genomes_dict, population_size, metadata_dict['genome_size'])
    distances_negative, unbiased_sigmasquared_sums_negative = slim_utils.calculate_unbiased_sigmasquared(mutation_frequency_dict_negative, individual_frequency_dict, genomes_dict, population_size, metadata_dict['genome_size'])

    distances_neutral_str =  ":".join([str(x) for x in distances_neutral])
    unbiased_sigmasquared_sums_neutral_str =  ":".join([str(x) for x in unbiased_sigmasquared_sums_neutral])

    distances_negative_str =  ":".join([str(x) for x in distances_negative])
    unbiased_sigmasquared_sums_negative_str =  ":".join([str(x) for x in unbiased_sigmasquared_sums_negative])

    record_strs.append( ", ".join([str(iteration), distances_neutral_str, unbiased_sigmasquared_sums_neutral_str, distances_negative_str, unbiased_sigmasquared_sums_negative_str]) )



sys.stderr.write("Writing intermediate LD file...\n")
intermediate_filename = intermediate_filename_template % (slim_data_directory, 'ld_%d' % simulation_number)
file = gzip.open(intermediate_filename,"w")
record_str = "\n".join(record_strs)
record_str = record_str.encode(encoding='UTF-8')
file.write(record_str)
file.close()
sys.stderr.write("Done!\n")
