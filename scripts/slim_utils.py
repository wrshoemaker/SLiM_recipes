import math
import numpy
import os

import config
from itertools import combinations

#file_path = '%stest_out2.txt' % config.data_directory

def get_slim_metadata(file_path):

    metadata_dict = {}

    to_int = ['population_size', 'genome_size', 'number_generations', 'tract_length']

    with open(file_path) as file:
        for line in file:
            line_split = line.strip().split(',')

            if line_split[0] in to_int:

                metadata_dict[line_split[0]] = int(line_split[1])

            else:

                metadata_dict[line_split[0]] = float(line_split[1])

    return metadata_dict



def read_slim_outputFull(file_path):

    populations = False
    mutations = False
    individuals = False
    genomes = False

    #population_sizes_dict = {}
    mutation_dict = {}
    individual_count_dict = {}
    genomes_dict = {}

    with open(file_path) as file:
        for line in file:
            line_split = line.strip().split(' ')

            if 'Populations:' in line_split:
                populations=True
                continue

            if 'Mutations:' in line_split:
                mutations=True
                continue

            if 'Individuals:' in line_split:
                individuals=True
                continue

            if 'Genomes:' in line_split:
                genomes=True
                continue

            if (populations == True) and (mutations == True) and (individuals == False) and (genomes == False):
                # position, fitness, time_origin, frequency
                # start with zero, add counts based on individuals and genomes data
                mutation_dict[line_split[0]] = [int(line_split[3]), float(line_split[4]),  int(line_split[7]), 0]


            if (populations == True) and (mutations == True) and (individuals == True) and (genomes == False):
                population_id = line_split[2].split(':')[0]
                if int(line_split[4]) == 0:
                    continue
                individual_count_dict[line_split[2]] = int(line_split[4])


            if (populations == True) and (mutations == True) and (individuals == True) and (genomes == True):

                # only get genomes that are present in the final sample

                if line_split[0] in individual_count_dict:

                    genome_mutations = [str(mutation) for mutation in line_split[2:] ]

                    for genome_mutation in genome_mutations:

                        mutation_dict[genome_mutation][-1] += individual_count_dict[line_split[0]]

                    genomes_dict[line_split[0]] = set( genome_mutations )


    population_size = sum(individual_count_dict.values())

    mutation_frequency_dict = {}
    individual_frequency_dict = {}

    #mutation_population_dict = {}

    for key in mutation_dict:
        key_copy = mutation_dict[key]
        key_copy[-1] = key_copy[-1]/population_size

        # add empty list, fill with population labels later
        key_copy.append(set())
        mutation_frequency_dict[key] = key_copy

        #mutation_population_dict[key] = set()


    for key, mutation_ids in individual_count_dict.items():
        individual_frequency_dict[key] = individual_count_dict[key]/population_size


    for genome_id, value in genomes_dict.items():

        for site_id in value:

            mutation_frequency_dict[str(site_id)][-1].add(genome_id)

    return mutation_frequency_dict, individual_frequency_dict, genomes_dict, population_size




#def filter_slim_dicts(mutation_frequency_dict,  selection_coefficient=float(0)):





def calculate_unbiased_sigmasquared(mutation_frequency_dict, individual_frequency_dict, genomes_dict, population_size, genome_size, synonymous=True, delta_l=1000):

    # identify pairs of sites within delta_l
    #sliding window of 0.2 log units?

    #genome_size = 100000
    gene_size = 1000

    gene_bins = numpy.logspace(2, numpy.log10(gene_size), num=50, base=10.0)
    gene_bins = numpy.floor(gene_bins)

    mid_gene_bins = numpy.logspace(1, 2, num=15, base=10.0, endpoint=False)
    mid_gene_bins = numpy.floor(mid_gene_bins)

    early_gene_bins = numpy.asarray([1, 6])

    gene_bins = numpy.insert(gene_bins, 0, mid_gene_bins, axis=0)
    gene_bins = numpy.insert(gene_bins, 0, early_gene_bins, axis=0)

    gene_bins_dict = {}

    for gene_bin_idx in range(0, len(gene_bins)-1):

        gene_bins_dict[gene_bins[gene_bin_idx]] = {}
        gene_bins_dict[gene_bins[gene_bin_idx]]['unbiased_sigmasquared_numerator'] = []
        gene_bins_dict[gene_bins[gene_bin_idx]]['unbiased_sigmasquared_denominator'] = []


    for range_i in range(0, genome_size, gene_size):

        gene_positions = range(range_i, range_i+gene_size)

        mutations_to_keep = [key for key, value in mutation_frequency_dict.items() if value[0] in gene_positions]

        # arbitrary minimum number of sites
        if len(mutations_to_keep) < 3:
            continue

        mutation_id_pairs = combinations(mutations_to_keep, 2)

        for mutation_id_pair in mutation_id_pairs:

            mutation_individuals_1 = mutation_frequency_dict[mutation_id_pair[0]][-1]
            mutation_individuals_2 = mutation_frequency_dict[mutation_id_pair[1]][-1]

            mutation_id_pair_distance = abs(mutation_frequency_dict[mutation_id_pair[0]][0] - mutation_frequency_dict[mutation_id_pair[1]][0])

            if mutation_id_pair_distance < early_gene_bins[0]:
                continue

            f_12 = len(mutation_individuals_1.intersection(mutation_individuals_2)) /population_size

            if f_12 ==0:
                continue

            f_1 = mutation_frequency_dict[mutation_id_pair[0]][-2]
            f_2 = mutation_frequency_dict[mutation_id_pair[1]][-2]

            unbiased_sigmasquared_numerator = (f_12 - (f_1 * f_2)) ** 2
            unbiased_sigmasquared_denominator = f_1 * (1-f_1) * f_2 * (1-f_2)

            # dont look at mutation pairs where both mutations occur at the same site
            if mutation_id_pair_distance == 0:
                continue

            sign_change_locations = (numpy.diff(numpy.sign( gene_bins  - mutation_id_pair_distance )) != 0)*1
            bin_location = numpy.where(sign_change_locations==1)[0][0]

            gene_bins_dict[gene_bins[bin_location]]['unbiased_sigmasquared_numerator'].append(unbiased_sigmasquared_numerator)
            gene_bins_dict[gene_bins[bin_location]]['unbiased_sigmasquared_denominator'].append(unbiased_sigmasquared_denominator)


            #if mutation_id_pair_distance == 0:
            #    print( mutation_frequency_dict[mutation_id_pair[0]][0] , mutation_frequency_dict[mutation_id_pair[1]][0] )
            #    print( mutation_id_pair, mutation_id_pair_distance,  f_12)

            #print(mutation_frequency_dict[mutation_pair[0]], mutation_frequency_dict[mutation_pair[1]])

    distances = []
    unbiased_sigmasquared_sums = []


    for key, value in gene_bins_dict.items():

        #if key < 13.0:
        #    print(key, value['unbiased_sigmasquared_denominator'])

        if sum(value['unbiased_sigmasquared_denominator']) == 0:
            continue

        unbiased_sigmasquared = sum(value['unbiased_sigmasquared_numerator']) / sum(value['unbiased_sigmasquared_denominator'])

        distances.append(key)
        unbiased_sigmasquared_sums.append(unbiased_sigmasquared)

    distances = numpy.asarray(distances)
    unbiased_sigmasquared_sums = numpy.asarray(unbiased_sigmasquared_sums)

    return distances, unbiased_sigmasquared_sums
