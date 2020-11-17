
import config

import slim_utils
import matplotlib.pyplot as plt

import scipy.stats

import numpy


def get_null(N, probability_hgt, tract_length, genome_size, distance_max=3000):

    distances = numpy.arange(1,distance_max+1)

    expected_sigmasquared = []

    for distance in distances:

        R = probability_hgt * (tract_length/genome_size) * scipy.stats.geom.cdf(distance, 1/tract_length)

        expected_sigmasquared_i = (10+(2*N*R)) / (22+ 26*N*R + 4*((N*R)**2) )

        expected_sigmasquared.append(expected_sigmasquared_i)

    expected_sigmasquared = numpy.asarray(expected_sigmasquared)


    return distances, expected_sigmasquared




def plot_null_ld():

    fig = plt.figure(figsize = (4, 4)) #

    row_idx = 0

    for number_idx, number in enumerate(list(range(200,209)) ):

        column_idx = number_idx %3

        if (column_idx == 0) and (number_idx >0):
            row_idx += 1

        data_filename = "%stest_HGT/test_out_nonWF_%d.txt" % ( config.data_directory, number)

        metadata_filename = "%stest_HGT/simulation_parameters_%d.csv" % ( config.data_directory, number)

        metadata_dict = slim_utils.get_slim_metadata(metadata_filename)

        probability_hgt = metadata_dict['genome_size'] * metadata_dict['per_base_recombination_rate']
        tract_length = metadata_dict['tract_length']

        mutation_frequency_dict, individual_frequency_dict, genomes_dict, population_size, = slim_utils.read_slim_outputFull(data_filename)

        print(mutation_frequency_dict)

        distances, unbiased_sigmasquared_sums = slim_utils.calculate_unbiased_sigmasquared(mutation_frequency_dict, individual_frequency_dict, genomes_dict, population_size, metadata_dict['genome_size'])

        ax_i = plt.subplot2grid((3, 3), (row_idx, column_idx))

        ax_i.scatter(distances, unbiased_sigmasquared_sums, alpha=0.5)

        distances, expected_sigmasquared = get_null(metadata_dict['population_size'], probability_hgt, tract_length, metadata_dict['genome_size'], distance_max=3000)

        ax_i.plot(distances, expected_sigmasquared, 'k--')

        ax_i.set_yscale("log")
        ax_i.set_xscale("log")
        ax_i.set_xlim(0.7, 3300)
        ax_i.set_ylim(0.01, 1.2)
        ax_i.tick_params(axis='x', labelsize=6)
        ax_i.tick_params(axis='y', labelsize=6)

        #ax_i.set_xlabel('Distance between SNVs, $\ell$', fontsize=8)
        #ax_i.set_ylabel('Linkage disequilibrium, $\sigma^2_d$', fontsize=8)

        ax_i.text(0.45,0.25, r'$P(HGT) = {{{}}}$'.format(str( round(probability_hgt, 3) )), fontsize=6, color='k', ha='center', va='center', transform=ax_i.transAxes  )
        ax_i.text(0.35,0.1, '$\ell_{r} = %d$kbp' % int(tract_length/1000), fontsize=6, color='k', ha='center', va='center', transform=ax_i.transAxes  )

    fig.text(0.5, 0.02, 'Distance between SNVs, $\ell$', va='center', ha='center', fontsize=14)
    fig.text(-0.02, 0.5, 'Linkage disequilibrium, $\sigma^2_d$', va='center', rotation='vertical', fontsize=14)

    fig.subplots_adjust(wspace=0.3)

    fig.savefig('%stest.png' % (config.analysis_directory) , bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
    plt.close()



def plot_ld_selection():

    fig = plt.figure(figsize = (4, 4)) #

    row_idx = 0

    for number_idx, number in enumerate(list(range(0,9)) ):

        column_idx = number_idx %3

        if (column_idx == 0) and (number_idx >0):
            row_idx += 1

        data_filename = "%stest_selection_HGT/test_out_nonWF_%d.txt" % ( config.data_directory, number)

        metadata_filename = "%stest_selection_HGT/simulation_parameters_%d.csv" % ( config.data_directory, number)

        metadata_dict = slim_utils.get_slim_metadata(metadata_filename)

        probability_hgt = metadata_dict['genome_size'] * metadata_dict['per_base_recombination_rate']
        s_deleterious = metadata_dict['s_deleterious']
        NS = s_deleterious * metadata_dict['population_size']


        mutation_frequency_dict, individual_frequency_dict, genomes_dict, population_size, = slim_utils.read_slim_outputFull(data_filename)

        mutation_frequency_dict_negative = {key: value for key, value in mutation_frequency_dict.items() if value[1] < 0}
        mutation_frequency_dict_neutral = {key: value for key, value in mutation_frequency_dict.items() if value[1] == 0}

        distances_neutral, unbiased_sigmasquared_sums_neutral = slim_utils.calculate_unbiased_sigmasquared(mutation_frequency_dict_neutral, individual_frequency_dict, genomes_dict, population_size, metadata_dict['genome_size'])
        distances_negative, unbiased_sigmasquared_sums_negative = slim_utils.calculate_unbiased_sigmasquared(mutation_frequency_dict_negative, individual_frequency_dict, genomes_dict, population_size, metadata_dict['genome_size'])


        ax_i = plt.subplot2grid((3, 3), (row_idx, column_idx))

        ax_i.scatter(distances_neutral, unbiased_sigmasquared_sums_neutral, alpha=0.5, c='b', label='Neutral', edgecolors='k', s=20)
        ax_i.scatter(distances_negative, unbiased_sigmasquared_sums_negative, alpha=0.5, c='r', label='$Ns=%0.2f$' % NS, edgecolors='k', s=20)


        ax_i.set_yscale("log")
        ax_i.set_xscale("log")
        ax_i.set_xlim(0.7, 3300)
        ax_i.set_ylim(0.00001, 1.2)
        ax_i.tick_params(axis='x', labelsize=6)
        ax_i.tick_params(axis='y', labelsize=6)

        ax_i.set_title(r'$P(HGT) = {{{}}}$'.format(str( round(probability_hgt, 3) )), fontsize=6 )

        #ax_i.text(0.45,0.85, r'$P(HGT) = {{{}}}$'.format(str( round(probability_hgt, 3) )), fontsize=6, color='k', ha='center', va='center', transform=ax_i.transAxes  )

        ax_i.legend(loc="lower left", fontsize=4)




    fig.text(0.5, 0.02, 'Distance between SNVs, $\ell$', va='center', ha='center', fontsize=14)
    fig.text(-0.02, 0.5, 'Linkage disequilibrium, $\sigma^2_d$', va='center', rotation='vertical', fontsize=14)

    fig.subplots_adjust(wspace=0.4, hspace=0.4)

    fig.savefig('%stest_selection.png' % (config.analysis_directory) , bbox_inches = "tight", pad_inches = 0.2, dpi = 600)
    plt.close()
