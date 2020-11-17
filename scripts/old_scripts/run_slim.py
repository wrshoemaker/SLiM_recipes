import os
import argparse
import config

import slim_utils


slim_script = 'test_recombination.slim'

directory = config.directory

parser = argparse.ArgumentParser()

parser.add_argument('--generations', required=True, type=int)
parser.add_argument('--genome_size', required=True, type=int)
parser.add_argument('--population_size', required=True, type=int)

parser.add_argument('--per_base_mutation_rate', required=False, type=float, default=1e-8)
parser.add_argument('--per_base_recombination_rate', required=False, type=float, default=1e-7)

parser.add_argument('--beneficial_selection_coefficient', required=False, type=float, default=0.001)
parser.add_argument('--deleterious_selection_coefficient', required=False, type=float, default=-0.001)

parser.add_argument('--proportion_beneficial_sites', required=False, type=float, default=0)
parser.add_argument('--proportion_deleterious_sites', required=False, type=float, default=0)
parser.add_argument('--proportion_neutral_sites', required=False, type=float, default=1)


args = parser.parse_args()

generations = args.generations
population_size = args.population_size
genome_size = args.genome_size
per_base_u = args.per_base_mutation_rate
per_base_r = args.per_base_recombination_rate

beneficial_s = args.beneficial_selection_coefficient
deleterious_s = args.deleterious_selection_coefficient

proportion_b = args.proportion_beneficial_sites
proportion_d = args.proportion_deleterious_sites
proportion_n = args.proportion_neutral_sites

# extremely annoying, but to use % syntax with slim you have to have the outfile  in double parentheses and then single parentheses
file_out = "'%stest_out3.txt'" % directory

#slim_line = 'slim -d generations=%d -d population_size=%d -d genome_size=%d -d per_base_u=%f -d per_base_r=%f -d beneficial_s=%f -d deleterious_s=%f -d proportion_b=%f -d proportion_d=%f -d proportion_n=%f -d "file_out=%s" %stest_recombination.slim' % (generations, population_size, genome_size, per_base_u, per_base_r, beneficial_s, deleterious_s, proportion_b, proportion_d, proportion_n, file_out, directory)
#slim_line = 'slim -d generations=%d -d population_size=%d -d genome_size=%d -d per_base_u=%f -d per_base_r=%f -d beneficial_s=%f -d deleterious_s=%f -d proportion_b=%f -d proportion_d=%f -d proportion_n=%f %stest_recombination.slim' % (generations, population_size, genome_size, per_base_u, per_base_r, beneficial_s, deleterious_s, proportion_b, proportion_d, proportion_n, directory)

slim_line = 'slim -d generations={} -d population_size={} -d genome_size={} -d  per_base_u={} -d per_base_r={} -d beneficial_s={} -d deleterious_s={} -d proportion_b={} -d proportion_d={} -d proportion_n={} -d file_out="{}"  {}test_recombination.slim'.format(generations, population_size, genome_size, per_base_u, per_base_r, beneficial_s, deleterious_s, proportion_b, proportion_d, proportion_n, file_out, directory)



#print(slim_line)

#slim_utils.

os.system(slim_line)

#for s in [args.s]:
#    f = '%stest_file.txt' % directory
    #os.system('slim -d L=10000 -d u={}  -d r={} -d N={} -d s_b={} -d "fname=\'{}\'" -d nbases={} ~/project-ngarud/code/slim/{}.slim'.
    #          format(mu, r, ne, s, f, nbases, args.f))

    #os.system('slim -d L=%d '  )

    #          '%s%s' % (directory, slim_script)




# mutations info to parse output

#0: within-file numeric identifier
#1: Mutation id property, a within-run unique identifier for mutations that does not change over time
#2: mutation type
#3: chromosome position, zero-based
#4: selection coefficient
#5: dominance coefficient
#6: subpopulation ID
#7 generation where it originated
#8: mutation prevalence=count of the number of times that it occurs in any genome in the population.


# Individuals section of file
#0: individual ID
#1: sex
#2: Genome specifiers, p1:0 (indicating the 0th genome in p1).
#3: second genome, since slim makes us simulate diploids.
#4: number of individuals

# genomes section of file
#0: genome specifier, such as p1:0
#1: type of genome: an A for an autosome
#2: list of within-file mutation identifiers,as given in the Mutations section
