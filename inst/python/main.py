# This is meant to be a first attempt at recreating the work done in Taina's project
# Immonen, Taina T., et al. "Recombination enhances HIV-1 envelope diversity by facilitating
# the survival of latent genomic fragments in the plasma virus population."
# PLoS computational biology 11.12 (2015): e1004625.
# November 12, 2021
# Author: Narmada Sambaturu

# Import python libraries
from pandas import read_csv
from Bio import SeqIO
from time import time
from sys import argv, exit
from os.path import dirname
from os import makedirs
from numpy.random import default_rng
from yaml import safe_load
from random import seed

# Import custom classes
import swish.agents as agents

# Functions to read in data

def read_gen_info(filename):
    gen_info = read_csv(filename, sep='\t', index_col=['generation'])
    #gen_info.columns = ['T']  # To simplify future access. To get count at t=3, cd4_counts.loc[3]['T']
    return gen_info

# Epitopes file has the following columns (tab-delimited)
# col 1 -> epi_start
# col 2 -> epi_end
# col_3 -> min_fitness; this is not used for b-cell epitopes
# col_4 -> max_fitness
def read_b_epitopes(filename):
    epitopes_df = read_csv(filename, sep='\t')
    epitopes = []
    for row in epitopes_df.itertuples():
        epitopes.append(agents.Epitope(int(row[1]), int(row[2]), float(row[4])))
    return epitopes


def get_sequences(filename):
    founder_virus_sequences = [str(fasta.seq).upper() for fasta in  SeqIO.parse(open(filename), 'fasta')]
    # The lengths of the founder sequences must be the same
    len_founder = len(founder_virus_sequences[0])
    assert False not in [len(i) == len_founder for i in founder_virus_sequences], "Founder virus sequences must be of the same length"
    return founder_virus_sequences


def get_conserved_sites(conserved_sites_filename):
    return set(read_csv(conserved_sites_filename, header=None)[0])


def get_nucleotide_substitution_probabilities(substitution_probabilities_filename):
    # Read nucleotide substitution probabilities
    substitution_probabilities_df = read_csv(substitution_probabilities_filename, sep="\t", index_col='sub')

    # Make sure that the row and column labels are the same
    new_nucleotides_order = tuple(substitution_probabilities_df.columns)
    old_nucleotides_order = tuple(substitution_probabilities_df.index)
    assert ''.join(old_nucleotides_order) == ''.join(
        new_nucleotides_order), "Nucleotide substitution matrix row and column labels are different"

    # Get probabilities as a list of tuples
    probabilities = list(substitution_probabilities_df.itertuples(index=False, name=None))

    return new_nucleotides_order, probabilities

# Run model

if __name__ == '__main__':
    if len(argv) != 4 and len(argv) != 5:
        print(len(argv))
        print("argv[1] = config file")
        print("argv[2] = output sequence prefix")
        print("argv[3] = output cell counts prefix")
        print("argv[4] = seed")
        exit(1)

    start = time()

    # Read in config file
    config = safe_load(open(argv[1]))

    # Read in information
    input_files = config['input_files']

    # Read generation information, including cd4 counts and sampling times and numbers
    gen_info = read_gen_info(input_files['gen_info'])

    # Read details of reference sequence
    reference_sequence = get_sequences(input_files['reference_seq'])[0]

    # Epitope start positions and max fitness cost. Epitopes are 30 nucleotides long, contiguous and may overlap.
    founder_virus_epitopes = read_b_epitopes(input_files['epitopes'])

    # Conserved sites are determined separately. They correspond to the positions which don't vary in >=99% of
    # high quality HIV full genomes in the LANL HIV database (irrespective of sub-type).
    conserved_sites = get_conserved_sites(input_files['conserved_sites'])
    
    # Nucleotide substitution probabilities
    nucleotides_order, substitution_probabilities = \
            get_nucleotide_substitution_probabilities(input_files['nt_subst_probs'])

    # Read details about founder virus 
    founder_virus_sequences = get_sequences(input_files['founder_seqs'])

    # Initialize parameters
    params = config['parameters']
    mu = params['mu']
    rho = params['rho']
    eta = params['eta']
    a_L = params['a_L']
    delta = params['delta']
    prolif = params['prolif']
    cost_per_mutation_in_conserved_site = params['psi']
    seroconversion_time = params['seroconversion_time']
    immune_response_proportion = params['f']
    time_to_full_potency = params['d']
    immune_fitness = params['immune_fitness']
    conserved_fitness = params['conserved_fitness']
    replicative_fitness = params['replicative_fitness']
    rf_exp = float(params['rf_exp'])

    # Create founder viruses
    founder_viruses = [agents.HIV(seq, reference_sequence, conserved_sites, replicative_fitness, rf_exp) for seq in founder_virus_sequences]

    # Create host environment and add to viral sequences counter
    host = agents.HostEnv(founder_viruses, int(gen_info.loc[0]['active_cell_count']))

    # Output files
    active_sequences_filename_prefix = argv[2] + 'viral_seqs_active_CD4_'
    latent_sequences_filename_prefix = argv[2] + 'viral_seqs_latent_CD4_'
    counts_filename = argv[3] + 'counts.txt'

    # Make directories if they don't exist
    seq_dirs = dirname(argv[2])
    count_dirs = dirname(argv[3])
    makedirs(seq_dirs, exist_ok=True)
    makedirs(count_dirs, exist_ok=True)

    # Open output file
    out_f = open(counts_filename, 'w')
    out_f.write('\t'.join(['generation', 'active_CD4', 'latent_CD4', 'active_turned_latent', 
        'latent_turned_active', 'latent_died', 'latent_proliferated', 
        'number_mutations', 'number_dual_inf',
        'mean_fitness_active', 'mean_conserved_cost_active', 'mean_immune_cost_active', 'mean_replicative_cost_active']) + '\n')
    # I THINK WE NEED TO MOVE THIS UP BECAUSE WHEN MULTIPLE FOUNDERS THAT'S STOCHASTIC WHEN CREATING HOST ENVIRONMENT
    if len(argv) == 5:
        # Set seed (need to set a seed for random and for numpy.random)
        seed(argv[4])
        generator = default_rng(int(argv[4]))
    else:
        generator = default_rng()

    # Active cell counts for each generation
    active_cell_counts = gen_info['active_cell_count']
    # Number of cells to sample for each generation
    n_to_samp = gen_info['n_sample_active']
    # Last sampled generation (don't have to continue simulation after this)
    last_sampled_gen = [index for index, item in enumerate(n_to_samp) if item != 0][-1]

    if n_to_samp[0] != 0:
        host.record_viral_sequences(active_sequences_filename_prefix + str(0) + ".fasta",
                                        latent_sequences_filename_prefix + str(0) + ".fasta",
                                        n_to_samp = n_to_samp[0])
        mean_fitness_active, mean_conserved_cost_active, mean_immune_cost_active, mean_replicative_cost_active = \
                    host.summarize_fitness()
        out_f.write('\t'.join([str(0), str(len(host.C)), str(len(host.L)), str(0),
                str(0), str(0), str(0), str(0), str(0),
                str(mean_fitness_active), str(mean_conserved_cost_active),
                str(mean_immune_cost_active), str(mean_replicative_cost_active)]) + '\n')

    # Looping through generations until we sample everything we want
    for t in range(1, last_sampled_gen + 1):
        # Latent reservoir dynamics
        num_to_make_latent, num_to_activate, num_to_die, num_to_proliferate = \
                host.get_next_gen_latent(eta, a_L, delta, prolif, generator)
        # Productively infected cell dynamics
        n_mut, number_recombination = \
                host.get_next_gen_active(mu, rho, active_cell_counts[t], t,
                        seroconversion_time, nucleotides_order, substitution_probabilities, conserved_sites,
                        time_to_full_potency, cost_per_mutation_in_conserved_site, reference_sequence, 
                        immune_response_proportion, founder_virus_epitopes, 
                        immune_fitness, conserved_fitness, replicative_fitness, rf_exp, generator)
        # Record events
        if n_to_samp[t] != 0:
            host.record_viral_sequences(active_sequences_filename_prefix + str(t) + ".fasta",
                                        latent_sequences_filename_prefix + str(t) + ".fasta",
                                        n_to_samp = n_to_samp[t])
            mean_fitness_active, mean_conserved_cost_active, mean_immune_cost_active, mean_replicative_cost_active = \
                    host.summarize_fitness()
            out_f.write('\t'.join([str(t), str(len(host.C)), str(len(host.L)), str(num_to_make_latent),
                str(num_to_activate), str(num_to_die), str(num_to_proliferate),
                str(n_mut), str(number_recombination),
                str(mean_fitness_active), str(mean_conserved_cost_active), 
                str(mean_immune_cost_active), str(mean_replicative_cost_active)]) + '\n')

    out_f.close()

    end = time()
    print("Execution time = ", end - start)
