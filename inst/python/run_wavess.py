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
from yaml import safe_load
from csv import writer
from math import exp
import numpy as np

# Import custom classes and functions
import agents

# Functions to read in data


def read_pop_samp(filename):
    pop_samp = read_csv(filename, index_col=["generation"])
    return pop_samp


# Epitopes file has the following columns (tab-delimited)
# col 1 -> epi_start
# col 2 -> epi_end
# col_3 -> max_fitness
def read_b_epitopes(filename):
    epitopes_df = read_csv(filename)
    epitopes = []
    for row in epitopes_df.itertuples():
        epitopes.append(agents.Epitope(
            int(row[1]), int(row[2]), float(row[3])))
    return epitopes


def get_sequences(filename):
    founder_virus_sequences = [
        str(fasta.seq).upper() for fasta in SeqIO.parse(open(filename), "fasta")
    ]
    # The lengths of the founder sequences must be the same
    len_founder = len(founder_virus_sequences[0])
    assert False not in [
        len(i) == len_founder for i in founder_virus_sequences
    ], "Founder virus sequences must be of the same length"
    return founder_virus_sequences


def get_conserved_sites(conserved_sites_filename):
    return read_csv(conserved_sites_filename,header=0).set_index("position")["nucleotide"].to_dict()


def get_nucleotide_substitution_probabilities(q_filename, mut_rate):
    # Read nucleotide substitution probabilities
    q = read_csv(q_filename, index_col="nt_from")
    nucleotides_order, substitution_probabilities = agents.calc_nt_sub_probs_from_q(q, mut_rate)
    return nucleotides_order, substitution_probabilities
    


# Run model

if __name__ == "__main__":
    if len(argv) != 3 and len(argv) != 4:
        print(len(argv))
        print("argv[1] = config file")
        print("argv[2] = output prefix")
        print("argv[3] = seed")
        exit(1)

    start = time()
    
    # error if no slash in output file
    if "/" not in argv[2]:
      exit('Please provide a directory name that contains a / (e.g. wavess_output/)')

    s = None
    if len(argv) == 4:
        s = argv[3]
    generator = agents.set_python_seed(s)

    # Read in config file
    config = safe_load(open(argv[1]))

    ### Read in information ###
    input_files = config["input_files"]
    params = config["parameters"]

    ## Required inputs ##

    # Generation information, including cd4 counts and sampling times and
    # numbers
    pop_samp = read_pop_samp(input_files["pop_samp"])

    # Founder virus
    founder_virus_sequences = get_sequences(input_files["founder_seqs"])
    assert (
        len(founder_virus_sequences) == pop_samp["active_cell_count"][0]
    ), "The number of founder sequences must equal the active cell count at generation 0."
    # Create founder viruses
    founder_viruses = {}
    for i, v in enumerate(founder_virus_sequences):
        founder_viruses["founder" + str(i)] = v

    # Nucleotide substitution probabilities
    nucleotides_order, substitution_probabilities = (
        get_nucleotide_substitution_probabilities(
            input_files["q"], params["mut_rate"])
    )

    ## Optional inputs related to selection ##

    # Conserved sites
    conserved_sites = {}
    if input_files["conserved_sites"] != "":
        conserved_sites = get_conserved_sites(input_files["conserved_sites"])

    # Reference sequence
    reference_sequence = ""
    if input_files["ref_seq"] != "":
        reference_sequence = get_sequences(input_files["ref_seq"])[0]
        if len(conserved_sites):
          # remove conserved sites that are different between the reference and any founder,
          # and mask any conserved sites with a - in the reference
          reference_sequence, conserved_sites = agents.prep_ref_conserved(founder_viruses, reference_sequence, conserved_sites)
          
    # Epitope start positions and max fitness cost.
    epitope_locations = None
    if input_files["epitope_locations"] != "":
        epitope_locations = read_b_epitopes(input_files["epitope_locations"])

    # Initialize parameters

    # Create host environment and add to viral sequences counter
    host = agents.create_host_env(
        founder_viruses,
        reference_sequence,
        float(params["replicative_cost"]),
        int(pop_samp.loc[0]["active_cell_count"]),
    )

    # Active cell counts for each generation
    active_cell_count = pop_samp["active_cell_count"]
    # Number of cells to sample for each generation
    n_to_samp = pop_samp["n_sample_active"]
    # Last sampled generation (don't have to continue simulation after this)
    last_sampled_gen = [index for index,
                        item in enumerate(n_to_samp) if item != 0][-1]

    # Loop through generations
    counts, fitness, seqs = host.loop_through_generations(
        active_cell_count,
        n_to_samp,
        last_sampled_gen,
        founder_viruses,
        nucleotides_order,
        substitution_probabilities,
        1 - exp(-params["mut_rate"]),
        1 - exp(-params["recomb_rate"]),
        params["prob_act_to_lat"],
        params["prob_lat_to_act"],
        params["prob_lat_die"],
        params["prob_lat_prolif"],
        conserved_sites,
        params["conserved_cost"],
        reference_sequence,
        float(params["replicative_cost"]),
        epitope_locations,
        params["seroconversion_time"],
        params["immune_response_proportion"],
        params["time_to_full_potency"],
        generator
    )

    # Make directories if they don't exist
    directory = dirname(argv[2])
    makedirs(directory, exist_ok=True)

    # Write counts to csv
    keys = [
        "generation",
        "active_cell_count",
        "latent_cell_count",
        "active_turned_latent",
        "latent_turned_active",
        "latent_died",
        "latent_proliferated",
        "number_mutations",
        "number_dual_inf",
        "mean_fitness_active",
        "mean_conserved_active",
        "mean_immune_active",
        "mean_replicative_active",
    ]
    with open(argv[2] + "counts.csv", "w") as outfile:
        w = writer(outfile)
        w.writerow(keys)
        w.writerows(zip(*[counts[key] for key in keys]))
    
    # Write fitness to csv
    keys = [
        "generation",
        "seq_id",
        "immune",
        "conserved",
        "replicative",
        "overall",
    ]
    with open(argv[2] + "fitness.csv", "w") as outfile:
        w = writer(outfile)
        w.writerow(keys)
        w.writerows(zip(*[fitness[key] for key in keys]))

    # Write sequences to fasta
    with open(argv[2] + "viral_seqs_active_CD4.fasta", "w") as outfile:
        for key, value in seqs.items():
            outfile.write(">" + key + "\n")
            outfile.write(value + "\n")

    end = time()
    print("Execution time = ", end - start)
