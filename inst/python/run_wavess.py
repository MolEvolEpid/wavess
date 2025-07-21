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
import pandas as pd
import os

# Import custom classes and functions
import model.__init__ as agents

# Functions to read in data

# Function to read in sequence data
def get_sequences(filename):
    founder_virus_sequences = [
        str(fasta.seq).upper() for fasta in SeqIO.parse(open(filename), "fasta")
    ]
    # The lengths of the founder sequences must be the same
    assert False not in [
        len(i) == len(founder_virus_sequences[0]) for i in founder_virus_sequences
    ], "Founder virus sequences must be of the same length"
    return founder_virus_sequences

def load_recombination_rate(recomb_rate, seq_len):
    if isinstance(recomb_rate, (float, int)):
        return recomb_rate
    elif isinstance(recomb_rate, str):
        if not os.path.isfile(recomb_rate):
            raise FileNotFoundError(f"Recombination rate file not found: {recomb_rate}")
        arr = np.loadtxt(recomb_rate)
        if arr.shape[0] != seq_len - 1:
            raise ValueError(f"Expected {seq_len-1} rates in {recomb_rate}, found {arr.shape[0]}")
        return arr
    else:
        raise TypeError("recomb_rate must be a float, int, or a filepath string")
      

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

    # Read in config file
    config = safe_load(open(argv[1]))

    ### Read in information ###
    input_files = config["input_files"]
    params = config["parameters"]

    ## Required inputs ##

    # Generation information, including cd4 counts and sampling times and
    # numbers
    inf_pop_size = read_csv(input_files["inf_pop_size"], index_col=["generation"])
    samp_scheme = read_csv(input_files["samp_scheme"])
    samp_scheme['generation'] = pd.to_numeric(round(samp_scheme['day']/params['generation_time']), downcast = 'integer')
    samp_scheme = samp_scheme.drop(['day'], axis = 1)
    samp_scheme = samp_scheme.set_index(['generation'])
    pop_samp = pd.concat([inf_pop_size, samp_scheme], axis = 1)
    pop_samp['n_sample_active'] = pop_samp['n_sample_active'].fillna(0)
    pop_samp['n_sample_latent'] = pop_samp['n_sample_latent'].fillna(0)
    pop_samp = pop_samp[['active_cell_count', 'n_sample_active', 'n_sample_latent']]

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
    q_matrix = read_csv(input_files["q"], index_col="nt_from")

    ## Optional inputs related to selection ##

    # Conserved sites
    conserved_sites = {}
    if input_files["conserved_sites"] != "":
        conserved_sites = read_csv(input_files["conserved_sites"], header=0).set_index("position")["nucleotide"].to_dict() 

    # Reference sequence
    reference_sequence = ""
    if input_files["ref_seq"] != "":
        reference_sequence = get_sequences(input_files["ref_seq"])[0]

    # Epitope start positions and max fitness cost.
    b_epitope_locations = None
    if input_files["b_epitope_locations"] != "":
        b_epitope_locations = read_csv(input_files["b_epitope_locations"])
        
    t_epitope_locations = None
    if input_files["t_epitope_locations"] != "":
        t_epitope_locations = read_csv(input_files["t_epitope_locations"])

    # Initialize parameters
    
    # Check that fitness parameters are in the range [0,1)
    if float(params["conserved_cost"]) < 0 or float(params["conserved_cost"]) >= 1:
      raise ValueError("conserved fitness cost must be in the range [0,1)")
    if float(params["replicative_cost"]) < 0 or float(params["replicative_cost"]) >= 1:
      raise ValueError("replicative fitness cost must be in the range [0,1)")

    # Prepare recombination rate (variable or constant)
    seq_len = len(reference_sequence)
    prob_recomb = load_recombination_rate(params["recomb_rate"], seq_len)
    if isinstance(prob_recomb, (list, np.ndarray)):
        # User provided an array for variable recombination rate
        prob_recombination = np.asarray(prob_recomb) #1 - exp(-params["recomb_rate"])
        if prob_recombination.shape[0] != seq_len - 1:
            raise ValueError("Variable recombination rate must have length seq_len-1")
    else:
        # Use constant rate
        prob_recombination = 1 - exp(-prob_recomb) 

    # Active cell counts for each generation
    active_cell_count = pop_samp["active_cell_count"]
    # Number of cells to sample for each generation
    n_to_samp_active = pop_samp["n_sample_active"]
    n_to_samp_latent = pop_samp["n_sample_latent"]
    # Last sampled generation (don't have to continue simulation after this)
    sampled_gens_active = [index for index, item in enumerate(n_to_samp_active) if item != 0]
    sampled_gens_latent = [index for index, item in enumerate(n_to_samp_latent) if item != 0]
    last_sampled_gen_active = 0
    if len(sampled_gens_active) > 0:
        last_sampled_gen_active = sampled_gens_active[-1]
    last_sampled_gen_latent = 0
    if len(sampled_gens_latent) > 0:
        last_sampled_gen_latent = sampled_gens_latent[-1]
    if last_sampled_gen_active == 0 and last_sampled_gen_latent == 0:
        raise ValueError("at least one sequence must be sampled from the active or latent reservoir")
    last_sampled_gen = max(last_sampled_gen_active, last_sampled_gen_latent)
    
    model_config = agents.Config(params["generation_time"], last_sampled_gen, founder_viruses, q_matrix,
                                 params["mut_rate"], params["recomb_rate"],
                                 params["act_to_lat"], params["lat_to_act"], params["lat_die"], params["lat_prolif"],
                                 conserved_sites, params["conserved_cost"],
                                 reference_sequence, float(params["replicative_cost"]),
                                 b_epitope_locations, params["b_immune_start_day"],
                                 params["b_n_for_imm"], params["b_days_to_full_potency"],
                                 t_epitope_locations, params["t_max_immune_cost"], params["t_days_to_full_potency"],
                                 s)

    # Loop through generations
    counts, fitness, seqs_active, seqs_latent = model_config.host.loop_through_generations(
        active_cell_count,
        n_to_samp_active,
        n_to_samp_latent,
        model_config
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
        "number_recombinations",
        "mean_fitness_active",
        "mean_conserved_active",
        "mean_b_immune_active",
        "mean_t_immune_active",
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
        "b_immune",
        "t_immune",
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
        for key, value in seqs_active.items():
            outfile.write(">" + key + "\n")
            outfile.write(value + "\n")
    with open(argv[2] + "viral_seqs_latent_CD4.fasta", "w") as outfile:
        for key, value in seqs_latent.items():
            outfile.write(">" + key + "\n")
            outfile.write(value + "\n")

    end = time()
    print("Execution time = ", end - start)
