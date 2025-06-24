from collections import namedtuple
from numpy.random import default_rng

from model.seed import set_python_seed


class Config:
    def __init__(self, last_sampled_gen, founder_seqs, 
        substitution_probabilities, prob_mut, prob_recomb,
        prob_act_to_lat, prob_lat_to_act, prob_lat_die, prob_lat_prolif,
        conserved_sites, conserved_cost, ref_seq, replicative_cost,
        epitope_locations, seroconversion_time, n_for_imm, gen_full_potency,
        seed
    ):
        self.last_sampled_gen = last_sampled_gen
        self.founder_seqs = founder_seqs
        self.seq_len = len(next(iter(founder_seqs.values())))
        
        if isinstance(substitution_probabilities, dict):
            # Create the namedtuple type
            SubProb = namedtuple("SubProb", substitution_probabilities.keys())
            # Convert dict to namedtuple
            substitution_probabilities = SubProb(**substitution_probabilities)
        
        self.substitution_probabilities = substitution_probabilities
        self.prob_mut = prob_mut
        self.prob_recomb = prob_recomb
        
        self.prob_act_to_lat = prob_act_to_lat
        self.prob_lat_to_act = prob_lat_to_act
        self.prob_lat_die = prob_lat_die
        self.prob_lat_prolif = prob_lat_prolif
        
        conserved_sites = {int(k): v.upper() for k, v in conserved_sites.items()}
        
        self.conserved_sites = conserved_sites
        self.conserved_cost = conserved_cost
        
        self.ref_seq = ref_seq
        self.replicative_cost = replicative_cost
        
        self.epitope_locations = epitope_locations
        self.seroconversion_time = seroconversion_time
        self.n_for_imm = n_for_imm
        self.gen_full_potency = gen_full_potency

        self.generator = set_python_seed(seed)
