from collections import namedtuple
from numpy.random import default_rng
import numpy as np

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
        
        # determine type of recombination breakpoint method to do
        recrate_is_sparse = True
        base_prob = prob_recomb
        if isinstance(prob_recomb, (float, int)):
            prob_recomb = np.full(self.seq_len - 1, prob_recomb)
        else:
            prob_recomb = np.array(prob_recomb)
            sparse_threshold = 0.05
            if prob_recomb.shape[0] != self.seq_len - 1:
                raise ValueError("Length of per-breakpoint recombination rate must be seq_len-1.")
            # Check for sparseness: is there a dominant value?
            probs, nprob = np.unique(prob_recomb, return_counts=True)
            L = self.seq_len - 1
            maxidx = np.argmax(nprob)
            base_prob = probs[maxidx]
            recrate_is_sparse = ((L - nprob[maxidx]) / L) < sparse_threshold
        
        self.prob_recomb = prob_recomb
        self.base_prob = base_prob
        self.recrate_is_sparse = recrate_is_sparse
        
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
