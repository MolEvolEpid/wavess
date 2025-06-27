from collections import namedtuple
from numpy.random import default_rng
import numpy as np

from model.seed import set_python_seed


class Config:
    def __init__(self, last_sampled_gen, founder_seqs, 
        substitution_probabilities, prob_mut, prob_recomb,
        prob_act_to_lat, prob_lat_to_act, prob_lat_die, prob_lat_prolif,
        conserved_sites, conserved_cost, ref_seq, replicative_cost,
        b_epitope_locations, seroconversion_time, n_for_imm, b_gen_full_potency,
        t_epitope_locations, t_max_imm, t_gen_full_potency,
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
        
        self.b_epitope_locations = b_epitope_locations
        self.b_seroconversion_time = seroconversion_time
        self.b_n_for_imm = n_for_imm
        self.b_gen_full_potency = b_gen_full_potency
        
        self.t_epitope_locations = t_epitope_locations
        self.t_max_imm = t_max_imm
        self.t_gen_full_potency = t_gen_full_potency
        
        t_epitope_mask = np.full(self.seq_len, -1, dtype=int)
        recognition_motifs = []
        epi_num = 0
        for epi in self.t_epitope_locations:
            for position in epi.escape_positions:
                t_epitope_mask[(epi.start + (position-1)*3):(epi.start + (position-1)*3+3)] = epi_num
            recognition_motifs.append(''.join([epi.escape_positions[pos] if pos in epi.escape_positions else 'N' for pos in range(1, max(epi.escape_positions.keys())+1)]))
            epi_num += 1
        self.t_epitope_mask = t_epitope_mask
        self.t_recognition_motifs = recognition_motifs

        self.generator = set_python_seed(seed)
