from collections import namedtuple
from numpy.random import default_rng
import numpy as np
from numpy import fill_diagonal
from numpy import sum as npsum
from scipy.linalg import expm
from math import exp

from model.virus import HIV
from model.epitope import BEpitope, TEpitope
from model.host import HostEnv
from model.seed import set_python_seed


class Config:
    def __init__(self, generation_time, last_sampled_gen, founder_seqs, 
        q, mut_rate, recomb_rate,
        act_to_lat_rate, lat_to_act_rate, lat_die_rate, lat_prolif_rate,
        conserved_sites, conserved_cost, ref_seq, replicative_cost,
        b_epitope_locations, b_immune_start_day, b_n_for_imm, b_days_full_potency,
        t_epitope_locations, t_max_imm, #t_days_full_potency,
        seed
    ):
        self.last_sampled_gen = last_sampled_gen
        self.founder_seqs = founder_seqs
        self.seq_len = len(next(iter(founder_seqs.values())))
        
        # get substitution probabilities
        self.q = q
        self.mut_rate = mut_rate
        self.calc_nt_sub_probs_from_q() 
        
        # convert rates to per-generation probabilities
        # these are in generations already
        self.prob_mut = self.rate_to_probability(mut_rate)
        if np.ndim(recomb_rate) != 0:
            recomb_rate = np.asarray(recomb_rate)
        self.prob_recomb = self.rate_to_probability(recomb_rate)
        # these are in days, so we convert to generations first
        self.prob_act_to_lat = self.rate_to_probability(act_to_lat_rate / generation_time)
        self.prob_lat_to_act = self.rate_to_probability(lat_to_act_rate / generation_time)
        self.prob_lat_die = self.rate_to_probability(lat_die_rate / generation_time)
        self.prob_lat_prolif = self.rate_to_probability(lat_prolif_rate / generation_time)

        # determine type of recombination breakpoint method to do
        self.recrate_is_sparse = True
        self.base_prob = self.prob_recomb

        if isinstance(self.prob_recomb, (float, int)):
            self.prob_recomb = np.full(self.seq_len - 1, self.prob_recomb)
        else:
            self.prob_recomb = np.array(self.prob_recomb)
            sparse_threshold = 0.05
            if self.prob_recomb.shape[0] != self.seq_len - 1:
                raise ValueError("Length of per-breakpoint recombination rate must be seq_len-1.")
            # Check for sparseness: is there a dominant value?
            probs, nprob = np.unique(self.prob_recomb, return_counts=True)
            L = self.seq_len - 1
            maxidx = np.argmax(nprob)

            self.base_prob = probs[maxidx]
            self.recrate_is_sparse = ((L - nprob[maxidx]) / L) < sparse_threshold
        
        self.b_epitope_locations = b_epitope_locations
        if b_epitope_locations is not None:
            self.b_epitope_locations = [BEpitope(int(epi[1]), int(epi[2]), float(epi[3])) for epi in b_epitope_locations.itertuples()]
        self.b_immune_start_gen = b_immune_start_day / generation_time
        self.b_n_for_imm = b_n_for_imm
        # convert days to generations
        self.b_gen_full_potency = b_days_full_potency / generation_time
        
        # epitope, start, escape position, recognized_aa
        starts = []
        t_days_to_full_potency = []
        positions = []
        aas = []
        for row in t_epitope_locations.itertuples():
            starts.append(int(row[1]))
            t_days_to_full_potency.append(int(row[2]))
            positions.append(int(row[3]))
            aas.append(row[4])
        epitopes = []
        for start in set(starts):
            epitopes.append(TEpitope(
                start, np.unique(np.array(t_days_to_full_potency)[np.array(starts) == start])[0] / generation_time, np.array(positions)[np.array(starts) == start], np.array(aas)[np.array(starts) == start]))
        self.t_epitope_locations = epitopes

        self.t_max_imm = t_max_imm
        # self.t_gen_full_potency = t_days_full_potency / generation_time
        
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
        
        self.conserved_sites = {int(k): v.upper() for k, v in conserved_sites.items()}
        self.conserved_cost = conserved_cost
        
        self.ref_seq = ref_seq
        self.replicative_cost = replicative_cost
        
        if(len(ref_seq) > 0):
            self.prep_ref_conserved() # MAKE A TEST TO BE SURE THIS WORKS
        
        # Get founder viruses
        founder_viruses = [HIV(seq, self) for seq in self.founder_seqs.values()]
        # Create host environment
        self.host = HostEnv(founder_viruses, len(founder_viruses))

        self.generator = set_python_seed(seed)


    # Probability that an odd number of events occur in the time period assuming a Poisson process
    # asked chatgpt then checked on wikipedia and with wolfram alpha
    # this also means that if there are n bases between two bases in the sequence 
    # (e.g. if you're modeling pol and env), then you can input the rate as (n+1)*base_rate 
    # (assuming independence)
    def rate_to_probability(self, rate, time = 1):
        return (1 - np.exp(-2 * rate * time))/2 
    
    def calc_nt_sub_probs_from_q(self):
        # Convert to probabilities
        sub_probs = expm(self.q * self.mut_rate)
        fill_diagonal(sub_probs, 0)
        sub_probs = sub_probs / npsum(sub_probs, axis=1, keepdims=True)
    
        # Make sure that the row and column labels are the same
        new_nucleotides_order = tuple(self.q.columns)
        assert "".join(tuple(self.q.index)) == "".join(
            new_nucleotides_order
        ), "Nucleotide substitution matrix row and column labels are different"
    
        # Get probabilities as a list of tuples
        probabilities = sub_probs.tolist()
        
        SubProb = namedtuple('SubProb', ['order', 'probs'])
    
        self.substitution_probabilities = SubProb(order=new_nucleotides_order, probs=probabilities) 
    

    def prep_ref_conserved(self):
        founder_virus_sequences = list(self.founder_seqs.values())
        assert len(founder_virus_sequences[0]) == len(self.ref_seq), "Founder sequence and reference sequence must be the same length"
        # remove any conserved sites that are variable in the founder sequence compared to the reference
        diff_sites = set({})
        for f in founder_virus_sequences:
            diff_sites.update({i for i, (left, right) in enumerate(
                zip(self.ref_seq, f)) if left != right})
        # pop these out of conserved if they exist in there
        [self.conserved_sites.pop(x, None) for x in diff_sites]
        # mask conserved sites so they aren't included in replicative fitness computation
        self.ref_seq = "".join(
            [x if i not in self.conserved_sites else "-" for i, x in enumerate(self.ref_seq)])
            
