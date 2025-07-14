from collections import namedtuple
from numpy.random import default_rng
from numpy import fill_diagonal
from numpy import sum as npsum
from scipy.linalg import expm
from math import exp

from model.virus import HIV
from model.epitope import Epitope
from model.host import HostEnv
from model.seed import set_python_seed


class Config:
    def __init__(self, generation_time, last_sampled_gen, founder_seqs, 
        q, mut_rate, recomb_rate,
        act_to_lat_rate, lat_to_act_rate, lat_die_rate, lat_prolif_rate,
        conserved_sites, conserved_cost, ref_seq, replicative_cost,
        epitope_locations, immune_start_day, n_for_imm, days_full_potency,
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
        self.prob_recomb = self.rate_to_probability(recomb_rate)
        # these are in days, so we convert to generations first
        self.prob_act_to_lat = self.rate_to_probability(act_to_lat_rate / generation_time)
        self.prob_lat_to_act = self.rate_to_probability(lat_to_act_rate / generation_time)
        self.prob_lat_die = self.rate_to_probability(lat_die_rate / generation_time)
        self.prob_lat_prolif = self.rate_to_probability(lat_prolif_rate / generation_time)
        
        self.conserved_sites = {int(k): v.upper() for k, v in conserved_sites.items()}
        self.conserved_cost = conserved_cost
        
        self.ref_seq = ref_seq
        self.replicative_cost = replicative_cost
        
        if(len(ref_seq) > 0):
            self.prep_ref_conserved() # MAKE A TEST TO BE SURE THIS WORKS
        
        self.epitope_locations = epitope_locations
        if epitope_locations is not None:
            self.epitope_locations = [Epitope(int(epi[1]), int(epi[2]), float(epi[3])) for epi in epitope_locations.itertuples()]
        self.n_for_imm = n_for_imm
         # convert days to generations
        self.gen_full_potency = days_full_potency / generation_time
        self.seroconversion_time = immune_start_day / generation_time

        # Get founder viruses
        founder_viruses = [HIV(seq, self) for seq in self.founder_seqs.values()]
        # Create host environment
        self.host = HostEnv(founder_viruses, len(founder_viruses))

        self.generator = set_python_seed(seed)


    # Probability that at least one event occurs in the time period assuming a Poisson process
    def rate_to_probability(self, rate, time = 1):
        return 1 - exp(-rate * time)
    
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
        [self.conserved_sites.pop(x, None) for x in diff_sites]
        # mask conserved sites so they aren't included in replicative fitness computation
        self.ref_seq = "".join(
            [x if i not in self.conserved_sites else "-" for i, x in enumerate(self.ref_seq)])
            
