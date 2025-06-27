from scipy.linalg import expm
from numpy import fill_diagonal
from numpy import sum as npsum
from collections import namedtuple

from model.virus import HIV
from model.host import HostEnv

def prep_ref_conserved(founder_viruses, reference_sequence, conserved_sites):
    founder_virus_sequences = list(founder_viruses.values())
    assert len(founder_virus_sequences[0]) == len(reference_sequence)
    # remove any conserved sites that are variable in the founder sequence compared to the reference
    diff_sites = set({})
    for f in founder_virus_sequences:
        diff_sites.update({i for i, (left, right) in enumerate(
            zip(reference_sequence, f)) if left != right})
    [conserved_sites.pop(x, None) for x in diff_sites]
    # mask conserved sites so they aren't included in replicative
    # fitness computation
    reference_sequence = "".join(
        [x if i not in conserved_sites else "-" for i, x in enumerate(reference_sequence)])
    # number of sites to be compared for replicative fitness
    # len_ref_compare = len([1 for x in range(len(reference_sequence)) if reference_sequence[x] in ["A", "T", "C", "G"]])
    return reference_sequence, conserved_sites  # , len_ref_compare


def calc_nt_sub_probs_from_q(q, mut_rate):
    # Convert to probabilities
    sub_probs = expm(q * mut_rate)
    fill_diagonal(sub_probs, 0)
    sub_probs = sub_probs / npsum(sub_probs, axis=1, keepdims=True)

    # Make sure that the row and column labels are the same
    new_nucleotides_order = tuple(q.columns)
    old_nucleotides_order = tuple(q.index)
    assert "".join(old_nucleotides_order) == "".join(
        new_nucleotides_order
    ), "Nucleotide substitution matrix row and column labels are different"

    # Get probabilities as a list of tuples
    probabilities = sub_probs.tolist()
    
    SubProb = namedtuple('SubProb', ['order', 'probs'])

    return SubProb(order=new_nucleotides_order, probs=probabilities) 


def create_host_env(config, initial_cell_count):
    founder_viruses = [HIV(seq, config) for seq in config.founder_seqs.values()]
    return HostEnv(founder_viruses, initial_cell_count)


def create_b_epitope(start, end, max_fc):
    return BEpitope(int(start), int(end), float(max_fc))
  
class BEpitope:
    def __init__(self, start, end, max_fitness):
        self.start = start
        self.end = end
        self.max_fitness = max_fitness

    def __repr__(self):
        return "(%s to %s, maxfit: %s)" % (
            self.start, self.end, self.max_fitness)
            
class TEpitope:
    def __init__(self, start, escape_positions, recognized_aa):
        self.start = start
        self.end = start + max(escape_positions)*3
        self.escape_positions = dict(zip(escape_positions, recognized_aa))
        
    def __repr__(self):
        return "(start: %s, escape positions: %s)" % (self.start, self.escape_positions)

