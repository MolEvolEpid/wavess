from numpy.random import default_rng
from random import sample
from random import choices
from collections import defaultdict

def get_substitution(old_nucleotide, new_nucleotides_order, probabilities):
    try:
        idx = new_nucleotides_order.index(old_nucleotide)
    except ValueError:
        raise Exception(f"Unknown nucleotide {old_nucleotide}")
    return choices(new_nucleotides_order, probabilities[idx])[0]


def get_recomb_breakpoints(seq_len, num_cells, prob_recomb, seed):
    # Get random number generator
    rng = default_rng(seed)
    # Number of potential breakpoints
    n_potential_bps = int(seq_len * num_cells - num_cells)
    # Number of cross-over events
    n_recomb = int(rng.binomial(n_potential_bps, prob_recomb))
    # Indices for cross-over events
    indices = sample(range(n_potential_bps), n_recomb)
    # Actual indices of cross-over events (have to skip numbers between cells)
    corrected_indices = [x + x // (seq_len - 1) + 1 for x in indices]
    # Cell in which cross-over event occurred
    recomb_cells = [x // seq_len for x in corrected_indices]
    # Recombination breakpoint in cell
    breakpoints = [x % seq_len for x in corrected_indices]
    # Total number of dually infected cells in which recombination occurred
    num_cells_recomb = len(set(recomb_cells))
    # Cross-over positions for each cell
    cell_breakpoints = defaultdict(list)
    for cell, breakpoint in zip(recomb_cells, breakpoints):
        cell_breakpoints[cell].append(breakpoint)
    return num_cells_recomb, cell_breakpoints.values()


def get_recombined_sequence(sequence1, sequence2, breakpoints):
    seq_len = len(sequence1)
    assert seq_len not in breakpoints, "Breakpoint at end of sequence"
    assert 0 not in breakpoints, "Breakpoint at beginning of sequence"
    cross_over_positions = sorted(breakpoints) + [seq_len]
    # Switch between copying first and second strand
    curr_pos = cross_over_positions.pop(0)
    # Start with sequence1 first, then alternate
    recombined_sequence = sequence1[:curr_pos]
    next_strand = 2
    while len(cross_over_positions) > 0:
        prev_pos = curr_pos
        curr_pos = cross_over_positions.pop(0)
        if next_strand == 2:
            recombined_sequence += sequence2[prev_pos:curr_pos]
            next_strand = 1
        else:
            recombined_sequence += sequence1[prev_pos:curr_pos]
            next_strand = 2
    return recombined_sequence, breakpoints


def muts_rel_ref(nuc_sequence, ref_sequence):
    return sum(1 for x, y in zip(nuc_sequence, ref_sequence) if x != y and y in {"A", "T", "C", "G"})


def get_conserved_sites_mutated(v1_muts, v2_muts, cross_over_positions, seq_len):
    cross_over_positions = cross_over_positions + [seq_len]
    n_muts_conserved_virus1 = len(v1_muts)
    n_muts_conserved_virus2 = len(v2_muts)
    muts_in_conserved = set()
    if n_muts_conserved_virus1 != 0 or n_muts_conserved_virus2 != 0:
        curr_pos = 0
        next_strand = 1
        while len(cross_over_positions) > 0:
            prev_pos = curr_pos
            curr_pos = cross_over_positions.pop(0)
            if next_strand == 1:
                muts = set(x for x in v1_muts if x >= prev_pos and x < curr_pos)
                for x in muts:
                    muts_in_conserved.add(x)
                next_strand = 2
            else:
                muts = set(x for x in v2_muts if x >= prev_pos and x < curr_pos)
                for x in muts:
                    muts_in_conserved.add(x)
                next_strand = 1

    return muts_in_conserved
