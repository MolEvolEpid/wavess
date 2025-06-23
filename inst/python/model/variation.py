from numpy.random import default_rng
from random import sample
from random import choices
from collections import defaultdict
import numpy as np

def get_substitution(old_nucleotide, new_nucleotides_order, probabilities):
    try:
        idx = new_nucleotides_order.index(old_nucleotide)
    except ValueError:
        raise Exception(f"Unknown nucleotide {old_nucleotide}")
    return choices(new_nucleotides_order, probabilities[idx])[0]


def get_recomb_breakpoints(seq_len, num_cells, prob_recomb, is_sparse, base_prob, seed):
    """
    Efficient recombination breakpoint sampling that automatically chooses the fastest method.

    If prob_recomb is:
      - a scalar: uses constant-rate binomial sampling.
      - an array where most entries are identical: uses a sparse-optimized method for the rare positions.
      - otherwise, uses a fully vectorized method.

    Args:
        seq_len (int): Length of the genome sequence.
        num_cells (int): Number of cells.
        prob_recomb (float or array-like): Recombination probability per potential breakpoint.
        seed (int): Seed for random number generation.
        sparse_threshold (float): Fraction of non-base-rate entries below which to use the sparse method.

    Returns:
        num_cells_recomb (int): Number of dually infected cells in which recombination occurred.
        cell_breakpoints (list): List of breakpoint lists for each cell.
    """
    # Get random number generator
    rng = default_rng(seed)

    if is_sparse: # --- Sparse-optimized method ---
        # Bulk positions
        bulk_positions = np.where(prob_recomb == base_prob)[0]
        n_bulk = int(len(bulk_positions) * num_cells)
        n_recomb_bulk = rng.binomial(n_bulk, base_prob)
        if n_recomb_bulk > 0:
            flat_bulk_indices = rng.choice(n_bulk, n_recomb_bulk, replace=False)
            bulk_breakpoints = bulk_positions[flat_bulk_indices % len(bulk_positions)]
            bulk_cells = flat_bulk_indices // len(bulk_positions)
        else:
            bulk_breakpoints = np.array([], dtype=int)
            bulk_cells = np.array([], dtype=int)

        # Special positions
        special_idx = np.where(prob_recomb != base_prob)[0]
        special_breakpoints = []
        special_cells = []
        for idx in special_idx:
            events = rng.random(num_cells) < prob_recomb[idx]
            cells = np.flatnonzero(events)
            special_breakpoints.extend([idx] * len(cells))
            special_cells.extend(cells)
        special_breakpoints = np.array(special_breakpoints, dtype=int)
        special_cells = np.array(special_cells, dtype=int)

        # Combine
        # need to add 1 because want breakpoint to be after base
        breakpoints = np.concatenate([bulk_breakpoints, special_breakpoints]) + 1
        recomb_cells = np.concatenate([bulk_cells, special_cells])

        # Convert to flat indices for downstream calculation
        num_cells_recomb = len(set(recomb_cells))

    else:
        # --- Fully vectorized method ---
        per_bp_probs = np.tile(prob_recomb, num_cells)
        events = rng.random(per_bp_probs.size) < per_bp_probs
        flat_indices = np.flatnonzero(events)
        corrected_indices = flat_indices + flat_indices // (seq_len - 1) + 1
        recomb_cells = corrected_indices // seq_len
        breakpoints = corrected_indices % seq_len
        num_cells_recomb = len(set(recomb_cells))

    # Cross-over positions for each cell
    cell_breakpoints = defaultdict(list)
    for cell, breakpoint in zip(recomb_cells, breakpoints):
        cell_breakpoints[cell].append(breakpoint)
    return num_cells_recomb, list(cell_breakpoints.values())


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
