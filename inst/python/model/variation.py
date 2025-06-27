from numpy.random import default_rng
from random import sample
from random import choices
from collections import defaultdict
import numpy as np


def get_substitution(old_nucleotide, subprob):
    try:
        idx = subprob.order.index(old_nucleotide)
    except ValueError:
        raise Exception(f"Unknown nucleotide {old_nucleotide}")
    return choices(subprob.order, subprob.probs[idx])[0]
  

def get_var_positions(num_cells, var_type, config, gen = None):
    """
    Efficient mutation position and recombination breakpoint sampling that automatically chooses the fastest method.

    If prob is:
      - a scalar: uses constant-rate binomial sampling.
      - an array where most entries are identical: uses a sparse-optimized method for the rare positions.
      - otherwise, uses a fully vectorized method.

    Args:
        num_cells (int): Number of cells.
        var_type: mutation or recombination
        config: config

    Returns:
        num_cells_recomb (int): Number of dually infected cells in which recombination occurred.
        cell_breakpoints (list): List of breakpoint lists for each cell.
    """
    # Get random number generator
    rng = default_rng(config.generator)
    
    assert var_type in ['mutation', 'recombination']
    
    if var_type == 'mutation':
        prob = config.prob_mut
        base_prob = config.base_prob_mut 
        is_sparse = config.mutrate_is_sparse 
        if config.mut_rate_scalar is not None:
            prob = prob * config.mut_rate_scalar[gen]
            base_prob = base_prob * config.mut_rate_scalar[gen]
    else: #if var_type == 'recombination':
        prob = config.prob_recomb
        base_prob = config.base_prob_recomb
        is_sparse = config.recrate_is_sparse 

    if is_sparse: # --- Sparse-optimized method ---
        # Bulk positions
        bulk_positions = np.where(prob == base_prob)[0]
        n_bulk = int(len(bulk_positions) * num_cells)
        n_var_bulk = rng.binomial(n_bulk, base_prob)
        flat_bulk_indices = rng.choice(n_bulk, n_var_bulk, replace=False)
        bulk_var_positions = bulk_positions[flat_bulk_indices % len(bulk_positions)]
        bulk_cells = flat_bulk_indices // len(bulk_positions)

        # Special positions
        special_idx = np.where(prob != base_prob)[0]
        special_positions = []
        special_cells = []
        for idx in special_idx:
            events = rng.random(num_cells) < prob[idx]
            cells = np.flatnonzero(events)
            special_positions.extend([idx] * len(cells))
            special_cells.extend(cells)
        special_positions = np.array(special_positions, dtype=int)
        special_cells = np.array(special_cells, dtype=int)

        # Combine
        positions = np.concatenate([bulk_var_positions, special_positions]) 
        # need to add 1 for recombinations because want breakpoint to be after base
        if var_type == 'recombination':
            positions += 1
        cells = np.concatenate([bulk_cells, special_cells])

        # Convert to flat indices for downstream calculation
        num_cells_var = len(set(cells))

    else:
        # --- Fully vectorized method ---
        per_bp_probs = np.tile(prob, num_cells)
        events = rng.random(per_bp_probs.size) < per_bp_probs
        flat_indices = np.flatnonzero(events)
        if var_type == 'recombination':
            corrected_indices = flat_indices + flat_indices // (config.seq_len - 1) + 1
        else:
            corrected_indices = flat_indices + flat_indices // (config.seq_len) 
        cells = corrected_indices // config.seq_len
        positions = corrected_indices % config.seq_len
        #cells, positions = divmod(
        #        corrected_indices, config.seq_len
        #    )  # Find cell and position to mutate
        num_cells_var = len(set(cells))

    # Var positions for each cell
    cell_positions = defaultdict(list)
    for cell, position in zip(cells, positions):
        cell_positions[cell].append(position)
    return num_cells_var, cell_positions


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
