from copy import deepcopy
from random import seed
from random import choice
from random import choices
from random import sample
from random import randrange
from numpy.random import default_rng
from numpy import arange
from numpy import where
from numpy import concatenate
from collections import defaultdict
from collections import Counter

def set_python_seed(s):
  if s != None:
    # Set seed (need to set a seed for random and for numpy.random)
    seed(s)
    generator = default_rng(int(s))
  else:
    generator = default_rng()
  return generator

def create_host_env(founder_seqs, ref_seq, conserved_sites, rep_fit, rf_exp, initial_cell_count):
  founder_viruses = [HIV(seq, ref_seq, conserved_sites, rep_fit, rf_exp) for seq in founder_seqs]
  return HostEnv(founder_viruses, initial_cell_count)

def create_epitope(start, end, max_fc):
  return Epitope(int(start), int(end), float(max_fc))

def get_substitution(old_nucleotide, new_nucleotides_order, probabilities):
    if old_nucleotide == 'A':
        return choices(new_nucleotides_order, probabilities[new_nucleotides_order.index('A')])[0]
    elif old_nucleotide == 'C':
        return choices(new_nucleotides_order, probabilities[new_nucleotides_order.index('C')])[0]
    elif old_nucleotide == 'G':
        return choices(new_nucleotides_order, probabilities[new_nucleotides_order.index('G')])[0]
    elif old_nucleotide == 'T':
        return choices(new_nucleotides_order, probabilities[new_nucleotides_order.index('T')])[0]
    else:
        raise Exception("Unknown nucleotide %s" % old_nucleotide)

def get_recomb_breakpoints(seq_len, num_cells, prob_recomb, seed):
    # Get random number generator
    rng = default_rng(seed)
    # Number of potential breakpoints
    n_potential_bps = int(seq_len*num_cells-num_cells)
    # Number of cross-over events
    n_recomb = int(rng.binomial(n_potential_bps, prob_recomb))
    # Indices for cross-over events
    indices = sample(range(n_potential_bps), n_recomb)
    # Actual indices of cross-over events (have to skip numbers between cells)
    corrected_indices = [x+x//(seq_len-1)+1 for x in indices]
    # Cell in which cross-over event occurred
    recomb_cells = [x//seq_len for x in corrected_indices]
    # Recombination breakpoint in cell
    breakpoints = [x%seq_len for x in corrected_indices]
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
    while(len(cross_over_positions) > 0):
        prev_pos = curr_pos
        curr_pos = cross_over_positions.pop(0)
        if next_strand == 2:
            recombined_sequence += sequence2[prev_pos:curr_pos]
            next_strand = 1
        else:
            recombined_sequence += sequence1[prev_pos:curr_pos]
            next_strand = 2
    return recombined_sequence, breakpoints


def conserved_fitness_cost(num_mutations_in_conserved_sites, cost_per_mutation_in_conserved_site):
    return 1 - pow((1 - cost_per_mutation_in_conserved_site), num_mutations_in_conserved_sites)

# UPDATE THIS TO NOT TAKE EXPONENT IF WE END UP KEEPING IT THE WAY IT IS NOW
def replicative_fitness_cost(s1, s2, conserved_sites, rf_exp):
    assert len(s1) == len(s2), "In replicative fitness, sequences are not of same length!"
    #return (sum(1 for x,y in zip(s1, s2) if x != y and y in ['A','T','C','G']))/len(s1)
    #return len([1 for x,y in zip(s1, s2) if x != y and y in ['A','T','C','G']])/len(s1) # slightly faster
    # don't include conserved sites in calculation
    s1 = [s1[x] for x in range(len(s1)) if x not in conserved_sites]
    s2 = [s2[x] for x in range(len(s2)) if x not in conserved_sites]
    n_compare = len([1 for x in range(len(s2)) if s2[x] in ['A','T','C','G']])
    return (len([1 for x,y in zip(s1, s2) if x != y and y in ['A','T','C','G']])/n_compare) #**rf_exp


def normalize(likelihoods):
    total = sum(likelihoods)
    return [likelihood / total for likelihood in likelihoods]

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


class Epitope:
    def __init__(self, start, end, max_fitness):
        self.start = start
        self.end = end
        self.max_fitness = max_fitness

    def __repr__(self):
        return "(%s to %s, maxfit: %s)" % (self.start, self.end, self.max_fitness)


class HIV:
    def __init__(self, nuc_seq, reference_sequence, conserved_sites, replicative_fitness, rf_exp):
        # Make sure values supplied are as expected
        assert isinstance(nuc_seq, str), "Nucleotide sequence needs to be a string"

        # Initialize instance variables
        self.nuc_sequence = nuc_seq
        self.conserved_sites_mutated = set()

        # Tracking different components of fitness
        self.immune_fitness_cost = 0
        self.conserved_fitness_cost = 0
        self.replicative_fitness_cost = 0
        if replicative_fitness:
            self.replicative_fitness_cost = replicative_fitness_cost(self.nuc_sequence, reference_sequence, conserved_sites, rf_exp)
        self.fitness = 1

    def __repr__(self):
        return_str = "HIV with sequence %s" % self.nuc_sequence
        return return_str

    def mutate(self, position_to_mutate, nucleotides_order, substitution_probabilities, conserved_sites, 
            cost_per_mutation_in_conserved_site, reference_sequence, conserved_fitness, replicative_fitness, rf_exp):
        # Figure out what the mutation is based on mutation probabilities
        old_nucleotide = self.nuc_sequence[position_to_mutate]
        new_nucleotide = get_substitution(old_nucleotide, nucleotides_order, substitution_probabilities)
        # Fastest way I can find to mutate - add part before to new nt to part after
        self.nuc_sequence = self.nuc_sequence[:position_to_mutate] + new_nucleotide + self.nuc_sequence[position_to_mutate+1:]

        if conserved_fitness:
            # If mutation is in a conserved site, update mutations in conserved sites and conserved sites fitness
            # Note that we are not tracking whether there are multiple mutations at a single conserved site
            if position_to_mutate in conserved_sites:
                self.conserved_sites_mutated.add(position_to_mutate)
                self.conserved_fitness_cost = conserved_fitness_cost(len(self.conserved_sites_mutated), cost_per_mutation_in_conserved_site)

        # Update replicative fitness cost only if needed
        if replicative_fitness:
            ref_base = reference_sequence[position_to_mutate]
            prev_comp = old_nucleotide == ref_base
            if prev_comp != (ref_base == new_nucleotide):
                self.replicative_fitness_cost = replicative_fitness_cost(self.nuc_sequence, reference_sequence, conserved_sites, rf_exp)


class InfectedCD4:
    def __init__(self, infecting_virus: HIV, active: bool):
        assert isinstance(active, bool), "Provide bool value for whether CD4 T cell is active (vs latent)"
        self.active = active
        self.infecting_virus = infecting_virus

    def __repr__(self):
        return "Infected CD4. Active: %s. %s" % (
            str(self.active), str(self.infecting_virus))

    def become_latent(self):
        self.active = False  # Set to False to indicate the cell is latent

    def become_active(self):
        self.active = True  # Set to True to indicate the cell is actively proliferating


class HostEnv:  # This is the 'compartment' where the model dynamics take place
    def __init__(self, founder_viruses: HIV, NC: int): 
        self.epitope_variants_translated = defaultdict(lambda: '') # store epitope translations 
        
        # Create NC actively proliferating infected cells
        self.C = [InfectedCD4(deepcopy(choice(founder_viruses)), True) for i in range(NC)]
        #self.C = [InfectedCD4(deepcopy(founder_viruses), True) for i in range(NC)]

        # Initiate latent infected cells
        self.L = list()

        # Track epitopes recognized by the immune system
        # key -> epitope sequence variant, value -> generation when it started being recognized
        self.epitopes_recognition_generation = defaultdict(lambda: 0)
        # Immune cost of epitopes
        self.cross_reactive_epitope_cost_frac = defaultdict(lambda: 0)

    def __repr__(self):
        return_str = "Host has %s active and %s latent infected cells\n" % (str(len(self.C)), str(len(self.L)))
        return_str += "%s epitopes recognized" % len(self.epitopes_recognition_generation)
        return return_str

    # Translate nucleotide sequence to amino acid sequence
    # Nucleotide sequence length must be a multiple of 3
    # Modified from https://www.geeksforgeeks.org/dna-protein-python-3/
    def translate(self, seq_to_translate):
        table = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }
        assert len(seq_to_translate)%3 == 0, "The nucleotide sequence length is not a multiple of 3"
        protein = ""
        for i in range(0, len(seq_to_translate), 3):
            codon = seq_to_translate[i:i + 3]
            protein += table[codon]
        return protein

    def update_epitopes_recognized(self, current_generation, epitopes, immune_response_proportion, time_to_full_potency, seed):
        # Get random number generator
        rng = default_rng(seed)
        # Immune strengths
        immune_strength_dict = {k:min((current_generation - v)/time_to_full_potency, 1) for k,v in self.epitopes_recognition_generation.items()}
        # Number of active infected cells
        num_active_cd4 = len(self.C)
        # Get unique sequences
        seqs = defaultdict(lambda: 0)
        for index in range(num_active_cd4):
            seqs[self.C[index].infecting_virus.nuc_sequence] += 1
        # Track frequency of epitope sequence variants in sample
        epitope_variants = defaultdict(lambda: 0)
        for seq, freq in seqs.items():
            for epi in epitopes:
                epitope_variants[seq[epi.start: epi.end]] += freq

        # Translate epitopes, get frequencies, and store translations
        aa_epitope_variants = defaultdict(lambda: 0)
        for variant, frequency in epitope_variants.items():
            # Translate and add to dictionary if new variant
            if variant not in self.epitope_variants_translated:
                self.epitope_variants_translated[variant] = self.translate(variant)
            # Add variant frequency to dictionary
            aa_epitope_variants[self.epitope_variants_translated[variant]] += frequency
        
        # Get cross-reactive strength for each new epitope
        for variant, frequency in aa_epitope_variants.items():
            if variant not in self.epitopes_recognition_generation and len(self.epitopes_recognition_generation) != 0:
                if variant not in self.cross_reactive_epitope_cost_frac:
                   most_similar = ''
                   n_diff_min = len(variant)
                   for v in self.epitopes_recognition_generation:
                       if variant != v and len(variant) == len(v):
                           n_diff = len([1 for x,y in zip(variant, v) if x != y])
                           if n_diff <= n_diff_min:
                               n_diff_min = n_diff
                               most_similar = v
                   # EVENTUALLY UPDATE THIS TO BE HOWEVER WE WANT CROSS-REACTIVITY TO ACTUALLY WORK
                   self.cross_reactive_epitope_cost_frac[variant] = immune_strength_dict[most_similar]*rng.beta(1, n_diff_min**2) #0.5/n_diff_min
        
        # Mark frequent epitopes as recognized
        for variant, frequency in aa_epitope_variants.items():
            if (frequency / num_active_cd4) >= immune_response_proportion:  # common enough to be recognized
                if variant not in self.epitopes_recognition_generation:  # if first time epitope is being recognized
                    if variant not in self.cross_reactive_epitope_cost_frac:
                        self.epitopes_recognition_generation[variant] = current_generation
                    else:
                       self.epitopes_recognition_generation[variant] = int(current_generation - self.cross_reactive_epitope_cost_frac[variant]*time_to_full_potency)

    def update_immune_fitness(self, epitopes, current_generation, time_to_full_potency):
        immune_strength_dict = {k:min((current_generation - v)/time_to_full_potency, 1) for k,v in self.epitopes_recognition_generation.items()}
        for CD4_index in range(len(self.C)):
            # Fitness cost due to immune response
            max_epitope_fitness_cost = 0 # by default the virus is super fit
            for epi in epitopes:
                # get translated epitope
                epitope_sequence = self.epitope_variants_translated[self.C[CD4_index].infecting_virus.nuc_sequence[epi.start:epi.end]]
                if epitope_sequence in self.epitopes_recognition_generation:
                    timed_fitness = epi.max_fitness * immune_strength_dict[epitope_sequence]
                    max_epitope_fitness_cost = max(max_epitope_fitness_cost, timed_fitness)
                elif epitope_sequence in self.cross_reactive_epitope_cost_frac:
                   timed_fitness = epi.max_fitness * self.cross_reactive_epitope_cost_frac[epitope_sequence]
                   max_epitope_fitness_cost = max(max_epitope_fitness_cost, timed_fitness)
            self.C[CD4_index].infecting_virus.immune_fitness_cost = max_epitope_fitness_cost 


    def get_fitness_of_infecting_virus(self, CD4_index, rf_exp):

        # Compute overall fitness
        viral_fitness = (1 - self.C[CD4_index].infecting_virus.immune_fitness_cost) * \
                (1 - self.C[CD4_index].infecting_virus.conserved_fitness_cost) * \
                (1 - self.C[CD4_index].infecting_virus.replicative_fitness_cost)**rf_exp

        # Update the virus object with the fitness costs
        self.C[CD4_index].infecting_virus.fitness = viral_fitness

        return viral_fitness

    def get_fitness(self, rf_exp):
        fitness = []
        for CD4_index in range(len(self.C)):
            fit = self.get_fitness_of_infecting_virus(CD4_index, rf_exp)
            fitness.append(fit)
        return fitness
    

    def singly_infect_cd4(self, next_singly_infecting_viruses):
        # Create a new list with the newly infected cells
        newly_infected = [None]*len(next_singly_infecting_viruses)

        n_added = 0
        # Count up the number of times each virus infects a cell in the next generation
        cd4_norecomb_counts = Counter(next_singly_infecting_viruses)
        # Create next generation of signly infected viruses
        for newly_infected_cd4 in cd4_norecomb_counts:
            # Assign the first instance to the original instance
            newly_infected[n_added] = self.C[newly_infected_cd4]
            n_added += 1
            # Deep copy other instances (if virus infects >1 cell)
            n_cells_to_infect = cd4_norecomb_counts[newly_infected_cd4]
            if n_cells_to_infect > 1:
                for i in range(n_cells_to_infect-1):
                    newly_infected[n_added] = deepcopy(self.C[newly_infected_cd4])
                    n_added += 1

        # Return singly infected active cells
        return newly_infected

    def dually_infect_cd4(self, next_dually_infecting_viruses, cell_breakpoints, seq_len, 
            cost_per_mutation_in_conserved_site, reference_sequence, conserved_sites,
            conserved_fitness, replicative_fitness, rf_exp):
        # Create a new list with the newly infected cells
        newly_infected = [None]*int(len(next_dually_infecting_viruses)/2)

        n_added = 0
        newly_infected_index = 0
        # Create next generation of dually infected viruses
        for breakpoints in cell_breakpoints:
            # Get index of cd4 with the infecting virus
            index_of_cd4_with_virus1 = next_dually_infecting_viruses[newly_infected_index]
            newly_infected_index += 1
            # Get second infecting virus
            index_of_cd4_with_virus2 = next_dually_infecting_viruses[newly_infected_index]
            newly_infected_index += 1

            # Append the cd4 with the infecting virus to newly_infected
            # Deepcopy all of these since we don't know if any of them also singly-infected a cell
            newly_infected[n_added] = deepcopy(self.C[index_of_cd4_with_virus1])

            # Update viral sequence of newly appended cd4 with recombined sequence
            newly_infected[n_added].infecting_virus.nuc_sequence, breakpoints  = \
                get_recombined_sequence(
                    self.C[index_of_cd4_with_virus1].infecting_virus.nuc_sequence,
                    self.C[index_of_cd4_with_virus2].infecting_virus.nuc_sequence,
                    breakpoints)

            if conserved_fitness:
                # Get mutations in conserved sites for recombined virus
                newly_infected[n_added].infecting_virus.conserved_sites_mutated = get_conserved_sites_mutated(
                        self.C[index_of_cd4_with_virus1].infecting_virus.conserved_sites_mutated,
                        self.C[index_of_cd4_with_virus2].infecting_virus.conserved_sites_mutated,
                        breakpoints, seq_len)
                # Update conserved sites fitness for recombined virus
                newly_infected[n_added].infecting_virus.conserved_fitness_cost = conserved_fitness_cost(
                        len(newly_infected[n_added].infecting_virus.conserved_sites_mutated),
                        cost_per_mutation_in_conserved_site)

            # Update replicative fitness for recombined virus
            if replicative_fitness:
                newly_infected[n_added].infecting_virus.replicative_fitness_cost = replicative_fitness_cost(
                        newly_infected[n_added].infecting_virus.nuc_sequence, reference_sequence, conserved_sites, rf_exp)

            n_added +=1

        # Return dually infected active cells
        return newly_infected


    def record_viral_sequences(self, active_filename, latent_filename, n_to_samp):
        with open(active_filename, 'w') as active_f:
            c_sub = sample(self.C, n_to_samp)
            for index, CD4 in enumerate(c_sub):
                name = str(index)
                name += "_ic_" + str("%.4f" % CD4.infecting_virus.immune_fitness_cost)
                name += "_cc_" + str("%.4f" % CD4.infecting_virus.conserved_fitness_cost)
                name += "_rc_" + str("%.4f" % CD4.infecting_virus.replicative_fitness_cost)
                name += "_f_" + str("%.4f" % CD4.infecting_virus.fitness)
                active_f.write(">" + name + "\n")
                active_f.write(CD4.infecting_virus.nuc_sequence + "\n")

        # with open(latent_filename, 'w') as latent_f:
        #    for index, CD4 in enumerate(self.L):
        #        name = str(index)
        #        name += "_ic_" + str("%.4f" % CD4.infecting_virus.immune_fitness_cost)
        #        name += "_cc_" + str("%.4f" % CD4.infecting_virus.conserved_fitness_cost)
        #        name += "_rc_" + str("%.4f" % CD4.infecting_virus.replicative_fitness_cost)
        #        name += "_f_" + str("%.4f" % CD4.infecting_virus.fitness)
        #        latent_f.write(">" + name + "\n")
        #        latent_f.write(CD4.infecting_virus.nuc_sequence + "\n")

    def proliferate_latent_CD4(self, index_to_proliferate):
        assert index_to_proliferate < len(self.L), "Index %s out of range for latent CD4 T cells" % index_to_proliferate
        self.L.append(deepcopy(self.L[index_to_proliferate]))

    def die_latent_CD4(self, index_to_die):
        assert index_to_die < len(self.L), "Index %s out of range for latent CD4 T cells" % index_to_die
        self.L.pop(index_to_die)

    def make_latent(self, index_to_make_latent):
        assert index_to_make_latent < len(self.C), "Index %s out of range for active CD4 T cells" % index_to_make_latent
        self.C[index_to_make_latent].become_latent()
        self.L.append(self.C.pop(index_to_make_latent))

    def make_active(self, index_to_make_active):
        assert index_to_make_active < len(self.L), "Index %s out of range for latent CD4 T cells" % index_to_make_active
        self.L[index_to_make_active].become_active()
        self.C.append(self.L.pop(index_to_make_active))

    def mutate_virus_in_productive_CD4(self, positions_to_mutate, viral_sequence_length, 
            nucleotides_order, substitution_probabilities, conserved_sites, cost_per_mutation_in_conserved_site, reference_sequence,
            conserved_fitness, replicative_fitness, rf_exp):
        # Positions to mutate are given assuming all viral sequences are concatenated
        # We need to identify the right cell number and position within each viral sequence, and initiate the mutation
        for pos in positions_to_mutate:
            cell_number, viral_seq_position = divmod(pos, viral_sequence_length)  # Find cell and position to mutate
            self.C[cell_number].infecting_virus.mutate(viral_seq_position, 
                    nucleotides_order, substitution_probabilities, 
                    conserved_sites, cost_per_mutation_in_conserved_site, reference_sequence,
                    conserved_fitness, replicative_fitness, rf_exp)  # Initiate mutation

    def get_next_gen_latent(self, prob_active_to_latent, prob_latent_to_active, prob_latent_die, prob_latent_proliferate, seed):
        # ***************************** Latent reservoir dynamics ***************************** 
        # Get random number generator
        rng = default_rng(seed)
        
        # Get number of active and latent cells
        num_active_cd4 = len(self.C)
        n_latent_cd4 = len(self.L)
        
        # Move some productively infected cells to the latent pool
        num_to_make_latent = min(rng.binomial(num_active_cd4, prob_active_to_latent), num_active_cd4-1) # make it so not all cells can become latent
        indices_to_make_latent = sample(range(num_active_cd4), num_to_make_latent)
        for i in range(num_to_make_latent):
            index_to_make_latent = indices_to_make_latent[i] - i
            self.make_latent(index_to_make_latent)

        # # Draw the number of cells that leave the latent reservoir by activation or death, or proliferate
        # prob_event = prob_latent_to_active + prob_latent_die + prob_latent_proliferate
        # num_latent_events = rng.binomial(n_latent_cd4, prob_event)
        # num_to_activate, num_to_die, num_to_proliferate = rng.multinomial(num_latent_events,
        #                                                     [prob_latent_to_active / prob_event,
        #                                                      prob_latent_die / prob_event,
        #                                                      prob_latent_proliferate / prob_event])
        # assert num_to_activate + num_to_die + num_to_proliferate == num_latent_events, "Activation, death, and proliferation" \
        #         "don't add up to number of cells leaving latent pool"
        # 
        # # Get indices of cells that activate, die, or proliferate
        # latent_event_indices = sample(range(n_latent_cd4), num_latent_events)
        # to_active = latent_event_indices[:num_to_activate]
        # to_die = latent_event_indices[num_to_activate:num_to_activate+num_to_die]
        # to_proliferate = latent_event_indices[num_to_activate+num_to_die:]
        
        latent_indices = range(n_latent_cd4)
        to_active = where(rng.binomial(1, prob_latent_to_active, n_latent_cd4))[0].tolist()
        to_die = where(rng.binomial(1, prob_latent_die, n_latent_cd4))[0].tolist()
        to_proliferate = where(rng.binomial(1, prob_latent_proliferate, n_latent_cd4))[0].tolist()
        latent_event_indices = set(to_active + to_die + to_proliferate)

        # Perform events
        # Sort so that the highest indices are popped first (doesn't influence earlier indices)
        n_to_active = 0
        n_to_die = 0
        n_to_proliferate = 0
        for i in sorted(latent_event_indices, reverse=True):
            if i in to_active:
                self.make_active(i)
                n_to_active += 1
            elif i in to_die:
                self.die_latent_CD4(i)
                n_to_die += 1
            else:
                self.proliferate_latent_CD4(i)
                n_to_proliferate += 1

        return len(latent_event_indices), n_to_active, n_to_die, n_to_proliferate

    def get_next_gen_active(self, prob_mutation, prob_recombination, n_active_next_gen, gen, seroconversion_time, 
            nucleotides_order, substitution_probabilities, conserved_sites, time_to_full_potency, cost_per_mutation_in_conserved_site, 
            reference_sequence, immune_response_proportion, epitopes, immune_fitness, conserved_fitness, replicative_fitness, rf_exp, seed):
        # ***************************** Productively infected cell dynamics ***************************** #
        # Get random number generator
        rng = default_rng(seed)
        # Mutate virus in productively infected cells
        # Identify number of mutations
        num_active_cd4 = len(self.C)
        seq_len = len(self.C[0].infecting_virus.nuc_sequence)
        total_nucleotides = num_active_cd4 * seq_len
        n_mut = rng.binomial(total_nucleotides, prob_mutation)
        # Distribute the mutations across all possible positions
        positions_to_mutate = sample(range(total_nucleotides), n_mut)
        self.mutate_virus_in_productive_CD4(positions_to_mutate, seq_len, nucleotides_order, substitution_probabilities, 
                conserved_sites, cost_per_mutation_in_conserved_site, reference_sequence, conserved_fitness, replicative_fitness, rf_exp)

        if immune_fitness:
            # Update immune response
            if gen >= seroconversion_time:
                # Epitope recognition
                self.update_epitopes_recognized(gen, epitopes, immune_response_proportion, time_to_full_potency, seed)
                # Virus immune fitness
                self.update_immune_fitness(epitopes, gen, time_to_full_potency)

        # Compute fitness
        fitness = self.get_fitness(rf_exp)

        # Determine number of dual infections with recombination
        num_cells_recomb, cell_breakpoints = get_recomb_breakpoints(seq_len, n_active_next_gen, prob_recombination, rng)
        
        # Sample num_active_cd4 + num_cells_recomb virus variants based on relative fitness
        # We normalize the fitness values. The next generation of CD4s will be infected by virions according to this
        # normalized fitness.
        next_infecting_viruses = choices(range(len(fitness)), normalize(fitness),
                                                  k=(int(n_active_next_gen + num_cells_recomb)))

        # Infect n_active_next_gen - num_cells_recomb cells with a single virus. Here we will be reusing the same data structure (host.C)
        # but replacing the virus. First n_active_next_gen - num_cells_recomb indices from the sampled next_infecting_virus will be used.
        self.C = (self.singly_infect_cd4(next_infecting_viruses[:int(n_active_next_gen-num_cells_recomb)]) +
        # Infect num_cells_recomb with two viruses and create recombinants
        self.dually_infect_cd4(next_infecting_viruses[int(n_active_next_gen-num_cells_recomb):], cell_breakpoints, seq_len, 
                        cost_per_mutation_in_conserved_site, reference_sequence, conserved_sites, conserved_fitness, replicative_fitness, rf_exp))

        return(n_mut, num_cells_recomb)

    
    def summarize_fitness(self):
        n_active = len(self.C)
        return sum([x.infecting_virus.fitness for x in self.C])/n_active, \
                sum([x.infecting_virus.conserved_fitness_cost for x in self.C])/n_active, \
                sum([x.infecting_virus.immune_fitness_cost for x in self.C])/n_active, \
                sum([x.infecting_virus.replicative_fitness_cost for x in self.C])/n_active # MAYBE UPDATE THIS IF WE KEEP REPLICATIVE FITNESS AS (1-RFC)^EXP
