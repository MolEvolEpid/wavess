from numpy.random import default_rng
from numpy import where
from collections import Counter
from copy import deepcopy
from collections import namedtuple

from model.virus import HIV
from model.cells import InfectedCD4
from model.fitness import *
from model.variation import *

class HostEnv:  # This is the 'compartment' where the model dynamics take place
    def __init__(self, founder_viruses: HIV, NC: int):

        assert NC == len(
            founder_viruses
        ), "Initial population size must equal the number of founder sequences"

        # Create NC actively proliferating infected cells
        self.C = [InfectedCD4(founder, True) for founder in founder_viruses]

        # Initiate latent infected cells
        self.L = list()

        # Store epitope translations
        self.epitope_variants_translated = defaultdict(lambda: "")
        # Track epitopes recognized by the immune system
        # key -> epitope sequence variant, value -> generation when it started
        # being recognized
        self.epitopes_recognition_generation = defaultdict(lambda: 0)
        # Immune cost of epitopes
        self.cross_reactive_epitope_cost_frac = defaultdict(lambda: 0)

    def __repr__(self):
        return_str = "Host has %s active and %s latent infected cells\n" % (
            str(len(self.C)),
            str(len(self.L)),
        )
        return_str += "%s epitopes recognized" % len(
            self.epitopes_recognition_generation
        )
        return return_str
      

    def update_b_epitopes_recognized(self, current_generation, config):
        # Get random number generator
        rng = default_rng(config.generator)
        # Immune strengths
        immune_strength_dict = {
            k: min((current_generation - v) / config.b_gen_full_potency, 1)
            for k, v in self.epitopes_recognition_generation.items()
        }
        # Number of active infected cells
        num_active_cd4 = len(self.C)
        # Get unique sequences
        seqs = defaultdict(lambda: 0)
        for index in range(num_active_cd4):
            seqs[self.C[index].infecting_virus.nuc_sequence] += 1
        # Track frequency of epitope sequence variants in sample
        epitope_variants = defaultdict(lambda: 0)
        for seq, freq in seqs.items():
            for epi in config.b_epitope_locations:
                epitope_variants[seq[epi.start: epi.end]] += freq

        # Translate epitopes, get frequencies, and store translations
        aa_epitope_variants = defaultdict(lambda: 0)
        for variant, frequency in epitope_variants.items():
            # Translate and add to dictionary if new variant
            if variant not in self.epitope_variants_translated:
                self.epitope_variants_translated[variant] = translate(
                    variant)
            # Add variant frequency to dictionary
            aa_epitope_variants[self.epitope_variants_translated[variant]] += frequency

        # Get cross-reactive strength for each new epitope
        for variant, frequency in aa_epitope_variants.items():
            if (
                variant not in self.epitopes_recognition_generation
                and len(self.epitopes_recognition_generation) != 0
            ):
                if variant not in self.cross_reactive_epitope_cost_frac:
                    most_similar = ""
                    n_diff_min = len(variant)
                    for v in self.epitopes_recognition_generation:
                        if variant != v and len(variant) == len(v):
                            n_diff = sum(1 for x, y in zip(
                                variant, v) if x != y)
                            if n_diff <= n_diff_min:
                                n_diff_min = n_diff
                                most_similar = v
                    self.cross_reactive_epitope_cost_frac[
                        variant
                    ] = immune_strength_dict[most_similar] * rng.beta(
                        1, n_diff_min**2
                    )
                    
        # Mark frequent epitopes as recognized
        for variant, frequency in aa_epitope_variants.items():
            if frequency >= config.b_n_for_imm:  # common enough to be recognized
                if (
                    variant not in self.epitopes_recognition_generation
                ):  # if first time epitope is being recognized
                    if variant not in self.cross_reactive_epitope_cost_frac:
                        self.epitopes_recognition_generation[variant] = (
                            current_generation
                        )
                    else:
                        self.epitopes_recognition_generation[variant] = int(
                            current_generation
                            - self.cross_reactive_epitope_cost_frac[variant]
                            * config.b_gen_full_potency
                        )

    def update_b_immune_fitness(self, current_generation, config):
        immune_strength_dict = {
            k: min((current_generation - v) / config.b_gen_full_potency, 1)
            for k, v in self.epitopes_recognition_generation.items()
        }
        for CD4_index in range(len(self.C)):
            # Fitness cost due to immune response
            max_epitope_fitness_cost = 0  # by default the virus is super fit
            for epi in config.b_epitope_locations:
                # get translated epitope
                epitope_sequence = self.epitope_variants_translated[
                    self.C[CD4_index].infecting_virus.nuc_sequence[epi.start: epi.end]
                ]
                if epitope_sequence in self.epitopes_recognition_generation:
                    timed_fitness = (
                        epi.max_fitness *
                        immune_strength_dict[epitope_sequence]
                    )
                    max_epitope_fitness_cost = max(
                        max_epitope_fitness_cost, timed_fitness
                    )
                elif epitope_sequence in self.cross_reactive_epitope_cost_frac:
                    timed_fitness = (
                        epi.max_fitness
                        * self.cross_reactive_epitope_cost_frac[epitope_sequence]
                    )
                    max_epitope_fitness_cost = max(
                        max_epitope_fitness_cost, timed_fitness
                    )
            self.C[CD4_index].infecting_virus.b_immune_fitness = 1 - \
                max_epitope_fitness_cost

    
    def get_fitness_of_infecting_virus(self, CD4_index):
        
        # Compute overall fitness
        viral_fitness = (
            self.C[CD4_index].infecting_virus.b_immune_fitness
            * self.C[CD4_index].infecting_virus.t_immune_fitness
            * self.C[CD4_index].infecting_virus.conserved_fitness
            * self.C[CD4_index].infecting_virus.replicative_fitness
        )

        # Update the virus object with the fitness costs
        self.C[CD4_index].infecting_virus.fitness = viral_fitness

        return viral_fitness

    def get_fitness(self):
        fitness = []
        for CD4_index in range(len(self.C)):
            fit = self.get_fitness_of_infecting_virus(CD4_index)
            fitness.append(fit)
        return fitness

    def singly_infect_cd4(self, next_singly_infecting_viruses):
        # Create a new list with the newly infected cells
        newly_infected = [None] * len(next_singly_infecting_viruses)

        n_added = 0
        # Count up the number of times each virus infects a cell in the next
        # generation
        cd4_norecomb_counts = Counter(next_singly_infecting_viruses)
        # Create next generation of signly infected viruses
        for newly_infected_cd4 in cd4_norecomb_counts:
            # Assign the first instance to the original instance
            newly_infected[n_added] = self.C[newly_infected_cd4]
            n_added += 1
            # Deep copy other instances (if virus infects >1 cell)
            n_cells_to_infect = cd4_norecomb_counts[newly_infected_cd4]
            if n_cells_to_infect > 1:
                for i in range(n_cells_to_infect - 1):
                    newly_infected[n_added] = deepcopy(
                        self.C[newly_infected_cd4])
                    n_added += 1

        # Return singly infected active cells
        return newly_infected

    def dually_infect_cd4(self, next_dually_infecting_viruses, cell_breakpoints, config):
        # Create a new list with the newly infected cells
        newly_infected = [None] * int(len(next_dually_infecting_viruses) / 2)

        n_added = 0
        newly_infected_index = 0
        # Create next generation of dually infected viruses
        for breakpoints in cell_breakpoints:
            # Get index of cd4 with the infecting virus
            index_of_cd4_with_virus1 = next_dually_infecting_viruses[
                newly_infected_index
            ]
            newly_infected_index += 1
            # Get second infecting virus
            index_of_cd4_with_virus2 = next_dually_infecting_viruses[
                newly_infected_index
            ]
            newly_infected_index += 1

            # Append the cd4 with the infecting virus to newly_infected
            # Deepcopy all of these since we don't know if any of them also
            # singly-infected a cell
            newly_infected[n_added] = deepcopy(
                self.C[index_of_cd4_with_virus1])

            # Update viral sequence of newly appended cd4 with recombined
            # sequence
            newly_infected[n_added].infecting_virus.nuc_sequence, breakpoints = (
                get_recombined_sequence(
                    self.C[index_of_cd4_with_virus1].infecting_virus.nuc_sequence,
                    self.C[index_of_cd4_with_virus2].infecting_virus.nuc_sequence,
                    breakpoints
                )
            )

            if len(config.conserved_sites):
                # Get mutations in conserved sites for recombined virus
                newly_infected[n_added].infecting_virus.conserved_sites_mutated = (
                    get_conserved_sites_mutated(
                        self.C[
                            index_of_cd4_with_virus1
                        ].infecting_virus.conserved_sites_mutated,
                        self.C[
                            index_of_cd4_with_virus2
                        ].infecting_virus.conserved_sites_mutated,
                        breakpoints, config.seq_len
                    )
                )
                # Update conserved sites fitness for recombined virus
                newly_infected[n_added].infecting_virus.conserved_fitness = (
                    calc_seq_fitness(
                        len(
                            newly_infected[
                                n_added
                            ].infecting_virus.conserved_sites_mutated
                        ),
                        config.conserved_cost
                    )
                )

            # Update replicative fitness for recombined virus
            if len(config.ref_seq):
                newly_infected[n_added].infecting_virus.replicative_fitness = (
                    calc_seq_fitness(muts_rel_ref(
                        newly_infected[n_added].infecting_virus.nuc_sequence, config.ref_seq), config.replicative_cost)
                )

            n_added += 1

        # Return dually infected active cells
        return newly_infected

    def proliferate_latent_CD4(self, index_to_proliferate):
        assert index_to_proliferate < len(self.L), (
            "Index %s out of range for latent CD4 T cells" % index_to_proliferate
        )
        self.L.append(deepcopy(self.L[index_to_proliferate]))

    def die_latent_CD4(self, index_to_die):
        assert index_to_die < len(self.L), (
            "Index %s out of range for latent CD4 T cells" % index_to_die
        )
        self.L.pop(index_to_die)

    def make_latent(self, index_to_make_latent):
        assert index_to_make_latent < len(self.C), (
            "Index %s out of range for active CD4 T cells" % index_to_make_latent
        )
        self.C[index_to_make_latent].become_latent()
        self.L.append(self.C.pop(index_to_make_latent))

    def make_active(self, index_to_make_active):
        assert index_to_make_active < len(self.L), (
            "Index %s out of range for latent CD4 T cells" % index_to_make_active
        )
        self.L[index_to_make_active].become_active()
        self.C.append(self.L.pop(index_to_make_active))

    def mutate_virus_in_productive_CD4(self, positions_to_mutate, config):
        # Positions to mutate are given assuming all viral sequences are concatenated
        # We need to identify the right cell number and position within each
        # viral sequence, and initiate the mutation
        for pos in positions_to_mutate:
            cell_number, viral_seq_position = divmod(
                pos, config.seq_len
            )  # Find cell and position to mutate
            self.C[cell_number].infecting_virus.mutate(viral_seq_position, config)  # Initiate mutation

    def get_next_gen_latent(self, config):
        # ***************************** Latent reservoir dynamics *************
        # Get random number generator
        rng = default_rng(config.generator)

        # Get number of active and latent cells
        num_active_cd4 = len(self.C)
        n_latent_cd4 = len(self.L)

        # Move some productively infected cells to the latent pool
        num_to_make_latent = min(
            rng.binomial(
                num_active_cd4, config.prob_act_to_lat), num_active_cd4 - 1
        )  # make it so not all cells can become latent
        indices_to_make_latent = sample(
            range(num_active_cd4), num_to_make_latent)
        for i in range(num_to_make_latent):
            index_to_make_latent = indices_to_make_latent[i] - i
            self.make_latent(index_to_make_latent)

        latent_indices = range(n_latent_cd4)
        to_active = where(rng.binomial(1, config.prob_lat_to_act, n_latent_cd4))[
            0
        ].tolist()
        to_die = where(
            rng.binomial(1, config.prob_lat_die, n_latent_cd4))[0].tolist()
        to_proliferate = where(rng.binomial(1, config.prob_lat_prolif, n_latent_cd4))[
            0
        ].tolist()
        latent_event_indices = set(to_active + to_die + to_proliferate)

        # Perform events
        # Sort so that the highest indices are popped first (doesn't influence
        # earlier indices)
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

        return num_to_make_latent, n_to_active, n_to_die, n_to_proliferate

    def get_next_gen_active(self, n_active_next_gen, gen, config):
        # ***************************** Productively infected cell dynamics ***************************** #
        # Get random number generator
        rng = default_rng(config.generator)
        # Mutate virus in productively infected cells
        # Identify number of mutations
        num_active_cd4 = len(self.C)
        #seq_len = len(self.C[0].infecting_virus.nuc_sequence)
        total_nucleotides = num_active_cd4 * config.seq_len
        n_mut = rng.binomial(total_nucleotides, config.prob_mut)
        # Distribute the mutations across all possible positions
        #positions_to_mutate = sample(range(total_nucleotides), n_mut)
        if n_mut > 0:
            positions_to_mutate = rng.choice(total_nucleotides, n_mut, replace=False)
        else:
            positions_to_mutate = []
        self.mutate_virus_in_productive_CD4(positions_to_mutate, config)
        
        if config.b_epitope_locations is not None:
            # Update immune response
            if gen >= config.b_immune_start_gen:
                # Epitope recognition
                self.update_b_epitopes_recognized(gen, config)
                # Virus immune fitness
                self.update_b_immune_fitness(gen, config)

        # t immune fitness depends on generation so have to update every time
        if config.t_epitope_locations is not None:
            for CD4_index in range(len(self.C)):
                # THIS IS A LOT SLOWER NOW BECAUSE EACH EPITOPE CAN HAVE A DIFFERENT TIME TO MAX POTENCY
                # SO CAN'T JUST UPDATE IT WHEN THERE'S A MUTATION IN AN EPITOPE
                # CAN PROBABLY MAKE IT FASTER BECAUSE AT SOME POINT THEY'LL ALL BE THE SAME
                self.C[CD4_index].infecting_virus.determine_t_recognition(gen, config)
            # self.C[CD4_index].infecting_virus.recognized_by_t
            # if gen >= config.t_gen_full_potency:
            #     self.C[CD4_index].infecting_virus.t_immune_fitness = (1 - config.t_max_imm) ** self.C[CD4_index].infecting_virus.recognized_by_t
            # else:
            #     self.C[CD4_index].infecting_virus.t_immune_fitness = (1 - config.t_max_imm * gen / config.t_gen_full_potency) ** self.C[CD4_index].infecting_virus.recognized_by_t

        # WHY IS IT NOT EXACTLY 2000? BECAUSE OF LATENT CELLS?
        # print(gen)
        # print(Counter([self.C[CD4_index].infecting_virus.t_immune_fitness for CD4_index in range(len(self.C))]))

        # Compute fitness
        fitness = self.get_fitness()

        # Determine number of dual infections with recombination
        num_cells_recomb, cell_breakpoints = get_recomb_breakpoints(int(n_active_next_gen), config)

        # Sample num_active_cd4 + num_cells_recomb virus variants based on relative fitness
        # We normalize the fitness values. The next generation of CD4s will be infected by virions according to this
        # normalized fitness.
        next_infecting_viruses = choices(
            range(len(fitness)),
            normalize(fitness),
            k=(int(n_active_next_gen + num_cells_recomb)),
        )

        # Infect n_active_next_gen - num_cells_recomb cells with a single virus. 
        # Here we will be reusing the same data structure (host.C)
        # but replacing the virus. First n_active_next_gen - num_cells_recomb
        # indices from the sampled next_infecting_virus will be used.
        self.C = (
            self.singly_infect_cd4(
                next_infecting_viruses[: int(
                    n_active_next_gen - num_cells_recomb)]
            )
            +
            # Infect num_cells_recomb with two viruses and create recombinants
            self.dually_infect_cd4(
                next_infecting_viruses[int(n_active_next_gen - num_cells_recomb):],
                cell_breakpoints,
                config
            )
        )

        return (n_mut, num_cells_recomb)

    def summarize_fitness(self):
        n_active = len(self.C)
        return (
            sum([x.infecting_virus.fitness for x in self.C]) / n_active,
            sum([x.infecting_virus.conserved_fitness for x in self.C]) / n_active,
            sum([x.infecting_virus.b_immune_fitness for x in self.C]) / n_active,
            sum([x.infecting_virus.t_immune_fitness for x in self.C]) / n_active,
            sum([x.infecting_virus.replicative_fitness for x in self.C])
            / n_active,
        )

    def record_counts(self, counts, generation,
                      latent_nums, var_nums, fitness):
        counts["generation"].append(generation)
        counts["active_cell_count"].append(len(self.C))
        counts["latent_cell_count"].append(len(self.L))
        # num_to_make_latent, num_to_activate, num_to_die, num_to_proliferate
        counts["active_turned_latent"].append(latent_nums[0]),
        counts["latent_turned_active"].append(latent_nums[1]),
        counts["latent_died"].append(latent_nums[2]),
        counts["latent_proliferated"].append(latent_nums[3]),
        # n_mut, number_recombinations
        counts["number_mutations"].append(var_nums[0]),
        counts["number_recombinations"].append(var_nums[1]),
        # mean_fitness_active, mean_conserved_active,
        # mean_immune_active, mean_replicative_active
        counts["mean_fitness_active"].append(fitness[0]),
        counts["mean_conserved_active"].append(fitness[1]),
        counts["mean_b_immune_active"].append(fitness[2]),
        counts["mean_t_immune_active"].append(fitness[3]),
        counts["mean_replicative_active"].append(fitness[4])
        return counts

    def sample_viral_sequences(self, seqs, fitness, generation, n_to_samp, cell_type = 'active'):
        assert cell_type in ['active', 'latent'], "sampling cell type must be active or latent"
        c_sub = sample(self.C, int(n_to_samp))
        if cell_type == 'active':
            c_sub = sample(self.C, int(n_to_samp))
        elif cell_type == 'latent':
            c_sub = sample(self.L, int(n_to_samp))
        for index, CD4 in enumerate(c_sub):
            name = "gen" + str(generation)
            name += "_" + cell_type
            name += "_" + str(index)
            seqs[name] = CD4.infecting_virus.nuc_sequence
            if cell_type == 'active':
                fitness["generation"].append(str(generation))
                fitness["seq_id"].append(name)
                fitness["b_immune"].append(float(CD4.infecting_virus.b_immune_fitness))
                fitness["t_immune"].append(float(CD4.infecting_virus.t_immune_fitness))
                fitness["conserved"].append(
                    float(CD4.infecting_virus.conserved_fitness))
                fitness["replicative"].append(
                    float(CD4.infecting_virus.replicative_fitness))
                fitness["overall"].append(float(CD4.infecting_virus.fitness))
        return seqs, fitness

    def loop_through_generations(self, active_cell_count, n_sample_active, n_sample_latent, config):
        
        # Initialize counts, fitness, and sequence objects
        counts = {
            "generation": [],
            "active_cell_count": [],
            "latent_cell_count": [],
            "active_turned_latent": [],
            "latent_turned_active": [],
            "latent_died": [],
            "latent_proliferated": [],
            "number_mutations": [],
            "number_recombinations": [],
            "mean_fitness_active": [],
            "mean_conserved_active": [],
            "mean_b_immune_active": [],
            "mean_t_immune_active": [],
            "mean_replicative_active": [],
        }

        fitness = {
            "generation": [],
            "seq_id": [],
            "b_immune": [],
            "t_immune": [],
            "conserved": [],
            "replicative": [],
            "overall": []
        }

        # put founders at top of fasta and fitness
        seqs_active = config.founder_seqs
        for fname, fseq in seqs_active.items():
            fitness["generation"].append('founder')
            fitness["seq_id"].append(fname)
            fitness["b_immune"].append(float(1))
            fitness["t_immune"].append(float(1)) # UPDATE!!!
            fitness["conserved"].append(float(1))
            if len(config.ref_seq):
                repfit = calc_seq_fitness(muts_rel_ref(
                    fseq, config.ref_seq), config.replicative_cost)
                fitness["replicative"].append(repfit)
                fitness["overall"].append(repfit)
            else:
                fitness["replicative"].append(float(1))
                fitness["overall"].append(float(1))
                
        seqs_latent = defaultdict(lambda: 0)

        if n_sample_active[0] != 0:
            # num_to_make_latent, num_to_activate, num_to_die, num_to_proliferate
            latent_nums = [0, 0, 0, 0]
            # n_mut, number_recombinations
            var_nums = [0, 0]
            # mean_fitness_active, mean_conserved_active, 
            # mean_b_immune_active, mean_t_immune_active
            # mean_replicative_active
            mean_fitness = self.summarize_fitness()
            counts = self.record_counts(
                counts, 0, latent_nums, var_nums, mean_fitness)
            seqs_active, fitness = self.sample_viral_sequences(
                seqs_active, fitness, 0, min(int(n_sample_active[0]), len(seqs_active)))

        # Looping through generations until we sample everything we want
        for t in range(1, int(config.last_sampled_gen) + 1):
            # Only get latent reservoir dynamics if modeling
            latent_nums = [0, 0, 0, 0]
            if config.prob_act_to_lat:
                # num_to_make_latent, num_to_activate, num_to_die, num_to_proliferate
                latent_nums = self.get_next_gen_latent(config)
            # Productively infected cell dynamics
            # n_mut, number_recombinations
            var_nums = self.get_next_gen_active(active_cell_count[t], t, config)
            # Record events
            if n_sample_active[t] != 0:
                mean_fitness = self.summarize_fitness()
                counts = self.record_counts(
                    counts, t, latent_nums, var_nums, mean_fitness)
                seqs_active, fitness = self.sample_viral_sequences(
                    seqs_active, fitness, t, min(int(n_sample_active[t]), active_cell_count[t]))
                    
            if n_sample_latent[t] != 0:
                seqs_latent, tmp = self.sample_viral_sequences(
                    seqs_latent, None, t, min(int(n_sample_latent[t]), len(self.L)), 'latent')

        return counts, fitness, seqs_active, dict(seqs_latent)
