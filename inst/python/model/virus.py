from model.fitness import calc_seq_fitness
from model.variation import muts_rel_ref, get_substitution, translate

class HIV:
    def __init__(self, nuc_seq, config):
        # Make sure values supplied are as expected
        assert isinstance(
            nuc_seq, str), "Nucleotide sequence needs to be a string"

        # Initialize instance variables
        self.nuc_sequence = nuc_seq
        self.conserved_sites_mutated = set()

        # Tracking different components of fitness
        self.b_immune_fitness = 1
        self.t_immune_fitness = 1 # SHOULD THIS START AT 1 IN GEN ZERO OR ALREADY LESS? COST WILL INCREASE QUICKLY REGARDLESS
        self.recognized_by_t = False
        if len(config.t_epitope_locations):
            self.determine_t_recognition(config)
        self.conserved_fitness = 1
        self.replicative_fitness = 1
        self.fitness = 1
        if len(config.ref_seq):
            self.replicative_fitness = calc_seq_fitness(muts_rel_ref(
                self.nuc_sequence, config.ref_seq), config.replicative_cost)
            self.fitness = self.replicative_fitness
            

    def __repr__(self):
        return_str = "HIV with sequence %s" % self.nuc_sequence
        return return_str
      
    # REMEMBER TO MAKE CONSERVED SITES AT EPITOPE LOCATIONS NOT CONSERVED IN PREP!!
    def determine_t_recognition(self, config):
        n_epitopes = len(config.t_epitope_locations)
        recognized_by_t = 0 #False
        for epi in config.t_epitope_locations:
            virus_epitope = list(translate(self.nuc_sequence[epi.start:epi.end]))
            virus_epitope = ''.join([virus_epitope[pos-1] if pos in epi.escape_positions else 'N' for pos in range(1, max(epi.escape_positions.keys())+1)])
            #print(config.t_recognition_motifs)
            #print(virus_epitope)
            if virus_epitope in config.t_recognition_motifs:
                recognized_by_t += 1 #= True
                #break
                 #t_immune_cost = epi.max_fitness_cost if current_generation >= config.t_gen_full_potency else epi.max_fitness_cost * current_generation / config.t_gen_full_potency
        # WE MIGHT NEED TO MAKE THE POPULATION SIZE BIGGER TO GET A MUTATION IN THE RIGHT SPOT OR IS THERE JUST A BUG IN THE CODE
        # OR MAKE IT SO THAT THE FITNESS COST IS PROPORTIONAL TO THE NUMBER OF EPITOPES RECOGNIZED
        #if not recognized_by_t: print(recognized_by_t)
        # print(recognized_by_t)
        self.recognized_by_t = recognized_by_t#/n_epitopes


    def mutate(self, position_to_mutate, config):
        # Figure out what the mutation is based on mutation probabilities
        old_nucleotide = self.nuc_sequence[position_to_mutate]
        new_nucleotide = get_substitution(
            old_nucleotide, 
            config.substitution_probabilities
        )
        # Fastest way I can find to mutate - add part before to new nt to part
        # after
        self.nuc_sequence = (
            self.nuc_sequence[:position_to_mutate]
            + new_nucleotide
            + self.nuc_sequence[position_to_mutate + 1:]
        )
        # If the mutation is at a mutated conserved site, check to see if it reverted back to conserved base
        if position_to_mutate in self.conserved_sites_mutated:
            if new_nucleotide == config.conserved_sites[position_to_mutate]:
                self.conserved_sites_mutated.remove(position_to_mutate)
                self.conserved_fitness = calc_seq_fitness(
                    len(self.conserved_sites_mutated), config.conserved_cost
                )
        # If mutation is in a conserved site, update mutations in conserved sites and conserved sites fitness
        elif position_to_mutate in config.conserved_sites.keys():
            self.conserved_sites_mutated.add(position_to_mutate)
            self.conserved_fitness = calc_seq_fitness(
                len(self.conserved_sites_mutated), config.conserved_cost
            )
            
        # if mutation is at a t-cell epitope recognition location, check to see if virus is still recognized
        if config.t_epitope_mask[position_to_mutate] >= 0:
            self.determine_t_recognition(config)

        # Update replicative fitness cost only if needed
        if len(config.ref_seq):
            ref_base = config.ref_seq[position_to_mutate]
            prev_comp = old_nucleotide == ref_base
            if prev_comp != (ref_base == new_nucleotide):
                self.replicative_fitness = calc_seq_fitness(muts_rel_ref(
                    self.nuc_sequence, config.ref_seq), config.replicative_cost)
                    
