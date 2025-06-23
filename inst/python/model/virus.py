from model.fitness import calc_seq_fitness
from model.variation import muts_rel_ref, get_substitution

class HIV:
    def __init__(self, nuc_seq, reference_sequence, replicative_cost):
        # Make sure values supplied are as expected
        assert isinstance(
            nuc_seq, str), "Nucleotide sequence needs to be a string"

        # Initialize instance variables
        self.nuc_sequence = nuc_seq
        self.conserved_sites_mutated = set()

        # Tracking different components of fitness
        self.immune_fitness = 1
        self.conserved_fitness = 1
        self.replicative_fitness = 1
        self.fitness = 1
        if len(reference_sequence):
            self.replicative_fitness = calc_seq_fitness(muts_rel_ref(
                self.nuc_sequence, reference_sequence), replicative_cost)
            self.fitness = self.replicative_fitness
       

    def __repr__(self):
        return_str = "HIV with sequence %s" % self.nuc_sequence
        return return_str

    def mutate(
        self,
        position_to_mutate,
        #nucleotides_order,
        substitution_probabilities,
        conserved_sites,
        cost_per_mutation_in_conserved_site,
        reference_sequence,
        replicative_cost,
    ):
        # Figure out what the mutation is based on mutation probabilities
        old_nucleotide = self.nuc_sequence[position_to_mutate]
        new_nucleotide = get_substitution(
            old_nucleotide, #nucleotides_order, 
            substitution_probabilities
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
            if new_nucleotide == conserved_sites[position_to_mutate]:
                self.conserved_sites_mutated.remove(position_to_mutate)
                self.conserved_fitness = calc_seq_fitness(
                    len(self.conserved_sites_mutated), cost_per_mutation_in_conserved_site
                )
        # If mutation is in a conserved site, update mutations in conserved sites and conserved sites fitness
        elif position_to_mutate in conserved_sites.keys():
            self.conserved_sites_mutated.add(position_to_mutate)
            self.conserved_fitness = calc_seq_fitness(
                len(self.conserved_sites_mutated), cost_per_mutation_in_conserved_site
            )

        # Update replicative fitness cost only if needed
        if len(reference_sequence):
            ref_base = reference_sequence[position_to_mutate]
            prev_comp = old_nucleotide == ref_base
            if prev_comp != (ref_base == new_nucleotide):
                self.replicative_fitness = calc_seq_fitness(muts_rel_ref(
                    self.nuc_sequence, reference_sequence), replicative_cost)
