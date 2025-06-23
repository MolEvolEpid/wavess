def calc_seq_fitness(n_muts, cost):
    return (1-cost)**n_muts
  
  
def normalize(likelihoods):
    total = sum(likelihoods)
    if total == 0:
        raise ValueError("All virus fitnesses are 0")
    norm = [likelihood / total for likelihood in likelihoods]
    return norm
