from random import seed
from numpy.random import default_rng

def set_python_seed(s):
    if s is not None:
        # Set seed (need to set a seed for random and for numpy.random)
        seed(s)
        generator = default_rng(int(s))
    else:
        generator = default_rng()
    return generator
