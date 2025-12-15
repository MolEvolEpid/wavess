import pytest
from random import seed
from numpy.random import default_rng
from random import sample
from collections import defaultdict
import pandas as pd

from run_wavess import *
from model.seed import set_python_seed
from model.epitope import BEpitope, TEpitope
from model.variation import get_substitution, get_recomb_breakpoints, get_recombined_sequence, muts_rel_ref, get_conserved_sites_mutated, translate
from model.fitness import calc_seq_fitness, normalize
from model.virus import HIV
from model.cells import InfectedCD4
from model.host import HostEnv
from model.config import Config


q = read_csv("../extdata/hiv_q_mat.csv", index_col="nt_from")
#pd.DataFrame({'A': {'A': -1.3935601992358562, 'C': 0.9167813093518908, 'G': 2.933888073407108, 'T': 0.5499598145717008}, 'C': {'A': 0.165021502856394, 'C': -3.2089352256667243, 'G': 0.0182889027942411, 'T': 1.8338005458527384}, 'G': {'A': 1.100214466648487, 'C': 0.0915971984019693, 'G': -3.318923580552964, 'T': 0.550113237779834}, 'T': {'A': 0.1283242296110112, 'C': 2.200556717904873, 'G': 0.366746604293706, 'T': -2.9338735982413424}})


# Functions outside of any class

def test_set_python_seed():
    assert isinstance(set_python_seed(None), np.random._generator.Generator)
    g = set_python_seed(1)
    assert sample(range(1000), 1) == [137]
    assert g.binomial(1000, 0.1) == 100
    

# def test_create_epitope():
#     epitope = create_epitope(10, 20, 0.4)
#     assert vars(epitope) == {'start': 10, 'end': 20, 'max_fitness': 0.4}


# Functions in Config

def test_prep_ref_conserved():
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 0, 0,
          0, 0, 0, 0,
          {}, 0, "AAA", 1,
          None, 0, 0, 0,
          None, 0, 0)
    assert (config.ref_seq, config.conserved_sites) == ('AAA', {})
    config.prep_ref_conserved() 
    assert (config.ref_seq, config.conserved_sites) == ('AAA', {})
    config.conserved_sites = {0: 'a'}
    config.prep_ref_conserved() 
    assert (config.ref_seq, config.conserved_sites) == ('-AA', {0: 'a'})
    config.founder_seqs = {'founder0': 'TAA'}
    config.ref_seq = 'AAA'
    config.prep_ref_conserved() 
    assert (config.ref_seq, config.conserved_sites) == ('AAA', {})
    config.founder_seqs = {'founder0': 'TAA', 'founder1': 'TAA'}
    config.conserved_sites = {0: 'a', 2: 'c'}
    config.prep_ref_conserved() 
    assert (config.ref_seq, config.conserved_sites) == ('AA-', {2: 'c'})
    config.founder_seqs = {'founder0': 'TAA', 'founder1': 'AAT'}
    config.ref_seq = 'AAA'
    config.conserved_sites = {0: 'a', 2: 'c'}
    config.prep_ref_conserved() 
    assert (config.ref_seq, config.conserved_sites) == ('AAA', {})


def test_create_host_env():
    config = Config(1, None, {"founder0": "AAA"}, 
        q, 3.5e-5, 0, 
        0, 0, 0, 0,
        {}, None, "AAA", 1,
        None, 0, 0, 0,
        None, 0, 0)
    assert config.host.C[0].active
    assert config.host.C[0].infecting_virus.nuc_sequence == "AAA"
    config = Config(1, None,  {"founder0": "AAA", "founder1": "GGG"}, 
        q, 3.5e-5, 0, 
        0, 0, 0, 0,
        {}, None, "AAA", 1,
        None, 0, 0, 0,
        None, 0, None)
    assert config.host.C[0].infecting_virus.nuc_sequence == "AAA"
    assert config.host.C[1].infecting_virus.nuc_sequence == "GGG"

    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 0, 0,
          0, 0, 0, 0,
          {}, 0, "AAA", 1,
          None, 0, 0, 0,
          None, 0, 0)
    assert config.host.C[0].active
    assert config.host.C[0].infecting_virus.nuc_sequence == "AAA"
    config = Config(1, 0, {"founder0": "AAA", "founder1": "GGG"}, 
          q, 0, 0,
          0, 0, 0, 0,
          {}, 0, "AAA", 1,
          None, 0, 0, 0,
          None, 0, 0)
    assert config.host.C[0].infecting_virus.nuc_sequence == "AAA"
    assert config.host.C[1].infecting_virus.nuc_sequence == "GGG"


def test_get_nucleotide_substitution_probabilities():
    config = Config(1, 0, {"founder0": "AAA", "founder1": "GGG"}, 
          q, 3.5e-5, 0,
          0, 0, 0, 0,
          {}, 0, "AAA", 1,
          None, 0, 0, 0,
          None, 0, 0)
    assert config.substitution_probabilities[0] == ("A", "C", "G", "T")
    assert config.substitution_probabilities[1][1] == [0.2857040450786617, 0.0, 0.02855551755714605, 0.6857404373641922]
    

def test_get_substitution():
    seed(123)
    config = Config(1, 0, {"founder0": "AAA", "founder1": "GGG"}, 
          q, 3.5e-5, 0,
          0, 0, 0, 0,
          {}, 0, "AAA", 1,
          None, 0, 0, 0,
          None, 0, 0)
    assert get_substitution("A", config.substitution_probabilities) != "A"
    assert get_substitution("C", config.substitution_probabilities) != "C"
    assert get_substitution("G", config.substitution_probabilities) != "G"
    assert get_substitution("T", config.substitution_probabilities) != "T"
    with pytest.raises(Exception):
        get_substitution("X", config.substitution_probabilities)
    with pytest.raises(Exception):
        get_substitution("A")


def test_get_recomb_breakpoints():
    config = Config(1, 0, {"founder0": "AAA", "founder1": "GGG"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {}, 0, "AAA", 1,
          None, 0, 0, 0,
          None, 0, 0)
    assert (config.prob_recomb == np.full(2, 0.5)).all() # probability of an odd number of recombinations
    config = Config(1, 0, {"founder0": "AAA", "founder1": "GGG"},
          q, 3.5e-5, np.array([0, 100]), 
          0, 0, 0, 0,
          {}, 0, "AAA", 1,
          None, 0, 0, 0,
          None, 0, 0)
    assert (config.prob_recomb == np.array([0, 0.5])).all()
    # force prob_recomb to be 1
    config.prob_recomb = np.full(2, 1)
    nc, bp = get_recomb_breakpoints(1, config)
    assert nc == 1
    assert list(bp) == [[1, 2]]
    nc, bp = get_recomb_breakpoints(2, config)
    assert nc == 2
    assert list(bp) == [[1, 2], [1, 2]]
    config.recrate_is_sparse = True
    config.prob_recomb = np.full(2, 0.5)
    config.base_prob = 0.5
    nc, bp = get_recomb_breakpoints(2, config)
    assert nc == 1
    assert list(bp) == [[1, 2]]
    config.prob_recomb = np.full(2, 0.25)
    config.base_prob = 0.25
    nc, bp = get_recomb_breakpoints(2, config)
    assert nc == 2
    assert list(bp) == [[2], [2]]
    config.prob_recomb = np.full(2, 0.5)
    config.base_prob = 0.5
    nc, bp = get_recomb_breakpoints(5, config)
    assert nc == 1
    assert list(bp) == [[2]]
    config.prob_recomb = np.array([1,1])
    config.base_prob = 1
    nc, bp = get_recomb_breakpoints(1, config)
    assert nc == 1
    assert list(bp) == [[1, 2]]
    
    config.prob_recomb = np.array([1, 0])
    nc, bp = get_recomb_breakpoints(1, config)
    assert nc == 1
    assert list(bp) == [[1]]
    config.prob_recomb = np.array([0, 1])
    nc, bp = get_recomb_breakpoints(1, config)
    assert nc == 1
    assert list(bp) == [[2]]
    config.recrate_is_sparse = False
    nc, bp = get_recomb_breakpoints(1, config)
    assert nc == 1
    assert list(bp) == [[2]]

    seq_len = 5
    num_cells = 3

    # All zeros: should be no recombination
    config.prob_recomb = np.zeros(seq_len - 1)
    config.recrate_is_sparse = True
    config.base_prob = 0
    nc, bp = get_recomb_breakpoints(num_cells, config)
    assert nc == 0
    assert bp == []

    # All ones: should always have maximum recombination
    config.prob_recomb = np.ones(seq_len - 1)
    nc, bp = get_recomb_breakpoints(num_cells, config)
    assert nc >= 1
    for bplist in bp:
        for b in bplist:
            assert 1 <= b <= seq_len - 1

    # All unique rates (dense, not all equal or sparse)
    config.prob_recomb = np.array([0.1, 0.2, 0.3, 0.4])
    config.recrate_is_sparse = False
    config.base_prob = 0
    nc, bp = get_recomb_breakpoints(num_cells, config)
    assert isinstance(bp, list)

    nc, bp = get_recomb_breakpoints(1, config)
    assert nc == 1
    assert list(bp) == [[2]]
    nc, bp = get_recomb_breakpoints(3, config)
    assert nc == 4
    assert list(bp) == [[1], [1], [1], [1]]
    config.prob_recomb = 0.5
    config.base_prob = 0.5
    nc, bp = get_recomb_breakpoints(6, config)
    assert nc == 3
    assert list(bp) == [[1, 2], [2], [2]]
    config.prob_recomb = 0
    config.base_prob = 0
    nc, bp = get_recomb_breakpoints(2, config)
    assert nc == 0
    assert list(bp) == []
    config.prob_recomb = 0.25
    config.base_prob = 0
    nc, bp = get_recomb_breakpoints(3, config)
    assert nc == 1
    assert list(bp) == [[1]]
    config.base_prob = 0
    nc, bp = get_recomb_breakpoints(3, config)
    assert nc == 0

    with pytest.raises(ValueError):
      Config(1, 0, {"founder0": "AAA", "founder1": "GGG"},
          q, 3.5e-5, np.full(0, 10), 
          0, 0, 0, 0,
          {}, 0, "AAA", 1,
          None, 0, 0, 0,
          None, 0, 0)


def test_get_recombined_sequence():
    seq1 = "AAA"
    seq2 = "TTT"
    assert get_recombined_sequence(seq1, seq2, [1]) == ("ATT", [1])
    assert get_recombined_sequence(seq1, seq2, [1, 2]) == ("ATA", [1, 2])
    with pytest.raises(Exception):
        get_recombined_sequence(seq1[0], seq2[0], [1])
    with pytest.raises(Exception):
        get_recombined_sequence(seq1, seq2, [0])
    with pytest.raises(Exception):
        get_recombined_sequence(seq1, seq2, 1)
        
    seq1 = "ACGTACGT"
    seq2 = "ACGTACGT"
    # Try with various breakpoints
    for breakpoints in ([], [1], [3, 5], [1, 3, 5, 7]):
        out_seq, out_bp = get_recombined_sequence(seq1, seq2, breakpoints)
        assert out_seq == seq1
        assert out_seq == seq2


def test_calc_seq_fitness():
    assert calc_seq_fitness(0, 0.99) == 1
    assert calc_seq_fitness(1, 0.99) == (1-0.99)
    assert calc_seq_fitness(2, 0.99) == (1-0.99)**2


def test_replicative_fitness():
    assert calc_seq_fitness(muts_rel_ref("AAA", "AAA"), 0.99) == 1
    assert calc_seq_fitness(muts_rel_ref("AAA", "AAT"), 0.99) == (1-0.99)
    assert calc_seq_fitness(muts_rel_ref("AAA", "AAN"), 0.99) == 1
    assert calc_seq_fitness(muts_rel_ref("AAA", "AA-"), 0.99) == 1


def test_normalize():
    nums = [1, 2, 3, 4]
    assert normalize(nums) == [0.1, 0.2, 0.3, 0.4]
    assert normalize([1]) == [1]
    with pytest.raises(ValueError):
        normalize([0])
    with pytest.raises(Exception):
        normalize(1)
    with pytest.raises(Exception):
        normalize(["a"])


def test_get_conserved_sites_mutated():
    assert get_conserved_sites_mutated(set(), set(), [1], 3) == set()
    assert get_conserved_sites_mutated(set([0]), set(), [1], 3) == set([0])
    assert get_conserved_sites_mutated(set(), set([0]), [1], 3) == set()
    assert get_conserved_sites_mutated(set([2]), set(), [1], 3) == set()
    assert get_conserved_sites_mutated(set(), set([2]), [1], 3) == set([2])
    assert get_conserved_sites_mutated(
        set([0]), set([2]), [1], 3) == set([0, 2])
    assert get_conserved_sites_mutated(set([2]), set([0]), [1], 3) == set()
    assert get_conserved_sites_mutated(set([2]), set(), [1, 2], 3) == set([2])
    assert get_conserved_sites_mutated(set(), set([1]), [1, 2], 3) == set([1])


# Epitopes class


def test_BEpitopes():
    epitope = BEpitope(0, 3, 0.3)
    assert repr(epitope) == "(0 to 3, maxfit: 0.3)"
    assert str(epitope) == "(0 to 3, maxfit: 0.3)"


# HIV class


def test_HIV():
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {}, 0, "AAA", 0.99,
          None, 0, 0, 0,
          None, 0, 1234)
    hiv = HIV("AAT", config)
    assert repr(hiv) == "HIV with sequence AAT"
    assert str(hiv) == "HIV with sequence AAT"
    assert hiv.conserved_sites_mutated == set()
    assert hiv.b_immune_fitness == 1
    assert hiv.t_immune_fitness == 1
    assert hiv.conserved_fitness == 1
    assert hiv.replicative_fitness == 1-0.99
    assert hiv.fitness == 1-0.99
    config.ref_seq = ""
    hiv = HIV("AAT", config)
    assert hiv.replicative_fitness == 1


def test_mutate():
    seed(123)
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {1: 'A'}, 0.00, "AAA", 0,
          None, 0, 0, 0,
          None, 0, 1234)
    hiv = HIV("AAA", config)
    hiv.mutate(0, config)
    assert hiv.nuc_sequence == "TAA"
    assert hiv.conserved_sites_mutated == set()
    assert hiv.replicative_fitness == 1
    config.ref_seq = "AAA"
    config.replicative_cost = 0.1
    hiv.mutate(1, config)
    assert hiv.nuc_sequence == "TGA"
    assert hiv.conserved_sites_mutated == set([1])
    assert (
        hiv.replicative_fitness == (1-0.1)**2
    )  # don't account for conserved sites here, it's done in advance
    hiv.mutate(0, config)
    assert hiv.nuc_sequence == "AGA"
    assert hiv.conserved_sites_mutated == set([1])
    assert hiv.replicative_fitness == (1-0.1)**1

    hiv.mutate(1, config)
    assert hiv.nuc_sequence == "ATA"
    assert hiv.conserved_sites_mutated == set(
        [1])  # now allow it to mutate back
    assert hiv.replicative_fitness == (1-0.1)

    hiv.mutate(1, config)
    assert hiv.nuc_sequence == "AGA"
    assert hiv.conserved_sites_mutated == set([1])
    assert hiv.replicative_fitness == (1-0.1)**1
    
    hiv.mutate(1, config) # mutate back
    assert hiv.nuc_sequence == "AAA"
    assert hiv.conserved_sites_mutated == set()
    
    with pytest.raises(Exception):
        hiv.mutate(10, config)
    with pytest.raises(Exception):
        hiv.mutate("A", config)
    with pytest.raises(Exception):
        hiv.mutate(1.4, config)


# InfectedCD4 class


def test_InfectedCD4():
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {1: 'A'}, 0.99, "AAA", 0,
          None, 0, 0, 0,
          None, 0, 1234)
    hiv = HIV("AAA", config)
    inf_cell = InfectedCD4(hiv, True)
    assert repr(inf_cell) == "Infected CD4. Active: True. HIV with sequence AAA"
    assert str(inf_cell) == "Infected CD4. Active: True. HIV with sequence AAA"
    inf_cell.become_latent()
    assert inf_cell.active == False
    inf_cell.become_active()
    assert inf_cell.active == True


# HostEnv class


def test_HostEnv():
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {1: 'A'}, 0.99, "AAA", 0,
          None, 0, 0, 0,
          None, 0, 1234)
    host = HostEnv([HIV(seq, config)
                   for seq in ["AAA"] * 10], 10)
    assert (
        repr(host)
        == "Host has 10 active and 0 latent infected cells\n0 epitopes recognized"
    )
    assert (
        str(host)
        == "Host has 10 active and 0 latent infected cells\n0 epitopes recognized"
    )
    assert host.C[0].active == True
    assert host.epitopes_recognition_generation == defaultdict(lambda: 0)
    assert host.epitope_variants_translated == defaultdict(lambda: "")
    host.C[0].infecting_virus.nuc_sequence = "ATA"
    assert host.C[1].infecting_virus.nuc_sequence == "AAA"
    host = HostEnv([HIV(seq, config)
                   for seq in ["AAA", "TTT"] * 5], 10)
    assert "AAA" in [host.C[i].infecting_virus.nuc_sequence for i in range(10)]
    assert "TTT" in [host.C[i].infecting_virus.nuc_sequence for i in range(10)]


def test_translate():
    config = Config(1, 0, {"founder0": "ATGATTGTGTAG"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {1: 'A'}, 0.99, "ATGATTGTGTAG", 0,
          None, 0, 0, 0,
          None, 0, 1234)
    hiv = [HIV("ATGATTGTGTAG", config)]
    host = HostEnv(hiv, 1)
    assert translate(hiv[0].nuc_sequence) == "MIV_"
    assert translate(hiv[0].nuc_sequence[3:9])
    with pytest.raises(Exception):
        translate(hiv[0].nuc_sequence[0:7])


def test_update_b_epitopes_recognized():
    config = Config(1, 0, {"founder0": "AAAAAAAAA"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {1: 'A'}, 0.99, "AAAAAAAAA", 0,
          pd.DataFrame({'start': [0], 'end': [3], 'maxepi': [0.3]}), 0, 2, 30,
          None, 0, 1234)
    # rng = default_rng(1234)
    # reference_sequence = "AAAAAAAAA"
    host = HostEnv([HIV(seq, config)
                   for seq in ["AAAAAAAAA"] * 11], 11)
    # immune_response_num = 2
    host.update_b_epitopes_recognized(
        10, config)
    assert host.epitopes_recognition_generation == {"K": 10}
    assert host.epitope_variants_translated == {"AAA": "K"}
    host.C[0].infecting_virus.nuc_sequence = "ATAAAAAAA"
    host.update_b_epitopes_recognized(
        15, config)
    assert host.epitopes_recognition_generation == {"K": 10}
    assert host.epitope_variants_translated == {"AAA": "K", "ATA": "I"}
    host.C[1].infecting_virus.nuc_sequence = "ATAAAAAAA"
    host.update_b_epitopes_recognized(
        20, config)
    assert host.epitopes_recognition_generation == {
        "K": 10,
        "I": 16,
    }  # because of cross-reactivity
    assert host.epitope_variants_translated == {"AAA": "K", "ATA": "I"}


def test_update_b_immune_fitness():
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {1: 'A'}, 0.99, "AAA", 0,
          pd.DataFrame({'start': [0, 4], 'end': [3, 7], 'maxepi': [0.2, 0.3]}),
          0, 0.1, 90,
          None, 0, 1234)
    host = HostEnv([HIV(seq, config)
                   for seq in ["AAAAAAAAA"] * 11], 11)
    host.update_b_epitopes_recognized(10, config)
    host.update_b_immune_fitness(40, config)
    assert (
        round(host.C[0].infecting_virus.b_immune_fitness, 1)
        == 1 - 0.3 * (40 - 10) / 90
    )


def test_get_fitness_of_infecting_virus():
    #seed(123)
    #rng = default_rng(1234)
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {1: 'A'}, 0.99, "AAA", 0,
          pd.DataFrame({'start': [0], 'end': [3], 'maxepi': [0.3]}), 0, 10, 30,
          None, 0, 1234)
    hiv = HIV("AAA", config)
    host = HostEnv([HIV(seq, config)
                   for seq in ["AAA"] * 11], 11)
    assert host.get_fitness_of_infecting_virus(0) == 1
    assert host.C[0].infecting_virus.fitness == 1
    host.C[0].infecting_virus.mutate(1, config)
    assert host.get_fitness_of_infecting_virus(
        0) == (1 - 0) * (1 - 0) * (1-0.99)
    assert host.C[0].infecting_virus.b_immune_fitness == 1
    assert host.C[0].infecting_virus.t_immune_fitness == 1
    assert host.C[0].infecting_virus.replicative_fitness == 1
    assert host.C[0].infecting_virus.conserved_fitness == (1-0.99)
    assert host.C[0].infecting_virus.fitness == (1 - 0) * (1 - 0) * (1-0.99)
    assert host.get_fitness_of_infecting_virus(1) == 1
    assert host.C[1].infecting_virus.b_immune_fitness == 1
    assert host.C[0].infecting_virus.t_immune_fitness == 1
    assert host.C[1].infecting_virus.replicative_fitness == 1
    assert host.C[1].infecting_virus.conserved_fitness == 1
    host.update_b_epitopes_recognized(1, config)
    host.update_b_immune_fitness(40, config)
    assert host.get_fitness_of_infecting_virus(
        0) == (1 - 0) * (1 - 0) * (1-0.99)
    assert (
        host.get_fitness_of_infecting_virus(1)
        == 1 - 0.3 * (40 - 10) / 30
    )
    with pytest.raises(Exception):
        host.get_fitness_of_infecting_virus(20)


def test_singly_infect_cd4():
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {1: 'A'}, 0.99, "AAA", 0,
          pd.DataFrame({'start': [0], 'end': [3], 'maxepi': [0.3]}), 0, 10, 30,
          None, 0, 123)
    # seed(123)
    
    # reference_sequence = "AAA"
    host = HostEnv([HIV(seq, config) for seq in ["AAA"] * 4], 4)
    host.C[0].infecting_virus.nuc_sequence = "GGG"
    host.C[1].infecting_virus.nuc_sequence = "TTT"
    newly_infected = host.singly_infect_cd4([0, 1, 2, 0])
    assert [newly_infected[i].infecting_virus.nuc_sequence for i in range(4)] == [
        "GGG",
        "GGG",
        "TTT",
        "AAA",
    ]
    newly_infected[0].infecting_virus.mutate(0, config)
    newly_infected[1].infecting_virus.mutate(1, config)
    assert [newly_infected[i].infecting_virus.nuc_sequence for i in range(2)] == [
        "AGG",
        "GAG",
    ]


def test_dually_infect_cd4():
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {1: 'A'}, 0.99, "AAA", 0.1,
          pd.DataFrame({'start': [0], 'end': [3], 'maxepi': [0.3]}), 0, 10, 30,
          None, 0, 123)
    host = HostEnv([HIV(seq, config) for seq in ["AAA"] * 3], 3)
    host.C[0].infecting_virus.nuc_sequence = "GGG"
    host.C[1].infecting_virus.nuc_sequence = "TTT"
    host.C[0].infecting_virus.conserved_sites_mutated = set([2])
    newly_infected = host.dually_infect_cd4([0, 1, 0, 2], [[1], [1, 2]], config)
    host.C[0].infecting_virus.mutate(0, config)
    assert [newly_infected[i].infecting_virus.nuc_sequence for i in range(2)] == [
        "GTT",
        "GAG",
    ]
    assert [
        newly_infected[i].infecting_virus.conserved_fitness for i in range(2)
    ] == [1, (1-0.99)]
    assert [
        newly_infected[i].infecting_virus.replicative_fitness for i in range(2)
    ] == [(1-0.1)**2, (1-0.1)**2]


def test_latent_active_CD4():
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {1: 'A'}, 0.99, "AAA", 0.1,
          pd.DataFrame({'start': [0], 'end': [3], 'maxepi': [0.3]}), 0, 10, 30,
          None, 0, 123)
    host = HostEnv([HIV(seq, config)
                   for seq in ["AAA"] * 10], 10)
    assert len(host.C) == 10
    assert len(host.L) == 0
    host.make_latent(0)
    assert len(host.L) == 1
    host.proliferate_latent_CD4(0)
    assert len(host.L) == 2
    host.die_latent_CD4(0)
    assert len(host.L) == 1
    assert len(host.C) == 9
    assert host.L[0].active == False
    host.make_active(0)
    assert len(host.L) == 0
    assert len(host.C) == 10
    assert host.C[9].active == True
    with pytest.raises(Exception):
        host.proliferate_latent_CD4(100)
    with pytest.raises(Exception):
        host.die_latent_CD4(100)
    with pytest.raises(Exception):
        host.make_latent(100)
    with pytest.raises(Exception):
        host.make_active(100)


def test_mutate_virus_in_productive_CD4():
    seed(123)
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {1: 'A'}, 0.99, "AAA", 0.1,
          pd.DataFrame({'start': [0], 'end': [3], 'maxepi': [0.3]}), 0, 10, 30,
          None, 0, 123)
    host = HostEnv([HIV(seq, config)
                   for seq in ["AAA"] * 10], 10)
    host.mutate_virus_in_productive_CD4([1], config)
    assert host.C[0].infecting_virus.nuc_sequence == "ACA"
    assert host.C[0].infecting_virus.conserved_sites_mutated == set([1])
    assert host.C[1].infecting_virus.nuc_sequence == "AAA"
    assert host.C[1].infecting_virus.conserved_sites_mutated == set()
    host.mutate_virus_in_productive_CD4([4, 9], config)
    assert host.C[0].infecting_virus.nuc_sequence == "ACA"
    assert host.C[1].infecting_virus.nuc_sequence == "ACA"
    assert host.C[2].infecting_virus.nuc_sequence == "AAA"
    assert host.C[3].infecting_virus.nuc_sequence == "GAA"
    with pytest.raises(Exception):
        host.mutate_virus_in_productive_CD4(1, config)
    with pytest.raises(Exception):
        host.mutate_virus_in_productive_CD4([300000], config)


def test_get_next_gen_latent():
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 100,
          100, 0.00001, 0, 0,
          {1: 'A'}, 0.99, "AAA", 0.1,
          pd.DataFrame({'start': [0], 'end': [3], 'maxepi': [0.3]}), 0, 10, 30,
          None, 0, 1234)
    host = HostEnv([HIV(seq, config) for seq in ["GGG"] * 3], 3)
    assert host.get_next_gen_latent(config) == (2, 0, 0, 0)
    assert host.L[0].active == False
    assert host.L[0].infecting_virus.nuc_sequence == "GGG"
    assert len(host.L) == 2
    assert len(host.C) == 1
    config.prob_act_to_lat = 0
    config.prob_lat_to_act = 1
    config.prob_lat_die = 0
    config.prob_lat_prolif = 0
    assert host.get_next_gen_latent(config) == (0, 2, 0, 0)
    assert len(host.L) == 0
    assert len(host.C) == 3
    config.prob_act_to_lat = 1
    config.prob_lat_to_act = 0.00001
    host.get_next_gen_latent(config) == (2, 0, 0, 0)
    config.prob_act_to_lat = 0
    config.prob_lat_to_act = 0
    config.prob_lat_prolif = 1
    assert host.get_next_gen_latent(config) == (0, 0, 0, 2)
    assert len(host.L) == 4
    assert len(host.C) == 1
    config.prob_act_to_lat = 0
    config.prob_lat_to_act = 0.1
    config.prob_lat_die = 0.1
    config.prob_lat_prolif = 0.1
    assert host.get_next_gen_latent(config) == (0, 1, 0, 0)
    config.prob_act_to_lat = 0
    config.prob_lat_to_act = 0
    config.prob_lat_die = 1
    config.prob_lat_prolif = 0
    assert host.get_next_gen_latent(config) == (0, 0, 3, 0)
    assert len(host.L) == 0
    assert len(host.C) == 2
    config.prob_act_to_lat = 0
    config.prob_lat_to_act = 0
    config.prob_lat_die = 0
    config.prob_lat_prolif = 2
    with pytest.raises(Exception):
        host.get_next_gen_latent(config)


def test_get_next_gen_active():
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 0,
          0, 0, 0, 0,
          {1: 'A'}, 0.99, "TTT", 0.1,
          pd.DataFrame({'start': [0], 'end': [3], 'maxepi': [0.3]}), 30, 0.2, 30,
          None, 0, 1234)
    host = HostEnv([HIV(seq, config) for seq in ["AAA"] * 3], 3)
    host.C[0].infecting_virus.mutate(0, config)
    host.C[1].infecting_virus.mutate(0, config)
    assert host.get_next_gen_active(2, 40, config) == (0, 0)
    assert [host.C[i].infecting_virus.nuc_sequence for i in range(2)] == [
        "TAA", "AAA"]
    assert [host.C[i].infecting_virus.fitness for i in range(2)] == [(1-0.1)**2, (1-0.1)**3]
    config.prob_mut = 0.5
    assert host.get_next_gen_active(10, 40, config) == (3, 0)
    assert [host.C[i].infecting_virus.nuc_sequence for i in range(1)] == [
        "CAA"]
    config.prob_mut = 0.1
    config.prob_recomb = 0.1

    config.base_prob = 0.1
    assert host.get_next_gen_active(10, 40, config) == (3, 0)
    assert "CAA" in [host.C[i].infecting_virus.nuc_sequence for i in range(10)]


def test_summarize_fitness():
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {1: 'A'}, 0.99, "AAA", 0.1,
          pd.DataFrame({'start': [0], 'end': [3], 'maxepi': [0.3]}), 30, 0.2, 30,
          None, 0, 1234)
    host = HostEnv([HIV(seq, config)
                   for seq in ["GGG", "AAA"] * 2], 4)
    assert [host.C[i].infecting_virus.nuc_sequence for i in range(len(host.C))] == [
        "GGG",
        "AAA",
        "GGG",
        "AAA",
    ]
    host.C[1].infecting_virus.mutate(1, config)
    assert [
        host.C[i].infecting_virus.replicative_fitness for i in range(len(host.C))
    ] == [(1-0.1)**2, 1, (1-0.1)**2, 1]
    assert [
        host.C[i].infecting_virus.conserved_fitness for i in range(len(host.C))
    ] == [1, (1-0.99), 1, 1]
    assert [
        host.C[i].infecting_virus.b_immune_fitness for i in range(len(host.C))
    ] == [1, 1, 1, 1]
    assert [
        host.C[i].infecting_virus.t_immune_fitness for i in range(len(host.C))
    ] == [1, 1, 1, 1]
    host.get_fitness()
    assert [host.C[i].infecting_virus.fitness for i in range(len(host.C))] == [
        (1-0.1)**2,
        (1-0.99),
        (1-0.1)**2,
        1.0,
    ]
    assert host.summarize_fitness() == (0.6575, 0.7525, 1.0, 1.0, 0.905)


def test_record_counts():
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
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {1: 'A'}, 0.99, "TTT", 0.1,
          pd.DataFrame({'start': [0], 'end': [3], 'maxepi': [0.3]}), 30, 0.2, 30,
          None, 0, 1234)
    assert list(
        config.host.record_counts(
            counts, 1, (1, 2, 3, 4), (5, 6), (7, 8, 9, 10, 11)).values()) == [
        [1], [1], [0], [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11]]
    assert list(
        config.host.record_counts(
            counts, 2, (11, 12, 13, 14), (15, 16), (17, 18, 19, 20, 21)).values()) == [
        [
            1, 2], [
            1, 1], [
            0, 0], [
            1, 11], [
            2, 12], [
            3, 13], [
            4, 14], [
            5, 15], [
            6, 16], [
            7, 17], [
            8, 18], [
            9, 19], [
            10, 20], [
            11, 21
            ]]


def test_sample_viral_sequences():
    config = Config(1, 0, {"founder0": "AAA"}, 
          q, 3.5e-5, 100,
          0, 0, 0, 0,
          {1: 'A'}, 0.99, "AAA", 0.1,
          pd.DataFrame({'start': [0], 'end': [3], 'maxepi': [0.3]}), 30, 0.2, 30,
          None, 0, 1234)
    fitness = {
        "generation": [],
        "seq_id": [],
        "b_immune": [],
        "t_immune": [],
        "conserved": [],
        "replicative": [],
        "overall": []
    }
    assert config.host.sample_viral_sequences(
        config.founder_seqs,
        fitness,
        1,
        1) == ({'founder0': 'AAA', 'gen1_active_0': 'AAA'}, {'generation': ['1'], 'seq_id': ['gen1_active_0'], 'b_immune': [1.0], 't_immune': [1.0], 'conserved': [1.0], 'replicative': [1.0], 'overall': [1.0]})


def test_loop_through_generations():
    config = Config(1, 2, {"founder0": "AAA"}, 
          q, 3.5e-5, 0,
          0.1, 0, 0, 0,
          {1: 'A'}, 0.99, "", 0.1,
          pd.DataFrame({'start': [0], 'end': [3], 'maxepi': [0.3]}), 30, 0.2, 30,
          None, 0, 1234)
    out = config.host.loop_through_generations([1, 2, 3], [1, 2, 3], [0, 1, 1], config)
    assert out[0] == {
        'generation': [
            0, 1, 2], 'active_cell_count': [
            1, 2, 3], 'latent_cell_count': [
                0, 0, 1], 'active_turned_latent': [
                    0, 0, 1], 'latent_turned_active': [
                        0, 0, 0], 'latent_died': [
                            0, 0, 0], 'latent_proliferated': [
                                0, 0, 0], 'number_mutations': [
                                    0, 0, 0], 'number_recombinations': [
                                        0, 0, 0], 'mean_fitness_active': [
                                            1.0, 1.0, 1.0], 'mean_conserved_active': [
                                                1.0, 1.0, 1.0], 'mean_b_immune_active': [
                                                    1.0, 1.0, 1.0], 'mean_t_immune_active': [
                                                    1.0, 1.0, 1.0], 'mean_replicative_active': [
                                                        1.0, 1.0, 1.0]}
    assert out[1] == {'generation': ['founder', '0', '1', '1', '2', '2', '2'], 'seq_id': ['founder0', 'gen0_active_0', 'gen1_active_0', 'gen1_active_1', 'gen2_active_0', 'gen2_active_1', 'gen2_active_2'], 
    'b_immune': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 't_immune': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 
    'conserved': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 'replicative': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 'overall': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]}
    assert out[2] == {'founder0': 'AAA', 'gen0_active_0': 'AAA', 'gen1_active_0': 'AAA',
                      'gen1_active_1': 'AAA', 'gen2_active_0': 'AAA', 'gen2_active_1': 'AAA', 'gen2_active_2': 'AAA'}
    
    config = Config(1, 2, {"founder0": "AAA"}, 
          q, 3.5e-5, 0,
          0.001, 0, 0, 0,
          {1: 'A'}, 0.99, "ACA", 0.1,
          pd.DataFrame({'start': [0], 'end': [3], 'maxepi': [0.3]}), 30, 0.2, 30,
          None, 0, 1234)
    out = config.host.loop_through_generations([1, 2, 3], [1, 2, 3], [0, 1, 1], config)
    assert out[0]['mean_replicative_active'] == [0.9,0.9,0.9]
