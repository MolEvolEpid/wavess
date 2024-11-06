import pytest
from random import seed
from numpy.random import default_rng

from run_wavess import *
from agents import *

# Functions outside of any class


def test_set_python_seed():
    assert isinstance(set_python_seed(None), np.random._generator.Generator)
    g = set_python_seed(1)
    assert sample(range(1000), 1) == [137]
    assert g.binomial(1000, 0.1) == 100


def test_prep_ref_conserved():
    assert prep_ref_conserved({'founder0': 'AAA'}, 'AAA', {}) == ('AAA', {})
    assert prep_ref_conserved({'founder0': 'AAA'}, 'AAA', {
                              0: 'a'}) == ('-AA', {0: 'a'})
    assert prep_ref_conserved({'founder0': 'TAA'}, 'AAA', {
                              0: 'a'}) == ('AAA', {})
    assert prep_ref_conserved({'founder0': 'TAA', 'founder1': 'TAA'}, 'AAA', {
                              0: 'a', 2: 'c'}) == ('AA-', {2: 'c'})
    assert prep_ref_conserved({'founder0': 'TAA', 'founder1': 'AAT'}, 'AAA', {
                              0: 'a', 2: 'c'}) == ('AAA', {})


def test_create_host_env():
    host = create_host_env({"founder0": "AAA"}, "AAA", 1, 1)
    assert host.C[0].active
    assert host.C[0].infecting_virus.nuc_sequence == "AAA"
    host = create_host_env({"founder0": "AAA", "founder1": "GGG"}, "AAA", 2, 2)
    assert host.C[0].infecting_virus.nuc_sequence == "AAA"
    assert host.C[1].infecting_virus.nuc_sequence == "GGG"
    with pytest.raises(Exception):
        create_host_env({"founder0": "AAA"}, "AAA", 1, 10)
    with pytest.raises(Exception):
        create_host_env({"founder0": "AAA", "founder2": "GGG"}, "AAA", 1, 1)


def test_create_epitope():
    epitope = create_epitope(10, 20, 0.4)
    assert vars(epitope) == {'start': 10, 'end': 20, 'max_fitness': 0.4}


def test_get_nucleotide_substitution_probabilities():
    nt_sub_probs = get_nucleotide_substitution_probabilities(
        "../extdata/hiv_q_mat.csv", 3.5e-5
    )
    assert nt_sub_probs[0] == ("A", "C", "G", "T")
    assert nt_sub_probs[1][1] == [0.28570164890645106, 0.0, 0.028552030139711436, 0.6857463209538375]
    with pytest.raises(Exception):
        get_nucleotide_substitution_probabilities("", 3.5e-5)
    with pytest.raises(Exception):
        get_nucleotide_substitution_probabilities(0, 3.5e-5)


def test_get_substitution():
    seed(123)
    new_nt_order, probs = get_nucleotide_substitution_probabilities(
        "../extdata/hiv_q_mat.csv", 3.5e-5
    )
    assert get_substitution("A", new_nt_order, probs) == "C"
    assert get_substitution("C", new_nt_order, probs) == "A"
    assert get_substitution("G", new_nt_order, probs) == "A"
    assert get_substitution("T", new_nt_order, probs) == "A"
    with pytest.raises(Exception):
        get_substitution("X", new_nt_order, probs)
    with pytest.raises(Exception):
        get_substitution("A")


def test_get_recomb_breakpoints():
    rng = default_rng(1234)
    nc, bp = get_recomb_breakpoints(3, 1, 1, rng)
    assert nc == 1
    assert list(bp) == [[1, 2]]
    nc, bp = get_recomb_breakpoints(3, 2, 1, rng)
    assert nc == 2
    assert list(bp) == [[1, 2], [2, 1]]
    nc, bp = get_recomb_breakpoints(3, 2, 0.5, rng)
    assert nc == 2
    assert list(bp) == [[2], [2, 1]]
    nc, bp = get_recomb_breakpoints(3, 2, 0.25, rng)
    assert nc == 0
    assert list(bp) == []
    nc, bp = get_recomb_breakpoints(3, 3, 0.25, rng)
    assert nc == 1
    assert list(bp) == [[2]]


def test_get_recombined_sequence():
    # seed(123)
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


def test_Epitopes():
    epitope = Epitope(0, 3, 0.3)
    assert repr(epitope) == "(0 to 3, maxfit: 0.3)"
    assert str(epitope) == "(0 to 3, maxfit: 0.3)"


# HIV class


def test_HIV():
    reference_sequence = "AAA"
    hiv = HIV("AAT", reference_sequence, 0.99)
    assert repr(hiv) == "HIV with sequence AAT"
    assert str(hiv) == "HIV with sequence AAT"
    assert hiv.conserved_sites_mutated == set()
    assert hiv.immune_fitness == 1
    assert hiv.conserved_fitness == 1
    assert hiv.replicative_fitness == 1-0.99
    assert hiv.fitness == 1-0.99
    hiv = HIV("AAT", "", 1)
    assert hiv.replicative_fitness == 1


def test_mutate():
    seed(123)
    new_nt_order, probs = get_nucleotide_substitution_probabilities(
        "../extdata/hiv_q_mat.csv", 3.5e-5
    )
    reference_sequence = ""
    hiv = HIV("AAA", reference_sequence, 1)
    hiv.mutate(0, new_nt_order, probs, {1: 'A'}, 0.99, reference_sequence, 1)
    assert hiv.nuc_sequence == "CAA"
    assert hiv.conserved_sites_mutated == set()
    assert hiv.replicative_fitness == 1
    reference_sequence = "AAA"
    hiv.mutate(1, new_nt_order, probs, {1: 'A'}, 0.99, reference_sequence, 0.1)
    assert hiv.nuc_sequence == "CCA"
    assert hiv.conserved_sites_mutated == set([1])
    assert (
        hiv.replicative_fitness == (1-0.1)**2
    )  # don't account for conserved sites here, it's done in advance
    hiv.mutate(0, new_nt_order, probs, {1: 'A'}, 0.99, reference_sequence, 0.1)
    assert hiv.nuc_sequence == "TCA"
    assert hiv.conserved_sites_mutated == set([1])
    assert hiv.replicative_fitness == (1-0.1)**2

    hiv.mutate(1, new_nt_order, probs, {1: 'A'}, 0.99, reference_sequence, 0.1)
    assert hiv.nuc_sequence == "TAA"
    assert hiv.conserved_sites_mutated == set(
        [])  # now allow it to mutate back
    assert hiv.replicative_fitness == (1-0.1)

    hiv.mutate(1, new_nt_order, probs, {1: 'A'}, 0.99, reference_sequence, 0.1)
    assert hiv.nuc_sequence == "TGA"
    assert hiv.conserved_sites_mutated == set([1])
    assert hiv.replicative_fitness == (1-0.1)**2
    with pytest.raises(Exception):
        hiv.mutate(10, 0.99, reference_sequence, 1)
    with pytest.raises(Exception):
        hiv.mutate("A", 0.99, reference_sequence, 1)
    with pytest.raises(Exception):
        hiv.mutate(1.4, 0.99, reference_sequence, 1)


# InfectedCD4 class


def test_InfectedCD4():
    reference_sequence = "AAA"
    hiv = HIV("AAA", reference_sequence, 1)
    inf_cell = InfectedCD4(hiv, True)
    assert repr(inf_cell) == "Infected CD4. Active: True. HIV with sequence AAA"
    assert str(inf_cell) == "Infected CD4. Active: True. HIV with sequence AAA"
    inf_cell.become_latent()
    assert inf_cell.active == False
    inf_cell.become_active()
    assert inf_cell.active == True


# HostEnv class


def test_HostEnv():
    reference_sequence = "AAA"
    host = HostEnv([HIV(seq, reference_sequence, 1)
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
    host = HostEnv([HIV(seq, reference_sequence, 1)
                   for seq in ["AAA", "TTT"] * 5], 10)
    assert "AAA" in [host.C[i].infecting_virus.nuc_sequence for i in range(10)]
    assert "TTT" in [host.C[i].infecting_virus.nuc_sequence for i in range(10)]


def test_translate():
    reference_sequence = "ATGATTGTGTAG"
    hiv = [HIV("ATGATTGTGTAG", reference_sequence, 1)]
    host = HostEnv(hiv, 1)
    assert host.translate(hiv[0].nuc_sequence) == "MIV_"
    assert host.translate(hiv[0].nuc_sequence[3:9])
    with pytest.raises(Exception):
        host.translate(hiv[0].nuc_sequence[0:7])


def test_update_epitopes_recognized():
    rng = default_rng(1234)
    reference_sequence = "AAAAAAAAA"
    host = HostEnv([HIV(seq, reference_sequence, 1)
                   for seq in ["AAAAAAAAA"] * 11], 11)
    epi = [Epitope(0, 3, 0.3)]
    immune_response_proportion = 0.1
    host.update_epitopes_recognized(
        10, epi, immune_response_proportion, 30, rng)
    assert host.epitopes_recognition_generation == {"K": 10}
    assert host.epitope_variants_translated == {"AAA": "K"}
    host.C[0].infecting_virus.nuc_sequence = "ATAAAAAAA"
    host.update_epitopes_recognized(
        15, epi, immune_response_proportion, 30, rng)
    assert host.epitopes_recognition_generation == {"K": 10}
    assert host.epitope_variants_translated == {"AAA": "K", "ATA": "I"}
    host.C[1].infecting_virus.nuc_sequence = "ATAAAAAAA"
    host.update_epitopes_recognized(
        20, epi, immune_response_proportion, 30, rng)
    assert host.epitopes_recognition_generation == {
        "K": 10,
        "I": 16,
    }  # because of cross-reactivity
    assert host.epitope_variants_translated == {"AAA": "K", "ATA": "I"}


def test_update_immune_fitness():
    rng = default_rng(1234)
    reference_sequence = "AAAAAAAAA"
    epi = [Epitope(0, 3, 0.2), Epitope(4, 7, 0.3)]
    time_to_full_potency = 90
    host = HostEnv([HIV(seq, reference_sequence, 1)
                   for seq in ["AAAAAAAAA"] * 11], 11)
    host.update_epitopes_recognized(10, epi, 0.1, time_to_full_potency, rng)
    host.update_immune_fitness(epi, 40, time_to_full_potency)
    assert (
        round(host.C[0].infecting_virus.immune_fitness, 1)
        == 1 - 0.3 * (40 - 10) / time_to_full_potency
    )


def test_get_fitness_of_infecting_virus():
    seed(123)
    rng = default_rng(1234)
    new_nt_order, probs = get_nucleotide_substitution_probabilities(
        "../extdata/hiv_q_mat.csv", 3.5e-5
    )
    reference_sequence = ""
    hiv = HIV("AAA", reference_sequence, 1)
    host = HostEnv([HIV(seq, reference_sequence, 1)
                   for seq in ["AAA"] * 11], 11)
    immune_response_proportion = 0.1
    time_to_full_potency = 30
    cost_per_mutation_in_conserved_site = 0.99
    epi = [Epitope(0, 3, 0.3)]
    assert host.get_fitness_of_infecting_virus(0) == 1
    assert host.C[0].infecting_virus.fitness == 1
    host.C[0].infecting_virus.mutate(
        1, new_nt_order, probs, {1: 'A'}, 0.99, reference_sequence, 1
    )
    assert host.get_fitness_of_infecting_virus(
        0) == (1 - 0) * (1 - 0) * (1-0.99)
    assert host.C[0].infecting_virus.immune_fitness == 1
    assert host.C[0].infecting_virus.replicative_fitness == 1
    assert host.C[0].infecting_virus.conserved_fitness == (1-0.99)
    assert host.C[0].infecting_virus.fitness == (1 - 0) * (1 - 0) * (1-0.99)
    assert host.get_fitness_of_infecting_virus(1) == 1
    assert host.C[1].infecting_virus.immune_fitness == 1
    assert host.C[1].infecting_virus.replicative_fitness == 1
    assert host.C[1].infecting_virus.conserved_fitness == 1
    host.update_epitopes_recognized(
        1, epi, immune_response_proportion, time_to_full_potency, rng
    )
    host.update_immune_fitness(epi, 40, time_to_full_potency)
    assert host.get_fitness_of_infecting_virus(
        0) == (1 - 0) * (1 - 0) * (1-0.99)
    assert (
        host.get_fitness_of_infecting_virus(1)
        == 1 - 0.3 * (40 - 10) / time_to_full_potency
    )
    with pytest.raises(Exception):
        host.get_fitness_of_infecting_virus(20)


def test_singly_infect_cd4():
    seed(123)
    new_nt_order, probs = get_nucleotide_substitution_probabilities(
        "../extdata/hiv_q_mat.csv", 3.5e-5
    )
    reference_sequence = "AAA"
    host = HostEnv([HIV(seq, reference_sequence, 1) for seq in ["AAA"] * 4], 4)
    host.C[0].infecting_virus.nuc_sequence = "GGG"
    host.C[1].infecting_virus.nuc_sequence = "TTT"
    newly_infected = host.singly_infect_cd4([0, 1, 2, 0])
    assert [newly_infected[i].infecting_virus.nuc_sequence for i in range(4)] == [
        "GGG",
        "GGG",
        "TTT",
        "AAA",
    ]
    newly_infected[0].infecting_virus.mutate(
        0, new_nt_order, probs, {1: 'A'}, 0.99, reference_sequence, 1
    )
    newly_infected[1].infecting_virus.mutate(
        1, new_nt_order, probs, {1: 'A'}, 0.99, reference_sequence, 1
    )
    assert [newly_infected[i].infecting_virus.nuc_sequence for i in range(2)] == [
        "AGG",
        "GAG",
    ]


def test_dually_infect_cd4():
    new_nt_order, probs = get_nucleotide_substitution_probabilities(
        "../extdata/hiv_q_mat.csv", 3.5e-5
    )
    reference_sequence = "AAA"
    host = HostEnv([HIV(seq, reference_sequence, 1) for seq in ["AAA"] * 3], 3)
    host.C[0].infecting_virus.nuc_sequence = "GGG"
    host.C[1].infecting_virus.nuc_sequence = "TTT"
    host.C[0].infecting_virus.conserved_sites_mutated = set([2])
    newly_infected = host.dually_infect_cd4(
        [0, 1, 0, 2], [[1], [1, 2]], 3, 0.99, reference_sequence, set([2]), 0.1
    )
    host.C[0].infecting_virus.mutate(
        0, new_nt_order, probs, {1: 'A'}, 0.99, reference_sequence, 0.1
    )
    assert [newly_infected[i].infecting_virus.nuc_sequence for i in range(2)] == [
        "GTT",
        "GAG",
    ]
    assert [
        newly_infected[i].infecting_virus.conserved_fitness for i in range(2)
    ] == [1, (1-0.99)]
    assert [
        newly_infected[i].infecting_virus.replicative_fitness for i in range(2)
    ] == [(1-0.1)**3, (1-0.1)**2]


def test_latent_active_CD4():
    reference_sequence = "AAA"
    host = HostEnv([HIV(seq, reference_sequence, 1)
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
    new_nt_order, probs = get_nucleotide_substitution_probabilities(
        "../extdata/hiv_q_mat.csv", 3.5e-5
    )
    reference_sequence = "AAA"
    host = HostEnv([HIV(seq, reference_sequence, 1)
                   for seq in ["AAA"] * 10], 10)
    host.mutate_virus_in_productive_CD4(
        [1], 3, new_nt_order, probs, {1: 'A'}, 0.99, reference_sequence, 1
    )
    assert host.C[0].infecting_virus.nuc_sequence == "ACA"
    assert host.C[0].infecting_virus.conserved_sites_mutated == set([1])
    assert host.C[1].infecting_virus.nuc_sequence == "AAA"
    assert host.C[1].infecting_virus.conserved_sites_mutated == set()
    host.mutate_virus_in_productive_CD4(
        [4, 9], 3, new_nt_order, probs, {1: 'A'}, 0.99, reference_sequence, 1
    )
    assert host.C[0].infecting_virus.nuc_sequence == "ACA"
    assert host.C[1].infecting_virus.nuc_sequence == "ACA"
    assert host.C[2].infecting_virus.nuc_sequence == "AAA"
    assert host.C[3].infecting_virus.nuc_sequence == "GAA"
    with pytest.raises(Exception):
        host.mutate_virus_in_productive_CD4(1, 3, new_nt_order, probs, [1], 1)
    with pytest.raises(Exception):
        host.mutate_virus_in_productive_CD4(
            [300000], 3, new_nt_order, probs, {1: 'A'}, 1)


def test_get_next_gen_latent():
    seed(1234)
    rng = default_rng(1234)
    reference_sequence = "AAA"
    epi = [Epitope(0, 3, 0.3)]
    host = HostEnv([HIV(seq, reference_sequence, 1) for seq in ["GGG"] * 3], 3)
    assert host.get_next_gen_latent(1, 0.00001, 0, 0, rng) == (2, 0, 0, 0)
    assert host.L[0].active == False
    assert host.L[0].infecting_virus.nuc_sequence == "GGG"
    assert len(host.L) == 2
    assert len(host.C) == 1
    assert host.get_next_gen_latent(0, 1, 0, 0, rng) == (0, 2, 0, 0)
    assert len(host.L) == 0
    assert len(host.C) == 3
    host.get_next_gen_latent(1, 0.00001, 0, 0, rng) == (2, 0, 0, 0)
    assert host.get_next_gen_latent(0, 0, 0, 1, rng) == (0, 0, 0, 2)
    assert len(host.L) == 4
    assert len(host.C) == 1
    assert host.get_next_gen_latent(0, 0.1, 0.1, 0.1, rng) == (0, 1, 0, 0)
    assert host.get_next_gen_latent(0, 0, 1, 0, rng) == (0, 0, 3, 0)
    assert len(host.L) == 0
    assert len(host.C) == 2
    with pytest.raises(Exception):
        host.get_next_gen_latent(0, 0, 0, 2, rng)


def test_get_next_gen_active():
    seed(1234)
    rng = default_rng(1234)
    new_nt_order, probs = get_nucleotide_substitution_probabilities(
        "../extdata/hiv_q_mat.csv", 3.5e-5
    )
    reference_sequence = "AAA"
    epi = [Epitope(0, 3, 0.3)]
    host = HostEnv([HIV(seq, reference_sequence, 1) for seq in ["GGG"] * 3], 3)
    host.C[0].infecting_virus.mutate(
        0, new_nt_order, probs, {1: 'A'}, 0.99, reference_sequence, 0.1
    )
    host.C[1].infecting_virus.mutate(
        0, new_nt_order, probs, {1: 'A'}, 0.99, reference_sequence, 0.1
    )
    assert host.get_next_gen_active(
        0,
        0,
        2,
        40,
        30,
        new_nt_order,
        probs,
        {1: 'A'},
        30,
        0.99,
        reference_sequence,
        0.2,
        epi,
        0.1,
        rng,
    ) == (0, 0)
    assert [host.C[i].infecting_virus.nuc_sequence for i in range(2)] == [
        "AGG", "AGG"]
    assert [
        host.C[i].infecting_virus.fitness for i in range(2)] == [(1-0.1)**2, (1-0.1)**2]
    assert host.get_next_gen_active(
        0.5,
        0,
        10,
        40,
        30,
        new_nt_order,
        probs,
        {1: 'A'},
        30,
        0.99,
        reference_sequence,
        0.2,
        epi,
        0.1,
        rng,
    ) == (5, 0)
    assert [host.C[i].infecting_virus.nuc_sequence for i in range(1)] == [
        "GAA"]
    assert host.get_next_gen_active(
        0.1,
        0.1,
        10,
        40,
        30,
        new_nt_order,
        probs,
        {1: 'A'},
        30,
        0.99,
        reference_sequence,
        0.2,
        epi,
        0.1,
        rng,
    ) == (2, 2)
    assert "GAA" in [host.C[i].infecting_virus.nuc_sequence for i in range(10)]


def test_summarize_fitness():
    seed(1234)
    new_nt_order, probs = get_nucleotide_substitution_probabilities(
        "../extdata/hiv_q_mat.csv", 3.5e-5
    )
    reference_sequence = "AAA"
    epi = [Epitope(0, 3, 0.3)]
    host = HostEnv([HIV(seq, reference_sequence, 0.1)
                   for seq in ["GGG", "AAA"] * 2], 4)
    assert [host.C[i].infecting_virus.nuc_sequence for i in range(len(host.C))] == [
        "GGG",
        "AAA",
        "GGG",
        "AAA",
    ]
    host.C[1].infecting_virus.mutate(
        1, new_nt_order, probs, {1: 'A'}, 0.99, reference_sequence, 0.1
    )
    assert [
        host.C[i].infecting_virus.replicative_fitness for i in range(len(host.C))
    ] == [(1-0.1)**3, (1-0.1)**1, (1-0.1)**3, 1]
    assert [
        host.C[i].infecting_virus.conserved_fitness for i in range(len(host.C))
    ] == [1, (1-0.99), 1, 1]
    assert [
        host.C[i].infecting_virus.immune_fitness for i in range(len(host.C))
    ] == [1, 1, 1, 1]
    host.get_fitness()
    assert [host.C[i].infecting_virus.fitness for i in range(len(host.C))] == [
        (1-0.1)**3,
        (1 - 0) * (1-0.1)**1 * (1-0.99),
        (1-0.1)**3,
        1.0,
    ]
    assert host.summarize_fitness() == (0.61675, 0.7525, 1.0, 0.8395)


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
        "number_dual_inf": [],
        "mean_fitness_active": [],
        "mean_conserved_active": [],
        "mean_immune_active": [],
        "mean_replicative_active": [],
    }
    host = create_host_env({"founder0": "AAA"}, "AAA", 1, 1)
    assert list(
        host.record_counts(
            counts, 1, (1, 2, 3, 4), (5, 6), (7, 8, 9, 10)).values()) == [
        [1], [1], [0], [1], [2], [3], [4], [5], [6], [7], [8], [9], [10]]
    assert list(
        host.record_counts(
            counts, 2, (11, 12, 13, 14), (15, 16), (17, 18, 19, 20)).values()) == [
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
            10, 20]]


def test_sample_viral_sequences():
    seqs = {"founder0": "AAA"}
    host = create_host_env({"founder0": "AAA"}, "AAA", 1, 1)
    fitness = {
        "generation": [],
        "seq_id": [],
        "immune": [],
        "conserved": [],
        "replicative": [],
        "overall": []
    }
    assert host.sample_viral_sequences(
        seqs,
        fitness,
        1,
        1) == ({'founder0': 'AAA', 'gen1_0': 'AAA'}, {'generation': ['1'], 'seq_id': ['gen1_0'], 'immune': [1.0], 'conserved': [1.0], 'replicative': [1.0], 'overall': [1.0]})


def test_loop_through_generations():
    g = set_python_seed(1)
    new_nt_order, probs = get_nucleotide_substitution_probabilities(
        "../extdata/hiv_q_mat.csv", 3.5e-5
    )
    host = create_host_env({"founder0": "AAA"}, "AAA", 1, 1)
    out = host.loop_through_generations([1, 2, 3],
                                        [1, 2, 3],
                                        2,
                                        {"founder0": "AAA"},
                                        new_nt_order,
                                        probs,
                                        0.1,
                                        0,
                                        0.1,
                                        0,
                                        0,
                                        0,
                                        {1: 'A'},
                                        0.99,
                                        "",
                                        1,
                                        None,
                                        30,
                                        0.1,
                                        90,
                                        g)
    assert out[0] == {
        'generation': [
            0, 1, 2], 'active_cell_count': [
            1, 2, 3], 'latent_cell_count': [
                0, 0, 0], 'active_turned_latent': [
                    0, 0, 0], 'latent_turned_active': [
                        0, 0, 0], 'latent_died': [
                            0, 0, 0], 'latent_proliferated': [
                                0, 0, 0], 'number_mutations': [
                                    0, 1, 2], 'number_dual_inf': [
                                        0, 0, 0], 'mean_fitness_active': [
                                            1.0, 1.0, 1.0], 'mean_conserved_active': [
                                                1.0, 1.0, 1.0], 'mean_immune_active': [
                                                    1.0, 1.0, 1.0], 'mean_replicative_active': [
                                                        1.0, 1.0, 1.0]}
    assert out[1] == {'generation': ['founder', '0', '1', '1', '2', '2', '2'], 'seq_id': ['founder0', 'gen0_0', 'gen1_0', 'gen1_1', 'gen2_0', 'gen2_1', 'gen2_2'], 'immune': [
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 'conserved': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 'replicative': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 'overall': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]}
    assert out[2] == {'founder0': 'AAA', 'gen0_0': 'AAA', 'gen1_0': 'AAG',
                      'gen1_1': 'AAG', 'gen2_0': 'AAG', 'gen2_1': 'AAG', 'gen2_2': 'GAA'}
