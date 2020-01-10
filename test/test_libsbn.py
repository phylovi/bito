"""Some basic testing and demo code for the libsbn module.

If you want to see the results of the print statements, use `pytest -s`.
"""

import json
import pandas as pd
import pprint
import pytest
import numpy as np
import libsbn
import libsbn.beagle_flags as beagle_flags

SIMPLE_SPECIFICATION = libsbn.PhyloModelSpecification(
    substitution="JC69", site="constant", clock="none"
)


def convert_dict_to_int(dictionary):
    """Change the values of a dict to ints."""
    return {k: int(v) for k, v in dictionary.items()}


def hello_demo():
    """Demonstrate basic phylogenetic likelihood calculation using the "hello"
    data set."""
    inst = libsbn.instance("charlie")
    inst.tree_collection = libsbn.TreeCollection(
        [libsbn.Tree.of_parent_id_vector([3, 3, 3])], ["mars", "saturn", "jupiter"]
    )
    inst.read_fasta_file("data/hello.fasta")
    inst.prepare_for_phylo_likelihood(
        SIMPLE_SPECIFICATION, 2, [beagle_flags.VECTOR_SSE]
    )
    branch_lengths = np.array(inst.tree_collection.trees[0].branch_lengths, copy=False)
    branch_lengths[:] = np.array([0.1, 0.1, 0.3, 0.0])
    print(inst.tree_collection.newick())
    print(np.array(inst.log_likelihoods()))
    branch_lengths[0] = 0.2
    print(inst.tree_collection.newick())
    print(np.array(inst.log_likelihoods()))


def sampling_and_indexers_demo():
    """Demonstrate sampling and indexers.

    Demonstrate loading trees, processing them to make the SBN data
    structures, and then sampling from the SBN with arbitrarily-set
    parameters.
    """
    inst = libsbn.instance("charlie")
    inst.read_newick_file("data/five_taxon.nwk")
    assert inst.tree_count() == 4
    # Showing off tree sampling.
    inst.process_loaded_trees()
    inst.train_expectation_maximization(0.0001, 1)
    # Note that this puts the trees into the instance object, replacing the trees loaded
    # from the file.
    inst.sample_trees(2)
    print("\ntaxon names:")
    print(inst.taxon_names)
    print("\nSBN indexing:")
    print(inst.get_indexer_representations())
    print("\nPSP indexing:")
    print(inst.get_psp_indexer_representations())
    print("\nPSP details:")
    print(inst.psp_indexer.details())
    print("\nSBN parameters:")
    print(inst.sbn_parameters)
    print()


def ds1_support_test():
    """Check the subplit support calculation on DS1."""
    inst = libsbn.instance("DS1")
    # Checking split supports
    inst.read_nexus_file("data/DS1.subsampled_10.t.reordered")
    inst.process_loaded_trees()
    [rootsplit_support, subsplit_support] = inst.split_counters()
    with open("data/DS1.subsampled_10.t_support.json") as json_fp:
        supports = json.load(json_fp)
        vbpi_rootsplit_supp_dict = convert_dict_to_int(supports["rootsplit_supp_dict"])
        vbpi_subsplit_supp_dict = {
            ss: convert_dict_to_int(d)
            for ss, d in supports["subsplit_supp_dict"].items()
        }
    # vbpi and libsbn differ a little concerning how they compute the values of
    # the subsplit support dictionaries. However, we primarily care about the
    # actual support, so that's compared here using a call to keys.
    assert rootsplit_support.keys() == vbpi_rootsplit_supp_dict.keys()
    assert subsplit_support.keys() == vbpi_subsplit_supp_dict.keys()
    return inst


def ds1_phylo_model_demo(inst):
    """Demonstrate how phylogenetic models and likelihoods work using DS1."""
    inst.read_fasta_file("data/DS1.fasta")
    # Just use the first tree for likelihood comparison.
    inst.tree_collection.erase(1, 10)
    branch_lengths = np.array(inst.tree_collection.trees[0].branch_lengths, copy=False)
    branch_lengths[:] = 0.1

    inst.prepare_for_phylo_likelihood(SIMPLE_SPECIFICATION, 2)
    jc69_likelihood = np.array(inst.log_likelihoods())

    # Showing off phylo_model_param_block_map.
    gtr_specification = libsbn.PhyloModelSpecification(
        substitution="GTR", site="constant", clock="none"
    )
    inst.prepare_for_phylo_likelihood(gtr_specification, 2)
    phylo_model_param_block_map = inst.get_phylo_model_param_block_map()
    phylo_model_param_block_map["GTR rates"][:] = 1.0
    phylo_model_param_block_map["frequencies"][:] = 0.25
    print("\nHere's a look at phylo_model_param_block_map:")
    pprint.pprint(phylo_model_param_block_map)
    print("\nWe can see that we are changing the phylo_model_params matrix:")
    print(inst.get_phylo_model_params(), "\n")
    assert jc69_likelihood == pytest.approx(np.array(inst.log_likelihoods()))


def rootings_indexer_test():
    """Ensure that various rootings of a given tree give the same indexer
    representations in a mathematical sense.

    In order to compare them, we sort their representations with respect
    to the order on the rootsplits. This makes sense because we want to
    make sure that for each virtual rooting (corresponding to each
    rootsplit) we have the same _set_ of PCSSs (the order within PCSS
    sets doesn't matter).
    """
    inst = libsbn.instance("rootings")
    inst.read_newick_file("data/many_rootings.nwk")
    inst.process_loaded_trees()
    # First we turn the PCSS sets into actual Python sets for unordered comparison.
    reps = [
        (rootsplits, [set(pcss_set) for pcss_set in pcss_set_list])
        for (rootsplits, pcss_set_list) in inst.get_indexer_representations()
    ]
    # Next we sort the representations with respect to the order on the rootsplits
    # (the first component, corresponding to the various virtual rootings).
    final_sorted = [zip(*sorted(zip(*rep))) for rep in reps]
    first_sorted_rep = list(final_sorted[0])
    for rep in final_sorted[1:]:
        assert first_sorted_rep == list(rep)


def sbn_training_test():
    """Test SBN training."""
    inst = libsbn.instance("sbn")
    inst.read_newick_file("data/DS1.100_topologies.nwk")
    inst.process_loaded_trees()
    inst.train_simple_average()
    probabilities = inst.calculate_sbn_probabilities()
    pd.Series(probabilities).to_csv(
        "_ignore/simple-average.csv", index=False, header=False
    )
    inst.train_expectation_maximization(0.0, 1)
    probabilities = inst.calculate_sbn_probabilities()
    pd.Series(probabilities).to_csv(
        "_ignore/em-alpha0-loop1.csv", index=False, header=False
    )
    inst.train_expectation_maximization(0.0, 23)
    probabilities = inst.calculate_sbn_probabilities()
    pd.Series(probabilities).to_csv(
        "_ignore/em-alpha0-loop23.csv", index=False, header=False
    )
    inst.train_expectation_maximization(0.3, 7)
    probabilities = inst.calculate_sbn_probabilities()
    pd.Series(probabilities).to_csv(
        "_ignore/em-alpha0.3-loop7.csv", index=False, header=False
    )
    inst.train_expectation_maximization(0.3, 1)
    probabilities = inst.calculate_sbn_probabilities()
    pd.Series(probabilities).to_csv(
        "_ignore/em-alpha0.3-loop1.csv", index=False, header=False
    )


def test_libsbn():
    """Test the libsbn instance."""

    sbn_training_test()
    hello_demo()
    sampling_and_indexers_demo()
    inst = ds1_support_test()
    ds1_phylo_model_demo(inst)
    rootings_indexer_test()
