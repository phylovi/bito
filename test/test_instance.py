import json
import numpy as np
import libsbn


def test_instance():
    """Test the libsbn instance.

    This is also a bit of a demo. If you want to see the results of the
    print statements, use `pytest -s`.
    """
    inst = libsbn.instance("charlie")
    inst.read_newick_file("data/five_taxon.nwk")
    assert inst.tree_count() == 4
    inst.process_loaded_trees()

    # Showing off tree sampling.
    [indexer, range_indexer] = inst.get_indexers()
    sbn_parameters = np.array(inst.sbn_parameters, copy=False)
    sbn_parameters[0] = 0.2
    # Note that this puts the trees into the instance object, replacing the trees loaded
    # from the file.
    inst.sample_trees(2)
    print("\nSBN indexing:")
    print(inst.get_indexer_representations())
    print("\nPSP indexing:")
    print(inst.get_psp_indexer_representations())
    print("\nPSP details:")
    print(inst.psp_indexer.details())
    print()

    # Checking split supports.
    def convert_dict_to_int(d):
        return {k: int(v) for k, v in d.items()}

    inst.read_nexus_file("data/DS1.subsampled_10.t")
    inst.process_loaded_trees()
    [rootsplit_support, subsplit_support] = inst.split_counters()
    with open("data/DS1.subsampled_10.t_support.json") as fp:
        supports = json.load(fp)
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

    inst.read_fasta_file("data/DS1.fasta")
    inst.sample_trees(1)
    branch_lengths = np.array(inst.tree_collection.trees[0].branch_lengths, copy=False)
    branch_lengths[:] = 0.1
    gtr_specification = libsbn.PhyloModelSpecification(
        substitution="GTR", site="constant", clock="strict"
    )
    inst.make_engine(gtr_specification, 2)
    parametrization = {
        "frequencies": np.array([0.2, 0.55, 0.1, 0.15]),
        "rates": np.array([1.0, 0.5, 1.0, 2.0, 1.0, 1.0]),
    }
    print("with parameterization:", np.array(inst.log_likelihoods([parametrization])))

    inst.tree_collection = libsbn.TreeCollection(
        [libsbn.Tree.of_parent_id_vector([3, 3, 3])], ["mars", "saturn", "jupiter"]
    )
    inst.read_fasta_file("data/hello.fasta")
    simple_specification = libsbn.PhyloModelSpecification(
        substitution="JC69", site="constant", clock="strict"
    )
    inst.make_engine(simple_specification, 2)
    branch_lengths = np.array(inst.tree_collection.trees[0].branch_lengths, copy=False)
    branch_lengths[:] = np.array([0.1, 0.1, 0.3, 0.0])
    print(inst.tree_collection.newick())
    print(np.array(inst.log_likelihoods([])))
    branch_lengths[0] = 0.2
    print(inst.tree_collection.newick())
    print(np.array(inst.log_likelihoods()))

    # Here we ensure that various rootings of a given tree give the same indexer
    # representations in a mathematical sense.
    # In order to compare them, we sort their representations with respect to the order on
    # the rootsplits. This makes sense because we want to make sure that for each virtual
    # rooting (corresponding to each rootsplit) we have the same _set_ of PCSSs (the
    # order within PCSS sets doesn't matter).
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
