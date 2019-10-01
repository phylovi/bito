import json
import numpy as np
import libsbn


def test_instance():
    inst = libsbn.instance('charlie')
    inst.read_newick_file('data/five_taxon.nwk')
    assert inst.tree_count() == 4
    inst.process_loaded_trees()

    # Showing off tree sampling.
    [indexer, range_indexer] = inst.get_indexers()
    sbn_parameters = np.array(inst.sbn_parameters, copy=False)
    sbn_parameters[0] = 0.2
    # Note that this puts the trees into the instance object, replacing the trees loaded from the file.
    inst.sample_trees(2)
    print("\nSBN indexing:")
    print(inst.get_indexer_representations())
    print("\nPSP indexing:")
    print(inst.get_psp_indexer_representations())
    print("\nPSP details:")
    print(inst.psp_indexer.details())
    print()

    # Checking split supports
    def convert_dict_to_int(d):
        return {k: int(v) for k, v in d.items()}

    inst.read_nexus_file('data/DS1.subsampled_10.t')
    inst.process_loaded_trees()
    [rootsplit_support, subsplit_support] = inst.split_counters()
    with open('data/DS1.subsampled_10.t_support.json') as fp:
        supports = json.load(fp)
        vbpi_rootsplit_supp_dict = convert_dict_to_int(
            supports["rootsplit_supp_dict"])
        vbpi_subsplit_supp_dict = {
            ss: convert_dict_to_int(d)
            for ss, d in supports["subsplit_supp_dict"].items()
        }
    # vbpi and libsbn differ a little concerning how they compute the values of
    # the subsplit support dictionaries. However, we primarily care about the
    # actual support, so that's compared here using a call to keys.
    assert rootsplit_support.keys() == vbpi_rootsplit_supp_dict.keys()
    assert subsplit_support.keys() == vbpi_subsplit_supp_dict.keys()

    inst.read_fasta_file('data/DS1.fasta')
    inst.make_beagle_instances(2)
    print(np.array(inst.log_likelihoods()))

    log_likelihoods, gradients = zip(*inst.branch_gradients())
    print(np.array(log_likelihoods))
    print(np.array(gradients[-1]))

    inst.tree_collection = libsbn.TreeCollection(
        [libsbn.Tree.of_parent_id_vector([3, 3, 3])],
        ["mars", "saturn", "jupiter"])
    inst.read_fasta_file('data/hello.fasta')
    inst.make_beagle_instances(2)
    branch_lengths = np.array(inst.tree_collection.trees[0].branch_lengths,
                              copy=False)
    branch_lengths[:] = np.array([0.1, 0.1, 0.3, 0.])
    print(inst.tree_collection.newick())
    print(np.array(inst.log_likelihoods()))
    branch_lengths[0] = 0.2
    print(inst.tree_collection.newick())
    print(np.array(inst.log_likelihoods()))
