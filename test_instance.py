import json
import numpy as np
import sbn


def test_instance():
    inst = sbn.instance('charlie')
    inst.read_newick_file('data/five_taxon.nwk')
    inst.print_status()
    assert inst.tree_count() == 4
    inst.process_loaded_trees()

    [indexer, range_indexer] = inst.get_indexers()

    sbn_probs = np.array(inst.sbn_probs, copy=False)
    sbn_probs[3] = 3.14159265359

    # print(sbn_probs)
    # print(inst.sbn_total_prob())

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
    gradients = [np.array(gradient) for gradient in inst.branch_gradients()]
    print(gradients[-1])

#    t = sbn.Tree.of_index_vector([5, 5, 4, 5, 5])
    t = sbn.Tree.of_index_vector([6, 5, 4, 4, 5, 6])
    branch_lengths = np.array(t, copy=False)
    branch_lengths[0] = 3.14159265359
#    print(t.newick())


