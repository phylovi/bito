import numpy as np
import sbn

import timeit

def test_instance():
    inst = sbn.instance('charlie')
    inst.read_newick_file('data/five_taxon.nwk')
    inst.print_status()
    assert inst.tree_count() == 4
    [rootsplit_support, subsplit_support] = inst.split_supports()

    with open('_build/support.txt', 'w') as fp:
        for support in [rootsplit_support, subsplit_support]:
            support_list = list(support.keys())
            support_list.sort()
            for support in support_list:
                fp.write(support+'\n')
            fp.write('\n')

    inst.read_nexus_file('data/DS1.subsampled_10.t')
    inst.read_fasta_file('data/DS1.fasta')
    inst.make_beagle_instances(2)
    print(inst.tree_log_likelihoods())

    v1_c = sbn.make_vector()
    v2_c = sbn.make_vector()
    v3_c = sbn.make_vector()
    a1_c = np.array(v1_c, copy=False)
    a2_c = np.array(v2_c, copy=False)
    a3_c = np.array(v2_c, copy=False)
