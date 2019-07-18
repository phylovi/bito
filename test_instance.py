import numpy as np
import sbn

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

    c_v = sbn.make_vector()
    v = np.array(c_v, copy=False)
    v[3] = 666
    print(sbn.total_vector(c_v))

    inst.build_indexer()
    sbn_probs = np.array(inst.sbn_probs, copy=False)
    sbn_probs[3] = 0.4
    print(inst.sbn_total_prob())
