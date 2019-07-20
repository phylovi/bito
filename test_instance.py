import json
import numpy as np
import sbn


def test_instance():
    inst = sbn.instance('charlie')
    inst.read_newick_file('data/five_taxon.nwk')
    inst.print_status()
    assert inst.tree_count() == 4
    inst.process_loaded_trees()

#    [indexer, range_indexer, index_to_child_str] = inst.get_indexers()
#    print(indexer)
#    print(range_indexer)
#
#    sbn_probs = np.array(inst.sbn_probs, copy=False)
#    sbn_probs[3] = 3.14159265359
#    print(sbn_probs)
#    print(inst.sbn_total_prob())


    def convert_dict_to_int(d):
        return {k:int(v) for k, v in d.items()}

    [rootsplit_support, subsplit_support] = inst.split_counters()
    with open('data/five_taxon_support.json') as fp:
        supports = json.load(fp)
        vbpi_rootsplit_supp_dict = convert_dict_to_int(supports["rootsplit_supp_dict"])
        vbpi_subsplit_supp_dict = {ss:convert_dict_to_int(d) for ss, d in supports["subsplit_supp_dict"].items()}
    assert rootsplit_support == vbpi_rootsplit_supp_dict
    assert subsplit_support == vbpi_subsplit_supp_dict

#    inst.read_nexus_file('data/DS1.subsampled_10.t')
#    inst.read_fasta_file('data/DS1.fasta')
#    inst.make_beagle_instances(2)
#    print(np.array(inst.tree_log_likelihoods()))
