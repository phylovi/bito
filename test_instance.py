import numpy as np
import sbn

def test_instance():
    inst = sbn.instance('charlie')
    inst.parse_file('data/five_taxon.tre')
    inst.print_status()
    # assert inst.tree_count() == 3
    sbn.f(np.array([3,4]))
    supports = inst.supports()
    supports_list = list(supports.keys())
    supports_list.sort()
    with open('_build/erick-support.txt', 'w') as fp:
        for support in supports_list:
            fp.write(support+'\n')
