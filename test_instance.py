import numpy as np
import sbn

def test_instance():
    inst = sbn.instance('charlie')
    inst.parse_file('data/four_taxon.tre')
    inst.print_status()
    # There are 5 input trees, but 2 pairs are identical.
    assert inst.tree_count() == 3
    sbn.f(np.array([3,4]))
    supports = inst.supports()
    supports_list = list(supports.keys())
    supports_list.sort()
    with open('/home/ematsen/erick-support.txt', 'w') as fp:
        for support in supports_list:
            fp.write(support+'\n')
