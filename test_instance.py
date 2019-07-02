import numpy as np
import sbn

def test_instance():
    inst = sbn.instance('charlie')
    inst.parse_file('data/four_taxon.tre')
    inst.print_status()
    # There are 5 input trees, but 2 pairs are identical.
    assert inst.tree_count() == 3
    inst.parse_file('data/one_four_taxon.tre')
    assert inst.tree_count() == 1
    sbn.f(np.array([3,4]))
    s = inst.supports()
    for k,v in s.items():
        print(k)
