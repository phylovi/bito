import sbn

def test_instance():
    inst = sbn.instance('charlie')
    inst.parse_file('data/four_taxon.tre')
    inst.print_status()
    assert inst.tree_count() == 5
    inst.init_indexer()
    print(inst.indexer)
    assert inst.indexer[4] == 2
