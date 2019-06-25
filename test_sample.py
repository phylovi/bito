import sbn

def test_instance():
    inst = sbn.instance('charlie')
    inst.parse_file('data/four_taxon.tre')
    inst.print_status()
    assert 5 == 5
