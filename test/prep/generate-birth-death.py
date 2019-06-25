import dendropy
from dendropy.model import birthdeath

taxon_count = 1000
tree_count = 1000

namespace = dendropy.TaxonNamespace(["x" + str(i) for i in range(taxon_count)])

def gen_bd():
    return dendropy.model.birthdeath.birth_death_tree(birth_rate = 1., death_rate=0., birth_rate_sd=0.0, death_rate_sd=0.0, num_extant_tips=taxon_count, taxon_namespace=namespace)


tl = dendropy.TreeList()

for i in range(tree_count):
    print(i)
    tl.append(gen_bd())

tl.write(path = f"../../_ignore/birth_death_{tree_count}_of_{taxon_count}_taxa.nwk", schema="newick", suppress_rooting=True)
