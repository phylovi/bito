import numpy as np
from utils import get_tree_list_raw, get_support_from_samples
import json

tree_dict_total, tree_names_total, tree_wts_total = get_tree_list_raw('data/five_taxon.nwk', hpd=1.)
taxa = sorted([leaf.name for leaf in tree_dict_total['tree_3'].iter_leaves()])
rootsplit_supp_dict, subsplit_supp_dict = get_support_from_samples(taxa, tree_dict_total, tree_names_total)

with open('data/five_taxon_support.json', 'w') as fp:
    fp.write(json.dumps({"rootsplit_supp_dict": rootsplit_supp_dict, "subsplit_supp_dict": subsplit_supp_dict}))
