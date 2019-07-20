import numpy as np
from utils import get_tree_list_raw, get_support_from_samples
import json
import os

def get_supports(path):
    tree_dict_total, tree_names_total, tree_wts_total = get_tree_list_raw(path, hpd=1.)
    taxa = sorted([leaf.name for leaf in tree_dict_total['tree_1'].iter_leaves()])
    rootsplit_supp_dict, subsplit_supp_dict = get_support_from_samples(taxa, tree_dict_total, tree_names_total)

    out_path = os.path.splitext(path)[0] + "_support.json"
    with open(out_path, 'w') as fp:
        fp.write(json.dumps({"rootsplit_supp_dict": rootsplit_supp_dict, "subsplit_supp_dict": subsplit_supp_dict}))


get_supports('data/five_taxon.nwk')
get_supports('data/DS1.subsampled_10.t.nwk')
