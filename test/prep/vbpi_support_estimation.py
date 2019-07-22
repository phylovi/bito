import numpy as np
from utils import get_tree_list_raw, get_support_from_samples
import json
import os

five_taxon_taxa_str = "x0 x1 x2 x3 x4"
DS1_taxa_str = """
Rattus_norvegicus Mus_musculus Ichthyophis_bannanicus Turdus_migratorius Latimeria_chalumnae Heterodon_platyrhinos Discoglossus_pictus Gallus_gallus Trachemys_scripta Grandisonia_alternans Oryctolagus_cuniculus Nesomantis_thomasseti Gastrophryne_carolinensis Siren_intermedia Homo_sapiens Sceloporus_undulatus Eleutherodactylus_cuneatus Typhlonectes_natans Plethodon_yonhalossee Xenopus_laevis Hypogeophis_rostratus Hyla_cinerea Bufo_valliceps Scaphiopus_holbrooki Ambystoma_mexicanum Amphiuma_tridactylum Alligator_mississippiensis
"""

def get_supports(path, taxa):
    tree_dict_total, tree_names_total, tree_wts_total = get_tree_list_raw(path, hpd=1.)
    rootsplit_supp_dict, subsplit_supp_dict = get_support_from_samples(taxa, tree_dict_total, tree_names_total)

    out_path = os.path.splitext(path)[0] + "_support.json"
    with open(out_path, 'w') as fp:
        fp.write(json.dumps({"rootsplit_supp_dict": rootsplit_supp_dict, "subsplit_supp_dict": subsplit_supp_dict}))


get_supports('data/five_taxon.nwk', five_taxon_taxa_str.split())
get_supports('data/DS1.subsampled_10.t.nwk', DS1_taxa_str.split())
