"""Some basic testing and demo code for the bito module.

If you want to see the results of the print statements, use `pytest -s`.
"""

import os
import json
import pprint
import pytest
import numpy as np
import bito
import bito.beagle_flags as beagle_flags
import bito.phylo_model_mapkeys as model_keys


def nni_search(fasta_file: text, newick_file: text):
    inst = bito.gp_instance()
    inst.read_newick_file(newick_file)
    inst.read_fasta_file(fasta_file)

    # initialize engines
    inst.make_engine("_ignore/gp_engine.data")
    inst.make_tp_engine("_ignore/tp_like.data", True, "_ignore/tp_pars.data")
    inst.make_nni_engine()
    nni_engine = inst.get_nni_engine()

    # run
    inst.sync_adjacent_nnis_with_dag()
    iter_count = 0
    while inst.adjacent_nni_count() > 0:
        for i in range(len(inst.adjacent_nni_count())):

            iter_count += 1

    pass


if __name__ == "__main__":
    print("Calling NNI Search from main...")

    if len(sys.argv) != 3:
        print("Usage: <fasta_file> <newick_file>")
        exit(1)

    fasta_file = sys.argv[1]
    newick_file = sys.argv[2]
    nni_search(fasta_file, newick_file)

