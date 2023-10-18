"""
Tests that NNI search finds NNIs in the correct order.
"""

import os
from nni_search import Program, Loader
import bito

golden_git_commit = "811b735"

script_path = os.path.dirname(__file__)
data_path = os.path.join(script_path, "../data/ds1/")
fasta_path = os.path.join(data_path, "ds1.fasta")
nwk_seed_path = os.path.join(data_path, "ds1.top1.nwk")
nwk_cred_path = os.path.join(data_path, "ds1.credible.with-branches.rerooted.nwk")
tree_pp_path = os.path.join(data_path, "ds1.mb-pp.csv")
pcsp_pp_path = os.path.join(data_path, "ds1.pcsp-pp.csv")
golden_run_path = os.path.join(data_path, "test", f"run.{golden_git_commit}.csv")
iter_max = 200  # should only take 104 iterations to find all credible edges
cmd_line_args = f"nni-search {fasta_path} {nwk_seed_path} {nwk_cred_path} {tree_pp_path} {pcsp_pp_path} --tp --iter-max {iter_max}"
cmd_line_args = cmd_line_args.split()

parsed_args = Program.parse_args(cmd_line_args)
nni_inst, results, timer = Program.run_program(parsed_args)

nni_list = results.data_['nni_hash']
golden_git_commit, golden_nni_list = Loader.load_nni_list(golden_run_path)

# compare nni_list result from run against golden run
print("\n=== COMPARE_RESULTS_TO_GOLDEN ===")
print(
    f"# GIT_COMMIT: TEST::{bito.git_commit()} GOLDEN::{golden_git_commit}")
test_passes = (nni_list == golden_nni_list)
if test_passes:
    print("# TEST_PASSED -- accepted NNIs matches golden run.")
else:
    print("# TEST_FAILED -- accepted NNIs do not match golden run.")
if not test_passes:
    print(f"# ITERS: TEST::{len(nni_list)} GOLDEN::{len(golden_nni_list)}")
    print(
        f"# {'ITER':<7} {'RESULT':<10} {'TEST_NNI':<20} {'GOLDEN_NNI':<20}")
    for i in range(max(len(golden_nni_list), len(nni_list))):
        nni = "None"
        if (i < len(nni_list)):
            nni = nni_list[i]
        golden_nni = "None"
        if (i < len(golden_nni_list)):
            golden_nni = golden_nni_list[i]
        if (nni == golden_nni):
            res = "pass"
        else:
            res = "fail"
        print(f"  {i:<7} {res:<10} {nni:<20} {golden_nni:<20}")
