"""
Tests that NNI search finds NNIs in the correct order.
"""

from nni_search import Program, Loader
import bito

data_path = "../data/ds1/"
fasta_path = "../data/ds1/ds1.fasta"
nwk_seed_path = "../data/ds1/ds1.top1.nwk"
nwk_cred_path = "../data/ds1/ds1.credible.with-branches.rerooted.nwk"
tree_pp_path = "../data/ds1/ds1.mb-pp.csv"
pcsp_pp_path = "../data/ds1/ds1.pcsp-pp.csv"
golden_nni_list_path = "../data/ds1/golden-results.csv"
cmd_line_args = f"nni-search {fasta_path} {nwk_seed_path} {nwk_cred_path} {tree_pp_path} {pcsp_pp_path} --tp --iter-max 200"
cmd_line_args = cmd_line_args.split()
print(f"cmd_line_args: {cmd_line_args}")

parsed_args = Program.parse_args(cmd_line_args)
nni_inst, results, times = Program.run_program(parsed_args)

nni_list = nni_inst.data_['nni_hash']
golden_git_commit, golden_nni_list = Loader.load_nni_list(golden_nni_list_path)

# compare nni_lists
print("\n=== TEST ===")
print(
    f"# GIT_COMMIT: TEST::{bito.git_commit()} GOLDEN::{golden_git_commit}")
if (nni_list == golden_nni_list):
    print("# TEST_PASSED -- accepted NNIs matches golden run.")
else:
    print("# TEST_FAILED -- accepted NNIs do not match golden run.")

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
