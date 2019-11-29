import libsbn
import pandas as pd

data_sets = {
    "andy": "data/DS1.100_topologies.nwk",
    "ds1": "../ds-experiments/DS1/DS1_out.t",
    "ds4": "../ds-experiments/DS4/DS4_out.t",
    "ds7": "../ds-experiments/DS7/DS7_out.t",
    "beast": "_ignore/H1N1_HA.trees",
    "ebola_beast": "_ignore/Ebola_BEAST.trees",
    "ebola_mb": "../../test-libsbn/mb/makona_out.t",
}

alpha = 0.0001
iter_count = 10
score_epsilon = 0
# score_epsilon = 1e-6
# data_set = "ebola_mb"
data_set = "ds7"

inst = libsbn.instance("charlie")

for data_set in ["ds1", "ds4", "ds7"]:
    print(data_set)
    print("reading trees...")
    inst.read_newick_file(data_sets[data_set])
    inst.tree_collection.drop_first(0.25)
    print("processing trees...")
    inst.process_loaded_trees()
    print("training sbn...")
    score_history = pd.Series(
        inst.train_expectation_maximization(alpha, iter_count, score_epsilon)
    )
    prefix = f"_ignore/{data_set}-alpha{alpha}-loop{iter_count}-eps{score_epsilon}"
    score_history.to_csv(prefix + "-score.csv", header=False)
    inst.sbn_parameters.tofile(prefix + "-parameters.csv", ",")

    indexer, parent_to_range = inst.get_indexers()
    single_child_count = 0
    for _, indexer_range in parent_to_range.items():
        (start, end) = indexer_range
        child_count = end - start
        if child_count == 1:
            single_child_count += 1

    print(f"{single_child_count} / {len(parent_to_range)} are single children")
