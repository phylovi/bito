"""

"""

from __future__ import print_function
import sys
import argparse
import numpy as np
import scipy.stats as ss
import pandas as pd
import cProfile
import random
import time
import pprint
import bito

### Manual Settings (some can be overridden by commandline arguments)
verbose = 1
# abbreviate hash length for readability (None to not truncate)
abbr = 5
sort_taxa = False
# rounding digits in printed output
digits = 5
# number of optimization iterations
# whether to rescore or reevaluate rejected nnis (use None for default behavior)
do_rescore_all_nnis = None
do_reeval_all_nnis = None
# nni evaluation settings
optimization_max_iteration = 5
do_optimize_new_edges = True
do_use_best_edge_map = True
do_init_proposed_bls_with_dag = True
do_fix_proposed_bls_from_dag = True
# terminate search when all credible edges found
do_end_when_all_creds_found = False
# print data
do_print_settings = True
do_print_setup_data = False
do_print_iter_data = True
do_print_times = True
do_print_dag_stats = True
do_print_scored_nnis = False
do_print_accepted_nnis = True
do_print_summary = True
# diagnostics: track changes to the DAG
do_run_tracker = False
do_print_tracker_summary = False
do_check_for_dag_changes = True
do_check_for_nni_score_changes = True
do_check_for_choice_map_changes = True
do_check_for_pcsp_map_changes = True
do_check_for_nni_map_changes = True
do_check_for_pv_map_changes = True
do_check_for_bl_map_changes = True
do_check_for_dag_score_changes = True

# specific subsplits to watch.
subsplit_watchlist = {}


def print_v(*args, v=verbose):
    ''' Optionally prints based on verbosity argument '''
    if verbose >= v:
        print(*args)


class Timer:
    default_print = True

    def __init__(self):
        self.times_prv = {}
        self.times = {}
        self.start()

    def start(self, timer_name="__default__"):
        self.times_prv[timer_name] = time.time()

    def stop(self, timer_name="__default__"):
        self.times_prv[timer_name] = self.elapsed(timer_name)

    def resume(self, timer_name="__default__"):
        self.times_prv[timer_name] = time.time() - self.times_prv[timer_name]

    def elapsed(self, timer_name="__default__"):
        return round((time.time() - self.times_prv[timer_name]), digits)

    def lap(self, timer_name="__default__", timer_prv_name="__default__"):
        time_elapsed = self.elapsed(timer_prv_name)
        if (timer_name not in self.times):
            self.times[timer_name] = []
        self.times[timer_name].append(time_elapsed)
        self.times_prv[timer_prv_name] = time.time()
        if do_print_times:
            print_v(f'  #TIME {timer_name}: {time_elapsed}')
        return time_elapsed

    def lap_next(self, timer_name="__default__"):
        return self.lap(timer_name, "__default__")

    def lap_name(self, timer_name="__default__"):
        return self.lap(timer_name, timer_name)

    def lap_last(self, timer_name="__default__"):
        return self.times[timer_name][len(self.times[timer_name]) - 1]

    def get_lap(self, timer_name="__default__", lap_count=0):
        return self.times[timer_name][lap_count]

    def to_summary(self):
        times_perc = {}
        for key in self.times:
            times_perc[key] = []
            for time, total in zip(self.times[key], self.times["full_iter"]):
                times_perc[key].append(round(time / total, digits))

        print_v("\n  == TIMES_REAL ==")
        for key in self.times:
            bin_len = min(len(self.times[key]), 10)
            print_v(
                f'  {key}:\n    {Utils.average_over_ranges(self.times[key], bin_len)}')
            x_s = range(len(self.times[key]))
            fit = np.polyfit(y=self.times[key], x=x_s, deg=2)
            print_v(
                f'    #POLYFIT: {round(fit[0], digits)}x^2 + {round(fit[1], digits)}x + {round(fit[2], digits)}')
            print_v(
                f'    #POLYFIT: {round(fit[0]/abs(fit[2]), digits)}x^2 + {round(fit[1]/abs(fit[2]), digits)}x + {round(fit[2]/fit[2], digits)}')

        print_v("\n  == TIMES_PERC ==")
        for key in times_perc:
            bin_len = min(len(self.times[key]), 10)
            print_v(
                f'  {key}:\n    {Utils.average_over_ranges(times_perc[key], bin_len)}')
            x_s = range(len(times_perc[key]))
            fit = np.polyfit(y=times_perc[key], x=x_s, deg=2)
            print_v(
                f'    #POLYFIT: {round(fit[0], digits)}x^2 + {round(fit[1], digits)}x + {round(fit[2], digits)}')
        pass


class Utils:
    @staticmethod
    def build_ranked_list(list):
        total_list = [len(list)] * len(list)
        rank_list = total_list - ss.rankdata(list, method='ordinal')
        return rank_list

    @staticmethod
    def average_over_ranges(data, buckets):
        step = int(len(data) / buckets)
        averages = []
        for i in range(0, len(data), step):
            start = i
            end = min(i + step, len(data))
            sublist = data[start:end]
            if sublist:  # Check if the sublist is not empty
                avg = sum(sublist) / len(sublist)
                averages.append(round(avg, digits))
        return averages

    @staticmethod
    def to_hash(obj, abbr=None):
        hash = f'0x{obj.__hash__():016X}'
        if (abbr != None):
            hash = hash[:abbr + 2]
        return hash

    @staticmethod
    def bitset_to_hash(bitset, abbr=16):
        return bitset.to_hash_string(abbr)

    @staticmethod
    def subsplit_to_hash(subsplit, abbr=16):
        # return subsplit.subsplit_to_hash_string(abbr)
        return subsplit.to_hash_string(abbr)

    @staticmethod
    def pcsp_to_hash(pcsp, abbr=16):
        # return pcsp.pcsp_to_hash_string(abbr)
        return pcsp.to_hash_string(abbr)

    @staticmethod
    def nni_to_hash(nni, abbr=16):
        return nni.to_hash_string(abbr)
        # return nni.get_central_edge_pcsp().pcsp_to_hash_string(abbr)


class Loader:
    @staticmethod
    def load_trprobs_and_t(trprobs_path, t_path, args):
        taxon_id_map = {}
        tree_nwk_map = {}
        tree_pp_map = {}
        tree_cpp_map = {}
        trprobs_fp = open(trprobs_path, 'r')
        tree_id = 1
        for line in trprobs_fp.readlines():
            fields = line.strip().split()
            if len(fields) == 0:
                continue
            if fields[0].isdigit():
                taxon_id = int(fields[0])
                taxon_name = fields[1]
                taxon_id_map[taxon_id] = taxon_name
            if fields[0] == "tree":
                # tree_id = int(fields[1].split("_")[1])
                tree_pp = float(fields[10].replace("]", ""))
                tree_cpp = float(fields[7].replace("]", ""))
                tree_nwk = fields[len(fields) - 1]
                tree_nwk_map[tree_id] = tree_nwk
                tree_pp_map[tree_id] = tree_pp
                tree_cpp_map[tree_id] = tree_cpp
                tree_id += 1
        trprobs_fp.close()
        return taxon_id_map, tree_nwk_map, tree_pp_map, tree_cpp_map

    @staticmethod
    def load_t(t_path, args):
        taxon_id_map = {}
        tree_nwk_map = {}
        t_fp = open(t_path, 'r')
        tree_id = 1
        for line in t_fp.readlines():
            fields = line.strip().split()
            if fields[0].isdigit():
                taxon_id = int(fields[0])
                taxon_name = fields[1]
                taxon_id_map[taxon_id] = taxon_name
            if fields[0] == "tree":
                tree_nwk = fields[len(fields) - 1]
                tree_nwk_map[tree_id] = tree_nwk
                tree_id += 1
        return taxon_id_map, tree_nwk_map

    @staticmethod
    def load_pcsp_weight_table(pcsp_weight_path):
        pcsp_weight_path = sys.argv[1]
        pcsp_weight_header = ["pcsp", "pp"]
        pcsp_wt_df = pd.read_csv(pcsp_weight_path, names=pcsp_weight_header)
        pcsp_wt_df = pcsp_wt_df[1:]

        parent_subsplits = []
        child_subsplits = []
        for pcsp in pcsp_wt_df["pcsp"]:
            pcsp_bitset = bito.bitset(pcsp)
            parent_subsplits.append(
                pcsp_bitset.pcsp_get_parent_subsplit().subsplit_to_string())
            child_subsplits.append(
                pcsp_bitset.pcsp_get_child_subsplit().subsplit_to_string())
        pcsp_wt_df["parent"] = parent_subsplits
        pcsp_wt_df["child"] = child_subsplits
        pcsp_wt_df = pcsp_wt_df.drop(["pcsp"], axis=1)
        return pcsp_wt_df

    @staticmethod
    def load_trees(fasta_path, newick_path):
        tree_inst = bito.rooted_instance("trees")
        tree_inst.read_fasta_file(fasta_path)
        # tree_inst.read_newick_file(newick_path)
        tree_inst.read_newick_file(newick_path, sort_taxa)
        trees = tree_inst.tree_collection.trees
        return tree_inst, trees

    @staticmethod
    def load_dag(fasta_path, newick_path):
        rand = random.randint(10000, 99999)
        dag_inst = bito.gp_instance(f"_ignore/mmap.${rand}.data")
        dag_inst.read_fasta_file(fasta_path)
        # dag_inst.read_newick_file(newick_path)
        dag_inst.read_newick_file(newick_path, sort_taxa)
        dag_inst.make_dag()
        dag = dag_inst.get_dag()
        return dag_inst, dag

    @staticmethod
    def load_pps(pp_path):
        pps = []
        with open(pp_path, 'r') as fp:
            for line in fp.readlines():
                pps.append(float(line))
        return pps

    @staticmethod
    def load_pcsp_pp_map(pcsp_pp_path):
        pcsp_pp_map = {}
        df = pd.read_csv(pcsp_pp_path)
        for index, row in df.iterrows():
            parent = row["parent"].split("|")
            parent = bito.subsplit(parent[0], parent[1])
            child = row["child"].split("|")
            child = bito.subsplit(child[0], child[1])
            pcsp = bito.pcsp(parent, child)
            pp = float(row["pcsp_pp"])
            pcsp_pp_map[pcsp] = pp
        return pcsp_pp_map

    @staticmethod
    def build_tree_id_map(trees):
        tree_id_map = {}
        for id, tree in enumerate(trees):
            tree_id_map[id] = tree
        return tree_id_map

    @staticmethod
    def build_tree_pp_map(tree_id_map, pps):
        tree_pp_map = {}
        for (tree_id, pp) in zip(tree_id_map, pps):
            tree_pp_map[tree_id] = pp
        return tree_pp_map

    @staticmethod
    def build_pcsp_pp_map(dag, tree_id_map, tree_pp_map):
        dag_pcsps = dag.build_set_of_edge_bitsets()
        pcsp_pp_map = {}
        for pcsp_count, pcsp in enumerate(dag_pcsps):
            print_v(f"# loading pcsp {pcsp_count} of {len(dag_pcsps)}...")
            pcsp_pp_map[pcsp] = 0.0
        for tree_id in tree_pp_map:
            tree = tree_id_map[tree_id]
            pp = tree_pp_map[tree_id]
            tree_pcsps = tree.build_set_of_pcsps()
            tree_pcsps = tree.build_set_of_pcsps()
            for pcsp in tree_pcsps:
                pcsp_pp_map[pcsp] += pp
        return pcsp_pp_map

    @staticmethod
    def load_nni_list(nni_list_path):
        git_commit = None
        nni_list = []
        nni_list_fp = open(nni_list_path, 'r')
        for line in nni_list_fp.readlines():
            fields = line.strip().split()
            if not fields[0].startswith("#"):
                nni_list.append(fields[0])
            elif fields[0].startswith("#GIT_COMMIT"):
                git_commit = fields[1]
        nni_list_fp.close()
        return git_commit, nni_list


class Results:
    def __init__(self, args):
        self.data_ = {}
        self.pre_data_ = {}
        self.df_ = None
        self.args = args
        pass

    def data_init(self):
        self.data_ = pd.DataFrame(columns=['iter', 'acc_nni_id', 'acc_nni_count', 'score', 'is_nni_cred', 'is_nni_new', 'tree_pp', 'pcsp_pp', 'pcsp_pp_rank', 'node_count', 'edge_count', 'cred_edge_count', 'tree_count', 'adj_nni_count', 'new_adj_nni_count', 'cred_adj_nni_count', 'llhs_computed', 'parent', 'child', 'pcsp', 'nni_hash'])
        pass

    def predata_init(self, dag, nni_engine, pp_maps):
        self.pre_data_['iter_count'] = 0
        self.pre_data_['llhs_computed'] = 0
        self.pre_data_['prev_adj_nni_count'] = 0
        self.pre_data_['tree_pp_total'] = pp_maps.get_tree_pp_total()
        self.pre_data_['cred_edge_total'] = pp_maps.get_credible_edge_total()
        pass

    def predata_begin_iter(self, iter_count, dag, nni_engine, pp_maps):
        self.pre_data_['iter_count'] = iter_count
        self.pre_data_['tree_pp'] = pp_maps.get_tree_pp(dag)
        pass

    def predata_mid_iter(self, iter_count, dag, nni_engine, pp_maps):
        new_adj_nni_count = len(
            nni_engine.adjacent_nnis()) - self.pre_data_['prev_adj_nni_count']
        self.pre_data_['new_adj_nni_count'] = new_adj_nni_count
        prev_adj_nni_count = len(nni_engine.adjacent_nnis())
        self.pre_data_['prev_adj_nni_count'] = prev_adj_nni_count
        llhs_computed = self.pre_data_['llhs_computed']
        llhs_computed += new_adj_nni_count
        self.pre_data_['llhs_computed'] = llhs_computed
        pass

    def predata_end_iter(self, iter_count, dag, nni_engine, pp_maps):
        best_nni_score = -np.inf
        scored_nnis = nni_engine.scored_nnis()
        for nni in nni_engine.accepted_nnis():
            best_nni_score = scored_nnis[nni]
            break
        self.pre_data_['best_nni_score'] = best_nni_score
        cred_edge_count, noncred_edge_count = pp_maps.get_credible_edge_count(
            dag)
        self.pre_data_['cred_edge_count'] = cred_edge_count
        self.pre_data_['non_cred_edge_count'] = noncred_edge_count
        cred_adj_nni_count = pp_maps.get_credible_adj_nni_count(
            nni_engine.adjacent_nnis())
        self.pre_data_['cred_adj_nni_count'] = cred_adj_nni_count
        if self.args.pcsp:
            self.pre_data_['acc_nni_count'] = 1
            self.pre_data_['score'] = -np.inf
        else:
            self.pre_data_['acc_nni_count'] = len(nni_engine.accepted_nnis())
            self.pre_data_['score'] = best_nni_score
        pcsp_pp = pp_maps.get_pcsp_pp(nni)
        self.pre_data_['pcsp_pp'] = pcsp_pp
        pcsp_pp_rank = pp_maps.get_pcsp_pp_rank(
            nni, nni_engine.adjacent_nnis())
        self.pre_data_['pcsp_pp_rank'] = pcsp_pp_rank
        is_nni_new = False
        self.pre_data_['is_nni_new'] = is_nni_new
        pass

    def add_entry(self, iter_count, dag, nni_engine, pp_maps):
        for (nni_id, nni) in enumerate(nni_engine.accepted_nnis()):
            new_row = {
              'iter': iter_count,
              'acc_nni_id': nni_id,
              'acc_nni_count': self.pre_data_['acc_nni_count'],
              'score': nni_engine.scored_nnis()[nni],
              'tree_pp': self.pre_data_['tree_pp'],
              'pcsp_pp': self.pre_data_['pcsp_pp'],
              'is_nni_cred': (self.pre_data_['pcsp_pp'] > 0.0),
              'is_nni_new': ("T" if self.pre_data_['is_nni_new'] else "F"),
              'pcsp_pp_rank': self.pre_data_['pcsp_pp_rank'],
              'node_count': dag.node_count(),
              'edge_count': dag.edge_count(),
              'cred_edge_count': self.pre_data_['cred_edge_count'],
              'tree_count': dag.topology_count(),
              'adj_nni_count': nni_engine.adjacent_nni_count(),
              'new_adj_nni_count': self.pre_data_['new_adj_nni_count'],
              'cred_adj_nni_count': self.pre_data_['cred_adj_nni_count'],
              'llhs_computed': self.pre_data_['llhs_computed'],
              'parent': nni.get_parent().subsplit_to_string(),
              'child': nni.get_child().subsplit_to_string(),
              'pcsp': nni.get_central_edge_pcsp().pcsp_to_string(),
              'nni_hash': f'{Utils.to_hash(nni)}'
            }
            self.data_.loc[len(self.data_.index)] = new_row
        pass

    def found_all_credible_edges(self):
        found_all = self.data_['cred_edge_count'][len(
            self.data_['cred_edge_count']) - 2] == self.pre_data_['cred_edge_total']
        return found_all

    def write_dataframe_to_file(self, file_path):
        df = self.data_
        df.to_csv(file_path)
        return df

    def to_summary(self, timer):
        df = self.data_
        prv_iter = -1
        total_time = 0
        for loop_iter in range(len(df)):
            iter = df['iter'][loop_iter]
            tree_pp = df['tree_pp'][loop_iter]
            score = df['score'][loop_iter]
            hash = df['nni_hash'][loop_iter]
            acc_count = df['acc_nni_count'][loop_iter]
            adj_nni_count = df['adj_nni_count'][loop_iter]
            pcsp_pp = df['pcsp_pp'][loop_iter]
            pcsp_rank = df['pcsp_pp_rank'][loop_iter]
            nni_new = df['is_nni_new'][loop_iter]
            node_count = df['node_count'][loop_iter]
            edge_count = df['edge_count'][loop_iter]
            cred_edge_count = df['cred_edge_count'][loop_iter]
            cred_edge_total = self.pre_data_['cred_edge_total']
            cred_adj_nni_count = df['cred_adj_nni_count'][loop_iter]
            curr_time = timer.get_lap("full_iter", iter - 1)
            if (pcsp_rank > cred_adj_nni_count):
                pcsp_rank = "X"
            if (prv_iter != iter):
                total_time += curr_time
            tree_pp = tree_pp / self.pre_data_['tree_pp_total']
            print_v(
                f'{loop_iter+1:<3} {iter:<3} | {acc_count:<3} {score:.{digits}f} {hash:<16} {nni_new:<2} | {tree_pp:.{digits}f} {pcsp_pp:.{digits}f} {pcsp_rank}/{cred_adj_nni_count} | {node_count:<4} {edge_count:<4} {adj_nni_count:<4} {cred_edge_count}/{cred_edge_total} | {curr_time:.4f}s {int(total_time/60)}m{(total_time%60):.4f}s')
            prv_iter = iter
        pass

    def write_nni_list(self, file_path):
        fp = open(file_path, 'w')
        fp.write(f"#GIT_COMMIT: {bito.git_commit()}\n")
        for i in range(len(self.data_['nni_hash'])):
            nni_hash = self.data_['nni_hash'][i]
            nni_parent = self.data_['parent'][i]
            nni_child = self.data_['child'][i]
            nni_pcsp = self.data_['pcsp'][i]
            fp.write(f"{nni_hash} {nni_pcsp}\n")
        fp.close()


class PosteriorProbabilityMaps:
    def __init__(self, args):
        self.args = args
        self.tree_inst = None
        self.trees = None
        self.pps = None
        self.tree_id_map = None
        self.tree_pp_map = None
        self.pcsp_pp_map = None
        self.init()
        pass

    def init(self):
        self.tree_inst, self.trees = Loader.load_trees(
            self.args.fasta, self.args.credible_newick)
        self.pps = Loader.load_pps(self.args.pp_csv)
        self.tree_id_map = Loader.build_tree_id_map(self.trees)
        self.tree_pp_map = Loader.build_tree_pp_map(
            self.tree_id_map, self.pps)
        self.pcsp_pp_map = Loader.load_pcsp_pp_map(self.args.pcsp_pp_csv)
        pass

    def get_tree_pp(self, dag):
        dag_pp = 0.0
        for tree_id in self.tree_id_map:
            tree = self.tree_id_map[tree_id]
            if (dag.contains_tree(tree)):
                pp = self.tree_pp_map[tree_id]
                dag_pp += pp
        return dag_pp

    def get_tree_pp_total(self):
        tree_pp_total = 0
        for tree_id in self.tree_pp_map:
            tree_pp_total += self.tree_pp_map[tree_id]
        return tree_pp_total

    def get_pcsp_pp(self, nni):
        if type(nni) == bito.nni_op:
            pcsp = nni.get_central_edge_pcsp()
        else:
            pcsp = nni
        if pcsp in self.pcsp_pp_map.keys():
            return self.pcsp_pp_map[pcsp]
        else:
            return 0.0

    def get_pcsp_pp_rank(self, best_nni, adj_nnis):
        best_pcsp_pp = self.get_pcsp_pp(best_nni)
        pcsp_pp_rank = 1
        for nni in adj_nnis:
            nni_pcsp = self.get_pcsp_pp(nni)
            if nni_pcsp > best_pcsp_pp:
                pcsp_pp_rank += 1
        return pcsp_pp_rank

    def get_credible_edge_count(self, dag):
        cred_edge_count = 0
        noncred_edge_count = 0
        pcsps = dag.build_set_of_edge_bitsets()
        for pcsp in pcsps:
            if pcsp in self.pcsp_pp_map:
                cred_edge_count += 1
            else:
                noncred_edge_count += 1
        return (cred_edge_count, noncred_edge_count)

    def get_credible_edge_total(self):
        return len(self.pcsp_pp_map)

    def get_credible_adj_nni_count(self, adj_nnis):
        cred_adj_nni_count = 0
        for nni in adj_nnis:
            pcsp_pp = self.get_pcsp_pp(nni)
            if pcsp_pp > 0.0:
                cred_adj_nni_count += 1
        return cred_adj_nni_count


class NNISearchInstance:
    def __init__(self, args):
        self.args = args
        self.dag_inst = None
        self.dag = None
        self.nni_engine = None
        self.pp_maps = None
        self.init()
        pass

    def args_info(self):
        pass

    def init(self):
        self.dag_inst, self.dag = Loader.load_dag(
            self.args.fasta, self.args.seed_newick)
        pass

    def set_pp_maps(self, pp_maps):
        self.pp_maps = pp_maps
        pass

    def init_engine_for_search(self):
        if self.args.tp:
            self.init_engine_for_tp_search()
        if self.args.gp:
            self.init_engine_for_gp_search()
        if self.args.pcsp:
            self.init_engine_for_pcsp_search()
        self.nni_engine = self.dag_inst.get_nni_engine()
        self.nni_engine.run_init()
        pass

    def init_engine_for_eval_search(self):
        pass

    def init_engine_for_gp_search(self):
        self.dag_inst.make_gp_engine()
        self.dag_inst.make_nni_engine()
        self.dag_inst.take_first_branch_length()
        nni_engine = self.dag_inst.get_nni_engine()
        nni_engine.set_include_rootsplits(args.include_rootsplits)
        nni_engine.set_gp_likelihood_cutoff_filtering_scheme(0.0)
        if (do_rescore_all_nnis != None):
            nni_engine.set_rescore_rejected_nnis(do_rescore_all_nnis)
        if (do_reeval_all_nnis != None):
            nni_engine.set_reevaluate_rejected_nnis(do_reeval_all_nnis)
        nni_engine.set_top_k_score_filtering_scheme(1)
        if self.args.use_cutoff:
            nni_engine.set_gp_likelihood_cutoff_filtering_scheme(
                self.args.threshold)
        if self.args.use_dropoff:
            nni_engine.set_gp_likelihood_drop_filtering_scheme(
                self.args.threshold)
        if self.args.use_top_k:
            nni_engine.set_top_k_score_filtering_scheme(self.args.top_k)
        self.dag_inst.estimate_branch_lengths(1e-5, 3, True)
        pass

    def init_engine_for_tp_search(self):
        self.dag_inst.make_tp_engine()
        self.dag_inst.make_nni_engine()
        self.dag_inst.tp_engine_set_branch_lengths_by_taking_first()
        self.dag_inst.tp_engine_set_choice_map_by_taking_first()
        nni_engine = self.dag_inst.get_nni_engine()
        nni_engine.set_include_rootsplits(self.args.include_rootsplits)
        nni_engine.set_tp_likelihood_cutoff_filtering_scheme(0.0)
        if (do_rescore_all_nnis != None):
            nni_engine.set_rescore_rejected_nnis(do_rescore_all_nnis)
        if (do_reeval_all_nnis != None):
            nni_engine.set_reevaluate_rejected_nnis(do_reeval_all_nnis)
        nni_engine.set_top_k_score_filtering_scheme(1)
        tp_engine = self.dag_inst.get_tp_engine()
        tp_engine.set_optimization_max_iteration(self.args.opt_max)
        tp_engine.set_optimize_new_edges(do_optimize_new_edges)
        tp_engine.set_use_best_edge_map(do_use_best_edge_map)
        tp_engine.set_init_proposed_branch_lengths_with_dag(do_init_proposed_bls_with_dag)
        tp_engine.set_fix_proposed_branch_lengths_from_dag(do_fix_proposed_bls_from_dag)
        if self.args.use_cutoff:
            nni_engine.set_tp_likelihood_cutoff_filtering_scheme(
                args.threshold)
        if self.args.use_dropoff:
            nni_engine.set_tp_likelihood_drop_filtering_scheme(
                self.args.threshold)
        if self.args.use_top_k:
            nni_engine.set_top_k_score_filtering_scheme(self.args.top_k)
        pass

    def init_engine_for_pcsp_search(self):
        self.dag_inst.make_nni_engine()
        nni_engine = self.dag_inst.get_nni_engine()

        def get_pscp_score(nni_engine, nni):
            return self.pp_maps.get_pcsp_pp(nni)
        pcsp_func = get_pscp_score

        nni_engine.set_filter_score_loop_function(pcsp_func)
        nni_engine.set_top_k_score_filtering_scheme(1)
        pass

    def get_best_nni_score(self):
        best_score = -np.inf
        scored_nnis = self.nni_engine.scored_nnis()
        for nni in scored_nnis:
            if best_score < scored_nnis[nni]:
                best_score = scored_nnis[nni]
        return best_score

    def print_scored_nnis(self):
        scored_nnis = self.nni_engine.scored_nnis()
        entries_per_line = 4
        print_v("# scored_nnis:", len(scored_nnis))
        sorted_scored_nnis = [(k, scored_nnis[k]) for k in sorted(
            scored_nnis, key=scored_nnis.get, reverse=True)]
        for i, (nni, score) in enumerate(sorted_scored_nnis):
            sys.stdout.write(
                f'    [{Utils.to_hash(nni)}:  {score:.5f}]')
            if (i + 1) % entries_per_line == 0 or i == len(sorted_scored_nnis) - 1:
                sys.stdout.write('\n')
        pass

    def print_accepted_nnis(self):
        scored_nnis = self.nni_engine.scored_nnis()
        accepted_nnis = self.nni_engine.accepted_nnis()
        print_v("# accepted_nnis:", len(accepted_nnis))
        for nni in accepted_nnis:
            print_v(
                f"    [{Utils.nni_to_hash(nni, abbr)}:  {scored_nnis[nni]:.5f}]")
        pass


class Tracker:
    def __init__(self, args, dag_inst, run_tracker=do_run_tracker):
        self.args = args
        self.dag_inst = dag_inst
        self.old_choice_map = {}
        self.old_proposed_scores = {}
        self.old_dag_scores = {}
        self.old_nni_map = {}
        self.old_pcsp_map = {}
        self.old_pv_hash_map = {}
        self.old_pv_value_map = {}
        self.old_bl_map = {}
        # tracks what iter pcsp/subsplit was added to the DAG.
        self.timedag_pcsps = {}
        self.timedag_subsplits = {}
        self.run_tracker = run_tracker
        self.tolerance = 1e-3

        # track mismatches
        self.choice_map_mismatches = 0
        self.proposed_score_mismatches = 0
        self.dag_score_mismatches = 0
        self.compare_score_mismatches = 0
        self.pv_hash_mismatches = 0
        self.pv_value_mismatches = 0
        self.bl_mismatches = 0
        return

    def is_pcsp_in_watchlist(self, pcsp):
        if (not self.run_tracker):
            return
        parent_hash = Utils.subsplit_to_hash(
            pcsp.pcsp_get_parent_subsplit(), abbr)
        child_hash = Utils.subsplit_to_hash(
            pcsp.pcsp_get_child_subsplit(), abbr)
        in_watchlist = (parent_hash in subsplit_watchlist) or (
            child_hash in subsplit_watchlist)
        return in_watchlist, parent_hash, child_hash

    def update_timedag_pcsps(self, iter_count):
        if (not self.run_tracker):
            return
        for pcsp in self.dag_inst.get_dag().build_set_of_edge_bitsets():
            if pcsp not in self.timedag_pcsps:
                self.timedag_pcsps[pcsp] = iter_count
        return

    def update_timedag_subsplits(self, iter_count):
        if (not self.run_tracker):
            return
        for pcsp in self.dag_inst.get_dag().build_set_of_node_bitsets():
            if pcsp not in self.timedag_subsplits:
                self.timedag_subsplits[pcsp] = iter_count
        return

    def check_for_choice_map_changes(self):
        if (not self.run_tracker):
            return
        tp_engine = self.dag_inst.get_tp_engine()
        new_choice_map = tp_engine.build_map_from_pcsp_to_edge_choice_pcsps()
        for pcsp in new_choice_map:
            if pcsp in self.old_choice_map:
                new_result = new_choice_map[pcsp]
                old_result = self.old_choice_map[pcsp]
                matches = (old_result == new_result)
                result = "MATCH" if matches else "MISMATCH"
                if not matches:
                    pcsp = Utils.pcsp_to_hash(pcsp, abbr)
                    print(f"#CHOICE_MAP_{result}: {pcsp} {new_result} {old_result}")
            else:
                # pcsp = Utils.pcsp_to_hash(pcsp, abbr)
                # print(f"#CHOICE_MAP: PCSP not found! {pcsp}")
                pass
        for pcsp in new_choice_map:
            self.old_choice_map[pcsp] = new_choice_map[pcsp]
        return self.old_choice_map

    def check_for_proposed_score_changes(self):
        if (not self.run_tracker):
            return
        nni_engine = self.dag_inst.get_nni_engine()
        new_proposed_scores = nni_engine.scored_nnis()
        new_nnis = nni_engine.new_adjacent_nnis()
        rescored_nnis = nni_engine.nnis_to_rescore()
        match_cnt = 0
        mismatch_cnt = 0
        new_cnt = 0
        for nni in rescored_nnis:
            pcsp = nni.get_central_edge_pcsp()
            if nni in new_proposed_scores:
                if pcsp in self.old_proposed_scores:
                    old_score = self.old_proposed_scores[pcsp]
                    new_score = new_proposed_scores[nni]
                    score_diff = old_score - new_score
                    score_improved = old_score <= new_score
                    matches = (abs(score_diff) < self.tolerance)
                    result = "MATCH" if matches else "MISMATCH"
                    if not matches:
                        mismatch_cnt += 1
                        pcsp = Utils.pcsp_to_hash(pcsp, abbr)
                        print(
                            f"  #PROP_SCORE_{result}: {pcsp} new::{round(new_score, digits)} old::{round(old_score, digits)} change::{round(score_diff, digits)} is_improved::{score_improved}")
                    else:
                        match_cnt += 1
                else:
                    new_cnt += 1
        self.proposed_score_mismatches += mismatch_cnt
        for nni in new_proposed_scores:
            pcsp = nni.get_central_edge_pcsp()
            self.old_proposed_scores[pcsp] = new_proposed_scores[nni]
        return self.old_proposed_scores

    def check_for_nni_map_changes(self):
        if (not self.run_tracker):
            return
        nni_engine = self.dag_inst.get_nni_engine()
        tp_engine = self.dag_inst.get_tp_engine()
        new_nni_map = tp_engine.build_map_of_proposed_nnis_to_best_pre_nnis(
            nni_engine.new_adjacent_nnis())
        for post_nni in new_nni_map:
            if post_nni in self.old_nni_map:
                new_pre_nni = new_nni_map[post_nni]
                old_pre_nni = self.old_nni_map[post_nni]
                matches = (old_pre_nni == new_pre_nni)
                result = "MATCH" if matches else "MISMATCH"
                if not matches:
                    post_nni = Utils.nni_to_hash(post_nni, abbr)
                    new_pre_nni = Utils.nni_to_hash(new_pre_nni, abbr)
                    old_pre_nni = Utils.nni_to_hash(old_pre_nni, abbr)
                    print(
                        f"  #NNI_MAP_{result}: {post_nni} => {new_pre_nni} {old_pre_nni}")
        for post_nni in new_nni_map:
            self.old_nni_map[post_nni] = new_nni_map[post_nni]
        return self.old_nni_map

    def check_for_pcsp_map_changes(self):
        if (not self.run_tracker):
            return
        nni_engine = self.dag_inst.get_nni_engine()
        tp_engine = self.dag_inst.get_tp_engine()
        new_pcsp_map = tp_engine.build_map_from_pcsp_to_edge_choice_pcsps()
        for pcsp in new_pcsp_map:
            if pcsp in self.old_pcsp_map:
                new_pcsps = new_pcsp_map[pcsp]
                old_pcsps = self.old_pcsp_map[pcsp]
                matches = (old_pcsps == new_pcsps)
                result = "MATCH" if matches else "MISMATCH"
                in_watchlist, parent_hash, child_hash = self.is_pcsp_in_watchlist(pcsp)
                if (not matches):
                    pcsp_hash = Utils.pcsp_to_hash(pcsp, abbr)
                    new_pcsps_hash = [Utils.pcsp_to_hash(x, abbr) for x in new_pcsps]
                    old_pcsps_hash = [Utils.pcsp_to_hash(x, abbr) for x in old_pcsps]
                    print(
                        f"  #PCSP_MAP_{result}: {pcsp_hash} \n    {new_pcsps_hash} \n    {old_pcsps_hash}")
        for focal_pcsp in new_pcsp_map:
            self.old_pcsp_map[focal_pcsp] = new_pcsp_map[focal_pcsp]
        return self.old_pcsp_map

    def check_for_pv_hash_map_changes(self, pcsps=None):
        if (not self.run_tracker):
            return
        nni_engine = self.dag_inst.get_nni_engine()
        tp_engine = self.dag_inst.get_tp_engine()
        new_pv_map = tp_engine.build_map_from_pcsp_to_pv_hashes()
        if pcsps == None:
            pcsps = new_pv_map.keys()
        for pcsp in pcsps:
            if pcsp in self.old_pv_hash_map:
                new_pv = new_pv_map[pcsp]
                old_pv = self.old_pv_hash_map[pcsp]
                matches = (old_pv == new_pv)
                result = "MATCH" if matches else "MISMATCH"
                pcsp_hash = Utils.pcsp_to_hash(pcsp, abbr)
                in_watchlist, parent_hash, child_hash = self.is_pcsp_in_watchlist(pcsp)
                if (not matches):
                    parent_is_root = pcsp.pcsp_get_parent_subsplit().subsplit_is_rootsplit()
                    parent_is_root = "IS_ROOT" if parent_is_root else "IS_NOT_ROOT"
                    child_is_leaf = pcsp.pcsp_get_child_subsplit().subsplit_is_leaf()
                    child_is_leaf = "IS_LEAF" if child_is_leaf else "IS_NOT_LEAF"
                    print(
                        f"  #PV_HASH_MAP_{result}: {pcsp_hash} ({parent_hash} {parent_is_root} {child_hash} {child_is_leaf}) => {pcsp.pcsp_to_string()}")
                    print(f"    NEW::{new_pv}")
                    print(f"    OLD::{old_pv}")
        for pcsp in pcsps:
            self.old_pv_hash_map[pcsp] = new_pv_map[pcsp]
        return self.old_pv_hash_map

    def check_for_pv_value_map_changes(self, pcsps=None):
        if (not self.run_tracker):
            return
        tp_engine = self.dag_inst.get_tp_engine()
        new_pv_map = tp_engine.build_map_from_pcsp_to_pv_values()
        if pcsps == None:
            pcsps = new_pv_map.keys()
        for pcsp in pcsps:
            if pcsp in self.old_pv_value_map:
                new_pv = new_pv_map[pcsp]
                old_pv = self.old_pv_value_map[pcsp]
                new_pv = [np.array(x) for x in new_pv]
                old_pv = [np.array(x) for x in old_pv]
                max_diffs = []
                for i in range(len(new_pv)):
                    max_diff = np.max(np.abs(np.subtract(new_pv[i], old_pv[i])))
                    max_diffs.append(round(max_diff, 5))
                matches = (max(max_diffs) < self.tolerance)
                result = "MATCH" if matches else "MISMATCH"
                if not matches:
                    pcsp = Utils.pcsp_to_hash(pcsp, abbr)
                    print(f"  #PV_VALUE_MAP_{result}: {pcsp} => {max_diffs}")
                    # for i in range(len(new_pv)):
                    #     print(f"### PV_[{i}]")
                    #     print("NEW:", new_pv[i])
                    #     print("OLD:", old_pv[i])
        for pcsp in pcsps:
            self.old_pv_value_map[pcsp] = new_pv_map[pcsp]
        return self.old_pv_value_map

    def check_for_branch_length_map_changes(self):
        if (not self.run_tracker):
            return
        tp_engine = self.dag_inst.get_tp_engine()
        new_bl_map = tp_engine.build_map_from_pcsp_to_branch_length()
        for pcsp in new_bl_map:
            if pcsp in self.old_bl_map:
                new_bl = new_bl_map[pcsp]
                old_bl = self.old_bl_map[pcsp]
                matches = (abs(new_bl - old_bl) < self.tolerance)
                result = "MATCH" if matches else "MISMATCH"
                in_watchlist, parent_hash, child_hash = self.is_pcsp_in_watchlist(pcsp)
                if (not matches):
                    pcsp = Utils.pcsp_to_hash(pcsp, abbr)
                    print(
                        f"  #BL_MAP_{result}: {pcsp} => {round(new_bl, digits)} {round(old_bl, digits)}")
        for pcsp in new_bl_map:
            self.old_bl_map[pcsp] = new_bl_map[pcsp]
        return self.old_bl_map

    def check_for_dag_score_changes(self):
        if (not self.run_tracker):
            return
        tp_engine = self.dag_inst.get_tp_engine()
        new_dag_scores = tp_engine.build_map_from_pcsp_to_score(True)
        # compare old dag_scores to new dag_scores
        match_cnt = 0
        mismatch_cnt = 0
        for pcsp in new_dag_scores:
            if pcsp in self.old_dag_scores:
                new_score = new_dag_scores[pcsp]
                old_score = self.old_dag_scores[pcsp]
                score_diff = abs(new_score - old_score)
                matches = (score_diff < self.tolerance)
                result = "MATCH" if matches else "MISMATCH"
                if (not matches):
                    mismatch_cnt += 1
                    pcsp_hash = Utils.pcsp_to_hash(pcsp, abbr)
                    print(
                        f"  #DAG_SCORE_{result}: {pcsp_hash} => {round(new_score, digits)} {round(old_score, digits)} | {round(score_diff, digits)}")
                    # print(
                    #     f"    {pcsp_hash} {pcsp.pcsp_get_parent_subsplit().subsplit_to_string()} {pcsp.pcsp_get_child_subsplit().subsplit_to_string()}")
                else:
                    match_cnt += 1
        self.dag_score_mismatches += mismatch_cnt
        # update old dag_scores
        for pcsp in new_dag_scores:
            self.old_dag_scores[pcsp] = new_dag_scores[pcsp]
        # compare proposed_scores and dag_scores
        match_cnt = 0
        mismatch_cnt = 0
        for pcsp in self.old_dag_scores:
            if pcsp in self.old_proposed_scores:
                dag_score = self.old_dag_scores[pcsp]
                proposed_score = self.old_proposed_scores[pcsp]
                score_diff = abs(dag_score - proposed_score)
                matches = (score_diff < self.tolerance)
                result = "MATCH" if matches else "MISMATCH"
                if (not matches):
                    mismatch_cnt += 1
                    pcsp_hash = Utils.pcsp_to_hash(pcsp, abbr)
                    print(
                        f"  #DAG_VS_PROP_SCORE_{result}: {pcsp_hash} => {round(dag_score, digits)} {round(proposed_score, digits)} | {round(score_diff, digits)}")
                else:
                    match_cnt += 1
        self.compare_score_mismatches += mismatch_cnt
        # remove proposed_score if it has been added to the dag_scores
        for pcsp in self.old_dag_scores:
            if pcsp in self.old_proposed_scores:
                del self.old_proposed_scores[pcsp]
        # print all dag_scores
        do_print = False
        if do_print:
            for pcsp in self.old_dag_scores:
                pcsp_hash = Utils.pcsp_to_hash(pcsp, abbr)
                dag_score = round(self.old_dag_scores[pcsp], digits)
                print(f"  {pcsp_hash} {dag_score}")
        pass


class Program:
    @staticmethod
    def parse_args(args):
        parser = argparse.ArgumentParser(
            description='Tools for performing NNI systematic search.')
        parser.add_argument('-v', '--verbose',
                            help='verbose', type=int, default=verbose)
        parser.add_argument('-p', '--profiler',
                            action='store_true', help='profile program')
        parser.add_argument('--output-basename',
                            help='set basename for all output files', type=str)

        # ### main programs
        subparsers = parser.add_subparsers(title='programs', dest='program')
        # ### nni search
        subparser1 = subparsers.add_parser(
            'nni-search', help='Perform systematic NNI search.')
        # positional arguments
        subparser1.add_argument('fasta', help='fasta file', type=str)
        subparser1.add_argument(
            'seed_newick', help='newick file for initial trees in DAG', type=str)
        subparser1.add_argument(
            'credible_newick', help='newick file for trees in credible posterior', type=str)
        subparser1.add_argument(
            'pp_csv', help='csv file containing the posterior weights of the trees corresponding to credible newick.', type=str)
        subparser1.add_argument(
            'pcsp_pp_csv', help='csv file containing the per-PCSP posterior weights', type=str)
        # search method group
        group = subparser1.add_mutually_exclusive_group(
            required=True)
        group.add_argument('--gp', action='store_true',
                           help='Selects best NNI according to Generalized Pruning.')
        group.add_argument('--tp', action='store_true',
                           help='Selects best NNI according to Top Pruning via Likelihood.')
        group.add_argument('--pcsp', action='store_true',
                           help='Selects best NNI according to per-PCSP cumulative posterior.')
        # search scheme group
        group = subparser1.add_mutually_exclusive_group(
            required=False)
        group.add_argument('--use-top-k', action='store_true',
                           help='Selects the top N scoring NNIs.')
        group.add_argument('--use-cutoff', action='store_true',
                           help='Selects all NNIs scoring above a given threshold.')
        group.add_argument('--use-dropoff', action='store_true',
                           help='Selects all NNIs scoring above a given dropoff below best NNI.')
        # search scheme arguments
        subparser1.add_argument(
            '--top-n', help='number of NNIs to accept per iteration', type=int, default=1)
        subparser1.add_argument(
            '--threshold', help='cutoff/dropoff threshold', type=float, default=0.0)
        subparser1.add_argument('--opt-init', action='store_true',
                                help='optimize branch lengths when initializing DAG')
        subparser1.add_argument('--opt-max', help='Maximum number of iterations of branch length optimization.',
                                type=int, default=optimization_max_iteration)
        # options
        subparser1.add_argument('-o', '--output', help='output csv file', type=str,
                                default='_out/results.nni_search.csv')
        subparser1.add_argument('--no-output', action='store_true',
                                help='Do not save csv results file.')
        subparser1.add_argument(
            '--iter-max', help='Number of NNI search iterations', type=int, default=10)
        subparser1.add_argument(
            '--test', help='Compare NNI results to a golden run.', type=str)
        subparser1.add_argument(
            '--save-test', help='Save NNI results as golden run.', type=str)
        subparser1.add_argument("--include-rootsplits", action='store_true',
                                help='Whether to include rootsplits in NNI search')

        # ### pcsp map builder
        subparser2 = subparsers.add_parser(
            'build-pcsp-map', help='Build per-PCSP map.')
        # positional arguments
        subparser2.add_argument('fasta', help='fasta file', type=str)
        subparser2.add_argument(
            'credible_newick', help='newick file for trees in credible posterior', type=str)
        subparser2.add_argument(
            'pp_csv', help='csv file containing the posterior weights of the trees from credible_trees', type=str)
        # options
        subparser2.add_argument('-o', '--output', help='output file', type=str,
                                default='results.pcsp_pp_map.csv')

        # ### build credible set
        subparser3 = subparsers.add_parser(
            'build-credible', help='Build credible tree and pp files from trprobs.')
        subparser3.add_argument('trprobs', help='trprobs file', type=str)
        # options
        subparser3.add_argument(
            '-c', '--credible', help='credible cutoff', type=float, default=0.95)
        subparser3.add_argument(
            '-A', '--accept-all', action='store_true', help='accept all trees from input')
        subparser3.add_argument('-p', '--pp-output', help='output pp file', type=str,
                                default='results.cred-pps.csv')
        subparser3.add_argument('-t', '--tree-output', help='output tree newick file', type=str,
                                default='results.cred-trees.nwk')

        # ### fix results
        subparser4 = subparsers.add_parser(
            'fix-results', help='Fix/update results using selected NNI from results.')
        # positional arguments
        subparser4.add_argument('fasta', help='fasta file', type=str)
        subparser4.add_argument(
            'seed_newick', help='newick file for initial trees in DAG', type=str)
        subparser4.add_argument(
            'credible_newick', help='newick file for trees in credible posterior', type=str)
        subparser4.add_argument(
            'pp_csv', help='csv file containing the posterior weights of the trees from credible_trees', type=str)
        subparser4.add_argument(
            'pcsp_pp_csv', help='csv file containing the per-PCSP posterior weights', type=str)
        subparser4.add_argument(
            'nni_search_csv', help='csv file from the per-iteration results of NNI search', type=str)
        # options
        subparser3.add_argument(
            '--update-tree-pp', action='store_true', help='Update tree pp using new pp csv.')

        # run parser
        parsed_args = parser.parse_args(args)
        args_dict = vars(parsed_args)
        return parsed_args

    @staticmethod
    def run_program(args):
        verbose = args.verbose
        if (args.program == 'nni-search'):
            return Program.nni_search(args)
        if (args.program == 'build-pcsp-map'):
            return Program.build_and_save_pcsp_pp_map(args)
        if (args.program == 'build-credible'):
            return Program.build_credible_set(args)
        if (args.program == 'fix-results'):
            return Program.fix_results(args)

    @staticmethod
    def nni_search(args):
        timer = Timer()
        results = Results(args)

        if do_print_settings:
            print_v("# NNI SEARCH SETTINGS:")
            print_v(args)

        if do_print_setup_data:
            print_v("# load nni_search_instance...")
        nni_inst = NNISearchInstance(args)
        if do_print_setup_data:
            print_v("# load pp_maps...")
        pp_maps = PosteriorProbabilityMaps(args)
        if do_print_setup_data:
            print_v("# all data loaded.")
        nni_inst.set_pp_maps(pp_maps)

        if do_print_setup_data:
            print_v("# initialize nni engine...")
        nni_inst.init_engine_for_search()

        dag = nni_inst.dag
        nni_engine = nni_inst.nni_engine
        if do_print_dag_stats:
            print_v(f"# init_dag: {dag.node_count()} {dag.edge_count()}")

        results.data_init()
        results.predata_init(dag, nni_engine, pp_maps)

        tracker = Tracker(args, nni_inst.dag_inst)
        subroutine_names = ["graft_adjacent_nnis_to_dag", "filter_pre_score", "filter_score_adjacent_nnis", "filter_post_score", "filter_evaluate_adjacent_nnis", "remove_all_graft_nnis_from_dag", "add_accepted_nnis_to_dag"]

        ### NNI SEARCH MAIN LOOP ###
        iter_count = 1
        while iter_count <= args.iter_max:
            if do_print_iter_data:
                print_v("--- + ---")
                print_v(f"# iter_count: {iter_count} of {args.iter_max}...")

            # capture pcsps before grafting for PV comparison.
            if do_check_for_pv_map_changes:
                pcsps = dag.build_set_of_edge_bitsets()

            # capture initial data.
            results.predata_begin_iter(iter_count, dag, nni_engine, pp_maps)
            tree_pp = results.pre_data_['tree_pp']
            tree_pp_total = results.pre_data_['tree_pp_total']
            if do_print_dag_stats:
                print_v("# dag:", dag.node_count(), dag.edge_count())
                print_v(f"# dag_tree_pp: {round(tree_pp, digits)} {round((tree_pp/tree_pp_total)*100, 2)}%")

            # nni_engine.sync_adjacent_nnis_with_dag()

            if do_print_dag_stats:
                print_v(f"# adjacent_nnis: {len(nni_engine.adjacent_nnis())}")
                print_v(f"# new_adjacent_nnis: {len(nni_engine.new_adjacent_nnis())}")
                if nni_inst.args.tp or nni_inst.args.gp:
                    bls = nni_engine.get_branch_lengths()
                    print_v(f"# branch_lengths: {len(bls)} {bls}")

            # begin inner loop
            timer.start()
            timer.start("full_iter")
            timer.start("inner_iter")

            timer.start()
            nni_engine.graft_adjacent_nnis_to_dag()
            timer.lap_next("graft_adjacent_nnis_to_dag")

            timer.start()
            nni_engine.filter_pre_score()
            timer.lap_next("filter_pre_score")

            timer.start()
            nni_engine.filter_score_adjacent_nnis()
            timer.lap_next("filter_score_adjacent_nnis")

            timer.start()
            nni_engine.filter_post_score()
            timer.lap_next("filter_post_score")

            timer.stop("full_iter")
            timer.stop("inner_iter")
            results.predata_mid_iter(iter_count, dag, nni_engine, pp_maps)

            if do_print_scored_nnis:
                nni_inst.print_scored_nnis()

            # if searching by pcsp and no nni has a pcsp pp over zero, end search.
            if args.pcsp and nni_inst.get_best_nni_score() <= 0.0:
                break
            timer.resume("full_iter")
            timer.resume("inner_iter")

            timer.start()
            nni_engine.filter_evaluate_adjacent_nnis()
            timer.lap_next("filter_evaluate_adjacent_nnis")

            timer.stop("inner_iter")
            if do_print_accepted_nnis:
                nni_inst.print_accepted_nnis()
            timer.resume("inner_iter")

            timer.start()
            nni_engine.remove_all_graft_nnis_from_dag()
            timer.lap_next("remove_all_graft_nnis_from_dag")

            timer.start()
            nni_engine.add_accepted_nnis_to_dag()
            timer.lap_next("add_accepted_nnis_to_dag")

            timer.stop("full_iter")
            timer.lap_name("inner_iter")

            # check that inner subroutine times sum to inner total time
            if True:
                inner_iter_total = 0.0
                for subroutine_name in subroutine_names:
                    inner_iter_total += timer.lap_last(subroutine_name)
                print(f"#COMPARE inner_iter_total: {timer.lap_last('inner_iter')} {inner_iter_total}")

            # check for changes to the PVs.
            if do_check_for_pv_map_changes:
                pcsps = dag.build_set_of_edge_bitsets()
                current_pv_map = tracker.check_for_pv_hash_map_changes(pcsps)
                current_pv_map = tracker.check_for_pv_value_map_changes(pcsps)
                pass

            ### LOG DATA ENTRY ###
            results.predata_end_iter(iter_count, dag, nni_engine, pp_maps)
            results.add_entry(iter_count, dag, nni_engine, pp_maps)
            results.write_dataframe_to_file(args.output)

            # check for NNI score changes.
            if do_check_for_nni_score_changes:
                current_proposed_score_map = tracker.check_for_proposed_score_changes()
                pass

            # check for changes to choice_map.
            if do_check_for_choice_map_changes:
                current_choice_map = tracker.check_for_choice_map_changes()
                pass

            # check for changes to the NNI map from proposed-NNI to best pre-NNI.
            if do_check_for_nni_map_changes:
                current_nni_map = tracker.check_for_nni_map_changes()
                pass

            # check for changes to the map from proposed-NNI PCSPs to best pre-NNI PCSPs.
            if do_check_for_pcsp_map_changes:
                current_pcsp_map = tracker.check_for_pcsp_map_changes()
                pass

            # check for changes to the branch lengths.
            if do_check_for_bl_map_changes:
                current_bl_map = tracker.check_for_branch_length_map_changes()
                pass

            # check for changes to DAG scores.
            if do_check_for_dag_score_changes:
                current_dag_score_map = tracker.check_for_dag_score_changes()
                pass

            # end iterations if no accepted nnis.
            if len(nni_engine.accepted_nnis()) == 0:
                print_v("# NO ACCEPTED NNIS")
                break

            timer.resume("full_iter")
            timer.start()
            nni_engine.run_post_loop()
            timer.lap_next("run_post_loop")
            timer.stop("full_iter")

            # track changes made to dag.
            if do_check_for_dag_changes:
                tracker.update_timedag_pcsps(iter_count)
                tracker.update_timedag_subsplits(iter_count)

            # bookmark by outputting top trees
            iter_count += 1
            timer.resume("full_iter")
            timer.lap_name("full_iter")
            # end search if all credible edges have been found
            if do_end_when_all_creds_found and results.found_all_credible_edges():
                break

        if do_print_iter_data:
            print_v("--- + ---")

        # write final results to file
        df = results.data_
        print(f"OUTER args.no_output: {args.no_output}")
        if not args.no_output:
            print(f"INNER args.no_output: {args.no_output}")
            results.write_dataframe_to_file(args.output)

        ### SUMMARY OUTPUT ###
        if do_print_tracker_summary:
            print_v("\n=== TRACKER SUMMARY ===")
            print_v(f"  proposed_score_mismatches: {tracker.proposed_score_mismatches}")
            print_v(f"  dag_score_mismatches: {tracker.dag_score_mismatches}")
            print_v(f"  compare_score_mismatches: {tracker.compare_score_mismatches}")

        if do_print_summary:
            print_v("\n=== FINAL DATAFRAME ===")
            pd.set_option('display.max_colwidth', None)
            print_v(df)
            print_v("\n=== TIMES_SUMMARY ===")
            timer.to_summary()
            print_v("\n=== RESULTS_SUMMARY ===")
            results.to_summary(timer)

        # compare results to golden run.
        if args.test:
            nni_list = results.data_["nni_hash"]
            golden_git_commit, golden_nni_list = Loader.load_nni_list(args.test)
            print_v("\n=== TEST ===")
            print_v(f"# GIT_COMMIT: TEST::{bito.git_commit()} GOLDEN::{golden_git_commit}")
            if (nni_list == golden_nni_list):
                print_v("# TEST_PASSED -- accepted NNIs matches golden run.")
            else:
                print_v("# TEST_FAILED -- accepted NNIs do not match golden run.")

            print_v(f"# {'ITER':<7} {'RESULT':<10} {'TEST_NNI':<20} {'GOLDEN_NNI':<20}")
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
                print_v(f"  {i:<7} {res:<10} {nni:<20} {golden_nni:<20}")
        # write results
        if args.save_test:
            print_v(f"# NNI list written to: {args.save_test}")
            results.write_nni_list(args.save_test)
        return nni_inst, results, timer

    @staticmethod
    def build_credible_set(args):
        taxon_id_map, tree_nwk_map, tree_pp_map, tree_cpp_map = Loader.load_trprobs(args.trprobs, args)
        print_v("Credible:")
        print_v(tree_cpp_map)
        fp_trees = open(args.tree_output, 'w')
        fp_pps = open(args.pp_output, 'w')
        for tree_id in tree_nwk_map:
            if (tree_cpp_map[tree_id] < args.credible) or args.accept_all:
                fp_trees.write(f"{tree_nwk_map[tree_id]}\n")
                fp_pps.write(f"{tree_pp_map[tree_id]}\n")
        fp_trees.close()
        fp_pps.close()
        return

    @staticmethod
    def build_and_save_pcsp_pp_map(args):
        print_v("# load dag...")
        dag_inst, dag = Loader.load_dag(args.fasta, args.credible_newick)
        print_v("# load trees...")
        tree_inst, trees = Loader.load_trees(args.fasta, args.credible_newick)
        print_v("# load pps...")
        pps = Loader.load_pps(args.pp_csv)
        print_v("# build maps...")
        tree_id_map = Loader.build_tree_id_map(trees)
        tree_pp_map = Loader.build_tree_pp_map(tree_id_map, pps)
        pcsp_pp_map = Loader.build_pcsp_pp_map(dag, tree_id_map, tree_pp_map)
        print_v("pcsp_pp_map:", len(pcsp_pp_map), pcsp_pp_map)
        my_dict = {
            'parent': [],
            'child': [],
            'pcsp_pp': []
        }
        for pcsp in pcsp_pp_map:
            parent = pcsp.pcsp_get_parent_subsplit().subsplit_to_string()
            child = pcsp.pcsp_get_child_subsplit().subsplit_to_string()
            pcsp_pp = get_pcsp_pp(pcsp, pcsp_pp_map)
            my_dict["parent"].append(parent)
            my_dict["child"].append(child)
            my_dict["pcsp_pp"].append(pcsp_pp)
        df = pd.DataFrame(my_dict)
        df.to_csv(args.output)
        return pcsp_pp_map

    @staticmethod
    def fix_results(args):
        print_v("nni_search")
        print_v("# load data...")
        df = pd.DataFrame(args.nni_search_csv)
        print_v("# load trees...")
        tree_inst, trees = load_trees(args.fasta, args.credible_newick)
        print_v("# load pps...")
        pps = load_pps(args.pp_csv)
        print_v("# build maps...")
        tree_id_map = build_tree_id_map(trees)
        tree_pp_map = build_tree_pp_map(tree_id_map, pps)
        pcsp_pp_map = load_pcsp_pp_map(args.pcsp_pp_csv)
        final_dict = final_data_init()

        print_v("# load dag...")
        dag_inst, _ = load_dag(args.fasta, args.seed_newick)
        dag = dag_inst.get_dag()
        print_v("# init engine...")
        dag_inst.make_nni_engine()
        nni_engine.get_nni_engine()

        iter_max = len(df)
        while iter_count < args.iter_max:
            parent = bito.bitset(df.iloc[iter_count]["parent"])
            child = bito.bitset(df.iloc[iter_count]["child"])
            nni = bito.nni_op(parent, child)
            dag.add_node_pair(nni.get_parent(), nni.get_child())

            iter_count += 1
        pass


############
### MAIN ###
############

if __name__ == "__main__":
    print_v(f"# begin... [commit={bito.git_commit()}]")
    args = Program.parse_args(sys.argv[1:])
    Program.run_program(args)
    print_v(f"# ...done [commit={bito.git_commit()}]")
