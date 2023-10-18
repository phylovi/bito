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


verbose = True
sort_taxa = False
# rounding digits
digits = 5
# number of optimization iterations
max_opt_count = 5
do_rescore_rejected_nnis = False
do_reeval_rejected_nnis = True
# terminate search when all credible edges found
do_end_when_all_creds_found = True
# print data
do_print_scored_nnis = True
do_print_accepted_nnis = True
# track if individual nni scores have changed between iterations.
do_track_nni_score_changes = True


def print(*args):
    __builtins__.print(
        *("%.5f" % a if isinstance(a, float) else a for a in args))


def print_v(*args):
    if verbose:
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

    def lap(self, timer_name="__default__", timer_prv_name="__default__", do_print=default_print):
        time_elapsed = self.elapsed(timer_prv_name)
        if (timer_name not in self.times):
            self.times[timer_name] = []
        self.times[timer_name].append(time_elapsed)
        self.times_prv[timer_prv_name] = time.time()
        if (do_print):
            print_v(f'  #TIME {timer_name}: {time_elapsed}')
        return time_elapsed

    def lap_next(self, timer_name="__default__", do_print=default_print):
        return self.lap(timer_name, "__default__", do_print)

    def lap_name(self, timer_name="__default__", do_print=default_print):
        return self.lap(timer_name, timer_name, do_print)

    def lap_last(self, timer_name="__default__"):
        return self.times[timer_name][len(self.times[timer_name])-1]

    def get_lap(self, timer_name="__default__", lap_count=0):
        return self.times[timer_name][lap_count]

    def to_summary(self):
        times_perc = {}
        for key in self.times:
            times_perc[key] = []
            for time, total in zip(self.times[key], self.times["full_iter"]):
                times_perc[key].append(round(time / total, digits))

        print_v("=== TIMES_REAL ===")
        for key in self.times:
            print_v(
                f'{key}:\n    {Utils.average_over_ranges(self.times[key], 10)}')
            x_s = range(len(self.times[key]))
            fit = np.polyfit(y=self.times[key], x=x_s, deg=2)
            print_v(
                f'    #POLYFIT: {round(fit[0], digits)}x^2 + {round(fit[1], digits)}x + {round(fit[2], digits)}')
            print_v(
                f'    #POLYFIT: {round(fit[0]/abs(fit[2]), digits)}x^2 + {round(fit[1]/abs(fit[2]), digits)}x + {round(fit[2]/fit[2], digits)}')

        print_v("\n=== TIMES_PERC ===")
        for key in times_perc:
            print_v(
                f'{key}:\n    {Utils.average_over_ranges(times_perc[key], 10)}')
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
        step = int(len(data)/buckets)
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
    def to_hash_string(obj):
        string = f'0x{obj.__hash__():016X}'
        return string


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
        # !CHANGE
        dag_pcsps = dag.build_set_of_edge_bitsets()
        # dag_pcsps = dag.build_vector_of_edge_bitsets()
        # dag_pcsps = dag.build_sorted_vector_of_edge_bitsets()
        pcsp_pp_map = {}
        for pcsp_count, pcsp in enumerate(dag_pcsps):
            print_v(f"# loading pcsp {pcsp_count} of {len(dag_pcsps)}...")
            pcsp_pp_map[pcsp] = 0.0
        for tree_id in tree_pp_map:
            tree = tree_id_map[tree_id]
            pp = tree_pp_map[tree_id]
            # !CHANGE
            tree_pcsps = tree.build_set_of_pcsps()
            # tree_pcsps = tree.build_vector_of_pcsps()
            # tree_pcsps = tree.build_sorted_vector_of_pcsps()
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
            elif fields[1] == "git_commit:":
                git_commit = fields[2]
        nni_list_fp.close()
        return git_commit, nni_list


class Results:
    def __init__(self, args):
        self.data_ = {}
        self.pre_data_ = {}
        self.df_ = None
        self.args_ = args
        pass

    def data_init(self):
        self.data_ = {
            'iter': [],
            'acc_nni_id': [],
            'acc_nni_count': [],
            'score': [],
            'is_nni_cred': [],
            'is_nni_new': [],
            'tree_pp': [],
            'pcsp_pp': [],
            'pcsp_pp_rank': [],
            'node_count': [],
            'edge_count': [],
            'cred_edge_count': [],
            'tree_count': [],
            'adj_nni_count': [],
            'new_adj_nni_count': [],
            "cred_adj_nni_count": [],
            'llhs_computed': [],
            'nni_hash': [],
            'parent': [],
            'child': [],
            'pcsp': []
        }
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
        # !CHANGE
        new_adj_nni_count = len(
            nni_engine.adjacent_nnis()) - self.pre_data_['prev_adj_nni_count']
        # new_adj_nni_count = len(nni_engine.new_adjacent_nnis())
        self.pre_data_['new_adj_nni_count'] = new_adj_nni_count
        prev_adj_nni_count = len(nni_engine.adjacent_nnis())
        self.pre_data_['prev_adj_nni_count'] = prev_adj_nni_count
        # !CHANGE
        llhs_computed = self.pre_data_['llhs_computed']
        # llhs_computed += len(nni_engine.nnis_to_rescore())
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
        if args.pcsp:
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
        # !CHANGE
        # is_nni_new = nni in nni_engine.new_adjacent_nnis()
        is_nni_new = False
        self.pre_data_['is_nni_new'] = is_nni_new
        pass

    def add_entry(self, iter_count, dag, nni_engine, pp_maps):
        for (nni_id, nni) in enumerate(nni_engine.accepted_nnis()):
            self.data_['iter'].append(iter_count + 1)
            self.data_['acc_nni_id'].append(nni_id)
            self.data_['acc_nni_count'].append(self.pre_data_['acc_nni_count'])
            self.data_['score'].append(nni_engine.scored_nnis()[nni])
            self.data_['tree_pp'].append(self.pre_data_['tree_pp'])
            self.data_['pcsp_pp'].append(self.pre_data_['pcsp_pp'])
            self.data_['is_nni_cred'].append((self.pre_data_['pcsp_pp'] > 0.0))
            self.data_['is_nni_new'].append(
                "T" if self.pre_data_['is_nni_new'] else "F")
            self.data_['pcsp_pp_rank'].append(self.pre_data_['pcsp_pp_rank'])
            self.data_['node_count'].append(dag.node_count())
            self.data_['edge_count'].append(dag.edge_count())
            self.data_['cred_edge_count'].append(
                self.pre_data_['cred_edge_count'])
            self.data_['tree_count'].append(dag.topology_count())
            self.data_['adj_nni_count'].append(nni_engine.adjacent_nni_count())
            self.data_['new_adj_nni_count'].append(
                self.pre_data_['new_adj_nni_count'])
            self.data_['cred_adj_nni_count'].append(
                self.pre_data_['cred_adj_nni_count'])
            self.data_['llhs_computed'].append(self.pre_data_['llhs_computed'])
            self.data_['nni_hash'].append(f'{Utils.to_hash_string(nni)}')
            self.data_['parent'].append(
                nni.get_parent().subsplit_to_string())
            self.data_['child'].append(
                nni.get_child().subsplit_to_string())
            self.data_['pcsp'].append(
                nni.get_central_edge_pcsp().pcsp_to_string())
        pass

    def found_all_credible_edges(self):
        found_all = self.data_['cred_edge_count'][len(
            self.data_['cred_edge_count'])-2] == self.pre_data_['cred_edge_total']
        return found_all

    def to_dataframe(self):
        self.df_ = pd.DataFrame(self.data_)
        return self.df_

    def write_dataframe_to_file(self, file_path):
        df = self.to_dataframe()
        df.to_csv(file_path)
        return df

    def to_summary(self, timer, file_path):
        fp = open(file_path, 'w')
        df = self.to_dataframe()
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
            tree_pp = tree_pp/self.pre_data_['tree_pp_total']
            print_v(
                f'{loop_iter+1:<3} {iter:<3} | {acc_count:<3} {score:.{digits}f} {hash:<16} {nni_new:<2} | {tree_pp:.{digits}f} {pcsp_pp:.{digits}f} {pcsp_rank}/{cred_adj_nni_count} | {node_count:<4} {edge_count:<4} {adj_nni_count:<4} {cred_edge_count}/{cred_edge_total} | {curr_time:.4f}s {int(total_time/60)}m{(total_time%60):.4f}s')
            fp.write(
                f'{tree_pp:.{digits}f},{score:.{digits}f},{hash},{adj_nni_count}\n')
            prv_iter = iter
        fp.close()

    def write_nni_list(self, file_path):
        fp = open(file_path, 'w')
        fp.write(f"# git_commit: {bito.git_commit()}\n")
        for i in range(len(self.data_['nni_hash'])):
            nni_hash = self.data_['nni_hash'][i]
            nni_parent = self.data_['parent'][i]
            nni_child = self.data_['child'][i]
            nni_pcsp = self.data_['pcsp'][i]
            fp.write(f"{nni_hash} {nni_pcsp}\n")
        fp.close()


class PosteriorProbabilityMaps:
    def __init__(self, args):
        self.args_ = args
        self.tree_inst_ = None
        self.trees_ = None
        self.pps_ = None
        self.tree_id_map_ = None
        self.tree_pp_map_ = None
        self.pcsp_pp_map_ = None
        self.init()
        pass

    def init(self):
        self.tree_inst_, self.trees_ = Loader.load_trees(
            self.args_.fasta, self.args_.credible_newick)
        self.pps_ = Loader.load_pps(self.args_.pp_csv)
        self.tree_id_map_ = Loader.build_tree_id_map(self.trees_)
        self.tree_pp_map_ = Loader.build_tree_pp_map(
            self.tree_id_map_, self.pps_)
        self.pcsp_pp_map_ = Loader.load_pcsp_pp_map(self.args_.pcsp_pp_csv)
        pass

    def get_tree_pp(self, dag):
        dag_pp = 0.0
        for tree_id in self.tree_id_map_:
            tree = self.tree_id_map_[tree_id]
            if (dag.contains_tree(tree)):
                pp = self.tree_pp_map_[tree_id]
                dag_pp += pp
        return dag_pp

    def get_tree_pp_total(self):
        tree_pp_total = 0
        for tree_id in self.tree_pp_map_:
            tree_pp_total += self.tree_pp_map_[tree_id]
        return tree_pp_total

    def get_pcsp_pp(self, nni):
        if type(nni) == bito.nni_op:
            pcsp = nni.get_central_edge_pcsp()
        else:
            pcsp = nni
        if pcsp in self.pcsp_pp_map_.keys():
            return self.pcsp_pp_map_[pcsp]
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
        # !CHANGE
        pcsps = dag.build_set_of_edge_bitsets()
        # pcsps = dag.build_vector_of_edge_bitsets()
        # pcsps = dag.build_sorted_vector_of_edge_bitsets()
        for pcsp in pcsps:
            if pcsp in self.pcsp_pp_map_:
                cred_edge_count += 1
            else:
                noncred_edge_count += 1
        return (cred_edge_count, noncred_edge_count)

    def get_credible_edge_total(self):
        return len(self.pcsp_pp_map_)

    def get_credible_adj_nni_count(self, adj_nnis):
        cred_adj_nni_count = 0
        for nni in adj_nnis:
            pcsp_pp = self.get_pcsp_pp(nni)
            if pcsp_pp > 0.0:
                cred_adj_nni_count += 1
        return cred_adj_nni_count


class NNISearchInstance:
    def __init__(self, args):
        self.args_ = args
        self.dag_inst_ = None
        self.dag_ = None
        self.nni_engine_ = None
        self.pp_maps_ = None
        self.init()
        pass

    def init(self):
        self.dag_inst_, self.dag_ = Loader.load_dag(
            self.args_.fasta, self.args_.seed_newick)
        pass

    def set_pp_maps(self, pp_maps):
        self.pp_maps_ = pp_maps
        pass

    def init_engine_for_search(self):
        if self.args_.tp:
            self.init_engine_for_tp_search()
        if self.args_.gp:
            self.init_engine_for_gp_search()
        if self.args_.pcsp:
            self.init_engine_for_pcsp_search()
        self.nni_engine_ = self.dag_inst_.get_nni_engine()
        self.nni_engine_.run_init(False)
        pass

    def init_engine_for_eval_search(self):
        pass

    def init_engine_for_gp_search(self):
        self.dag_inst_.make_gp_engine()
        self.dag_inst_.make_nni_engine()
        self.dag_inst_.take_first_branch_length()
        nni_engine = self.dag_inst_.get_nni_engine()
        nni_engine.set_include_rootsplits(args.include_rootsplits)
        nni_engine.set_gp_likelihood_cutoff_filtering_scheme(0.0)
        # !CHANGE
        # nni_engine.set_top_n_score_filtering_scheme(1)
        nni_engine.set_top_k_score_filtering_scheme(1)
        if self.args_.use_cutoff:
            nni_engine.set_gp_likelihood_cutoff_filtering_scheme(
                self.args_.threshold)
        if self.args_.use_dropoff:
            nni_engine.set_gp_likelihood_drop_filtering_scheme(
                self.args_.threshold)
        if self.args_.use_top_k:
            # !CHANGE
            # nni_engine.set_top_n_score_filtering_scheme(args.top_k)
            nni_engine.set_top_k_score_filtering_scheme(self.args_.top_k)
        self.dag_inst_.estimate_branch_lengths(1e-5, 3, True)
        pass

    def init_engine_for_tp_search(self):
        self.dag_inst_.make_tp_engine()
        self.dag_inst_.make_nni_engine()
        self.dag_inst_.tp_engine_set_branch_lengths_by_taking_first()
        self.dag_inst_.tp_engine_set_choice_map_by_taking_first()
        nni_engine = self.dag_inst_.get_nni_engine()
        nni_engine.set_include_rootsplits(self.args_.include_rootsplits)
        nni_engine.set_tp_likelihood_cutoff_filtering_scheme(0.0)
        # !CHANGE
        nni_engine.set_rescore_rejected_nnis(do_rescore_rejected_nnis)
        nni_engine.set_reevaluate_rejected_nnis(do_reeval_rejected_nnis)
        # !CHANGE
        # nni_engine.set_top_n_score_filtering_scheme(1)
        nni_engine.set_top_k_score_filtering_scheme(1)
        # !CHANGE
        self.dag_inst_.get_tp_engine().set_optimization_max_iteration(max_opt_count)
        if self.args_.use_cutoff:
            nni_engine.set_tp_likelihood_cutoff_filtering_scheme(
                args.threshold)
        if self.args_.use_dropoff:
            nni_engine.set_tp_likelihood_drop_filtering_scheme(
                self.args_.threshold)
        if self.args_.use_top_k:
            # !CHANGE
            # nni_engine.set_top_n_score_filtering_scheme(args.top_k)
            nni_engine.set_top_k_score_filtering_scheme(self.args_.top_k)
        pass

    def init_engine_for_pcsp_search(self):
        self.dag_inst_.make_nni_engine()
        nni_engine = self.dag_inst_.get_nni_engine()

        def get_pscp_score(nni_engine, nni):
            return self.pp_maps_.get_pcsp_pp(nni)
        pcsp_func = get_pscp_score

        nni_engine.set_filter_score_loop_function(pcsp_func)
        nni_engine.set_top_k_score_filtering_scheme(1)
        pass

    def get_best_nni_score(self):
        best_score = -np.inf
        scored_nnis = self.nni_engine_.scored_nnis()
        for nni in scored_nnis:
            if best_score < scored_nnis[nni]:
                best_score = scored_nnis[nni]
        return best_score

    def print_scored_nnis(self):
        scored_nnis = self.nni_engine_.scored_nnis()
        entries_per_line = 4
        print_v("# scored_nnis:", len(scored_nnis))
        sorted_scored_nnis = [(k, scored_nnis[k]) for k in sorted(
            scored_nnis, key=scored_nnis.get, reverse=True)]
        for i, (nni, score) in enumerate(sorted_scored_nnis):
            sys.stdout.write(
                f'    [{Utils.to_hash_string(nni)}:  {score:.5f}]')
            if (i+1) % entries_per_line == 0 or i == len(sorted_scored_nnis)-1:
                sys.stdout.write('\n')
        pass

    def print_accepted_nnis(self):
        scored_nnis = self.nni_engine_.scored_nnis()
        accepted_nnis = self.nni_engine_.accepted_nnis()
        print_v("# accepted_nnis:", len(accepted_nnis))
        for nni in accepted_nnis:
            print_v(
                f"    [{Utils.to_hash_string(nni)}:  {scored_nnis[nni]:.5f}]")
        pass


class Program:
    @staticmethod
    def parse_args(args):
        parser = argparse.ArgumentParser(
            description='Tools for performing NNI systematic search.')
        parser.add_argument('-v', '--verbose',
                            help='verbose', type=int, default=1)
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
        # positional arguments
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
        subparser1.add_argument('--init-optimize', action='store_true',
                                help='optimize branch lengths when initializing DAG')
        # options
        subparser1.add_argument('-o', '--output', help='output file', type=str,
                                default='results.nni_search.csv')
        subparser1.add_argument(
            '--iter-max', help='number of NNI search iterations', type=int, default=10)
        subparser1.add_argument(
            '--test', help='Compare NNI results to a golden run.', type=str)
        subparser1.add_argument(
            '--save-test', help='Save NNI results as golden run.', type=str)
        subparser1.add_argument(
            '--nni-info', help='Give additional NNI info per iteration', type=str)
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
        verbose = (args.verbose > 0)
        if (args.program == 'nni-search'):
            Program.nni_search(args)
        if (args.program == 'build-pcsp-map'):
            Program.build_and_save_pcsp_pp_map(args)
        if (args.program == 'build-credible'):
            Program.build_credible_set(args)
        if (args.program == 'fix-results'):
            Program.fix_results(args)

    @staticmethod
    def nni_search(args):
        timer = Timer()
        results = Results(args)

        print_v("# load nni_search_instance...")
        nni_inst = NNISearchInstance(args)
        print_v("# load pp_maps...")
        pp_maps = PosteriorProbabilityMaps(args)
        print_v("# all data loaded.")
        nni_inst.set_pp_maps(pp_maps)

        print_v("# initialize nni engine...")
        nni_inst.init_engine_for_search()

        dag = nni_inst.dag_
        nni_engine = nni_inst.nni_engine_
        print_v(f"# init_dag: {dag.node_count()} {dag.edge_count()}")

        print_v("# initialize results...")
        results.data_init()
        results.predata_init(dag, nni_engine, pp_maps)

        iter_count = 0
        while iter_count < args.iter_max:
            print_v("--- ---")
            # run iteration of search
            results.predata_begin_iter(
                iter_count, dag, nni_engine, pp_maps)
            print_v(f"# iter_count: {iter_count + 1} of {args.iter_max}...")
            print_v("# dag:", dag.node_count(), dag.edge_count())
            tree_pp = results.pre_data_['tree_pp']
            tree_pp_total = results.pre_data_['tree_pp_total']
            print_v(f"# dag_tree_pp: {tree_pp} {tree_pp/tree_pp_total}")

            timer.start()
            timer.start("full_iter")
            timer.start("inner_iter")

            ### NNI_SEARCH MAIN_LOOP ###
            # nni_engine.sync_adjacent_nnis_with_dag()
            print_v("# adjacent_nnis:", nni_engine.adjacent_nni_count())
            timer.start()

            nni_engine.graft_adjacent_nnis_to_dag()
            timer.lap_next("graft_adjacent_nnis_to_dag")

            nni_engine.filter_pre_score()
            timer.lap_next("filter_pre_score")

            # !CHANGE
            # nni_engine.filter_eval_adjacent_nnis()
            nni_engine.filter_score_adjacent_nnis()
            timer.lap_next("filter_score_adjacent_nnis")

            nni_engine.filter_post_score()
            timer.lap_next("filter_post_score")

            results.predata_mid_iter(
                iter_count, dag, nni_engine, pp_maps)

            if do_print_scored_nnis:
                nni_inst.print_scored_nnis()

            # if searching by pcsp and no nni has a pcsp pp over zero, end search.
            if args.pcsp and nni_inst.get_best_nni_score() <= 0.0:
                break

            timer.start()
            # !CHANGE
            # nni_engine.filter_process_adjacent_nnis()
            nni_engine.filter_evaluate_adjacent_nnis()
            timer.lap_next("filter_evaluate_adjacent_nnis")

            if do_print_accepted_nnis:
                nni_inst.print_accepted_nnis()

            timer.start()
            nni_engine.remove_all_graft_nnis_from_dag()
            timer.lap_next("remove_all_graft_nnis_from_dag")

            nni_engine.add_accepted_nnis_to_dag(False)
            timer.lap_next("add_accepted_nnis_to_dag")

            timer.lap_name("inner_iter")
            timer.stop("full_iter")

            ### LOG DATA ENTRY ###
            results.predata_end_iter(
                iter_count, dag, nni_engine, pp_maps)
            results.add_entry(iter_count, dag, nni_engine, pp_maps)
            results.write_dataframe_to_file(args.output)

            if len(nni_engine.accepted_nnis()) == 0:
                print_v("# NO ACCEPTED NNIS")
                break
            timer.resume("full_iter")
            timer.start()
            nni_engine.run_post_loop()
            timer.lap_next("run_post_loop")

            # bookmark by outputting top trees
            iter_count += 1
            timer.lap_name("full_iter")
            # end search if all credible edges have been found
            if do_end_when_all_creds_found and results.found_all_credible_edges():
                break

        print_v("\n=== FINAL DATAFRAME ===")
        df = results.write_dataframe_to_file(args.output)
        pd.set_option('display.max_colwidth', None)
        print_v(df)

        ### SUMMARY OUTPUT ###
        print_v("\n=== TIMES_SUMMARY ===")
        timer.to_summary()
        print_v("\n=== RESULTS_SUMMARY ===")
        results.to_summary(timer, "result.csv")

        # compare results to golden run.
        if args.test:
            nni_list = results.data_["nni_hash"]
            golden_git_commit, golden_nni_list = Loader.load_nni_list(
                args.test)
            print_v("\n=== TEST ===")
            print_v(
                f"# GIT_COMMIT: TEST::{bito.git_commit()} GOLDEN::{golden_git_commit}")
            if (nni_list == golden_nni_list):
                print_v("# TEST_PASSED -- accepted NNIs matches golden run.")
            else:
                print_v("# TEST_FAILED -- accepted NNIs do not match golden run.")

            print_v(
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
                print_v(f"  {i:<7} {res:<10} {nni:<20} {golden_nni:<20}")
        # write results
        if args.save_test:
            print_v(f"# NNI list written to: {args.save_test}")
            results.write_nni_list(args.save_test)

        return nni_inst, results, timer

    @staticmethod
    def build_credible_set(args):
        taxon_id_map, tree_nwk_map, tree_pp_map, tree_cpp_map = Loader.load_trprobs(
            args.trprobs, args)
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
    print_v("# begin...")
    args = Program.parse_args(sys.argv[1:])
    Program.run_program(args)
    print_v("# ...done")
