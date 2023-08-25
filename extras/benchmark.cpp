#include <iostream>
#include <fstream>

#include "stopwatch.hpp"
#include "gp_instance.hpp"
#include "rooted_sbn_instance.hpp"
#include "unrooted_sbn_instance.hpp"
#include "rooted_tree_collection.hpp"

// This is just a place to muck around, and check out performance.

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

auto GetArgVec(int argc, char* argv[]) {
  std::string current_exec_name = argv[0];  // Name of the current exec program
  std::vector<std::string> all_args;
  if (argc > 1) {
    all_args.assign(argv + 1, argv + argc);
  }
  return all_args;
}

auto MakeDAGInstanceFromFiles(const std::string& fasta_path,
                              const std::string& newick_path) {
  GPInstance inst("_ignore/mmap.data");
  inst.ReadFastaFile(fasta_path);
  inst.ReadNewickFile(newick_path);
  return inst;
}

void OutputNewickToFile(const std::string& file_path, const std::string& newick_str) {
  std::ofstream file_out;
  file_out.open(file_path);
  file_out << newick_str << std::endl;
  file_out.close();
}

namespace SubsplitSetBuilder {
void BuildAllSubsplitsRecurse(std::vector<int>& subsplit_assign, size_t n,
                              std::set<Bitset>& results) {
  if (n > 0) {
    for (size_t i = 0; i < 3; i++) {
      subsplit_assign[subsplit_assign.size() - 1 - n] = i;
      BuildAllSubsplitsRecurse(subsplit_assign, n - 1, results);
    }
    return;
  }
  if (n == 0) {
    Bitset clade_left(subsplit_assign.size(), false);
    Bitset clade_right(subsplit_assign.size(), false);
    for (size_t i = 0; i < subsplit_assign.size(); i++) {
      if (subsplit_assign[i] == 0) {
        clade_left.set(i);
      }
      if (subsplit_assign[i] == 1) {
        clade_right.set(i);
      }
    }
    if ((clade_left.Count() == 0) && (clade_right.Count() == 0)) {
      return;
    }
    if ((clade_left.Count() == 0) && (clade_right.Count() > 1)) {
      return;
    }
    if ((clade_left.Count() > 1) && (clade_right.Count() == 0)) {
      return;
    }
    Bitset subsplit = Bitset::Subsplit(clade_left, clade_right);
    results.insert(subsplit);
    return;
  }
}

std::set<Bitset> BuildAllSubsplits(size_t n) {
  std::vector<int> subsplit_assign(n);
  std::set<Bitset> all_subsplits;
  BuildAllSubsplitsRecurse(subsplit_assign, n, all_subsplits);
  std::cout << "all_subsplits: " << n << " " << all_subsplits.size() << " "
            << all_subsplits << std::endl;
  // return all_subsplits;
}
}  // namespace SubsplitSetBuilder

int main(int argc, char* argv[]) {
  Stopwatch timer(true, Stopwatch::TimeScale::SecondScale);

  // Parse commandline args.
  auto args = GetArgVec(argc, argv);
  if (args.size() < 2 or args.size() > 3) {
    std::cout << "usage: <fasta_path> <newick_path> <search_type_opt>" << std::endl;
    exit(0);
  }
  auto fasta_path = args[0];
  auto newick_path = args[1];

  std::cout << "Fasta: " << fasta_path << std::endl;
  std::cout << "Newick: " << newick_path << std::endl;

  // Build DAG and Engines.
  auto inst = MakeDAGInstanceFromFiles(fasta_path, newick_path);
  inst.MakeDAG();
  auto& dag = inst.GetDAG();
  inst.MakeGPEngine();
  inst.MakeTPEngine();
  inst.MakeNNIEngine();
  auto& nni_engine = inst.GetNNIEngine();
  auto& gp_engine = inst.GetGPEngine();
  auto& tp_engine = inst.GetTPEngine();
  auto& graft_dag = nni_engine.GetGraftDAG();

  nni_engine.SetTPLikelihoodCutoffFilteringScheme(0.0);
  nni_engine.SetTopNScoreFilteringScheme(5);

  std::cout << "# Build DAG and NNIEngine: " << timer.Lap() << " sec" << std::endl;

  tp_engine.OptimizeBranchLengths();
  std::cout << "# Optimize Branch Lengths: " << timer.Lap() << " sec" << std::endl;

  nni_engine.RunInit();
  std::cout << "# nni_engine.RunInit(): " << timer.Lap() << " sec" << std::endl;

  std::cout << "DAG_COUNTS: " << dag.NodeCount() << " "
            << dag.EdgeCountWithLeafSubsplits() << std::endl;
  std::cout << "TOPO_SORT: " << dag.LeafwardNodeTraversalTrace(true) << std::endl;
  std::cout << "BRANCH_LENGTHS: " << tp_engine.GetBranchLengths() << std::endl;

  Stopwatch iter_timer(true, Stopwatch::TimeScale::SecondScale);
  size_t max_iter = 5;
  for (size_t iter = 0; iter < max_iter; iter++) {
    nni_engine.SyncAdjacentNNIsWithDAG();
    std::cout << "### Iteration " << iter << " of " << max_iter << "..." << std::endl;
    std::cout << "ADJACENT_NNIs: " << nni_engine.GetAdjacentNNIs().size() << std::endl;
    // Main Loop
    nni_engine.GraftAdjacentNNIsToDAG();
    std::cout << "# nni_engine.GraftAdjacentNNIsToDAG(): " << timer.Lap() << " sec"
              << std::endl;
    std::cout << "DAG_COUNTS (INNER): " << dag.NodeCount() << " "
              << dag.EdgeCountWithLeafSubsplits() << std::endl;
    std::cout << "GRAFT_COUNTS (INNER): " << graft_dag.HostNodeCount() << " "
              << graft_dag.HostEdgeCount() << std::endl;
    std::cout << "GRAFT_COUNTS (INNER): " << graft_dag.GraftNodeCount() << " "
              << graft_dag.GraftEdgeCount() << std::endl;
    nni_engine.FilterPreUpdate();
    std::cout << "# nni_engine.FilterPreUpdate(): " << timer.Lap() << " sec"
              << std::endl;
    nni_engine.FilterEvaluateAdjacentNNIs();
    std::cout << "# nni_engine.FilterEvaluateAdjacentNNIs(): " << timer.Lap()
              << std::endl;
    nni_engine.FilterPostUpdate();
    std::cout << "# nni_engine.FilterPostUpdate(): " << timer.Lap() << " sec"
              << std::endl;
    nni_engine.FilterProcessAdjacentNNIs();
    std::cout << "# nni_engine.FilterProcessAdjacentNNIs(): " << timer.Lap() << " sec"
              << std::endl;
    nni_engine.RemoveAllGraftedNNIsFromDAG();
    std::cout << "# nni_engine.RemoveAllGraftedNNIsFromDAG(): " << timer.Lap() << " sec"
              << std::endl;
    nni_engine.AddAcceptedNNIsToDAG(false);
    std::cout << "# nni_engine.AddAcceptedNNIsToDAG(): " << timer.Lap() << " sec"
              << std::endl;
    std::cout << "NNI_SCORES: " << nni_engine.GetScoredNNIs().size() << " "
              << nni_engine.GetScoredNNIs() << std::endl;

    // Post Loop
    nni_engine.UpdateAdjacentNNIs();
    std::cout << "# nni_engine.UpdateAdjacentNNIs(): " << timer.Lap() << " sec"
              << std::endl;
    nni_engine.UpdateAcceptedNNIs();
    std::cout << "# nni_engine.UpdateAcceptedNNIs(): " << timer.Lap() << " sec"
              << std::endl;
    nni_engine.UpdateRejectedNNIs();
    std::cout << "# nni_engine.UpdateRejectedNNIs(): " << timer.Lap() << " sec"
              << std::endl;
    nni_engine.UpdateScoredNNIs();
    std::cout << "# nni_engine.UpdateScoredNNIs(): " << timer.Lap() << " sec"
              << std::endl;

    // Iteration details
    std::cout << "DAG_COUNTS (POST): " << dag.NodeCount() << " "
              << dag.EdgeCountWithLeafSubsplits() << std::endl;
    std::cout << "GRAFT_COUNTS (POST): " << graft_dag.HostNodeCount() << " "
              << graft_dag.HostEdgeCount() << std::endl;
    std::cout << "GRAFT_COUNTS (POST): " << graft_dag.GraftNodeCount() << " "
              << graft_dag.GraftEdgeCount() << std::endl;
    std::cout << "### iter_time " << iter << ": " << iter_timer.Lap() << " sec"
              << std::endl;
  }
}
