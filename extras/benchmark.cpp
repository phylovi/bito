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

void VerifyAndOutputDAGToTrees(GPInstance& inst, bool verify = false) {
  Stopwatch timer(true, Stopwatch::TimeScale::SecondScale);
  std::string all_newick_path = "_ignore/all_dag.nwk";
  std::string spanning_newick_path = "_ignore/spanning_dag.nwk";
  std::string tp_newick_path = "_ignore/tp_dag.nwk";

  auto& dag = inst.GetDAG();
  auto& tp_engine = inst.GetTPEngine();
  OutputNewickToFile(all_newick_path, dag.ToNewickOfAllTopologies());
  OutputNewickToFile(spanning_newick_path, dag.ToNewickOfSpanningTopologies());
  OutputNewickToFile(tp_newick_path, tp_engine.ToNewickOfTopTopologies());

  std::cout << "# building spanning topologies..." << std::endl;
  auto spanning_topologies = dag.GenerateSpanningTopologies();
  std::cout << timer.Lap() << std::endl;
  std::cout << "SpanningTopologies: " << spanning_topologies.size() << std::endl;
  std::cout << "# building all topologies..." << std::endl;
  auto all_topologies = dag.GenerateAllTopologies();
  std::cout << timer.Lap() << std::endl;
  std::cout << "AllTopologies: " << all_topologies.size() << std::endl;
  std::cout << "# building top topologies..." << std::endl;
  auto top_topologies = tp_engine.BuildMapOfTreeIdToTopTopologies();
  std::cout << timer.Lap() << std::endl;
  size_t top_topologies_size = 0;
  for (const auto& [tree_id, top_topo_vec] : top_topologies) {
    top_topologies_size += top_topo_vec.size();
  }
  std::cout << "TopTopologies: " << top_topologies_size << std::endl;
  std::cout << std::endl;
  // std::cout << "==> DAG::AllTopologies:" << std::endl
  //           << dag.ToNewickOfAllTopologies() << std::endl;
  // std::cout << "==> DAG::SpanningTopologies:" << std::endl
  //           << dag.ToNewickOfSpanningTopologies() << std::endl;
  // std::cout << "==> TPEngine::TopTopologies:" << std::endl
  //           << tp_engine.ToNewickOfTopTopologies() << std::endl;

  if (!verify) return;
  for (const auto& newick_path :
       {all_newick_path, spanning_newick_path, tp_newick_path}) {
    auto tmp_inst = MakeDAGInstanceFromFiles(inst.GetFastaSourcePath(), newick_path);
    tmp_inst.MakeDAG();
    auto& tmp_dag = tmp_inst.GetDAG();

    bool dag_equal = (dag.Compare(tmp_dag, false) == 0);
    std::cout << "DAG Compare " << (dag_equal ? "PASS" : "FAIL") << ": " << newick_path
              << std::endl;
    if (!dag_equal) {
      std::cout << "node_counts: " << dag.NodeCount() << " " << tmp_dag.NodeCount()
                << std::endl;
      std::cout << "edge_counts: " << dag.EdgeCountWithLeafSubsplits() << " "
                << tmp_dag.EdgeCountWithLeafSubsplits() << std::endl;
      std::cout << "DAG_A Taxon Map: " << dag.GetTaxonMap() << std::endl;
      std::cout << "DAG_B Taxon Map: " << tmp_dag.GetTaxonMap() << std::endl;
      std::cout << "Taxon Translation Map: "
                << SubsplitDAG::BuildTaxonTranslationMap(dag, tmp_dag) << std::endl;
    }
    if (!dag_equal) {
      auto [common, not_in_rhs, not_in_lhs] =
          SubsplitDAG::CompareSubsplits(dag, tmp_dag);
      std::cout << "common: " << common.size() << std::endl;
      std::cout << "lhs_not_in_rhs: " << not_in_rhs.size() << " " << not_in_rhs
                << std::endl;
      std::cout << "rhs_not_in_lhs: " << not_in_lhs.size() << " " << not_in_lhs
                << std::endl;
    }
    if (!dag_equal) {
      auto [common, not_in_rhs, not_in_lhs] = SubsplitDAG::ComparePCSPs(dag, tmp_dag);
      std::cout << "common: " << common.size() << std::endl;
      std::cout << "lhs_not_in_rhs: " << not_in_rhs.size() << " " << not_in_rhs
                << std::endl;
      std::cout << "rhs_not_in_lhs: " << not_in_lhs.size() << " " << not_in_lhs
                << std::endl;
    }
  }
  std::cout << std::endl;
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
  auto args = GetArgVec(argc, argv);
  if (args.size() < 2 or args.size() > 3) {
    std::cout << "usage: <fasta_path> <newick_path> <search_type_opt>" << std::endl;
    exit(0);
  }

  auto fasta_path = args[0];
  auto newick_path = args[1];

  bool use_gp = false;
  if (args.size() == 3) {
    if (args[2] == "gp") {
      use_gp = true;
    }
    if (args[2] == "tp") {
      use_gp = false;
    }
  }

  std::cout << "Fasta: " << fasta_path << std::endl;
  std::cout << "Newick: " << newick_path << std::endl;

  auto inst = MakeDAGInstanceFromFiles(fasta_path, newick_path);
  inst.MakeDAG();
  auto& dag = inst.GetDAG();
  inst.MakeGPEngine();
  inst.MakeTPEngine();
  inst.MakeNNIEngine();
  auto& nni_engine = inst.GetNNIEngine();
  auto& gp_engine = inst.GetGPEngine();
  auto& tp_engine = inst.GetTPEngine();

  // inst.TakeFirstBranchLength();
  auto& gp_bls = gp_engine.GetBranchLengthHandler();
  auto& tp_bls = tp_engine.GetLikelihoodEvalEngine().GetDAGBranchHandler();

  std::cout << "Initial DAG: " << dag.TaxonCount() << " " << dag.NodeCount() << " "
            << dag.EdgeCountWithLeafSubsplits() << std::endl;
  std::cout << "[before optim]" << std::endl;
  std::cout << "TPEngine::BranchLengths: " << tp_bls.GetBranchLengthData() << std::endl;
  std::cout << "GPEngine::BranchLengths: " << gp_bls.GetBranchLengthData() << std::endl;

  size_t optimization_count = 2;
  inst.EstimateBranchLengths(1e-5, optimization_count, true);
  tp_engine.GetLikelihoodEvalEngine().SetOptimizationMaxIteration(optimization_count);
  tp_engine.GetLikelihoodEvalEngine().BranchLengthOptimization(false);
  if (use_gp) {
    nni_engine.SetGPLikelihoodCutoffFilteringScheme(0.0);
  } else {
    nni_engine.SetTPLikelihoodCutoffFilteringScheme(0.0);
  }
  nni_engine.SetTopNScoreFilteringScheme(5);
  nni_engine.SetReevaluateRejectedNNIs(true);
  nni_engine.RunInit(false);

  std::cout << "[after optim]" << std::endl;
  std::cout << "TPEngine::BranchLengths: " << tp_bls.GetBranchLengthData() << std::endl;
  std::cout << "GPEngine::BranchLengths: " << gp_bls.GetBranchLengthData() << std::endl;

  for (EdgeId edge_id{0}; edge_id < dag.EdgeCountWithLeafSubsplits(); edge_id++) {
    std::cout << "Edge" << edge_id << " " << dag.GetDAGEdgeBitset(edge_id) << ": "
              << tp_bls(edge_id) << " " << gp_bls(edge_id) << " "
              << dag.IsEdgeLeaf(edge_id) << std::endl;
  }

  std::vector<NNIOperation> added_nnis;
  size_t iter_max = 2;
  for (size_t iter = 0; iter < iter_max; iter++) {
    std::cout << "Iteration: " << iter << " of " << iter_max << std::endl;

    std::cout << "Newick: " << tp_engine.ToNewickOfTopTrees() << std::endl;
    std::ofstream fp;
    fp.open("_ignore/test_dag.nwk");
    fp << tp_engine.ToNewickOfTopTrees() << std::endl;
    fp.close();

    std::cout << "ChoiceMap: " << tp_engine.GetChoiceMap().ToString() << std::endl;

    // std::cout << "DAG Counts: " << dag.NodeCount() << " "
    //           << dag.EdgeCountWithLeafSubsplits() << std::endl;
    // if (!use_gp) {
    //   std::cout << "TP Engine Top Trees: " << std::endl;
    //   std::cout << tp_engine.ToNewickOfTopTrees() << std::endl;
    // }
    // std::cout << "DAG Spanning Trees: " << std::endl;
    // const auto trees = dag.GenerateSpanningTrees(gp_engine.GetBranchLengths());
    // for (const auto tree : trees) {
    //   std::cout << tree.Newick() << std::endl;
    // }

    // if (use_gp) {
    //   inst.EstimateBranchLengths(1e-5, 1, false);
    // }
    // auto& pvs = gp_engine.GetPLVHandler();

    // auto NodeData = [&dag, &pvs](NodeId node_id) {
    //   const auto& node = dag.GetDAGNode(node_id);
    //   std::cout << "Subsplit: " << dag.GetDAGNodeBitset(node_id).SubsplitToString()
    //             << std::endl;
    //   for (const auto dir : DirectionEnum::Iterator()) {
    //     std::cout << ((dir == Direction::Rootward) ? "Parent" : "Child") << ": ";
    //     for (const auto clade : SubsplitCladeEnum::Iterator()) {
    //       const auto neighbors_view = node.GetNeighbors(dir, clade);
    //       std::set<NodeId> neighbors{neighbors_view.begin(), neighbors_view.end()};
    //       std::cout << neighbors << " ";
    //     }
    //   }
    //   std::cout << std::endl;
    //   std::cout << "Node" << node_id << std::endl;
    //   for (const auto pv_type : PLVTypeEnum::Iterator()) {
    //     auto range = pvs.ValueRange(pvs.GetPVIndex(pv_type, node_id));
    //     std::cout << PLVTypeEnum::ToString(pv_type) << "_" << range << " ";
    //   }
    //   std::cout << std::endl;
    // };

    nni_engine.RunMainLoop(false);
    std::cout << "Scored NNIs: " << nni_engine.GetScoredNNIs().size() << std::endl;
    for (const auto& [nni, score] : nni_engine.GetScoredNNIs()) {
      std::cout << "\t" << nni << ": " << score << std::endl;
    }
    for (const auto& nni : nni_engine.GetAcceptedNNIs()) {
      added_nnis.push_back(nni);
    }
    nni_engine.RunPostLoop(false);
    std::cout << "Added NNIs: " << added_nnis.size() << std::endl;
    for (const auto& nni : added_nnis) {
      std::cout << "\t" << nni << std::endl;
    }
  }
}
