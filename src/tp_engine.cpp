// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//

#include "tp_engine.hpp"
#include "tp_evaluation_engine.hpp"
#include "gp_engine.hpp"
#include "sbn_maps.hpp"
#include "optimization.hpp"
#include "stopwatch.hpp"

TPEngine::TPEngine(GPDAG &dag, SitePattern &site_pattern,
                   std::optional<std::string> mmap_likelihood_path,
                   std::optional<std::string> mmap_parsimony_path,
                   std::optional<const RootedTreeCollection> tree_collection,
                   std::optional<const BitsetSizeMap> edge_indexer)
    : choice_map_(dag), site_pattern_(site_pattern), dag_(&dag) {
  // Initialize site pattern-based data.
  auto weights = site_pattern_.GetWeights();
  site_pattern_weights_ = EigenVectorXdOfStdVectorDouble(weights);
  // Initialize node-based data
  GrowNodeData(GetDAG().NodeCount(), std::nullopt, std::nullopt, true);
  // Initialize edge-based data
  GrowEdgeData(GetDAG().EdgeCountWithLeafSubsplits(), std::nullopt, std::nullopt, true);
  // Initialize Tree Source.
  if (tree_collection.has_value() && edge_indexer.has_value()) {
    SetTreeSourceByTakingFirst(tree_collection.value(), edge_indexer.value());
  } else {
    std::fill(GetTreeSource().begin(), GetTreeSource().end(), TreeId{1});
  }
  // Initialize Choice Map.
  InitializeChoiceMap();
  // Initialize Eval Engines
  if (mmap_likelihood_path.has_value()) {
    MakeLikelihoodEvalEngine(mmap_likelihood_path.value());
  }
  if (mmap_parsimony_path.has_value()) {
    MakeParsimonyEvalEngine(mmap_parsimony_path.value());
  }
  InitializeScores();
}

// ** Comparators

int TPEngine::Compare(const TPEngine &lhs, const TPEngine &rhs, const bool is_quiet) {
  std::stringstream dev_null;
  auto &os = is_quiet ? dev_null : std::cout;

  if (lhs.GetDAG() != rhs.GetDAG()) {
    os << "TPEngine::Compare: DAGs are not equal." << std::endl;
    return SubsplitDAG::Compare(lhs.GetDAG(), rhs.GetDAG(), is_quiet);
  }
  const auto node_map =
      SubsplitDAG::BuildNodeIdMapBetweenDAGs(lhs.GetDAG(), rhs.GetDAG());
  if (node_map.size() != lhs.GetDAG().NodeCount()) {
    os << "TPEngine::Compare: node_maps do not cover entire node set: "
       << node_map.size() << " " << lhs.GetDAG().NodeCount() << " "
       << rhs.GetDAG().NodeCount() << std::endl;
    return node_map.size() - lhs.GetDAG().NodeCount();
  }
  const auto edge_map =
      SubsplitDAG::BuildEdgeIdMapBetweenDAGs(lhs.GetDAG(), rhs.GetDAG());
  if (edge_map.size() != lhs.GetDAG().EdgeCountWithLeafSubsplits()) {
    os << "TPEngine::Compare: edge_maps do not cover entire edge set: "
       << edge_map.size() << " " << lhs.GetDAG().EdgeCountWithLeafSubsplits() << " "
       << rhs.GetDAG().EdgeCountWithLeafSubsplits() << std::endl;
    return edge_map.size() - lhs.GetDAG().EdgeCountWithLeafSubsplits();
  }

  auto TranslateEdgeChoice =
      [&edge_map](const TPChoiceMap::EdgeChoice &pre_edge_choice) {
        TPChoiceMap::EdgeChoice edge_choice(pre_edge_choice);
        if (pre_edge_choice.parent_edge_id != NoId) {
          edge_choice.parent_edge_id =
              edge_map.find(pre_edge_choice.parent_edge_id)->second;
        }
        if (pre_edge_choice.sister_edge_id != NoId) {
          edge_choice.sister_edge_id =
              edge_map.find(pre_edge_choice.sister_edge_id)->second;
        }
        if (pre_edge_choice.left_child_edge_id != NoId) {
          edge_choice.left_child_edge_id =
              edge_map.find(pre_edge_choice.left_child_edge_id)->second;
        }
        if (pre_edge_choice.right_child_edge_id != NoId) {
          edge_choice.right_child_edge_id =
              edge_map.find(pre_edge_choice.right_child_edge_id)->second;
        }

        return edge_choice;
      };

  auto CompareEdgeChoice = [](const TPChoiceMap::EdgeChoice &lhs,
                              const TPChoiceMap::EdgeChoice &rhs) {
    if (lhs.parent_edge_id != rhs.parent_edge_id) {
      return int(lhs.parent_edge_id.value_ - rhs.parent_edge_id.value_);
    }
    if (lhs.sister_edge_id != rhs.sister_edge_id) {
      return int(lhs.sister_edge_id.value_ - rhs.sister_edge_id.value_);
    }
    if ((lhs.left_child_edge_id != rhs.left_child_edge_id) or
        (lhs.right_child_edge_id != rhs.left_child_edge_id)) {
      // Right and left child can be swapped due to taxon bitset ordering.
      if ((lhs.left_child_edge_id == rhs.right_child_edge_id) and
          (lhs.right_child_edge_id == rhs.left_child_edge_id)) {
        return 0;
      }
      if (lhs.left_child_edge_id != rhs.left_child_edge_id) {
        return int(lhs.left_child_edge_id.value_ - rhs.left_child_edge_id.value_);
      }
      return int(lhs.right_child_edge_id.value_ - rhs.right_child_edge_id.value_);
    }
    return 0;
  };

  int final_diff = 0;
  for (const auto [lhs_edge_id, rhs_edge_id] : edge_map) {
    const auto &lhs_choice = lhs.GetChoiceMap(lhs_edge_id);
    const auto &rhs_choice = rhs.GetChoiceMap(rhs_edge_id);
    auto trans_choice = TranslateEdgeChoice(lhs_choice);
    auto choice_diff = CompareEdgeChoice(trans_choice, rhs_choice);
    if (choice_diff != 0) {
      os << "TPEngine::Compare: edge_choice does not match at -- " << lhs_edge_id << " "
         << rhs_edge_id << std::endl;
      os << "lhs_choice_map: " << lhs_edge_id << " " << lhs_choice << std::endl;
      os << "trans_choice_map: " << lhs_edge_id << " " << trans_choice << std::endl;
      os << "rhs_choice_map: " << rhs_edge_id << " " << rhs_choice << std::endl;

      for (const auto dir : DirectionEnum::Iterator()) {
        auto &hs = lhs;
        auto &hs2 = rhs;
        auto edge_id = lhs_edge_id;
        const auto node_id = (dir == Direction::Leafward)
                                 ? hs.GetDAG().GetDAGEdge(edge_id).GetChild()
                                 : hs.GetDAG().GetDAGEdge(edge_id).GetParent();
        const auto &node = hs.GetDAG().GetDAGNode(node_id);
        for (const auto clade : SubsplitCladeEnum::Iterator()) {
          const auto &adj_node_ids = node.GetNeighbors(dir, clade);
          os << "hs::" << DirectionEnum::ToString(dir) << ":"
             << SubsplitCladeEnum::ToString(clade) << ": ";
          for (const auto adj_node_id : adj_node_ids) {
            auto edge_id = (dir == Direction::Leafward)
                               ? hs.GetDAG().GetEdgeIdx(node.Id(), adj_node_id)
                               : hs.GetDAG().GetEdgeIdx(adj_node_id, node.Id());
            os << "Edge" << edge_id << "->"
               << "Tree" << hs.GetTreeSource(edge_id) << " => ";
            os << "Edge" << edge_map.find(edge_id)->second << "->"
               << "Tree" << hs2.GetTreeSource(edge_map.find(edge_id)->second) << ", ";
          }
          os << std::endl;
        }
      }

      final_diff = choice_diff;
    }
  }
  if (lhs.HasLikelihoodEvalEngine() && rhs.HasLikelihoodEvalEngine()) {
    auto branch_diff =
        DAGBranchHandler::Compare(lhs.GetLikelihoodEvalEngine().GetDAGBranchHandler(),
                                  rhs.GetLikelihoodEvalEngine().GetDAGBranchHandler());
    if (branch_diff != 0) {
      os << "TPEngine::Compare: Branch lengths do not match" << std::endl;
      final_diff = branch_diff;
    }
  }

  return final_diff;
}

bool operator==(const TPEngine &lhs, const TPEngine &rhs) {
  return TPEngine::Compare(lhs, rhs) == 0;
}

// ** Maintenance

void TPEngine::Initialize() {
  // Initialize node-based data
  GrowNodeData(GetDAG().NodeCount(), std::nullopt, std::nullopt, true);
  // Initialize edge-based data
  GrowEdgeData(GetDAG().EdgeCountWithLeafSubsplits(), std::nullopt, std::nullopt, true);
  // Initialize scores.
  InitializeChoiceMap();
}

void TPEngine::UpdateAfterModifyingDAG(
    const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni,
    const size_t prev_node_count, const Reindexer &node_reindexer,
    const size_t prev_edge_count, const Reindexer &edge_reindexer) {
  Stopwatch timer(true, Stopwatch::TimeScale::SecondScale);
  UpdateChoiceMapAfterModifyingDAG(nni_to_pre_nni, prev_node_count, node_reindexer,
                                   prev_edge_count, edge_reindexer);
  std::cout << "UpdateAfterModifying::ChoiceMap: " << timer.Lap() << std::endl;
  if (IsEvalEngineInUse(TPEvalEngineType::LikelihoodEvalEngine)) {
    GetLikelihoodEvalEngine().UpdateEngineAfterModifyingDAG(
        nni_to_pre_nni, prev_node_count, node_reindexer, prev_edge_count,
        edge_reindexer);
    std::cout << "UpdateAfterModifying::LikelihoodEvalEngine: " << timer.Lap()
              << std::endl;
  }
  if (IsEvalEngineInUse(TPEvalEngineType::ParsimonyEvalEngine)) {
    GetParsimonyEvalEngine().UpdateEngineAfterModifyingDAG(
        nni_to_pre_nni, prev_node_count, node_reindexer, prev_edge_count,
        edge_reindexer);
    std::cout << "UpdateAfterModifying::ParsimonyEvalEngine: " << timer.Lap()
              << std::endl;
  }
}

// ** Tree/Topology Builder

Node::Topology TPEngine::GetTopTopologyWithEdge(const EdgeId edge_id) const {
  auto topology = GetChoiceMap().ExtractTopology(edge_id);
  return topology;
}

RootedTree TPEngine::GetTopTreeWithEdge(const EdgeId edge_id) const {
  Assert(HasLikelihoodEvalEngine(),
         "Must MakeLikelihoodEvalEngine before getting top tree.");
  auto topology = GetTopTopologyWithEdge(edge_id);
  return BuildTreeFromTopologyInDAG(topology);
}

RootedTree TPEngine::BuildTreeFromTopologyInDAG(const Node::Topology &topology) const {
  auto rooted_tree = GetDAG().BuildTreeFromTopology(
      topology, GetDAGBranchHandler().GetBranchLengthData());
  return rooted_tree;
}

std::set<EdgeId> TPEngine::BuildSetOfEdgesRepresentingTopology(
    const Node::Topology &topology) const {
  std::set<EdgeId> edges;
  topology->Preorder([this, &edges](const Node *node) {
    if (!node->IsLeaf()) {
      for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
        const auto edge_id = GetDAG().GetEdgeIdx(node->BuildPCSP(clade));
        edges.insert(edge_id);
      }
    }
  });
  return edges;
}

std::set<TreeId> TPEngine::FindTreeIdsInTreeEdgeVector(
    const std::set<EdgeId> edge_ids) const {
  Assert(!edge_ids.empty(), "EdgeVector representation of tree cannot be empty.");
  std::set<TreeId> tree_ids;
  for (const auto edge_id : edge_ids) {
    tree_ids.insert(GetTreeSource()[edge_id.value_]);
  }
  return tree_ids;
}

TPEngine::EdgeIdTopologyMap TPEngine::BuildMapOfEdgeIdToTopTopologies() const {
  EdgeIdTopologyMap topology_map;
  BoolVector visited_edges(GetEdgeCount(), false);
  // Ignore rootsplit edges.
  for (const auto edge_id : GetDAG().GetRootsplitEdgeIds()) {
    visited_edges[edge_id.value_] = true;
  }
  // Build trees as we encounter edges not yet assigned to a tree.
  for (EdgeId edge_id(0); edge_id < GetDAG().EdgeCountWithLeafSubsplits(); edge_id++) {
    if (visited_edges[edge_id.value_]) {
      continue;
    }
    const auto top_topology = GetTopTopologyWithEdge(edge_id);
    const auto tree_visited_edges = BuildSetOfEdgesRepresentingTopology(top_topology);
    std::set<EdgeId> new_visited_edges;
    // For edges contained in the tree, add them if edge has not already been assigned
    // to previous tree.
    for (const auto edge_id : tree_visited_edges) {
      if (!visited_edges[edge_id.value_]) {
        new_visited_edges.insert(edge_id);
      }
      visited_edges[edge_id.value_] = true;
    }
    topology_map.push_back({new_visited_edges, top_topology});
  }
  if (!std::all_of(visited_edges.begin(), visited_edges.end(),
                   [](bool n) { return n; })) {
    Failwith("Top Topologies does not cover entire DAG.");
  }
  return topology_map;
}

TPEngine::TreeIdTopologyMap TPEngine::BuildMapOfTreeIdToTopTopologies() const {
  TreeIdTopologyMap topology_map;
  const auto edge_topology_map = BuildMapOfEdgeIdToTopTopologies();
  for (const auto &[edge_ids, topology] : edge_topology_map) {
    const auto tree_ids = FindTreeIdsInTreeEdgeVector(edge_ids);
    const auto tree_id = *std::min_element(tree_ids.begin(), tree_ids.end());
    if (topology_map.find(tree_id) == topology_map.end()) {
      topology_map[tree_id] = std::vector<Node::Topology>();
    }
    topology_map[tree_id].push_back(topology);
  }
  return topology_map;
}

TPEngine::TreeIdTreeMap TPEngine::BuildMapOfTreeIdToTopTrees() const {
  TreeIdTreeMap tree_map;
  const auto topology_map = BuildMapOfTreeIdToTopTopologies();
  for (const auto &[tree_id, topology_vec] : topology_map) {
    tree_map[tree_id] = std::vector<RootedTree>();
    for (const auto topology : topology_vec) {
      RootedTree tree = BuildTreeFromTopologyInDAG(topology);
      tree_map[tree_id].push_back(tree);
    }
  }
  return tree_map;
}

std::string TPEngine::ToNewickOfTopTopologies() const {
  std::stringstream str;
  const auto topology_map = BuildMapOfTreeIdToTopTopologies();
  for (TreeId tree_id(0); tree_id < GetMaxTreeId(); tree_id++) {
    if (topology_map.find(tree_id) == topology_map.end()) continue;
    const auto &topology_vec = topology_map.find(tree_id)->second;
    for (const auto topology : topology_vec) {
      str << topology->Newick(std::nullopt, GetDAG().GetTagTaxonMap()) << std::endl;
    }
  }
  return str.str();
}

std::string TPEngine::ToNewickOfTopTrees() const {
  std::stringstream str;
  const auto tree_map = BuildMapOfTreeIdToTopTrees();
  for (TreeId tree_id(0); tree_id < GetMaxTreeId(); tree_id++) {
    if (tree_map.find(tree_id) == tree_map.end()) continue;
    const auto &tree_vec = tree_map.find(tree_id)->second;
    for (const auto tree : tree_vec) {
      str << tree.Newick(GetDAG().GetTagTaxonMap()) << std::endl;
    }
  }
  return str.str();
}

TPChoiceMap::EdgeChoicePCSPs TPEngine::BuildAdjacentPCSPsToProposedNNI(
    const NNIOperation &nni, const TPChoiceMap::EdgeChoiceNodeIds &adj_node_ids) const {
  TPChoiceMap::EdgeChoicePCSPs adj_pcsps;
  const auto &parent_subsplit = GetDAG().GetDAGNodeBitset(adj_node_ids.parent_node_id);
  adj_pcsps.parent_pcsp = Bitset::PCSP(parent_subsplit, nni.GetParent());
  adj_pcsps.focal_pcsp = Bitset::PCSP(nni.GetParent(), nni.GetChild());
  const auto &sister_subsplit = GetDAG().GetDAGNodeBitset(adj_node_ids.sister_node_id);
  adj_pcsps.sister_pcsp = Bitset::PCSP(nni.GetParent(), sister_subsplit);
  const auto &leftchild_subsplit =
      GetDAG().GetDAGNodeBitset(adj_node_ids.left_child_node_id);
  adj_pcsps.left_child_pcsp = Bitset::PCSP(nni.GetChild(), leftchild_subsplit);
  const auto &rightchild_subsplit =
      GetDAG().GetDAGNodeBitset(adj_node_ids.right_child_node_id);
  adj_pcsps.right_child_pcsp = Bitset::PCSP(nni.GetChild(), rightchild_subsplit);
  return adj_pcsps;
}

// ** Choice Map / Tree Source

void TPEngine::InitializeChoiceMap() {
  for (EdgeId edge_id = EdgeId(0); edge_id < GetEdgeCount(); edge_id++) {
    UpdateEdgeChoiceByTakingHighestPriorityTree(edge_id);
  }
}

void TPEngine::UpdateChoiceMapAfterModifyingDAG(
    const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni,
    const size_t prev_node_count, const Reindexer &node_reindexer,
    const size_t prev_edge_count, const Reindexer &edge_reindexer) {
  // std::cout << "UpdateChoiceMap [begin]" << std::endl;
  Stopwatch full_timer(true, Stopwatch::TimeScale::SecondScale);
  Stopwatch timer(true, Stopwatch::TimeScale::SecondScale);

  const size_t new_edge_count = edge_reindexer.size();
  // Tree counter for assigning a unique tree id for each NNI added to the DAG.
  tree_counter_++;
  const TreeId min_tree_id(tree_counter_);
  const TreeId max_tree_id(tree_counter_ + nni_to_pre_nni.size());
  TreeId nni_tree_id(tree_counter_);

  // Edges that have been newly added to the DAG.
  std::set<EdgeId> new_edges, nni_edges;
  for (size_t i = prev_edge_count; i < new_edge_count; i++) {
    const EdgeId edge_id = EdgeId(edge_reindexer.GetNewIndexByOldIndex(i));
    new_edges.insert(edge_id);
  }
  // Edges that need their choicemaps initialized.
  std::set<EdgeId> edges_to_init_choicemap(new_edges);

  // Initialize all new tree sources to default tree source.
  for (const auto edge_id : new_edges) {
    GetTreeSource(edge_id) = max_tree_id;
    GetChoiceMap().ResetEdgeChoice(edge_id);
    GetLikelihoodEvalEngine().GetDAGBranchHandler()(edge_id) =
        GetLikelihoodEvalEngine().GetDAGBranchHandler().GetDefaultBranchLength();
  }

  // Build map of best reference edges to update.
  NNISet nnis;
  for (const auto &[post_nni, pre_nni] : nni_to_pre_nni) {
    std::ignore = pre_nni;
    nnis.insert(post_nni);
  }
  const auto best_pcsp_edge_map =
      BuildBestEdgeMapOverNNIs(nnis, prev_edge_count, edge_reindexer);
  std::unordered_map<EdgeId, EdgeId> best_edge_map;
  for (const auto &[pcsp, pre_edge_id] : best_pcsp_edge_map) {
    if (!GetDAG().ContainsEdge(pcsp)) {
      std::cout << "PCSP not found in DAG: " << pcsp.PCSPToString() << std::endl;
      std::cout << "NNIs: " << std::endl;
      for (const auto &[post_nni, pre_nni] : nni_to_pre_nni) {
        std::cout << post_nni.GetCentralEdgePCSP().PCSPToString() << std::endl;
      }
    }
    const auto post_edge_id = GetDAG().GetEdgeIdx(pcsp);
    best_edge_map[post_edge_id] = pre_edge_id;
  }

  // Update branch lengths from best reference edges.
  if (HasLikelihoodEvalEngine()) {
    auto &branch_handler = GetLikelihoodEvalEngine().GetDAGBranchHandler();
    for (const auto &[post_edge_id, pre_edge_id] : best_edge_map) {
      branch_handler(post_edge_id) = branch_handler(pre_edge_id);
    }
  }

  // Iterate over newly added NNIs to initialize choice map and tree source.
  for (const auto &[post_nni, pre_nni] : nni_to_pre_nni) {
    std::ignore = pre_nni;
    // Get edge mapping from pre-NNI to post-NNI.
    const auto post_edge_id = GetDAG().GetEdgeIdx(post_nni);
    nni_edges.insert(post_edge_id);
    edges_to_init_choicemap.erase(post_edge_id);
    const auto mapped_post_choice =
        GetRemappedEdgeChoiceFromPreNNIToPostNNI(pre_nni, post_nni);

    // Use edge mapping to update each edge.
    auto UpdateEdge = [this, &nni_tree_id, &new_edges](const EdgeId post_edge_id) {
      // Update tree source.
      if (GetTreeSource(post_edge_id) > nni_tree_id) {
        GetTreeSource(post_edge_id) = nni_tree_id;
      }
    };
    UpdateEdge(post_edge_id);
    UpdateEdge(mapped_post_choice.parent_edge_id);
    UpdateEdge(mapped_post_choice.sister_edge_id);
    UpdateEdge(mapped_post_choice.left_child_edge_id);
    UpdateEdge(mapped_post_choice.right_child_edge_id);
    // Use edge mapping to update each choice map.
    GetChoiceMap(post_edge_id) = mapped_post_choice;

    // Increment source tree priority.
    nni_tree_id++;
    tree_counter_++;
  }

  // For non-central edges, initialize their choicemap.
  for (const auto edge_id : edges_to_init_choicemap) {
    UpdateEdgeChoiceByTakingHighestPriorityTree(edge_id);
    GetTreeSource(edge_id) = nni_tree_id;

    // Increment source tree priority.
    nni_tree_id++;
    tree_counter_++;
  }
  // Update adjacent edge choice.
  for (const auto &[post_nni, pre_nni] : nni_to_pre_nni) {
    std::ignore = pre_nni;
    using AdjacentEdge = TPChoiceMap::AdjacentEdge;
    // Get edge mapping from pre-NNI to post-NNI.
    const auto post_edge_id = GetDAG().GetEdgeIdx(post_nni);
    // Update given choice with given adj_edge_id if edge is new.
    auto UpdateChoice = [this, &new_edges](EdgeId choice_edge_id,
                                           AdjacentEdge edge_type, EdgeId adj_edge_id) {
      if (new_edges.find(choice_edge_id) != new_edges.end()) {
        GetChoiceMap().SetEdgeChoice(choice_edge_id, edge_type, adj_edge_id);
      }
    };
    const auto &choice = GetChoiceMap(post_edge_id);
    const auto focal_clade = GetDAG().GetFocalClade(post_edge_id);
    if (focal_clade == SubsplitClade::Left) {
      UpdateChoice(choice.parent_edge_id, AdjacentEdge::LeftChild, post_edge_id);
    } else {
      UpdateChoice(choice.parent_edge_id, AdjacentEdge::RightChild, post_edge_id);
    }
    UpdateChoice(choice.sister_edge_id, AdjacentEdge::Sister, post_edge_id);
    UpdateChoice(choice.left_child_edge_id, AdjacentEdge::Parent, post_edge_id);
    UpdateChoice(choice.right_child_edge_id, AdjacentEdge::Parent, post_edge_id);
  }

  // std::cout << "Update ChoiceMap [end]: " << timer.Lap() << std::endl;
}

void TPEngine::UpdateEdgeChoiceByTakingFirstTree(const EdgeId edge_id) {
  using AdjacentEdge = TPChoiceMap::AdjacentEdge;
  const auto edge = GetDAG().GetDAGEdge(EdgeId(edge_id));
  GetChoiceMap().ResetEdgeChoice(edge_id);

  auto GetFirstEdgeId = [this](const NodeId node_id, const Direction direction,
                               const SubsplitClade clade) {
    auto adj_edge_id = EdgeId(NoId);
    const auto &node = GetDAG().GetDAGNode(node_id);
    for (const auto adj_node_id : node.GetNeighbors(direction, clade)) {
      const auto parent_node_id =
          (direction == Direction::Rootward) ? adj_node_id : node_id;
      const auto child_node_id =
          (direction == Direction::Rootward) ? node_id : adj_node_id;
      adj_edge_id = GetDAG().GetEdgeIdx(parent_node_id, child_node_id);
      return adj_edge_id;
    }
    return adj_edge_id;
  };

  // Select parent.
  auto first_edge_id = EdgeId(NoId);
  for (const auto clade : SubsplitCladeEnum::Iterator()) {
    first_edge_id = GetFirstEdgeId(edge.GetParent(), Direction::Rootward, clade);
    if (first_edge_id != NoId) {
      break;
    }
  }
  GetChoiceMap().SetEdgeChoice(edge_id, AdjacentEdge::Parent, first_edge_id);
  // Select sister.
  first_edge_id = GetFirstEdgeId(edge.GetParent(), Direction::Leafward,
                                 Bitset::Opposite(edge.GetSubsplitClade()));
  GetChoiceMap().SetEdgeChoice(edge_id, AdjacentEdge::Sister, first_edge_id);
  // Select left child.
  first_edge_id =
      GetFirstEdgeId(edge.GetChild(), Direction::Leafward, SubsplitClade::Left);
  GetChoiceMap().SetEdgeChoice(edge_id, AdjacentEdge::LeftChild, first_edge_id);
  // Select right child.
  first_edge_id =
      GetFirstEdgeId(edge.GetChild(), Direction::Leafward, SubsplitClade::Right);
  GetChoiceMap().SetEdgeChoice(edge_id, AdjacentEdge::RightChild, first_edge_id);
}

void TPEngine::UpdateEdgeChoiceByTakingHighestPriorityTree(const EdgeId edge_id) {
  using AdjacentEdge = TPChoiceMap::AdjacentEdge;
  const auto edge = GetDAG().GetDAGEdge(EdgeId(edge_id));
  // GetChoiceMap().ResetEdgeChoice(edge_id);

  auto GetBestEdgeIdByHighestPriorityTree =
      [this, edge_id](const NodeId node_id, const Direction direction,
                      const SubsplitClade clade, TreeId *opt_tree_id = nullptr) {
        auto best_tree_id = TreeId(NoId);
        auto best_edge_id = EdgeId(NoId);
        if (!GetDAG().ContainsNode(node_id)) {
          std::cout << "Encountered invalid node id during edge choice update: Edge"
                    << edge_id << " Node" << node_id << std::endl;
        }
        const auto &node = GetDAG().GetDAGNode(node_id);
        bool has_first_edge = false;
        for (const auto adj_node_id : node.GetNeighbors(direction, clade)) {
          const auto parent_node_id =
              (direction == Direction::Rootward) ? adj_node_id : node_id;
          const auto child_node_id =
              (direction == Direction::Rootward) ? node_id : adj_node_id;
          const auto adj_edge_id = GetDAG().GetEdgeIdx(parent_node_id, child_node_id);
          auto adj_tree_id = GetTreeSource()[adj_edge_id.value_];
          if ((best_tree_id > adj_tree_id) || !has_first_edge) {
            best_tree_id = adj_tree_id;
            best_edge_id = adj_edge_id;
            has_first_edge = true;
          }
        }
        if (opt_tree_id) {
          *opt_tree_id = best_tree_id;
        }
        return best_edge_id;
      };

  auto best_edge_id = EdgeId(NoId);
  auto best_tree_id = TreeId(NoId);
  // Select parent.
  // std::cout << "select_parent" << std::endl;
  for (const auto clade : SubsplitCladeEnum::Iterator()) {
    TreeId clade_tree_id = TreeId(NoId);
    const auto clade_edge_id = GetBestEdgeIdByHighestPriorityTree(
        edge.GetParent(), Direction::Rootward, clade, &clade_tree_id);
    if ((best_edge_id == NoId) || (best_tree_id > clade_tree_id)) {
      best_edge_id = clade_edge_id;
      best_tree_id = clade_tree_id;
    }
  }
  GetChoiceMap().SetEdgeChoice(edge_id, AdjacentEdge::Parent, best_edge_id);
  // Select sister.
  // std::cout << "select_sister" << std::endl;
  best_edge_id = GetBestEdgeIdByHighestPriorityTree(
      edge.GetParent(), Direction::Leafward, Bitset::Opposite(edge.GetSubsplitClade()));
  GetChoiceMap().SetEdgeChoice(edge_id, AdjacentEdge::Sister, best_edge_id);
  // Select left child.
  // std::cout << "select_leftchild" << std::endl;
  best_edge_id = GetBestEdgeIdByHighestPriorityTree(
      edge.GetChild(), Direction::Leafward, SubsplitClade::Left);
  GetChoiceMap().SetEdgeChoice(edge_id, AdjacentEdge::LeftChild, best_edge_id);
  // Select right child.
  // std::cout << "select_rightchild" << std::endl;
  best_edge_id = GetBestEdgeIdByHighestPriorityTree(
      edge.GetChild(), Direction::Leafward, SubsplitClade::Right);
  GetChoiceMap().SetEdgeChoice(edge_id, AdjacentEdge::RightChild, best_edge_id);
}

void TPEngine::UpdateEdgeChoiceByTakingHighestScoringTree(const EdgeId edge_id) {}

void TPEngine::SetTreeSourceByTakingFirst(const RootedTreeCollection &tree_collection,
                                          const BitsetSizeMap &edge_indexer) {
  input_tree_count_ = tree_collection.TreeCount();
  tree_counter_ = input_tree_count_ + 1;
  const auto tree_id_max = TreeId(input_tree_count_ + 1);
  GetTreeSource().resize(GetEdgeCount(), tree_id_max);
  std::fill(GetTreeSource().begin(), GetTreeSource().end(), tree_id_max);
  // Set tree source map for each edge in DAG.
  auto set_tree_source = [tree_id_max, this](
                             const EdgeId edge_id, const Bitset &edge_bitset,
                             const RootedTree &tree, const size_t tree_id,
                             const Node *focal_node) {
    if (GetTreeSource(edge_id) == tree_id_max) {
      GetTreeSource(edge_id) = TreeId(tree_id + 1);
    }
  };
  RootedSBNMaps::FunctionOverRootedTreeCollection(set_tree_source, tree_collection,
                                                  edge_indexer, NoId);
  // Set tree source map for rootsplit edges from the best tree.
  const auto root_node = GetDAG().GetDAGNode(GetDAG().GetDAGRootNodeId());
  for (const auto rootsplit_node_id : GetDAG().GetRootsplitNodeIds()) {
    const auto rootsplit_node = GetDAG().GetDAGNode(rootsplit_node_id);
    const auto rootsplit_edge_id =
        GetDAG().GetEdgeIdx(root_node.Id(), rootsplit_node.Id());
    TreeId best_tree_source = tree_id_max;
    GetDAG().IterateOverLeafwardEdges(
        rootsplit_node, [this, &rootsplit_node, &rootsplit_edge_id, &best_tree_source](
                            bool is_edge_on_left, SubsplitDAGNode child_node) {
          const auto edge_id =
              GetDAG().GetEdgeIdx(rootsplit_node.Id(), child_node.Id());
          if (best_tree_source > GetTreeSource(edge_id)) {
            best_tree_source = GetTreeSource()[edge_id.value_];
            GetTreeSource(rootsplit_edge_id) = best_tree_source;
          }
        });
  }
}

void TPEngine::SetChoiceMapByTakingFirst(const RootedTreeCollection &tree_collection,
                                         const BitsetSizeMap &edge_indexer,
                                         const bool use_subsplit_method) {
  // First, update tree sources over tree by taking first occurance of tree.
  SetTreeSourceByTakingFirst(tree_collection, edge_indexer);
  // Use Subsplit Heuristic.
  if (use_subsplit_method) {
    // Assign choice map by choosing adjacent edges from tree source.
    for (EdgeId edge_id = 0; edge_id < GetEdgeCount(); edge_id++) {
      UpdateEdgeChoiceByTakingHighestPriorityTree(edge_id);
    }
  }
  // Use PCSP Heuristic.
  else {
    // Build maps so we can traverse trees rootward.
    std::unordered_map<size_t, std::unordered_map<size_t, const Node *>>
        all_parent_maps;
    size_t tree_id = 0;
    for (const auto &tree : tree_collection) {
      auto parent_map = tree.Topology()->BuildParentNodeMap();
      all_parent_maps[tree_id] = parent_map;
      tree_id++;
    }

    auto SetEdgeChoice = [this](const Node *grandparent, const Node *parent,
                                const Node *sister, const Node *child,
                                const Node *grandchild0, const Node *grandchild1) {
      // Find central edge and get
      const auto parent_node_id = GetDAG().GetDAGNodeId(parent->BuildSubsplit());
      const auto child_node_id = GetDAG().GetDAGNodeId(child->BuildSubsplit());
      const auto central_edge_id = GetDAG().GetEdgeIdx(parent_node_id, child_node_id);
      auto &edge_choice = GetChoiceMap(central_edge_id);
      // Assign sister edge.
      const auto sister_node_id = GetDAG().GetDAGNodeId(sister->BuildSubsplit());
      const auto sister_edge_id = GetDAG().GetEdgeIdx(parent_node_id, sister_node_id);
      edge_choice.sister_edge_id = sister_edge_id;
      // Assign parent edge if parent node is not root.
      if (grandparent) {
        const auto grandparent_node_id =
            GetDAG().GetDAGNodeId(grandparent->BuildSubsplit());
        const auto parent_edge_id =
            GetDAG().GetEdgeIdx(grandparent_node_id, parent_node_id);
        edge_choice.parent_edge_id = parent_edge_id;
      }
      // Assign child edges if child node is not leaf.
      if (grandchild0) {
        const auto grandchild0_node_id =
            GetDAG().GetDAGNodeId(grandchild0->BuildSubsplit());
        const auto grandchild1_node_id =
            GetDAG().GetDAGNodeId(grandchild1->BuildSubsplit());
        const auto child0_edge_id =
            GetDAG().GetEdgeIdx(child_node_id, grandchild0_node_id);
        const auto child1_edge_id =
            GetDAG().GetEdgeIdx(child_node_id, grandchild1_node_id);
        const auto clade0 = GetDAG().GetFocalClade(child0_edge_id);
        edge_choice.left_child_edge_id =
            (clade0 == SubsplitClade::Left) ? child0_edge_id : child1_edge_id;
        edge_choice.right_child_edge_id =
            (clade0 == SubsplitClade::Left) ? child1_edge_id : child0_edge_id;
      }
    };

    auto FuncOnNeighboringNodes = [this, &all_parent_maps, &SetEdgeChoice](
                                      const EdgeId edge_id, const Bitset &edge_bitset,
                                      const RootedTree &tree, const size_t tree_id,
                                      const Node *node) {
      // Only set edge choices if this is the tree source assigned to edge.
      const Node *grandparent_node = nullptr;
      const Node *parent_node = nullptr;
      const Node *sister_node = nullptr;
      const Node *child_node = node;
      const Node *grandchild0_node = nullptr;
      const Node *grandchild1_node = nullptr;
      auto &parent_map = all_parent_maps[tree_id];
      // Only continue if tree source of central edge is from tree.
      if (GetTreeSource(edge_id) != tree_id + 1) {
        return;
      }
      // Find parent node. Ignore current node if child is a root node.
      if (parent_map.find(child_node->Id()) == parent_map.end()) {
        return;
      }
      parent_node = parent_map[child_node->Id()];
      // Find grandparent node if exists.
      if (parent_map.find(parent_node->Id()) != parent_map.end()) {
        grandparent_node = parent_map[parent_node->Id()];
      }
      // Find sister node.
      for (const auto &childx_node : parent_node->Children()) {
        if (childx_node->Id() != child_node->Id()) {
          sister_node = childx_node.get();
        }
      }
      // Find grandchild nodes.
      grandchild0_node = child_node->Children()[0].get();
      grandchild1_node = child_node->Children()[1].get();

      SetEdgeChoice(grandparent_node, parent_node, sister_node, child_node,
                    grandchild0_node, grandchild1_node);
    };

    RootedSBNMaps::FunctionOverRootedTreeCollection(
        FuncOnNeighboringNodes, tree_collection, edge_indexer, NoId);
  }
}

NNIOperation TPEngine::FindHighestPriorityNeighborNNIInDAG(
    const NNIOperation &nni) const {
  // Select pre-NNI by taking one with the highest priority in ChoiceMap.
  NNIOperation best_pre_nni;
  TreeId best_tree_id = TreeId(NoId);
  bool nni_found = false;
  auto pre_nnis = GetDAG().FindAllNNINeighborsInDAG(nni);
  for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
    const auto pre_nni = pre_nnis[clade];
    if (pre_nni.has_value()) {
      nni_found = true;
      const auto edge_id = GetDAG().GetEdgeIdx(pre_nni.value());
      const auto tree_id = GetTreeSource(edge_id);
      if ((best_tree_id == NoId) || (best_tree_id.value_ > tree_id.value_)) {
        best_tree_id = tree_id;
        best_pre_nni = pre_nni.value();
      }
    }
  }
  Assert(nni_found, "DAG does not contain a neighboring NNI to given NNI.");
  return best_pre_nni;
}

std::unordered_map<EdgeId, EdgeId> TPEngine::BuildAdjacentEdgeMapFromPostNNIToPreNNI(
    const NNIOperation &pre_nni, const NNIOperation &post_nni) const {
  using NNIClade = NNIOperation::NNIClade;
  using Adj = std::pair<NodeId, Direction>;
  using AdjMap = NNIOperation::NNICladeEnum::Array<Adj>;

  std::unordered_map<EdgeId, EdgeId> edge_map;
  const auto pre_edge_id = GetDAG().GetEdgeIdx(pre_nni);
  const auto &pre_edge = GetDAG().GetDAGEdge(pre_edge_id);
  const auto post_edge_id = GetDAG().GetEdgeIdx(post_nni);
  const auto &post_edge = GetDAG().GetDAGEdge(post_edge_id);
  AdjMap pre_adj_map_1, pre_adj_map_2, post_adj_map;

  // Base maps for adjacent edges.
  for (const auto &edge_id : {pre_edge_id, post_edge_id}) {
    const auto &edge = (edge_id == pre_edge_id) ? pre_edge : post_edge;
    auto &adj_map = (edge_id == pre_edge_id) ? pre_adj_map_1 : post_adj_map;
    adj_map[NNIClade::ParentFocal] = {edge.GetParent(), Direction::Rootward};
    adj_map[NNIClade::ParentSister] = {edge.GetParent(), Direction::Leafward};
    adj_map[NNIClade::ChildLeft] = {edge.GetChild(), Direction::Leafward};
    adj_map[NNIClade::ChildRight] = {edge.GetChild(), Direction::Leafward};
  }

  // Map central edge.
  edge_map[post_edge_id] = pre_edge_id;
  // Map adjacent edges.
  for (const auto post_node_dir : {Direction::Rootward, Direction::Leafward}) {
    const auto post_node_id = (post_node_dir == Direction::Rootward)
                                  ? post_edge.GetParent()
                                  : post_edge.GetChild();
    const auto &post_node = GetDAG().GetDAGNode(post_node_id);
    for (const auto pre_node_dir : {Direction::Rootward, Direction::Leafward}) {
      const auto pre_node_id = (pre_node_dir == Direction::Rootward)
                                   ? pre_edge.GetParent()
                                   : pre_edge.GetChild();
      // Iterate over all edges adjacent to post-NNI nodes.
      for (const auto adj_node_dir : {Direction::Rootward, Direction::Leafward}) {
        for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
          for (const auto adj_node_id : post_node.GetNeighbors(adj_node_dir, clade)) {
            // Get adjacent post-NNI edge id.
            const auto post_parent_node_id =
                (adj_node_dir == Direction::Rootward) ? adj_node_id : post_node_id;
            const auto post_child_node_id =
                (adj_node_dir == Direction::Rootward) ? post_node_id : adj_node_id;
            const auto adj_post_edge_id =
                GetDAG().GetEdgeIdx(post_parent_node_id, post_child_node_id);
            // Get matching adjacent pre-NNI edge id, if a match exists.
            const auto pre_parent_node_id =
                (adj_node_dir == Direction::Rootward) ? adj_node_id : pre_node_id;
            const auto pre_child_node_id =
                (adj_node_dir == Direction::Rootward) ? pre_node_id : adj_node_id;
            // If match exists, add it to the map.
            if (GetDAG().ContainsEdge(pre_parent_node_id, pre_child_node_id)) {
              const auto adj_pre_edge_id =
                  GetDAG().GetEdgeIdx(pre_parent_node_id, pre_child_node_id);
              edge_map[adj_post_edge_id] = adj_pre_edge_id;
            }
          }
        }
      }
    }
  }
  return edge_map;
}

TPChoiceMap::EdgeChoice TPEngine::RemapEdgeChoiceFromPreNNIToPostNNI(
    const TPChoiceMap::EdgeChoice &in_choice,
    const NNIOperation::NNICladeArray &clade_map) const {
  using NNIClade = NNIOperation::NNIClade;
  using NNICladeEnum = NNIOperation::NNICladeEnum;
  TPChoiceMap::EdgeChoice out_choice;
  NNICladeEnum::Array<EdgeId> pre_clade_map, post_clade_map;
  pre_clade_map[NNIClade::ParentFocal] = in_choice.parent_edge_id;
  pre_clade_map[NNIClade::ParentSister] = in_choice.sister_edge_id;
  pre_clade_map[NNIClade::ChildLeft] = in_choice.left_child_edge_id;
  pre_clade_map[NNIClade::ChildRight] = in_choice.right_child_edge_id;
  for (const auto clade : NNICladeEnum::Iterator()) {
    post_clade_map[clade] = pre_clade_map[clade_map[clade]];
  }
  out_choice.parent_edge_id = post_clade_map[NNIClade::ParentFocal];
  out_choice.sister_edge_id = post_clade_map[NNIClade::ParentSister];
  out_choice.left_child_edge_id = post_clade_map[NNIClade::ChildLeft];
  out_choice.right_child_edge_id = post_clade_map[NNIClade::ChildRight];
  return out_choice;
}

TPChoiceMap::EdgeChoice TPEngine::GetRemappedEdgeChoiceFromPreNNIToPostNNI(
    const NNIOperation &pre_nni, const NNIOperation &post_nni) const {
  TPChoiceMap::EdgeChoice mapped_choice;
  const auto pre_edge_id = GetDAG().GetEdgeIdx(pre_nni);
  const auto post_edge_id = GetDAG().GetEdgeIdx(post_nni);
  const auto &post_edge = GetDAG().GetDAGEdge(post_edge_id);

  const auto pre_choice = GetChoiceMap().GetEdgeChoice(pre_edge_id);
  const auto clade_map =
      NNIOperation::BuildNNICladeMapFromPreNNIToNNI(post_nni, pre_nni);
  const auto post_choice = RemapEdgeChoiceFromPreNNIToPostNNI(pre_choice, clade_map);
  const auto node_choice = GetChoiceMap().GetNodeIdsFromEdgeChoice(post_choice);

  // Find common nodes between pre-NNI and post-NNI.  Then use them to
  // find edges that go to the common nodes in post-NNI.
  auto GetEdgeId = [this](const NodeId parent_node_id, const NodeId child_node_id) {
    if ((parent_node_id == NoId) || (child_node_id == NoId)) {
      return EdgeId(NoId);
    }
    return GetDAG().GetEdgeIdx(parent_node_id, child_node_id);
  };
  mapped_choice.parent_edge_id =
      GetEdgeId(node_choice.parent_node_id, post_edge.GetParent());
  mapped_choice.sister_edge_id =
      GetEdgeId(post_edge.GetParent(), node_choice.sister_node_id);
  mapped_choice.left_child_edge_id =
      GetEdgeId(post_edge.GetChild(), node_choice.left_child_node_id);
  mapped_choice.right_child_edge_id =
      GetEdgeId(post_edge.GetChild(), node_choice.right_child_node_id);
  return mapped_choice;
}

double TPEngine::GetAvgLengthOfAdjEdges(
    const NodeId parent_node_id, const NodeId child_node_id,
    const std::optional<size_t> prev_node_count,
    const std::optional<Reindexer> node_reindexer,
    const std::optional<size_t> prev_edge_count,
    const std::optional<Reindexer> edge_reindexer) const {
  double total_branch_length = 0.0;
  double branch_count = 0;
  // Checks if edge was in previous version of DAG.
  auto IsEdgeOld = [this, &prev_edge_count, &edge_reindexer](const EdgeId edge_id) {
    if (!prev_edge_count.has_value()) {
      return true;
    }
    if (edge_reindexer->GetOldIndexByNewIndex(edge_id.value_) <
        prev_edge_count.value()) {
      return true;
    }
    return false;
  };

  // If an old edge exists, then use that branch length.
  if (GetDAG().ContainsEdge(parent_node_id, child_node_id)) {
    const auto edge_id = GetDAG().GetEdgeIdx(parent_node_id, child_node_id);
    if (IsEdgeOld(edge_id)) {
      return GetLikelihoodEvalEngine().GetDAGBranchHandler()(edge_id);
    }
  }
  // Iterate over alternate edges from parent and child and take average of all
  // neighboring branches.
  for (const auto node_dir : {Direction::Rootward, Direction::Leafward}) {
    const auto node_id =
        (node_dir == Direction::Rootward) ? parent_node_id : child_node_id;
    const auto other_node_id =
        (node_dir == Direction::Rootward) ? child_node_id : parent_node_id;
    if (!GetDAG().ContainsNode(node_id)) {
      continue;
    }
    const auto &node = GetDAG().GetDAGNode(node_id);
    const auto adj_node_dir =
        (node_dir == Direction::Rootward) ? Direction::Leafward : Direction::Rootward;
    for (const auto clade : {SubsplitClade::Left, SubsplitClade::Right}) {
      for (const auto adj_node_id : node.GetNeighbors(adj_node_dir, clade)) {
        if (adj_node_id == other_node_id) {
          continue;
        }
        const auto parent_id =
            (adj_node_dir == Direction::Rootward) ? adj_node_id : node_id;
        const auto child_id =
            (adj_node_dir == Direction::Rootward) ? node_id : adj_node_id;
        const auto edge_id = GetDAG().GetEdgeIdx(parent_id, child_id);
        if (IsEdgeOld(edge_id)) {
          branch_count += 1;
          total_branch_length +=
              GetLikelihoodEvalEngine().GetDAGBranchHandler()(edge_id);
        }
      }
    }
  }
  if (branch_count == 0) {
    return GetLikelihoodEvalEngine().GetDAGBranchHandler().GetDefaultBranchLength();
  }
  return total_branch_length / branch_count;
}

BitsetEdgeIdMap TPEngine::BuildBestEdgeMapOverNNIs(
    const NNISet &nnis, std::optional<const size_t> prev_edge_count,
    std::optional<const Reindexer> edge_reindexer) const {
  // Stores the best pre-NNI reference edge.
  std::unordered_map<Bitset, EdgeId> best_edge_ids;
  // Stores the best edge's tree id.
  std::unordered_map<Bitset, TreeId> best_tree_ids;
  // Checks if edge is a most recent addition.
  auto IsEdgeOld = [this, &prev_edge_count, &edge_reindexer](const EdgeId edge_id) {
    if (!prev_edge_count.has_value()) {
      return true;
    }
    if (edge_reindexer->GetOldIndexByNewIndex(edge_id.value_) <
        prev_edge_count.value()) {
      return true;
    }
    return false;
  };
  // Update map if better edge is found.
  auto AssignBestEdge = [this, &best_edge_ids, &best_tree_ids, &IsEdgeOld](
                            const Bitset &pcsp, const EdgeId proposed_edge_id) {
    // If edge was not recently added, ref itself with the highest priority,
    if (GetDAG().ContainsEdge(pcsp)) {
      const auto edge_id = GetDAG().GetEdgeIdx(pcsp);
      if (IsEdgeOld(edge_id)) {
        best_edge_ids[pcsp] = edge_id;
        best_tree_ids[pcsp] = TreeId(0);
      }
    }
    // If edge_id has not been seen yet or current reference edge has a higher
    // tree_id, then update.
    if ((best_edge_ids.find(pcsp) == best_edge_ids.end()) ||
        (best_tree_ids[pcsp] > GetTreeSource(proposed_edge_id))) {
      best_edge_ids[pcsp] = proposed_edge_id;
      best_tree_ids[pcsp] = GetTreeSource(proposed_edge_id);
    }
  };

  // Iterate over NNIs to find best branches.
  for (const auto &post_nni : nnis) {
    // Get pre_nni and build map between them.
    const NNIOperation pre_nni = FindHighestPriorityNeighborNNIInDAG(post_nni);
    const auto pre_edge_id = GetDAG().GetEdgeIdx(pre_nni);
    const auto rev_clade_map =
        NNIOperation::BuildNNICladeMapFromPreNNIToNNI(post_nni, pre_nni);
    // Get edge choice map from pre-NNI in DAG, then remap according to post-NNI.
    const auto pre_choice = GetChoiceMap(pre_edge_id);
    const auto mapped_pre_choice =
        RemapEdgeChoiceFromPreNNIToPostNNI(pre_choice, rev_clade_map);
    const auto adj_node_ids =
        GetChoiceMap().GetNodeIdsFromEdgeChoice(mapped_pre_choice);
    if ((adj_node_ids.parent_node_id == NoId) &&
        (adj_node_ids.sister_node_id == NoId) &&
        (adj_node_ids.left_child_node_id == NoId) &&
        (adj_node_ids.right_child_node_id == NoId)) {
      std::cout << "WARNING: Skipped NNI -- " << post_nni << std::endl;
      continue;
    }
    // Parent edge.
    const auto &parent_subsplit =
        GetDAG().GetDAGNodeBitset(adj_node_ids.parent_node_id);
    const auto parent_pcsp = Bitset::PCSP(parent_subsplit, post_nni.GetParent());
    AssignBestEdge(parent_pcsp, mapped_pre_choice.parent_edge_id);
    // Central edge.
    const auto central_pcsp = Bitset::PCSP(post_nni.GetParent(), post_nni.GetChild());
    AssignBestEdge(central_pcsp, pre_edge_id);
    // Sister edge.
    const auto &sister_subsplit =
        GetDAG().GetDAGNodeBitset(adj_node_ids.sister_node_id);
    const auto sister_pcsp = Bitset::PCSP(post_nni.GetParent(), sister_subsplit);
    AssignBestEdge(sister_pcsp, mapped_pre_choice.sister_edge_id);
    // LeftChild edge.
    const auto &leftchild_subsplit =
        GetDAG().GetDAGNodeBitset(adj_node_ids.left_child_node_id);
    const auto leftchild_pcsp = Bitset::PCSP(post_nni.GetChild(), leftchild_subsplit);
    AssignBestEdge(leftchild_pcsp, mapped_pre_choice.left_child_edge_id);
    // RightChild edge.
    const auto &rightchild_subsplit =
        GetDAG().GetDAGNodeBitset(adj_node_ids.right_child_node_id);
    const auto rightchild_pcsp = Bitset::PCSP(post_nni.GetChild(), rightchild_subsplit);
    AssignBestEdge(rightchild_pcsp, mapped_pre_choice.right_child_edge_id);
  }
  return best_edge_ids;
}

// ** Branch Lengths

void TPEngine::SetBranchLengthsToDefault() {
  auto &dag_branch_lengths = GetLikelihoodEvalEngine().GetDAGBranchHandler();
  for (EdgeId edge_idx = 0; edge_idx < GetEdgeCount(); edge_idx++) {
    dag_branch_lengths(edge_idx) = dag_branch_lengths.GetDefaultBranchLength();
  }
}

void TPEngine::SetBranchLengthsByTakingFirst(
    const RootedTreeCollection &tree_collection, const BitsetSizeMap &edge_indexer,
    const bool set_uninitialized_to_default) {
  auto &dag_branch_lengths = GetLikelihoodEvalEngine().GetDAGBranchHandler();
  if (set_uninitialized_to_default) {
    SetBranchLengthsToDefault();
  }
  // Unique edges in collection should be the same as the number of total edges in
  // DAG created from collection.
  EigenVectorXi observed_edge_counts = EigenVectorXi::Zero(GetEdgeCount());
  // Set branch lengths on first occurance.
  auto set_first_branch_length = [this, &observed_edge_counts, &dag_branch_lengths](
                                     const EdgeId edge_idx, const Bitset &edge_bitset,
                                     const RootedTree &tree, const size_t tree_id,
                                     const Node *focal_node) {
    if (observed_edge_counts(edge_idx.value_) == 0) {
      dag_branch_lengths(edge_idx) = tree.BranchLength(focal_node);
      observed_edge_counts(edge_idx.value_)++;
    }
  };
  RootedSBNMaps::FunctionOverRootedTreeCollection(
      set_first_branch_length, tree_collection, edge_indexer,
      dag_branch_lengths.GetBranchLengthData().size());
}

void TPEngine::OptimizeBranchLengths(std::optional<bool> check_branch_convergence) {
  GetLikelihoodEvalEngine().BranchLengthOptimization(check_branch_convergence);
}

// ** TP Evaluation Engine

void TPEngine::MakeLikelihoodEvalEngine(const std::string &mmap_likelihood_path) {
  likelihood_engine_ =
      std::make_unique<TPEvalEngineViaLikelihood>(*this, mmap_likelihood_path);
  eval_engine_ = likelihood_engine_.get();
}

void TPEngine::MakeParsimonyEvalEngine(const std::string &mmap_parsimony_path) {
  parsimony_engine_ =
      std::make_unique<TPEvalEngineViaParsimony>(*this, mmap_parsimony_path);
  eval_engine_ = parsimony_engine_.get();
}

void TPEngine::ClearEvalEngineInUse() {
  for (auto eval_engine_type : TPEvalEngineTypeEnum::Iterator()) {
    eval_engine_in_use_[eval_engine_type] = false;
  }
}

void TPEngine::SelectEvalEngine(const TPEvalEngineType eval_engine_type) {
  switch (eval_engine_type) {
    case TPEvalEngineType::LikelihoodEvalEngine:
      SelectLikelihoodEvalEngine();
      break;
    case TPEvalEngineType::ParsimonyEvalEngine:
      SelectParsimonyEvalEngine();
      break;
    default:
      Failwith("Invalid TPEvalEngineType.");
  }
}

void TPEngine::SelectLikelihoodEvalEngine() {
  Assert(HasLikelihoodEvalEngine(),
         "Must MakeLikelihoodEvalEngine before selecting it.");
  ClearEvalEngineInUse();
  eval_engine_in_use_[TPEvalEngineType::LikelihoodEvalEngine] = true;
  eval_engine_ = likelihood_engine_.get();
}

void TPEngine::SelectParsimonyEvalEngine() {
  Assert(HasParsimonyEvalEngine(), "Must MakeParimonyEvalEngine before selecting it.");
  ClearEvalEngineInUse();
  eval_engine_in_use_[TPEvalEngineType::ParsimonyEvalEngine] = true;
  eval_engine_ = parsimony_engine_.get();
}

void TPEngine::UpdateEvalEngineAfterModifyingDAG(
    const std::map<NNIOperation, NNIOperation> &nni_to_pre_nni,
    const size_t prev_node_count, const Reindexer &node_reindexer,
    const size_t prev_edge_count, const Reindexer &edge_reindexer) {
  UpdateChoiceMapAfterModifyingDAG(nni_to_pre_nni, prev_node_count, node_reindexer,
                                   prev_edge_count, edge_reindexer);
  if (HasLikelihoodEvalEngine()) {
    GetLikelihoodEvalEngine().UpdateEngineAfterModifyingDAG(
        nni_to_pre_nni, prev_node_count, node_reindexer, prev_edge_count,
        edge_reindexer);
  }
  if (HasParsimonyEvalEngine()) {
    GetParsimonyEvalEngine().UpdateEngineAfterModifyingDAG(
        nni_to_pre_nni, prev_node_count, node_reindexer, prev_edge_count,
        edge_reindexer);
  }
}

// ** Scoring

double TPEngine::GetTopTreeScoreWithProposedNNI(const NNIOperation &post_nni,
                                                const NNIOperation &pre_nni,
                                                const size_t spare_offset) {
  return GetEvalEngine().GetTopTreeScoreWithProposedNNI(post_nni, pre_nni,
                                                        spare_offset);
}

// ** Maintenance

void TPEngine::InitializeScores() {
  if (HasLikelihoodEvalEngine()) {
    GetLikelihoodEvalEngine().Initialize();
  }
  if (HasParsimonyEvalEngine()) {
    GetParsimonyEvalEngine().Initialize();
  }
}

void TPEngine::ComputeScores() {
  if (HasLikelihoodEvalEngine()) {
    GetLikelihoodEvalEngine().ComputeScores();
  }
  if (HasParsimonyEvalEngine()) {
    GetParsimonyEvalEngine().ComputeScores();
  }
}

void TPEngine::UpdateScoresAfterDAGAddNodePair(const NNIOperation &post_nni,
                                               const NNIOperation &pre_nni,
                                               std::optional<size_t> new_tree_id) {
  if (HasLikelihoodEvalEngine()) {
    GetLikelihoodEvalEngine().UpdateEngineAfterDAGAddNodePair(post_nni, pre_nni,
                                                              new_tree_id);
  }
  if (HasParsimonyEvalEngine()) {
    GetParsimonyEvalEngine().UpdateEngineAfterDAGAddNodePair(post_nni, pre_nni,
                                                             new_tree_id);
  }
}

EdgeId TPEngine::FindHighestPriorityEdgeAdjacentToNode(
    const NodeId node_id, const Direction direction) const {
  const auto node = GetDAG().GetDAGNode(node_id);
  EdgeId best_edge_id = EdgeId(NoId);
  auto best_tree_source = GetNextTreeId();
  for (const auto focal_clade : {SubsplitClade::Left, SubsplitClade::Right}) {
    for (const auto rootward_node_id : node.GetNeighbors(direction, focal_clade)) {
      const auto edge_id = GetDAG().GetEdgeIdx(rootward_node_id, node_id);
      if (best_tree_source > GetTreeSource(edge_id)) {
        best_tree_source = GetTreeSource(edge_id);
        best_edge_id = edge_id;
      }
    }
  }
  return best_edge_id;
}

EdgeId TPEngine::FindHighestPriorityEdgeAdjacentToNode(
    const NodeId node_id, const Direction direction, const SubsplitClade clade) const {
  const auto node = GetDAG().GetDAGNode(node_id);
  EdgeId best_edge_id = EdgeId(NoId);
  auto best_tree_source = GetNextTreeId();
  for (const auto rootward_node_id : node.GetNeighbors(direction, clade)) {
    const auto edge_id = GetDAG().GetEdgeIdx(rootward_node_id, node_id);
    if (best_tree_source > GetTreeSource(edge_id)) {
      best_tree_source = GetTreeSource(edge_id);
      best_edge_id = edge_id;
    }
  }
  return best_edge_id;
}

void TPEngine::CopyOverEdgeDataFromPreNNIToPostNNI(const NNIOperation &post_nni,
                                                   const NNIOperation &pre_nni,
                                                   CopyEdgeDataFunc copy_data_func,
                                                   std::optional<size_t> new_tree_id) {
  new_tree_id =
      (new_tree_id.has_value()) ? new_tree_id.value() : GetNextTreeId().value_;
  input_tree_count_ = new_tree_id.value() + 1;
  // Copy over all adjacent branch lengths from those
  auto CopyBranchLengthFromCommonAdjacentNodes =
      [this, &copy_data_func, new_tree_id](
          const NodeId pre_node_id, const NodeId post_node_id,
          const Direction direction, const SubsplitClade clade) {
        const auto pre_node = GetDAG().GetDAGNode(pre_node_id);
        for (const auto parent_id : pre_node.GetNeighbors(direction, clade)) {
          const auto pre_edge_id = GetDAG().GetEdgeIdx(parent_id, pre_node_id);
          const auto post_edge_id = GetDAG().GetEdgeIdx(parent_id, post_node_id);
          copy_data_func(pre_edge_id, post_edge_id);
          // GetTreeSource(post_edge_id) = TreeId(new_tree_id.value());
        }
      };
  const auto pre_parent_id = GetDAG().GetDAGNodeId(pre_nni.GetParent());
  const auto pre_child_id = GetDAG().GetDAGNodeId(pre_nni.GetChild());
  const auto pre_edge_id = GetDAG().GetEdgeIdx(pre_parent_id, pre_child_id);
  const auto pre_edge = GetDAG().GetDAGEdge(pre_edge_id);
  const auto post_parent_id = GetDAG().GetDAGNodeId(post_nni.GetParent());
  const auto post_child_id = GetDAG().GetDAGNodeId(post_nni.GetChild());
  const auto post_edge_id = GetDAG().GetEdgeIdx(post_parent_id, post_child_id);
  // Copy over central edge.
  copy_data_func(pre_edge_id, post_edge_id);
  // Copy over parent and sister edges.
  CopyBranchLengthFromCommonAdjacentNodes(pre_parent_id, post_parent_id,
                                          Direction::Rootward, SubsplitClade::Left);
  CopyBranchLengthFromCommonAdjacentNodes(pre_parent_id, post_parent_id,
                                          Direction::Rootward, SubsplitClade::Right);
  CopyBranchLengthFromCommonAdjacentNodes(
      pre_parent_id, post_child_id, Direction::Leafward,
      Bitset::Opposite(pre_edge.GetSubsplitClade()));
  // Copy over left and right edges.
  NodeId post_leftchild_id;
  NodeId post_rightchild_id;
  if (pre_nni.GetSisterClade() == post_nni.GetLeftChildClade()) {
    // If post_nni swapped pre_nni sister with pre_nni left child.
    post_leftchild_id = post_parent_id;
    post_rightchild_id = post_child_id;
  } else {
    // If post_nni swapped pre_nni sister swapped with pre_nni right child.
    post_leftchild_id = post_child_id;
    post_rightchild_id = post_parent_id;
  }
  CopyBranchLengthFromCommonAdjacentNodes(pre_child_id, post_leftchild_id,
                                          Direction::Leafward, SubsplitClade::Left);
  CopyBranchLengthFromCommonAdjacentNodes(pre_child_id, post_rightchild_id,
                                          Direction::Leafward, SubsplitClade::Right);
}

// ** Parameter Data

void TPEngine::GrowNodeData(const size_t new_node_count,
                            std::optional<const Reindexer> node_reindexer,
                            std::optional<const size_t> explicit_alloc,
                            const bool on_init) {
  if (HasLikelihoodEvalEngine()) {
    GetLikelihoodEvalEngine().GrowNodeData(new_node_count, node_reindexer,
                                           explicit_alloc, on_init);
  }
  if (HasParsimonyEvalEngine()) {
    GetParsimonyEvalEngine().GrowNodeData(new_node_count, node_reindexer,
                                          explicit_alloc, on_init);
  }
  const size_t old_node_count = GetNodeCount();
  SetNodeCount(new_node_count);
  // Reallocate more space if needed.
  if ((GetPaddedNodeCount() > GetAllocatedNodeCount()) || explicit_alloc.has_value()) {
    SetAllocatedNodeCount(
        size_t(ceil(double(GetPaddedNodeCount()) * resizing_factor_)));
    if (explicit_alloc.has_value()) {
      Assert(explicit_alloc.value() >= GetNodeCount(),
             "Attempted to reallocate space smaller than node_count.");
      SetAllocatedNodeCount(explicit_alloc.value() + GetSpareNodeCount());
    }
  }
  // Reindex work space to realign with DAG.
  if (node_reindexer.has_value()) {
    ReindexNodeData(node_reindexer.value(), old_node_count);
  }
}

void TPEngine::GrowEdgeData(const size_t new_edge_count,
                            std::optional<const Reindexer> edge_reindexer,
                            std::optional<const size_t> explicit_alloc,
                            const bool on_init) {
  GetChoiceMap().GrowEdgeData(new_edge_count, edge_reindexer, explicit_alloc, on_init);
  if (HasLikelihoodEvalEngine()) {
    GetLikelihoodEvalEngine().GrowEdgeData(new_edge_count, edge_reindexer,
                                           explicit_alloc, on_init);
  }
  if (HasParsimonyEvalEngine()) {
    GetParsimonyEvalEngine().GrowEdgeData(new_edge_count, edge_reindexer,
                                          explicit_alloc, on_init);
  }
  const size_t old_edge_count = GetEdgeCount();
  SetEdgeCount(new_edge_count);
  // Reallocate more space if needed.
  if ((GetPaddedEdgeCount() > GetAllocatedEdgeCount()) || explicit_alloc.has_value()) {
    SetAllocatedEdgeCount(
        size_t(ceil(double(GetPaddedEdgeCount()) * resizing_factor_)));
    if (explicit_alloc.has_value()) {
      Assert(explicit_alloc.value() >= GetNodeCount(),
             "Attempted to reallocate space smaller than node_count.");
      SetAllocatedEdgeCount(explicit_alloc.value() + GetSpareEdgeCount());
    }
    GetTreeSource().reserve(GetAllocatedEdgeCount());
  }

  GetTreeSource().resize(GetPaddedEdgeCount());
  // Initialize data.
  tree_counter_++;
  for (EdgeId i(old_edge_count); i < new_edge_count; i++) {
    GetTreeSource(i) = TreeId(NoId);
  }
  // Reindex work space to realign with DAG.
  if (edge_reindexer.has_value()) {
    ReindexEdgeData(edge_reindexer.value(), old_edge_count);
  }
}

void TPEngine::ReindexNodeData(const Reindexer &node_reindexer,
                               const size_t old_node_count) {
  Assert(node_reindexer.size() == GetNodeCount(),
         "Node Reindexer is the wrong size for TPEngine.");
  Assert(node_reindexer.IsValid(GetNodeCount()), "Node Reindexer is not valid.");
}

void TPEngine::ReindexEdgeData(const Reindexer &edge_reindexer,
                               const size_t old_edge_count) {
  Assert(edge_reindexer.size() == GetEdgeCount(),
         "Edge Reindexer is the wrong size for TPEngine.");
  Assert(edge_reindexer.IsValid(GetEdgeCount()),
         "Edge Reindexer is not valid for TPEngine size.");
  GetTreeSource() = Reindexer::Reindex(GetTreeSource(), edge_reindexer, GetEdgeCount());
}

void TPEngine::GrowSpareNodeData(const size_t new_node_spare_count) {
  if (new_node_spare_count > GetSpareNodeCount()) {
    SetSpareNodeCount(new_node_spare_count);
    GrowNodeData(GetNodeCount());
  }
}

void TPEngine::GrowSpareEdgeData(const size_t new_edge_spare_count) {
  if (new_edge_spare_count > GetSpareEdgeCount()) {
    SetSpareEdgeCount(new_edge_spare_count);
    GrowEdgeData(GetEdgeCount());
  }
}

// ** I/O

std::string TPEngine::LikelihoodPVToString(const PVId pv_id) const {
  return GetLikelihoodEvalEngine().GetPVs().ToString(pv_id);
}

std::string TPEngine::LogLikelihoodMatrixToString() const {
  std::stringstream out;
  auto &log_likelihoods = GetLikelihoodEvalEngine().GetMatrix();
  for (Eigen::Index i = 0; i < log_likelihoods.rows(); i++) {
    for (Eigen::Index j = 0; j < log_likelihoods.cols(); j++) {
      out << "[" << i << "," << j << "]: " << log_likelihoods(i, j) << "\t";
    }
    out << std::endl;
  }
  return out.str();
}

std::string TPEngine::ParsimonyPVToString(const PVId pv_id) const {
  return GetParsimonyEvalEngine().GetPVs().ToString(pv_id);
}

std::string TPEngine::TreeSourceToString() const {
  std::stringstream out;
  out << "TreeSource: { ";
  for (EdgeId i(0); i < GetEdgeCount(); i++) {
    out << "[Edge" << i << "]: Tree" << GetTreeSource(i) << ", ";
  }
  out << " }";
  return out.str();
}
