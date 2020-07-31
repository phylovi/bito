// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_dag.hpp"

#include <iostream>

#include "numerical_utils.hpp"

using namespace GPOperations;

GPDAG::GPDAG() : taxon_count_(0), rootsplit_and_pcsp_count_(0) {}

GPDAG::GPDAG(const RootedTreeCollection &tree_collection) :
tree_collection_(tree_collection)
{
  ProcessTrees(tree_collection_);
  BuildNodes();
  BuildEdges();
  SetPCSPIndexerEncodingToFullSubsplits();
  ExpandPCSPIndexerAndSubsplitToRange();
  CountTrees();
}

void GPDAG::CountTrees() {
  tree_count_below_ = EigenVectorXd::Ones(NodeCount());
  for (const auto &node_id : RootwardPassTraversal()) {
    const auto &node = GetDagNode(node_id);
    if (!node->IsLeaf()) {
      for (const bool rotated : {false, true}) {
        double per_rotated_count = 0.;
        // Sum options across the possible children.
        for (const auto &child_id : node->GetLeafward(rotated)) {
          per_rotated_count += tree_count_below_[child_id];
        }
        // Take the product across the number of options for the left and right branches
        // of the tree.
        tree_count_below_[node_id] *= per_rotated_count;
      }
    }
  }
  tree_count_ = 0;
  IterateOverRootsplitIds([this](size_t root_id) {
    tree_count_ += tree_count_below_[root_id];
  });
}

double GPDAG::TreeCount() const { return tree_count_; }

size_t GPDAG::NodeCount() const { return dag_nodes_.size(); }

size_t GPDAG::RootsplitAndPCSPCount() const { return rootsplit_and_pcsp_count_; }

size_t GPDAG::GeneralizedPCSPCount() const {
  size_t fake_subsplit_parameter_count = 0;
  for (size_t taxon_idx = 0; taxon_idx < taxon_count_; taxon_idx++) {
    fake_subsplit_parameter_count += dag_nodes_[taxon_idx]->GetRootwardRotated().size();
    fake_subsplit_parameter_count += dag_nodes_[taxon_idx]->GetRootwardSorted().size();
  }
  return RootsplitAndPCSPCount() + fake_subsplit_parameter_count;
}

void GPDAG::Print() const {
  for (size_t i = 0; i < dag_nodes_.size(); i++) {
    std::cout << dag_nodes_[i]->ToString() << std::endl;
  }
}

void GPDAG::PrintGPCSPIndexer() const {
  for (const auto &[pcsp, idx] : gpcsp_indexer_) {
    std::cout << pcsp.SubsplitToString() << ", " << idx << std::endl;
  }
}

GPDAGNode *GPDAG::GetDagNode(const size_t node_id) const {
  return dag_nodes_.at(node_id).get();
}

size_t GPDAG::GetPLVIndexStatic(PLVType plv_type, size_t node_count, size_t src_idx) {
  switch (plv_type) {
    case PLVType::P:
      return src_idx;
    case PLVType::P_HAT:
      return node_count + src_idx;
    case PLVType::P_HAT_TILDE:
      return 2 * node_count + src_idx;
    case PLVType::R_HAT:
      return 3 * node_count + src_idx;
    case PLVType::R:
      return 4 * node_count + src_idx;
    case PLVType::R_TILDE:
      return 5 * node_count + src_idx;
    default:
      Failwith("Invalid PLV index requested.");
  }
}

size_t GPDAG::GetPLVIndex(PLVType plv_type, size_t src_idx) const {
  return GetPLVIndexStatic(plv_type, dag_nodes_.size(), src_idx);
}

EigenVectorXd GPDAG::BuildUniformQ() const {
  EigenVectorXd q = EigenVectorXd::Ones(GeneralizedPCSPCount());
  q.segment(0, rootsplits_.size()).array() = 1. / rootsplits_.size();
  for (const auto &[_, range] : subsplit_to_range_) {
    const auto num_child_subsplits = range.second - range.first;
    double val = 1. / num_child_subsplits;
    q.segment(range.first, num_child_subsplits).array() = val;
  }
  return q;
}

Node::NodePtrVec GPDAG::GenerateAllTrees() const {
  std::vector<Node::NodePtrVec> nodes_below(NodeCount());
    
  auto GetSubnodes = [&nodes_below](const GPDAGNode *node,
                                    Node::NodePtrVec &rotated_subtrees,
                                    Node::NodePtrVec &ordered_subtrees) {
    for (const bool rotated : {false, true}) {
      for (const auto &child_id : node->GetLeafward(rotated)) {
        for (const auto &sub_tree : nodes_below.at(child_id)) {
          rotated ? rotated_subtrees.push_back(sub_tree) :
          ordered_subtrees.push_back(sub_tree);
        }
      }
    }
  };
  
  auto MergeNodes = [](Node::NodePtrVec &trees,
                       Node::NodePtrVec &rotated_subtrees,
                       Node::NodePtrVec &ordered_subtrees) {
    for (auto &rotated_subtree : rotated_subtrees) {
      for (auto &ordered_subtree : ordered_subtrees) {
        Node::NodePtr new_tree = Node::Join(rotated_subtree, ordered_subtree);
        trees.push_back(new_tree);
      }
    }
  };

  
  for (const auto &node_id : RootwardPassTraversal()) {
    const auto &node = GetDagNode(node_id);

    if (node->IsLeaf()) {
      nodes_below.at(node_id).push_back(Node::Leaf(node_id));
    } else {
      Node::NodePtrVec rotated_nodes, ordered_nodes;
      GetSubnodes(node, rotated_nodes, ordered_nodes);
      MergeNodes(nodes_below.at(node_id), rotated_nodes, ordered_nodes);
    }
  }

  Node::NodePtrVec nodes;
  IterateOverRootsplitIds([this, &nodes, &GetSubnodes, &MergeNodes](size_t root_id) {
    Node::NodePtrVec rotated_nodes, ordered_nodes;
    GetSubnodes(GetDagNode(root_id), rotated_nodes, ordered_nodes);
    MergeNodes(nodes, rotated_nodes, ordered_nodes);
  });
  
  for (auto &node : nodes) {
    node->Polish();
    std::cout << node->Newick() << std::endl;
  }
  
  return nodes;
}

EigenVectorXd GPDAG::BuildUniformPrior() const {
  EigenVectorXd q = EigenVectorXd::Ones(GeneralizedPCSPCount());

  for (const auto &node_id : RootwardPassTraversal()) {
    const auto &node = GetDagNode(node_id);
    if (!node->IsLeaf()) {
      for (const bool rotated : {false, true}) {
        double per_rotated_count = 0.;
        for (const auto &child_id : node->GetLeafward(rotated)) {
          per_rotated_count += tree_count_below_[child_id];
        }
        for (const auto &child_id : node->GetLeafward(rotated)) {
          Bitset gpcsp = node->GetBitset(rotated) + GetDagNode(child_id)->GetBitset();
          size_t gpcsp_idx = gpcsp_indexer_.at(gpcsp);
          q[gpcsp_idx] = tree_count_below_[child_id] / per_rotated_count;
        }
      }
    }
  }

  IterateOverRootsplitIds([this, &q](size_t root_id) {
    auto node = GetDagNode(root_id);
    auto gpcsp_idx = gpcsp_indexer_.at(node->GetBitset());
    q[gpcsp_idx] = tree_count_below_[root_id] / tree_count_;
  });

  return q;
}

GPOperationVector GPDAG::BranchLengthOptimization() const {
  GPOperationVector operations;
  std::unordered_set<size_t> visited_nodes;
  IterateOverRootsplitIds([this, &operations, &visited_nodes](size_t rootsplit_id) {
    ScheduleBranchLengthOptimization(rootsplit_id, visited_nodes, operations);
  });
  return operations;
}

GPOperationVector GPDAG::ComputeLikelihoods() const {
  GPOperationVector operations;
  IterateOverRealNodes([this, &operations](const GPDAGNode *node) {
    IterateOverLeafwardEdges(
        node,
        [this, node, &operations](const bool rotated, const GPDAGNode *child_node) {
          const auto gpcsp_idx =
              gpcsp_indexer_.at(node->GetBitset(rotated) + child_node->GetBitset());
          operations.push_back(Likelihood{gpcsp_idx,
                                          GetPLVIndex(RPLVType(rotated), node->Id()),
                                          GetPLVIndex(PLVType::P, child_node->Id())});
        });
  });

  const auto marginal_likelihood_operations = MarginalLikelihood();
  operations.insert(operations.end(), marginal_likelihood_operations.begin(),
                    marginal_likelihood_operations.end());
  return operations;
}

GPOperationVector GPDAG::LeafwardPass() const {
  return LeafwardPass(LeafwardPassTraversal());
}

GPOperationVector GPDAG::MarginalLikelihood() const {
  GPOperationVector operations;
  for (size_t rootsplit_idx = 0; rootsplit_idx < rootsplits_.size(); rootsplit_idx++) {
    const auto rootsplit = rootsplits_[rootsplit_idx];
    const auto root_subsplit = rootsplit + ~rootsplit;
    size_t root_idx = subsplit_to_id_.at(root_subsplit);
    operations.push_back(GPOperations::IncrementMarginalLikelihood{
        GetPLVIndex(PLVType::R_HAT, root_idx), rootsplit_idx,
        GetPLVIndex(PLVType::P, root_idx)});
  }
  return operations;
}

GPOperationVector GPDAG::RootwardPass() const {
  return RootwardPass(RootwardPassTraversal());
}

GPOperationVector GPDAG::OptimizeSBNParameters() const {
  GPOperationVector operations;
  std::unordered_set<size_t> visited_nodes;
  for (size_t &id : LeafwardPassTraversal()) {
    const auto node = GetDagNode(id);
    OptimizeSBNParametersForASubsplit(node->GetBitset(), operations);
    OptimizeSBNParametersForASubsplit(node->GetBitset().RotateSubsplit(), operations);
  }
  operations.push_back(UpdateSBNProbabilities{0, rootsplits_.size()});
  return operations;
}

GPOperationVector GPDAG::SetLeafwardZero() const {
  GPOperationVector operations;
  const auto node_count = dag_nodes_.size();
  for (size_t i = 0; i < node_count; i++) {
    operations.push_back(Zero{GetPLVIndex(PLVType::R_HAT, i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::R, i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::R_TILDE, i)});
  }
  return operations;
}

GPOperationVector GPDAG::SetRhatToStationary() const {
  GPOperationVector operations;
  IterateOverRootsplitIds([this, &operations](size_t rootsplit_id) {
    operations.push_back(
        SetToStationaryDistribution{GetPLVIndex(PLVType::R_HAT, rootsplit_id)});
  });
  return operations;
}

GPOperationVector GPDAG::SetRootwardZero() const {
  GPOperationVector operations;
  const auto node_count = dag_nodes_.size();
  for (size_t i = taxon_count_; i < node_count; i++) {
    operations.push_back(Zero{GetPLVIndex(PLVType::P, i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::P_HAT, i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::P_HAT_TILDE, i)});
  }
  return operations;
}

void GPDAG::IterateOverRealNodes(NodeLambda f) const {
  Assert(taxon_count_ < dag_nodes_.size(), "No real DAG nodes!");
  for (auto it = dag_nodes_.cbegin() + taxon_count_; it < dag_nodes_.cend(); it++) {
    f((*it).get());
  }
}

void GPDAG::IterateOverLeafwardEdges(const GPDAGNode *node,
                                     EdgeDestinationLambda f) const {
  for (bool rotated : {false, true}) {
    for (const size_t child_idx : node->GetLeafward(rotated)) {
      f(rotated, GetDagNode(child_idx));
    }
  }
}

void GPDAG::IterateOverRootwardEdges(const GPDAGNode *node,
                                     EdgeDestinationLambda f) const {
  for (bool rotated : {false, true}) {
    for (const size_t parent_idx : node->GetRootward(rotated)) {
      f(rotated, GetDagNode(parent_idx));
    }
  }
}

void GPDAG::IterateOverRootsplitIds(std::function<void(size_t)> f) const {
  for (const auto &rootsplit : rootsplits_) {
    const auto subsplit = rootsplit + ~rootsplit;
    f(subsplit_to_id_.at(rootsplit + ~rootsplit));
  }
}

std::vector<Bitset> GPDAG::GetChildrenSubsplits(const Bitset &subsplit,
                                                bool include_fake_subsplits) {
  std::vector<Bitset> children_subsplits;
  if (subsplit_to_range_.count(subsplit)) {
    const auto [start, stop] = subsplit_to_range_.at(subsplit);
    for (auto idx = start; idx < stop; idx++) {
      const auto child_subsplit = index_to_child_.at(idx);
      children_subsplits.push_back(child_subsplit);
    }
  } else if (include_fake_subsplits) {
    // In the case where second chunk of the subsplit is just a single taxon,
    // the subsplit will not map to any value (parent_to_range_[subsplit] doesn't
    // exist). But we still need to create and connect to fake subsplits in the DAG. A
    // subsplit has a fake subsplit as a child if the first chunk is non-zero and the
    // second chunk has exactly one bit set to 1.
    if (subsplit.SplitChunk(0).Any() &&
        subsplit.SplitChunk(1).SingletonOption().has_value()) {
      // The fake subsplit corresponds to the second chunk of subsplit.
      // Prepend it by 0's.
      Bitset zero(subsplit.size() / 2);
      Bitset fake_subsplit = zero + subsplit.SplitChunk(1);
      children_subsplits.push_back(fake_subsplit);
    }
  }

  return children_subsplits;
}

void GPDAG::ProcessTrees(const RootedTreeCollection &tree_collection) {
  taxon_count_ = tree_collection.TaxonCount();
  const auto topology_counter = tree_collection.TopologyCounter();

  std::tie(rootsplits_, gpcsp_indexer_, index_to_child_, subsplit_to_range_,
           rootsplit_and_pcsp_count_) =
      SBNMaps::BuildIndexerBundle(RootedSBNMaps::RootsplitCounterOf(topology_counter),
                                  RootedSBNMaps::PCSPCounterOf(topology_counter));
}

void GPDAG::CreateAndInsertNode(const Bitset &subsplit) {
  size_t id = dag_nodes_.size();
  SafeInsert(subsplit_to_id_, subsplit, id);
  dag_nodes_.push_back(std::make_unique<GPDAGNode>(id, subsplit));
}

void GPDAG::ConnectNodes(size_t idx, bool rotated) {
  const auto node = dag_nodes_[idx].get();
  // Retrieve children subsplits, set edge relation.
  const Bitset subsplit = node->GetBitset(rotated);
  const auto children = GetChildrenSubsplits(subsplit, true);
  for (const auto &child_subsplit : children) {
    const auto child_node = GetDagNode(subsplit_to_id_.at(child_subsplit));
    if (rotated) {
      node->AddLeafwardRotated(child_node->Id());
      child_node->AddRootwardRotated(node->Id());
    } else {
      node->AddLeafwardSorted(child_node->Id());
      child_node->AddRootwardSorted(node->Id());
    }
  }
}

void GPDAG::BuildNodesDepthFirst(const Bitset &subsplit,
                                 std::unordered_set<Bitset> &visited_subsplits) {
  visited_subsplits.insert(subsplit);
  for (bool rotated : {false, true}) {
    for (const auto &child_subsplit :
         GetChildrenSubsplits(PerhapsRotateSubsplit(subsplit, rotated), false)) {
      if (!visited_subsplits.count(child_subsplit)) {
        BuildNodesDepthFirst(child_subsplit, visited_subsplits);
      }
    }
  }
  CreateAndInsertNode(subsplit);
}

void GPDAG::BuildNodes() {
  std::unordered_set<Bitset> visited_subsplits;

  // We will create fake subsplits and insert to dag_nodes_.
  // These nodes will take IDs in [0, taxon_count_).
  Bitset zero(taxon_count_);
  for (size_t taxon_idx = 0; taxon_idx < taxon_count_; taxon_idx++) {
    Bitset fake(taxon_count_);
    fake.set(taxon_idx);
    Bitset fake_subsplit = zero + fake;
    CreateAndInsertNode(fake_subsplit);
  }
  // We are going to add the remaining nodes.
  // The root splits will take on the higher IDs compared to the non-rootsplits.
  for (const auto &rootsplit : rootsplits_) {
    const auto subsplit = rootsplit + ~rootsplit;
    BuildNodesDepthFirst(subsplit, visited_subsplits);
  }
}

void GPDAG::BuildEdges() {
  for (size_t i = taxon_count_; i < dag_nodes_.size(); i++) {
    ConnectNodes(i, false);
    ConnectNodes(i, true);
  }
}

void GPDAG::SetPCSPIndexerEncodingToFullSubsplits() {
  BitsetSizeMap indexer;
  for (auto it = gpcsp_indexer_.cbegin(); it != gpcsp_indexer_.cend(); ++it) {
    if (it->first.size() == taxon_count_) {
      SafeInsert(indexer, it->first + ~it->first, it->second);
    } else {
      SafeInsert(indexer, it->first.PCSPParent() + it->first.PCSPChildSubsplit(),
                 it->second);
    }
  }
  gpcsp_indexer_ = indexer;
}

void GPDAG::ExpandPCSPIndexerAndSubsplitToRange() {
  // Add fake subsplits to expand gpcsp_indexer_ and subsplit_to_range_.
  for (size_t i = 0; i < taxon_count_; i++) {
    const auto current_bitset = dag_nodes_[i]->GetBitset();
    IterateOverRootwardEdges(
        GetDagNode(i),
        [this, current_bitset](const bool rotated, const GPDAGNode *node) {
          Bitset gpcsp = node->GetBitset(rotated) + current_bitset;
          SafeInsert(subsplit_to_range_, node->GetBitset(rotated),
                     {gpcsp_indexer_.size(), gpcsp_indexer_.size() + 1});
          SafeInsert(gpcsp_indexer_, gpcsp, gpcsp_indexer_.size());
        });
  }
}

void RootwardDepthFirst(size_t id,
                        const std::vector<std::unique_ptr<GPDAGNode>> &dag_nodes,
                        SizeVector &visit_order,
                        std::unordered_set<size_t> &visited_nodes) {
  SafeInsert(visited_nodes, id);
  for (size_t child_id : dag_nodes.at(id)->GetRootwardSorted()) {
    if (!visited_nodes.count(child_id)) {
      RootwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  for (size_t child_id : dag_nodes.at(id)->GetRootwardRotated()) {
    if (!visited_nodes.count(child_id)) {
      RootwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  visit_order.push_back(id);
}

void LeafwardDepthFirst(size_t id,
                        const std::vector<std::unique_ptr<GPDAGNode>> &dag_nodes,
                        SizeVector &visit_order,
                        std::unordered_set<size_t> &visited_nodes) {
  SafeInsert(visited_nodes, id);
  for (size_t child_id : dag_nodes.at(id)->GetLeafwardSorted()) {
    if (!visited_nodes.count(child_id)) {
      LeafwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  for (size_t child_id : dag_nodes.at(id)->GetLeafwardRotated()) {
    if (!visited_nodes.count(child_id)) {
      LeafwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  visit_order.push_back(id);
}

GPOperationVector GPDAG::LeafwardPass(SizeVector visit_order) const {
  GPOperationVector operations;
  for (const size_t node_id : visit_order) {
    const auto node = GetDagNode(node_id);

    // Build rhat(s) via rhat(s) += \sum_t q(s|t) P'(s|t) r(t)
    AddRhatOperations(node, operations);
    // Multiply to get r(s) = rhat(s) \circ phat(s_tilde).
    operations.push_back(Multiply{GetPLVIndex(PLVType::R, node_id),
                                  GetPLVIndex(PLVType::R_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
    // Multiply to get r(s_tilde) = rhat(s) \circ phat(s).
    operations.push_back(Multiply{GetPLVIndex(PLVType::R_TILDE, node_id),
                                  GetPLVIndex(PLVType::R_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT, node_id)});
  }
  return operations;
}

GPOperationVector GPDAG::RootwardPass(SizeVector visit_order) const {
  GPOperationVector operations;
  for (const size_t node_id : visit_order) {
    const auto node = GetDagNode(node_id);
    if (node->IsLeaf()) {
      continue;
    }
    // Build phat(s).
    AddPhatOperations(node, false, operations);
    // Build phat(s_tilde).
    AddPhatOperations(node, true, operations);
    // Multiply to get p(s) = phat(s) \circ phat(s_tilde).
    operations.push_back(Multiply{node_id, GetPLVIndex(PLVType::P_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
  }
  return operations;
}

SizeVector GPDAG::LeafwardPassTraversal() const {
  SizeVector visit_order;
  std::unordered_set<size_t> visited_nodes;
  for (size_t leaf_id = 0; leaf_id < taxon_count_; leaf_id++) {
    RootwardDepthFirst(leaf_id, dag_nodes_, visit_order, visited_nodes);
  }
  return visit_order;
}

SizeVector GPDAG::RootwardPassTraversal() const {
  SizeVector visit_order;
  std::unordered_set<size_t> visited_nodes;
  IterateOverRootsplitIds([this, &visit_order, &visited_nodes](size_t root_id) {
    LeafwardDepthFirst(root_id, dag_nodes_, visit_order, visited_nodes);
  });
  return visit_order;
}

// Take in some new operations, determine an appropriate PrepForMarginalization for
// them, then append the PrepForMarginalization and the new operations to `operations`
// (in that order).
void AppendOperationsAfterPrepForMarginalization(
    GPOperationVector &operations, const GPOperationVector &new_operations) {
  if (!new_operations.empty()) {
    operations.push_back(PrepForMarginalizationOfOperations(new_operations));
    operations.insert(operations.end(), new_operations.begin(), new_operations.end());
  }
}

void GPDAG::AddPhatOperations(const GPDAGNode *node, bool rotated,
                              GPOperationVector &operations) const {
  PLVType plv_type = rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT;
  const auto parent_subsplit = node->GetBitset(rotated);
  const size_t dest_idx = GetPLVIndex(plv_type, node->Id());
  GPOperationVector new_operations;
  for (const size_t &child_idx : node->GetLeafward(rotated)) {
    const auto child_node = GetDagNode(child_idx);
    const auto child_subsplit = child_node->GetBitset();
    const auto pcsp = parent_subsplit + child_subsplit;
    Assert(gpcsp_indexer_.count(pcsp), "Non-existent PCSP index.");
    const auto gpcsp_idx = gpcsp_indexer_.at(pcsp);
    new_operations.push_back(IncrementWithWeightedEvolvedPLV{
        dest_idx, gpcsp_idx, GetPLVIndex(PLVType::P, child_idx)});
  }
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::AddRhatOperations(const GPDAGNode *node,
                              GPOperationVector &operations) const {
  const auto subsplit = node->GetBitset();
  GPOperationVector new_operations;
  IterateOverRootwardEdges(node, [this, node, &new_operations, subsplit](
                                     const bool rotated, const GPDAGNode *parent_node) {
    const auto parent_subsplit = parent_node->GetBitset(rotated);
    const auto gpcsp_idx = gpcsp_indexer_.at(parent_subsplit + subsplit);

    new_operations.push_back(IncrementWithWeightedEvolvedPLV{
        GetPLVIndex(PLVType::R_HAT, node->Id()), gpcsp_idx,
        GetPLVIndex(RPLVType(rotated), parent_node->Id())});
  });
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::OptimizeSBNParametersForASubsplit(const Bitset &subsplit,
                                              GPOperationVector &operations) const {
  if (subsplit_to_range_.count(subsplit)) {
    const auto param_range = subsplit_to_range_.at(subsplit);
    if (param_range.second - param_range.first > 1) {
      operations.push_back(
          UpdateSBNProbabilities{param_range.first, param_range.second});
    }
  }
}

void GPDAG::ScheduleBranchLengthOptimization(size_t node_id,
                                             std::unordered_set<size_t> &visited_nodes,
                                             GPOperationVector &operations) const {
  visited_nodes.insert(node_id);
  const auto node = GetDagNode(node_id);

  if (!node->IsRoot()) {
    // Compute R_HAT(s) = \sum_{t : s < t} q(s|t) P(s|t) r(t).
    operations.push_back(Zero{GetPLVIndex(PLVType::R_HAT, node_id)});
    UpdateRHat(node_id, false, operations);
    UpdateRHat(node_id, true, operations);

    // Update r(s) and r_tilde(s)
    operations.push_back(Multiply{GetPLVIndex(PLVType::R, node_id),
                                  GetPLVIndex(PLVType::R_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
    operations.push_back(Multiply{GetPLVIndex(PLVType::R_TILDE, node_id),
                                  GetPLVIndex(PLVType::R_HAT, node_id),
                                  GetPLVIndex(PLVType::P_HAT, node_id)});
  }

  if (node->IsLeaf()) return;

  operations.push_back(Zero{GetPLVIndex(PLVType::P_HAT, node_id)});
  for (size_t child_id : dag_nodes_.at(node_id)->GetLeafwardSorted()) {
    if (!visited_nodes.count(child_id)) {
      ScheduleBranchLengthOptimization(child_id, visited_nodes, operations);
    }
    OptimizeBranchLengthUpdatePHat(node_id, child_id, false, operations);
  }
  // Update r_tilde(t) = r_hat(t) \circ p_hat(t).
  operations.push_back(Multiply{GetPLVIndex(PLVType::R_TILDE, node_id),
                                GetPLVIndex(PLVType::R_HAT, node_id),
                                GetPLVIndex(PLVType::P_HAT, node_id)});

  operations.push_back(Zero{GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
  for (size_t child_id : dag_nodes_.at(node_id)->GetLeafwardRotated()) {
    if (!visited_nodes.count(child_id)) {
      ScheduleBranchLengthOptimization(child_id, visited_nodes, operations);
    }
    OptimizeBranchLengthUpdatePHat(node_id, child_id, true, operations);
  }
  // Update r(t) = r_hat(t) \circ p_hat_tilde(t).
  operations.push_back(Multiply{GetPLVIndex(PLVType::R, node_id),
                                GetPLVIndex(PLVType::R_HAT, node_id),
                                GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});

  // Update p(t).
  operations.push_back(Multiply{GetPLVIndex(PLVType::P, node_id),
                                GetPLVIndex(PLVType::P_HAT, node_id),
                                GetPLVIndex(PLVType::P_HAT_TILDE, node_id)});
};

void GPDAG::UpdateRHat(size_t node_id, bool rotated,
                       GPOperationVector &operations) const {
  const auto node = GetDagNode(node_id);
  PLVType src_plv_type = rotated ? PLVType::R_TILDE : PLVType::R;
  const auto parent_nodes =
      rotated ? node->GetRootwardRotated() : node->GetRootwardSorted();
  GPOperationVector new_operations;
  for (size_t parent_id : parent_nodes) {
    const auto parent_node = GetDagNode(parent_id);
    auto pcsp =
        rotated ? parent_node->GetBitset().RotateSubsplit() : parent_node->GetBitset();
    pcsp = pcsp + node->GetBitset();
    size_t gpcsp_idx = gpcsp_indexer_.at(pcsp);
    new_operations.push_back(
        IncrementWithWeightedEvolvedPLV{GetPLVIndex(PLVType::R_HAT, node_id), gpcsp_idx,
                                        GetPLVIndex(src_plv_type, parent_id)});
  }
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::UpdatePHatComputeLikelihood(size_t node_id, size_t child_node_id,
                                        bool rotated,
                                        GPOperationVector &operations) const {
  const auto node = GetDagNode(node_id);
  const auto child_node = GetDagNode(child_node_id);
  auto pcsp = rotated ? node->GetBitset().RotateSubsplit() : node->GetBitset();
  pcsp = pcsp + child_node->GetBitset();
  size_t gpcsp_idx = gpcsp_indexer_.at(pcsp);
  // Update p_hat(s)
  GPOperationVector new_operations;
  new_operations.push_back(IncrementWithWeightedEvolvedPLV{
      GetPLVIndex(rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT, node_id),
      gpcsp_idx,
      GetPLVIndex(PLVType::P, child_node_id),
  });
  new_operations.push_back(Likelihood{gpcsp_idx,
                                      GetPLVIndex(RPLVType(rotated), node->Id()),
                                      GetPLVIndex(PLVType::P, child_node->Id())});
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

void GPDAG::OptimizeBranchLengthUpdatePHat(size_t node_id, size_t child_node_id,
                                           bool rotated,
                                           GPOperationVector &operations) const {
  const auto node = GetDagNode(node_id);
  const auto child_node = GetDagNode(child_node_id);
  auto pcsp = rotated ? node->GetBitset().RotateSubsplit() : node->GetBitset();
  pcsp = pcsp + child_node->GetBitset();
  size_t gpcsp_idx = gpcsp_indexer_.at(pcsp);
  operations.push_back(OptimizeBranchLength{GetPLVIndex(PLVType::P, child_node_id),
                                            GetPLVIndex(RPLVType(rotated), node_id),
                                            gpcsp_idx});
  // Update p_hat(s)
  GPOperationVector new_operations;
  new_operations.push_back(IncrementWithWeightedEvolvedPLV{
      GetPLVIndex(rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT, node_id),
      gpcsp_idx,
      GetPLVIndex(PLVType::P, child_node_id),
  });
  AppendOperationsAfterPrepForMarginalization(operations, new_operations);
}

Bitset GPDAG::PerhapsRotateSubsplit(const Bitset &subsplit, bool rotated) {
  return rotated ? subsplit.RotateSubsplit() : subsplit;
}
