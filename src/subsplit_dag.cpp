// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "subsplit_dag.hpp"
#include "combinatorics.hpp"
#include "numerical_utils.hpp"
#include "sbn_probability.hpp"

SubsplitDAG::SubsplitDAG()
    : taxon_count_(0), gpcsp_count_without_fake_subsplits_(0), topology_count_(0.) {}

SubsplitDAG::SubsplitDAG(size_t taxon_count,
                         const Node::TopologyCounter &topology_counter)
    : taxon_count_(taxon_count) {
  Assert(topology_counter.size() > 0, "Empty topology counter given to SubsplitDAG.");
  Assert(topology_counter.begin()->first->LeafCount() == taxon_count,
         "Taxon count mismatch in SubsplitDAG constructor.");
  auto [gpcsp_indexer, index_to_child, rootsplits] =
      ProcessTopologyCounter(topology_counter);
  BuildNodes(index_to_child, rootsplits);
  BuildEdges(index_to_child);
  BuildDAGEdgesFromGPCSPIndexer(gpcsp_indexer);
  AddFakeSubsplitsToDAGEdgesAndParentToRange();
  CountTopologies();
}

SubsplitDAG::SubsplitDAG(const RootedTreeCollection &tree_collection)
    : SubsplitDAG(tree_collection.TaxonCount(), tree_collection.TopologyCounter()) {}

void SubsplitDAG::CountTopologies() {
  topology_count_below_ = EigenVectorXd::Ones(NodeCount());

  for (const auto &node_id : RootwardPassTraversal(true)) {
    const auto &node = GetDAGNode(node_id);
    for (const bool rotated : {true, false}) {
      // When there are no leafward nodes in the `rotated` direction, we set the number
      // of topologies for the rotation of the node to be 1.
      double per_rotated_count = node->GetLeafward(rotated).empty() ? 1. : 0.;
      // Sum options across the possible children.
      for (const auto &child_id : node->GetLeafward(rotated)) {
        per_rotated_count += topology_count_below_[child_id];
      }
      // Take the product across the number of options for the left and right branches
      // of the tree.
      topology_count_below_[node_id] *= per_rotated_count;
    }
  }
  topology_count_ = topology_count_below_[DAGRootNodeId()];
}

size_t SubsplitDAG::NodeCount() const { return dag_nodes_.size(); }

size_t SubsplitDAG::NodeCountWithoutDAGRoot() const { return NodeCount() - 1; }

double SubsplitDAG::TopologyCount() const { return topology_count_; }

size_t SubsplitDAG::RootsplitCount() const { return RootsplitIds().size(); }

size_t SubsplitDAG::GPCSPCount() const { return gpcsp_count_without_fake_subsplits_; }

size_t SubsplitDAG::GPCSPCountWithFakeSubsplits() const { return dag_edges_.size(); }

void SubsplitDAG::Print() const {
  for (const auto &dag_node : dag_nodes_) {
    std::cout << dag_node->ToString() << std::endl;
  }
}

void SubsplitDAG::PrintGPCSPIndexer() const {
  for (const auto &[pcsp, idx] : BuildGPCSPIndexer()) {
    std::cout << pcsp.PCSPToString() << ", " << idx << std::endl;
  }
}

void SubsplitDAG::PrintDAGEdges() const {
  for (const auto &[parent_child_id, gpcsp_idx] : dag_edges_) {
    const auto &[parent_id, child_id] = parent_child_id;
    std::cout << "[" << parent_id << ", " << child_id << "], " << gpcsp_idx
              << std::endl;
  }
}

std::string SubsplitDAG::ToDot(bool show_index_labels) const {
  std::stringstream string_stream;
  string_stream << "digraph g {\n";
  string_stream << "node [shape=record];\n";
  string_stream << "edge [colorscheme=dark23];\n";
  DepthFirstWithAction(
      true, SubsplitDAGTraversalAction(
                // BeforeNode
                [this, &string_stream, &show_index_labels](size_t node_id) {
                  auto bs = GetDAGNode(node_id)->GetBitset();
                  string_stream << node_id << " [label=\"<f0>"
                                << bs.SubsplitChunk(0).ToIndexSetString() << "|<f1>";
                  if (show_index_labels) {
                    string_stream << node_id;
                  }
                  string_stream << "|<f2>" << bs.SubsplitChunk(1).ToIndexSetString()
                                << "\"]\n";
                },
                // AfterNode
                [](size_t node_id) {},
                // BeforeNodeClade
                [](size_t node_id, bool rotated) {},
                // VisitEdge
                [this, &string_stream, &show_index_labels](
                    size_t node_id, size_t child_id, bool rotated) {
                  if (GetDAGNode(child_id)->IsLeaf()) {
                    string_stream << child_id << " [label=\"<f1>" << child_id
                                  << "\"]\n";
                  }
                  string_stream << "\"" << node_id << "\":";
                  string_stream << (rotated ? "f0" : "f2");
                  string_stream << "->\"";
                  string_stream << child_id << "\":f1";
                  if (show_index_labels) {
                    string_stream << " [label=\"" << GPCSPIndexOfIds(node_id, child_id);
                    if (rotated) {
                      string_stream << "\", color=1, fontcolor=1]";
                    } else {
                      string_stream << "\", color=3, fontcolor=3]";
                    }
                  }
                  string_stream << "\n";
                }));
  string_stream << "}";
  return string_stream.str();
}

BitsetSizeMap SubsplitDAG::BuildGPCSPIndexer() const {
  auto gpcsp_indexer = BitsetSizeMap();
  ReversePostorderIndexTraversal([this, &gpcsp_indexer](size_t parent_id, bool rotated,
                                                        size_t child_id,
                                                        size_t gpcsp_idx) {
    const auto parent_subsplit = GetDAGNode(parent_id)->GetBitset(rotated);
    const auto child_subsplit = GetDAGNode(child_id)->GetBitset();
    SafeInsert(gpcsp_indexer, Bitset::PCSPOfPair(parent_subsplit, child_subsplit),
               gpcsp_idx);
  });
  return gpcsp_indexer;
}

SubsplitDAGNode *SubsplitDAG::GetDAGNode(const size_t node_id) const {
  return dag_nodes_.at(node_id).get();
}

size_t SubsplitDAG::DAGRootNodeId() const { return NodeCount() - 1; }

const SizeVector &SubsplitDAG::RootsplitIds() const {
  return GetDAGNode(DAGRootNodeId())->GetLeafwardRotated();
}

size_t SubsplitDAG::GetGPCSPIndex(const Bitset &parent_subsplit,
                                  const Bitset &child_subsplit) const {
  return GPCSPIndexOfIds(subsplit_to_id_.at(parent_subsplit),
                         subsplit_to_id_.at(child_subsplit));
}

size_t SubsplitDAG::GPCSPIndexOfIds(size_t parent_id, size_t child_id) const {
  return dag_edges_.at({parent_id, child_id});
}

EigenVectorXd SubsplitDAG::BuildUniformOnTopologicalSupportPrior() const {
  EigenVectorXd q = EigenVectorXd::Ones(GPCSPCountWithFakeSubsplits());

  for (const auto &node_id : RootwardPassTraversal(true)) {
    const auto &node = GetDAGNode(node_id);
    for (const bool rotated : {false, true}) {
      if (!node->GetLeafward(rotated).empty()) {
        double per_rotated_count = 0.;
        for (const auto &child_id : node->GetLeafward(rotated)) {
          per_rotated_count += topology_count_below_[child_id];
        }
        for (const auto &child_id : node->GetLeafward(rotated)) {
          size_t gpcsp_idx = GPCSPIndexOfIds(node->Id(), child_id);
          q(gpcsp_idx) = topology_count_below_(child_id) / per_rotated_count;
        }
      }
    }
  }
  return q;
}

Node::NodePtrVec SubsplitDAG::GenerateAllTopologies() const {
  std::vector<Node::NodePtrVec> topology_below(NodeCount());

  auto GetSubtopologies = [&topology_below](const SubsplitDAGNode *node) {
    Node::NodePtrVec rotated_subtopologies, sorted_subtopologies;
    for (const bool rotated : {false, true}) {
      for (const auto &child_id : node->GetLeafward(rotated)) {
        for (const auto &subtopology : topology_below.at(child_id)) {
          rotated ? rotated_subtopologies.push_back(subtopology)
                  : sorted_subtopologies.push_back(subtopology);
        }
      }
    }
    return std::make_pair(rotated_subtopologies, sorted_subtopologies);
  };

  auto MergeTopologies = [this](size_t node_id, Node::NodePtrVec &rotated_subtopologies,
                                Node::NodePtrVec &sorted_subtopologies) {
    Node::NodePtrVec topologies;
    for (const auto &rotated_subtopology : rotated_subtopologies) {
      for (const auto &sorted_subtopology : sorted_subtopologies) {
        Node::NodePtr new_topology =
            Node::Join(sorted_subtopology, rotated_subtopology, node_id);
        topologies.push_back(new_topology);
      }
    }
    if (node_id == DAGRootNodeId()) {
      // DAG root node has no `sorted_subtopologies`, so loop above yields empty
      // `topologies` vector.
      return rotated_subtopologies;
    }
    return topologies;
  };

  for (const auto &node_id : RootwardPassTraversal(true)) {
    const auto &node = GetDAGNode(node_id);
    if (node->IsLeaf()) {
      topology_below.at(node_id).push_back(Node::Leaf(node_id));
    } else {
      auto [rotated_topologies, sorted_topologies] = GetSubtopologies(node);
      topology_below[node_id] =
          MergeTopologies(node_id, rotated_topologies, sorted_topologies);
    }
  }

  const auto &topologies = topology_below.at(DAGRootNodeId());
  Assert(topologies.size() == TopologyCount(),
         "The realized number of topologies does not match the expected count.");

  // We return a deep copy of every Polished topology to avoid loops in the pointer
  // structure. Such loops can create problems when we Polish the topologies one at a
  // time: polishing a second topology can change the numbering of a previous topology.
  // This is checked for in the "GPInstance: GenerateCompleteRootedTreeCollection" test.
  Node::NodePtrVec final_topologies;
  final_topologies.reserve(topologies.size());
  for (auto &topology : topologies) {
    topology->Polish();
    final_topologies.push_back(topology->DeepCopy());
  }

  return final_topologies;
}

EigenVectorXd SubsplitDAG::BuildUniformOnAllTopologiesPrior() const {
  EigenVectorXd result = EigenVectorXd::Zero(GPCSPCountWithFakeSubsplits());
  for (const auto &[parent_child_id, gpcsp_idx] : dag_edges_) {
    const auto &[parent_id, child_id] = parent_child_id;
    size_t child0_taxon_count =
        GetDAGNode(child_id)->GetBitset().SubsplitGetChild(0).Count();
    size_t child1_taxon_count =
        GetDAGNode(child_id)->GetBitset().SubsplitGetChild(1).Count();
    result(gpcsp_idx) = Combinatorics::LogChildSubsplitCountRatio(child0_taxon_count,
                                                                  child1_taxon_count);
  }
  NumericalUtils::Exponentiate(result);
  return result;
}

void SubsplitDAG::IterateOverRealNodes(const NodeLambda &f) const {
  Assert(taxon_count_ < NodeCount(), "No real DAG nodes!");
  for (auto it = dag_nodes_.cbegin() + taxon_count_; it < dag_nodes_.cend() - 1; it++) {
    f((*it).get());
  }
}

void SubsplitDAG::IterateOverLeafwardEdges(const SubsplitDAGNode *node, bool rotated,
                                           const NodeLambda &f) const {
  for (const size_t child_id : node->GetLeafward(rotated)) {
    f(GetDAGNode(child_id));
  }
}

void SubsplitDAG::IterateOverLeafwardEdges(const SubsplitDAGNode *node,
                                           const EdgeDestinationLambda &f) const {
  for (bool rotated : {false, true}) {
    for (const size_t child_id : node->GetLeafward(rotated)) {
      f(rotated, GetDAGNode(child_id));
    }
  }
}

void SubsplitDAG::IterateOverLeafwardEdgesAndChildren(
    const SubsplitDAGNode *node, const EdgeAndNodeLambda &f) const {
  IterateOverLeafwardEdges(
      node, [this, &node, &f](bool rotated, const SubsplitDAGNode *child) {
        f(GPCSPIndexOfIds(node->Id(), child->Id()), rotated, child->Id());
      });
}

void SubsplitDAG::IterateOverRootwardEdges(const SubsplitDAGNode *node,
                                           const EdgeDestinationLambda &f) const {
  if (not node->IsRootsplit()) {
    for (bool rotated : {false, true}) {
      for (const size_t parent_id : node->GetRootward(rotated)) {
        f(rotated, GetDAGNode(parent_id));
      }
    }
  }
}

void SubsplitDAG::IterateOverRootwardEdgesAndParents(const SubsplitDAGNode *node,
                                                     const EdgeAndNodeLambda &f) const {
  IterateOverRootwardEdges(
      node, [this, &node, &f](bool rotated, const SubsplitDAGNode *parent) {
        f(GPCSPIndexOfIds(parent->Id(), node->Id()), rotated, parent->Id());
      });
}

RootedIndexerRepresentation SubsplitDAG::IndexerRepresentationOf(
    const BitsetSizeMap &indexer, const Node::NodePtr &topology,
    size_t default_index) const {
  return RootedSBNMaps::IndexerRepresentationOf(indexer, topology, default_index);
}

EigenVectorXd SubsplitDAG::UnconditionalNodeProbabilities(
    EigenConstVectorXdRef normalized_sbn_parameters) const {
  EigenVectorXd node_probabilities(NodeCount());
  node_probabilities.setZero();
  node_probabilities[DAGRootNodeId()] = 1.;

  ReversePostorderIndexTraversal([this, &node_probabilities,
                                  &normalized_sbn_parameters](
                                     const size_t parent_id, const bool,
                                     const size_t child_id, const size_t gpcsp_idx) {
    const double child_probability_given_parent = normalized_sbn_parameters[gpcsp_idx];
    Assert(child_probability_given_parent >= 0. && child_probability_given_parent <= 1.,
           "UnconditionalNodeProbabilities: got an out-of-range probability. Are these "
           "normalized and in linear space?");
    const double parent_probability = node_probabilities[parent_id];
    node_probabilities[child_id] += parent_probability * child_probability_given_parent;
  });

  return node_probabilities;
}

BitsetDoubleMap SubsplitDAG::UnconditionalSubsplitProbabilities(
    EigenConstVectorXdRef normalized_sbn_parameters) const {
  auto node_probabilities = UnconditionalNodeProbabilities(normalized_sbn_parameters);
  BitsetDoubleMap subsplit_probability_map;
  for (size_t node_id = 0; node_id < node_probabilities.size(); node_id++) {
    const auto &subsplit_bitset = GetDAGNode(node_id)->GetBitset();
    if (!subsplit_bitset.SubsplitIsFake()) {
      SafeInsert(subsplit_probability_map, subsplit_bitset,
                 node_probabilities[node_id]);
    }
  }
  return subsplit_probability_map;
}

EigenVectorXd SubsplitDAG::InvertedGPCSPProbabilities(
    EigenConstVectorXdRef normalized_sbn_parameters,
    EigenConstVectorXdRef node_probabilities) const {
  EigenVectorXd inverted_probabilities =
      EigenVectorXd(normalized_sbn_parameters.size());
  inverted_probabilities.setOnes();
  ReversePostorderIndexTraversal(
      [this, &node_probabilities, &normalized_sbn_parameters, &inverted_probabilities](
          const size_t parent_id, const bool, const size_t child_id,
          const size_t gpcsp_idx) {
        // The traversal doesn't set the rootsplit probabilities, but those are always 1
        // (there is only one "parent" of a rootsplit).
        if (parent_id != DAGRootNodeId()) {
          // For a PCSP t -> s:
          inverted_probabilities[gpcsp_idx] =         // P(t|s)
              node_probabilities[parent_id] *         // P(t)
              normalized_sbn_parameters[gpcsp_idx] /  // P(s|t)
              node_probabilities[child_id];           // P(s)
        }
      });
  return inverted_probabilities;
}

std::vector<Bitset> SubsplitDAG::GetChildSubsplits(const SizeBitsetMap &index_to_child,
                                                   const Bitset &subsplit,
                                                   bool include_fake_subsplits) {
  std::vector<Bitset> children_subsplits;
  if (parent_to_range_.count(subsplit) > 0) {
    const auto [start, stop] = parent_to_range_.at(subsplit);
    for (auto idx = start; idx < stop; idx++) {
      children_subsplits.push_back(index_to_child.at(idx));
    }
  } else if (include_fake_subsplits) {
    // This method is designed to be called before calling
    // AddFakeSubsplitsToDAGEdgesAndParentToRange. In that case, if the second chunk
    // of the subsplit is just a single taxon, the subsplit will not map to any value in
    // parent_to_range_.
    //
    // But we still need to create and connect to fake subsplits in the DAG. So, here we
    // make a fake child subsplit.
    children_subsplits.push_back(Bitset::FakeChildSubsplit(subsplit));
  }

  return children_subsplits;
}

std::tuple<BitsetSizeMap, SizeBitsetMap, BitsetVector>
SubsplitDAG::ProcessTopologyCounter(const Node::TopologyCounter &topology_counter) {
  BitsetSizeMap gpcsp_indexer;
  SizeBitsetMap index_to_child;
  BitsetVector rootsplits;
  std::tie(rootsplits, gpcsp_indexer, index_to_child, parent_to_range_,
           gpcsp_count_without_fake_subsplits_) =
      SBNMaps::BuildIndexerBundle(RootedSBNMaps::RootsplitCounterOf(topology_counter),
                                  RootedSBNMaps::PCSPCounterOf(topology_counter));
  return {gpcsp_indexer, index_to_child, rootsplits};
}

void SubsplitDAG::CreateAndInsertNode(const Bitset &subsplit) {
  size_t id = NodeCount();
  SafeInsert(subsplit_to_id_, subsplit, id);
  SafeInsert(subsplit_to_id_, subsplit.RotateSubsplit(), id);
  dag_nodes_.push_back(std::make_unique<SubsplitDAGNode>(id, subsplit));
}

void SubsplitDAG::ConnectNodes(const SizeBitsetMap &index_to_child, size_t id,
                               bool rotated) {
  const auto node = GetDAGNode(id);
  // Retrieve children subsplits, set edge relation.
  const Bitset subsplit = node->GetBitset(rotated);
  const auto children = GetChildSubsplits(index_to_child, subsplit, true);
  for (const auto &child_subsplit : children) {
    const auto child_node = GetDAGNode(subsplit_to_id_.at(child_subsplit));
    if (rotated) {
      node->AddLeafwardRotated(child_node->Id());
      child_node->AddRootwardRotated(node->Id());
    } else {
      node->AddLeafwardSorted(child_node->Id());
      child_node->AddRootwardSorted(node->Id());
    }
  }
}

void SubsplitDAG::BuildNodesDepthFirst(const SizeBitsetMap &index_to_child,
                                       const Bitset &subsplit,
                                       std::unordered_set<Bitset> &visited_subsplits) {
  visited_subsplits.insert(subsplit);
  for (bool rotated : {false, true}) {
    for (const auto &child_subsplit : GetChildSubsplits(
             index_to_child, PerhapsRotateSubsplit(subsplit, rotated), false)) {
      if (visited_subsplits.count(child_subsplit) == 0) {
        BuildNodesDepthFirst(index_to_child, child_subsplit, visited_subsplits);
      }
    }
  }
  CreateAndInsertNode(subsplit);
}

void SubsplitDAG::BuildNodes(const SizeBitsetMap &index_to_child,
                             const BitsetVector &rootsplits) {
  std::unordered_set<Bitset> visited_subsplits;
  // We will create fake subsplits and insert to dag_nodes_.
  // Fake subsplits (i.e. leaf nodes) will take IDs in [0, taxon_count_).
  for (size_t taxon_idx = 0; taxon_idx < taxon_count_; taxon_idx++) {
    CreateAndInsertNode(
        Bitset::FakeSubsplit(Bitset::Singleton(taxon_count_, taxon_idx)));
  }
  // Then we add the remaining nodes.
  // The root splits will take on the higher IDs compared to the non-rootsplits.
  for (const auto &rootsplit : rootsplits) {
    BuildNodesDepthFirst(index_to_child, rootsplit, visited_subsplits);
  }

  // Finally, we add the DAG root node.
  CreateAndInsertNode(Bitset::DAGRootSubsplitOfTaxonCount(taxon_count_));
}

void SubsplitDAG::BuildEdges(const SizeBitsetMap &index_to_child) {
  // Connect every node except for the DAG root node.
  for (size_t i = taxon_count_; i < DAGRootNodeId(); i++) {
    ConnectNodes(index_to_child, i, false);
    ConnectNodes(index_to_child, i, true);
  }
  // Connect the DAG root node.
  ConnectNodes(index_to_child, DAGRootNodeId(), true);
}

void SubsplitDAG::BuildDAGEdgesFromGPCSPIndexer(BitsetSizeMap &gpcsp_indexer) {
  for (const auto &[gpcsp, index] : gpcsp_indexer) {
    Assert(gpcsp.size() == 3 * taxon_count_,
           "All edges should be bitsets with size 3 times taxon_count_.");
    const auto parent_id = subsplit_to_id_.at(gpcsp.PCSPParentSubsplit());
    const auto child_id = subsplit_to_id_.at(gpcsp.PCSPChildSubsplit());
    SafeInsert(dag_edges_, {parent_id, child_id}, index);
  }
}

void SubsplitDAG::AddFakeSubsplitsToDAGEdgesAndParentToRange() {
  for (size_t i = 0; i < taxon_count_; i++) {
    const auto current_bitset = dag_nodes_.at(i)->GetBitset();
    IterateOverRootwardEdges(
        GetDAGNode(i),
        [this, current_bitset](const bool rotated, const SubsplitDAGNode *node) {
          SafeInsert(parent_to_range_, node->GetBitset(rotated),
                     {dag_edges_.size(), dag_edges_.size() + 1});
          SafeInsert(dag_edges_, {node->Id(), subsplit_to_id_.at(current_bitset)},
                     dag_edges_.size());
        });
  }
}

void SubsplitDAG::RootwardDepthFirst(size_t id, SizeVector &visit_order,
                                     std::unordered_set<size_t> &visited_nodes) const {
  SafeInsert(visited_nodes, id);
  const auto &node = GetDAGNode(id);
  for (size_t parent_id : node->GetRootwardSorted()) {
    if (visited_nodes.count(parent_id) == 0) {
      RootwardDepthFirst(parent_id, visit_order, visited_nodes);
    }
  }
  for (size_t parent_id : node->GetRootwardRotated()) {
    if (visited_nodes.count(parent_id) == 0) {
      RootwardDepthFirst(parent_id, visit_order, visited_nodes);
    }
  }
  visit_order.push_back(id);
}

void SubsplitDAG::LeafwardDepthFirst(size_t id, SizeVector &visit_order,
                                     std::unordered_set<size_t> &visited_nodes) const {
  SafeInsert(visited_nodes, id);
  for (size_t child_id : GetDAGNode(id)->GetLeafwardSorted()) {
    if (visited_nodes.count(child_id) == 0) {
      LeafwardDepthFirst(child_id, visit_order, visited_nodes);
    }
  }
  for (size_t child_id : GetDAGNode(id)->GetLeafwardRotated()) {
    if (visited_nodes.count(child_id) == 0) {
      LeafwardDepthFirst(child_id, visit_order, visited_nodes);
    }
  }
  visit_order.push_back(id);
}

SizeVector SubsplitDAG::LeafwardPassTraversal(bool include_dag_root_node) const {
  SizeVector visit_order;
  std::unordered_set<size_t> visited_nodes;
  if (!include_dag_root_node) {
    SafeInsert(visited_nodes, DAGRootNodeId());
  }
  for (size_t leaf_id = 0; leaf_id < taxon_count_; leaf_id++) {
    RootwardDepthFirst(leaf_id, visit_order, visited_nodes);
  }
  return visit_order;
}

SizeVector SubsplitDAG::RootwardPassTraversal(bool include_dag_root_node) const {
  SizeVector visit_order;
  std::unordered_set<size_t> visited_nodes;
  for (const auto &rootsplit_id : RootsplitIds()) {
    LeafwardDepthFirst(rootsplit_id, visit_order, visited_nodes);
  }
  if (include_dag_root_node) {
    visit_order.push_back(DAGRootNodeId());
  }
  return visit_order;
}

SizeVector SubsplitDAG::ReversePostorderTraversal() const {
  auto visit_order = RootwardPassTraversal(true);
  std::reverse(visit_order.begin(), visit_order.end());
  return visit_order;
}

void SubsplitDAG::ReversePostorderIndexTraversal(
    ParentRotationChildEdgeLambda f) const {
  for (const auto node_id : ReversePostorderTraversal()) {
    IterateOverLeafwardEdgesAndChildren(
        GetDAGNode(node_id), [&f, &node_id](const size_t gpcsp_idx, const bool rotated,
                                            const size_t child_id) {
          f(node_id, rotated, child_id, gpcsp_idx);
        });
  }
}

Bitset SubsplitDAG::PerhapsRotateSubsplit(const Bitset &subsplit, bool rotated) {
  return rotated ? subsplit.RotateSubsplit() : subsplit;
}
