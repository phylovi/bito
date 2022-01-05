// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

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

void SubsplitDAG::PrintParentToRange() const {
  for (const auto &[subsplit, pcsp_range] : parent_to_range_) {
    std::cout << subsplit.SubsplitToString() << ": [" << pcsp_range.first << ", "
              << pcsp_range.second << "]" << std::endl;
  }
}

void SubsplitDAG::ToDot(const std::string file_path, bool show_index_labels) const {
  std::ofstream out_file(file_path);
  out_file << ToDot(show_index_labels);
  out_file.close();
}

std::string SubsplitDAG::ToDot(bool show_index_labels) const {
  std::stringstream string_stream;
  string_stream << "digraph g {\n";
  string_stream << "node [shape=record];\n";
  string_stream << "edge [colorscheme=dark23];\n";
  DepthFirstWithAction(
      {DAGRootNodeId()},
      SubsplitDAGTraversalAction(
          // BeforeNode
          [this, &string_stream, &show_index_labels](size_t node_id) {
            const auto &node = GetDAGNode(node_id);
            if (node->IsDAGRootNode()) {
              string_stream << node_id << " [label=\"<f0>&rho;\"]\n";
              return;
            }
            auto bs = node->GetBitset();
            string_stream << node_id << " [label=\"<f0>"
                          << bs.SubsplitGetClade(0).ToVectorOfSetBitsAsString()
                          << "|<f1>";
            if (show_index_labels) {
              string_stream << node_id;
            }
            string_stream << "|<f2>"
                          << bs.SubsplitGetClade(1).ToVectorOfSetBitsAsString()
                          << "\"]\n";
          },
          // AfterNode
          [](size_t node_id) {},
          // BeforeNodeClade
          [](size_t node_id, bool rotated) {},
          // VisitEdge
          [this, &string_stream, &show_index_labels](size_t node_id, size_t child_id,
                                                     bool rotated) {
            if (GetDAGNode(child_id)->IsLeaf()) {
              string_stream << child_id << " [label=\"<f1>" << child_id << "\"]\n";
            }
            string_stream << "\"" << node_id << "\":";
            string_stream << (rotated ? "f0" : "f2");
            string_stream << "->\"";
            string_stream << child_id << "\":f1";
            if (show_index_labels) {
              string_stream << " [label=\"" << GPCSPIndexOfIds(node_id, child_id);
              if (rotated) {
                string_stream << "\", color=1, fontcolor=1";
              } else {
                string_stream << "\", color=3, fontcolor=3";
              }
              if (GetDAGNode(node_id)->IsDAGRootNode()) {
                string_stream << ",style=dashed]";
              } else {
                string_stream << "]";
              }
            } else {
              if (GetDAGNode(node_id)->IsDAGRootNode()) {
                string_stream << "[style=dashed]";
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
    SafeInsert(gpcsp_indexer, Bitset::PCSP(parent_subsplit, child_subsplit), gpcsp_idx);
  });
  return gpcsp_indexer;
}

SubsplitDAGNode *SubsplitDAG::GetDAGNode(const size_t node_id) const {
  Assert(ContainsNode(node_id),
         "Node with the given node_id does not exist in SubsplitDAG::GetDAGNode.");
  return dag_nodes_.at(node_id).get();
}

size_t SubsplitDAG::GetDAGNodeId(const Bitset &subsplit) const {
  Assert(ContainsNode(subsplit),
         "Node with the given subsplit does not exist in SubsplitDAG::GetDAGNodeId.");
  return subsplit_to_id_.at(subsplit);
}

size_t SubsplitDAG::DAGRootNodeId() const { return NodeCount() - 1; }

const SizeVector &SubsplitDAG::RootsplitIds() const {
  return GetDAGNode(DAGRootNodeId())->GetLeafwardRotated();
}

size_t SubsplitDAG::GetGPCSPIndex(const Bitset &parent_subsplit,
                                  const Bitset &child_subsplit) const {
  return GPCSPIndexOfIds(GetDAGNodeId(parent_subsplit), GetDAGNodeId(child_subsplit));
}

size_t SubsplitDAG::GPCSPIndexOfIds(size_t parent_id, size_t child_id) const {
  return dag_edges_.at({parent_id, child_id});
}

SizePair SubsplitDAG::GetEdgeRange(const Bitset &subsplit, const bool rotated) const {
  Assert(ContainsNode(subsplit),
         "Node with the given subsplit does not exist in SubsplitDAG::GetEdgeRange.");
  return parent_to_range_.at(PerhapsSubsplitRotate(subsplit, rotated));
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
    std::ignore = parent_id;
    size_t child0_taxon_count =
        GetDAGNode(child_id)->GetBitset().SubsplitGetCladeByBinaryOrder(0).Count();
    size_t child1_taxon_count =
        GetDAGNode(child_id)->GetBitset().SubsplitGetCladeByBinaryOrder(1).Count();
    //
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
  if (!node->IsRootsplit()) {
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

void SubsplitDAG::IterateOverParentAndChildAndLeafwardEdges(
    const SubsplitDAGNode *node, const ParentRotationChildEdgeLambda &f) const {
  IterateOverLeafwardEdges(
      node, [this, &node, &f](bool rotated, const SubsplitDAGNode *child) {
        f(node->Id(), rotated, child->Id(), GPCSPIndexOfIds(node->Id(), child->Id()));
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

  ReversePostorderIndexTraversal([&node_probabilities, &normalized_sbn_parameters](
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
  for (size_t node_id = 0;
       static_cast<Eigen::Index>(node_id) < node_probabilities.size(); node_id++) {
    const auto &subsplit_bitset = GetDAGNode(node_id)->GetBitset();
    if (node_id != DAGRootNodeId() && !subsplit_bitset.SubsplitIsLeaf()) {
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
    // AddFakeSubsplitsToDAGEdgesAndParentToRange. In that case, if the second clade
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
  size_t node_id = NodeCount();
  dag_nodes_.push_back(std::make_unique<SubsplitDAGNode>(node_id, subsplit));
  SafeInsert(subsplit_to_id_, subsplit, node_id);
}

void SubsplitDAG::CreateAndInsertEdge(const size_t parent_id, const size_t child_id,
                                      bool rotated) {
  Assert(ContainsNode(parent_id), "Node with the given parent_id does not exist.");
  Assert(ContainsNode(child_id), "Node with the given child_id does not exist.");
  ConnectGivenNodes(parent_id, child_id, rotated);
  SafeInsert(dag_edges_, {parent_id, child_id}, GPCSPCountWithFakeSubsplits());
}

void SubsplitDAG::ConnectGivenNodes(const size_t parent_id, const size_t child_id,
                                    bool rotated) {
  const auto parent_node = GetDAGNode(parent_id);
  const auto child_node = GetDAGNode(child_id);
  if (rotated) {
    parent_node->AddLeafwardRotated(child_node->Id());
    child_node->AddRootwardRotated(parent_node->Id());
  } else {
    parent_node->AddLeafwardSorted(child_node->Id());
    child_node->AddRootwardSorted(parent_node->Id());
  }
}

void SubsplitDAG::ConnectNodes(const SizeBitsetMap &index_to_child, size_t node_id,
                               bool rotated) {
  // Retrieve children subsplits, set edge relation.
  const Bitset subsplit = GetDAGNode(node_id)->GetBitset(rotated);
  const auto children = GetChildSubsplits(index_to_child, subsplit, true);
  for (const auto &child_subsplit : children) {
    ConnectGivenNodes(node_id, GetDAGNodeId(child_subsplit), rotated);
  }
}

void SubsplitDAG::BuildNodesDepthFirst(const SizeBitsetMap &index_to_child,
                                       const Bitset &subsplit,
                                       std::unordered_set<Bitset> &visited_subsplits) {
  visited_subsplits.insert(subsplit);
  for (bool rotated : {false, true}) {
    for (const auto &child_subsplit : GetChildSubsplits(
             index_to_child, PerhapsSubsplitRotate(subsplit, rotated), false)) {
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
  for (size_t node_id = taxon_count_; node_id < DAGRootNodeId(); node_id++) {
    ConnectNodes(index_to_child, node_id, true);
    ConnectNodes(index_to_child, node_id, false);
  }
  // Connect the DAG root node.
  ConnectNodes(index_to_child, DAGRootNodeId(), true);
}

void SubsplitDAG::BuildDAGEdgesFromGPCSPIndexer(BitsetSizeMap &gpcsp_indexer) {
  for (const auto &[gpcsp, index] : gpcsp_indexer) {
    Assert(gpcsp.size() == 3 * taxon_count_,
           "All edges should be bitsets with size 3 times taxon_count_.");
    const auto parent_id = GetDAGNodeId(gpcsp.PCSPGetParentSubsplit());
    const auto child_id = GetDAGNodeId(gpcsp.PCSPGetChildSubsplit());
    SafeInsert(dag_edges_, {parent_id, child_id}, index);
  }
}

void SubsplitDAG::AddFakeSubsplitsToDAGEdgesAndParentToRange() {
  for (size_t node_id = 0; node_id < taxon_count_; node_id++) {
    const auto current_bitset(GetDAGNode(node_id)->GetBitset());
    IterateOverRootwardEdges(GetDAGNode(node_id), [this, current_bitset](
                                                      const bool rotated,
                                                      const SubsplitDAGNode *node) {
      SafeInsert(parent_to_range_, node->GetBitset(rotated),
                 {GPCSPCountWithFakeSubsplits(), GPCSPCountWithFakeSubsplits() + 1});
      SafeInsert(dag_edges_, {node->Id(), GetDAGNodeId(current_bitset)},
                 GPCSPCountWithFakeSubsplits());
    });
  }
}

void SubsplitDAG::RootwardDepthFirst(size_t node_id, SizeVector &visit_order,
                                     std::unordered_set<size_t> &visited_nodes) const {
  SafeInsert(visited_nodes, node_id);
  const auto &node = GetDAGNode(node_id);
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
  visit_order.push_back(node_id);
}

void SubsplitDAG::LeafwardDepthFirst(size_t node_id, SizeVector &visit_order,
                                     std::unordered_set<size_t> &visited_nodes) const {
  SafeInsert(visited_nodes, node_id);
  for (size_t child_id : GetDAGNode(node_id)->GetLeafwardSorted()) {
    if (visited_nodes.count(child_id) == 0) {
      LeafwardDepthFirst(child_id, visit_order, visited_nodes);
    }
  }
  for (size_t child_id : GetDAGNode(node_id)->GetLeafwardRotated()) {
    if (visited_nodes.count(child_id) == 0) {
      LeafwardDepthFirst(child_id, visit_order, visited_nodes);
    }
  }
  visit_order.push_back(node_id);
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

Bitset SubsplitDAG::PerhapsSubsplitRotate(const Bitset &subsplit, bool rotated) const {
  return rotated ? subsplit.SubsplitRotate() : subsplit;
}

bool SubsplitDAG::ContainsNode(const Bitset &subsplit) const {
  return subsplit_to_id_.find(subsplit) != subsplit_to_id_.end();
}

bool SubsplitDAG::ContainsNode(const size_t node_id) const {
  return node_id < NodeCount();
}

bool SubsplitDAG::ContainsEdge(const size_t parent_id, const size_t child_id) const {
  return dag_edges_.find({parent_id, child_id}) != dag_edges_.end();
}

std::pair<SizeVector, SizeVector> SubsplitDAG::BuildParentIdVector(
    const Bitset &subsplit) const {
  SizeVector rotated_parents, sorted_parents;
  for (const auto &[potential_parent_subsplit, node_id] : subsplit_to_id_) {
    if (subsplit.SubsplitIsRotatedChildOf(potential_parent_subsplit)) {
      rotated_parents.push_back(node_id);
    } else if (subsplit.SubsplitIsSortedChildOf(potential_parent_subsplit)) {
      sorted_parents.push_back(node_id);
    }
  }
  return {rotated_parents, sorted_parents};
}

std::pair<SizeVector, SizeVector> SubsplitDAG::BuildChildIdVector(
    const Bitset &subsplit) const {
  SizeVector rotated_children, sorted_children;
  for (const auto &[potential_child_subsplit, node_id] : subsplit_to_id_) {
    if (potential_child_subsplit.SubsplitIsRotatedChildOf(subsplit)) {
      rotated_children.push_back(node_id);
    } else if (potential_child_subsplit.SubsplitIsSortedChildOf(subsplit)) {
      sorted_children.push_back(node_id);
    }
  }
  return {rotated_children, sorted_children};
}

bool SubsplitDAG::IsValidNewNodePair(const Bitset &parent_subsplit,
                                     const Bitset &child_subsplit) const {
  const auto [rotated_parents_of_parent, sorted_parents_of_parent] =
      BuildParentIdVector(parent_subsplit);
  const auto [rotated_children_of_parent, sorted_children_of_parent] =
      BuildChildIdVector(parent_subsplit);
  const auto [rotated_children_of_child, sorted_children_of_child] =
      BuildChildIdVector(child_subsplit);
  return parent_subsplit.size() == 2 * taxon_count_ &&
         child_subsplit.size() == 2 * taxon_count_ &&
         !(rotated_parents_of_parent.empty() && sorted_parents_of_parent.empty()) &&
         ((child_subsplit.SubsplitIsRotatedChildOf(parent_subsplit) &&
           !sorted_children_of_parent.empty()) ||
          (child_subsplit.SubsplitIsSortedChildOf(parent_subsplit) &&
           !rotated_children_of_parent.empty())) &&
         !rotated_children_of_child.empty() && !sorted_children_of_child.empty();
}

void SubsplitDAG::ConnectChildToAllChildren(const Bitset &child_subsplit,
                                            SizeVector &new_edge_idxs) {
  const auto [rotated_children_of_child, sorted_children_of_child] =
      BuildChildIdVector(child_subsplit);
  for (const auto &[children_of_child, rotated] :
       std::vector<std::pair<SizeVector, bool>>{{rotated_children_of_child, true},
                                                {sorted_children_of_child, false}}) {
    SafeInsert(parent_to_range_, PerhapsSubsplitRotate(child_subsplit, rotated),
               {GPCSPCountWithFakeSubsplits(),
                GPCSPCountWithFakeSubsplits() + children_of_child.size()});
    for (const size_t child_of_child_id : children_of_child) {
      new_edge_idxs.push_back(GPCSPCountWithFakeSubsplits());
      CreateAndInsertEdge(GetDAGNodeId(child_subsplit), child_of_child_id, rotated);
    }
  }
}

void SubsplitDAG::ConnectParentToAllChildrenExcept(const Bitset &parent_subsplit,
                                                   const Bitset &child_subsplit,
                                                   SizeVector &new_edge_idxs) {
  const auto [rotated_children_of_parent, sorted_children_of_parent] =
      BuildChildIdVector(parent_subsplit);
  for (const auto &[children_of_parent, rotated] :
       std::vector<std::pair<SizeVector, bool>>{{rotated_children_of_parent, true},
                                                {sorted_children_of_parent, false}}) {
    SafeInsert(parent_to_range_, PerhapsSubsplitRotate(parent_subsplit, rotated),
               {GPCSPCountWithFakeSubsplits(),
                GPCSPCountWithFakeSubsplits() + children_of_parent.size()});
    for (const size_t child_of_parent_id : children_of_parent) {
      if (child_of_parent_id != GetDAGNodeId(child_subsplit)) {
        new_edge_idxs.push_back(GPCSPCountWithFakeSubsplits());
        CreateAndInsertEdge(GetDAGNodeId(parent_subsplit), child_of_parent_id, rotated);
      }
    }
  }
}

void SubsplitDAG::ConnectChildToAllParentsExcept(const Bitset &parent_subsplit,
                                                 const Bitset &child_subsplit,
                                                 SizeVector &new_edge_idxs) {
  const auto [rotated_parents_of_child, sorted_parents_of_child] =
      BuildParentIdVector(child_subsplit);
  for (const auto &[parents_of_child, rotated] :
       std::vector<std::pair<SizeVector, bool>>{{rotated_parents_of_child, true},
                                                {sorted_parents_of_child, false}}) {
    for (const size_t parent_of_child_id : parents_of_child) {
      if (parent_of_child_id != GetDAGNodeId(parent_subsplit)) {
        new_edge_idxs.push_back(GPCSPCountWithFakeSubsplits());
        CreateAndInsertEdge(parent_of_child_id, GetDAGNodeId(child_subsplit), rotated);
      }
    }
  }
}

void SubsplitDAG::ConnectParentToAllParents(const Bitset &parent_subsplit,
                                            SizeVector &new_edge_idxs) {
  const auto [rotated_parents_of_parent, sorted_parents_of_parent] =
      BuildParentIdVector(parent_subsplit);
  for (const auto &[parents_of_parent, rotated] :
       std::vector<std::pair<SizeVector, bool>>{{rotated_parents_of_parent, true},
                                                {sorted_parents_of_parent, false}}) {
    for (const size_t parent_of_parent_id : parents_of_parent) {
      new_edge_idxs.push_back(GPCSPCountWithFakeSubsplits());
      CreateAndInsertEdge(parent_of_parent_id, GetDAGNodeId(parent_subsplit), rotated);
    }
  }
}

SubsplitDAG::NodeAdditionResult SubsplitDAG::AddNodePair(const Bitset &parent_subsplit,
                                                         const Bitset &child_subsplit) {
  Assert(
      IsValidNewNodePair(parent_subsplit, child_subsplit),
      "The given pair of nodes is incompatible with DAG in SubsplitDAG::AddNodePair.");
  SizeVector new_node_ids, new_edge_idxs, node_reindexer, edge_reindexer;
  const bool parent_is_new = !ContainsNode(parent_subsplit);
  const bool child_is_new = !ContainsNode(child_subsplit);
  // If both the parent and child exists, return new_node_ids and  new_edge_idxs as
  // empty, and node_reindexer and edge_reindexer as default reindexers.
  if (!parent_is_new && !child_is_new) {
    // Return default reindexers if both nodes already exist.
    node_reindexer = Reindexer::IdentityReindexer(NodeCount());
    edge_reindexer = Reindexer::IdentityReindexer(GPCSPCountWithFakeSubsplits());
    return {new_node_ids, new_edge_idxs, node_reindexer, edge_reindexer};
  }
  // Note: `prev_node_count` is a marker so that we know what the DAG root node id is
  // (`prev_node_count - 1`).
  const size_t prev_node_count = NodeCount();
  if (child_is_new) {
    CreateAndInsertNode(child_subsplit);
    new_node_ids.push_back(GetDAGNodeId(child_subsplit));
    // Don't reindex these edges.
    ConnectChildToAllChildren(child_subsplit, new_edge_idxs);
  }
  if (parent_is_new) {
    CreateAndInsertNode(parent_subsplit);
    new_node_ids.push_back(GetDAGNodeId(parent_subsplit));
    // Don't reindex these edges.
    ConnectParentToAllChildrenExcept(parent_subsplit, child_subsplit, new_edge_idxs);
  }
  // Note: `prev_edge_count` is a marker conveying where we need to start
  // reindexing edge idxs.
  // Edges are only reindexed if the parent node already existed in the DAG
  // (so as to ensure that edges descending from the same node clade have
  // contiguous idxs).
  size_t prev_edge_count = GPCSPCountWithFakeSubsplits();
  // Connect the given parent node to the given child node.
  new_edge_idxs.push_back(GPCSPCountWithFakeSubsplits());
  CreateAndInsertEdge(GetDAGNodeId(parent_subsplit), GetDAGNodeId(child_subsplit),
                      child_subsplit.SubsplitIsRotatedChildOf(parent_subsplit));
  // Don't reindex the edge between the given parent and child if the parent is new.
  if (parent_is_new) {
    prev_edge_count = GPCSPCountWithFakeSubsplits();
  }
  if (child_is_new) {
    // Reindex these edges.
    ConnectChildToAllParentsExcept(parent_subsplit, child_subsplit, new_edge_idxs);
  }
  if (parent_is_new) {
    // Reindex these edges.
    ConnectParentToAllParents(parent_subsplit, new_edge_idxs);
  }
  // Create reindexers.
  node_reindexer = BuildNodeReindexer(prev_node_count);
  edge_reindexer = BuildEdgeReindexer(prev_edge_count);
  // Update the ids in new_node_ids and new_edge_idxs according to the reindexers.
  Reindexer::RemapIdVector(new_node_ids, node_reindexer);
  Reindexer::RemapIdVector(new_edge_idxs, edge_reindexer);
  // Update fields in the Subsplit DAG according to the reindexers.
  RemapNodeIds(node_reindexer);
  RemapEdgeIdxs(edge_reindexer);
  // Recount topologies to update `topology_count_below_`.
  CountTopologies();
  return {new_node_ids, new_edge_idxs, node_reindexer, edge_reindexer};
}

SizeVector SubsplitDAG::BuildNodeReindexer(const size_t prev_node_count) {
  SizeVector node_reindexer = Reindexer::IdentityReindexer(NodeCount());
  size_t running_traversal_idx = taxon_count_;
  size_t dag_root_node_id = prev_node_count - 1;
  // Build node_reindexer by traversing entire DAG and assigning new ids after child ids
  // are assigned.
  DepthFirstWithAction({dag_root_node_id},
                       SubsplitDAGTraversalAction(
                           // BeforeNode
                           [](size_t node_id) {},
                           // AfterNode
                           [&node_reindexer, &running_traversal_idx](size_t node_id) {
                             node_reindexer.at(node_id) = running_traversal_idx;
                             running_traversal_idx++;
                           },
                           // BeforeNodeClade
                           [](size_t node_id, bool rotated) {},
                           // VisitEdge
                           [](size_t node_id, size_t child_id, bool rotated) {}));
  return node_reindexer;
}

SizeVector SubsplitDAG::BuildEdgeReindexer(const size_t prev_edge_count) {
  SizeVector edge_reindexer =
      Reindexer::IdentityReindexer(GPCSPCountWithFakeSubsplits());
  // Only edges from an existing parent node to a new child node need to be reindexed.
  // See SubsplitDAG::AddNodePair().
  for (size_t edge_idx = prev_edge_count; edge_idx < GPCSPCountWithFakeSubsplits();
       edge_idx++) {
    const auto &element =
        std::find_if(dag_edges_.begin(), dag_edges_.end(),
                     [edge_idx](auto pair) { return pair.second == edge_idx; });
    Assert(element != dag_edges_.end(),
           "An edge with given edge_idx did not exist in "
           "SubsplitDAG::BuildEdgeReindexer.");
    const auto &[node_pair, idx] = *element;
    std::ignore = idx;
    const auto &[parent_id, child_id] = node_pair;
    const Bitset parent_subsplit = GetDAGNode(parent_id)->GetBitset();
    const Bitset child_subsplit = GetDAGNode(child_id)->GetBitset();
    const auto idx_range = GetEdgeRange(
        parent_subsplit, child_subsplit.SubsplitIsRotatedChildOf(parent_subsplit));
    // New edge is added to the end of the range.
    const size_t new_idx = edge_reindexer.at(idx_range.second);
    Reindexer::ReassignAndShift(edge_reindexer, edge_idx, new_idx);
  }
  return edge_reindexer;
}

void SubsplitDAG::RemapNodeIds(const SizeVector &node_reindexer) {
  // Update `dag_nodes_`.
  auto dag_nodes_copy = Reindexer::Reindex(dag_nodes_, node_reindexer);
  dag_nodes_.swap(dag_nodes_copy);
  // Update each node's id and leafward/rootward ids.
  for (size_t node_id = 0; node_id < NodeCount(); node_id++) {
    GetDAGNode(node_id)->RemapNodeIds(node_reindexer);
  }
  // Update `subsplit_to_id_`.
  for (const auto &[subsplit, node_id] : subsplit_to_id_) {
    subsplit_to_id_.at(subsplit) = node_reindexer.at(node_id);
  }
  // Update `dag_edges_`.
  std::map<SizePair, size_t> dag_edges_copy;
  for (const auto &[node_pair, edge_idx] : dag_edges_) {
    SafeInsert(
        dag_edges_copy,
        {node_reindexer.at(node_pair.first), node_reindexer.at(node_pair.second)},
        edge_idx);
  }
  dag_edges_.swap(dag_edges_copy);
}

void SubsplitDAG::RemapEdgeIdxs(const SizeVector &edge_reindexer) {
  // Update `dag_edges_`.
  for (const auto &[node_pair, edge_idx] : dag_edges_) {
    dag_edges_.at(node_pair) = edge_reindexer.at(edge_idx);
  }
  // Update `parent_to_range_`.
  for (const auto &[subsplit, idx_range] : parent_to_range_) {
    parent_to_range_.at(subsplit) = {edge_reindexer.at(idx_range.first),
                                     edge_reindexer.at(idx_range.second)};
  }
}
