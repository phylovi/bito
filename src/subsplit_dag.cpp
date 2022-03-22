// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "subsplit_dag.hpp"

#include "combinatorics.hpp"
#include "numerical_utils.hpp"
#include "sbn_probability.hpp"

// ** Constructor methods:

SubsplitDAG::SubsplitDAG()
    : taxon_count_(0), edge_count_without_leaf_subsplits_(0), topology_count_(0.) {}

SubsplitDAG::SubsplitDAG(const RootedTreeCollection &tree_collection)
    : SubsplitDAG(tree_collection.TaxonCount(), tree_collection.TopologyCounter(),
                  tree_collection.TagTaxonMap()) {}

SubsplitDAG::SubsplitDAG(size_t taxon_count,
                         const Node::TopologyCounter &topology_counter,
                         const TagStringMap &tag_taxon_map)
    : dag_taxa_(), taxon_count_(taxon_count) {
  Assert(topology_counter.size() > 0, "Empty topology counter given to SubsplitDAG.");
  Assert(topology_counter.begin()->first->LeafCount() == taxon_count,
         "Taxon count mismatch in SubsplitDAG constructor.");
  auto [edge_indexer, index_to_child, rootsplits] =
      ProcessTopologyCounter(topology_counter);
  BuildTaxonMap(tag_taxon_map);
  BuildNodes(index_to_child, rootsplits);
  BuildEdges(index_to_child);
  BuildDAGEdgesFromEdgeIndexer(edge_indexer);
  AddLeafSubsplitsToDAGEdgesAndParentToRange();
  StoreEdgeIds();
  CountTopologies();

  Assert(IsValidTaxonMap(),
         "SubsplitDAG::SubsplitDAG(): Given Taxon Map does not represent a valid map. "
         "Taxon ids do not cover all indices from 0 to taxon_count.");
}

// ** Comparator

int SubsplitDAG::Compare(const SubsplitDAG &other) {
  return SubsplitDAG::Compare(*this, other);
}

int SubsplitDAG::Compare(const SubsplitDAG &lhs, const SubsplitDAG &rhs) {
  // (1) Compare Taxon Sizes.
  int taxon_diff = lhs.TaxonCount() - rhs.TaxonCount();
  if (taxon_diff != 0) {
    // #350 let's talk about -100 vs -200 here.
    return taxon_diff;
  }
  // Create translation map (lhs->rhs) for bitset clades.
  auto taxon_map = SubsplitDAG::BuildTaxonTranslationMap(lhs, rhs);
  // (2) Compare Subsplit Nodes.
  auto lhs_nodes = lhs.GetSortedVectorOfNodeBitsets();
  auto rhs_nodes = rhs.GetSortedVectorOfNodeBitsets();
  // Translate to account for different Taxon mappings and sort output.
  for (size_t i = 0; i < lhs_nodes.size(); i++) {
    lhs_nodes[i] =
        SubsplitDAG::BitsetTranslateViaTaxonTranslationMap(lhs_nodes[i], taxon_map);
    lhs_nodes[i] = lhs_nodes[i].SubsplitSortClades();
  }
  std::sort(lhs_nodes.begin(), lhs_nodes.end());
  if (lhs_nodes != rhs_nodes) {
    return (lhs_nodes < rhs_nodes) ? -1 : 1;
  }
  // (3) Compare PCSP Edges.
  auto lhs_edges = lhs.GetSortedVectorOfEdgeBitsets();
  auto rhs_edges = rhs.GetSortedVectorOfEdgeBitsets();
  // Translate to account for different Taxon mappings and sort output.
  for (size_t i = 0; i < lhs_edges.size(); i++) {
    lhs_edges[i] =
        SubsplitDAG::BitsetTranslateViaTaxonTranslationMap(lhs_edges[i], taxon_map);
    lhs_edges[i] = lhs_edges[i].PCSPSortClades();
  }
  std::sort(lhs_edges.begin(), lhs_edges.end());
  if (lhs_edges != rhs_edges) {
    return (lhs_edges < rhs_edges) ? -1 : 1;
  }
  return 0;
}

bool operator==(const SubsplitDAG &lhs, const SubsplitDAG &rhs) {
  return (SubsplitDAG::Compare(lhs, rhs) == 0);
}

bool operator!=(const SubsplitDAG &lhs, const SubsplitDAG &rhs) {
  return (SubsplitDAG::Compare(lhs, rhs) != 0);
}

// ** Count methods:

void SubsplitDAG::CountTopologies() {
  topology_count_below_ = EigenVectorXd::Ones(NodeCount());

  for (const auto &node_id : RootwardEdgeTraversalTrace(true)) {
    const auto &node = GetDAGNode(node_id);
    for (const bool rotated : {true, false}) {
      // When there are no leafward nodes in the `rotated` direction, we set the number
      // of topologies for the rotation of the node to be 1.
      double per_rotated_count = node.GetLeafward(rotated).empty() ? 1. : 0.;
      // Sum options across the possible children.
      for (const auto &child_id : node.GetLeafward(rotated)) {
        per_rotated_count += topology_count_below_[child_id];
      }
      // Take the product across the number of options for the left and right branches
      // of the tree.
      topology_count_below_[node_id] *= per_rotated_count;
    }
  }
  topology_count_ = topology_count_below_[GetDAGRootNodeId()];
}

size_t SubsplitDAG::TaxonCount() const { return dag_taxa_.size(); }

size_t SubsplitDAG::NodeCount() const { return storage_.GetVertices().size(); }

size_t SubsplitDAG::NodeCountWithoutDAGRoot() const { return NodeCount() - 1; }

double SubsplitDAG::TopologyCount() const { return topology_count_; }

size_t SubsplitDAG::RootsplitCount() const { return GetRootsplitNodeIds().size(); }

size_t SubsplitDAG::EdgeCount() const { return edge_count_without_leaf_subsplits_; }

size_t SubsplitDAG::EdgeCountWithLeafSubsplits() const {
  return storage_.GetLines().size();
}

// ** Output methods:

StringSizeMap SubsplitDAG::SummaryStatistics() const {
  return {{"node_count", NodeCount()}, {"edge_count", EdgeCountWithLeafSubsplits()}};
}

void SubsplitDAG::Print() const {
  for (auto dag_node : storage_.GetVertices()) {
    std::cout << dag_node.ToString() << std::endl;
  }
}

void SubsplitDAG::PrintNodes() const {
  for (const auto &dag_node : storage_.GetVertices()) {
    std::cout << dag_node.Id() << ": " << dag_node.GetBitset().SubsplitToString()
              << std::endl;
  }
}

void SubsplitDAG::PrintEdgeIndexer() const {
  for (const auto &[edge, idx] : BuildEdgeIndexer()) {
    std::cout << idx << ": " << edge.PCSPToString() << std::endl;
  }
}

void SubsplitDAG::PrintDAGEdges() const {
  for (const auto &[parent_child_id, edge_idx] : storage_.GetLines()) {
    const auto &[parent_id, child_id] = parent_child_id;
    std::cout << "[" << parent_id << ", " << child_id << "], " << edge_idx << std::endl;
  }
}

void SubsplitDAG::PrintParentToRange() const {
  for (const auto &[subsplit, edge_range] : parent_to_child_range_) {
    std::cout << subsplit.SubsplitToString() << ": [" << edge_range.first << ", "
              << edge_range.second << "]" << std::endl;
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
      {GetDAGRootNodeId()},
      SubsplitDAGTraversalAction(
          // BeforeNode
          [this, &string_stream, &show_index_labels](size_t node_id) {
            const auto &node = GetDAGNode(node_id);
            if (node.IsDAGRootNode()) {
              string_stream << node_id << " [label=\"<f0>&rho;\"]\n";
              return;
            }
            auto bs = node.GetBitset();
            string_stream << node_id << " [label=\"<f0>"
                          << bs.SubsplitGetClade(Bitset::SubsplitClade::Left)
                                 .ToVectorOfSetBitsAsString()
                          << "|<f1>";
            if (show_index_labels) {
              string_stream << node_id;
            }
            string_stream << "|<f2>"
                          << bs.SubsplitGetClade(Bitset::SubsplitClade::Right)
                                 .ToVectorOfSetBitsAsString()
                          << "\"]\n";
          },
          // AfterNode
          [](size_t node_id) {},
          // BeforeNodeClade
          [](size_t node_id, bool rotated) {},
          // VisitEdge
          [this, &string_stream, &show_index_labels](size_t node_id, size_t child_id,
                                                     bool rotated) {
            if (GetDAGNode(child_id).IsLeaf()) {
              string_stream << child_id << " [label=\"<f1>" << child_id << "\"]\n";
            }
            string_stream << "\"" << node_id << "\":";
            string_stream << (rotated ? "f0" : "f2");
            string_stream << "->\"";
            string_stream << child_id << "\":f1";
            if (show_index_labels) {
              string_stream << " [label=\"" << GetEdgeIdx(node_id, child_id);
              if (rotated) {
                string_stream << "\", color=1, fontcolor=1";
              } else {
                string_stream << "\", color=3, fontcolor=3";
              }
              if (GetDAGNode(node_id).IsDAGRootNode()) {
                string_stream << ",style=dashed]";
              } else {
                string_stream << "]";
              }
            } else {
              if (GetDAGNode(node_id).IsDAGRootNode()) {
                string_stream << "[style=dashed]";
              }
            }
            string_stream << "\n";
          }));
  string_stream << "}";
  return string_stream.str();
}

// ** Build output indexes/vectors methods:

BitsetSizeMap SubsplitDAG::BuildEdgeIndexer() const {
  auto edge_indexer = BitsetSizeMap();
  TopologicalEdgeTraversal([this, &edge_indexer](size_t parent_id, bool rotated,
                                                 size_t child_id, size_t edge_idx) {
    const auto parent_subsplit = GetDAGNode(parent_id).GetBitset(rotated);
    const auto child_subsplit = GetDAGNode(child_id).GetBitset();
    SafeInsert(edge_indexer, Bitset::PCSP(parent_subsplit, child_subsplit), edge_idx);
  });
  return edge_indexer;
}

SizeBoolVectorMap SubsplitDAG::BuildEdgeIdxToPCSPMap() const {
  auto edge_to_pcsp_map = SizeBoolVectorMap();
  TopologicalEdgeTraversal([this, &edge_to_pcsp_map](size_t parent_id, bool rotated,
                                                     size_t child_id, size_t edge_idx) {
    const auto parent_subsplit = GetDAGNode(parent_id).GetBitset(rotated);
    const auto child_subsplit = GetDAGNode(child_id).GetBitset();
    SafeInsert(edge_to_pcsp_map, edge_idx,
               Bitset::PCSP(parent_subsplit, child_subsplit).GetData());
  });
  return edge_to_pcsp_map;
}

size_t SubsplitDAG::GetTaxonId(const std::string &name) const {
  Assert(
      ContainsTaxon(name),
      "SubsplitDAG::GetTaxonId(): Taxon with given name is not contained in the DAG.");
  return dag_taxa_.find(name)->second;
}

SubsplitDAGNode SubsplitDAG::GetDAGNode(const size_t node_id) const {
  Assert(ContainsNode(node_id), "Node with the given node_id does not exist in DAG.");
  return storage_.GetVertices().at(node_id);
}

MutableSubsplitDAGNode SubsplitDAG::GetDAGNode(const size_t node_id) {
  Assert(ContainsNode(node_id), "Node with the given node_id does not exist in DAG.");
  return storage_.GetVertices().at(node_id);
}

size_t SubsplitDAG::GetDAGNodeId(const Bitset &subsplit) const {
  Assert(ContainsNode(subsplit), "Node with the given subsplit does not exist in DAG.");
  return subsplit_to_id_.at(subsplit);
}

size_t SubsplitDAG::GetDAGRootNodeId() const { return NodeCount() - 1; }

ConstNeighborsView SubsplitDAG::GetRootsplitNodeIds() const {
  return GetDAGNode(GetDAGRootNodeId()).GetLeftLeafward();
}

ConstLineView SubsplitDAG::GetDAGEdge(const size_t edge_id) const {
  Assert(ContainsEdge(edge_id), "Node with the given node_id does not exist in DAG.");
  return *storage_.GetLine(edge_id);
}

size_t SubsplitDAG::GetEdgeIdx(const Bitset &parent_subsplit,
                               const Bitset &child_subsplit) const {
  return GetEdgeIdx(GetDAGNodeId(parent_subsplit), GetDAGNodeId(child_subsplit));
}

size_t SubsplitDAG::GetEdgeIdx(size_t parent_id, size_t child_id) const {
  auto edge = storage_.GetLine(parent_id, child_id);
  if (!edge.has_value()) Failwith("Edge not found");
  return edge.value().GetId();
}

SizePair SubsplitDAG::GetChildEdgeRange(const Bitset &subsplit,
                                        const bool rotated) const {
  Assert(
      ContainsNode(subsplit),
      "Node with the given subsplit does not exist in SubsplitDAG::GetChildEdgeRange.");
  return parent_to_child_range_.at(SubsplitToSortedOrder(subsplit, rotated));
}

std::vector<std::string> SubsplitDAG::GetSortedVectorOfTaxonNames() const {
  std::vector<std::string> taxa;
  for (const auto &name_id : dag_taxa_) {
    taxa.push_back(name_id.first);
  }
  std::sort(taxa.begin(), taxa.end());
  return taxa;
}

std::vector<Bitset> SubsplitDAG::GetSortedVectorOfNodeBitsets() const {
  std::vector<Bitset> nodes;
  for (size_t i = 0; i < NodeCount(); i++) {
    Bitset node_bitset = GetDAGNode(i).GetBitset();
    nodes.push_back(node_bitset);
  }
  std::sort(nodes.begin(), nodes.end());
  return nodes;
}

std::vector<Bitset> SubsplitDAG::GetSortedVectorOfEdgeBitsets() const {
  std::vector<Bitset> edges;
  for (auto i : storage_.GetLines()) {
    auto parent_bitset = GetDAGNode(i.GetParent()).GetBitset();
    auto child_bitset = GetDAGNode(i.GetChild()).GetBitset();
    Bitset edge_bitset = Bitset::PCSP(parent_bitset, child_bitset);
    edges.push_back(edge_bitset);
  }
  std::sort(edges.begin(), edges.end());
  return edges;
}

const std::map<std::string, size_t> &SubsplitDAG::GetTaxonMap() const {
  return dag_taxa_;
}

const BitsetSizeMap &SubsplitDAG::GetSubsplitToIdMap() const { return subsplit_to_id_; }

EigenVectorXd SubsplitDAG::BuildUniformOnTopologicalSupportPrior() const {
  EigenVectorXd q = EigenVectorXd::Ones(EdgeCountWithLeafSubsplits());

  for (const auto &node_id : RootwardEdgeTraversalTrace(true)) {
    const auto &node = GetDAGNode(node_id);
    for (const bool rotated : {false, true}) {
      if (!node.GetLeafward(rotated).empty()) {
        double per_rotated_count = 0.;
        for (auto child_id : node.GetLeafward(rotated)) {
          per_rotated_count += topology_count_below_[child_id];
        }
        for (auto child_id : node.GetLeafward(rotated)) {
          size_t edge_idx = GetEdgeIdx(node.Id(), child_id);
          q(edge_idx) = topology_count_below_(child_id) / per_rotated_count;
        }
      }
    }
  }
  return q;
}

Node::NodePtrVec SubsplitDAG::GenerateAllTopologies() const {
  std::vector<Node::NodePtrVec> topology_below(NodeCount());

  auto GetSubtopologies = [&topology_below](SubsplitDAGNode node) {
    Node::NodePtrVec rotated_subtopologies, sorted_subtopologies;
    for (const bool rotated : {false, true}) {
      for (const auto &child_id : node.GetLeafward(rotated)) {
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
    if (node_id == GetDAGRootNodeId()) {
      // DAG root node has no `sorted_subtopologies`, so loop above yields empty
      // `topologies` vector.
      return rotated_subtopologies;
    }
    return topologies;
  };

  for (const auto &node_id : RootwardEdgeTraversalTrace(true)) {
    const auto &node = GetDAGNode(node_id);
    if (node.IsLeaf()) {
      topology_below.at(node_id).push_back(Node::Leaf(node_id));
    } else {
      auto [rotated_topologies, sorted_topologies] = GetSubtopologies(node);
      topology_below[node_id] =
          MergeTopologies(node_id, rotated_topologies, sorted_topologies);
    }
  }

  const auto &topologies = topology_below.at(GetDAGRootNodeId());
  Assert(topologies.size() == TopologyCount(),
         "The realized number of topologies does not match the expected count.");

  // We return a deep copy of every Polished topology to avoid loops in the pointer
  // structure. Such loops can create problems when we Polish the topologies one at a
  // time: polishing a second topology can change the numbering of a} previous topology.
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
  EigenVectorXd result = EigenVectorXd::Zero(EdgeCountWithLeafSubsplits());
  for (const auto &[parent_child_id, edge_idx] : storage_.GetLines()) {
    const auto &[parent_id, child_id] = parent_child_id;
    std::ignore = parent_id;
    // If child is a leaf and subsplit is sorted, then child0 will have a zero taxon
    // count.
    size_t child0_taxon_count = GetDAGNode(child_id)
                                    .GetBitset()
                                    .SubsplitGetClade(Bitset::SubsplitClade::Left)
                                    .Count();
    // As long as subsplit is sorted and nonempty, then child1 will have a nonzero taxon
    // count.
    size_t child1_taxon_count = GetDAGNode(child_id)
                                    .GetBitset()
                                    .SubsplitGetClade(Bitset::SubsplitClade::Right)
                                    .Count();
    // The ordering of this subsplit is flipped so that this ratio will be nonzero in
    // the denominator in the case of root & leaves.
    result(edge_idx) = Combinatorics::LogChildSubsplitCountRatio(child1_taxon_count,
                                                                 child0_taxon_count);
  }
  NumericalUtils::Exponentiate(result);
  return result;
}

// ** DAG Lambda Iterator methods:

void SubsplitDAG::IterateOverRealNodes(const NodeLambda &f) const {
  Assert(taxon_count_ < NodeCount(), "No real DAG nodes!");
  for (auto it = storage_.GetVertices().cbegin() + taxon_count_;
       it < storage_.GetVertices().cend() - 1; it++) {
    f(*it);
  }
}

void SubsplitDAG::IterateOverLeafwardEdges(SubsplitDAGNode node, bool rotated,
                                           const NodeLambda &f) const {
  for (const size_t child_id : node.GetLeafward(rotated)) {
    f(GetDAGNode(child_id));
  }
}

void SubsplitDAG::IterateOverLeafwardEdges(SubsplitDAGNode node,
                                           const EdgeDestinationLambda &f) const {
  for (bool rotated : {false, true}) {
    for (const size_t child_id : node.GetLeafward(rotated)) {
      f(rotated, GetDAGNode(child_id));
    }
  }
}

void SubsplitDAG::IterateOverLeafwardEdgesAndChildren(
    SubsplitDAGNode node, const EdgeAndNodeLambda &f) const {
  IterateOverLeafwardEdges(node,
                           [this, &node, &f](bool rotated, SubsplitDAGNode child) {
                             f(GetEdgeIdx(node.Id(), child.Id()), rotated, child.Id());
                           });
}

void SubsplitDAG::IterateOverRootwardEdges(SubsplitDAGNode node,
                                           const EdgeDestinationLambda &f) const {
  if (!node.IsRootsplit()) {
    for (bool rotated : {false, true}) {
      for (const size_t parent_id : node.GetRootward(rotated)) {
        f(rotated, GetDAGNode(parent_id));
      }
    }
  }
}

void SubsplitDAG::IterateOverRootwardEdgesAndParents(SubsplitDAGNode node,
                                                     const EdgeAndNodeLambda &f) const {
  IterateOverRootwardEdges(
      node, [this, &node, &f](bool rotated, SubsplitDAGNode parent) {
        f(GetEdgeIdx(parent.Id(), node.Id()), rotated, parent.Id());
      });
}

void SubsplitDAG::IterateOverParentAndChildAndLeafwardEdges(
    SubsplitDAGNode node, const ParentRotationChildEdgeLambda &f) const {
  IterateOverLeafwardEdges(
      node, [this, &node, &f](bool rotated, SubsplitDAGNode child) {
        f(node.Id(), rotated, child.Id(), GetEdgeIdx(node.Id(), child.Id()));
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
  node_probabilities[GetDAGRootNodeId()] = 1.;

  TopologicalEdgeTraversal([&node_probabilities, &normalized_sbn_parameters](
                               const size_t parent_id, const bool,
                               const size_t child_id, const size_t edge_idx) {
    const double child_probability_given_parent = normalized_sbn_parameters[edge_idx];
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
    const auto &subsplit_bitset = GetDAGNode(node_id).GetBitset();
    if (node_id != GetDAGRootNodeId() && !subsplit_bitset.SubsplitIsLeaf()) {
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
  TopologicalEdgeTraversal([this, &node_probabilities, &normalized_sbn_parameters,
                            &inverted_probabilities](const size_t parent_id, const bool,
                                                     const size_t child_id,
                                                     const size_t edge_idx) {
    // The traversal doesn't set the rootsplit probabilities, but those are always 1
    // (there is only one "parent" of a rootsplit).
    if (parent_id != GetDAGRootNodeId()) {
      // For a PCSP t -> s:
      inverted_probabilities[edge_idx] =         // P(t|s)
          node_probabilities[parent_id] *        // P(t)
          normalized_sbn_parameters[edge_idx] /  // P(s|t)
          node_probabilities[child_id];          // P(s)
    }
  });
  return inverted_probabilities;
}

std::vector<Bitset> SubsplitDAG::GetChildSubsplits(const SizeBitsetMap &index_to_child,
                                                   const Bitset &parent_subsplit,
                                                   bool include_leaf_subsplits) {
  std::vector<Bitset> children_subsplits;
  // Add all child subsplit bitsets to vector.
  if (parent_to_child_range_.count(parent_subsplit) > 0) {
    const auto [start, stop] = parent_to_child_range_.at(parent_subsplit);
    for (auto idx = start; idx < stop; idx++) {
      children_subsplits.push_back(index_to_child.at(idx));
    }
  } else if (include_leaf_subsplits) {
    // This method is designed to be called before calling
    // AddLeafSubsplitsToDAGEdgesAndParentToRange. In that case, if the second clade
    // of the subsplit is just a single taxon, the subsplit will not map to any value in
    // parent_to_child_range_.
    //
    // But we still need to create and connect to leaf subsplits in the DAG. So, here we
    // make a leaf child subsplit.
    children_subsplits.push_back(Bitset::LeafSubsplitOfParentSubsplit(parent_subsplit));
  }

  return children_subsplits;
}

// ** Modify DAG methods:

std::tuple<BitsetSizeMap, SizeBitsetMap, BitsetVector>
SubsplitDAG::ProcessTopologyCounter(const Node::TopologyCounter &topology_counter) {
  BitsetSizeMap edge_indexer;
  SizeBitsetMap index_to_child;
  BitsetVector rootsplits;
  std::tie(rootsplits, edge_indexer, index_to_child, parent_to_child_range_,
           edge_count_without_leaf_subsplits_) =
      SBNMaps::BuildIndexerBundle(RootedSBNMaps::RootsplitCounterOf(topology_counter),
                                  RootedSBNMaps::PCSPCounterOf(topology_counter));
  return {edge_indexer, index_to_child, rootsplits};
}

void SubsplitDAG::BuildTaxonMap(const TagStringMap &tag_taxon_map) {
  // Insert all taxa from tree_collections's map to SubsplitDAG map.
  for (const auto &[tag, name] : tag_taxon_map) {
    // The "tag" key of the tree_collection's taxon_map is 2 bitpacked ints: [id,
    // topology count]. We only care about the id.
    size_t id = static_cast<size_t>(UnpackFirstInt(tag));
    dag_taxa_.insert(std::make_pair(name, id));
  }
}

size_t SubsplitDAG::CreateAndInsertNode(const Bitset &subsplit) {
  size_t node_id = NodeCount();
  storage_.AddVertex({node_id, subsplit});
  SafeInsert(subsplit_to_id_, subsplit, node_id);
  return node_id;
}

size_t SubsplitDAG::CreateAndInsertEdge(const size_t parent_id, const size_t child_id,
                                        bool rotated) {
  Assert(ContainsNode(parent_id), "Node with the given parent_id does not exist.");
  Assert(ContainsNode(child_id), "Node with the given child_id does not exist.");
  // Insert edge between parent and child.
  size_t edge_idx = EdgeCountWithLeafSubsplits();
  ConnectGivenNodes(parent_id, child_id, rotated, edge_idx);
  storage_.AddLine(
      {edge_idx, parent_id, child_id, rotated ? Clade::Left : Clade::Right});
  return edge_idx;
}

void SubsplitDAG::ConnectGivenNodes(const size_t parent_id, const size_t child_id,
                                    bool rotated, size_t edge_id) {
  auto parent_node = GetDAGNode(parent_id);
  auto child_node = GetDAGNode(child_id);

  Clade which_parent_clade = (rotated ? Clade::Left : Clade::Right);
  parent_node.AddEdge(child_node.Id(), edge_id, Direction::Leafward,
                      which_parent_clade);
  child_node.AddEdge(parent_node.Id(), edge_id, Direction::Rootward,
                     which_parent_clade);
}

void SubsplitDAG::ConnectNodes(const SizeBitsetMap &index_to_child, size_t node_id,
                               bool rotated) {
  // #350 look for "rotation / rotated" in the code and fix.
  // Get bitset of parent node according to its rotation.
  const Bitset subsplit = GetDAGNode(node_id).GetBitset(rotated);
  // Build vector of child node's subsplits.
  const auto children = GetChildSubsplits(index_to_child, subsplit, true);
  // Connect parent node to all child nodes.
  for (const auto &child_subsplit : children) {
    ConnectGivenNodes(node_id, GetDAGNodeId(child_subsplit), rotated, NoId);
  }
}

void SubsplitDAG::BuildNodes(const SizeBitsetMap &index_to_child,
                             const BitsetVector &rootsplits) {
  std::unordered_set<Bitset> visited_subsplits;
  // We will create leaf subsplits and insert to storage_ nodes.
  // Leaf subsplits (i.e. leaf nodes) will take IDs in [0, taxon_count_).
  for (size_t taxon_idx = 0; taxon_idx < taxon_count_; taxon_idx++) {
    CreateAndInsertNode(Bitset::LeafSubsplitOfNonemptyClade(
        Bitset::Singleton(taxon_count_, taxon_idx)));
  }
  // Then we add the remaining nodes.
  // The root splits will take on the higher IDs compared to the non-rootsplits.
  for (const auto &rootsplit : rootsplits) {
    BuildNodesDepthFirst(index_to_child, rootsplit, visited_subsplits);
  }
  // #350 look for "DAG root" and fix
  // Finally, we add the DAG root node.
  CreateAndInsertNode(Bitset::UCASubsplitOfTaxonCount(taxon_count_));
}

void SubsplitDAG::BuildNodesDepthFirst(const SizeBitsetMap &index_to_child,
                                       const Bitset &subsplit,
                                       std::unordered_set<Bitset> &visited_subsplits) {
  visited_subsplits.insert(subsplit);
  for (bool rotated : {false, true}) {
    for (const auto &child_subsplit : GetChildSubsplits(
             index_to_child, SubsplitToSortedOrder(subsplit, rotated), false)) {
      if (visited_subsplits.count(child_subsplit) == 0) {
        BuildNodesDepthFirst(index_to_child, child_subsplit, visited_subsplits);
      }
    }
  }
  CreateAndInsertNode(subsplit);
}

void SubsplitDAG::BuildEdges(const SizeBitsetMap &index_to_child) {
  // Connect every node except for the DAG root node.
  for (size_t node_id = taxon_count_; node_id < GetDAGRootNodeId(); node_id++) {
    ConnectNodes(index_to_child, node_id, true);
    ConnectNodes(index_to_child, node_id, false);
  }
  // Connect the DAG root node.
  ConnectNodes(index_to_child, GetDAGRootNodeId(), true);
}

void SubsplitDAG::BuildDAGEdgesFromEdgeIndexer(BitsetSizeMap &edge_indexer) {
  for (const auto &[edge, index] : edge_indexer) {
    Assert(edge.size() == 3 * taxon_count_,
           "All edges should be bitsets with size 3 times taxon_count_.");
    const auto parent_id = GetDAGNodeId(edge.PCSPGetParentSubsplit());
    const auto child_id = GetDAGNodeId(edge.PCSPGetChildSubsplit());
    auto left = GetDAGNode(parent_id).GetLeftLeafward();
    auto clade = std::find(left.begin(), left.end(), child_id) != left.end()
                     ? Clade::Left
                     : Clade::Right;
    if (clade == Clade::Right) {
      auto right = GetDAGNode(parent_id).GetRightLeafward();
      Assert(std::find(right.begin(), right.end(), child_id) != right.end(),
             "Parent has no connection to child");
    }
    storage_.AddLine({index, parent_id, child_id, clade});
  }
}

void SubsplitDAG::AddLeafSubsplitsToDAGEdgesAndParentToRange() {
  for (size_t node_id = 0; node_id < taxon_count_; node_id++) {
    const auto current_bitset(GetDAGNode(node_id).GetBitset());
    IterateOverRootwardEdges(
        GetDAGNode(node_id),
        [this, current_bitset](const bool rotated, SubsplitDAGNode node) {
          SafeInsert(parent_to_child_range_, node.GetBitset(rotated),
                     {EdgeCountWithLeafSubsplits(), EdgeCountWithLeafSubsplits() + 1});
          storage_.AddLine({EdgeCountWithLeafSubsplits(), node.Id(),
                            GetDAGNodeId(current_bitset),
                            rotated ? Clade::Left : Clade::Right});
        });
  }
}

void SubsplitDAG::StoreEdgeIds() {
  for (auto edge : storage_.GetLines()) {
    auto parent = storage_.GetVertices().at(edge.GetParent());
    auto child = storage_.GetVertices().at(edge.GetChild());
    parent.SetLineId(edge.GetChild(), edge.GetId());
    child.SetLineId(edge.GetParent(), edge.GetId());
  }
}

void SubsplitDAG::RootwardDepthFirst(size_t node_id, SizeVector &visit_order,
                                     std::unordered_set<size_t> &visited_nodes) const {
  // Add to set of all visited nodes.
  SafeInsert(visited_nodes, node_id);
  // #350: Look for sorted/rotated and update.
  // Recurse on sorted children.
  const auto &node = GetDAGNode(node_id);
  for (size_t parent_id : node.GetRightRootward()) {
    if (visited_nodes.count(parent_id) == 0) {
      RootwardDepthFirst(parent_id, visit_order, visited_nodes);
    }
  }
  // Recurse on rotated children.
  for (size_t parent_id : node.GetLeftRootward()) {
    if (visited_nodes.count(parent_id) == 0) {
      RootwardDepthFirst(parent_id, visit_order, visited_nodes);
    }
  }
  // Append to vector post-order (after all children have been visited).
  visit_order.push_back(node_id);
}

void SubsplitDAG::LeafwardDepthFirst(size_t node_id, SizeVector &visit_order,
                                     std::unordered_set<size_t> &visited_nodes) const {
  // Add to set of all visited nodes.
  SafeInsert(visited_nodes, node_id);
  // Recurse on right/sorted children.
  for (size_t child_id : GetDAGNode(node_id).GetRightLeafward()) {
    if (visited_nodes.count(child_id) == 0) {
      LeafwardDepthFirst(child_id, visit_order, visited_nodes);
    }
  }
  // Recurse on left/rotated children.
  for (size_t child_id : GetDAGNode(node_id).GetLeftLeafward()) {
    if (visited_nodes.count(child_id) == 0) {
      LeafwardDepthFirst(child_id, visit_order, visited_nodes);
    }
  }
  // Append to vector post-order (after all children have been visited).
  visit_order.push_back(node_id);
}

SizeVector SubsplitDAG::LeafwardEdgeTraversalTrace(bool include_dag_root_node) const {
  SizeVector visit_order;
  std::unordered_set<size_t> visited_nodes;
  if (!include_dag_root_node) {
    SafeInsert(visited_nodes, GetDAGRootNodeId());
  }
  for (size_t leaf_id = 0; leaf_id < taxon_count_; leaf_id++) {
    RootwardDepthFirst(leaf_id, visit_order, visited_nodes);
  }
  return visit_order;
}

SizeVector SubsplitDAG::RootwardEdgeTraversalTrace(bool include_dag_root_node) const {
  SizeVector visit_order;
  std::unordered_set<size_t> visited_nodes;
  for (const auto &rootsplit_id : GetRootsplitNodeIds()) {
    LeafwardDepthFirst(rootsplit_id, visit_order, visited_nodes);
  }
  if (include_dag_root_node) {
    visit_order.push_back(GetDAGRootNodeId());
  }
  return visit_order;
}

SizeVector SubsplitDAG::TopologicalEdgeTraversalTrace() const {
  auto visit_order = RootwardEdgeTraversalTrace(true);
  std::reverse(visit_order.begin(), visit_order.end());
  return visit_order;
}

void SubsplitDAG::TopologicalEdgeTraversal(ParentRotationChildEdgeLambda f) const {
  for (const auto node_id : TopologicalEdgeTraversalTrace()) {
    IterateOverLeafwardEdgesAndChildren(
        GetDAGNode(node_id), [&f, &node_id](const size_t edge_idx, const bool rotated,
                                            const size_t child_id) {
          f(node_id, rotated, child_id, edge_idx);
        });
  }
}

// ** Miscellaneous

Bitset SubsplitDAG::SubsplitToSortedOrder(const Bitset &subsplit, bool rotated) const {
  return rotated ? subsplit.SubsplitRotate() : subsplit;
}

SizeVector SubsplitDAG::BuildTaxonTranslationMap(const SubsplitDAG &dag_a,
                                                 const SubsplitDAG &dag_b) {
  auto names_a = dag_a.GetSortedVectorOfTaxonNames();
  auto names_b = dag_b.GetSortedVectorOfTaxonNames();
  Assert(names_a == names_b,
         "SubsplitDAG::BuildTaxonTranslationMap(): SubsplitDAGs do not cover the same "
         "taxon set.");
  Assert(names_a.size() == dag_a.TaxonCount(),
         "SubsplitDAG::BuildTaxonTranslationMap(): Number of taxon names does not "
         "match the number of taxa in the DAG.");

  SizeVector taxon_map(names_a.size());
  for (const auto &name : names_a) {
    taxon_map[dag_a.GetTaxonId(name)] = dag_b.GetTaxonId(name);
  }
  return taxon_map;
}

int SubsplitDAG::BitsetCompareViaTaxonTranslationMap(const Bitset &bitset_a,
                                                     const Bitset &bitset_b,
                                                     const SizeVector &taxon_map) {
  Bitset bitset_a_translated_to_b =
      BitsetTranslateViaTaxonTranslationMap(bitset_a, taxon_map);
  return Bitset::Compare(bitset_a_translated_to_b, bitset_b);
}

Bitset SubsplitDAG::BitsetTranslateViaTaxonTranslationMap(
    const Bitset &bitset, const SizeVector &taxon_map, const bool forward_translate) {
  Bitset translated_bitset(bitset.size());
  size_t clade_size = taxon_map.size();
  size_t clade_count = bitset.size() / taxon_map.size();
  // Remap each clade individually.
  for (size_t i = 0; i < clade_count; i++) {
    // Remap each bit in clade according to a new position according to the taxon map.
    for (size_t j = 0; j < clade_size; j++) {
      size_t input_offset = (i * clade_size) + j;
      size_t output_offset = (i * clade_size) + taxon_map[j];
      if (forward_translate) {
        translated_bitset.set(input_offset, bitset[output_offset]);
      } else {
        translated_bitset.set(output_offset, bitset[input_offset]);
      }
    }
  }
  return translated_bitset;
}

TagStringMap SubsplitDAG::BuildDummyTagTaxonMap(const size_t taxon_count) {
  TagStringMap tag_taxon_map;
  for (size_t i = 0; i < taxon_count; i++) {
    std::string name = std::string("x" + std::to_string(i));
    tag_taxon_map.insert(std::make_pair(PackInts(i, 0), name));
  }
  return tag_taxon_map;
}

// ** Contains

bool SubsplitDAG::ContainsTaxon(const std::string &name) const {
  return dag_taxa_.find(name) != dag_taxa_.end();
}

bool SubsplitDAG::ContainsNode(const Bitset &subsplit) const {
  return subsplit_to_id_.find(subsplit) != subsplit_to_id_.end();
}

bool SubsplitDAG::ContainsNode(const size_t node_id) const {
  return node_id < NodeCount();
}

bool SubsplitDAG::ContainsEdge(const size_t parent_id, const size_t child_id) const {
  return storage_.GetLine(parent_id, child_id).has_value();
}

bool SubsplitDAG::ContainsEdge(const size_t edge_id) const {
  return storage_.GetLine(edge_id).has_value();
}

// ** Build output indexes/vectors

std::pair<SizeVector, SizeVector> SubsplitDAG::BuildParentIdVectors(
    const Bitset &subsplit) const {
  SizeVector left_parents, right_parents;
  // Linear search for all parents.
  for (const auto &[potential_parent_subsplit, node_id] : subsplit_to_id_) {
    if (subsplit.SubsplitIsLeftChildOf(potential_parent_subsplit)) {
      left_parents.push_back(node_id);
    } else if (subsplit.SubsplitIsRightChildOf(potential_parent_subsplit)) {
      right_parents.push_back(node_id);
    }
  }
  return {left_parents, right_parents};
}

std::pair<SizeVector, SizeVector> SubsplitDAG::BuildChildIdVectors(
    const Bitset &subsplit) const {
  SizeVector left_children, right_children;
  // Linear search for all children.
  for (const auto &[potential_child_subsplit, node_id] : subsplit_to_id_) {
    if (potential_child_subsplit.SubsplitIsLeftChildOf(subsplit)) {
      left_children.push_back(node_id);
    } else if (potential_child_subsplit.SubsplitIsRightChildOf(subsplit)) {
      right_children.push_back(node_id);
    }
  }
  return {left_children, right_children};
}

void SubsplitDAG::ConnectChildToAllChildren(const Bitset &child_subsplit,
                                            SizeVector &added_edge_idxs) {
  const auto [left_leafward_of_child, right_leafward_of_child] =
      BuildChildIdVectors(child_subsplit);

  for (const auto &[children_of_child, rotated] :
       std::vector<std::pair<SizeVector, bool>>{{left_leafward_of_child, true},
                                                {right_leafward_of_child, false}}) {
    SafeInsert(parent_to_child_range_, SubsplitToSortedOrder(child_subsplit, rotated),
               {EdgeCountWithLeafSubsplits(),
                EdgeCountWithLeafSubsplits() + children_of_child.size()});

    for (const size_t child_of_child_id : children_of_child) {
      const auto new_edge_idx =
          CreateAndInsertEdge(GetDAGNodeId(child_subsplit), child_of_child_id, rotated);
      added_edge_idxs.push_back(new_edge_idx);
    }
  }
}

void SubsplitDAG::ConnectParentToAllChildrenExcept(const Bitset &parent_subsplit,
                                                   const Bitset &child_subsplit,
                                                   SizeVector &added_edge_idxs) {
  const auto [left_leafward_of_parent, right_leafward_of_parent] =
      BuildChildIdVectors(parent_subsplit);

  for (const auto &[children_of_parent, rotated] :
       std::vector<std::pair<SizeVector, bool>>{{left_leafward_of_parent, true},
                                                {right_leafward_of_parent, false}}) {
    SafeInsert(parent_to_child_range_, SubsplitToSortedOrder(parent_subsplit, rotated),
               {EdgeCountWithLeafSubsplits(),
                EdgeCountWithLeafSubsplits() + children_of_parent.size()});

    for (const size_t child_of_parent_id : children_of_parent) {
      if (child_of_parent_id != GetDAGNodeId(child_subsplit)) {
        const auto new_edge_idx = CreateAndInsertEdge(GetDAGNodeId(parent_subsplit),
                                                      child_of_parent_id, rotated);
        added_edge_idxs.push_back(new_edge_idx);
      }
    }
  }
}

void SubsplitDAG::ConnectChildToAllParentsExcept(const Bitset &parent_subsplit,
                                                 const Bitset &child_subsplit,
                                                 SizeVector &added_edge_idxs) {
  const auto [left_rootward_of_child, right_rootward_of_child] =
      BuildParentIdVectors(child_subsplit);

  for (const auto &[parents_of_child, rotated] :
       std::vector<std::pair<SizeVector, bool>>{{left_rootward_of_child, true},
                                                {right_rootward_of_child, false}}) {
    for (const size_t parent_of_child_id : parents_of_child) {
      if (parent_of_child_id != GetDAGNodeId(parent_subsplit)) {
        const auto new_edge_idx = CreateAndInsertEdge(
            parent_of_child_id, GetDAGNodeId(child_subsplit), rotated);
        added_edge_idxs.push_back(new_edge_idx);
      }
    }
  }
}

void SubsplitDAG::ConnectParentToAllParents(const Bitset &parent_subsplit,
                                            SizeVector &added_edge_idxs) {
  const auto [left_rootward_of_parent, right_rootward_of_parent] =
      BuildParentIdVectors(parent_subsplit);

  for (const auto &[parents_of_parent, rotated] :
       std::vector<std::pair<SizeVector, bool>>{{left_rootward_of_parent, true},
                                                {right_rootward_of_parent, false}}) {
    for (const size_t parent_of_parent_id : parents_of_parent) {
      const auto new_edge_idx = CreateAndInsertEdge(
          parent_of_parent_id, GetDAGNodeId(parent_subsplit), rotated);
      added_edge_idxs.push_back(new_edge_idx);
    }
  }
}

SubsplitDAG::ModificationResult SubsplitDAG::AddNodePair(const NNIOperation &nni) {
  return AddNodePair(nni.parent_, nni.child_);
}

SubsplitDAG::ModificationResult SubsplitDAG::AddNodePair(const Bitset &parent_subsplit,
                                                         const Bitset &child_subsplit) {
  // Check that node pair will create a valid SubsplitDAG.
  Assert(
      IsValidAddNodePair(parent_subsplit, child_subsplit),
      "The given pair of nodes is incompatible with DAG in SubsplitDAG::AddNodePair.");
  // Initialize output vectors.
  SizeVector added_node_ids, added_edge_idxs, node_reindexer, edge_reindexer;
  // Check if either parent or child don't already exist in the DAG.
  const bool parent_is_new = !ContainsNode(parent_subsplit);
  const bool child_is_new = !ContainsNode(child_subsplit);
  // Soft assert: This allows for parent-child pair to exist in the DAG, but no work is
  // done. If both the parent and child already exist in DAG, return added_node_ids and
  // added_edge_idxs as empty, and node_reindexer and edge_reindexer as identity
  // reindexers.
  if (!parent_is_new && !child_is_new) {
    // Return default reindexers if both nodes already exist.
    node_reindexer = Reindexer::IdentityReindexer(NodeCount());
    edge_reindexer = Reindexer::IdentityReindexer(EdgeCountWithLeafSubsplits());
    return {added_node_ids, added_edge_idxs, node_reindexer, edge_reindexer};
  }

  // Note: `prev_node_count` acts as a place marker. We know what the DAG root node id
  // is
  // (`prev_node_count - 1`).
  const size_t prev_node_count = NodeCount();

  // Add parent/child nodes and connect them to their children
  // If child node is new, add node and connect it to all its children.
  if (child_is_new) {
    CreateAndInsertNode(child_subsplit);
    added_node_ids.push_back(GetDAGNodeId(child_subsplit));
    // Don't reindex these edges.
    ConnectChildToAllChildren(child_subsplit, added_edge_idxs);
  }
  // If parent node is new, add node it to all its children (except )
  if (parent_is_new) {
    CreateAndInsertNode(parent_subsplit);
    added_node_ids.push_back(GetDAGNodeId(parent_subsplit));
    // Don't reindex these edges.
    ConnectParentToAllChildrenExcept(parent_subsplit, child_subsplit, added_edge_idxs);
  }
  // Note: `prev_edge_count` is a marker conveying where we need to start
  // reindexing edge idxs.
  // Edges are only reindexed if the parent node already existed in the DAG
  // (so as to ensure that edges descending from the same node clade have
  // contiguous idxs).
  size_t prev_edge_count = EdgeCountWithLeafSubsplits();
  // Connect the given parent node to the given child node.
  added_edge_idxs.push_back(EdgeCountWithLeafSubsplits());
  CreateAndInsertEdge(GetDAGNodeId(parent_subsplit), GetDAGNodeId(child_subsplit),
                      child_subsplit.SubsplitIsLeftChildOf(parent_subsplit));
  // Don't reindex the edge between the given parent and child if the parent is new.
  if (parent_is_new) {
    prev_edge_count = EdgeCountWithLeafSubsplits();
  }
  if (child_is_new) {
    // Reindex these edges.
    ConnectChildToAllParentsExcept(parent_subsplit, child_subsplit, added_edge_idxs);
  }
  if (parent_is_new) {
    // Reindex these edges.
    ConnectParentToAllParents(parent_subsplit, added_edge_idxs);
  }
  // Create reindexers.
  node_reindexer = BuildNodeReindexer(prev_node_count);
  edge_reindexer = BuildEdgeReindexer(prev_edge_count);
  // Update the ids in added_node_ids and added_edge_idxs according to the reindexers.
  Reindexer::RemapIdVector(added_node_ids, node_reindexer);
  Reindexer::RemapIdVector(added_edge_idxs, edge_reindexer);
  // Update fields in the Subsplit DAG according to the reindexers.
  RemapNodeIds(node_reindexer);
  RemapEdgeIdxs(edge_reindexer);
  // Recount topologies to update `topology_count_below_`.
  CountTopologies();

  return {added_node_ids, added_edge_idxs, node_reindexer, edge_reindexer};
}

// ** Validation Test methods:

bool SubsplitDAG::IsConsistent() const {
  Failwith("SubsplitDAG::IsConsistent() is not yet implemented.");
  return false;
}

bool SubsplitDAG::IsValid() const {
  size_t correct_id = 0;
  for (auto node : storage_.GetVertices()) {
    if (correct_id++ != node.Id()) {
      return false;
    }
    if (!node.IsValid()) {
      return false;
    }
  }
  return true;
}

bool SubsplitDAG::IsValidAddNodePair(const Bitset &parent_subsplit,
                                     const Bitset &child_subsplit) const {
  // Get the number of adjacent nodes in the given direction.
  auto GetNodeCounts = [this](const Bitset &subsplit, const bool is_rootward) {
    const auto [left, right] =
        is_rootward ? BuildParentIdVectors(subsplit) : BuildChildIdVectors(subsplit);
    SizePair sides = {left.size(), right.size()};
    return sides;
  };
  // Get all adjacent nodes, not including parent and child.
  auto [left_rootward_of_parent, right_rootward_of_parent] =
      GetNodeCounts(parent_subsplit, true);
  auto [left_leafward_of_parent, right_leafward_of_parent] =
      GetNodeCounts(parent_subsplit, false);
  auto [left_leafward_of_child, right_leafward_of_child] =
      GetNodeCounts(child_subsplit, false);
  // Add child to parent's adjacent nodes.
  const bool is_left_child = child_subsplit.SubsplitIsLeftChildOf(parent_subsplit);
  if (is_left_child) {
    left_leafward_of_parent++;
  } else {
    right_leafward_of_parent++;
  }

  // (1) Added nodes are parent/child pair.
  if (Bitset::SubsplitIsParentChildPair(parent_subsplit, child_subsplit) == false) {
    return false;
  }
  // (2) Nodes do not add or remove taxa.
  if ((parent_subsplit.size() != 2 * taxon_count_) ||
      (child_subsplit.size() != 2 * taxon_count_)) {
    return false;
  }
  // (3) The parent node has at least one parent, and at least one left and right
  // child (including the added child node).
  const bool parent_has_parent =
      (left_rootward_of_parent > 0) || (right_rootward_of_parent > 0);
  const bool parent_has_children =
      (left_leafward_of_parent > 0) && (right_leafward_of_parent > 0);
  if ((parent_has_parent && parent_has_children) == false) {
    return false;
  }
  // (4) The child node has at least one parent, and at least one rotated and sorted
  // child.
  // (* we know child node has a parent node, so only need to check children)
  const bool child_has_children =
      (left_leafward_of_child > 0) && (right_leafward_of_child > 0);
  if (child_has_children == false) {
    return false;
  }
  return true;
}

bool SubsplitDAG::IsValidTaxonMap() const {
  std::vector<bool> id_exists(TaxonCount());
  // Get all ids from map.
  for (const auto &name_id : dag_taxa_) {
    size_t taxon_id = name_id.second;
    if (taxon_id > id_exists.size()) {
      return false;
    }
    if (id_exists[taxon_id]) {
      return false;
    } else {
      id_exists[taxon_id] = true;
    }
  }
  // Taxon map should cover the all ids from 0 to taxon_count-1.
  for (size_t i = 0; i < id_exists.size(); i++) {
    if (!id_exists[i]) {
      return false;
    }
  }
  return true;
}

// ** Reindexer methods:

// NOTE: To be performed *after* DAG modification.
SizeVector SubsplitDAG::BuildNodeReindexer(const size_t prev_node_count) {
  SizeVector node_reindexer = Reindexer::IdentityReindexer(NodeCount());
  // Begin reindex values at taxon count to account for ...leaves?
  size_t running_traversal_idx = taxon_count_;
  size_t dag_root_node_id = prev_node_count - 1;
  // Build node_reindexer by using post-order traversal (topological sort) of entire DAG
  // to assign new ids, where the index is the "before" node_id (stored in the node
  // object), and the value is the "after" node_id.
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
      Reindexer::IdentityReindexer(EdgeCountWithLeafSubsplits());
  // Only edges from an existing parent node to a new child node need to be reindexed.
  // See SubsplitDAG::AddNodePair().
  for (size_t edge_idx = prev_edge_count; edge_idx < EdgeCountWithLeafSubsplits();
       edge_idx++) {
    // Find edge with given idx.
    auto element = storage_.GetLine(edge_idx);
    Assert(element.has_value(),
           "An edge with given edge_idx did not exist in "
           "SubsplitDAG::BuildEdgeReindexer.");
    //
    auto [node_pair, idx] = element.value();
    std::ignore = idx;
    const auto &[parent_id, child_id] = node_pair;
    const Bitset parent_subsplit = GetDAGNode(parent_id).GetBitset();
    const Bitset child_subsplit = GetDAGNode(child_id).GetBitset();
    const auto idx_range = GetChildEdgeRange(
        parent_subsplit, child_subsplit.SubsplitIsLeftChildOf(parent_subsplit));
    // New edge is added to the end of the range.
    const size_t new_idx = edge_reindexer.at(idx_range.second);
    Reindexer::ReassignAndShift(edge_reindexer, edge_idx, new_idx);
  }
  return edge_reindexer;
}

void SubsplitDAG::RemapNodeIds(const SizeVector &node_reindexer) {
  std::vector<DAGVertex> nodes = {storage_.GetVertices().begin(),
                                  storage_.GetVertices().end()};
  std::vector<DAGVertex> nodes_copy = Reindexer::Reindex(nodes, node_reindexer);
  storage_.SetVertices(nodes_copy);

  // Update each node's id and leafward/rootward ids.
  for (size_t node_id = 0; node_id < NodeCount(); node_id++) {
    GetDAGNode(node_id).RemapNodeIds(node_reindexer);
  }
  // Update `subsplit_to_id_`.
  for (const auto &[subsplit, node_id] : subsplit_to_id_) {
    subsplit_to_id_.at(subsplit) = node_reindexer.at(node_id);
  }
  // Update edges.
  for (auto i : storage_.GetLines()) {
    storage_.ReindexLine(i.GetId(), node_reindexer.at(i.GetParent()),
                         node_reindexer.at(i.GetChild()));
  }
}

void SubsplitDAG::RemapEdgeIdxs(const SizeVector &edge_reindexer) {
  // Update edges.
  std::vector<DAGLineStorage> edges_copy(storage_.GetLines().size());
  for (auto i : storage_.GetLines()) {
    LineId new_idx = edge_reindexer.at(i.GetId());
    edges_copy[new_idx] = i;
    edges_copy[new_idx].SetId(new_idx);
  }
  storage_.SetLines(edges_copy);
  // Update `parent_to_child_range_`.
  for (const auto &[subsplit, idx_range] : parent_to_child_range_) {
    parent_to_child_range_.at(subsplit) = {edge_reindexer.at(idx_range.first),
                                           edge_reindexer.at(idx_range.second)};
  }
}
