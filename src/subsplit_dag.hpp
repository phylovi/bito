// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The purpose of this class is to hold a DAG built from the parent-child relationships
// of the subsplits. We wish to have information associated with both the nodes and
// edges of the DAG. Our strategy for doing that is via non-negative integer indices:
// nodes have a unique size_t `Id`, and each edge has a unique `gpcsp_idx`. We can then
// store arbitrary associated information in other data structures associated with these
// indices.
//
// The DAG has a well-defined notion of rootward and leafward.
//
// The data structure for this DAG is as follows:
// - The nodes of the DAG have vectors of indices representing their edges to other
// nodes. These include edges to the children (leafward edges) and also the edges that
// are connecting to that node (rootward edges). Furthermore, these sets of edges are
// separated into two classes: those that split apart the left component of the subsplit
// ("rotated" edges) and those that split apart the right component of the subsplit
// ("sorted" edges). The nodes of the DAG also include bitsets describing the subsplit
// that they represent.
// - The edges of the DAG are indexed separately, and there is a map (`dag_edges_`)
// which maps from pairs of node Ids to this edge index. These DAG edges are indexed
// such that all of the sorted edges descending from a given node have a contiguous set
// of indices, as do all of the rotated indices. The range of indices for such a set of
// edges is given by the `parent_to_range_` map.
// There is also a `subsplit_to_id_` map that maps from the subsplit bitset of the DAG
// node (and its rotated version) to the node Id.
// - To clarify terminology, the "DAG root node" refers to the universal ancestor and is
// representated by a bitset where the first half is all ones and the second half is all
// zeros (e.g. for taxon_count_ = 4, 1111|0000). Note that this implies that the DAG
// root node only has rotated children. Children of the DAG root node are called
// "rootsplits" and partition the whole taxon set.

#pragma once

#include "reindexer.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"
#include "subsplit_dag_action.hpp"
#include "subsplit_dag_nni.hpp"
#include "subsplit_dag_node.hpp"
#include "subsplit_dag_storage.hpp"

class SubsplitDAG {
 public:
  // NodeLambda is for iterating over nodes.
  using NodeLambda = std::function<void(const SubsplitDAGNode *)>;
  // EdgeDestinationLambda takes in a rotation status (true is rotated, false is not)
  // and a "destination" node. For iterating over DAG edges with a rotation status.
  using EdgeDestinationLambda = std::function<void(bool, const SubsplitDAGNode *)>;
  // EdgeAndNodeLambda takes a GPCSP index of an edge, its rotation status, and an index
  // of the node on the other side of the edge.
  using EdgeAndNodeLambda = std::function<void(const size_t, const bool, const size_t)>;
  // ParentEdgeChildLambda takes: the parent id in the DAG, the rotation status of the
  // edge, the child id, and the GCPSP index of the edge.
  using ParentRotationChildEdgeLambda =
      std::function<void(const size_t, const bool, const size_t, const size_t)>;

  // NodeAdditionResult is the return value of SubsplitDAG::AddNodePair.
  struct NodeAdditionResult {
    SizeVector new_node_ids, new_edge_idxs, node_reindexer, edge_reindexer;
  };

  SubsplitDAG();
  explicit SubsplitDAG(const RootedTreeCollection &tree_collection);

  // The total node count (including the DAG root node).
  size_t NodeCount() const;
  size_t NodeCountWithoutDAGRoot() const;
  // How many topologies can be expressed by the subsplit DAG? Expressed as a double
  // because this number can be big.
  double TopologyCount() const;
  size_t RootsplitCount() const;
  size_t GPCSPCount() const;
  size_t GPCSPCountWithFakeSubsplits() const;
  StringSizeMap SummaryStatistics() const;

  void Print() const;
  void PrintGPCSPIndexer() const;
  void PrintDAGEdges() const;
  void PrintParentToRange() const;
  void ToDot(const std::string file_path, bool show_index_labels = true) const;
  std::string ToDot(bool show_index_labels = true) const;

  // Create a GPCSPIndexer representing the DAG.
  // The gpcsp indexer is "expanded" meaning it contains fake PCSPs and rootsplit
  // bitsets are formatted as subsplits: 1110|0001.
  BitsetSizeMap BuildGPCSPIndexer() const;
  SubsplitDAGNode *GetDAGNode(size_t node_id) const;
  size_t GetDAGNodeId(const Bitset &subsplit) const;
  size_t DAGRootNodeId() const;
  // Return the node ids corresponding to the rootsplits.
  const SizeVector &RootsplitIds() const;
  const DAGTraits::Lines& GetEdges() const { return storage_.GetLines(); }
  const BitsetSizePairMap& GetParentToRange() const { return parent_to_range_; }
  // Access the GPCSP index from a parent-child pair of DAG nodes.
  size_t GetGPCSPIndex(const Bitset &parent_subsplit,
                       const Bitset &child_subsplit) const;
  // Get the GPCSP index from a parent-child pair of DAG nodes using the dag_edges_.
  size_t GPCSPIndexOfIds(size_t parent_id, size_t child_id) const;
  // Get the range of outgoing idxs from the given clade of a subsplit.
  SizePair GetEdgeRange(const Bitset &subsplit, const bool rotated) const;
  // Iterate over the "real" nodes, i.e. those that do not correspond to
  // fake subsplits or the DAG root node.
  void IterateOverRealNodes(const NodeLambda &f) const;
  // Iterate over the all leafward edges, rotated and sorted, of node using an
  // EdgeDestinationLambda.
  void IterateOverLeafwardEdges(const SubsplitDAGNode *node,
                                const EdgeDestinationLambda &f) const;
  // Iterate over only the rotated/sorted leafward edges of node using a NodeLambda.
  void IterateOverLeafwardEdges(const SubsplitDAGNode *node, bool rotated,
                                const NodeLambda &f) const;
  // Iterate over the leafward edges, supplying both the a GPCSP index of the edge and
  // the SubsplitDAGNode of the corresponding child.
  void IterateOverLeafwardEdgesAndChildren(const SubsplitDAGNode *node,
                                           const EdgeAndNodeLambda &f) const;
  // Iterate over the rootward edges of node using an EdgeDestinationLambda.
  // Excludes edges to the DAG root node.
  void IterateOverRootwardEdges(const SubsplitDAGNode *node,
                                const EdgeDestinationLambda &f) const;
  // Iterate over the rootward edges, supplying both the a GPCSP index of the edge and
  // the SubsplitDAGNode of the corresponding child.
  void IterateOverRootwardEdgesAndParents(const SubsplitDAGNode *node,
                                          const EdgeAndNodeLambda &f) const;
  // Iterate over the leafward edges, supplying the parent node id, child node id,
  // rotation of child, and the GPCSP index of the rootward edge connecting the two.
  void IterateOverParentAndChildAndLeafwardEdges(
      const SubsplitDAGNode *node, const ParentRotationChildEdgeLambda &f) const;

  // Each node in a topology is constructed with SubsplitDAGNode ID as Node ID.
  Node::NodePtrVec GenerateAllTopologies() const;

  // Apply an Action via a depth first traversal. Do not visit leaf nodes.
  // Applied to a given node, we:
  // * Apply BeforeNode
  // * For each of the clades of the node, we:
  //     * Apply BeforeNodeClade
  //     * For each edge descending from that clade, we:
  //         * Recur into the child node of the clade if it is not a leaf
  //         * Apply VisitEdge to the edge
  // * Apply AfterNode
  template <typename TraversalActionT>
  void DepthFirstWithAction(const SizeVector &starting_nodes,
                            const TraversalActionT &action) const {
    std::unordered_set<size_t> visited_nodes;
    for (const auto &node_id : starting_nodes) {
      DepthFirstWithActionForNode(action, node_id, visited_nodes);
    }
  };

  // The portion of the traversal that is below a given node.
  template <typename TraversalActionT>
  void DepthFirstWithActionForNode(const TraversalActionT &action, size_t node_id,
                                   std::unordered_set<size_t> &visited_nodes) const {
    action.BeforeNode(node_id);
    DepthFirstWithActionForNodeClade(action, node_id, false, visited_nodes);
    DepthFirstWithActionForNodeClade(action, node_id, true, visited_nodes);
    action.AfterNode(node_id);
  };

  // The portion of the traversal that is below a given clade of a given node.
  // Do not recur into leaf nodes.
  template <typename TraversalActionT>
  void DepthFirstWithActionForNodeClade(
      const TraversalActionT &action, size_t node_id, bool rotated,
      std::unordered_set<size_t> &visited_nodes) const {
    action.BeforeNodeClade(node_id, rotated);
    const auto node = GetDAGNode(node_id);
    for (const size_t child_id : node->GetLeafward(rotated)) {
      if (visited_nodes.count(child_id) == 0) {
        visited_nodes.insert(child_id);
        if (!GetDAGNode(child_id)->IsLeaf()) {
          DepthFirstWithActionForNode(action, child_id, visited_nodes);
        }
      }
      action.VisitEdge(node_id, child_id, rotated);
    }
  };

  // #288: consider renaming these re leafward and rootward.
  // Also, note that they could be called "TraversalTrace" to signify that they are
  // recording the trace of a traversal.
  [[nodiscard]] SizeVector LeafwardPassTraversal(bool include_dag_root_node) const;
  [[nodiscard]] SizeVector RootwardPassTraversal(bool include_dag_root_node) const;
  [[nodiscard]] SizeVector ReversePostorderTraversal() const;

  // Do a reverse postorder traversal on the edges of the DAG, including edges from the
  // DAG root node to the rootsplits, supplying the relevant indices to a lambda.
  void ReversePostorderIndexTraversal(ParentRotationChildEdgeLambda f) const;

  // Discrete uniform distribution over each subsplit.
  [[nodiscard]] EigenVectorXd BuildUniformQ() const;
  // Uniform prior over all topologies in the subsplit support.
  [[nodiscard]] EigenVectorXd BuildUniformOnTopologicalSupportPrior() const;
  // Uniform prior over all topologies, whether or not they are in the support.
  // Thus, this will only be a normalized probability distribution for each subsplit if
  // all topologies are in the support.
  [[nodiscard]] EigenVectorXd BuildUniformOnAllTopologiesPrior() const;

  RootedIndexerRepresentation IndexerRepresentationOf(const BitsetSizeMap &indexer,
                                                      const Node::NodePtr &topology,
                                                      size_t default_index) const;

  // Get a vector from each DAG node index to the probability of sampling that DAG node
  // with the supplied SBN parameters. These SBN parameters must be indexed in a
  // compatible way as the dag_edges_ of the subsplit DAG.
  EigenVectorXd UnconditionalNodeProbabilities(
      EigenConstVectorXdRef normalized_sbn_parameters) const;
  // Get a map from each non-fake subsplit to the probability of observing that
  // subsplit with the supplied SBN parameters. See
  // UnconditionalSubsplitProbabilityVector for notes.
  BitsetDoubleMap UnconditionalSubsplitProbabilities(
      EigenConstVectorXdRef normalized_sbn_parameters) const;
  // Get a vector from each GPCSP index to the Bayes-inverted probability of sampling
  // the parent given the child.
  EigenVectorXd InvertedGPCSPProbabilities(
      EigenConstVectorXdRef normalized_sbn_parameters,
      EigenConstVectorXdRef node_probabilities) const;

  // Does a node with the given subsplit exist?
  bool ContainsNode(const Bitset &subsplit) const;
  // Does a node with the given id exist?
  bool ContainsNode(const size_t node_id) const;
  // Does an edge that connects the two nodes exist?
  bool ContainsEdge(const size_t parent_id, const size_t child_id) const;
  // Get the rotated and sorted parents of the node with the given subsplit.
  std::pair<SizeVector, SizeVector> BuildParentIdVector(const Bitset &subsplit) const;
  // Get the rotated and sorted children of the node with the given subsplit.
  std::pair<SizeVector, SizeVector> BuildChildIdVector(const Bitset &subsplit) const;
  // Check if it is valid to add the node pair. More specifically, check that:
  // - The nodes are adjacent.
  // - The nodes do not add/remove any taxa.
  // - The parent node has at least 1 parent.
  // - Including the child node, each clade of the parent node has at least 1 child.
  // - Each clade of the child node has at least 1 child.
  bool IsValidNewNodePair(const Bitset &parent_subsplit,
                          const Bitset &child_subsplit) const;
  // Add an adjacent node pair to the DAG
  // and maintain the following indexing conventions for nodes and edges:
  // - Tips have ids 0 to taxon_count_.
  // - Parents have higher ids than their children.
  // - The DAG root node has highest id.
  // - Edges descending from the same node clade have contiguous idxs.
  // - There are no gaps in node ids or edge idxs.
  // Get the new_node_ids, new_edge_idxs, node_reindexer, and edge_reindexer.
  // Note: if both nodes already existed in the DAG, then new_node_ids and
  // new_edge_idxs will be empty.
  NodeAdditionResult AddNodePair(const Bitset &parent_subsplit,
                                 const Bitset &child_subsplit);

  // Get the reindexer for node ids.
  SizeVector BuildNodeReindexer(const size_t prev_node_count);
  // Get the reindexer for edges idxs.
  SizeVector BuildEdgeReindexer(const size_t prev_edge_count);
  // Remap all node ids according to the node_reindexer.
  void RemapNodeIds(const SizeVector &node_reindexer);
  // Remap all edge idxs according to the edge_reindexer.
  void RemapEdgeIdxs(const SizeVector &edge_reindexer);

 protected:
  size_t taxon_count_;
  size_t gpcsp_count_without_fake_subsplits_;
  // This indexer is an expanded version of parent_to_range_ in sbn_instance:
  // it includes single element range for fake subsplits.
  BitsetSizePairMap parent_to_range_;

  // A map from node id pairs to gpcsp idxs.
  DAGStorage storage_;
  //std::map<SizePair, size_t> dag_edges_;
  
  // We will call the index of DAG nodes "ids" to distinguish them from GPCSP indexes.
  // This corresponds to the analogous concept for topologies.
  //
  // A map from Bitset to the corresponding index in dag_nodes_.
  // The first entries are reserved for fake subsplits.
  // The last entries are reserved for rootsplits.
  // The DAG root node has the highest node id.
  BitsetSizeMap subsplit_to_id_;
  std::vector<std::unique_ptr<SubsplitDAGNode>> dag_nodes_;

  // Total number of topologies spanned by the DAG.
  double topology_count_;
  // Storage for the number of topologies below for each node.
  EigenVectorXd topology_count_below_;

  SubsplitDAG(size_t taxon_count, const Node::TopologyCounter &topology_counter);

  // Gives the child subsplits of a given parent subsplit, optionally including fake
  // subsplits.
  std::vector<Bitset> GetChildSubsplits(const SizeBitsetMap &index_to_child,
                                        const Bitset &subsplit,
                                        bool include_fake_subsplits = false);

  std::tuple<BitsetSizeMap, SizeBitsetMap, BitsetVector> ProcessTopologyCounter(
      const Node::TopologyCounter &topology_counter);
  void CreateAndInsertNode(const Bitset &subsplit);
  void CreateAndInsertEdge(const size_t parent_id, const size_t child_id, bool rotated);
  void ConnectGivenNodes(const size_t parent_id, const size_t child_id, bool rotated);
  void ConnectNodes(const SizeBitsetMap &index_to_child, size_t node_id, bool rotated);
  void BuildNodesDepthFirst(const SizeBitsetMap &index_to_child, const Bitset &subsplit,
                            std::unordered_set<Bitset> &visited_subsplits);
  void BuildNodes(const SizeBitsetMap &index_to_child, const BitsetVector &rootsplits);
  void BuildEdges(const SizeBitsetMap &index_to_child);
  void BuildDAGEdgesFromGPCSPIndexer(BitsetSizeMap &gpcsp_indexer);
  void CountTopologies();
  // Expand dag_edges_ and parent_to_range_ with fake subsplits at the end.
  void AddFakeSubsplitsToDAGEdgesAndParentToRange();
  void LeafwardDepthFirst(size_t node_id, SizeVector &visit_order,
                          std::unordered_set<size_t> &visited_nodes) const;
  void RootwardDepthFirst(size_t node_id, SizeVector &visit_order,
                          std::unordered_set<size_t> &visited_nodes) const;
  // Rotates subsplit only if it is out of sorted order (rotated).
  Bitset PerhapsSubsplitRotate(const Bitset &subsplit, bool rotated) const;

  // Note: the below methods are helper funcions for AddNodePair.
  // Connect the child to all of its children and update new_edge_idxs.
  void ConnectChildToAllChildren(const Bitset &child_subsplit,
                                 SizeVector &new_edge_idxs);
  // Connect the parent to all of its children except for the given child node and
  // update new_edge_idxs.
  void ConnectParentToAllChildrenExcept(const Bitset &parent_subsplit,
                                        const Bitset &child_subsplit,
                                        SizeVector &new_edge_idxs);
  // Connect the child to all of its parents except for the given parent node and update
  // new_edge_idxs.
  void ConnectChildToAllParentsExcept(const Bitset &parent_subsplit,
                                      const Bitset &child_subsplit,
                                      SizeVector &new_edge_idxs);
  // Connect the parent to all of its parents and update new_edge_idxs.
  void ConnectParentToAllParents(const Bitset &parent_subsplit,
                                 SizeVector &new_edge_idxs);
};
