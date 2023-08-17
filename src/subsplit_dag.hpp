// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The purpose of this class is to hold a DAG built from the parent-child relationships
// of the subsplits. We wish to have information associated with both the nodes and
// edges of the DAG. Our strategy for doing that is via non-negative integer indices:
// nodes have a unique size_t `Id`, and each edge has a unique `edge_idx`. We can then
// store arbitrary associated information in other data structures associated with these
// indices.
//
// The DAG has a well-defined notion of rootward and leafward.
//
// The data structure for this DAG is as follows:
// - The nodes of the DAG are vectors of indices representing their edges to other
// nodes. These include edges to the children (leafward edges) and also the edges that
// are connecting to that node (rootward edges). Furthermore, these sets of edges are
// separated into two subsets: those that descend from the left component of the
// subsplit
// ("rotated" edges) and those that descend from the right component of the subsplit
// ("sorted" edges). The nodes of the DAG also include bitsets describing the taxa they
// contain.
// - The edges of the DAG are indexed separately, and there is a map (`dag_edges_`)
// which maps from pairs of node Ids to its edge index. These DAG edges are indexed
// such that all of the sorted edges descending from a given node have a contiguous set
// of indices, as do all of the rotated indices. The range of indices for such a set of
// edges is given by the `parent_to_child_range_` map.
// There is also a `subsplit_to_id_` map that maps from the subsplit bitset of the DAG
// node (and its rotated version) to the node Id.
// - To clarify terminology, the "DAG root node" refers to the universal ancestor and is
// representated by a bitset where the first half is all ones and the second half is all
// zeros (e.g. for taxon_count_ = 4, 1111|0000). Note that this implies that the DAG
// root node only has rotated children. Children of the DAG root node are called
// "rootsplits" and partition the whole taxon set.

#pragma once

#include "resizer.hpp"
#include "reindexer.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"
#include "subsplit_dag_action.hpp"
#include "nni_operation.hpp"
#include "subsplit_dag_node.hpp"
#include "node.hpp"
#include "stopwatch.hpp"

class SubsplitDAG {
 public:
  // ** Constructor

  SubsplitDAG(const SubsplitDAG &) = default;
  // Build empty Subsplit DAG with no topologies and no taxa.
  SubsplitDAG();
  // Build a Subsplit DAG expressing all tree topologies from tree_collection.
  explicit SubsplitDAG(const RootedTreeCollection &tree_collection);

  // ** Comparator

  // This compare ensures that both DAGs have the same topology according to their
  // set of node and edge bitsets.  However, it does ensure that DAGs have the same id
  // and idxs for their respective nodes and edges, only that they contain the
  // same set of nodes and edges (as long as taxon positions in the
  // clades have the same mapping).
  int Compare(const SubsplitDAG &other, const bool quiet = true) const;
  static int Compare(const SubsplitDAG &lhs, const SubsplitDAG &rhs,
                     const bool quiet = true);
  friend bool operator==(const SubsplitDAG &lhs, const SubsplitDAG &rhs);
  friend bool operator!=(const SubsplitDAG &lhs, const SubsplitDAG &rhs);
  // Compares the subsplits between two DAGs. Returns the subsplits from lhs_dag not in
  // rhs_dag, and subsplits from the rhs_dag not in lhs_dag.
  static std::tuple<std::set<Bitset>, std::set<Bitset>, std::set<Bitset>>
  CompareSubsplits(const SubsplitDAG &lhs, const SubsplitDAG &rhs);
  // Compares the subsplits between two DAGs. Returns the subsplits from lhs_dag not in
  // rhs_dag, and subsplits from the rhs_dag not in lhs_dag.
  static std::tuple<std::set<Bitset>, std::set<Bitset>, std::set<Bitset>> ComparePCSPs(
      const SubsplitDAG &lhs, const SubsplitDAG &rhs);

  // ** Count

  // The total number of individual taxa in the DAG.
  size_t TaxonCount() const;
  // The total number of nodes in the DAG (including the root and leaves).
  size_t NodeCount() const;
  // The total number of nodes in the DAG (excluding the root, but including the
  // leaves).
  size_t NodeCountWithoutDAGRoot() const;
  // The current minimum and maximum node Id values.
  NodeIdPair NodeIdRange() const;
  // The total number of rootsplits in DAG. These count all direct descendants of the
  // root (also, the union of each rootsplits clades cover the set of all taxa in the
  // DAG).
  size_t RootsplitCount() const;
  // The total number of edges in the DAG (excluding edges which terminate at a root or
  // leaf node).
  size_t EdgeCount() const;
  // The total number of edges in the DAG (including edges which terminat at a root of
  // leaf node).
  size_t EdgeCountWithLeafSubsplits() const;
  // The current minimum and maximum edge Idx values.
  EdgeIdPair EdgeIdxRange() const;
  // The total number of tree topologies expressable by the DAG.
  double TopologyCount() const;
  // Checks how many valid neighbors nodes of given type exist in the DAG for the
  // specified subsplit.
  SizePair GetSubsplitNodeNeighborCounts(const Bitset &subsplit,
                                         const Direction direction) const;

  // ** I/O

  // Print all nodes:
  // - One line is given for (node_id | subsplit_bitset).
  // - One line is given for each combination of leafward/rootward, sorted/rotated, for
  // all adjacent node_ids.
  void Print() const;
  // Print all nodes as (node_id | node_bitsets) pairs, one-per-line.
  void PrintNodes() const;
  // Print all edges/PCSP, as (pcsp_bitset | edge_idx) pairs, one-per-line.
  void PrintEdgeIndexer() const;
  // Print all edges/PCSP, as (parent_node_id | child_node_id) pairs, one-per-line.
  void PrintDAGEdges() const;
  // Print all nodes, as (bitset | range_begin | range_end)
  void PrintParentToRange() const;
  // Print DAG Storage contents.
  void PrintStorage() const { storage_.Print(); }
  // Output DOT format graph of DAG to file.
  void ToDot(const std::string file_path, bool show_index_labels = true) const;
  // Output DOT format graph of DAG to a string.
  std::string ToDot(bool show_index_labels = true) const;

  // ** Build Indexers/Vectors

  // Create a EdgeIndexer representing the DAG.
  // The EdgeIndexer is a map (edge/PCSP bitset -> edge index).
  // The edge/PCSP indexer contains leafs and rootsplits.
  BitsetSizeMap BuildEdgeIndexer() const;
  // Builds inverse of EdgeIndexer map: (edge index -> edge/PCSP bitset).
  SizeBitsetMap BuildInverseEdgeIndexer() const;

  // Get the rotated and sorted parents of the node with the given subsplit.
  NodeIdVectorPair FindParentNodeIds(const Bitset &subsplit) const;
  NodeIdVectorPair FindParentNodeIdsViaMap(const Bitset &subsplit) const;
  NodeIdVectorPair FindParentNodeIdsViaScan(const Bitset &subsplit) const;
  // Get the rotated and sorted children of the node with the given subsplit.
  NodeIdVectorPair FindChildNodeIds(const Bitset &subsplit) const;
  NodeIdVectorPair FindChildNodeIdsViaMap(const Bitset &subsplit) const;
  NodeIdVectorPair FindChildNodeIdsViaScan(const Bitset &subsplit) const;

  // Output RootedIndexerRepresentation of DAG (from RootedSBNMaps).
  // RootedIndexerRepresentation is a vector of edge idxs in topological preorder.
  RootedIndexerRepresentation IndexerRepresentationOf(const BitsetSizeMap &indexer,
                                                      const Node::NodePtr &topology,
                                                      size_t default_index) const;

  // ** Access

  // Get Taxon's bitset clade positional id.
  TaxonId GetTaxonId(const std::string &name) const;
  // Get vector of all taxon ids.
  TaxonIdVector GetTaxonIds() const;
  // Get node based on node id.
  SubsplitDAGNode GetDAGNode(const NodeId node_id) const;
  MutableSubsplitDAGNode GetDAGNode(const NodeId node_id);
  // Get the subsplit bitset for the given node.
  Bitset GetDAGNodeBitset(const NodeId node_id) const;
  // Get the node id based on the subsplit bitset.
  NodeId GetDAGNodeId(const Bitset &subsplit) const;
  // Gets the node id of the DAG root.
  NodeId GetDAGRootNodeId() const;
  // Get the node ids corresponding to the rootsplits.
  ConstNeighborsView GetRootsplitNodeIds() const;
  // Get the edge ids corresponding to the rootsplits.
  EdgeIdVector GetRootsplitEdgeIds() const;
  // Get the node ids corresponding to the leaves for all taxa.
  NodeIdVector GetLeafNodeIds() const;
  // Get the node id corresponding to the given taxon.
  NodeId GetLeafNodeId(const TaxonId taxon_id) const;
  // Get leaf edge ids for given taxon.
  EdgeIdVector GetLeafEdgeIds(const TaxonId taxon_id) const;
  // Get edge based on edge id.
  ConstLineView GetDAGEdge(const EdgeId edge_id) const;
  // Get the PCSP bitset for the given edge.
  Bitset GetDAGEdgeBitset(const EdgeId edge_id) const;
  // Get the edge id by for given PCSP.
  EdgeId GetEdgeIdx(const Bitset &parent_subsplit, const Bitset &child_subsplit) const;
  EdgeId GetEdgeIdx(const NodeId parent_id, const NodeId child_id) const;
  EdgeId GetEdgeIdx(const Bitset &edge_pcsp) const;
  EdgeId GetEdgeIdx(const NNIOperation &nni) const;
  // Get focal clade of edge.
  SubsplitClade GetFocalClade(const EdgeId edge_id) const;
  // Get sister clade of edge.
  SubsplitClade GetSisterClade(const EdgeId edge_id) const;

  // Get NNI from edge index.
  NNIOperation GetNNI(const EdgeId edge_id) const;
  // Find a neighboring NNI. Returns the first NNI discovered.
  NNIOperation FindNNINeighborInDAG(const NNIOperation &nni) const;
  // Finds all the neighboring NNIs in the DAG. Array accessed by the child clade
  // swapped with the sister to produce the neighbor NNI.
  SubsplitCladeEnum::Array<std::optional<NNIOperation>> FindAllNNINeighborsInDAG(
      const NNIOperation &nni) const;

  // Get the range of outgoing idxs from the given clade of a subsplit.
  EdgeIdPair GetChildEdgeRange(const Bitset &subsplit,
                               const bool is_edge_on_left) const;
  // Get set of all taxon names.
  StringVector BuildSetOfTaxonNames() const;
  // Get set of all node Subsplit bitsets.
  std::set<Bitset> BuildSetOfNodeBitsets() const;
  // Get set of all edge PCSP bitsets.
  std::set<Bitset> BuildSetOfEdgeBitsets() const;
  // Get reference to taxon map.
  const StringTaxonIdMap &GetTaxonMap() const;
  // Get reference to the tag taxon map.
  const TagStringMapOption GetTagTaxonMap() const;
  // Get reference to subsplit -> node_id map.
  const BitsetNodeIdMap &GetSubsplitToIdMap() const;
  // Get reference to parent_node -> child_edge_range map.
  const NodeIdEdgeIdPairMap &GetParentNodeToChildEdgeRangeMap() const;

  // ** DAG Lambda Iterators
  // These methods iterate over the nodes and take lambda functions with arguments
  // relative to current node.

  // NodeLambda is for iterating over nodes.
  using NodeLambda = std::function<void(SubsplitDAGNode)>;
  // EdgeDestinationLambda takes in a rotation status (true is rotated, false is not)
  // and a "destination" node. For iterating over DAG edges with a rotation status.
  using EdgeDestinationLambda = std::function<void(bool, SubsplitDAGNode)>;
  // EdgeAndNodeLambda takes a PCSP index of an edge, its rotation status, and an index
  // of the node on the other side of the edge.
  using EdgeAndNodeLambda = std::function<void(const EdgeId, const bool, const NodeId)>;
  // ParentEdgeChildLambda takes: the parent id in the DAG, the rotation status of the
  // edge, the child id, and the GCPSP index of the edge.
  using ParentRotationChildEdgeLambda =
      std::function<void(const NodeId, const bool, const NodeId, const EdgeId)>;

  // Iterate over the "real" nodes, i.e. those that do not correspond to
  // leaf subsplits or the DAG root node.
  void IterateOverRealNodes(const NodeLambda &f) const;
  // Iterate over the all leafward edges, rotated and sorted, of node using an
  // EdgeDestinationLambda.
  void IterateOverLeafwardEdges(SubsplitDAGNode node,
                                const EdgeDestinationLambda &f) const;
  // Iterate over only the rotated/sorted leafward edges of node using a NodeLambda.
  void IterateOverLeafwardEdges(SubsplitDAGNode node, bool rotated,
                                const NodeLambda &f) const;
  // Iterate over the leafward edges, supplying both the index of the edge and
  // the SubsplitDAGNode of the corresponding child.
  void IterateOverLeafwardEdgesAndChildren(SubsplitDAGNode node,
                                           const EdgeAndNodeLambda &f) const;
  // Iterate over the rootward edges of node using an EdgeDestinationLambda.
  // Excludes edges to the DAG root node.
  void IterateOverRootwardEdges(SubsplitDAGNode node,
                                const EdgeDestinationLambda &f) const;
  // Iterate over the rootward edges, supplying both the a PCSP index of the edge and
  // the SubsplitDAGNode of the corresponding child.
  void IterateOverRootwardEdgesAndParents(SubsplitDAGNode node,
                                          const EdgeAndNodeLambda &f) const;
  // Iterate over the leafward edges, supplying the parent node id, child node id,
  // rotation of child, and the PCSP index of the rootward edge connecting the two.
  void IterateOverParentAndChildAndLeafwardEdges(
      SubsplitDAGNode node, const ParentRotationChildEdgeLambda &f) const;

  // ** DAG Traversals
  // These perform a traversal of the DAG and takes an argument TraversalActionT.
  // TraversalActionT passes four functions of the following form and order:
  // (1)  void BeforeNode(size_t node_id)
  // (2)  void AfterNode(size_t node_id)
  // (3)  void BeforeNodeClade(size_t node_id, bool rotated)
  // (4)  void VisitEdge(size_t node_id, size_t child_id, bool rotated)
  // See below for order of usage.

  // Apply a TraversalAction via a depth first traversal. Do not visit leaf nodes.
  // Applied to a given node, we:
  // - Apply BeforeNode()
  // - For each of the clades of the node, we:
  //     - Apply BeforeNodeClade()
  //     - For each edge descending from that clade, we:
  //         - Recur into the child node of the clade if it is not a leaf
  //         - Apply VisitEdge() to the edge
  // - Apply AfterNode()
  template <typename TraversalActionT>
  void DepthFirstWithAction(const NodeIdVector &starting_nodes,
                            const TraversalActionT &action) const {
    std::unordered_set<NodeId> visited_nodes;
    for (const auto &node_id : starting_nodes) {
      DepthFirstWithActionForNode(action, NodeId(node_id), visited_nodes);
    }
  };
  // The portion of the traversal that is below a given node.
  template <typename TraversalActionT>
  void DepthFirstWithActionForNode(const TraversalActionT &action, NodeId node_id,
                                   std::unordered_set<NodeId> &visited_nodes) const {
    action.BeforeNode(node_id);
    DepthFirstWithActionForNodeClade(action, node_id, false, visited_nodes);
    DepthFirstWithActionForNodeClade(action, node_id, true, visited_nodes);
    action.AfterNode(node_id);
  };
  // The portion of the traversal that is below a given clade of a given node.
  // Does not recurse into leaf nodes.
  template <typename TraversalActionT>
  void DepthFirstWithActionForNodeClade(
      const TraversalActionT &action, NodeId node_id, bool is_edge_on_left,
      std::unordered_set<NodeId> &visited_nodes) const {
    action.BeforeNodeClade(node_id, is_edge_on_left);
    const auto node = GetDAGNode(node_id);
    for (const auto child_id : node.GetLeafward(is_edge_on_left)) {
      if (visited_nodes.count(NodeId(child_id)) == 0) {
        visited_nodes.insert(NodeId(child_id));
        if (!GetDAGNode(NodeId(child_id)).IsLeaf()) {
          DepthFirstWithActionForNode(action, NodeId(child_id), visited_nodes);
        }
      }
      action.VisitEdge(node_id, NodeId(child_id), is_edge_on_left);
    }
  };

  // ** DAG Node Traversals
  // These functions produce a vector of node IDs representing an ordered traversal
  // of the DAG.

  // Creates a vector of node IDs representing a leafward DFS postorder traversal of
  // the DAG.
  [[nodiscard]] NodeIdVector LeafwardNodeTraversalTrace(
      bool include_dag_root_node) const;
  // Creates a vector of node IDs representing a rootward DFS postorder traversal of
  // the DAG.
  [[nodiscard]] NodeIdVector RootwardNodeTraversalTrace(
      bool include_dag_root_node) const;
  // Creates a vector of node IDs representing a reverse DFS postorder, leafward
  // traversal of the DAG. NOTE: A reverse postorder traversal represents a leafward
  // topological sort.
  [[nodiscard]] NodeIdVector TopologicalNodeTraversalTrace() const;

  // ** DAG Edge Traversals

  // Creates a vector of edge IDs representing a leafward DFS postorder traversal of
  // the DAG.
  [[nodiscard]] EdgeIdVector LeafwardEdgeTraversalTrace(
      bool include_dag_root_node) const;
  // Creates a vector of edge IDs representing a rootward DFS postorder traversal of
  // the DAG.
  [[nodiscard]] EdgeIdVector RootwardEdgeTraversalTrace(
      bool include_dag_root_node) const;
  // Creates a vector of edge IDs representing a reverse DFS postorder, leafward
  // traversal of the DAG. NOTE: A reverse postorder traversal represents a leafward
  // topological sort.
  [[nodiscard]] EdgeIdVector TopologicalEdgeTraversalTrace(
      bool include_dag_root_node) const;
  // Do a topological traversal on the edges of the DAG, including edges from the
  // DAG root node to the rootsplits, supplying the relevant indices to a lambda.
  void TopologicalEdgeTraversal(ParentRotationChildEdgeLambda f) const;

  // ** Priors

  // Discrete uniform distribution over each subsplit.
  [[nodiscard]] EigenVectorXd BuildUniformQ() const;
  // Uniform prior over all topologies in the subsplit support.
  [[nodiscard]] EigenVectorXd BuildUniformOnTopologicalSupportPrior() const;
  // Uniform prior over all topologies, whether or not they are in the support.
  // Thus, this will only be a normalized probability distribution for each subsplit if
  // all topologies are in the support.
  [[nodiscard]] EigenVectorXd BuildUniformOnAllTopologiesPrior() const;

  // Get a vector from each DAG node index to the probability of sampling that DAG node
  // with the supplied SBN parameters. These SBN parameters must be indexed in a
  // compatible way as the dag_edges_ of the subsplit DAG.
  EigenVectorXd UnconditionalNodeProbabilities(
      EigenConstVectorXdRef normalized_sbn_parameters) const;
  // Get a map from each non-leaf subsplit to the probability of observing that
  // subsplit with the supplied SBN parameters. See
  // UnconditionalSubsplitProbabilityVector for notes.
  BitsetDoubleMap UnconditionalSubsplitProbabilities(
      EigenConstVectorXdRef normalized_sbn_parameters) const;
  // Get a vector from each PCSP index to the Bayes-inverted probability of sampling
  // the parent given the child.
  EigenVectorXd InvertedGPCSPProbabilities(
      EigenConstVectorXdRef normalized_sbn_parameters,
      EigenConstVectorXdRef node_probabilities) const;

  // ** Query DAG

  // Does a taxon with the given name exist in DAG?
  bool ContainsTaxon(const std::string &name) const;
  // Does a node with the given subsplit exist in DAG?
  bool ContainsNode(const Bitset &subsplit) const;
  bool ContainsNode(const NodeId node_id) const;
  // Does an edge that connects the two nodes exist in DAG?
  bool ContainsEdge(const Bitset &parent_subsplit, const Bitset &child_subsplit) const;
  bool ContainsEdge(const NodeId parent_id, const NodeId child_id) const;
  bool ContainsEdge(const Bitset &edge_subsplit) const;
  bool ContainsEdge(const EdgeId edge_id) const;
  // Does the node pair of an NNI exist in DAG?
  bool ContainsNNI(const NNIOperation &nni) const;
  // Does the given tree/topology exist in DAG?
  bool ContainsTree(const RootedTree &tree, const bool is_quiet = true) const;
  bool ContainsTopology(const Node::Topology &topology,
                        const bool is_quiet = true) const;

  // Is node the root?
  bool IsNodeRoot(const NodeId node_id) const;
  // Is node a leaf?
  bool IsNodeLeaf(const NodeId node_id) const;
  // Does edge connect to the root node?
  bool IsEdgeRoot(const EdgeId edge_id) const;
  // Does edge connect to a leaf node?
  bool IsEdgeLeaf(const EdgeId edge_id) const;

  // ** Tree/Topologies

  std::unordered_map<NodeId, const Node *> BuildDAGNodeIdToTreeNodeMapFromTopology(
      const Node::Topology &topology) const;
  // Build map from NodeId in DAG to NodeId in topology.
  std::unordered_map<NodeId, size_t> BuildNodeIdMapFromTopology(
      const Node::Topology &topology) const;
  // Build map from EdgeId in DAG to edge's child NodeId in topology.
  std::unordered_map<EdgeId, SizePair> BuildEdgeIdMapFromTopology(
      const Node::Topology &topology) const;

  using ParentToChildNodeIdMap =
      std::unordered_map<NodeId, SubsplitCladeEnum::Array<NodeId>>;
  using ChildToParentNodeIdMap = std::unordered_map<NodeId, NodeId>;
  // Generate topology for a Node Id vector.
  Node::Topology BuildTopologyFromNodeIdMap(ParentToChildNodeIdMap &tree_map,
                                            NodeId rootsplit_id) const;
  // Build Tree from a given topology using given DAG branch lengths.
  RootedTree BuildTreeFromTopology(const Node::Topology &topology,
                                   const EigenVectorXd &dag_branch_lengths) const;

  // Output given tree to newick topology string using DAG's leaf labels (without branch
  // lengths).
  std::string TreeToNewickTree(const RootedTree &tree) const;
  // Output given tree to newick tree string using DAG's leaf labels (with branch
  // lengths).
  std::string TreeToNewickTopology(const RootedTree &tree) const;
  // Output given topology to newick topology string using DAG's leaf labels (with
  // branch lengths).
  std::string TopologyToNewickTopology(const Node::Topology &topology) const;

  // Generates all possible topologies contained in DAG.
  // Each node in a topology is constructed with SubsplitDAGNode ID as Node ID.
  Node::NodePtrVec GenerateAllTopologies() const;
  std::vector<RootedTree> GenerateAllTrees(
      const EigenVectorXd &dag_branch_lengths) const;
  std::string ToNewickOfAllTopologies() const;

  // Generate a set of tree topologies that span all nodes and edges in the DAG.
  Node::NodePtrVec GenerateCoveringTopologies() const;
  std::vector<RootedTree> GenerateCoveringTrees(
      const EigenVectorXd &dag_branch_length) const;
  std::string ToNewickOfCoveringTopologies() const;

  // ** Modify DAG
  // These methods are for directly modifying the DAG by adding or removing nodes and
  // edges. All these methods promise that all data owned by the SubsplitDAG will be
  // valid and up-to-date at their conclusion. This includes maintaining the following
  // DAG invariants:
  // - Tips have ids 0 to taxon_count_.
  // - Parents have higher ids than their children.
  // - The DAG root node has highest id.
  // - Edges descending from the same node clade have contiguous idxs.
  // - There are no gaps in node ids or edge idxs.
  // Returned is the ModificationResult, which will specify changes and reordering to
  // data that occurred during modification.

  // ModificationResult is the return value of Modification methods.
  // Contains the output needed to update all related data to reflect modifications to
  // DAG.
  struct ModificationResult {
    // Nodes that were added or removed by modification.
    NodeIdVector added_node_ids;
    // Edges that were added or removed by modification.
    EdgeIdVector added_edge_idxs;
    // New ordering of node ids relative to their ordering before DAG modification.
    Reindexer node_reindexer;
    // New ordering of edge idxs relative to their ordering before DAG modification.
    Reindexer edge_reindexer;
    // Previous node count before modification.
    size_t prv_node_count;
    // Previous edge count before modification.
    size_t prv_edge_count;
    // Current node count after modification.
    size_t cur_node_count;
    // Current edge count after modification.
    size_t cur_edge_count;
  };

  // Add an adjacent node pair to the DAG.
  virtual ModificationResult AddNodePair(const NNIOperation &nni);
  virtual ModificationResult AddNodePair(const Bitset &parent_subsplit,
                                         const Bitset &child_subsplit);
  // Add multiple nodes to the DAG.
  ModificationResult AddNodes(
      const std::vector<std::pair<Bitset, Bitset>> &node_subsplit_pairs);
  // Add multiple edges to the DAG.
  ModificationResult AddEdges(const std::vector<Bitset> &edge_subsplits);

  // Add all pontential edges to DAG. Building DAGs from a collection of trees can
  // result in a DAG that is not fully connected, in which one or more potentially
  // adjacent nodes in the DAG do not have an edge between them.
  ModificationResult FullyConnect();

  // Add all tree topologies in topology_counter to DAG.
  std::tuple<BitsetSizeMap, SizeBitsetMap, BitsetVector> ProcessTopologyCounter(
      const Node::TopologyCounter &topology_counter);

  // ** Reindexers
  // These methods are for building and applying reindexers.
  // A reindexer describes a remapping/transform of identifiers (i.e. node ids, edge
  // idxs) before and after a modification to the DAG. The index is before the
  // modification, the value is after the modification. Reindexers are used to perform a
  // transform on data vectors,

  // Uses a depth first DAG traversal to build a reindexer for node_ids.
  Reindexer BuildNodeReindexer(const size_t prev_node_count);
  // Get the reindexer for edges idxs.
  Reindexer BuildEdgeReindexer(const size_t prev_edge_count);
  // Remap all node ids according to the node_reindexer.
  void RemapNodeIds(const Reindexer &node_reindexer);
  // Remap all edge idxs according to the edge_reindexer.
  void RemapEdgeIdxs(const Reindexer &edge_reindexer);

  // ** Validation Tests
  // These methods are used to assert that a DAG is in a valid state or that given
  // operation will result in a valid DAG.

  // Checks if SubsplitDAG's corresponding data is consistent and up-to-date.
  // Specifically, checks that:
  // - subsplit_to_id_ map consistent with nodes in dag_nodes_.
  // - parent_to_child_range_ map consistent with each parent and child node's
  // neighbors.
  // - parent_to_child_range_ map consistent with each parent's edges.
  bool IsConsistent() const;
  // Checks if SubsplitDAG is in a valid state (assumes that DAG is consistent).
  // Specifically, checks that:
  // - Each node is valid. That is, either:
  //    - Node has zero parents and zero children.
  //    - Node has 1+ parents, 1+ sorted children, and 1+ rotated children.
  bool IsValid() const;
  // Check if it is valid to add given node pair.
  // Specifically, check that:
  // - The nodes are adjacent.
  // - The nodes do not add/remove any taxa.
  // - The parent node has at least one parent.
  // - Including the child node, each clade of the parent node has at least one child.
  // - Each clade of the child node has at least 1 child.
  bool IsValidAddNodePair(const Bitset &parent_subsplit,
                          const Bitset &child_subsplit) const;
  // Check if the taxon map is valid. Specifically, check that:
  // - No duplicate ids.
  // - Ids cover all clade bits from 0 to map_size.
  bool IsValidTaxonMap() const;

  // ** Miscellaneous

  // Creates a dictionary of summary statistics for DAG.
  // Includes node_count and edge_count.
  StringSizeMap SummaryStatistics() const;
  // Rotates Node Subsplit only if it is out of sorted order (rotated).
  Bitset SubsplitToSortedOrder(const Bitset &subsplit, bool rotated) const;

  // Builds a map from the nodes of one DAG onto another DAG.  Nodes are considered
  // matching when they represent the same Subsplit bitset. Map only contains the edges
  // common to both DAGs.
  static std::unordered_map<NodeId, NodeId> BuildNodeIdMapBetweenDAGs(
      const SubsplitDAG &dag_a, const SubsplitDAG &dag_b);
  // Builds a map from the edges of one DAG onto another DAG.  Edges are considered
  // matching when they represent the same PCSP bitset. Map only contains the edges
  // common to both DAGs.
  static std::unordered_map<EdgeId, EdgeId> BuildEdgeIdMapBetweenDAGs(
      const SubsplitDAG &dag_a, const SubsplitDAG &dag_b);

  // Build vector between from SubsplitDAGs dag_a to dag_b corresponding to their taxon
  // ids. Can be treated as a "map" with indices representing keys. Requires that both
  // SubsplitDAGs use the same taxon set.
  // - Map: [ index: dag_a's taxon_id => value: dag_b's taxon_id ]
  static SizeVector BuildTaxonTranslationMap(const SubsplitDAG &dag_a,
                                             const SubsplitDAG &dag_b);
  // Compare two bitsets (subsplits or PCSPs) from two different DAGs using a taxon
  // map.
  static int BitsetCompareViaTaxonTranslationMap(const Bitset &bitset_a,
                                                 const Bitset &bitset_b,
                                                 const SizeVector &taxon_map);
  // Translate bitset using the taxon translation map output positions.
  // NOTE:  forward_translate uses index in map as input bit locations, values in map as
  // output bit locations.
  //        If false, map is used vice versa.
  static Bitset BitsetTranslateViaTaxonTranslationMap(
      const Bitset &bitset, const SizeVector &taxon_map,
      const bool forward_translate = true);
  // Build a default taxon map for constructor with dummy taxon names:
  // E.g. {{0, "x0"}, {1, "x1"}, ...}
  static TagStringMap BuildDummyTagTaxonMap(const size_t taxon_count);

  // ** GraftDAG Helper

  // Remove Graft nodes and edges from Host DAG.
  void ResetHostDAG(SubsplitDAG &host_dag);

 protected:
  // ** GraftDAG Helper

  explicit SubsplitDAG(SubsplitDAG &host_dag, HostDispatchTag);

  // ** Count

  // Traverses the DAG and refreshes count of the number of topologies contained in the
  // DAG. Updates topology_count_ and topology_count_below_, which contains the count of
  void CountTopologies();
  // Update edge count without leaf subsplits.
  void CountEdgesWithoutLeafSubsplits();

  // ** DAG Traversal
  // NOTES: - visit_order is the output vector of node ids in specified traversal order.
  //        - visited_nodes is a set of all node ids reached by traversal.

  // Creates vector of node ids in leafward depth first post-order traversal of DAG.
  void LeafwardDepthFirst(NodeId node_id, NodeIdVector &visit_order,
                          std::unordered_set<NodeId> &visited_nodes) const;
  // Create vector of node ids in rootward depth first post-order traversal of DAG.
  void RootwardDepthFirst(NodeId node_id, NodeIdVector &visit_order,
                          std::unordered_set<NodeId> &visited_nodes) const;

  // ** Modify DAG Helpers
  // These methods help calling functions to modify the DAG, but do NOT ensure a valid
  // state at the end of the function. Do not necessarily handle remapping ids and idxs,
  // removing old references to deleted objects, etc.

  // Add taxon map to DAG.
  void BuildTaxonMap(const TagStringMap &tag_taxon_map);
  // Create Node and insert it into the DAG.  Returns ID of created node.
  NodeId CreateAndInsertNode(const Bitset &subsplit);
  // Create Edge between given nodes and insert it into the DAG. Returns ID of created
  // edge.
  EdgeId CreateAndInsertEdge(const NodeId parent_id, const NodeId child_id,
                             const bool is_edge_on_left);
  // Add edge between given parent and child nodes to the DAG.
  void ConnectGivenNodes(const NodeId parent_id, const NodeId child_id,
                         const bool is_edge_on_left, const EdgeId edge_id);
  // Add edges between node_id and all children in map.
  void ConnectNodes(const SizeBitsetMap &index_to_child, const NodeId node_id,
                    const bool is_edge_on_left);
  // Add nodes for all children in map.
  void BuildNodes(const SizeBitsetMap &index_to_child, const BitsetVector &rootsplits);
  // Add nodes in depth first ordering for children in map.
  void BuildNodesDepthFirst(const SizeBitsetMap &index_to_child, const Bitset &subsplit,
                            std::unordered_set<Bitset> &visited_subsplits);
  // Add edges from all parent nodes according to child nodes in map.
  void BuildEdges(const SizeBitsetMap &index_to_child);
  // Add edges to DAG according to node_id pairs in edge indexer.
  void BuildDAGEdgesFromEdgeIndexer(BitsetSizeMap &edge_indexer);
  // Connect the child to all of its children. Push all new edges to
  // added_edge_idxs.
  void ConnectChildToAllChildren(const Bitset &child_subsplit,
                                 EdgeIdVector &added_edge_idxs);
  // Connect the parent to all of its children except for the given child node. Insert
  // all new edges to added_edge_idxs vector.
  void ConnectParentToAllChildrenExcept(const Bitset &parent_subsplit,
                                        const Bitset &child_subsplit,
                                        EdgeIdVector &added_edge_idxs);
  // Connect the child to all of its parents except for the given parent node. Insert
  // all new edge to added_edge_idxs vector.
  void ConnectChildToAllParentsExcept(const Bitset &parent_subsplit,
                                      const Bitset &child_subsplit,
                                      EdgeIdVector &added_edge_idxs);
  // Connect the parent to all of its parents. Insert all new edges to
  // added_edge_idxs vector.
  void ConnectParentToAllParents(const Bitset &parent_subsplit,
                                 EdgeIdVector &added_edge_idxs);
  // Expand dag_edges_ and parent_to_child_range_ with leaf subsplits at the end.
  void AddLeafSubsplitsToDAGEdgesAndParentToRange();
  // Builds a vector of subsplits of all children , optionally including leaf nodes.
  BitsetVector GetChildSubsplits(const SizeBitsetMap &index_to_child,
                                 const Bitset &subsplit,
                                 bool include_leaf_subsplits = false);

  // Internal logic helper for adding node pair.  Assumes that validation check
  // has already been performed.
  ModificationResult AddNodePairInternals(const Bitset &parent_subsplit,
                                          const Bitset &child_subsplit);
  ModificationResult AddNodePairInternals(
      const std::vector<std::pair<Bitset, Bitset>> &node_subsplit_pairs);
  // Internal logic helper that inserts node pair without reindexing. Just appends and
  // adds nodes to modification result.
  bool AddNodePairInternalsWithoutReindexing(
      const std::vector<std::pair<Bitset, Bitset>> &node_subsplit_pairs,
      ModificationResult &mods);

  // Check if a node would have at least one valid neighboring parent exist in
  // the DAG for given subsplit.
  bool SubsplitNodeHasParent(const Bitset &node_subsplit) const;
  // Check if a node would have at least one valid neighboring right and one left child
  // exist in the DAG for given subsplit.
  bool SubsplitNodeHasLeftAndRightChild(const Bitset &node_subsplit) const;

  // ** Constructor Helpers

  // Build a Subsplit DAG on given number of taxa, expressing all tree topologies from
  // tree_collection, with trees on the given taxa names/labels.
  SubsplitDAG(size_t taxon_count, const Node::TopologyCounter &topology_counter,
              const TagStringMap &tag_taxon_map);

 private:
  void StoreEdgeIds();

 protected:
  // Underlying data containing nodes and edges.
  SubsplitDAGStorage storage_;
  // NOTE: When using unique identifiers, for DAG nodes (aka Subsplits) we use the term
  // "ids", and for edges (aka PCSPs) we use the term index or "idx", to more easily
  // distinguish the two. This corresponds to the analogous concept for topologies.

  // - Map of Taxon Names
  //    - [ Taxon Name => Taxon Id (position of the "on" bit in the clades) ]
  StringTaxonIdMap dag_taxa_;
  const TagStringMap *tag_taxon_map_ = nullptr;

  // - Map of all DAG Nodes:
  //    - [ Node Subsplit (Bitset) => Node Id ]
  // A node's id is equivalent to its index in dag_nodes_. The first entries are
  // reserved for leaf subsplits. The last entries are reserved for rootsplits. The DAG
  // root node has the highest node id.
  BitsetNodeIdMap subsplit_to_id_;

  BitsetNodeIdSetMap subsplit_union_;
  BitsetNodeIdSetMap subsplit_clade_;

  // - Map of all DAG Nodes:
  //    - [ Node Subsplit (Bitset) => (begin, end) Range of Child Ids ]
  // This indexer is an expanded version of parent_to_child_range_ in sbn_instance:
  // It includes single element range for leaf subsplits.
  BitsetEdgeIdPairMap parent_to_child_range_;

  // The number of taxa in the DAG. This is equivalent to the size of the clades in each
  // subsplit. Also equivalent to the number of leaf nodes in the DAG.
  size_t taxon_count_;
  // Number of internal edges in the DAG (excludes all edges that go to a root or
  // leaf).
  size_t edge_count_without_leaf_subsplits_;
  // Total number of tree topologies spanned by the DAG.
  double topology_count_;
  // Storage for the number of topologies below for each node. Each index maps to the
  // count for the corresponding node_id.
  EigenVectorXd topology_count_below_;
};
