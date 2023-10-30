// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "tp_choice_map.hpp"

// ** Access

EdgeId TPChoiceMap::GetEdgeChoice(const EdgeId edge_id,
                                  EdgeAdjacent edge_choice_type) const {
  switch (edge_choice_type) {
    case EdgeAdjacent::Parent:
      return edge_choice_vector_[edge_id.value_].parent;
    case EdgeAdjacent::Sister:
      return edge_choice_vector_[edge_id.value_].sister;
    case EdgeAdjacent::LeftChild:
      return edge_choice_vector_[edge_id.value_].left_child;
    case EdgeAdjacent::RightChild:
      return edge_choice_vector_[edge_id.value_].right_child;
    default:
      Failwith("Invalid edge choice type.");
  }
}

void TPChoiceMap::SetEdgeChoice(const EdgeId edge_id,
                                const EdgeAdjacent edge_choice_type,
                                const EdgeId new_edge_choice) {
  switch (edge_choice_type) {
    case EdgeAdjacent::Parent:
      edge_choice_vector_[edge_id.value_].parent = new_edge_choice;
      break;
    case EdgeAdjacent::Sister:
      edge_choice_vector_[edge_id.value_].sister = new_edge_choice;
      break;
    case EdgeAdjacent::LeftChild:
      edge_choice_vector_[edge_id.value_].left_child = new_edge_choice;
      break;
    case EdgeAdjacent::RightChild:
      edge_choice_vector_[edge_id.value_].right_child = new_edge_choice;
      break;
    default:
      Failwith("Invalid edge choice type.");
  }
}

void TPChoiceMap::ResetEdgeChoice(const EdgeId edge_id) {
  edge_choice_vector_[edge_id.value_] = {EdgeId(NoId), EdgeId(NoId), EdgeId(NoId),
                                         EdgeId(NoId)};
}

TPChoiceMap::EdgeChoiceNodeIds TPChoiceMap::GetEdgeChoiceNodeIds(
    const EdgeId edge_id) const {
  const auto &edge_choice = GetEdgeChoice(edge_id);
  return GetEdgeChoiceNodeIds(edge_choice);
}

TPChoiceMap::EdgeChoiceNodeIds TPChoiceMap::GetEdgeChoiceNodeIds(
    const EdgeChoice &edge_choice) const {
  auto GetNodeFromEdge = [this](const EdgeId edge_id, const Direction direction) {
    if (edge_id == NoId) {
      return NodeId(NoId);
    }
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    return (direction == Direction::Rootward) ? edge.GetParent() : edge.GetChild();
  };

  EdgeChoiceNodeIds node_choice;
  node_choice.parent = GetNodeFromEdge(edge_choice.parent, Direction::Rootward);
  node_choice.sister = GetNodeFromEdge(edge_choice.sister, Direction::Leafward);
  node_choice.left_child = GetNodeFromEdge(edge_choice.left_child, Direction::Leafward);
  node_choice.right_child =
      GetNodeFromEdge(edge_choice.right_child, Direction::Leafward);
  return node_choice;
}

TPChoiceMap::EdgeChoicePCSPs TPChoiceMap::GetEdgeChoicePCSPs(
    const EdgeId edge_id) const {
  EdgeChoicePCSPs adj_pcsps;
  const auto &edge_choice = GetEdgeChoice(edge_id);
  auto SetPCSPIfExists = [this](Bitset &adj_edge_bitset, const EdgeId adj_edge_id) {
    if (adj_edge_id != NoId) {
      adj_edge_bitset = GetDAG().GetDAGEdgeBitset(adj_edge_id);
    }
  };
  SetPCSPIfExists(adj_pcsps.parent, edge_choice.parent);
  SetPCSPIfExists(adj_pcsps.focal, edge_id);
  SetPCSPIfExists(adj_pcsps.sister, edge_choice.sister);
  SetPCSPIfExists(adj_pcsps.left_child, edge_choice.left_child);
  SetPCSPIfExists(adj_pcsps.right_child, edge_choice.right_child);

  return adj_pcsps;
}

// ** Maintenance

// Grow and reindex data to fit new DAG. Initialize new choice map to first edge.
void TPChoiceMap::GrowEdgeData(const size_t new_edge_count,
                               std::optional<const Reindexer> edge_reindexer,
                               std::optional<const size_t> explicit_alloc,
                               const bool on_init) {
  edge_choice_vector_.resize(new_edge_count);
  if (edge_reindexer.has_value()) {
    auto &reindexer = edge_reindexer.value();
    // Remap edge choices.
    for (EdgeId edge_id(0); edge_id < new_edge_count; edge_id++) {
      for (const auto edge_choice_type : EdgeAdjacentEnum::Iterator()) {
        if (GetEdgeChoice(edge_id, edge_choice_type) != NoId) {
          SetEdgeChoice(edge_id, edge_choice_type,
                        EdgeId(reindexer.GetNewIndexByOldIndex(
                            size_t(GetEdgeChoice(edge_id, edge_choice_type)))));
        }
      }
    }
    // Reindex ordering of choice map.
    Reindexer::ReindexInPlace<EdgeChoiceVector, EdgeChoice>(edge_choice_vector_,
                                                            reindexer, new_edge_count);
  }
}

// ** Selectors

void TPChoiceMap::SelectFirstEdge() {
  for (EdgeId edge_id = EdgeId(0); edge_id < GetDAG().EdgeCountWithLeafSubsplits();
       edge_id++) {
    SelectFirstEdge(edge_id);
  }
}

void TPChoiceMap::SelectFirstEdge(const EdgeId edge_id) {
  const auto edge = GetDAG().GetDAGEdge(edge_id);
  const auto focal_clade = edge.GetSubsplitClade();
  const auto parent_node = GetDAG().GetDAGNode(NodeId(edge.GetParent()));
  const auto child_node = GetDAG().GetDAGNode(NodeId(edge.GetChild()));
  ResetEdgeChoice(edge_id);
  auto &edge_choice = GetEdgeChoice(edge_id);

  // Query neighbor nodes.
  const SizeVector &left_parents =
      parent_node.GetNeighbors(Direction::Rootward, SubsplitClade::Left);
  const SizeVector &right_parents =
      parent_node.GetNeighbors(Direction::Rootward, SubsplitClade::Right);
  const SizeVector &sisters =
      parent_node.GetNeighbors(Direction::Leafward, Bitset::Opposite(focal_clade));
  const SizeVector &left_children =
      child_node.GetNeighbors(Direction::Leafward, SubsplitClade::Left);
  const SizeVector &right_children =
      child_node.GetNeighbors(Direction::Leafward, SubsplitClade::Right);

  // If neighbor lists are non-empty, get first edge from list.
  if (!left_parents.empty()) {
    edge_choice.parent = GetDAG().GetEdgeIdx(NodeId(left_parents[0]), parent_node.Id());
  }
  if (!right_parents.empty()) {
    edge_choice.parent =
        GetDAG().GetEdgeIdx(NodeId(right_parents[0]), parent_node.Id());
  }
  if (!sisters.empty()) {
    edge_choice.sister = GetDAG().GetEdgeIdx(parent_node.Id(), NodeId(sisters[0]));
  }
  if (!left_children.empty()) {
    edge_choice.left_child =
        GetDAG().GetEdgeIdx(child_node.Id(), NodeId(left_children[0]));
  }
  if (!right_children.empty()) {
    edge_choice.right_child =
        GetDAG().GetEdgeIdx(child_node.Id(), NodeId(right_children[0]));
  }
}

bool TPChoiceMap::SelectionIsValid(const bool is_quiet) const {
  std::stringstream dev_null;
  std::ostream &os = (is_quiet ? dev_null : std::cerr);
  bool is_valid = true;

  EdgeId edge_max_id = GetDAG().EdgeIdxRange().second;
  for (EdgeId edge_id = EdgeId(0); edge_id < edge_choice_vector_.size(); edge_id++) {
    const auto &edge_choice = edge_choice_vector_[edge_id.value_];
    if ((edge_choice.parent == NoId) && (edge_choice.sister == NoId) &&
        (edge_choice.left_child == NoId) && (edge_choice.right_child == NoId)) {
      os << "Invalid Selection: Edge Choice is empty." << std::endl;
      is_valid = false;
    }
    // If edge id is outside valid range.
    if ((edge_choice.parent > edge_max_id) || (edge_choice.sister > edge_max_id)) {
      // If they are not NoId, then it is an invalid edge_id.
      if ((edge_choice.parent != NoId) || (edge_choice.sister != NoId)) {
        os << "Invalid Selection: Parent or Sister has invalid edge_id." << std::endl;
        is_valid = false;
      }
      // NoId is valid only if edge goes to a root.
      if (!GetDAG().IsEdgeRoot(edge_id)) {
        os << "Invalid Selection: Parent or Sister has NoId when edge is not a root."
           << std::endl;
        is_valid = false;
      }
    }
    for (const auto &child_edge_id :
         EdgeIdVector({edge_choice.left_child, edge_choice.right_child})) {
      // If edge id is outside valid range.
      if (child_edge_id > edge_max_id) {
        // If they are not NoId, then it is an invalid edge_id.
        if (child_edge_id != NoId) {
          os << "Invalid Selection: Child has invalid edge id." << std::endl;
          is_valid = false;
        }
        // NoId is valid only if edge goes to a leaf.
        if (!GetDAG().IsEdgeLeaf(edge_id)) {
          os << "Invalid Selection: Child has NoId when edge is not a leaf."
             << std::endl;
          is_valid = false;
        }
      }
    }
    if (!is_valid) {
      os << "Failed at Edge" << edge_id << ": IsLeaf? " << GetDAG().IsEdgeLeaf(edge_id)
         << ", IsRoot? " << GetDAG().IsEdgeRoot(edge_id) << std::endl;
      os << "EdgeChoice: " << EdgeChoiceToString(edge_id) << std::endl;
      os << std::endl;
      break;
    }
  }
  return is_valid;
}

// ** TreeMask

// - Makes two passes:
//   - The first pass goes up along the chosen edges of the DAG to the root, adding
//   each edge it encounters.
//   - The second pass goes leafward, descending along the chosen edges to the leaf
//   edges from the sister of each edge in the rootward pass and the child edges from
//   the focal edge.
TPChoiceMap::ExpandedTreeMask TPChoiceMap::ExtractExpandedTreeMask(
    const TreeMask &tree_mask) const {
  ExpandedTreeMask tree_mask_ext;
  for (const auto edge_id : tree_mask) {
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    const auto focal_clade = edge.GetSubsplitClade();
    const NodeId parent_id = edge.GetParent();
    const NodeId child_id = edge.GetChild();
    // Add nodes to map if they don't already exist.
    for (const auto &node_id : {parent_id, child_id}) {
      if (tree_mask_ext.find(node_id) == tree_mask_ext.end()) {
        NodeAdjacentArray<NodeId> adj_nodes;
        adj_nodes.fill(NodeId(NoId));
        tree_mask_ext.insert({node_id, adj_nodes});
      }
    }
    // Add adjacent nodes to map.
    const auto which_node = (focal_clade == SubsplitClade::Left)
                                ? NodeAdjacent::LeftChild
                                : NodeAdjacent::RightChild;
    Assert(tree_mask_ext[parent_id][which_node] == NoId,
           "Invalid TreeMask: Cannot reassign adjacent child node.");
    tree_mask_ext[parent_id][which_node] = child_id;
    Assert(tree_mask_ext[child_id][NodeAdjacent::Parent] == NoId,
           "Invalid TreeMask: Cannot reassign adjacent parent node.");
    tree_mask_ext[child_id][NodeAdjacent::Parent] = parent_id;
  }
  return tree_mask_ext;
}

// Convert stack to a vector.
template <typename T>
static std::set<T> StackToVector(std::stack<T> &stack) {
  T *end = &stack.top() + 1;
  T *begin = end - stack.size();
  std::vector<T> stack_contents(begin, end);
  return stack_contents;
}

// Push Id to Stack if it is not NoId.
template <typename T, typename IdType>
static void StackPushIfValidId(std::stack<T> &stack, IdType id) {
  if (id != NoId) {
    stack.push(id);
  }
}

TPChoiceMap::TreeMask TPChoiceMap::ExtractTreeMask(const EdgeId initial_edge_id) const {
  TreeMask tree_mask;
  std::stack<EdgeId> rootward_stack, leafward_stack;
  const auto edge_max_id = GetDAG().EdgeIdxRange().second;

  // Rootward Pass: Capture parent and sister edges above focal edge.
  // For focal edge, add children to stack for leafward pass.
  auto focal_edge_id = initial_edge_id;
  const auto &focal_choices = edge_choice_vector_.at(focal_edge_id.value_);
  StackPushIfValidId(rootward_stack, focal_choices.left_child);
  StackPushIfValidId(rootward_stack, focal_choices.right_child);
  // Follow parentage upward until root.
  bool at_root = false;
  while (!at_root) {
    if (focal_edge_id >= edge_max_id) {
      std::cout << "FOCAL EDGE OUT-OF-RANGE: " << focal_edge_id << " " << edge_max_id
                << std::endl;
    }
    Assert(focal_edge_id < edge_max_id, "edge_id is outside valid edge range.");
    tree_mask.insert(focal_edge_id);
    const auto &focal_choices = edge_choice_vector_.at(focal_edge_id.value_);
    // End upward pass if we are at the root.
    if (GetDAG().IsEdgeRoot(focal_edge_id)) {
      at_root = true;
    } else {
      // If not at root, add sister for leafward pass.
      focal_edge_id = focal_choices.parent;
      rootward_stack.push(focal_choices.sister);
    }
  }

  // Leafward Pass: Capture all children from parentage.
  while (!rootward_stack.empty()) {
    const auto parent_edge_id = rootward_stack.top();
    rootward_stack.pop();
    leafward_stack.push(parent_edge_id);
  }
  while (!leafward_stack.empty()) {
    const auto edge_id = leafward_stack.top();
    leafward_stack.pop();
    tree_mask.insert(edge_id);
    const auto edge_choice = edge_choice_vector_.at(edge_id.value_);
    StackPushIfValidId(leafward_stack, edge_choice.left_child);
    StackPushIfValidId(leafward_stack, edge_choice.right_child);
  }

  return tree_mask;
}

TPChoiceMap::ExpandedTreeMask TPChoiceMap::ExtractExpandedTreeMask(
    const EdgeId initial_edge_id) const {
  return ExtractExpandedTreeMask(ExtractTreeMask(initial_edge_id));
}

bool TPChoiceMap::TreeMaskIsValid(const TreeMask &tree_mask,
                                  const bool is_quiet) const {
  std::stringstream dev_null;
  std::ostream &os = (is_quiet ? dev_null : std::cerr);

  SizeVector node_ids;
  bool root_check = false;
  BoolVector leaf_check(GetDAG().TaxonCount(), false);
  // Node map for checking node connectivity.
  using NodeMap = std::map<NodeId, NodeAdjacentArray<bool>>;
  NodeMap nodemap_check;
  // Check each edge in tree mask.
  for (const auto edge_id : tree_mask) {
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    const auto &parent_node = GetDAG().GetDAGNode(NodeId(edge.GetParent()));
    const auto &child_node = GetDAG().GetDAGNode(NodeId(edge.GetChild()));
    // Check if edge goes to root exactly once.
    if (GetDAG().IsNodeRoot(parent_node.Id())) {
      if (root_check) {
        os << "Invalid TreeMask: Multiple edges going to tree root." << std::endl;
        return false;
      }
      root_check = true;
    }
    // Check if edge goes to each leaf exactly once.
    if (GetDAG().IsNodeLeaf(child_node.Id())) {
      const auto taxon_id = child_node.Id();
      if (leaf_check.at(taxon_id.value_)) {
        os << "Invalid TreeMask: Multiple edges going to a tree leaf." << std::endl;
        return false;
      }
      leaf_check.at(taxon_id.value_) = true;
    }
    // Update node map. If a node already has parent or child, it is an invalid
    // tree.
    for (const NodeId node_id : {parent_node.Id(), child_node.Id()}) {
      if (nodemap_check.find(node_id) == nodemap_check.end()) {
        nodemap_check.insert({node_id, {{false, false, false}}});
      }
    }
    const auto which_child = (edge.GetSubsplitClade() == SubsplitClade::Left)
                                 ? NodeAdjacent::LeftChild
                                 : NodeAdjacent::RightChild;
    if (nodemap_check[parent_node.Id()][which_child] == true) {
      os << "Invalid TreeMask: Node has muliple parents." << std::endl;
      return false;
    }
    nodemap_check[parent_node.Id()][which_child] = true;
    if (nodemap_check[child_node.Id()][NodeAdjacent::Parent] == true) {
      os << "Invalid TreeMask: Node has multiple children." << std::endl;
      return false;
    }
    nodemap_check[child_node.Id()][NodeAdjacent::Parent] = true;
  }
  // Check if all nodes are fully connected.
  for (const auto &[node_id, connections] : nodemap_check) {
    // Check children.
    if (!(connections[NodeAdjacent::LeftChild] ||
          connections[NodeAdjacent::RightChild])) {
      // If one child is not connected, neither should be.
      if (connections[NodeAdjacent::LeftChild] ^
          connections[NodeAdjacent::RightChild]) {
        os << "Invalid TreeMask: Node has only one child." << std::endl;
        return false;
      }
      if (!GetDAG().IsNodeLeaf(node_id)) {
        os << "Invalid TreeMask: Non-leaf node has no children." << std::endl;
        return false;
      }
    }
    // Check parent.
    if (!connections[NodeAdjacent::Parent]) {
      if (!GetDAG().IsNodeRoot(node_id)) {
        os << "Invalid TreeMask: Non-root node has no parent." << std::endl;
        return false;
      }
    }
  }
  // Check if spans root.
  if (!root_check) {
    os << "Invalid TreeMask: Tree does not span root." << std::endl;
    return false;
  }
  // Check if spans all leaf nodes.
  for (size_t i = 0; i < leaf_check.size(); i++) {
    if (!leaf_check[i]) {
      os << "Invalid TreeMask: Tree does not span all leaves." << std::endl;
      return false;
    }
  }

  return true;
}

std::string TPChoiceMap::TreeMaskToString(const TreeMask &tree_mask) const {
  std::stringstream os;
  os << "[ " << std::endl;
  for (const auto edge_id : tree_mask) {
    const auto &edge = GetDAG().GetDAGEdge(edge_id);
    os << "\t" << edge_id << ":(" << edge.GetParent() << "->" << edge.GetChild()
       << "), " << std::endl;
  }
  os << "]";
  return os.str();
}

std::string TPChoiceMap::ExpandedTreeMaskToString(
    const ExpandedTreeMask &tree_mask) const {
  std::stringstream os;
  os << "[ " << std::endl;
  for (const auto [node_id, adj_node_ids] : tree_mask) {
    os << "\t" << node_id << ":(" << adj_node_ids[NodeAdjacent::Parent] << ", "
       << adj_node_ids[NodeAdjacent::LeftChild] << ", "
       << adj_node_ids[NodeAdjacent::RightChild] << "), " << std::endl;
  }
  os << "]";
  return os.str();
}

// ** Topology

// - Makes two passes:
//   - The first pass goes up along the chosen edges of the DAG to the root, adding
//   each edge it encounters.
//   - The second pass goes leafward, descending along the chosen edges to the leaf
//   edges from the sister of each edge in the rootward pass and the child edges from
//   the focal edge.
Node::Topology TPChoiceMap::ExtractTopology(const EdgeId focal_edge_id) const {
  return ExtractTopology(ExtractTreeMask(focal_edge_id));
}

Node::Topology TPChoiceMap::ExtractTopology(const TreeMask &tree_mask) const {
  ExpandedTreeMask tree_mask_ext = ExtractExpandedTreeMask(tree_mask);
  return ExtractTopology(tree_mask_ext);
}

Node::Topology TPChoiceMap::ExtractTopology(ExpandedTreeMask &tree_mask_ext) const {
  const auto dag_root_id = GetDAG().GetDAGRootNodeId();
  Assert(tree_mask_ext.find(dag_root_id) != tree_mask_ext.end(),
         "DAG Root Id does not exist in ExpandedTreeMask map.");
  Assert(tree_mask_ext[dag_root_id][NodeAdjacent::LeftChild] != NoId,
         "DAG Root Id has no children in ExpandedTreeMask map.");
  const auto dag_rootsplit_id = tree_mask_ext[dag_root_id][NodeAdjacent::LeftChild];

  std::unordered_map<NodeId, bool> visited_left;
  std::unordered_map<NodeId, bool> visited_right;
  std::unordered_map<NodeId, Node::NodePtr> nodes;

  // Build tree skeleton.  Set all nodes to unvisited.
  for (const auto &[node_id, adj_node_ids] : tree_mask_ext) {
    std::ignore = adj_node_ids;
    visited_left[node_id] = false;
    visited_right[node_id] = false;
  }

  size_t node_id_counter = GetDAG().TaxonCount();
  NodeId next_node_id = NodeId(NoId);
  NodeId current_node_id = dag_root_id;
  // Continue until left and right children of rootsplit node have been visited.
  nodes[dag_rootsplit_id] = nullptr;
  while (nodes[dag_rootsplit_id] == nullptr) {
    // If right branch (and left branch) already visited, join child nodes and return
    // up the tree.
    if (visited_right[current_node_id]) {
      const auto left_child_id =
          tree_mask_ext[current_node_id][NodeAdjacent::LeftChild];
      const auto right_child_id =
          tree_mask_ext[current_node_id][NodeAdjacent::RightChild];
      nodes[current_node_id] = Node::Join(nodes.at(left_child_id),
                                          nodes.at(right_child_id), node_id_counter);
      node_id_counter++;
      next_node_id = tree_mask_ext[current_node_id][NodeAdjacent::Parent];
    }
    // If left branch already visited, go down the right branch.
    else if (visited_left[current_node_id]) {
      visited_right[current_node_id] = true;
      next_node_id = tree_mask_ext[current_node_id][NodeAdjacent::RightChild];
    }
    // If node is a leaf, create leaf node and return up the the tree
    else if (GetDAG().IsNodeLeaf(current_node_id)) {
      nodes[current_node_id] =
          Node::Leaf(current_node_id.value_, GetDAG().TaxonCount());
      next_node_id = tree_mask_ext[current_node_id][NodeAdjacent::Parent];
    }
    // If neither left or right child has been visited, go down the left branch.
    else {
      visited_left[current_node_id] = true;
      next_node_id = tree_mask_ext[current_node_id][NodeAdjacent::LeftChild];
    }
    Assert(next_node_id != current_node_id, "Node cannot be adjacent to itself.");
    current_node_id = NodeId(next_node_id);
  }

  Node::NodePtr topology = nodes[dag_rootsplit_id];
  Assert((nodes.size() + 1) == tree_mask_ext.size(),
         "Invalid TreeMask-to-Tree: Topology did not span every node in "
         "the TreeMask.");

  return topology;
}

// ** I/O

std::string TPChoiceMap::EdgeChoiceToString(const EdgeId edge_id) const {
  const auto &edge_choice = GetEdgeChoice(edge_id);
  std::stringstream os;

  auto PrintEdge = [this, &os](const std::string &name, const EdgeId edge_id) {
    NodeId parent_id, child_id;
    if (edge_id != NoId) {
      const auto &edge = GetDAG().GetDAGEdge(edge_id);
      parent_id = edge.GetParent();
      child_id = edge.GetChild();
    }
    os << name << ": " << edge_id << " -> (" << parent_id << "," << child_id << "), ";
  };

  os << "{ ";
  PrintEdge("focal", edge_id);
  PrintEdge("parent", edge_choice.parent);
  PrintEdge("sister", edge_choice.sister);
  PrintEdge("left_child", edge_choice.left_child);
  PrintEdge("right_child", edge_choice.right_child);
  os << " }";
  return os.str();
}

std::string TPChoiceMap::EdgeChoiceToString(
    const TPChoiceMap::EdgeChoice &edge_choice) {
  std::stringstream os;
  os << "{ ";
  os << "parent: " << edge_choice.parent << ", ";
  os << "sister: " << edge_choice.sister << ", ";
  os << "left_child: " << edge_choice.left_child << ", ";
  os << "right_child: " << edge_choice.right_child;
  os << " }";
  return os.str();
}

std::string TPChoiceMap::ToString() const {
  std::stringstream os;
  os << "[ " << std::endl;
  for (EdgeId edge_id = EdgeId(0); edge_id < size(); edge_id++) {
    os << "\t" << edge_id << ": ";
    os << EdgeChoiceToString(GetEdgeChoice(edge_id));
    if (edge_id.value_ + 1 < size()) {
      os << ", ";
    }
    os << std::endl;
  }
  os << "]";
  return os.str();
}

std::map<Bitset, std::vector<Bitset>> TPChoiceMap::BuildPCSPMap() const {
  std::map<Bitset, std::vector<Bitset>> pcsp_map;
  for (EdgeId edge_id(0); edge_id < GetDAG().EdgeCountWithLeafSubsplits(); edge_id++) {
    const auto pcsps = GetEdgeChoicePCSPs(edge_id);
    pcsp_map[pcsps.focal] = {
        {pcsps.parent, pcsps.sister, pcsps.left_child, pcsps.right_child}};
  }
  return pcsp_map;
}

std::ostream &operator<<(std::ostream &os, const TPChoiceMap::EdgeChoice &edge_choice) {
  os << TPChoiceMap::EdgeChoiceToString(edge_choice);
  return os;
}

std::ostream &operator<<(std::ostream &os, const TPChoiceMap &choice_map) {
  os << choice_map.ToString();
  return os;
}
