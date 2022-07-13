// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "choice_map.hpp"

// ** Selectors

void ChoiceMap::SelectFirstEdge() {
  for (size_t edge_idx = 0; edge_idx < dag_.EdgeCountWithLeafSubsplits(); edge_idx++) {
    const auto edge = dag_.GetDAGEdge(EdgeId(edge_idx));
    const auto focal_clade = edge.GetSubsplitClade();
    const auto parent_node = dag_.GetDAGNode(NodeId(edge.GetParent()));
    const auto child_node = dag_.GetDAGNode(NodeId(edge.GetChild()));
    edge_choice_vector_[edge_idx] = EdgeChoice();
    auto &edge_choice = edge_choice_vector_[edge_idx];

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
      edge_choice.parent_edge_id =
          dag_.GetEdgeIdx(NodeId(left_parents[0]), parent_node.Id());
    }
    if (!right_parents.empty()) {
      edge_choice.parent_edge_id =
          dag_.GetEdgeIdx(NodeId(right_parents[0]), parent_node.Id());
    }
    if (!sisters.empty()) {
      edge_choice.sister_edge_id =
          dag_.GetEdgeIdx(parent_node.Id(), NodeId(sisters[0]));
    }
    if (!left_children.empty()) {
      edge_choice.left_child_edge_id =
          dag_.GetEdgeIdx(child_node.Id(), NodeId(left_children[0]));
    }
    if (!right_children.empty()) {
      edge_choice.right_child_edge_id =
          dag_.GetEdgeIdx(child_node.Id(), NodeId(right_children[0]));
    }
  }
}

bool ChoiceMap::SelectionIsValid(const bool is_quiet) const {
  std::stringstream dev_null;
  std::ostream &os = (is_quiet ? dev_null : std::cerr);

  EdgeId edge_max_id = dag_.EdgeIdxRange().second;
  for (EdgeId edge_idx = EdgeId(0); edge_idx < edge_choice_vector_.size(); edge_idx++) {
    const auto &edge_choice = edge_choice_vector_[edge_idx];
    // If edge id is outside valid range.
    if ((edge_choice.parent_edge_id > edge_max_id) ||
        (edge_choice.sister_edge_id > edge_max_id)) {
      // If they are not NoId, then it is an invalid edge_id.
      if ((edge_choice.parent_edge_id != NoId) ||
          (edge_choice.sister_edge_id != NoId)) {
        os << "Parent or Sister has invalid edge id." << std::endl;
        return false;
      }
      // NoId is valid only if edge goes to a root.
      if (!dag_.IsEdgeRoot(edge_idx)) {
        os << "Parent or Sister has NoId when edge is not a root." << std::endl;
        return false;
      }
    }
    for (const auto &child_edge_id :
         {edge_choice.left_child_edge_id, edge_choice.right_child_edge_id}) {
      // If edge id is outside valid range.
      if (child_edge_id > edge_max_id) {
        // If they are not NoId, then it is an invalid edge_id.
        if (child_edge_id != NoId) {
          os << "Child has invalid edge id." << std::endl;
          return false;
        }
        // NoId is valid only if edge goes to a leaf.
        if (!dag_.IsEdgeLeaf(edge_idx)) {
          os << "Child has NoId when edge is not a leaf." << std::endl;
          return false;
        }
      }
    }
  }
  return true;
}

// ** TreeMask

// - Makes two passes:
//   - The first pass goes up along the chosen edges of the DAG to the root, adding
//   each edge it encounters.
//   - The second pass goes leafward, descending along the chosen edges to the leaf
//   edges from the sister of each edge in the rootward pass and the child edges from
//   the central edge.
ChoiceMap::ExpandedTreeMask ChoiceMap::ExtractExpandedTreeMask(
    const TreeMask &tree_mask) const {
  ExpandedTreeMask tree_mask_ext;
  for (const auto edge_id : tree_mask) {
    const auto &edge = dag_.GetDAGEdge(edge_id);
    const auto focal_clade = edge.GetSubsplitClade();
    const NodeId parent_id = NodeId(edge.GetParent());
    const NodeId child_id = NodeId(edge.GetChild());
    // Add nodes to map if they don't already exist.
    for (const auto &node_id : {parent_id, child_id}) {
      if (tree_mask_ext.find(node_id) == tree_mask_ext.end()) {
        tree_mask_ext.insert({node_id, AdjacentNodeArray<NodeId>(NodeId(NoId))});
      }
    }
    // Add adjacent nodes to map.
    const auto which_node = (focal_clade == SubsplitClade::Left)
                                ? AdjacentNode::LeftChild
                                : AdjacentNode::RightChild;
    Assert(tree_mask_ext[parent_id][which_node] == NoId,
           "Invalid TreeMask: Cannot reassign adjacent node.");
    tree_mask_ext[parent_id][which_node] = child_id;
    Assert(tree_mask_ext[child_id][AdjacentNode::Parent] == NoId,
           "Invalid TreeMask: Cannot reassign adjacent node.");
    tree_mask_ext[child_id][AdjacentNode::Parent] = parent_id;
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

ChoiceMap::TreeMask ChoiceMap::ExtractTreeMask(const EdgeId central_edge_id) const {
  TreeMask tree_mask;
  std::stack<EdgeId> rootward_stack, leafward_stack;
  const auto root_node_id = dag_.GetDAGRootNodeId();
  const auto edge_max_id = dag_.EdgeIdxRange().second;

  // Rootward Pass: Capture parent and sister edges above focal edge.
  // For central edge, add children to stack for leafward pass.
  auto focal_edge_id = central_edge_id;
  const auto &focal_choices = edge_choice_vector_.at(focal_edge_id);
  StackPushIfValidId(rootward_stack, focal_choices.left_child_edge_id);
  StackPushIfValidId(rootward_stack, focal_choices.right_child_edge_id);
  // Follow parentage upward until root.
  bool at_root = false;
  while (!at_root) {
    Assert(focal_edge_id < edge_max_id,
           "Focal edge idx is outside valid edge idx range.");
    tree_mask.insert(focal_edge_id);
    const auto &focal_choices = edge_choice_vector_.at(focal_edge_id);
    // End upward pass if we are at the root.
    if (dag_.GetDAGEdge(focal_edge_id).GetParent() == root_node_id) {
      at_root = true;
    } else {
      // If not at root, add sister for leafward pass.
      focal_edge_id = focal_choices.parent_edge_id;
      rootward_stack.push(focal_choices.sister_edge_id);
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
    const auto edge_choice = edge_choice_vector_.at(edge_id);
    StackPushIfValidId(leafward_stack, edge_choice.left_child_edge_id);
    StackPushIfValidId(leafward_stack, edge_choice.right_child_edge_id);
  }

  return tree_mask;
}

ChoiceMap::ExpandedTreeMask ChoiceMap::ExtractExpandedTreeMask(
    const EdgeId central_edge_id) const {
  return ExtractExpandedTreeMask(ExtractTreeMask(central_edge_id));
}

bool ChoiceMap::TreeMaskIsValid(const TreeMask &tree_mask, const bool is_quiet) const {
  std::stringstream dev_null;
  std::ostream &os = (is_quiet ? dev_null : std::cerr);

  SizeVector node_ids;
  bool root_check = false;
  BoolVector leaf_check(dag_.TaxonCount(), false);
  // Node map for checking node connectivity.
  using NodeMap = std::map<NodeId, AdjacentNodeArray<bool>>;
  NodeMap nodemap_check;

  for (const auto edge_id : tree_mask) {
    const auto &edge = dag_.GetDAGEdge(edge_id);
    const auto &parent_node = dag_.GetDAGNode(NodeId(edge.GetParent()));
    const auto &child_node = dag_.GetDAGNode(NodeId(edge.GetChild()));
    // Check if edge goes to root exactly once.
    if (dag_.IsNodeRoot(parent_node.Id())) {
      if (root_check) {
        os << "Invalid TreeMask: Multiple edges going to tree root." << std::endl;
        return false;
      }
      root_check = true;
    }
    // Check if edge goes to each leaf exactly once.
    if (dag_.IsNodeLeaf(child_node.Id())) {
      const auto taxon_id = child_node.Id();
      if (leaf_check.at(taxon_id)) {
        os << "Invalid TreeMask: Multiple edges going to a tree leaf." << std::endl;
        return false;
      }
      leaf_check.at(taxon_id) = true;
    }
    // Update node map. If a node already has parent or child, it is an invalid
    // tree.
    for (const NodeId node_id : {parent_node.Id(), child_node.Id()}) {
      if (nodemap_check.find(node_id) == nodemap_check.end()) {
        nodemap_check.insert({node_id, {{false, false, false}}});
      }
    }
    const auto which_child = (edge.GetSubsplitClade() == SubsplitClade::Left)
                                 ? AdjacentNode::LeftChild
                                 : AdjacentNode::RightChild;
    if (nodemap_check[parent_node.Id()][which_child] == true) {
      os << "Invalid TreeMask: Node has muliple parents." << std::endl;
      return false;
    }
    nodemap_check[parent_node.Id()][which_child] = true;
    if (nodemap_check[child_node.Id()][AdjacentNode::Parent] == true) {
      os << "Invalid TreeMask: Node has multiple children." << std::endl;
      return false;
    }
    nodemap_check[child_node.Id()][AdjacentNode::Parent] = true;
  }
  // Check if all nodes are fully connected.
  for (const auto &[node_id, connections] : nodemap_check) {
    // Check children.
    if (!(connections[AdjacentNode::LeftChild] ||
          connections[AdjacentNode::RightChild])) {
      // If one child is not connected, neither should be.
      if (connections[AdjacentNode::LeftChild] ^
          connections[AdjacentNode::RightChild]) {
        os << "Invalid TreeMask: Node has only one child." << std::endl;
        return false;
      }
      if (!dag_.IsNodeLeaf(node_id)) {
        os << "Invalid TreeMask: Non-leaf node has no children." << std::endl;
        return false;
      }
    }
    // Check parent.
    if (!connections[AdjacentNode::Parent]) {
      if (!dag_.IsNodeRoot(node_id)) {
        os << "Invalid TreeMask: Non-root node has no parent." << std::endl;
        return false;
      }
    }
  }
  // Check if spans root and all leaf nodes.
  if (!root_check) {
    os << "Invalid TreeMask: Tree does not span root." << std::endl;
    return false;
  }
  for (size_t i = 0; i < leaf_check.size(); i++) {
    if (!leaf_check[i]) {
      os << "Invalid TreeMask: Tree does not span all leaves." << std::endl;
      return false;
    }
  }

  return true;
}

std::string ChoiceMap::TreeMaskToString(const TreeMask &tree_mask) const {
  std::stringstream os;
  os << "[ " << std::endl;
  for (const auto edge_id : tree_mask) {
    const auto &edge = dag_.GetDAGEdge(edge_id);
    os << "\t" << edge_id << ":(" << edge.GetParent() << "->" << edge.GetChild()
       << "), " << std::endl;
  }
  os << "]";
  return os.str();
}

std::string ChoiceMap::ExpandedTreeMaskToString(
    const ExpandedTreeMask &tree_mask) const {
  std::stringstream os;
  os << "[ " << std::endl;
  for (const auto [node_id, adj_node_ids] : tree_mask) {
    os << "\t" << node_id << ":(" << adj_node_ids[AdjacentNode::Parent] << ", "
       << adj_node_ids[AdjacentNode::LeftChild] << ", "
       << adj_node_ids[AdjacentNode::RightChild] << "), " << std::endl;
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
//   the central edge.
Node::NodePtr ChoiceMap::ExtractTopology(const EdgeId central_edge_id) const {
  return ExtractTopology(ExtractTreeMask(central_edge_id));
}

Node::NodePtr ChoiceMap::ExtractTopology(const TreeMask &tree_mask) const {
  ExpandedTreeMask tree_mask_ext = ExtractExpandedTreeMask(tree_mask);
  return ExtractTopology(tree_mask_ext);
}

Node::NodePtr ChoiceMap::ExtractTopology(ExpandedTreeMask &tree_mask_ext) const {
  const auto dag_root_id = dag_.GetDAGRootNodeId();
  Assert(tree_mask_ext.find(dag_root_id) != tree_mask_ext.end(),
         "DAG Root Id does not exist in ExpandedTreeMask map.");
  Assert(tree_mask_ext[dag_root_id][AdjacentNode::LeftChild] != NoId,
         "DAG Root Id has no children in ExpandedTreeMask map.");
  const auto dag_rootsplit_id = tree_mask_ext[dag_root_id][AdjacentNode::LeftChild];

  std::unordered_map<NodeId, bool> visited_left;
  std::unordered_map<NodeId, bool> visited_right;
  std::unordered_map<NodeId, Node::NodePtr> nodes;

  for (const auto &[node_id, adj_node_ids] : tree_mask_ext) {
    std::ignore = adj_node_ids;
    visited_left[node_id] = false;
    visited_right[node_id] = false;
  }

  size_t node_id_counter = dag_.TaxonCount();
  auto next_node_id = NoId;
  auto current_node_id = dag_root_id;
  // Continue until left and right children of rootsplit node have been visited.
  nodes[dag_rootsplit_id] = nullptr;
  while (nodes[dag_rootsplit_id] == nullptr) {
    // If right branch (and left branch) already visited, join child nodes and return
    // up the tree.
    if (visited_right[current_node_id]) {
      const auto left_child_id =
          tree_mask_ext[current_node_id][AdjacentNode::LeftChild];
      const auto right_child_id =
          tree_mask_ext[current_node_id][AdjacentNode::RightChild];
      nodes[current_node_id] = Node::Join(nodes.at(left_child_id),
                                          nodes.at(right_child_id), node_id_counter);
      node_id_counter++;
      next_node_id = tree_mask_ext[current_node_id][AdjacentNode::Parent];
    }
    // If left branch already visited, go down the right branch.
    else if (visited_left[current_node_id]) {
      visited_right[current_node_id] = true;
      next_node_id = tree_mask_ext[current_node_id][AdjacentNode::RightChild];
    }
    // If node is a leaf, return up the the tree
    else if (dag_.IsNodeLeaf(current_node_id)) {
      nodes[current_node_id] = Node::Leaf(current_node_id, dag_.TaxonCount());
      next_node_id = tree_mask_ext[current_node_id][AdjacentNode::Parent];
    }
    // If neither left or right child has been visited, go down the left branch.
    else {
      visited_left[current_node_id] = true;
      next_node_id = tree_mask_ext[current_node_id][AdjacentNode::LeftChild];
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

std::ostream &operator<<(std::ostream &os, const ChoiceMap::EdgeChoice &edge_choice) {
  os << "{ ";
  os << "parent: " << edge_choice.parent_edge_id << ", ";
  os << "sister: " << edge_choice.sister_edge_id << ", ";
  os << "left_child: " << edge_choice.left_child_edge_id << ", ";
  os << "right_child: " << edge_choice.right_child_edge_id;
  os << " }";
  return os;
}

std::ostream &operator<<(std::ostream &os, const ChoiceMap &choice_map) {
  os << choice_map.edge_choice_vector_;
  return os;
}
