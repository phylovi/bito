// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// A ChoiceMap is a per-edge map of the best adjacent edges applied to a SubsplitDAG
// for Top-Pruning. Used for selecting, updating, and extracting the top tree from the
// DAG. A ChoiceMap can generate a TreeMask, which is a list of edge ids which express a
// single, complete tree embedded in the DAG.

#pragma once

#include <stack>
#include "sugar.hpp"
#include "gp_dag.hpp"
#include "tree.hpp"
#include "tree_collection.hpp"

class ChoiceMap {
 public:
  // Per-edge choices of best adjacent edges.
  struct EdgeChoice {
    size_t parent_edge_id = NoId;
    size_t sister_edge_id = NoId;
    size_t left_child_edge_id = NoId;
    size_t right_child_edge_id = NoId;
    double tree_likelihood = -INFINITY;
  };
  using EdgeChoiceVector = std::vector<EdgeChoice>;

  ChoiceMap(GPDAG &dag)
      : dag_(dag), edge_choice_vector_(dag.EdgeCountWithLeafSubsplits()){};

  // ** Selectors

  // Naive choice selector. Chooses the first edge from each list of candidates.
  void SelectFirstEdge() {
    for (size_t edge_idx = 0; edge_idx < dag_.EdgeCountWithLeafSubsplits();
         edge_idx++) {
      const auto edge = dag_.GetDAGEdge(edge_idx);
      const auto focal_clade = edge.GetSubsplitClade();
      const auto parent_node = dag_.GetDAGNode(edge.GetParent());
      const auto child_node = dag_.GetDAGNode(edge.GetChild());
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
        edge_choice.parent_edge_id = dag_.GetEdgeIdx(left_parents[0], parent_node.Id());
      } else if (!right_parents.empty()) {
        edge_choice.parent_edge_id =
            dag_.GetEdgeIdx(right_parents[0], parent_node.Id());
      }
      if (!sisters.empty()) {
        edge_choice.sister_edge_id = dag_.GetEdgeIdx(parent_node.Id(), sisters[0]);
      }
      if (!left_children.empty()) {
        edge_choice.left_child_edge_id =
            dag_.GetEdgeIdx(child_node.Id(), left_children[0]);
      }
      if (!right_children.empty()) {
        edge_choice.right_child_edge_id =
            dag_.GetEdgeIdx(child_node.Id(), right_children[0]);
      }
    }
  }

  // Check if choice selection is valid.
  // Specifically, checks that:
  // - Every edge choice vector has a valid id for all options, unless...
  //   - Edge goes to root (NoId for sister and parent).
  //   - Edge goes to leaf (NoId for left and right child).
  // - Edges span every leaf and root node.
  bool SelectionIsValid(const bool is_quiet = true) const {
    size_t edge_max_id = dag_.EdgeIdxRange().second;
    for (size_t edge_idx = 0; edge_idx < edge_choice_vector_.size(); edge_idx++) {
      const auto &edge_choice = edge_choice_vector_[edge_idx];
      // If edge id is outside valid range.
      if ((edge_choice.parent_edge_id > edge_max_id) ||
          (edge_choice.sister_edge_id > edge_max_id)) {
        // If they are not NoId, then it is an invalid edge_id.
        if ((edge_choice.parent_edge_id != NoId) ||
            (edge_choice.sister_edge_id != NoId)) {
          if (!is_quiet) {
            std::cerr << "Parent or Sister has invalid edge id." << std::endl;
          }
          return false;
        }
        // NoId is valid only if edge goes to a root.
        if (!dag_.IsEdgeRoot(edge_idx)) {
          if (!is_quiet) {
            std::cerr << "Parent or Sister has NoId when edge is not a root."
                      << std::endl;
          }
          return false;
        }
      }
      for (const auto &child_edge_id :
           {edge_choice.left_child_edge_id, edge_choice.right_child_edge_id}) {
        // If edge id is outside valid range.
        if (child_edge_id > edge_max_id) {
          // If they are not NoId, then it is an invalid edge_id.
          if (child_edge_id != NoId) {
            if (!is_quiet) {
              std::cerr << "Child has invalid edge id." << std::endl;
            }
            return false;
          }
          // NoId is valid only if edge goes to a leaf.
          if (!dag_.IsEdgeLeaf(edge_idx)) {
            if (!is_quiet) {
              std::cerr << "Child has NoId when edge is not a leaf." << std::endl;
            }
            return false;
          }
        }
      }
    }
    return true;
  }

  // ** TreeMask
  // A TreeMask is an unordered vector of edge Ids, which represent a tree contained in
  // the DAG, from the selected subset of DAG edges.
  using TreeMask = std::vector<size_t>;

  // Convert stack to a vector.
  template <typename T>
  static std::vector<T> StackToVector(std::stack<T> &stack) {
    T *end = &stack.top() + 1;
    T *begin = end - stack.size();
    std::vector<T> stack_contents(begin, end);
    return stack_contents;
  }

  // Push Id to Stack if it is not NoId.
  template <typename T>
  static void StackPushIfValidId(std::stack<T> &stack, size_t id) {
    if (id != NoId) {
      stack.push(id);
    }
  }

  // Extract TreeMask from DAG based on edge choices to find best tree with given
  // central edge.
  // - Makes two passes:
  //   - The first pass goes up along the chosen edges of the DAG to the root, adding
  //   each edge it encounters.
  //   - The second pass goes leafward, descending along the chosen edges to the leaf
  //   edges from the sister of each edge in the rootward pass and the child edges from
  //   the central edge.
  TreeMask ExtractTreeMask(size_t central_edge_id) const {
    TreeMask tree_mask;
    std::stack<size_t> rootward_stack, leafward_stack;
    const size_t root_node_id = dag_.GetDAGRootNodeId();
    const auto edge_max_id = dag_.EdgeIdxRange().second;

    // Rootward Pass: Capture parent and sister edges above focal edge.
    // For central edge, add children to stack for leafward pass.
    size_t focal_edge_id = central_edge_id;
    const auto &focal_choices = edge_choice_vector_.at(focal_edge_id);
    StackPushIfValidId(rootward_stack, focal_choices.left_child_edge_id);
    StackPushIfValidId(rootward_stack, focal_choices.right_child_edge_id);
    // Follow parentage upward until root.
    bool at_root = false;
    while (!at_root) {
      Assert(focal_edge_id < edge_max_id,
             "Focal edge idx is outside valid edge idx range.");
      tree_mask.push_back(focal_edge_id);
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
      tree_mask.push_back(edge_id);
      const auto edge_choice = edge_choice_vector_.at(edge_id);
      StackPushIfValidId(leafward_stack, edge_choice.left_child_edge_id);
      StackPushIfValidId(leafward_stack, edge_choice.right_child_edge_id);
    }

    return tree_mask;
  }

  // Checks that TreeMask represents a valid, complete tree in the DAG.
  // Specifically, checks that:
  // - There is a single edge that goes to the root.
  // - There is a single edge that goes to each leaf.
  // - For each node in mask, there is a single parent, left and right child.
  //   - Unless node is root (no parent) or leaf (no children).
  bool TreeMaskIsValid(const TreeMask &tree_mask, const bool is_quiet = true) const {
    SizeVector node_ids;
    bool root_check = false;
    BoolVector leaf_check(dag_.TaxonCount(), false);
    // Node map for checking node connectivity: Array is ordered as [left_child,
    // right_child, parent].
    using NodeMap = std::map<size_t, std::array<bool, 3>>;
    NodeMap nodemap_check;

    for (size_t i = 0; i < tree_mask.size(); i++) {
      const auto edge_idx = tree_mask[i];
      const auto &edge = dag_.GetDAGEdge(edge_idx);
      const auto &parent_node = dag_.GetDAGNode(edge.GetParent());
      const auto &child_node = dag_.GetDAGNode(edge.GetChild());
      // Check if edge goes to root.
      if (dag_.IsNodeRoot(parent_node.Id())) {
        if (root_check == true) {
          return false;
        }
        root_check = true;
      }
      // Check if edge goes to leaf.
      if (dag_.IsNodeLeaf(child_node.Id())) {
        const auto taxon_id = child_node.Id();
        if (leaf_check.at(taxon_id) == true) {
        }
        leaf_check.at(taxon_id) = true;
      }
      // Update node map. If a node already has parent or child, it is an invalid tree.
      for (const auto node_id : {parent_node.Id(), child_node.Id()}) {
        if (nodemap_check.find(node_id) == nodemap_check.end()) {
          nodemap_check.insert({node_id, {false, false, false}});
        }
      }
      const size_t which_child =
          (edge.GetSubsplitClade() == SubsplitClade::Left) ? 0 : 1;
      if (nodemap_check[parent_node.Id()][which_child] == true) {
        return false;
      }
      nodemap_check[parent_node.Id()][which_child] = true;
      if (nodemap_check[child_node.Id()][2] == true) {
        return false;
      }
      nodemap_check[child_node.Id()][2] = true;
    }
    // Check if all nodes are fully connected.
    for (const auto &[node_id, connections] : nodemap_check) {
      // Check children.
      if (!(connections[0] || connections[1])) {
        // If one child is not connected, neither should be.
        if (connections[0] ^ connections[1]) {
          return false;
        }
        if (!dag_.IsNodeLeaf(node_id)) {
          return false;
        }
      }
      // Check parent.
      if (!connections[2]) {
        if (!dag_.IsNodeRoot(node_id)) {
          return false;
        }
      }
    }
    // Check if spans root and all leaf nodes.
    if (!root_check) {
      return false;
    }
    for (size_t i = 0; i < leaf_check.size(); i++) {
      if (!leaf_check[i]) {
        return false;
      }
    }

    return true;
  }

  // ** I/O

  // Output edge choice to string.
  static std::string EdgeChoiceToString(const EdgeChoice &edge_choice) {
    std::stringstream os;
    os << "{ ";
    os << "parent: " << edge_choice.parent_edge_id << ", ";
    os << "sister: " << edge_choice.sister_edge_id << ", ";
    os << "left_child: " << edge_choice.left_child_edge_id << ", ";
    os << "right_child: " << edge_choice.right_child_edge_id;
    os << " }";
    return os.str();
  }

  // Output edge choice vector to string.
  static std::string EdgeChoiceVectorToString(
      const EdgeChoiceVector &edge_choice_vector) {
    std::stringstream os;
    os << "[ " << std::endl;
    for (size_t i = 0; i < edge_choice_vector.size(); i++) {
      const auto &edge_choice = edge_choice_vector[i];
      os << "\t" << EdgeChoiceToString(edge_choice) << ", " << std::endl;
    }
    os << "]";
    return os.str();
  }

  // Output edge choice map to iostream.
  friend std::ostream &operator<<(std::ostream &os, const ChoiceMap &choice_map) {
    os << EdgeChoiceVectorToString(choice_map.edge_choice_vector_);
    return os;
  }

 private:
  // Un-owned reference DAG.
  GPDAG &dag_;
  // A vector that stores a map of each edge's best adjacent edges.
  EdgeChoiceVector edge_choice_vector_;
};
