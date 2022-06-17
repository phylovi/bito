bool ChoiceMap::TreeMaskIsValid(const TreeMask &tree_mask, const bool is_quiet) const {
  auto tree_mask_ext = ExtractExpandedTreeMask(tree_mask);
  return TreeMaskIsValid(tree_mask_ext, is_quiet);
}

bool ChoiceMap::TreeMaskIsValid(ExpandedTreeMask &tree_mask_ext,
                                const bool is_quiet) const {
  std::ofstream devnull("/dev/null");
  std::ostream &os = (is_quiet ? devnull : std::cerr);

  auto PrintAdjacentNodes = [this, &tree_mask_ext, &os](const size_t node_id) {
    os << "{ node_id: " << node_id;
    os << ", parent: " << tree_mask_ext[node_id][AdjacentNode::Parent];
    os << ", left_child: " << tree_mask_ext[node_id][AdjacentNode::Parent];
    os << ", right_child: " << tree_mask_ext[node_id][AdjacentNode::Parent];
    os << "} " << std::endl;
  };

  // Check all root and leaf nodes exist
  if (tree_mask_ext.find(dag_.GetDAGRootNodeId()) == tree_mask_ext.end()) {
    os << "Invalid TreeMask: Does not span root node." << std::endl;
    return false;
  }
  for (size_t i = 0; i < dag_.TaxonCount(); i++) {
    if (tree_mask_ext.find(i) == tree_mask_ext.end()) {
      os << "Invalid TreeMask: Does not span all leaf nodes." << std::endl;
      return false;
    }
  }

  // Check all nodes have valid connections.
  for (const auto &[node_id, adj_node_ids] : tree_mask_ext) {
    if (dag_.IsNodeRoot(node_id)) {
      if ((adj_node_ids[AdjacentNode::Parent] != NoId) ||
          (adj_node_ids[AdjacentNode::LeftChild] == NoId) ||
          (adj_node_ids[AdjacentNode::RightChild] == NoId)) {
        os << "Invalid TreeMask: The root_node has invalid adjacent nodes."
           << std::endl;
        PrintAdjacentNodes(node_id);
        return false;
      }
    } else if (dag_.IsNodeLeaf(node_id)) {
      if ((adj_node_ids[AdjacentNode::Parent] == NoId) ||
          (adj_node_ids[AdjacentNode::LeftChild] != NoId) ||
          (adj_node_ids[AdjacentNode::RightChild] != NoId)) {
        os << "Invalid TreeMask: A leaf_node has invalid adjacent nodes." << std::endl;
        PrintAdjacentNodes(node_id);
        return false;
      }
    } else {
      if ((adj_node_ids[AdjacentNode::Parent] == NoId) ||
          (adj_node_ids[AdjacentNode::LeftChild] == NoId) ||
          (adj_node_ids[AdjacentNode::RightChild] == NoId)) {
        os << "Invalid TreeMask: An internal_node has invalid adjacent nodes."
           << std::endl;
        PrintAdjacentNodes(node_id);
        return false;
      }
    }
  }
  return true;
}
