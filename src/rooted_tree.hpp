// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// This is a rooted tree class that had the extra parameters required to do node height
// gradients. In fact, because RootedTree also has branch lengths (inherited from Tree)
// all of the extra state in this object is redundant other than the tip dates.
//
// In the terminology here, imagine that the tree is rooted at the top and hangs down.
// The "height" of a node is how far we have to go back in time to that divergence event
// from the present.
//
// The most important parameterization here is in terms of node height ratios, which are
// of the form n/d, where
//
// n = time difference between this node's height and that of its earliest descendant E
// d = time difference between the parent's height and that of E.

#ifndef SRC_ROOTED_TREE_HPP_
#define SRC_ROOTED_TREE_HPP_

#include "eigen_sugar.hpp"
#include "tree.hpp"

class RootedTree : public Tree {
 public:
  typedef std::vector<RootedTree> RootedTreeVector;

  explicit RootedTree(const Tree& tree);

  bool operator==(const Tree& other) const = delete;
  bool operator==(const RootedTree& other) const;

  // Initialize the parameters of this tree using tip dates (and the branch lengths that
  // already exist in the tree). Note that these dates are the amount of time elapsed
  // between the sampling date and the present. Thus, older times have larger dates.
  void InitializeParameters(const TagDoubleMap& tag_date_map);

  // Set node_heights_ so that they have the given height ratios.
  void SetNodeHeightsViaHeightRatios(EigenConstVectorXdRef height_ratios);

  TagDoubleMap TagDateMapOfDateVector(EigenVectorXd leaf_date_vector);

  // This vector is of length equal to the number of internal nodes, and (except for the
  // last entry) has the node height ratios. The last entry is the root height.
  // The indexing is set up so that the ith entry has the node height ratio for the
  // (i+leaf_count)th node for all i except for the last.
  EigenVectorXd height_ratios_;
  // The actual node heights for all nodes.
  EigenVectorXd node_heights_;
  // The lower bound for the height of each node, which is the maximum of the tip dates
  // across all of the descendants of the node.
  EigenVectorXd node_bounds_;
  // The per-branch substitution rates.
  EigenVectorXd rates_;
  // Number of substitution rates (e.g. 1 rate for strict clock)
  size_t rate_count_;

  // The tree depicted in
  // https://github.com/phylovi/libsbn/issues/187#issuecomment-618421183
  static RootedTree Example();
};

inline bool operator!=(const RootedTree& lhs, const RootedTree& rhs) {
  return !(lhs == rhs);
}

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("RootedTree") {
  // To understand this test, please see
  // https://github.com/phylovi/libsbn/issues/187#issuecomment-618421183
  auto tree = RootedTree::Example();
  EigenVectorXd correct_height_ratios(3);
  correct_height_ratios << 1. / 3.5, 1.5 / 4., 7.;
  for (size_t i = 0; i < correct_height_ratios.size(); ++i) {
    CHECK_EQ(correct_height_ratios[i], tree.height_ratios_(i));
  }
  EigenVectorXd correct_node_heights(7);
  correct_node_heights << 5., 3., 0., 1., 2., 4.5, 7.;
  EigenVectorXd correct_node_bounds(7);
  correct_node_bounds << 5., 3., 0., 1., 1., 3., 5.;
  for (size_t i = 0; i < correct_node_heights.size(); ++i) {
    CHECK_EQ(correct_node_heights[i], tree.node_heights_(i));
    CHECK_EQ(correct_node_bounds[i], tree.node_bounds_(i));
  }
  // Test ratios to heights.
  const double arbitrary_dummy_number = -5.;
  std::fill(tree.LeafCount() + tree.node_heights_.begin(),  // First internal node.
            tree.node_heights_.end(), arbitrary_dummy_number);
  tree.SetNodeHeightsViaHeightRatios(correct_height_ratios);
  for (size_t i = 0; i < correct_node_heights.size(); ++i) {
    CHECK_EQ(correct_node_heights[i], tree.node_heights_(i));
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_ROOTED_TREE_HPP_
