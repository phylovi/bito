// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// This is a rooted tree class that had the extra parameters required to do node height
// gradients. In fact, because RootedTree also has branch lengths (inherited from Tree)
// all of the extra state in this object is redundant other than the tip dates.
//
// In the terminology here, imagine that the tree is rooted at the top and hangs down.
// Node "heights" is how far we have to go back in time to that divergence event from
// the present.

#ifndef SRC_ROOTED_TREE_HPP_
#define SRC_ROOTED_TREE_HPP_

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

  TagDoubleMap TagDateMapOfDateVector(std::vector<double> leaf_date_vector);

  // TODO and make unit test.
  std::vector<double> height_ratios_;
  // The actual node heights.
  std::vector<double> node_heights_;
  // The lower bound for the height of this node. This is set to be the maximum of the
  // current height of the two children.
  std::vector<double> node_bounds_;
};

inline bool operator!=(const RootedTree& lhs, const RootedTree& rhs) {
  return !(lhs == rhs);
}

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("RootedTree") {
  // (0,(1,(2,3)))
  auto topology = Node::ExampleTopologies()[3];
  RootedTree tree(Tree(topology, {2., 1.5, 2., 1., 2.5, 2.5, 0.}));
  std::vector<double> date_vector({5., 3., 0., 1.});
  auto tag_date_map = tree.TagDateMapOfDateVector(date_vector);
  tree.InitializeParameters(tag_date_map);
  std::vector<double> correct_node_heights({5., 3., 0., 1., 2., 4.5, 7.});
  std::vector<double> correct_node_bounds({5., 3., 0., 1., 1., 3., 5.});
  for (size_t i = 0; i < correct_node_heights.size(); ++i) {
    CHECK_EQ(correct_node_heights[i], tree.node_heights_[i]);
    CHECK_EQ(correct_node_bounds[i], tree.node_bounds_[i]);
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_ROOTED_TREE_HPP_
