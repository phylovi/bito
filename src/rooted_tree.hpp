// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ROOTED_TREE_HPP_
#define SRC_ROOTED_TREE_HPP_

#include "tree.hpp"

// TODO let's think about how we want the inheritance hierarchy to go. It makes a lot of
// sense to me to have UnrootedTree and RootedTree descend from a common, more general
// object. Then Tree will just be a general tree, which has branch lengths and a
// topology, and for which the rooting is unspecified.
// I assume then that we would have LogLikelihood, etc, for Unrooted and Rooted trees
// separately.
class RootedTree : public Tree {
 public:
  typedef std::vector<RootedTree> RootedTreeVector;

  explicit RootedTree(const Tree& tree);

  bool operator==(const RootedTree& other) const;

 private:
  // Contains ratios and the root height.
  // TODO I'd prefer a more descriptive name.
  std::vector<double> parameters_;

  // Node heights and node bounds are not parameters of the model but they will be
  // needed for calculating the gradient.
  std::vector<double> node_heights_;
  std::vector<double> node_bounds_;

  // TODO won't this be in the RootedTreeCollection?
  std::unordered_map<size_t, double> taxon_date_map_;
};

inline bool operator!=(const RootedTree& lhs, const RootedTree& rhs) {
  return !(lhs == rhs);
}

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("RootedTree") {
  // TODO
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_ROOTED_TREE_HPP_
