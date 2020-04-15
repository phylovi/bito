// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ROOTED_TREE_HPP_
#define SRC_ROOTED_TREE_HPP_

#include "tree.hpp"

class RootedTree : public Tree {
 public:
  typedef std::vector<RootedTree> RootedTreeVector;

  explicit RootedTree(const Tree& tree);

  bool operator==(const Tree& other) const = delete;
  bool operator==(const RootedTree& other) const;

  void InitializeParameters(const std::unordered_map<Tag, double>& tag_date_map);

  // Contains ratios and the root height.
  std::vector<double> parameters_;

  // Node heights and node bounds are not parameters of the model but they will be
  // needed for calculating the gradient.
  std::vector<double> node_heights_;
  std::vector<double> node_bounds_;
};

inline bool operator!=(const RootedTree& lhs, const RootedTree& rhs) {
  return !(lhs == rhs);
}

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("RootedTree") {
  // Issue 187: Unit tests here?
  // TODO or cut
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_ROOTED_TREE_HPP_
