// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ROOTED_TREE_HPP_
#define SRC_ROOTED_TREE_HPP_

#include "tree.hpp"

class RootedTree : public Tree {
 public:
  typedef std::vector<RootedTree> RootedTreeVector;

  explicit RootedTree(const Tree& tree);

  bool operator==(const RootedTree& other) const;
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
