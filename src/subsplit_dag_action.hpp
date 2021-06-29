// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A class to hold actions that we can perform on the subsplit DAG.

#ifndef SRC_SUBSPLIT_DAG_ACTION_HPP_
#define SRC_SUBSPLIT_DAG_ACTION_HPP_

#include <tuple>

template <size_t I, typename... Args>
using TypeOf = std::tuple_element_t<I, std::tuple<Args...>>;

// An action to be performed as part of a traversal of the subsplit DAG.
template <typename... Args>
auto SubsplitDAGTraversalAction(Args&&... args) {
  struct Impl {
    // Applied just before visiting a node.
    TypeOf<0, Args...> BeforeNode;
    // Applied after visiting a node.
    TypeOf<1, Args...> AfterNode;
    // Applied before visiting the set of edges below a (node, clade) pair.
    TypeOf<2, Args...> BeforeNodeClade;
    // Applied for each edge.
    TypeOf<3, Args...> VisitEdge;
  };
  return Impl{std::forward<Args>(args)...};
}

// An action to be performed as part of a traversal of a Tidy subsplit DAG.
template <typename... Args>
auto TidySubsplitDAGTraversalAction(Args&&... args) {
  struct Impl {
    // Applied just before visiting a node.
    TypeOf<0, Args...> BeforeNode;
    // Applied after visiting a node.
    TypeOf<1, Args...> AfterNode;
    // Applied before visiting the set of edges below a (node, clade) pair.
    TypeOf<2, Args...> BeforeNodeClade;
    // Applied for each edge, and "dirties" all of the nodes above it.
    TypeOf<3, Args...> ModifyEdge;
    // Applying this function "cleans" the node just above an edge that has been
    // "dirtied" by ModifyEdge assuming the node just below the edge is clean.
    // (Traversals using this action ensure that children are visited before parents.)
    TypeOf<4, Args...> UpdateEdge;
  };
  return Impl{std::forward<Args>(args)...};
}

#endif  // SRC_SUBSPLIT_DAG_ACTION_HPP_
