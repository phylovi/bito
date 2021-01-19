// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A class to hold actions that we can perform on the subsplit DAG.

#ifndef SRC_SUBSPLIT_DAG_ACTION_HPP_
#define SRC_SUBSPLIT_DAG_ACTION_HPP_

#include <tuple>

template <size_t I, typename... Args>
// std::get returns by reference, but adding decay_t will provide the actual type.
using TypeOf = std::decay_t<decltype(std::get<I>(std::declval<std::tuple<Args...>>()))>;

template <typename... Args>
auto SubsplitDAGAction(Args&&... args) {
  struct Impl {
    TypeOf<0, Args...> BeforeNode;
    TypeOf<1, Args...> AfterNode;
    TypeOf<2, Args...> BeforeNodeClade;
    TypeOf<3, Args...> AfterNodeClade;
    TypeOf<4, Args...> VisitEdge;
  };
  return Impl{std::forward<Args>(args)...};
}

#endif  // SRC_SUBSPLIT_DAG_ACTION_HPP_
