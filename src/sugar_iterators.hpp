// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Templates for building different types of iterators.

#include <iostream>
#include <vector>
#include <tuple>

// SuperIterator takes any number of iterators of the value_type and
// iterates through them one-after-another.
template <typename ValueType, typename... Iterators>
class SuperIterator {
 public:
  SuperIterator(Iterators... iterators) : iterators_(std::make_tuple(iterators...)) {}

  // Define the necessary iterator operations
  bool operator!=(const SuperIterator& other) const {
    return std::apply(
        [this](const auto&... iters) {
          return (... && (iters != std::get<0>(iterators_)));
        },
        iterators_);
  }

  void operator++() {
    std::apply([](auto&... iters) { (..., (++iters)); }, iterators_);
  }

  // Return the value_type directly since they are guaranteed to be the same
  ValueType operator*() const {
    return std::get<0>(
        std::apply([this](const auto&... iters) { return std::make_tuple(*iters...); },
                   iterators_));
  }

 private:
  std::tuple<Iterators...> iterators_;
};

// Template function to create the SuperIterator
template <typename ValueType, typename... Containers>
auto MakeSuperIterator(Containers&... containers) {
  return SuperIterator<ValueType, typename Containers::iterator...>(
      containers.begin()...);
}
