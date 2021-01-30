// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_DEFAULT_DICT_HPP_
#define SRC_DEFAULT_DICT_HPP_

#include <iostream>
#include <unordered_map>
#include <utility>

#include "sugar.hpp"

template <class Key, class T>
class DefaultDict {
 public:
  using iterator = typename std::unordered_map<Key, T>::iterator;
  using const_iterator = typename std::unordered_map<Key, T>::const_iterator;

  explicit DefaultDict(T default_value) : default_value_(default_value) {}

  size_t size() const { return map_.size(); }
  iterator begin() { return map_.begin(); }
  iterator end() { return map_.end(); }
  // Range-based for loops use const begin rather than cbegin.
  // https://stackoverflow.com/a/45732500/467327
  const_iterator begin() const { return map_.begin(); }
  const_iterator end() const { return map_.end(); }

  std::unordered_map<Key, T> Map() const { return map_; }

  T at(const Key &key) { return AtWithDefault(map_, key, default_value_); }

  bool contains(const Key &key) const { return (map_.find(key) != map_.end()); }

  void increment(const Key &key, const T &value) {
    auto search = map_.find(key);
    if (search == map_.end()) {
      SafeInsert(map_, key, value);
    } else {
      search->second += value;
    }
  }

  void increment(Key &&key, T value) {
    auto search = map_.find(key);
    if (search == map_.end()) {
      SafeInsert(map_, key, value);
    } else {
      search->second += value;
    }
  }

  void print() const {
    std::cout << "Default value: " << default_value_ << std::endl;
    for (const auto &iter : map_) {
      std::cout << std::to_string(iter.first) << " " << std::to_string(iter.second)
                << std::endl;
    }
  }

 private:
  const T default_value_;
  std::unordered_map<Key, T> map_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("DefaultDict") {
  auto d = DefaultDict<int, int>(0);
  CHECK_EQ(d.at(4), 0);
  d.increment(4, 5);
  CHECK_EQ(d.at(4), 5);
  d.increment(4, 2);
  CHECK_EQ(d.at(4), 7);
}

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_DEFAULT_DICT_HPP_
