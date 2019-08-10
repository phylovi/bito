// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_DEFAULT_DICT_HPP_
#define SRC_DEFAULT_DICT_HPP_

#include <cassert>
#include <iostream>
#include <unordered_map>
#include <utility>

template <class Key, class T>
class DefaultDict {
 private:
  const T default_value_;
  std::unordered_map<Key, T> map_;

 public:
  explicit DefaultDict(T default_value) : default_value_(default_value) {}

  size_t size() const { return map_.size(); }
  std::unordered_map<Key, T> Map() const { return map_; }
  // TODO
  // std::iterator<std::pair<Key, T>> begin() const { return map_.begin(); }

  T at(const Key &key) {
    auto search = map_.find(key);
    if (search == map_.end()) {
      return default_value_;
    }
    return search->second;
  }

  bool contains(const Key &key) { return (map_.find(key) != map_.end()); }

  void increment(const Key &key, const T &value) {
    auto search = map_.find(key);
    if (search == map_.end()) {
      SafeInsert(map_, key, value);
    } else {
      search->second += value;
    }
  }

  void print() {
    std::cout << "Default value: " << default_value_ << std::endl;
    for (const auto &iter : map_) {
      std::cout << std::to_string(iter.first) << " "
                << std::to_string(iter.second) << std::endl;
    }
  }
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
