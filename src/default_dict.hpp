// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_DEFAULT_DICT_HPP_
#define SRC_DEFAULT_DICT_HPP_

#include <cassert>
#include <iostream>
#include <unordered_map>
#include <utility>

// Inheriting from STL containers is frowned upon, but we do it here for fun!
// We could definitely implement using containment rather than inheritance if we
// wanted.
// Private inheritance mitigates the stated problems in our case, in particular
// it's not possible to have a pointer to the base class.
// https://stackoverflow.com/a/19355666/467327
// https://stackoverflow.com/a/2035044/467327

template <class Key, class T>
class DefaultDict : private std::unordered_map<Key, T> {
 private:
  const T default_value_;

 public:
  explicit DefaultDict(T default_value) : default_value_(default_value) {}

  // Inherit some functions as-is.
  using std::unordered_map<Key, T>::size;
  using std::unordered_map<Key, T>::begin;
  using std::unordered_map<Key, T>::end;

  T at(const Key &key) {
    auto search = this->find(key);
    if (search == this->end()) {
      return default_value_;
    }
    return search->second;
  }

  bool contains(const Key &key) { return (this->find(key) != this->end()); }

  void increment(const Key &key, const T &value) {
    auto search = this->find(key);
    if (search == this->end()) {
      SafeInsert(*this, key, value);
    } else {
      search->second += value;
    }
  }

  void print() {
    std::cout << "Default value: " << default_value_ << std::endl;
    for (const auto &iter : *this) {
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
