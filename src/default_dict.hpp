#ifndef SRC_DEFAULT_DICT_HPP_
#define SRC_DEFAULT_DICT_HPP_

#include <iostream>
#include <unordered_map>

// A special class for
//
// Inheriting from STL containers is frowned upon.
// However, private inheritance mitigates the problems, in particular it's not
// possible to have a pointer to the base class.
// https://stackoverflow.com/a/19355666/467327
// https://stackoverflow.com/a/2035044/467327
// I would have just had a polymorphic function to do things, but we need to
// store a default value.

template <class Key, class T>
class DefaultDict : private std::unordered_map<Key, T> {
 private:
  const T default_value_;

 public:
  explicit DefaultDict(T default_value) : default_value_(default_value){};

  // Inherit some functions as-is.
  using std::unordered_map<Key, T>::size;
  using std::unordered_map<Key, T>::begin;
  using std::unordered_map<Key, T>::end;

  T at(const Key &key) {
    auto search = this->find(key);
    if (search == this->end()) {
      return default_value_;
    }
    return std::unordered_map<Key, T>::at(key);
  }

  void increment(const Key &key, const T &value) {
    auto search = this->find(key);
    if (search == this->end()) {
      // std::pair doesn't make a copy, in contrast to std::make_pair.
      assert(this->insert(std::pair<Key, T>(key, value)).second);
    } else {
      search->second += value;
    }
  };

  void print() {
    std::cout << "Default value: " << default_value_ << std::endl;
    for (auto iter = this->begin(); iter != this->end(); ++iter) {
      std::cout << std::to_string(iter->first) << " "
                << std::to_string(iter->second) << std::endl;
    }
  }
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("DefaultDict") {
  auto d = DefaultDict<int, int>(0);
  d.increment(0, 5);
  d.increment(0, 2);
  CHECK(d.at(0) == 7);
}

#endif  // DOCTEST_LIBRARY_INCLUDED


#endif  // SRC_DEFAULT_DICT_HPP_
