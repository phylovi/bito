#ifndef SRC_BITSET_HPP_
#define SRC_BITSET_HPP_

#include <string>
#include <vector>

class Bitset {
 public:
  explicit Bitset(std::vector<bool> value);
  explicit Bitset(size_t n, bool initial_value = false);
  explicit Bitset(std::string);

  bool operator[](size_t i) const;
  void set(size_t i);
  void reset(size_t i);
  size_t size(void) const;
  size_t hash(void) const;

  bool operator==(const Bitset &bs) const;
  bool operator!=(const Bitset &bs) const;
  bool operator<(const Bitset &bs) const;
  bool operator<=(const Bitset &bs) const;
  bool operator>(const Bitset &bs) const;
  bool operator>=(const Bitset &bs) const;

  Bitset operator&(const Bitset &bs) const;
  Bitset operator|(const Bitset &bs) const;
  Bitset operator^(const Bitset &bs) const;
  Bitset operator~() const;

  std::string ToString();


 private:
  std::vector<bool> value_;
};


#endif  // SRC_BITSET_HPP_
