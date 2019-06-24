#ifndef RbBitset_H
#define RbBitset_H

#include <ostream>
#include <vector>

/**
 * RevBayes class for bit sets.
 *
 * This class is essentially a wrapper of a std::vector<bool> with some
 * additional functionality.
 *
 * @copyright Copyright 2009-
 * @author The RevBayes Development Core Team (Sebastian Hoehna)
 * @since 2016-09-06, version 1.0
 */
class Bitset {
 public:
  Bitset(void);
  Bitset(size_t n, bool def = false);
  Bitset(std::string);
  virtual ~Bitset(void) {}

  bool operator[](size_t i) const;

  bool operator==(const Bitset &bs) const;
  bool operator!=(const Bitset &bs) const;
  bool operator<(const Bitset &bs) const;
  bool operator<=(const Bitset &bs) const;
  bool operator>(const Bitset &bs) const;
  bool operator>=(const Bitset &bs) const;

  Bitset operator&(const Bitset &bs) const;
  Bitset operator|(const Bitset &bs) const;
  Bitset operator^(const Bitset &bs) const;
  Bitset &operator~();
  Bitset &operator&=(const Bitset &bs);
  Bitset &operator|=(const Bitset &bs);


  void clear(void);
  bool empty(void) const;
  void flip();
  void flip(size_t i);
  size_t getNumberSetBits(void) const;  //!< Get the number of bits set.
  size_t getFirstSetBit(void) const;    //!< Get the number of bits set.
  bool isSet(size_t i) const;
  void resize(size_t size);
  void set(size_t i);
  size_t size(void) const;
  void unset(size_t i);


 private:
  std::vector<bool> value_;
  size_t num_set_bits;


    };

    std::ostream &operator<<(std::ostream &o, const Bitset &x);


#endif
