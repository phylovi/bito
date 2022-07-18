// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// This is for laying out phylo alphabets for

#include <bitset>

#include "sugar.hpp"

namespace Alphabet {
// Alphabet for DNA
enum class DNA { A, C, G, T, X /* unknown */ };
static inline const size_t DNACount = 5;

class DNAEnum : public EnumWrapper<DNA, char, DNACount, DNA::A, DNA::T> {
 public:
  static inline const Array<std::string> Labels = {{"A", "C", "G", "T"}};

  friend std::ostream &operator<<(std::ostream &os, const Type e) {
    os << "DNA::" << Labels[e];
    return os;
  };
};

// Alphabet for Protein / Amino Acids
enum class Amino : size_t {
  A,
  R,
  N,
  D,
  C,
  Q,
  G,
  E,
  H,
  I,
  L,
  K,
  M,
  F,
  P,
  S,
  T,
  W,
  Y,
  V,
  X /* unknown */
};
static inline const size_t AminoCount = 21;

class AminoEnum : public EnumWrapper<Amino, size_t, AminoCount, Amino::A, Amino::V> {
 public:
  static inline const Array<std::string> Labels = {{"A", "R", "N", "D", "C", "Q", "G",
                                                    "E", "H", "I", "K", "M", "F", "P",
                                                    "S", "T", "W", "Y", "V", "X"}};

  friend std::ostream &operator<<(std::ostream &os, const Type e) {
    os << "Amino::" << Labels[e];
    return os;
  };
};

}  // namespace Alphabet
