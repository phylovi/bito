#include <cassert>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "bitset.hpp"
#include "driver.hpp"
#include "tree.hpp"

typedef std::unordered_map<uint64_t, Bitset> TagBitsetMap;
typedef std::unordered_set<Bitset> BitsetSet;
typedef std::unordered_map<Bitset, int> BitsetIndexer;
typedef std::unordered_map<Bitset, float> BitsetFloatMap;

// Using insert and at avoids needing to make a default constructor.
// https://stackoverflow.com/questions/17172080/insert-vs-emplace-vs-operator-in-c-map

TagBitsetMap TagBitsetMapOf(Node::NodePtr t) {
  TagBitsetMap m;
  auto leaf_count = t->LeafCount();
  t->PostOrder([&m, leaf_count](Node* n) {
    Bitset x((size_t)leaf_count);
    if (n->IsLeaf()) {
      x.set(n->MaxLeafID());
    } else {
      // Take the union of the children below.
      for (auto child : n->Children()) {
        x |= m.at(child->Tag());
      }
    }
    assert(m.insert(std::make_pair(n->Tag(), x)).second);
  });
  return m;
}

void PrintTagBitsetMap(TagBitsetMap m) {
  for (auto iter = m.begin(); iter != m.end(); ++iter) {
    std::cout << StringOfPackedInt(iter->first) << " "
              << iter->second.ToString() << std::endl;
  }
}

BitsetFloatMap RootsplitFrequencyOf(Node::NodePtr t) {
  BitsetFloatMap rootsplit_frequency;
  auto tag_to_bitset = TagBitsetMapOf(t);

  auto Aux = [&rootsplit_frequency, &tag_to_bitset](Node* n) {
    // TODO actually add frequency
    rootsplit_frequency.insert(std::make_pair(tag_to_bitset.at(n->Tag()), 1.));
  };

  for (auto child : t->Children()) {
    child->PreOrder(Aux);
  }
  return rootsplit_frequency;
}

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Build") {
  Driver driver;

  auto t = driver.ParseString("((0,1),(2,(3,4)));");

  auto m = TagBitsetMapOf(t);

  std::cout << t->Newick() << std::endl;
  // TODO Add an actual test.
  PrintTagBitsetMap(m);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
