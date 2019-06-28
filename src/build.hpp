#include <cassert>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "bitset.hpp"
#include "driver.hpp"
#include "tree.hpp"

typedef std::unordered_map<uint64_t, Bitset> TagToBitsetMap;
typedef std::unordered_set<std::string> BitsetSet;
// typedef std::unordered_set<Bitset, std::hash<Bitset>, std::equal_to<Bitset>>
// BitsetSet;
typedef std::unordered_map<Bitset, int> ParamIndexer;

// Using insert and at avoids needing to make a default constructor.
// https://stackoverflow.com/questions/17172080/insert-vs-emplace-vs-operator-in-c-map

TagToBitsetMap MakeTagToBitsetMap(Node::NodePtr t) {
  TagToBitsetMap m;
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

void PrintTagToBitsetMap(TagToBitsetMap m) {
  for (auto iter = m.begin(); iter != m.end(); ++iter) {
    std::cout << StringOfPackedInt(iter->first) << " "
              << iter->second.ToString() << std::endl;
  }
}

BitsetSet RootsplitSet(Node::NodePtr t) {
  BitsetSet s;
  // ParamIndexer indexer;
  auto m = MakeTagToBitsetMap(t);

  auto Aux = [&s, &m](Node* n) {
    Bitset x = m.at(n->Tag());
    s.insert(x.ToString());
    // indexer.insert(std::make_pair(x, 1));
  };

  for (auto child : t->Children()) {
    child->PreOrder(Aux);
  }
  return s;
}

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Build") {
  Driver driver;

  auto t = driver.ParseString("((0,1),(2,(3,4)));");

  auto m = MakeTagToBitsetMap(t);

  std::cout << t->Newick() << std::endl;
  PrintTagToBitsetMap(m);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
