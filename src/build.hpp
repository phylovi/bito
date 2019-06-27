#include <cassert>
#include <unordered_map>
#include "bitset.hpp"
#include "driver.hpp"
#include "tree.hpp"

typedef std::unordered_map<uint64_t, Bitset> TagToBitsetMap;
// typedef std::unordered_map<Bitset, std::pair<int, int>> ParamIndexer;

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

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Build") {
  Driver driver;

  auto t = driver.ParseString("((0,1),2);");

  auto m = MakeTagToBitsetMap(t);

  std::cout << t->Newick() << std::endl;
  PrintTagToBitsetMap(m);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
