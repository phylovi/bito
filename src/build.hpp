#include <cassert>
#include <unordered_map>
#include "bitset.hpp"
#include "tree.hpp"

typedef std::unordered_map<uint64_t, Bitset> TagToBitsetMap;
// typedef std::unordered_map<Bitset, std::pair<int, int>> ParamIndexer;

// Using insert and at avoids needing to make a default constructor.
// https://stackoverflow.com/questions/17172080/insert-vs-emplace-vs-operator-in-c-map

TagToBitsetMap MakeTagToBitsetMap(Node::NodePtr t) {
  TagToBitsetMap m;
  auto leaf_count = t->LeafCount();
  t->PreOrder([&m, leaf_count](Node* n) {
    Bitset x((size_t)leaf_count);
    for (auto child : n->Children()) {
      Bitset::AndWith(x, m.at(child->Tag()));
      assert(m.insert(std::make_pair(n->Tag(), x)).second);
    }
  });
  return m;
}
