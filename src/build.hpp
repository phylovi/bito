#include <unordered_map>
#include "bitset.hpp"
#include "tree.hpp"

// typedef std::unordered_map<uint64_t, Bitset> TagToBitsetMap;
// typedef std::unordered_map<Bitset, std::pair<int, int>> ParamIndexer;
//
// TagToBitsetMap MakeTagToBitsetMap(Node::NodePtr t) {
//   TagToBitsetMap m;
//   unsigned int leaf_count = t->LeafCount();
//   t->PreOrder([&m, leaf_count](Node* n) {
//     Bitset x(leaf_count, false);
//     for (auto child : n->Children()) {
//       Bitset::AndWith(x, m[child->Tag()]);
//       m[n->Tag()] = x;
//     }
//   });
//   return m;
// }
