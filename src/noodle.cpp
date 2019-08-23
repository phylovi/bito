#include "libsbn.hpp"

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

int main() {
  uint32_t leaf_count = 5000;

  auto t_start = now();

  Node::NodePtr topology = Node::Ladder(leaf_count);

  // std::cout << t->Newick() << std::endl;

  std::vector<size_t> ids;
  ids.reserve(1 + 2 * leaf_count);

  topology->Polish();

  /*
  TagSizeMap tag_id_map;
  size_t next_id = leaf_count;
  topology->MutablePostOrder([&tag_id_map, &next_id, &leaf_count](Node* node) {
    if (node->IsLeaf()) {
      node->id_ = node->MaxLeafID();
      node->leaves_ = Bitset::Singleton(leaf_count, node->id_);
    } else {
      node->id_ = next_id;
      next_id++;
      node->leaves_ = Node::LeavesOf(node->Children());
    }
    SafeInsert(tag_id_map, node->Tag(), node->id_);
  });

  TagSizeMap tag_id_map;
  size_t next_id = leaf_count;
  topology->PostOrder([&tag_id_map, &next_id, &leaf_count](const Node* node) {
    if (node->IsLeaf()) {
      Bitset::Singleton(leaf_count, node->Id());
    } else {
      next_id++;
      const auto& children = node->Children();
      Bitset leaves(children[0]->Leaves());
      for (size_t i = 1; i < children.size(); i++) {
        leaves |= children[i]->Leaves();
      }
    }
    SafeInsert(tag_id_map, node->Tag(), node->Id());
  });
  */
  std::chrono::duration<double> duration = now() - t_start;
  std::cout << "time: " << duration.count() << " seconds\n";
}
