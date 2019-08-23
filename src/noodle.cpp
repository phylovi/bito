#include "libsbn.hpp"
int main() {
  uint32_t leaf_count = 10000;

  auto now = std::chrono::high_resolution_clock::now;
  auto t_start = now();

  Node::NodePtr t = Node::Ladder(leaf_count);

  // std::cout << t->Newick() << std::endl;

  std::vector<size_t> ids;
  ids.reserve(1 + 2 * leaf_count);

  t->PostOrder([&ids](const Node* node) { ids.push_back(node->Id()); });
  std::chrono::duration<double> duration = now() - t_start;
  std::cout << "time: " << duration.count() << " seconds\n";
}
