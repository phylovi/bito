#include <stack>
#include "libsbn.hpp"

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

int main() {
  uint32_t leaf_count = 10000;

  Node::NodePtr topology = Node::Ladder(leaf_count);

  std::vector<size_t> ids;
  ids.reserve(1 + 2 * leaf_count);

  auto t_start = now();
  for (int i = 0; i < 1000; i++) {
    ids.clear();
    topology->PreOrder([&ids](const Node* node) { ids.push_back(node->Id()); });
  }

  std::chrono::duration<double> duration = now() - t_start;
  std::cout << "time: " << duration.count() << " seconds\n";

}
