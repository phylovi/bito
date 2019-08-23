#include <stack>
#include "libsbn.hpp"

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

/*
template <typename T>
void PostOrder(Node::NodePtr topology, std::function<void(T)> f) {
  // The stack records the nodes and whether they have been visited or not.
  std::stack<std::pair<T, bool>> stack;
  stack.push({topology.get(), false});
  const Node* node;
  bool visited;
  while (stack.size()) {
    std::tie(node, visited) = stack.top();
    stack.pop();
    if (visited) {
      // If we've already visited this node then we are on our way back.
      f(node);
    } else {
      // If not then we need to push ourself back on the stack (noting that
      // we've been visited)...
      stack.push({node, true});
      // And all of our children, which have not.
      const auto& children = node->Children();
      for (auto iter = children.rbegin(); iter != children.rend(); ++iter) {
        stack.push({(*iter).get(), false});
      }
    }
  }
}
*/

int main() {
  uint32_t leaf_count = 100000;

  Node::NodePtr topology = Node::Ladder(leaf_count);

  // topology->Polish();

  std::vector<size_t> ids;
  ids.reserve(1 + 2 * leaf_count);

  auto t_start = now();
  for (int i = 0; i < 100; i++) {
    ids.clear();
    topology->PostOrder(
        [&ids](const Node* node) { ids.push_back(node->Id()); });

    // topology->PreOrder([&ids](const Node* node) { ids.push_back(node->Id());
    // });
  }

  std::chrono::duration<double> duration = now() - t_start;
  std::cout << "time: " << duration.count() << " seconds\n";

}
