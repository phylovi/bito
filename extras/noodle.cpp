#include "gp_instance.hpp"
#include "rooted_sbn_instance.hpp"
#include "unrooted_sbn_instance.hpp"
#include "stopwatch.hpp"

// This is just a place to muck around, and check out performance.

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

int main() {
  uint32_t leaf_count = 10000;

  Node::NodePtr topology = Node::Ladder(leaf_count);

  std::vector<size_t> ids;
  ids.reserve(1 + 2 * leaf_count);

  Stopwatch timer;
  auto t_start = now();
  timer.Start();
  for (int i = 0; i < 100; i++) {
    ids.clear();
    topology->Preorder([&ids](const Node* node) { ids.push_back(node->Id()); });
  }

  double watch_duration = timer.Stop();
  std::chrono::duration<double> duration = now() - t_start;
  std::cout << "time: " << duration.count() << " seconds\n";
}

void MyTest() {
  Stopwatch timer;
  timer.Start();

  timer.Lap();

  double time = timer.Stop();
  DoubleVector laps = timer.GetLaps();
}
