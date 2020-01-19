#include "libsbn.hpp"

// This is just a place to muck around, and check out performance.

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

int main() {
  auto t_start = now();

  SBNInstance inst("charlie");
  std::cout << "reading trees...\n";
  inst.ReadNewickFile("data/DS1.100_topologies.nwk");
  // inst.ReadNewickFile("_ignore/makona_out.t");
  std::cout << "processing trees...\n";
  inst.ProcessLoadedTrees();
  std::cout << "training sbn...\n";
  auto score_history = inst.TrainExpectationMaximization(0.0001, 10000);
  EigenToCSV("_ignore/makona-alpha0.0001-loop10000-score.csv", score_history);
  std::chrono::duration<double> duration = now() - t_start;
  std::cout << "time to train SBN: " << duration.count() << " seconds\n";
}
