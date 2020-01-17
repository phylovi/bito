#include "libsbn.hpp"

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

  auto t_start = now();
  for (int i = 0; i < 100; i++) {
    ids.clear();
    topology->PreOrder([&ids](const Node* node) { ids.push_back(node->Id()); });
  }

  std::chrono::duration<double> duration = now() - t_start;
  std::cout << "time: " << duration.count() << " seconds\n";

  // Playing with counting.
  SBNInstance inst("charlie");
  inst.ReadNewickFile("data/five_taxon.nwk");
  // inst.ReadNewickFile("data/DS1.100_topologies.nwk");
  inst.ProcessLoadedTrees();
  inst.TrainSimpleAverage();

  RootedIndexerRepresentationSizeDict counter_from_file(0);
  for (const auto& indexer_representation : inst.MakeIndexerRepresentations()) {
    SBNMaps::IncrementRootedIndexerRepresentationSizeDict(counter_from_file,
                                                          indexer_representation);
  }
  std::cout << counter_from_file << std::endl;

  RootedIndexerRepresentationSizeDict counter_from_sampling(0);
  for (size_t sample_idx = 0; sample_idx < 10000; ++sample_idx) {
    const auto topology = inst.SampleTopology(true);
    SBNMaps::IncrementRootedIndexerRepresentationSizeDict(
        counter_from_sampling,
        SBNMaps::RootedIndexerRepresentationOf(inst.indexer_, topology, 99999999));
  }
  std::cout << counter_from_sampling << std::endl;
}
