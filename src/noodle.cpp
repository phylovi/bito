#include "unrooted_sbn_instance.hpp"

// This is just a place to muck around, and check out performance.

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

void Print(Node::NodePtr node, size_t ident = 0) {
  for (size_t i = 0; i < ident; ++i) std::cout << "  ";
  std::cout << node->Id() << "\n";
  for (auto&& i : node->Children()) Print(i, ident + 1);
}

void testSampling() {
  const size_t out_of_sample_index = 99999999;

  UnrootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/five_taxon_unrooted.nwk");
  inst.ProcessLoadedTrees();
  inst.TrainSimpleAverage();
  // Count the frequencies of rooted trees in a file.
  size_t rooted_tree_count_from_file = 0;
  RootedIndexerRepresentationSizeDict counter_from_file(0);
  for (const auto &indexer_representation : inst.MakeIndexerRepresentations()) {
    RootedSBNMaps::IncrementRootedIndexerRepresentationSizeDict(counter_from_file,
                                                                indexer_representation);
    rooted_tree_count_from_file += indexer_representation.size();
  }
  // Count the frequencies of trees when we sample after training with
  // SimpleAverage.
  size_t sampled_tree_count = 1'000'000;
  RootedIndexerRepresentationSizeDict counter_from_sampling(0);
  ProgressBar progress_bar(sampled_tree_count / 1000);
  [[maybe_unused]] TopologySampler sampler;
  for (size_t sample_idx = 0; sample_idx < sampled_tree_count; ++sample_idx) {
    const auto rooted_topology = 
      //inst.SampleTopology(true);
      sampler.SampleTopology(inst);
    Print(rooted_topology);
    RootedSBNMaps::IncrementRootedIndexerRepresentationSizeDict(
        counter_from_sampling,
        RootedSBNMaps::IndexerRepresentationOf(inst.SBNSupport().Indexer(),
                                               rooted_topology, out_of_sample_index));
    if (sample_idx % 1000 == 0) {
      ++progress_bar;
      progress_bar.display();
    }
  }
  // These should be equal in the limit when we're training with SA.
  for (const auto &[key, _] : counter_from_file) {
    std::ignore = _;
    double observed =
        static_cast<double>(counter_from_sampling.at(key)) / sampled_tree_count;
    double expected =
        static_cast<double>(counter_from_file.at(key)) / rooted_tree_count_from_file;
    Assert(fabs(observed - expected) < 5e-3, "Sampling failed");
  }
  progress_bar.done();
}

int main() {

  testSampling();
  return 0;

  uint32_t leaf_count = 10000;

  Node::NodePtr topology = Node::Ladder(leaf_count);

  std::vector<size_t> ids;
  ids.reserve(1 + 2 * leaf_count);

  auto t_start = now();
  for (int i = 0; i < 100; i++) {
    ids.clear();
    topology->Preorder([&ids](const Node* node) { ids.push_back(node->Id()); });
  }

  std::chrono::duration<double> duration = now() - t_start;
  std::cout << "time: " << duration.count() << " seconds\n";
}
