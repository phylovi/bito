#include <string>

#include "unrooted_sbn_instance.hpp"
#include "gp_instance.hpp"

// This is just a place to muck around, and check out performance.

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

int main() {
//  uint32_t leaf_count = 10000;
//
//  Node::NodePtr topology = Node::Ladder(leaf_count);
//
//  std::vector<size_t> ids;
//  ids.reserve(1 + 2 * leaf_count);
//
//  auto t_start = now();
//  for (int i = 0; i < 100; i++) {
//    ids.clear();
//    topology->Preorder([&ids](const Node* node) { ids.push_back(node->Id()); });
//  }
//
//  std::chrono::duration<double> duration = now() - t_start;
//  std::cout << "time: " << duration.count() << " seconds\n";

    // Step over GP estimation to understand where the small branch length estimation takes place.
    std::string fasta_path = "/Users/sjun2/gp-experiments/tmp/blexp2/rep1/run.fasta";
    std::string rerooted_newick_path = "/Users/sjun2/gp-experiments/tmp/blexp2/rep1/mcmc_trees_rerooted.nwk";

    GPInstance inst("mmaped_plv.dat");
    inst.ReadFastaFile(fasta_path);
    inst.ReadNewickFile(rerooted_newick_path);
    inst.MakeEngine();

    inst.EstimateSBNParameters();
    inst.EstimateBranchLengths(1e-6, 100);


}
