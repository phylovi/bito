#include "unrooted_sbn_instance.hpp"
#include "rooted_sbn_instance.hpp"
#include "gp_instance.hpp"

// This is just a place to muck around, and check out performance.

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

int main(int argc, char *argv[]) {
  Assert(argc == 4,
         "We need exactly 4 arguments: fasta, unrooted_nwk, rooted_nwk, and out_path.");
  std::string fasta_path = argv[0];
  std::string unrooted_nwk_path = argv[1];
  std::string rooted_nwk_path = argv[2];
  std::string out_path = argv[3];

  GPInstance gp_inst("_ignore/mmapped_plv.data");
  gp_inst.ReadNewickFile(rooted_nwk_path);
  gp_inst.ReadFastaFile(fasta_path);
  gp_inst.MakeEngine();
  auto tree_collection = gp_inst.GenerateCompleteRootedTreeCollection();

  auto indexer = gp_inst.GetDAG().BuildGPCSPIndexer();

  std::vector<RootedIndexerRepresentation> indexer_representations;
  for (const auto &tree : tree_collection.Trees()) {
    indexer_representations.push_back(
        RootedSBNMaps::IndexerRepresentationOf(indexer, tree.Topology(), SIZE_MAX));
  }
  std::cout << indexer_representations << std::endl;
  UnrootedSBNInstance ur_inst("charlie");
  ur_inst.ReadNewickFile(unrooted_nwk_path);
  ur_inst.ReadFastaFile(fasta_path);
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  ur_inst.PrepareForPhyloLikelihood(simple_specification, 2, {}, true,
                                    tree_collection.TreeCount());
  const auto log_likelihoods = ur_inst.UnrootedLogLikelihoods(tree_collection);

  std::ofstream out_stream(out_path);
  for (size_t which_tree = 0; which_tree < tree_collection.TreeCount(); which_tree++) {
    for (const auto &idx : indexer_representations.at(which_tree)) {
      out_stream << idx << ",";
    }
    out_stream << log_likelihoods.at(which_tree) << "\n";
  }
  out_stream.close();
}
