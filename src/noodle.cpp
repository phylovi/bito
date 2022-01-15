#include "unrooted_sbn_instance.hpp"
#include "rooted_sbn_instance.hpp"
#include "gp_instance.hpp"

// This is just a place to muck around, and check out performance.

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

int main() {
  GPInstance gp_inst("_ignore/mmapped_plv.data");
  gp_inst.ReadNewickFile("data/five_taxon_rooted.nwk");
  gp_inst.ReadFastaFile("data/five_taxon.fasta");
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
  ur_inst.ReadNewickFile("data/five_taxon_rooted.nwk");
  ur_inst.ReadFastaFile("data/five_taxon.fasta");
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  ur_inst.PrepareForPhyloLikelihood(simple_specification, 2);
  std::cout << ur_inst.UnrootedLogLikelihoods(tree_collection);
}
