#include "unrooted_sbn_instance.hpp"
#include "rooted_sbn_instance.hpp"
#include "gp_instance.hpp"

#include <thread>

// This is just a place to muck around, and check out performance.

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

// StringDoubleMap StringDoubleMapOfStringDoubleVector(StringDoubleVector vect) {
//   StringDoubleMap m;
//   for (const auto& [str, x] : vect) {
//     SafeInsert(m, str, x);
//   }
//   return m;
// }



int main(int argc, char *argv[]) {
  //std::cout << argc << std::endl;
  if (argc != 8) {
    std::cout << "We need exactly 7 arguments: fasta, rooted_nwk, credible_rooted_nwk,"
                "pp_rooted_nwk, repr_out_path, credible_repr_out_path, and pp_repr_out_path"
              << std::endl;
    abort();
  }
  std::string fasta_path = argv[1];
  std::string rooted_nwk_path = argv[2];
  std::string credible_rooted_nwk_path = argv[3];
  std::string pp_rooted_nwk_path = argv[4];
  std::string out_path = argv[5];
  std::string credible_out_path = argv[6];
  std::string pp_out_path = argv[7];
  auto thread_count = std::thread::hardware_concurrency();

  GPInstance gp_inst("mmapped_plv.data");
  gp_inst.ReadNewickFile(rooted_nwk_path);
  gp_inst.ReadFastaFile(fasta_path);
  gp_inst.MakeEngine();
  gp_inst.TakeFirstBranchLength();
  auto tree_collection = gp_inst.GenerateCompleteRootedTreeCollection();

  auto indexer = gp_inst.GetDAG().BuildGPCSPIndexer();
  std::vector<RootedIndexerRepresentation> indexer_representations;
  for (const auto &tree : tree_collection.Trees()) {
    indexer_representations.push_back(
        RootedSBNMaps::IndexerRepresentationOf(indexer, tree.Topology(), SIZE_MAX));
  }
 
  UnrootedSBNInstance ur_inst("charlie");
  ur_inst.ReadNewickFile(rooted_nwk_path);
  ur_inst.ReadFastaFile(fasta_path);
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  ur_inst.PrepareForPhyloLikelihood(simple_specification, thread_count, {}, true, 
                                    tree_collection.TreeCount());
  const auto log_likelihoods = ur_inst.UnrootedLogLikelihoods(tree_collection);

  std::ofstream out_stream(out_path);
  out_stream << std::setprecision(12);
  for (size_t which_tree = 0; which_tree < tree_collection.TreeCount(); which_tree++) {
    for (const auto &idx : indexer_representations.at(which_tree)) {
      out_stream << idx << ",";
    }
    out_stream << log_likelihoods.at(which_tree) << "\n";
  }
  out_stream.close();

  //We also need the representations of the credible_set topologies.
  RootedSBNInstance r_inst("just_a_name");
  r_inst.ReadNewickFile(credible_rooted_nwk_path);
  tree_collection = r_inst.tree_collection_;

  std::vector<RootedIndexerRepresentation> cred_indexer_representations;
  for (const auto &tree : tree_collection.Trees()) {
    cred_indexer_representations.push_back(
        RootedSBNMaps::IndexerRepresentationOf(indexer, tree.Topology(), SIZE_MAX));
  }
  out_stream.open(credible_out_path);
  for (size_t which_tree = 0; which_tree < tree_collection.TreeCount(); which_tree++) {
    for (const auto &idx : cred_indexer_representations.at(which_tree)) {
      out_stream << idx << ",";
    }
    out_stream << std::endl;
  }
  out_stream.close();



  //We also need the representations of the topologies with pp values.
  RootedSBNInstance r_inst1("again");
  r_inst1.ReadNewickFile(pp_rooted_nwk_path);
  tree_collection = r_inst1.tree_collection_;

  std::vector<RootedIndexerRepresentation> cred_indexer_representations1;
  for (const auto &tree : tree_collection.Trees()) {
    cred_indexer_representations1.push_back(
        RootedSBNMaps::IndexerRepresentationOf(indexer, tree.Topology(), SIZE_MAX));
  }
  out_stream.open(pp_out_path);
  for (size_t which_tree = 0; which_tree < tree_collection.TreeCount(); which_tree++) {
    for (const auto &idx : cred_indexer_representations1.at(which_tree)) {
      out_stream << idx << ",";
    }
    out_stream << std::endl;
  }
  out_stream.close();




}
