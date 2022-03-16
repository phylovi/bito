// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The purpose of this code is documented in the header.

#include "reps_and_likelihoods.hpp"

int main(int argc, char *argv[]) {
  if (argc != 8) {
    std::cout
        << "We need exactly 7 arguments: fasta, rooted_nwk, credible_rooted_nwk,"
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

  GPInstance all_trees_gp_inst("mmapped_plv.data");
  RootedSBNInstance cred_r_inst("cred_trees");
  RootedSBNInstance pp_r_inst("pp_trees");
  all_trees_gp_inst.ReadNewickFile(rooted_nwk_path);
  all_trees_gp_inst.ReadFastaFile(fasta_path);
  all_trees_gp_inst.MakeEngine();
  all_trees_gp_inst.TakeFirstBranchLength();
  cred_r_inst.ReadNewickFile(credible_rooted_nwk_path);
  pp_r_inst.ReadNewickFile(pp_rooted_nwk_path);

  auto all_trees = all_trees_gp_inst.GenerateCompleteRootedTreeCollection();
  auto cred_trees = cred_r_inst.tree_collection_;
  auto pp_trees = pp_r_inst.tree_collection_;
  auto indexer = all_trees_gp_inst.GetDAG().BuildGPCSPIndexer();
  //std::vector<RootedIndexerRepresentation> all_representations;
  //std::vector<RootedIndexerRepresentation> cred_representations;
  //std::vector<RootedIndexerRepresentation> pp_representations;
  //LoadIndexerRepresentations(all_representations, all_trees, indexer);
  //LoadIndexerRepresentations(cred_representations, cred_trees, indexer);
  //LoadIndexerRepresentations(pp_representations, pp_trees, indexer);
  auto all_representations = GetIndexerRepresentations(all_trees, indexer);
  auto cred_representations = GetIndexerRepresentations(cred_trees, indexer);
  auto pp_representations = GetIndexerRepresentations(pp_trees, indexer);

  UnrootedSBNInstance ur_inst("charlie");
  ur_inst.ReadNewickFile(rooted_nwk_path);
  ur_inst.ReadFastaFile(fasta_path);
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  ur_inst.PrepareForPhyloLikelihood(simple_specification, thread_count, {}, true,
                                    all_trees.TreeCount());
  const auto log_likelihoods = ur_inst.UnrootedLogLikelihoods(all_trees);

  WriteTreesToFile(out_path, all_representations, log_likelihoods);
  WriteTreesToFile(credible_out_path, cred_representations);
  WriteTreesToFile(pp_out_path, pp_representations);
}


std::vector<RootedIndexerRepresentation> GetIndexerRepresentations(
    PreRootedTreeCollection &trees, BitsetSizeMap &indexer) {
  std::vector<RootedIndexerRepresentation> indexer_representations;
  for (const auto &tree : trees.Trees()) {
    indexer_representations.push_back(
        RootedSBNMaps::IndexerRepresentationOf(indexer, tree.Topology(), SIZE_MAX));
  }
  return indexer_representations;
}


void WriteTreesToFile(
    const std::string &out_path,
    const std::vector<RootedIndexerRepresentation> &representations,
    const std::vector<double> &log_likelihoods) {
  std::ofstream out_stream(out_path);
  out_stream << std::setprecision(12);
  const auto write_likelihood = !log_likelihoods.empty();
  for (size_t which_tree = 0; which_tree < representations.size(); which_tree++) {
    for (const auto &idx : representations.at(which_tree)) {
      out_stream << idx << ",";
    }
    if (write_likelihood) {
      out_stream << log_likelihoods.at(which_tree);
    }
    out_stream << "\n";
  }
  out_stream.close();
}

