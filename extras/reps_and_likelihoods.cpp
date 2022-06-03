// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The purpose of this code is documented in the header.

#include "reps_and_likelihoods.hpp"

int main(int argc, char *argv[]) {
  if (argc <= 4 || argc % 2 !=1) {
    std::cout
        << "We need at least 4 arguments: fasta, rooted_nwk, repr_out_path, "
        << "and nwk_out_path." 
        << std::endl
        << "Additional arguments must come in pairs: extra_rooted_nwk, "
        << "extra_repr_out_path." << std::endl;
    abort();
  }
  std::string fasta_path = argv[1];
  std::string rooted_nwk_path = argv[2];
  std::string out_path = argv[3];
  std::string nwk_out_path = argv[4];
  std::vector<std::string> extra_nwk_paths;
  std::vector<std::string> extra_out_paths;
  for (int arg_index = 5; arg_index < argc; arg_index+=2) {
    extra_nwk_paths.push_back(argv[arg_index]);
    extra_out_paths.push_back(argv[arg_index + 1]);
  }
  const auto extras_count = extra_nwk_paths.size();
  auto thread_count = std::thread::hardware_concurrency();

  GPInstance all_trees_gp_inst("mmapped_plv.data");
  all_trees_gp_inst.ReadNewickFile(rooted_nwk_path);
  all_trees_gp_inst.ReadFastaFile(fasta_path);
  all_trees_gp_inst.MakeEngine();
  all_trees_gp_inst.TakeFirstBranchLength();
  std::vector<RootedSBNInstance> extra_r_insts;
  for (size_t i = 0; i < extras_count; i++) {
    extra_r_insts.push_back(RootedSBNInstance("extra_trees" + std::to_string(i)));
    extra_r_insts.at(i).ReadNewickFile(extra_nwk_paths.at(i));
  }

  const auto taxa_order = all_trees_gp_inst.GetTaxonNames();
  for (const auto &r_inst : extra_r_insts) {
    if (r_inst.tree_collection_.TaxonNames() != taxa_order) {
      std::cout << "The first tree of each newick file must have the taxa appearing in "
                << "the same order. Insert a dummy tree if needed." << std::endl;
      abort();
    }
  }

  auto all_trees = all_trees_gp_inst.GenerateCompleteRootedTreeCollection();
  auto indexer = all_trees_gp_inst.GetDAG().BuildEdgeIndexer();
  auto all_representations = GetIndexerRepresentations(all_trees, indexer);
  UnrootedSBNInstance ur_inst("charlie");
  ur_inst.ReadNewickFile(rooted_nwk_path);
  ur_inst.ReadFastaFile(fasta_path);
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  ur_inst.PrepareForPhyloLikelihood(simple_specification, thread_count, {}, true,
                                    all_trees.TreeCount());
  const auto log_likelihoods = ur_inst.UnrootedLogLikelihoods(all_trees);
  WriteTreesToFile(out_path, all_representations, log_likelihoods);
  WriteNewickToFile(nwk_out_path, all_trees);
  
  for (size_t i = 0; i < extras_count; i++) {
    WriteTreesToFile(
        extra_out_paths.at(i),
        GetIndexerRepresentations(extra_r_insts.at(i).tree_collection_, indexer));
  }
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

void WriteTreesToFile(const std::string &out_path,
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

void WriteNewickToFile(const std::string &out_path, const RootedTreeCollection &trees) {
  const auto node_labels = trees.TagTaxonMap();
  std::ofstream out_stream(out_path);
  for (const auto &tree : trees.Trees()) {
    out_stream << tree.NewickTopology(node_labels) << std::endl;
  }
  out_stream.close();
}

