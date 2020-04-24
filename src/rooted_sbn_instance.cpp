// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "rooted_sbn_instance.hpp"

void RootedSBNInstance::ReadNewickFile(std::string fname) {
  UnrootedSBNInstance::ReadNewickFile(fname);
  rooted_tree_collection_ = RootedTreeCollection::OfTreeCollection(tree_collection_);
  rooted_tree_collection_.ParseDatesFromTaxonNames();
  rooted_tree_collection_.InitializeParameters();
}

void RootedSBNInstance::ReadNexusFile(std::string fname) {
  UnrootedSBNInstance::ReadNexusFile(fname);
  rooted_tree_collection_ = RootedTreeCollection::OfTreeCollection(tree_collection_);
  rooted_tree_collection_.ParseDatesFromTaxonNames();
  rooted_tree_collection_.InitializeParameters();
}

std::vector<double> RootedSBNInstance::LogLikelihoods() {
  return GetEngine()->LogLikelihoods(rooted_tree_collection_, phylo_model_params_,
                                     rescaling_);
}

std::vector<std::pair<double, std::vector<double>>>
RootedSBNInstance::BranchGradients() {
  return GetEngine()->BranchGradients(rooted_tree_collection_, phylo_model_params_,
                                      rescaling_);
}
