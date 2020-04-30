// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_instance.hpp"
#include "driver.hpp"

void GPInstance::ReadFastaFile(std::string fname) {
  alignment_ = Alignment::ReadFasta(fname);
}

void GPInstance::ReadNewickFile(std::string fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(fname));
}

void GPInstance::ReadNexusFile(std::string fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNexusFile(fname));
}

void GPInstance::CheckSequencesAndTreesLoaded() const {
  if (alignment_.SequenceCount() == 0) {
    Failwith(
        "Load an alignment into your SBNInstance on which you wish to "
        "calculate phylogenetic likelihoods.");
  }
  // if (TreeCount() == 0) {
  //   Failwith(
  //       "Load some trees into your SBNInstance on which you wish to "
  //       "calculate phylogenetic likelihoods.");
  // }
}

void GPInstance::MakeEngine() {
  CheckSequencesAndTreesLoaded();
  SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());
  engine_ = std::make_unique<GPEngine>(site_pattern);
}

GPEngine *GPInstance::GetEngine() const {
  if (engine_ != nullptr) {
    return engine_.get();
  }
  // else
  Failwith(
      "Engine not available. Call PrepareForPhyloLikelihood to make an "
      "engine for phylogenetic likelihood computation computation.");
}
