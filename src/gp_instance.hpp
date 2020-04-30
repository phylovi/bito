// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_GP_INSTANCE_HPP_
#define SRC_GP_INSTANCE_HPP_

#include "gp_engine.hpp"
#include "rooted_tree_collection.hpp"
#include "site_pattern.hpp"

class GPInstance {
 public:
  GPInstance();

  void ReadFastaFile(std::string fname);
  void ReadNewickFile(std::string fname);
  void ReadNexusFile(std::string fname);

 private:
  Alignment alignment_;
  GPEngine engine_;
  RootedTreeCollection tree_collection_;

  void CheckSequencesAndTreesLoaded() const;
  void MakeEngine();
};

#endif  // SRC_GP_INSTANCE_HPP_
