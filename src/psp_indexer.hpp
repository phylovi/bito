// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_PSP_INDEXER_HPP_
#define SRC_PSP_INDEXER_HPP_

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "sbn_maps.hpp"
#include "tree.hpp"

class PSPIndexer {
 public:
  PSPIndexer(BitsetVector rootsplits, BitsetUInt32Map in_indexer);

  SizeVectorVector RepresentationOf(const Node::NodePtr& topology);

 private:
  BitsetUInt32Map indexer_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("PSPIndexer") {}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_PSP_INDEXER_HPP_
