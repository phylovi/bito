// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_UNROOTED_TREE_COLLECTION_HPP_
#define SRC_UNROOTED_TREE_COLLECTION_HPP_

#include "generic_tree_collection.hpp"
#include "unrooted_tree.hpp"

using UnrootedTreeCollection = GenericTreeCollection<UnrootedTree>;

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("UnrootedTreeCollection") {}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_ROOTED_TREE_COLLECTION_HPP_
