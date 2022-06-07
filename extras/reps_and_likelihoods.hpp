// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The purpose of this executable is to determine likelihoods and PCSP indexer
// representations of various trees as part of a proof of concept for nni sdag support.
// Specifically, this takes in a fasta file and three newick files of rooted trees with
// branch lengths. A subsplit DAG is built on the first set of trees and this is used
// for the indexer representations of all three sets of trees. When a subsplit appears
// in the second or third set of trees, but not the subsplit DAG, the subsplit has
// SIZE_MAX as the index.

#pragma once

#include "unrooted_sbn_instance.hpp"
#include "rooted_sbn_instance.hpp"
#include "gp_instance.hpp"
#include <thread>

std::vector<RootedIndexerRepresentation> GetIndexerRepresentations(
    PreRootedTreeCollection &trees, BitsetSizeMap &indexer);

void WriteTreesToFile(const std::string &out_path,
                      const std::vector<RootedIndexerRepresentation> &representations,
                      const std::vector<double> &log_likelihoods = {});

void WriteNewickToFile(const std::string &out_path, const RootedTreeCollection &trees);
