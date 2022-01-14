// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "psp_indexer.hpp"

#include <algorithm>

#include "sugar.hpp"

PSPIndexer::PSPIndexer(BitsetVector rootsplits, BitsetSizeMap in_indexer) {
  size_t index = 0;
  // First the rootsplits.
  for (const auto& rootsplit : rootsplits) {
    SafeInsert(indexer_, rootsplit, index);
    index++;
  }
  after_rootsplits_index_ = index;
  // Now onto the PCSPs.
  for (const auto& iter : in_indexer) {
    const auto& pcsp = iter.first;
    // The first condition allows us to skip the rootsplits. We only want the
    // PCSPs here. The second condition is because the "primary" part of Primary
    // Subsplit Pair means that the parent split is a rootsplit.
    if (iter.second >= rootsplits.size() && pcsp.PCSPIsParentRootsplit()) {
      SafeInsert(indexer_, pcsp.PCSPGetChildSubsplit(), index);
      index++;
    }
  }
  first_empty_index_ = index;
}

StringVector PSPIndexer::ToStringVector() const {
  std::vector<std::string> reversed_indexer(indexer_.size() + 1);
  for (const auto& iter : indexer_) {
    reversed_indexer[iter.second] = iter.first.SubsplitToString();
  }
  // Add an extra entry at the end for the split that doesn't exist.
  reversed_indexer[indexer_.size()] = "";
  return reversed_indexer;
}

SizeVectorVector PSPIndexer::RepresentationOf(const Node::NodePtr& topology) const {
  Assert(first_empty_index_ > 0, "This PSPIndexer is uninitialized.");
  SizeVector rootsplit_result(topology->Id(), first_empty_index_);
  SizeVector psp_result_down(topology->Id(), first_empty_index_);
  SizeVector psp_result_up(topology->Id(), first_empty_index_);
  auto rootsplit_index = [&indexer = this->indexer_](const Node* node) {
    return indexer.at(Bitset::RootsplitOfHalf(node->Leaves()));
  };
  // Here we use the terminology in the 2019 ICLR paper (screenshotted in
  // https://github.com/phylovi/bito/issues/95) looking at the right-hand case
  // in blue. The primary subsplit pair has Z_1 and Z_2 splitting apart Z. Here
  // we use analogous notation.
  auto psp_index = [&indexer = this->indexer_](const Bitset& z1, const Bitset& z2) {
    return indexer.at(Bitset::Subsplit(z1, z2));
  };
  topology->TriplePreorder(
      // f_rootsplit
      [&rootsplit_result, &psp_result_up, &rootsplit_index, &psp_index](
          const Node* node0, const Node* node1, const Node* node2) {
        rootsplit_result[node0->Id()] = rootsplit_index(node0);
        psp_result_up[node0->Id()] = psp_index(node1->Leaves(), node2->Leaves());
      },
      // f_internal
      [&rootsplit_result, &psp_result_up, &psp_result_down, &rootsplit_index,
       &psp_index](const Node* node, const Node* sister, const Node* parent) {
        rootsplit_result[node->Id()] = rootsplit_index(node);
        psp_result_up[node->Id()] = psp_index(~parent->Leaves(), sister->Leaves());
        psp_result_down[parent->Id()] = psp_index(node->Leaves(), sister->Leaves());
      });
  return {rootsplit_result, psp_result_down, psp_result_up};
}

StringVectorVector PSPIndexer::StringRepresentationOf(
    const Node::NodePtr& topology) const {
  StringVector reversed_indexer = ToStringVector();
  StringVectorVector result;
  for (const auto& partial_representation : RepresentationOf(topology)) {
    StringVector str_partial_representation;
    for (const auto& index : partial_representation) {
      str_partial_representation.push_back(reversed_indexer.at(index));
    }
    result.push_back(std::move(str_partial_representation));
  }
  return result;
}

DoubleVectorVector PSPIndexer::SplitLengths(
    const UnrootedTreeCollection& tree_collection) const {
  DoubleVectorVector result(after_rootsplits_index_);
  auto tree_count = tree_collection.TreeCount();
  for (size_t tree_index = 0; tree_index < tree_count; tree_index++) {
    const auto& tree = tree_collection.GetTree(tree_index);
    // The 0th part of the PSP representation is the rootsplit vector.
    auto split_indices = RepresentationOf(tree.Topology())[0];
    const auto branch_lengths = tree.BranchLengths();
    for (size_t edge_index = 0; edge_index < split_indices.size(); edge_index++) {
      result[split_indices[edge_index]].push_back(branch_lengths[edge_index]);
    }
  }
  return result;
}
