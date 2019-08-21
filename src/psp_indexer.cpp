// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

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
  // Now onto the PCSSs.
  for (const auto& iter : in_indexer) {
    const auto& pcss = iter.first;
    // The first condition allows us to skip the rootsplits. We only want the
    // PCSSs here. The second condition is because the "primary" part of Primary
    // Subsplit Pair means that the parent split is a rootsplit.
    if (iter.second >= rootsplits.size() && pcss.PCSSIsRootsplit()) {
      SafeInsert(indexer_, pcss.PCSSWithoutParent(), index);
      index++;
    }
  }
  first_empty_index_ = index;
}

StringVector PSPIndexer::ToStringVector() const {
  std::vector<std::string> reversed_indexer(indexer_.size() + 1);
  for (const auto& iter : indexer_) {
    if (iter.second < after_rootsplits_index_) {
      reversed_indexer[iter.second] = iter.first.ToString();
    } else {
      reversed_indexer[iter.second] = iter.first.SubsplitToString();
    }
  }
  // Add an extra at the end for the split that doesn't exist.
  reversed_indexer[indexer_.size()] = "";
  return reversed_indexer;
}

SizeVectorVector PSPIndexer::RepresentationOf(
    const Node::NodePtr& topology) const {
  Assert(first_empty_index_ > 0, "This PSPIndexer is uninitialized.");
  SizeVector rootsplit_result(topology->Id(), first_empty_index_);
  SizeVector psp_result_down(topology->Id(), first_empty_index_);
  SizeVector psp_result_up(topology->Id(), first_empty_index_);
  const auto leaf_count = topology->LeafCount();
  Bitset rootsplit_bitset(leaf_count);
  Bitset psp_bitset(2 * leaf_count);
  auto rootsplit_index =
      [&rootsplit_bitset, &indexer = this->indexer_ ](const Node* node) {
    rootsplit_bitset.Zero();
    rootsplit_bitset.CopyFrom(node->Leaves(), 0, false);
    rootsplit_bitset.Minorize();
    return indexer.at(rootsplit_bitset);
  };
  // Here we use the terminology in the 2018 ICLR paper (screenshotted in
  // https://github.com/phylovi/libsbn/issues/95) looking at the right-hand case
  // in blue. The primary subsplit pair has Z_1 and Z_2 splitting apart Z. Here
  // we use analogous notation.
  auto psp_index = [&psp_bitset, &leaf_count, &indexer = this->indexer_ ](
      const Bitset& z1, const Bitset& z2, const Bitset& z) {
    psp_bitset.Zero();
    // Set Z in the middle location.
    psp_bitset.CopyFrom(z, 0, false);
    psp_bitset.CopyFrom(std::min(z1, z2), leaf_count, false);
    return indexer.at(psp_bitset);
  };
  topology->TriplePreOrder(
      // f_root
      [&rootsplit_result, &psp_result_up, &rootsplit_index, &psp_index](
          const Node* node0, const Node* node1, const Node* node2) {
        rootsplit_result[node0->Id()] = rootsplit_index(node0);
        psp_result_up[node0->Id()] =
            psp_index(node1->Leaves(), node2->Leaves(), ~node0->Leaves());
      },
      // f_internal
      [&rootsplit_result, &psp_result_up, &psp_result_down, &rootsplit_index,
       &psp_index](const Node* node, const Node* sister, const Node* parent) {
        rootsplit_result[node->Id()] = rootsplit_index(node);
        psp_result_up[node->Id()] =
            psp_index(~parent->Leaves(), sister->Leaves(), ~node->Leaves());
        psp_result_down[parent->Id()] =
            psp_index(node->Leaves(), sister->Leaves(), parent->Leaves());
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
