// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "psp_indexer.hpp"
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
};

SizeVectorVector PSPIndexer::RepresentationOf(const Node::NodePtr& topology) {
  Assert(first_empty_index_ > 0, "This PSPIndexer is uninitialized.");
  SizeVector psp_result_down(topology->Id(), first_empty_index_);
  SizeVector psp_result_up(topology->Id(), first_empty_index_);
  const auto leaf_count = topology->LeafCount();
  Bitset bitset(2 * leaf_count);
  // Here we use the terminology in the 2018 ICLR paper (screenshotted in
  // https://github.com/phylovi/libsbn/issues/95) looking at the right-hand case
  // in blue. The primary subsplit pair has Z_1 and Z_2 splitting apart Z. Here
  // we use analogous notation, but for the corresponding edges of the tree.
  auto psp_index = [&bitset, &leaf_count, &indexer = this->indexer_ ](
      const Bitset& z1_bitset, const Node* z2, const Node* z, bool up) {
    bitset.Zero();
    // Set Z in the middle location.
    bitset.CopyFrom(z->Leaves(), 0, up);
    bitset.CopyFrom(std::min(z1_bitset, z2->Leaves()), leaf_count, false);
    return indexer.at(bitset);
  };
  std::cout << topology->Newick() << std::endl;
  topology->TriplePreOrder(
      // f_root
      [&psp_result_up, &psp_index](const Node* node0, const Node* node1,
                                   const Node* node2) {
        // std::cout << node0->Id() << ", true\n";
        psp_result_up[node0->Id()] =
            psp_index(node1->Leaves(), node2, node0, true);
      },
      // f_internal
      [&psp_result_up, &psp_result_down, &psp_index](
          const Node* node, const Node* sister, const Node* parent) {
        // std::cout << parent->Id() << ", false\n";
        // std::cout << node->Id() << ", true\n";
        psp_result_up[node->Id()] =
            psp_index(~parent->Leaves(), sister, node, true);
        psp_result_down[parent->Id()] =
            psp_index(node->Leaves(), sister, parent, false);
      });
  return {psp_result_up, psp_result_down};
}
