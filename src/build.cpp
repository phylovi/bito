// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "build.hpp"
#include <algorithm>
#include <cassert>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>

// Make a map from Tags to the bitset representing the leaves below the Tag.
TagBitsetMap TagLeafSetMapOf(Node::NodePtr topology) {
  TagBitsetMap map;
  auto leaf_count = topology->LeafCount();
  topology->PostOrder([&map, leaf_count](const Node* node) {
    Bitset bitset((size_t)leaf_count);
    if (node->IsLeaf()) {
      bitset.set(node->MaxLeafID());
    } else {
      // Take the union of the children below.
      for (const auto& child : node->Children()) {
        bitset |= map.at(child->Tag());
      }
    }
    assert(map.insert({node->Tag(), std::move(bitset)}).second);
  });
  return map;
}

// Make a map from Tags to the bitset representing the indices below the Tag.
SizeBitsetMap IndexIndexSetMapOf(Node::NodePtr topology) {
  SizeBitsetMap map;
  auto index_count = topology->Index() + 1;
  topology->PostOrder([&map, index_count](const Node* node) {
    Bitset bitset(static_cast<size_t>(index_count));
    if (node->Index() >= index_count) {
      std::cerr << "Malformed indices in TagIndexsetMapOf.\n";
      abort();
    }
    // Set the bit for the index of the current edge.
    bitset.set(node->Index());
    // Take the union of the children below.
    for (const auto& child : node->Children()) {
      bitset |= map.at(child->Index());
    }
    assert(map.insert({node->Index(), std::move(bitset)}).second);
  });
  return map;
}

BitsetUInt32Dict RootsplitCounterOf(const Node::TopologyCounter& topologies) {
  BitsetUInt32Dict rootsplit_counter(0);
  for (const auto& iter : topologies) {
    auto topology = iter.first;
    auto count = iter.second;
    auto tag_to_leafset = TagLeafSetMapOf(topology);
    auto Aux = [&rootsplit_counter, &tag_to_leafset, &count](const Node* n) {
      auto split = tag_to_leafset.at(n->Tag()).copy();
      split.Minorize();
      rootsplit_counter.increment(std::move(split), count);
    };
    for (const auto& child : topology->Children()) {
      child->PreOrder(Aux);
    }
  }
  return rootsplit_counter;
}

PCSSDict PCSSCounterOf(const Node::TopologyCounter& topologies) {
  PCSSDict pcss_dict;
  for (const auto& iter : topologies) {
    auto topology = iter.first;
    auto count = iter.second;
    auto tag_to_leafset = TagLeafSetMapOf(topology);
    auto leaf_count = topology->LeafCount();
    if (topology->Children().size() != 3) {
      std::cerr << "PCSSCounterOf was expecting a tree with a trifurcation at "
                   "the root!";
      abort();
    }
    topology->PCSSPreOrder([&pcss_dict, &tag_to_leafset, &count, &leaf_count](
                               const Node* sister_node, bool sister_direction,
                               const Node* focal_node, bool focal_direction,  //
                               const Node* child0_node,
                               bool child0_direction,  //
                               const Node* child1_node, bool child1_direction,
                               const Node*  // ignore virtual root clade
                           ) {
      Bitset parent(2 * leaf_count, false);
      // The first chunk is for the sister node.
      parent.CopyFrom(tag_to_leafset.at(sister_node->Tag()), 0,
                      sister_direction);
      // The second chunk is for the focal node.
      parent.CopyFrom(tag_to_leafset.at(focal_node->Tag()), leaf_count,
                      focal_direction);
      // Now we build the child bitset.
      auto child0 = tag_to_leafset.at(child0_node->Tag());
      if (child0_direction) {
        child0.flip();
      }
      auto child1 = tag_to_leafset.at(child1_node->Tag());
      if (child1_direction) {
        child1.flip();
      }
      auto child = std::min(child0, child1);
      // Insert the parent-child pair into the map.
      auto search = pcss_dict.find(parent);
      if (search == pcss_dict.end()) {
        // The first time we have seen this parent.
        BitsetUInt32Dict child_singleton(0);
        child_singleton.increment(std::move(child), count);
        assert(pcss_dict.insert({parent, child_singleton}).second);
      } else {
        search->second.increment(std::move(child), count);
      }
    });
  }
  return pcss_dict;
}

IndexerRepresentation IndexerRepresentationOf(const BitsetUInt32Map& indexer,
                                              const Node::NodePtr& topology) {
  auto tag_to_leafset = TagLeafSetMapOf(topology);
  auto leaf_count = topology->LeafCount();
  // First, the rootsplits.
  SizeVector rootsplit_result(topology->Index());
  topology->PreOrder([&topology, &rootsplit_result, &tag_to_leafset,
                      &indexer](const Node* node) {
    // Skip the root.
    if (node != topology.get()) {
      Bitset rootsplit = tag_to_leafset.at(node->Tag());
      rootsplit.Minorize();
      rootsplit_result[node->Index()] = indexer.at(rootsplit);
    }
  });
  // Next, the pcss_result.
  SizeVectorVector pcss_result(topology->Index());
  topology->PCSSPreOrder([&indexer, &tag_to_leafset, &leaf_count, &pcss_result,
                          &topology](
                             const Node* sister_node, bool sister_direction,
                             const Node* focal_node,
                             bool focal_direction,  //
                             const Node* child0_node,
                             bool child0_direction,  //
                             const Node* child1_node, bool child1_direction,
                             const Node* virtual_root_clade) {
    // Start by making the bitset representation of this PCSS.
    Bitset bitset(3 * leaf_count, false);
    bitset.CopyFrom(tag_to_leafset.at(sister_node->Tag()), 0, sister_direction);
    bitset.CopyFrom(tag_to_leafset.at(focal_node->Tag()), leaf_count,
                    focal_direction);
    auto child0_bitset = tag_to_leafset.at(child0_node->Tag());
    if (child0_direction) child0_bitset.flip();
    auto child1_bitset = tag_to_leafset.at(child1_node->Tag());
    if (child1_direction) child1_bitset.flip();
    bitset.CopyFrom(std::min(child0_bitset, child1_bitset), 2 * leaf_count,
                    false);
    auto indexer_position = indexer.at(bitset);
    const auto& focal_index = focal_node->Index();
    if (sister_node == focal_node) {
      // We are in the bidirectional edge situation.
      assert(focal_index < pcss_result.size());
      // Rooting at the present edge will indeed lead to the given PCSS.
      pcss_result[focal_index].push_back(indexer_position);
    } else {
      // The only time the virtual root clade should be nullptr should be when
      // sister_node == focal_node, but we check anyhow.
      assert(virtual_root_clade != nullptr);
      // Virtual-rooting on every edge in the sister will also lead to this
      // PCSS, because then the root will be "above" the PCSS.
      virtual_root_clade->ConditionalPreOrder([&pcss_result, &indexer_position,
                                               &sister_node, &focal_node,
                                               &topology](const Node* node) {
        if (node == sister_node || node == focal_node) {
          // Don't enter the sister or focal clades. This is only
          // activated in the second case on the bottom row of pcss.svg:
          // we want to add everything in the clade above the focal node,
          // but nothing else.
          return false;
        }  // else
        // Add all of the edges of the virtual rooting clade, except for the
        // root of the topology.
        if (node != topology.get()) {
          assert(node->Index() < pcss_result.size());
          pcss_result[node->Index()].push_back(indexer_position);
        }
        return true;
      });
    }
  });
  return std::pair<SizeVector, SizeVectorVector>(rootsplit_result, pcss_result);
}
