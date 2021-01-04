// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "sbn_maps.hpp"

#include <algorithm>
#include <memory>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>

SizeBitsetMap SBNMaps::IdIdSetMapOf(const Node::NodePtr& topology) {
  SizeBitsetMap map;
  auto id_count = topology->Id() + 1;
  topology->PostOrder([&map, id_count](const Node* node) {
    Bitset bitset(static_cast<size_t>(id_count));
    Assert(node->Id() < id_count, "Malformed ids in IdIdSetMapOf.");
    // Set the bit for the id of the current edge.
    bitset.set(node->Id());
    // Take the union of the children below.
    for (const auto& child : node->Children()) {
      bitset |= map.at(child->Id());
    }
    SafeInsert(map, node->Id(), std::move(bitset));
  });
  return map;
}

SizeVector SBNMaps::SplitIndicesOf(const BitsetSizeMap& indexer,
                                   const Node::NodePtr& topology) {
  SizeVector split_result(topology->Id());
  topology->PreOrder([&topology, &split_result, &indexer](const Node* node) {
    // Skip the root.
    if (node != topology.get()) {
      Bitset rootsplit = node->Leaves();
      rootsplit.Minorize();
      split_result[node->Id()] = indexer.at(rootsplit);
    }
  });
  return split_result;
}

StringPCSPMap SBNMaps::StringPCSPMapOf(PCSPCounter d) {
  StringPCSPMap d_str;
  for (const auto& [parent, child_dict] : d) {
    d_str[parent.ToString()] = StringifyMap(child_dict.Map());
  }
  return d_str;
}

StringDoubleVector SBNMaps::StringDoubleVectorOf(BitsetDoubleMap m) {
  StringDoubleVector result;
  result.reserve(m.size());
  for (const auto& [k, v] : m) {
    result.push_back({k.ToString(), v});
  }
  std::sort(result.begin(), result.end());
  return result;
}

Bitset SBNMaps::PCSPBitsetOf(const size_t leaf_count,  //
                             const Node* sister_node, bool sister_direction,
                             const Node* focal_node, bool focal_direction,
                             const Node* child0_node, bool child0_direction,
                             const Node* child1_node, bool child1_direction) {
  Bitset bitset(3 * leaf_count, false);
  bitset.CopyFrom(sister_node->Leaves(), 0, sister_direction);
  bitset.CopyFrom(focal_node->Leaves(), leaf_count, focal_direction);
  auto child0_bitset = child0_node->Leaves();
  if (child0_direction) {
    child0_bitset.flip();
  }
  auto child1_bitset = child1_node->Leaves();
  if (child1_direction) {
    child1_bitset.flip();
  }
  bitset.CopyFrom(std::min(child0_bitset, child1_bitset), 2 * leaf_count, false);
  return bitset;
}

IndexerBundle SBNMaps::BuildIndexerBundle(const BitsetSizeDict& rootsplit_counter,
                                          const PCSPCounter& pcsp_counter) {
  BitsetVector rootsplits;
  BitsetSizeMap indexer;
  SizeBitsetMap index_to_child;
  BitsetSizePairMap parent_to_range;
  size_t index = 0;
  // Start by adding the rootsplits.
  for (const auto& iter : rootsplit_counter) {
    SafeInsert(indexer, iter.first, index);
    rootsplits.push_back(iter.first);
    index++;
  }
  // Now add the PCSPs.
  for (const auto& [parent, child_counter] : pcsp_counter) {
    SafeInsert(parent_to_range, parent, {index, index + child_counter.size()});
    for (const auto& child_iter : child_counter) {
      const auto& child = child_iter.first;
      SafeInsert(indexer, parent + child, index);
      SafeInsert(index_to_child, index, Bitset::ChildSubsplit(parent, child));
      index++;
    }
  }
  return {rootsplits, indexer, index_to_child, parent_to_range, index};
}

BitsetSizeDict UnrootedSBNMaps::RootsplitCounterOf(
    const Node::TopologyCounter& topologies) {
  BitsetSizeDict rootsplit_counter(0);
  for (const auto& [topology, topology_count] : topologies) {
    auto Aux = [&rootsplit_counter, &topology_count = topology_count](const Node* n) {
      auto split = n->Leaves();
      split.Minorize();
      rootsplit_counter.increment(std::move(split), topology_count);
    };
    for (const auto& child : topology->Children()) {
      child->PreOrder(Aux);
    }
  }
  return rootsplit_counter;
}

// See functions below or the comments above the definition of UnrootedPCSPFun to
// understand the collection of arguments starting with `sister_node`.
void AddToPCSPCounter(PCSPCounter& pcsp_dict, const size_t topology_count,
                      const size_t leaf_count, const Node* sister_node,
                      bool sister_direction, const Node* focal_node,
                      bool focal_direction, const Node* child0_node,
                      bool child0_direction, const Node* child1_node,
                      bool child1_direction) {
  Bitset parent(2 * leaf_count, false);
  // The first chunk is for the sister node.
  parent.CopyFrom(sister_node->Leaves(), 0, sister_direction);
  // The second chunk is for the focal node.
  parent.CopyFrom(focal_node->Leaves(), leaf_count, focal_direction);
  // Now we build the child bitset.
  auto child0 = child0_node->Leaves();
  if (child0_direction) {
    child0.flip();
  }
  auto child1 = child1_node->Leaves();
  if (child1_direction) {
    child1.flip();
  }
  auto child = std::min(child0, child1);
  // Insert the parent-child pair into the map.
  auto search = pcsp_dict.find(parent);
  if (search == pcsp_dict.end()) {
    // The first time we have seen this parent.
    BitsetSizeDict child_singleton(0);
    child_singleton.increment(std::move(child), topology_count);
    SafeInsert(pcsp_dict, std::move(parent), std::move(child_singleton));
  } else {
    search->second.increment(std::move(child), topology_count);
  }
}

PCSPCounter UnrootedSBNMaps::PCSPCounterOf(const Node::TopologyCounter& topologies) {
  PCSPCounter pcsp_dict;
  for (const auto& [topology, topology_count] : topologies) {
    auto leaf_count = topology->LeafCount();
    Assert(topology->Children().size() == 3,
           "UnrootedSBNMaps::PCSPCounterOf was expecting a tree with a trifurcation at "
           "the root!");
    topology->UnrootedPCSPPreOrder(
        [&pcsp_dict, &topology_count = topology_count, &leaf_count](
            const Node* sister_node, bool sister_direction, const Node* focal_node,
            bool focal_direction,  //
            const Node* child0_node,
            bool child0_direction,  //
            const Node* child1_node, bool child1_direction,
            const Node*  // ignore virtual root clade
        ) {
          AddToPCSPCounter(pcsp_dict, topology_count, leaf_count, sister_node,
                           sister_direction, focal_node, focal_direction, child0_node,
                           child0_direction, child1_node, child1_direction);
        });
  }
  return pcsp_dict;
}

// Return the rootsplit of a rooted bifurcating topology.
Bitset Rootsplit(const Node* rooted_topology) {
  Assert(rooted_topology->Children().size() == 2,
         "Rootsplit expects a bifurcating tree.");
  Bitset subsplit = rooted_topology->Children()[0]->Leaves();
  subsplit.Minorize();
  return subsplit;
}

UnrootedIndexerRepresentation UnrootedSBNMaps::IndexerRepresentationOf(
    const BitsetSizeMap& indexer, const Node::NodePtr& topology,
    const size_t default_index) {
  const auto leaf_count = topology->LeafCount();
  // First, the rootsplits.
  SizeVector rootsplit_result = SBNMaps::SplitIndicesOf(indexer, topology);
  // We initialize each vector with the rootsplit index.
  SizeVectorVector result(topology->Id());
  std::transform(rootsplit_result.begin(), rootsplit_result.end(), result.begin(),
                 [&topology](const auto rootsplit) {
                   SizeVector v = {rootsplit};
                   // The number of PCSPs is less than number of internal nodes/2.
                   v.reserve(topology->Id() / 2);
                   return v;
                 });
  // Now we append the PCSPs.
  topology->UnrootedPCSPPreOrder(
      [&indexer, &default_index, &leaf_count, &result, &topology](
          const Node* sister_node, bool sister_direction, const Node* focal_node,
          bool focal_direction, const Node* child0_node, bool child0_direction,
          const Node* child1_node, bool child1_direction,
          const Node* virtual_root_clade) {
        const auto bitset = SBNMaps::PCSPBitsetOf(
            leaf_count, sister_node, sister_direction, focal_node, focal_direction,
            child0_node, child0_direction, child1_node, child1_direction);
        const auto indexer_position = AtWithDefault(indexer, bitset, default_index);
        const auto& focal_index = focal_node->Id();
        if (sister_node == focal_node) {
          // We are in the bidirectional edge situation.
          Assert(focal_index < result.size(), "focal_index out of range.");
          // Rooting at the present edge will indeed lead to the given PCSP.
          result[focal_index].push_back(indexer_position);
        } else {
          // The only time the virtual root clade should be nullptr should be when
          // sister_node == focal_node, but we check anyhow.
          Assert(virtual_root_clade != nullptr, "virtual_root_clade is null.");
          // Virtual-rooting on every edge in the virtual rooting clade will also
          // lead to this PCSP, because then the root will be "above" the PCSP.
          virtual_root_clade->ConditionalPreOrder([&result, &indexer_position,
                                                   &sister_node, &focal_node,
                                                   &topology](const Node* node) {
            if (node == sister_node || node == focal_node) {
              // Don't enter the sister or focal clades. This is only
              // activated in the second case on the bottom row of pcsp.svg:
              // we want to add everything in the clade above the focal node,
              // but nothing else.
              return false;
            }  // else
            // Add all of the edges of the virtual rooting clade, except for the
            // root of the topology.
            if (node != topology.get()) {
              Assert(node->Id() < result.size(), "node's root Id is out of range.");
              result[node->Id()].push_back(indexer_position);
            }
            return true;
          });
        }
      });
  return result;
}

UnrootedIndexerRepresentationCounter UnrootedSBNMaps::IndexerRepresentationCounterOf(
    const BitsetSizeMap& indexer, const Node::TopologyCounter& topology_counter,
    const size_t default_index) {
  UnrootedIndexerRepresentationCounter counter;
  counter.reserve(topology_counter.size());
  for (const auto& [topology, topology_count] : topology_counter) {
    counter.push_back(
        {UnrootedSBNMaps::IndexerRepresentationOf(indexer, topology, default_index),
         topology_count});
  }
  return counter;
}

StringSetVector UnrootedSBNMaps::StringIndexerRepresentationOf(
    const StringVector& reversed_indexer,
    const UnrootedIndexerRepresentation& indexer_representation) {
  StringSetVector string_sets;
  for (const auto& rooted_representation : indexer_representation) {
    string_sets.push_back(RootedSBNMaps::StringIndexerRepresentationOf(
        reversed_indexer, rooted_representation));
  }
  return string_sets;
}

BitsetSizeDict RootedSBNMaps::RootsplitCounterOf(
    const Node::TopologyCounter& topologies) {
  BitsetSizeDict rootsplit_counter(0);
  for (const auto& [topology, topology_count] : topologies) {
    rootsplit_counter.increment(Rootsplit(topology.get()), topology_count);
  }
  return rootsplit_counter;
}

PCSPCounter RootedSBNMaps::PCSPCounterOf(const Node::TopologyCounter& topologies) {
  PCSPCounter pcsp_dict;
  for (const auto& [topology, topology_count] : topologies) {
    auto leaf_count = topology->LeafCount();
    Assert(topology->Children().size() == 2,
           "RootedSBNMaps::PCSPCounterOf was expecting a bifurcating tree!");
    topology->RootedPCSPPreOrder(
        [&pcsp_dict, &topology_count = topology_count, &leaf_count](
            const Node* sister_node, const Node* focal_node, const Node* child0_node,
            const Node* child1_node) {
          AddToPCSPCounter(pcsp_dict, topology_count, leaf_count, sister_node, false,
                           focal_node, false, child0_node, false, child1_node, false);
        });
  }
  return pcsp_dict;
}

SizeVector RootedSBNMaps::IndexerRepresentationOf(const BitsetSizeMap& indexer,
                                                  const Node::NodePtr& topology,
                                                  const size_t default_index) {
  const auto leaf_count = topology->LeafCount();
  SizeVector result;
  // Start with the rootsplit.
  Bitset rootsplit = Rootsplit(topology.get());
  result.push_back(AtWithDefault(indexer, rootsplit, default_index));
  // Now add the PCSPs.
  topology->RootedPCSPPreOrder([&leaf_count, &indexer, &default_index, &result](
                                   const Node* sister_node, const Node* focal_node,
                                   const Node* child0_node, const Node* child1_node) {
    Bitset pcsp_bitset =
        SBNMaps::PCSPBitsetOf(leaf_count, sister_node, false, focal_node, false,
                              child0_node, false, child1_node, false);
    result.push_back(AtWithDefault(indexer, pcsp_bitset, default_index));
  });
  return result;
}

StringSet RootedSBNMaps::StringIndexerRepresentationOf(
    const StringVector& reversed_indexer,
    const RootedIndexerRepresentation& indexer_representation) {
  StringSet string_set;
  for (const auto index : indexer_representation) {
    SafeInsert(string_set, reversed_indexer.at(index));
  }
  return string_set;
}

RootedIndexerRepresentationCounter RootedSBNMaps::IndexerRepresentationCounterOf(
    const BitsetSizeMap& indexer, const Node::TopologyCounter& topology_counter,
    const size_t default_index) {
  RootedIndexerRepresentationCounter counter;
  counter.reserve(topology_counter.size());
  for (const auto& [topology, topology_count] : topology_counter) {
    counter.push_back(
        {RootedSBNMaps::IndexerRepresentationOf(indexer, topology, default_index),
         topology_count});
  }
  return counter;
}

void RootedSBNMaps::IncrementRootedIndexerRepresentationSizeDict(
    RootedIndexerRepresentationSizeDict& dict,
    RootedIndexerRepresentation rooted_indexer_representation) {
  Assert(rooted_indexer_representation.size() > 1,
         "Rooted indexer representation is too small in "
         "IncrementRootedIndexerRepresentationSizeDict!");
  std::sort(rooted_indexer_representation.begin() + 1,
            rooted_indexer_representation.end());
  dict.increment(rooted_indexer_representation, 1);
}

void RootedSBNMaps::IncrementRootedIndexerRepresentationSizeDict(
    RootedIndexerRepresentationSizeDict& dict,
    const UnrootedIndexerRepresentation& indexer_representation) {
  for (const auto& rooted_indexer_representation : indexer_representation) {
    IncrementRootedIndexerRepresentationSizeDict(dict, rooted_indexer_representation);
  }
}
