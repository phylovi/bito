// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// In our formulation of variational Bayes phylogenetics, we do inference of continuous
// parameters with respect to a "subsplit support", which is the collection of
// rootsplits and PCSPs that are allowed to be nonzero.
//
// We implement this concept as an SBNSupport here. We store the support as a collection
// of indexing maps that map from Bitset representations of rootsplits and PCSPs to some
// indexing scheme. Details differ if we are looking at rooted or unrooted trees, hence
// this class gets subclassed by RootedSBNSupport and UnrootedSBNSupport.
//
// To learn more about these indexer maps, see the unit tests in rooted_sbn_instance.hpp
// and unrooted_sbn_instance.hpp.

#pragma once

#include "psp_indexer.hpp"
#include "sbn_probability.hpp"
#include "subsplit_dag.hpp"

class SBNSupport {
 public:
  explicit SBNSupport(SubsplitDAG* dag) : dag_{dag} {}
  inline size_t GPCSPCount() const { return dag_->EdgeCountWithLeafSubsplits(); }
  inline bool Empty() const { return GPCSPCount() == 0; }
  inline size_t TaxonCount() const { return dag_->TaxonCount(); }
  inline StringVector TaxonNames() const { return dag_->GetSortedVectorOfTaxonNames(); }
  inline size_t RootsplitCount() const { return dag_->RootsplitCount(); }
  const Bitset &RootsplitsAt(size_t rootsplit_idx) const {
    auto iter = dag_->GetRootsplitNodeIds().begin();
    std::advance(iter, rootsplit_idx);
    return dag_->GetDAGNode(*iter).GetBitset();
  }
  const size_t &IndexerAt(const Bitset &bitset) const { return dag_->GetSubsplitToIdMap().at(bitset); }
  inline bool ParentInSupport(const Bitset &parent) const {
    return dag_->GetParentToChildRange().count(parent) > 0;
  }
  inline const SizePair &ParentToRangeAt(const Bitset &parent) const {
    return dag_->GetParentToChildRange().at(parent);
  }
  inline const Bitset &IndexToChildAt(size_t child_idx) const {
    auto edge = dag_->GetDAGEdge(child_idx);
    return dag_->GetDAGNode(edge.GetChild()).GetBitset();
  }

  const BitsetSizePairMap &ParentToRange() const { return dag_->GetParentToChildRange(); }

  BitsetSizeMap Indexer() const {
    BitsetSizeMap result;
    for (size_t edge_id = 0; edge_id < dag_->EdgeCount(); ++edge_id) {
      auto edge = dag_->GetDAGEdge(edge_id);
      auto parent = dag_->GetDAGNode(edge.GetParent());
      auto child = dag_->GetDAGNode(edge.GetChild());
      auto pcsp = Bitset::PCSP(parent.GetBitset(), child.GetBitset());
      result[pcsp] = edge_id;
    }
    return result;
  }

  PSPIndexer BuildPSPIndexer() const {
    BitsetVector rootsplits;
    for (auto id : dag_->GetRootsplitNodeIds()) {
      rootsplits.push_back(dag_->GetDAGNode(id).GetBitset());
    }
    return PSPIndexer(rootsplits, Indexer());
  }

  // "Pretty" string representation of the indexer.
  StringVector PrettyIndexer() const;
  void PrettyPrintIndexer() const;

  // Return indexer_ and parent_to_child_range_ converted into string-keyed maps.
  std::tuple<StringSizeMap, StringSizePairMap> GetIndexers() const;

  // Get the indexer, but reversed and with bitsets appropriately converted to
  // strings.
  StringVector StringReversedIndexer() const;

  void ProbabilityNormalizeSBNParametersInLog(EigenVectorXdRef sbn_parameters) const;

 protected:
  SubsplitDAG* dag_;
};
