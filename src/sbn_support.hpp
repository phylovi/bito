// Copyright 2019-2021 bito project contributors.
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

#ifndef SRC_SBN_SUPPORT_HPP_
#define SRC_SBN_SUPPORT_HPP_

#include "psp_indexer.hpp"
#include "sbn_probability.hpp"

class SBNSupport {
 public:
  explicit SBNSupport(StringVector taxon_names)
      : taxon_names_(std::move(taxon_names)){};
  inline size_t GPCSPCount() const { return gpcsp_count_; }
  inline bool Empty() const { return GPCSPCount() == 0; }
  inline size_t TaxonCount() const { return taxon_names_.size(); }
  inline const StringVector &TaxonNames() const { return taxon_names_; }
  inline size_t RootsplitCount() const { return rootsplits_.size(); }
  const Bitset &RootsplitsAt(size_t rootsplit_idx) const {
    return rootsplits_.at(rootsplit_idx);
  }
  const size_t &IndexerAt(const Bitset &bitset) const { return indexer_.at(bitset); }
  inline bool ParentInSupport(const Bitset &parent) const {
    return parent_to_range_.count(parent) > 0;
  }
  inline const SizePair &ParentToRangeAt(const Bitset &parent) const {
    return parent_to_range_.at(parent);
  }
  inline const Bitset &IndexToChildAt(size_t child_idx) const {
    return index_to_child_.at(child_idx);
  }

  const BitsetSizePairMap &ParentToRange() const { return parent_to_range_; }
  const BitsetSizeMap &Indexer() const { return indexer_; }

  PSPIndexer BuildPSPIndexer() const { return PSPIndexer(rootsplits_, indexer_); }

  // "Pretty" string representation of the indexer.
  StringVector PrettyIndexer() const;
  void PrettyPrintIndexer() const;

  // Return indexer_ and parent_to_range_ converted into string-keyed maps.
  std::tuple<StringSizeMap, StringSizePairMap> GetIndexers() const;

  // Get the indexer, but reversed and with bitsets appropriately converted to
  // strings.
  StringVector StringReversedIndexer() const;

  void ProbabilityNormalizeSBNParametersInLog(EigenVectorXdRef sbn_parameters) const;

 protected:
  // A vector of the taxon names.
  StringVector taxon_names_;
  // The total number of rootsplits and PCSPs.
  size_t gpcsp_count_ = 0;
  // The master indexer for SBN parameters.
  BitsetSizeMap indexer_;
  // The collection of rootsplits, with the same indexing as in the indexer_.
  BitsetVector rootsplits_;
  // A map going from the index of a PCSP to its child.
  SizeBitsetMap index_to_child_;
  // A map going from a parent subsplit to the range of indices in
  // sbn_parameters_ with its children. See the definition of Range for the indexing
  // convention.
  BitsetSizePairMap parent_to_range_;
};

#endif  // SRC_SBN_SUPPORT_HPP_
