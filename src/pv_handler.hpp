// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// PLVHandler is used for storing and manipulating Partial Vectors.  Partial Vector are
// intermediate computations, such as likelihoods or other cost matrices, used for
// performing dynamic programming on a tree or DAG.

#pragma once

#include "sugar.hpp"
#include "mmapped_plv.hpp"
#include "site_pattern.hpp"
#include "reindexer.hpp"

// Enumerated Types for Partial Vectors.
namespace PartialVectorType {
// PLV: Partial Likelihood Vectors
enum class PLVType : size_t {
  P,          // p(s)
  PHatRight,  // phat(s_right)
  PHatLeft,   // phat(s_left)
  RHat,       // rhat(s_right) = rhat(s_left)
  RRight,     // r(s_right)
  RLeft,      // r(s_left)
};
static inline const size_t PLVCount = 6;
class PLVTypeEnum
    : public EnumWrapper<PLVType, size_t, PLVCount, PLVType::P, PLVType::RLeft> {
 public:
  static inline const Array<std::string> Labels = {
      {"P", "PHatRight", "PHatLeft", "RHat", "RRight", "RLeft"}};

  friend std::ostream &operator<<(std::ostream &os, const Type e) {
    os << "PLV::" << Labels[e];
    return os;
  };
};

// PSV: Partial Sankoff Vectors
enum class PSVType : size_t {
  PRight,  // p(s_right)
  PLeft,   // p(s_left)
  Q        // q(s)
};
static inline const size_t PSVCount = 3;
class PSVTypeEnum
    : public EnumWrapper<PSVType, size_t, PSVCount, PSVType::PRight, PSVType::Q> {
 public:
  static inline const Array<std::string> Labels = {{"PRight", "PLeft", "Q"}};

  friend std::ostream &operator<<(std::ostream &os, const Type e) {
    os << "PSV::" << Labels[e];
    return os;
  };
};
};  // namespace PartialVectorType

using PLVType = PartialVectorType::PLVType;
using PLVTypeEnum = PartialVectorType::PLVTypeEnum;
using PSVType = PartialVectorType::PSVType;
using PSVTypeEnum = PartialVectorType::PSVTypeEnum;
using PVId = GenericId<struct PVIdTag>;

template <typename PVType, typename PVTypeEnum>
class PartialVectorHandler {
 public:
  using Type = PVType;
  using TypeEnum = PVTypeEnum;

  PartialVectorHandler(const std::string &mmap_file_path, const size_t node_count,
                       const size_t pattern_count, const double resizing_factor)
      : node_count_(node_count),
        pattern_count_(pattern_count),
        resizing_factor_(resizing_factor),
        mmap_file_path_(mmap_file_path),
        mmapped_master_pvs_(mmap_file_path_,
                            (node_count + node_spare_count_) * pv_count_per_node_ *
                                size_t(resizing_factor_) * pattern_count) {}

  // ** Counts

  double GetByteCount() const { return mmapped_master_pvs_.ByteCount(); }
  size_t GetPVCountPerNode() const { return pv_count_per_node_; }
  size_t GetSitePatternCount() const { return pattern_count_; }
  size_t GetNodeCount() const { return node_count_; }
  size_t GetSpareNodeCount() const { return node_spare_count_; }
  size_t GetAllocatedNodeCount() const { return node_alloc_; }
  size_t GetPaddedNodeCount() const { return GetNodeCount() + GetSpareNodeCount(); }
  size_t GetPVCount() const { return GetNodeCount() * GetPVCountPerNode(); }
  size_t GetSparePVCount() const { return GetSpareNodeCount() * GetPVCountPerNode(); }
  size_t GetPaddedPVCount() const { return GetPaddedNodeCount() * GetPVCountPerNode(); }
  size_t GetAllocatedPVCount() const {
    return GetAllocatedNodeCount() * GetPVCountPerNode();
  }

  void SetNodeCount(const size_t node_count) { node_count_ = node_count; }
  void SetSpareNodeCount(const size_t node_spare_count) {
    node_spare_count_ = node_spare_count;
  }
  void SetAllocatedNodeCount(const size_t node_alloc) { node_alloc_ = node_alloc; }

  // ** Resize

  // Resize PVHandler to accomodate DAG with given number of nodes.
  void Resize(const size_t new_node_count, const size_t new_node_alloc);
  // Reindex PV according to pv_reindexer.
  void Reindex(const Reindexer pv_reindexer);
  // Expand node_reindexer into pv_reindexer.
  Reindexer BuildPVReindexer(const Reindexer &node_reindexer,
                             const size_t old_node_count, const size_t new_node_count);

  // ** Access

  // Get vector of all Partial Vectors.
  NucleotidePLVRefVector &GetPVs() { return pvs_; }
  const NucleotidePLVRefVector &GetPVs() const { return pvs_; }
  // Get PV by absolute index from the vector of Partial Vectors.
  NucleotidePLVRef &GetPV(const size_t pv_idx) { return pvs_.at(pv_idx); };
  const NucleotidePLVRef &GetPV(const size_t pv_idx) const { return pvs_.at(pv_idx); };
  NucleotidePLVRef &operator()(const size_t pv_idx) { return GetPV(pv_idx); };
  const NucleotidePLVRef &operator()(const size_t pv_idx) const {
    return GetPV(pv_idx);
  };
  // Get PV by PV type and node index from the vector of Partial Vectors.
  NucleotidePLVRef &GetPV(const PVType pv_type, const size_t node_idx) {
    return pvs_.at(GetPVIndex(pv_type, node_idx));
  };
  const NucleotidePLVRef &GetPV(const PVType pv_type, const size_t node_idx) const {
    return pvs_.at(GetPVIndex(pv_type, node_idx));
  };
  NucleotidePLVRef &operator()(const PVType pv_type, const size_t node_idx) {
    return GetPV(GetPVIndex(pv_type, node_idx));
  };
  const NucleotidePLVRef &operator()(const PVType pv_type,
                                     const size_t node_idx) const {
    return GetPV(GetPVIndex(pv_type, node_idx));
  };
  // Get Spare PV by index from the vector of Partial Vectors.
  NucleotidePLVRef &GetSparePV(const size_t pv_idx) {
    return pvs_.at(GetSparePVIndex(pv_idx));
  };
  const NucleotidePLVRef &GetSparePV(const size_t pv_idx) const {
    return pvs_.at(GetSparePVIndex(pv_idx));
  };

  // Get total offset into PVs, indexed based on underlying DAG.
  static size_t GetPVIndex(const PVType pv_type, const size_t node_idx,
                           const size_t node_count) {
    // return GetPVIndex(PVTypeEnum::GetIndex(pv_type), node_idx, node_count);
    return GetPVIndex(GetPVTypeIndex(pv_type), node_idx, node_count);
  };
  size_t GetPVIndex(const PVType pv_type, const size_t node_idx) const {
    return GetPVIndex(pv_type, node_idx, GetNodeCount());
  }

  // Get total offset into temporary PVs, indexed based on underlying grafted DAG.
  size_t GetSparePVIndex(const size_t pv_idx) const {
    const size_t pv_scratch_size = GetPaddedPVCount() - GetPVCount();
    Assert(pv_idx < pv_scratch_size,
           "Requested temporary pv_idx outside of allocated scratch space.");
    return pv_idx + GetPVCount();
  }

  // Get vector of all node ids for given node.
  static SizeVector GetPVIndexVectorForNodeId(const size_t node_idx,
                                              const size_t node_count) {
    SizeVector pv_idxs;
    for (const auto pv_type : typename TypeEnum::Iterator()) {
      pv_idxs.push_back(GetPVIndex(pv_type, node_idx, node_count));
    }
    return pv_idxs;
  };

 protected:
  // Get total offset into PVs.
  static size_t GetPVIndex(const size_t pv_type_idx, const size_t node_idx,
                           const size_t node_count) {
    return (pv_type_idx * node_count) + node_idx;
  };
  // Get index for given PV enum.
  static size_t GetPVTypeIndex(const PVType pv_type) {
    return TypeEnum::GetIndex(pv_type);
  };

  // ** Data Sizing
  // "Count" is the currently occupied by data.
  // "Padding" is the amount of free working space added to end of occupied space.
  // "Alloc" is the total current memory allocation.
  // "Resizing factor" is the amount of extra storage allocated for when resizing.

  // Number of nodes in DAG.
  size_t node_count_ = 0;
  // Number of nodes of additional space for temporary graft nodes in DAG.
  size_t node_spare_count_ = 2;
  // Number of nodes allocated for in PVHandler.
  size_t node_alloc_ = 0;
  // Size of Site Pattern.
  size_t pattern_count_ = 0;
  // Number of PVs for each node in DAG.
  const size_t pv_count_per_node_ = PVTypeEnum::Count;
  // When size exceeds current allocation, ratio to grow new allocation.
  double resizing_factor_ = 2.0;

  // File path to data map.
  std::string mmap_file_path_;
  // Master PV: Large data block of virtual memory for Partial Likelihood Vectors.
  // Subdivided into sections for pvs_.
  MmappedNucleotidePLV mmapped_master_pvs_;
  // Partial Vectors.
  // Divides mmapped_master_pvs_.
  // For example, GP PLVs are divided as follows:
  // - [0, num_nodes): p(s).
  // - [num_nodes, 2*num_nodes): phat(s_right).
  // - [2*num_nodes, 3*num_nodes): phat(s_left).
  // - [3*num_nodes, 4*num_nodes): rhat(s_right) = rhat(s_left).
  // - [4*num_nodes, 5*num_nodes): r(s_right).
  // - [5*num_nodes, 6*num_nodes): r(s_left).
  NucleotidePLVRefVector pvs_;
};

// PLVHandler: Partial Likelihood Vector Handler
class PLVHandler : public PartialVectorHandler<PartialVectorType::PLVType,
                                               PartialVectorType::PLVTypeEnum> {
 public:
  using PLVTypeEnum = PartialVectorType::PLVTypeEnum;
  using PLVType = PLVTypeEnum::Type;
  using PLVTypeIterator = PLVTypeEnum::Iterator;
  static const inline size_t plv_count_ = PLVTypeEnum::Count;

  PLVHandler(const std::string &mmap_file_path, const size_t node_count,
             const size_t pattern_count, const double resizing_factor)
      : PartialVectorHandler(mmap_file_path, node_count, pattern_count,
                             resizing_factor){};

  static Type RPLVType(const bool is_on_left) {
    return is_on_left ? PLVType::RLeft : PLVType::RRight;
  };

  static Type PPLVType(const bool is_on_left) {
    return is_on_left ? PLVType::PHatLeft : PLVType::PHatRight;
  };
};

// PSVHandler: Partial Sankoff Vector Handler
class PSVHandler : public PartialVectorHandler<PartialVectorType::PSVType,
                                               PartialVectorType::PSVTypeEnum> {
 public:
  using PSVTypeEnum = PartialVectorType::PSVTypeEnum;
  using PSVType = PSVTypeEnum::Type;
  using PSVTypeIterator = PSVTypeEnum::Iterator;
  static const inline size_t psv_count_ = PartialVectorType::PSVTypeEnum::Count;

  PSVHandler(const std::string &mmap_file_path, const size_t node_count,
             const size_t pattern_count, const double resizing_factor)
      : PartialVectorHandler(mmap_file_path, node_count, pattern_count,
                             resizing_factor){};

  static Type PPLVType(const bool is_on_left) {
    return is_on_left ? PSVType::PLeft : PSVType::PRight;
  }
};

#ifdef DOCTEST_LIBRARY_INCLUDED

// Check that PLV iterator iterates over all PLVs exactly once.
TEST_CASE("PLVHandler: EnumIterator") {
  using namespace PartialVectorType;
  const auto plv_types = PLVTypeEnum::TypeArray();
  std::map<PLVType, size_t> plv_visited_map;
  // Iterate using vector.
  for (const PLVType plv_type : plv_types) {
    plv_visited_map.insert({plv_type, 0});
  }
  // Iterate using EnumIterator.
  for (const PLVType plv_type : PLVTypeEnum::Iterator()) {
    CHECK_MESSAGE(plv_visited_map.find(plv_type) != plv_visited_map.end(),
                  "Iterator has PLV not in plv_vector.");
    plv_visited_map.at(plv_type) += 1;
  }
  // Check that each was visited only once.
  for (const auto [plv_type, visit_count] : plv_visited_map) {
    std::ignore = plv_type;
    CHECK_FALSE_MESSAGE(visit_count < 1, "One or more PLVs skipped by EnumIterator.");
    CHECK_FALSE_MESSAGE(visit_count > 1,
                        "One or more PLVs in visited more than once by EnumIterator.");
  }
}

#endif  // DOCTEST_LIBRARY_INCLUDED
