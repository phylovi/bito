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
#include "subsplit_dag_storage.hpp"

// Helper Enumerated Types for Partial Vectors.
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
  static PLVType PPLVType(const bool is_on_left) {
    return is_on_left ? PLVType::PHatLeft : PLVType::PHatRight;
  };
  static PLVType RPLVType(const bool is_on_left) {
    return is_on_left ? PLVType::RLeft : PLVType::RRight;
  };

  static inline const Array<std::string> Labels = {
      {"P", "PHatRight", "PHatLeft", "RHat", "RRight", "RLeft"}};

  static std::string ToString(const PLVType e) {
    std::stringstream ss;
    ss << "PLV::" << Labels[e];
    return ss.str();
  }
  friend std::ostream &operator<<(std::ostream &os, const PLVType e) {
    os << ToString(e);
    return os;
  }
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

  static std::string ToString(const PSVType e) {
    std::stringstream ss;
    ss << "PSV::" << Labels[e];
    return ss.str();
  }
  friend std::ostream &operator<<(std::ostream &os, const PSVType e) {
    os << ToString(e);
    return os;
  }
};
};  // namespace PartialVectorType

// using PVId = GenericId<struct PVIdTag>;
using PVId = size_t;
using PLVType = PartialVectorType::PLVType;
using PLVTypeEnum = PartialVectorType::PLVTypeEnum;
using PSVType = PartialVectorType::PSVType;
using PSVTypeEnum = PartialVectorType::PSVTypeEnum;

template <typename PVType, typename PVTypeEnum, typename DAGElementId>
class PartialVectorHandler {
 public:
  using Type = PVType;
  using TypeEnum = PVTypeEnum;

  PartialVectorHandler(const std::string &mmap_file_path, const DAGElementId pv_count,
                       const size_t pattern_count, const double resizing_factor)
      : pv_count_(pv_count),
        pattern_count_(pattern_count),
        resizing_factor_(resizing_factor),
        mmap_file_path_(mmap_file_path),
        mmapped_master_pvs_(mmap_file_path_,
                            (pv_count_ + pv_spare_count_) * PVTypeEnum::Count *
                                size_t(resizing_factor_) * pattern_count) {}

  // ** Counts

  double GetByteCount() const { return mmapped_master_pvs_.ByteCount(); }
  size_t GetSitePatternCount() const { return pattern_count_; }
  size_t GetPVCount() const { return pv_count_; }
  size_t GetSparePVCount() const { return pv_spare_count_; }
  size_t GetPaddedPVCount() const { return GetPVCount() * GetSparePVCount(); }
  size_t GetAllocatedPVCount() const { return pv_alloc_count_; }

  void SetPVCount(const size_t pv_count) { pv_count_ = pv_count; }
  void SetSparePVCount(const size_t spare_pv_count) {
    pv_spare_count_ = spare_pv_count;
  }
  void SetAllocatedPVCount(const size_t alloc_pv_count) {
    pv_alloc_count_ = alloc_pv_count;
  }

  // ** Resize

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
  // NucleotidePLVRef &GetPV(const PVId pv_idx) { return pvs_.at(pv_idx.value_); }
  // const NucleotidePLVRef &GetPV(const PVId pv_idx) const {
  //   return pvs_.at(pv_idx.value_);
  // }
  NucleotidePLVRef &GetPV(const size_t pv_idx) { return pvs_.at(pv_idx); }
  const NucleotidePLVRef &GetPV(const size_t pv_idx) const { return pvs_.at(pv_idx); }
  NucleotidePLVRef &operator()(const PVId pv_idx) { return GetPV(pv_idx); }
  const NucleotidePLVRef &operator()(const PVId pv_idx) const { return GetPV(pv_idx); }
  // Get Spare PV by index from the vector of Partial Vectors.
  NucleotidePLVRef &GetSparePV(const PVId pv_idx) {
    return pvs_.at(GetSparePVIndex(pv_idx));
  }
  const NucleotidePLVRef &GetSparePV(const PVId pv_idx) const {
    return pvs_.at(GetSparePVIndex(pv_idx));
  }

  // Get total offset into temporary PVs, indexed based on underlying grafted DAG.
  PVId GetSparePVIndex(const PVId pv_idx) const {
    const size_t pv_scratch_size = GetPaddedPVCount() - GetPVCount();
    Assert(pv_idx < pv_scratch_size,
           "Requested temporary pv_idx outside of allocated scratch space.");
    return pv_idx + PVId(GetPVCount());
  }

  // Get vector of all node ids for given node.
  static SizeVector GetPVIndexVectorForDAGElementId(const DAGElementId node_idx,
                                                    const size_t node_count) {
    SizeVector pv_idxs;
    for (const auto pv_type : typename TypeEnum::Iterator()) {
      pv_idxs.push_back(GetPVIndex(pv_type, node_idx, node_count));
    }
    return pv_idxs;
  }

  // ** I/O

  // Output data to string.
  std::string ToString(const PVId pv_idx, const bool show_labels = false) const {
    std::stringstream out;
    out << "PV[" << pv_idx << "]: " << std::endl;
    for (auto &&row : GetPV(pv_idx).rowwise()) {
      out << row << std::endl;
    }
    return out.str();
  }

 protected:
  // Get index for given PV enum.
  static size_t GetPVTypeIndex(const PVType pv_type) {
    return TypeEnum::GetIndex(pv_type);
  }

  // ** Data Sizing
  // "Count" is the currently occupied by data.
  // "Padding" is the amount of free working space added to end of occupied space.
  // "Alloc" is the total current memory allocation.
  // "Resizing factor" is the amount of extra storage allocated for when resizing.

  // Number of active PVs.
  size_t pv_count_ = 0;
  // Number of nodes of additional space for temporary PVs.
  size_t pv_spare_count_ = 2;
  // Number of PVs allocated for in PVHandler.
  size_t pv_alloc_count_ = 0;

  // Size of Site Pattern.
  size_t pattern_count_ = 0;
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

  // ** DAG Element
 public:
  // ** Counts

  size_t GetPVCountPer() const { return pv_count_per_element_; }
  size_t GetCount() const { return element_count_; }
  size_t GetSpareCount() const { return element_spare_count_; }
  size_t GetAllocatedCount() const { return element_alloc_count_; }
  size_t GetPaddedCount() const { return GetCount() + GetSpareCount(); }

  void SetCount(const size_t element_count) {
    element_count_ = element_count;
    pv_count_ = element_count_ * pv_count_per_element_;
  }
  void SetSpareCount(const size_t element_spare_count) {
    element_spare_count_ = element_spare_count;
    pv_spare_count_ = element_spare_count_ * pv_count_per_element_;
  }
  void SetAllocatedCount(const size_t element_alloc_count) {
    element_alloc_count_ = element_alloc_count;
    pv_alloc_count_ = element_alloc_count_ * pv_count_per_element_;
  }

  // ** Access

  // Get PV by PV type and element index from the vector of Partial Vectors.
  NucleotidePLVRef &operator()(const PVType pv_type, const DAGElementId element_idx) {
    return GetPV(GetPVIndex(pv_type, element_idx));
  }
  const NucleotidePLVRef &operator()(const PVType pv_type,
                                     const DAGElementId element_idx) const {
    return GetPV(GetPVIndex(pv_type, element_idx));
  }
  NucleotidePLVRef &GetPV(const PVType pv_type, const DAGElementId element_idx) {
    return pvs_.at(GetPVIndex(pv_type, element_idx));
  }
  const NucleotidePLVRef &GetPV(const PVType pv_type,
                                const DAGElementId element_idx) const {
    return pvs_.at(GetPVIndex(pv_type, element_idx));
  }

  // Get total offset into PVs, indexed based on underlying DAG.
  static PVId GetPVIndex(const PVType pv_type, const DAGElementId element_idx,
                         const size_t element_count) {
    return GetPVIndex(GetPVTypeIndex(pv_type), element_idx, element_count);
  }
  PVId GetPVIndex(const PVType pv_type, const DAGElementId element_idx) const {
    Assert(element_idx.value_ < GetCount(), "Requested element_idx is out-of-range.");
    return GetPVIndex(pv_type, element_idx, GetCount());
  }

  // ** Resize

  // Resize PVHandler to accomodate DAG with given number of elements.
  void Resize(const size_t new_element_count, const size_t new_element_alloc);

  // ** I/O

  std::string ToString(const PVType pv_type, const DAGElementId element_idx,
                       const bool show_labels = false) const {
    std::stringstream out;
    out << "PV[" << PVTypeEnum::ToString(pv_type) << ", DAGElementId" << element_idx
        << ", PVId::" << GetPVIndex(pv_type, element_idx) << "]: " << std::endl;
    for (auto &&row : GetPV(pv_type, element_idx).rowwise()) {
      out << row << std::endl;
    }
    return out.str();
  }
  std::string AllPVsToString(const bool show_labels = false) {
    std::stringstream out;
    for (const auto pv_type : typename PVTypeEnum::Iterator()) {
      for (DAGElementId element_id = 0; element_id < GetCount(); element_id++) {
        out << ToString(pv_type, element_id, show_labels);
      }
    }
    return out.str();
  }

 protected:
  // Get total offset into PVs.
  static PVId GetPVIndex(const size_t pv_type_idx, const DAGElementId element_idx,
                         const size_t element_count) {
    return (pv_type_idx * element_count) + element_idx.value_;
  }
  PVId GetPVIndex(const size_t pv_type_idx, const DAGElementId element_idx) {
    return (pv_type_idx * element_count_) + element_idx.value_;
  }

  // ** Data Sizing
  // "Count" is the currently occupied by data.
  // "Padding" is the amount of free working space added to end of occupied space.
  // "Alloc" is the total current memory allocation.
  // "Resizing factor" is the amount of extra storage allocated for when resizing.

  // Number of PVs for each element in DAG.
  const size_t pv_count_per_element_ = PVTypeEnum::Count;
  // Number of elements in DAG.
  size_t element_count_ = 0;
  // Number of elements of additional space for temporary graft elements in DAG.
  size_t element_spare_count_ = 2;
  // Number of elements allocated for in PVHandler.
  size_t element_alloc_count_ = 0;
};

// PLVHandler: Partial Likelihood Vector Handler
template <typename DAGElementId>
class PLVHandler : public PartialVectorHandler<PLVType, PLVTypeEnum, DAGElementId> {
 public:
  using PLVTypeEnum = PartialVectorType::PLVTypeEnum;
  using PLVType = PLVTypeEnum::Type;
  using PLVTypeIterator = PLVTypeEnum::Iterator;
  static const inline size_t plv_count_ = PLVTypeEnum::Count;

  PLVHandler(const std::string &mmap_file_path, const DAGElementId pv_count,
             const size_t pattern_count, const double resizing_factor)
      : PartialVectorHandler<PLVType, PLVTypeEnum, DAGElementId>(
            mmap_file_path, pv_count, pattern_count, resizing_factor) {}

  static PLVType PPLVType(const bool is_on_left) {
    return is_on_left ? PLVType::PHatLeft : PLVType::PHatRight;
  }
  static PLVType RPLVType(const bool is_on_left) {
    return is_on_left ? PLVType::RLeft : PLVType::RRight;
  }
  static PLVType PPLVType(const SubsplitClade clade) {
    return (clade == SubsplitClade::Left) ? PLVType::PHatLeft : PLVType::PHatRight;
  }
  static PLVType RPLVType(const SubsplitClade clade) {
    return (clade == SubsplitClade::Left) ? PLVType::RLeft : PLVType::RRight;
  }
};

// PSVHandler: Partial Sankoff Vector Handler
template <typename DAGElementId>
class PSVHandler : public PartialVectorHandler<PSVType, PSVTypeEnum, DAGElementId> {
 public:
  using PSVTypeEnum = PartialVectorType::PSVTypeEnum;
  using PSVType = PSVTypeEnum::Type;
  using PSVTypeIterator = PSVTypeEnum::Iterator;
  static const inline size_t psv_count_ = PartialVectorType::PSVTypeEnum::Count;

  PSVHandler(const std::string &mmap_file_path, const DAGElementId pv_count,
             const size_t pattern_count, const double resizing_factor)
      : PartialVectorHandler<PSVType, PSVTypeEnum, DAGElementId>(
            mmap_file_path, pv_count, pattern_count, resizing_factor) {}

  static PSVType PPLVType(const bool is_on_left) {
    return is_on_left ? PSVType::PLeft : PSVType::PRight;
  }
};

// ** Explicit Instantiation
using PLVNodeHandler = PLVHandler<NodeId>;
using PLVEdgeHandler = PLVHandler<EdgeId>;
using PSVNodeHandler = PSVHandler<NodeId>;
using PSVEdgeHandler = PSVHandler<EdgeId>;

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
