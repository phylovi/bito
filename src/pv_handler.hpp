// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// PVHandler is used for storing and manipulating Partial Vectors.  Partial Vectors are
// intermediate computations, such as in likelihoods or parsimonies, used for performing
// dynamic programming on a tree or DAG. Partial Vectors can be "stored on" and indexed
// according to different elements of the DAG: either by the edges or the nodes.
//
// PSVHandler is used to perform the Sankoff algorithm. There are 3 partial vectors:
// PLeft, PRight, and Q. PLeft and Pright corresponds to Sankoff vectors for the left
// and right child respectively, and Q corresponds to the value of the partial vector
// pointing leaf-ward.

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
  static PSVType PPSVType(const bool is_on_left) {
    return is_on_left ? PSVType::PLeft : PSVType::PRight;
  }
  static PSVType PPSVType(const SubsplitClade clade) {
    return (clade == SubsplitClade::Left) ? PSVType::PLeft : PSVType::PRight;
  }

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

using PVId = GenericId<struct PVIdTag>;
using PVIdVector = std::vector<PVId>;

using PLVType = PartialVectorType::PLVType;
using PLVTypeEnum = PartialVectorType::PLVTypeEnum;
using PSVType = PartialVectorType::PSVType;
using PSVTypeEnum = PartialVectorType::PSVTypeEnum;

// PVTypeEnum determines which PV types need to be stored on each element of the
// PVHandler (e.g. P-PVs, Q-PVs, R-PVs). DAGElementId decides whether indexing PVs
// according to DAG's nodes or edges.
template <class PVTypeEnum, class DAGElementId>
class PartialVectorHandler {
 public:
  using TypeEnum = PVTypeEnum;
  using PVType = typename TypeEnum::Type;

  PartialVectorHandler(const std::string &mmap_file_path, const size_t element_count,
                       const size_t pattern_count, const double resizing_factor)
      : element_count_(element_count),
        pattern_count_(pattern_count),
        resizing_factor_(resizing_factor),
        mmap_file_path_(mmap_file_path),
        mmapped_master_pvs_(mmap_file_path_, (element_count + element_spare_count_) *
                                                 pv_count_per_element_ *
                                                 size_t(resizing_factor_) *
                                                 pattern_count) {}

  // ** Counts

  double GetByteCount() const { return mmapped_master_pvs_.ByteCount(); }
  size_t GetPVCountPer() const { return pv_count_per_element_; }
  size_t GetSitePatternCount() const { return pattern_count_; }
  // DAG element counts.
  size_t GetCount() const { return element_count_; }
  size_t GetSpareCount() const { return element_spare_count_; }
  size_t GetAllocatedCount() const { return element_alloc_; }
  size_t GetPaddedCount() const { return GetCount() + GetSpareCount(); }
  // PV counts.
  size_t GetPVCount() const { return GetCount() * GetPVCountPer(); }
  size_t GetSparePVCount() const { return GetSpareCount() * GetPVCountPer(); }
  size_t GetPaddedPVCount() const { return GetPaddedCount() * GetPVCountPer(); }
  size_t GetAllocatedPVCount() const { return GetAllocatedCount() * GetPVCountPer(); }

  void SetCount(const size_t element_count) { element_count_ = element_count; }
  void SetSpareCount(const size_t element_spare_count) {
    element_spare_count_ = element_spare_count;
  }
  void SetAllocatedCount(const size_t element_alloc) { element_alloc_ = element_alloc; }

  // ** Resize

  // Resize PVHandler to accomodate DAG with given number of nodes.
  void Resize(const size_t new_element_count, const size_t new_element_alloc);
  // Reindex PV according to pv_reindexer.
  void Reindex(const Reindexer pv_reindexer);
  // Expand element_reindexer into pv_reindexer.
  Reindexer BuildPVReindexer(const Reindexer &element_reindexer,
                             const size_t old_element_count,
                             const size_t new_element_count);

  // ** Access

  // Get vector of all Partial Vectors.
  NucleotidePLVRefVector &GetPVs() { return pvs_; }
  const NucleotidePLVRefVector &GetPVs() const { return pvs_; }
  // Get PV by absolute index from the vector of Partial Vectors.
  NucleotidePLVRef &GetPV(const PVId pv_idx) { return pvs_.at(pv_idx.value_); }
  const NucleotidePLVRef &GetPV(const PVId pv_idx) const {
    return pvs_.at(pv_idx.value_);
  }
  NucleotidePLVRef &operator()(const PVId pv_idx) { return GetPV(pv_idx); }
  const NucleotidePLVRef &operator()(const PVId pv_idx) const { return GetPV(pv_idx); }
  // Get PV by PV type and node index from the vector of Partial Vectors.
  NucleotidePLVRef &GetPV(const PVType pv_type, const DAGElementId element_idx) {
    return GetPV(GetPVIndex(pv_type, element_idx));
  }
  const NucleotidePLVRef &GetPV(const PVType pv_type,
                                const DAGElementId element_idx) const {
    return GetPV(GetPVIndex(pv_type, element_idx));
  }
  NucleotidePLVRef &operator()(const PVType pv_type, const DAGElementId element_idx) {
    return GetPV(GetPVIndex(pv_type, element_idx));
  }
  const NucleotidePLVRef &operator()(const PVType pv_type,
                                     const DAGElementId element_idx) const {
    return GetPV(GetPVIndex(pv_type, element_idx));
  }
  // Get Spare PV by index from the vector of Partial Vectors.
  NucleotidePLVRef &GetSparePV(const PVId pv_idx) {
    return GetPV(GetSparePVIndex(pv_idx));
  }
  const NucleotidePLVRef &GetSparePV(const PVId pv_idx) const {
    return GetPV(GetSparePVIndex(pv_idx));
  }

  // Get total offset into PVs, indexed based on underlying DAG.
  static PVId GetPVIndex(const PVType pv_type, const DAGElementId element_idx,
                         const size_t element_count) {
    return GetPVIndex(TypeEnum::GetIndex(pv_type), element_idx, element_count);
  }
  PVId GetPVIndex(const PVType pv_type, const DAGElementId element_idx) const {
    Assert(element_idx.value_ < GetCount(), "Requested element_idx is out-of-range.");
    return GetPVIndex(pv_type, element_idx, GetCount());
  }
  // Get total offset into temporary PVs, indexed based on underlying grafted DAG.
  PVId GetSparePVIndex(const PVId pv_idx) const {
    const size_t pv_scratch_size = GetPaddedPVCount() - GetPVCount();
    Assert(pv_idx < pv_scratch_size,
           "Requested temporary pv_idx outside of allocated scratch space.");
    return PVId(pv_idx.value_ + GetPVCount());
  }
  // Get vector of all node ids for given node.
  static PVIdVector GetPVIndexVectorForElementId(const DAGElementId element_idx,
                                                 const size_t element_count) {
    PVIdVector pv_idxs;
    for (const auto pv_type : typename TypeEnum::Iterator()) {
      pv_idxs.push_back(GetPVIndex(pv_type, element_idx, element_count));
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
  std::string ToString(const PVType pv_type, const DAGElementId element_idx,
                       const bool show_labels = false) const {
    std::stringstream out;
    out << "PV[" << PVTypeEnum::ToString(pv_type) << ", Element" << element_idx
        << ", PVId::" << GetPVIndex(pv_type, element_idx) << "]: " << std::endl;
    for (auto &&row : GetPV(pv_type, element_idx).rowwise()) {
      out << row << std::endl;
    }
    return out.str();
  }
  std::string AllPVsToString(const bool show_labels = false) const {
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
  // Get index for given PV enum.
  static size_t GetPVTypeIndex(const PVType pv_type) {
    return TypeEnum::GetIndex(pv_type);
  }

  // ** Data Sizing
  // "Count" is the currently occupied by data.
  // "Padding" is the amount of free working space added to end of occupied space.
  // "Alloc" is the total current memory allocation.
  // "Resizing factor" is the amount of extra storage allocated for when resizing.

  // Number of nodes in DAG.
  size_t element_count_ = 0;
  // Number of nodes of additional space for temporary graft nodes in DAG.
  size_t element_spare_count_ = 2;
  // Number of nodes allocated for in PVHandler.
  size_t element_alloc_ = 0;
  // Size of Site Pattern.
  size_t pattern_count_ = 0;
  // Number of PVs for each node in DAG.
  const size_t pv_count_per_element_ = PVTypeEnum::Count;
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
template <class DAGElementId>
class PLVHandler
    : public PartialVectorHandler<PartialVectorType::PLVTypeEnum, DAGElementId> {
 public:
  using PLVType = PartialVectorType::PLVType;
  using PLVTypeEnum = PartialVectorType::PLVTypeEnum;
  using PLVTypeIterator = PLVTypeEnum::Iterator;
  static const inline size_t plv_count_ = PLVTypeEnum::Count;

  PLVHandler(const std::string &mmap_file_path, const size_t element_count,
             const size_t pattern_count, const double resizing_factor)
      : PartialVectorHandler<PLVTypeEnum, DAGElementId>(
            mmap_file_path, element_count, pattern_count, resizing_factor) {}
};

// PSVHandler: Partial Sankoff Vector Handler
template <class DAGElementId>
class PSVHandler
    : public PartialVectorHandler<PartialVectorType::PSVTypeEnum, DAGElementId> {
 public:
  using PSVType = PartialVectorType::PSVType;
  using PSVTypeEnum = PartialVectorType::PSVTypeEnum;
  using PSVTypeIterator = PSVTypeEnum::Iterator;
  static const inline size_t psv_count_ = PartialVectorType::PSVTypeEnum::Count;

  PSVHandler(const std::string &mmap_file_path, const size_t element_count,
             const size_t pattern_count, const double resizing_factor)
      : PartialVectorHandler<PSVTypeEnum, DAGElementId>(
            mmap_file_path, element_count, pattern_count, resizing_factor) {}
};

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
