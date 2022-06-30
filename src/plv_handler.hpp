// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// PLVHandler is used by GPOperations to get index location of PLVs stored on the
// GPEngine.

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
  static inline const Array<std::string> Labels = {{"PLV::P", "PLV::PHatRight",
                                                    "PLV::PHatLeft", "PLV::RHat",
                                                    "PLV::RRight", "PLV::RLeft"}};

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
  static inline const Array<std::string> Labels = {
      {"PSV::PRight", "PSV::PLeft", "PSV::Q"}};

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

template <class PVType, class PVTypeEnum>
class PVHandler {
 public:
  using Type = PVType;
  using TypeEnum = PVTypeEnum;

  PVHandler(const std::string &mmap_file_path, const size_t node_count,
            const size_t pattern_count, const size_t plv_count_per_node,
            const double resizing_factor)
      : node_count_(node_count),
        pattern_count_(pattern_count),
        resizing_factor_(resizing_factor),
        mmap_file_path_(mmap_file_path),
        mmapped_master_plvs_(mmap_file_path_,
                             (node_count + node_padding_) * plv_count_per_node_ *
                                 size_t(resizing_factor_) * pattern_count) {}

  // ** Counts

  double GetByteCount() const { return mmapped_master_plvs_.ByteCount(); }
  size_t GetPLVCountPerNode() const { return plv_count_per_node_; }
  size_t GetSitePatternCount() const { return pattern_count_; }
  size_t GetNodeCount() const { return node_count_; }
  size_t GetTempNodeCount() const { return node_padding_; }
  size_t GetAllocatedNodeCount() const { return node_alloc_; }
  size_t GetPaddedNodeCount() const { return GetNodeCount() + GetTempNodeCount(); }
  size_t GetPLVCount() const { return GetNodeCount() * GetPLVCountPerNode(); }
  size_t GetTempPLVCount() const { return GetTempNodeCount() * GetPLVCountPerNode(); }
  size_t GetPaddedPLVCount() const {
    return GetPaddedNodeCount() * GetPLVCountPerNode();
  }
  size_t GetAllocatedPLVCount() const {
    return GetAllocatedNodeCount() * GetPLVCountPerNode();
  }

  void SetNodeCount(const size_t node_count) { node_count_ = node_count; }
  void SetTempNodeCount(const size_t node_padding) { node_padding_ = node_padding; }
  void SetAllocatedNodeCount(const size_t node_alloc) { node_alloc_ = node_alloc; }

  // ** Resize

  // Resize PVHandler to accomodate DAG with given number of nodes.
  void Resize(const size_t new_node_count, const size_t new_node_alloc) {
    const size_t old_plv_count = GetPLVCount();
    node_count_ = new_node_count;
    node_alloc_ = new_node_alloc;
    // Allocate mmapped data block.
    mmapped_master_plvs_.Resize(GetAllocatedPLVCount() * pattern_count_);
    // Subdivide mmapped data in individual PLVs.
    plvs_ = mmapped_master_plvs_.Subdivide(GetAllocatedPLVCount());
    // Initialize new work space.
    Assert((plvs_.back().rows() == MmappedNucleotidePLV::base_count_) &&
               (plvs_.back().cols() == static_cast<Eigen::Index>(pattern_count_)) &&
               (size_t(plvs_.size()) == GetAllocatedPLVCount()),
           "Didn't get the right shape of PLVs out of Subdivide.");
    for (size_t i = old_plv_count; i < GetPaddedPLVCount(); i++) {
      plvs_.at(i).setZero();
    }
  }

  // Reindex PV according to plv_reindexer.
  void Reindex(const Reindexer plv_reindexer) {
    Reindexer::ReindexInPlace(plvs_, plv_reindexer, GetPLVCount(),
                              GetPLV(GetPLVCount()), GetPLV(GetPLVCount() + 1));
  }

  // Expand node_reindexer into plv_reindexer.
  Reindexer BuildPLVReindexer(const Reindexer &node_reindexer,
                              const size_t old_node_count,
                              const size_t new_node_count) {
    node_count_ = new_node_count;
    Reindexer plv_reindexer(new_node_count * plv_count_per_node_);
    size_t new_plvs_idx = old_node_count * plv_count_per_node_;
    for (size_t i = 0; i < new_node_count; i++) {
      const size_t old_node_idx = i;
      const size_t new_node_idx = node_reindexer.GetNewIndexByOldIndex(old_node_idx);
      for (const auto plv_type : typename PVTypeEnum::Iterator()) {
        // Either get input plv_index from old plvs, or get new plv_index (new data is
        // irrelevant, so just get next available index).
        size_t old_plv_idx;
        if (old_node_idx < old_node_count) {
          old_plv_idx = GetPLVIndex(plv_type, old_node_idx, old_node_count);
        } else {
          old_plv_idx = new_plvs_idx;
          new_plvs_idx++;
        }
        const size_t new_plv_idx = GetPLVIndex(plv_type, new_node_idx, new_node_count);
        plv_reindexer.SetReindex(old_plv_idx, new_plv_idx);
      }
    }
    Assert(plv_reindexer.IsValid(GetPLVCount()), "PLV Reindexer is not valid.");
    return plv_reindexer;
  }

  // ** Access

  // Get vector of Partial Vectors.
  NucleotidePLVRefVector &GetPLVs() { return plvs_; }
  const NucleotidePLVRefVector &GetPLVs() const { return plvs_; }
  // Get PLV by index from the vector of Partial Vectors.
  NucleotidePLVRef &GetPLV(const size_t plv_id) { return plvs_.at(plv_id); };
  const NucleotidePLVRef &GetPLV(const size_t plv_id) const {
    return plvs_.at(plv_id);
  };
  // Get Temporary PLV by index from the vector of Partial Vectors.
  NucleotidePLVRef &GetTempPLV(const size_t plv_id) {
    return plvs_.at(GetTempPLVIndex(plv_id));
  };
  const NucleotidePLVRef &GetTempPLV(const size_t plv_id) const {
    return plvs_.at(GetTempPLVIndex(plv_id));
  };

  // Get total offset into PLVs, indexed based on size of underlying DAG.
  static size_t GetPLVIndex(const PVType plv_type, const size_t node_idx,
                            const size_t node_count) {
    // return GetPLVIndex(PVTypeEnum::GetIndex(plv_type), node_idx, node_count);
    return GetPLVIndex(GetPLVTypeIndex(plv_type), node_idx, node_count);
  };

  size_t GetTempPLVIndex(const size_t plv_id) const {
    const size_t plv_scratch_size = GetPaddedPLVCount() - GetPLVCount();
    Assert(plv_id < plv_scratch_size,
           "Requested temporary plv_id outside of allocated scratch space.");
    return plv_id + GetPLVCount();
  }

  // Get vector of all node ids for given node.
  static SizeVector GetPLVIndexVectorForNodeId(const size_t node_idx,
                                               const size_t node_count) {
    SizeVector plv_idxs;
    for (const auto plv_type : PVTypeEnum::Iterator()) {
      plv_idxs.push_back(GetPLVIndex(plv_type, node_idx, node_count));
    }
    return plv_idxs;
  };

 protected:
  // Get total offset into PLVs.
  static size_t GetPLVIndex(const size_t plv_type_idx, const size_t node_idx,
                            const size_t node_count) {
    return (plv_type_idx * node_count) + node_idx;
  };
  // Get index for given PLV enum.
  static size_t GetPLVTypeIndex(const PLVType plv_type) {
    return static_cast<typename std::underlying_type<PVType>::type>(plv_type);
  };

  // ** Data Sizing
  // "Count" is the currently occupied by data.
  // "Padding" is the amount of free working space added to end of occupied space.
  // "Alloc" is the total current memory allocation.
  // "Resizing factor" is the amount of extra storage allocated for when resizing.

  // Number of nodes in DAG.
  size_t node_count_ = 0;
  // Number of nodes of additional padding for temporary graft nodes in DAG.
  size_t node_padding_ = 2;
  // Number of nodes allocated for in PLVHandler.
  size_t node_alloc_ = 0;
  // Size of Site Pattern.
  size_t pattern_count_ = 0;
  // Number of PLVs for each node in DAG.
  const size_t plv_count_per_node_ = PVTypeEnum::Count;
  // When size exceeds current allocation, ratio to grow new allocation.
  double resizing_factor_ = 2.0;

  // File path to data map.
  std::string mmap_file_path_;
  // Master PLV: Large data block of virtual memory for Partial Likelihood Vectors.
  // Subdivided into sections for plvs_.
  MmappedNucleotidePLV mmapped_master_plvs_;
  // Partial Likelihood Vectors.
  // Divides mmapped_master_plvs_.
  // For example, GP PLVs are divided in the following:
  // - [0, num_nodes): p(s).
  // - [num_nodes, 2*num_nodes): phat(s_right).
  // - [2*num_nodes, 3*num_nodes): phat(s_left).
  // - [3*num_nodes, 4*num_nodes): rhat(s_right) = rhat(s_left).
  // - [4*num_nodes, 5*num_nodes): r(s_right).
  // - [5*num_nodes, 6*num_nodes): r(s_left).
  NucleotidePLVRefVector plvs_;
};

// PLVHandler: Partial Likelihood Vector Handler
class PLVHandler
    : public PVHandler<PartialVectorType::PLVType, PartialVectorType::PLVTypeEnum> {
 public:
  using PLVTypeEnum = PartialVectorType::PLVTypeEnum;
  using PLVType = PLVTypeEnum::Type;
  using PLVTypeIterator = PLVTypeEnum::Iterator;
  static const inline size_t plv_count_ = PLVTypeEnum::Count;

  PLVHandler(const std::string &mmap_file_path, const size_t node_count,
             const size_t pattern_count, const double resizing_factor)
      : PVHandler<PartialVectorType::PLVType, PartialVectorType::PLVTypeEnum>(
            mmap_file_path, node_count, pattern_count, PLVTypeEnum::Count,
            resizing_factor){};

  static Type RPLVType(const bool is_on_left) {
    return is_on_left ? PLVType::RLeft : PLVType::RRight;
  };

  static Type PPLVType(const bool is_on_left) {
    return is_on_left ? PLVType::PHatLeft : PLVType::PHatRight;
  };
};

// PSVHandler: Partial Sankoff Vector Handler
class PSVHandler
    : public PVHandler<PartialVectorType::PSVType, PartialVectorType::PSVTypeEnum> {
 public:
  using PSVTypeEnum = PartialVectorType::PSVTypeEnum;
  using PSVType = PSVTypeEnum::Type;
  using PSVTypeIterator = PSVTypeEnum::Iterator;
  static const inline size_t psv_count_ = PartialVectorType::PSVTypeEnum::Count;

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
