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
class PLVTypeEnum : public EnumWrapper<PLVType, PLVCount, PLVType::P, PLVType::RLeft> {
 public:
  static inline const Array<std::string> Labels = {{"PLV::P", "PLV::PHatRight",
                                                    "PLV::PHatLeft", "PLV::RHat",
                                                    "PLV::RRight", "PLV::RLeft"}};

  static std::string ToString(const Type e) {
    std::stringstream os;
    os << "PLV::" << Labels[e];
    return os.str();
  };

  friend std::ostream &operator<<(std::ostream &os, const Type e) {
    os << ToString(e);
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
class PSVTypeEnum : public EnumWrapper<PSVType, PSVCount, PSVType::PRight, PSVType::Q> {
 public:
  static inline const Array<std::string> Labels = {
      {"PSV::PRight", "PSV::PLeft", "PSV::Q"}};

  static std::string ToString(const Type e) {
    std::stringstream os;
    os << "PSV::" << Labels[e];
    return os.str();
  };

  friend std::ostream &operator<<(std::ostream &os, const Type e) {
    os << ToString(e);
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
      : mmap_file_path_(mmap_file_path),
        pattern_count_(pattern_count),
        node_count_(node_count),
        plv_count_per_node_(plv_count_per_node),
        resizing_factor_(resizing_factor),
        mmapped_master_plv_(mmap_file_path_,
                            (node_count + node_padding_) * plv_count_per_node_ *
                                size_t(resizing_factor_) * pattern_count) {}

  // ** Counts

  double PLVByteCount() const { return mmapped_master_plv_.ByteCount(); };
  size_t GetSitePatternCount() const { return pattern_count_; };
  size_t GetNodeCount() const { return node_count_; };
  size_t GetTempNodeCount() const { return node_padding_; }
  size_t GetPaddedNodeCount() const { return node_count_ + node_padding_; };
  size_t GetAllocatedNodeCount() const { return node_alloc_; }
  size_t GetPLVCount() const { return GetNodeCount() * plv_count_per_node_; };
  size_t GetTempPLVCount() const { return GetTempNodeCount() * plv_count_per_node_; };
  size_t GetPaddedPLVCount() const {
    return GetPaddedNodeCount() * plv_count_per_node_;
  };
  size_t GetAllocatedPLVCount() const {
    return GetAllocatedNodeCount() * plv_count_per_node_;
  }

  // ** Resize

  // Resize PVHandler to accomodate DAG with given number of nodes.
  void Resize(const size_t node_alloc) {
    const size_t old_node_count = GetNodeCount();
    const size_t old_plv_count = GetPLVCount();
    node_count_ = new_node_count;
    mmapped_master_plv_.Resize(GetAllocatedPLVCount() * pattern_count_);
    plvs_ = mmapped_master_plv_.Subdivide(GetAllocatedPLVCount());
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
  void GPEngine::Reindex(const Reindexer plv_reindexer, const size_t old_node_count) {
    Reindexer::ReindexInPlace(plvs_, plv_reindexer, GetPLVCount(),
                              plvs_.at(GetPLVCount()), plvs_.at(GetPLVCount() + 1));
  }

  // Expand node_reindexer into plv_reindexer.
  Reindexer BuildPLVReindexer(const Reindexer &node_reindexer) {
    Reindexer plv_reindexer(node_count_ * plv_count_per_node_);
    size_t new_data_idx = old_node_count * plv_count_per_node_;
    for (size_t i = 0; i < node_count_; i++) {
      const size_t old_node_idx = i;
      const size_t new_node_idx = node_reindexer.GetNewIndexByOldIndex(old_node_idx);
      for (const auto plv_type : PVTypeIterator()) {
        // Either get input plv_index from old plvs, or get new plv_index (new data is
        // irrelevant, so just get next available index).
        size_t old_plv_idx;
        if (old_node_idx < old_node_count) {
          old_plv_idx = PLVHandler::GetPLVIndex(plv_type, old_node_idx, old_node_count);
        } else {
          old_plv_idx = new_data_idx;
          new_data_idx++;
        }
        const size_t new_plv_idx =
            PLVHandler::GetPLVIndex(plv_type, new_node_idx, node_count_);
        plv_reindexer.SetReindex(old_plv_idx, new_plv_idx);
      }
    }
    Assert(plv_reindexer.IsValid(GetPLVCount()), "PLV Reindexer is not valid.");
    return plv_reindexer;
  }

  // ** Access

  //
  NucleotidePLVRefVector &Get() { return data_; }
  //
  NucleotidePLVRef &Get(const size_t offset) { return data_.at(offset); };

  // Get total offset into PLVs, indexed based on size of underlying DAG.
  static size_t GetPLVIndex(const PVType plv_type, const size_t node_idx,
                            const size_t node_count) {
    return GetPLVIndex(PVTypeEnum::GetIndex(plv_type), node_idx, node_count);
  };

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

  // File path to data map.
  std::string mmap_file_path_;
  // Master PLV: Large data block of virtual memory for Partial Likelihood Vectors.
  // Subdivided into sections for plvs_.
  MmappedNucleotidePLV mmapped_master_data_;
  // Partial Likelihood Vectors.
  // plvs_ store the following:
  // [0, num_nodes): p(s).
  // [num_nodes, 2*num_nodes): phat(s_right).
  // [2*num_nodes, 3*num_nodes): phat(s_left).
  // [3*num_nodes, 4*num_nodes): rhat(s_right) = rhat(s_left).
  // [4*num_nodes, 5*num_nodes): r(s_right).
  // [5*num_nodes, 6*num_nodes): r(s_left).
  NucleotidePLVRefVector data_;

  // Size of Site Pattern.
  size_t pattern_count_ = 0;
  // Number of nodes in DAG.
  size_t node_count_ = 0;
  // Number of nodes allocated for in PLVHandler.
  size_t node_alloc_ = 0;
  // Number of nodes of additional padding for temporary graft nodes in DAG.
  size_t node_padding_ = 2;
  // Number of PLVs for each node in DAG.
  const size_t plv_count_per_node_ = 6;
  // When size exceeds current allocation, ratio to grow new allocation.
  double resizing_factor_ = 2.0;
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
            mmap_file_path, node_count, pattern_count, PLVTypeEnum::Count){};

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
    return is_on_left ? PSV::PLeft : PSV::PRight;
  }
};

#ifdef DOCTEST_LIBRARY_INCLUDED

// Check that PLV iterator iterates over all PLVs exactly once.
TEST_CASE("PLVHandler: EnumIterator") {
  using namespace PartialVectorType;
  const std::array<PLVType, PLVTypeEnum::Count> plv_types = {
      {PLVType::P, PLVType::PHatRight, PLVType::PHatLeft, PLVType::RHat,
       PLVType::RRight, PLVType::RLeft}};
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

TEST_CASE("PartialVector") {
  std::cout << "PARTIAL_VECTOR" << std::endl;
  using namespace PartialVectorType;

  std::cout << "PLV_TYPE" << std::endl;
  for (const auto plv : PLVTypeEnum::Iterator()) {
    std::cout << "PLVType=>" << PLVTypeEnum::ToString(plv) << ", ";
  }
  std::cout << "total=" << PLVTypeEnum::Count << std::endl;

  std::cout << "PSV_TYPE" << std::endl;
  for (const auto psv : PSVTypeEnum::Iterator()) {
    std::cout << "PSVType=>" << PSVTypeEnum::ToString(psv) << ", ";
  }
  std::cout << "total=" << PSVTypeEnum::Count << std::endl;
}

#endif  // DOCTEST_LIBRARY_INCLUDED
