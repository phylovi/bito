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

  static inline const std::string Prefix = "PLV";
  static inline const Array<std::string> Labels = {
      {"P", "PHatRight", "PHatLeft", "RHat", "RRight", "RLeft"}};

  static std::string ToString(const PLVType e) {
    std::stringstream ss;
    ss << Prefix << "::" << Labels[e];
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

  static inline const std::string Prefix = "PSV";
  static inline const Array<std::string> Labels = {{"PRight", "PLeft", "Q"}};

  static std::string ToString(const PSVType e) {
    std::stringstream ss;
    ss << Prefix << "::" << Labels[e];
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
  using PVIdArray = typename TypeEnum::template Array<PVId>;

  PartialVectorHandler(const std::string &mmap_file_path, const size_t elem_count,
                       const size_t pattern_count, const double resizing_factor = 2.0)
      : element_count_(elem_count),
        pattern_count_(pattern_count),
        resizing_factor_(resizing_factor),
        mmap_file_path_(mmap_file_path),
        mmapped_master_pvs_(mmap_file_path_, (elem_count + element_spare_count_) *
                                                 pv_count_per_element_ *
                                                 size_t(resizing_factor_) *
                                                 pattern_count) {
    pv_reindexer_ = Reindexer::IdentityReindexer(GetPaddedPVCount());
    reindexer_init_size_ = pv_reindexer_.size();
  }

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

  void SetCount(const size_t elem_count) { element_count_ = elem_count; }
  void SetSpareCount(const size_t element_spare_count) {
    element_spare_count_ = element_spare_count;
  }
  void SetAllocatedCount(const size_t element_alloc) { element_alloc_ = element_alloc; }

  // ** Resize

  // Resize PVHandler to accomodate DAG with given number of nodes.
  void Resize(const size_t new_elem_count, const size_t new_element_alloc,
              std::optional<size_t> new_element_spare = std::nullopt);
  // Reindex PV according to pv_reindexer.
  void Reindex(const Reindexer pv_reindexer);
  // Reindex PVs by moving data to align with reindexer by copying.
  void ReindexViaMoveCopy(const Reindexer pv_reindexer);
  // Reindex PVs by updating the map from pv_id to data index.
  void ReindexViaRemap(const Reindexer pv_reindexer);
  // Expand element_reindexer into pv_reindexer.
  Reindexer BuildPVReindexer(const Reindexer &element_reindexer,
                             const size_t old_elem_count, const size_t new_elem_count);

  // ** Access

  // Get vector of all Partial Vectors.
  NucleotidePLVRefVector &GetPVs() { return pvs_; }
  const NucleotidePLVRefVector &GetPVs() const { return pvs_; }
  // Get PV by absolute index from the vector of Partial Vectors.
  NucleotidePLVRef &GetPV(const PVId pv_id) {
    auto &pv = pvs_.at(pv_reindexer_.GetOldIndexByNewIndex(pv_id.value_));
    Assert(pv_id < GetAllocatedPVCount(), "pv_id outside valid range.");
    return pv;
  }
  const NucleotidePLVRef &GetPV(const PVId pv_id) const {
    auto &pv = pvs_.at(pv_reindexer_.GetOldIndexByNewIndex(pv_id.value_));
    Assert(pv_id < GetAllocatedPVCount(), "pv_id outside valid range.");
    return pv;
  }
  NucleotidePLVRef &operator()(const PVId pv_id) { return GetPV(pv_id); }
  const NucleotidePLVRef &operator()(const PVId pv_id) const { return GetPV(pv_id); }
  // Get PV by PV type and node index from the vector of Partial Vectors.
  NucleotidePLVRef &GetPV(const PVType pv_type, const DAGElementId elem_id) {
    return GetPV(GetPVIndex(pv_type, elem_id));
  }
  const NucleotidePLVRef &GetPV(const PVType pv_type,
                                const DAGElementId elem_id) const {
    return GetPV(GetPVIndex(pv_type, elem_id));
  }
  NucleotidePLVRef &operator()(const PVType pv_type, const DAGElementId elem_id) {
    return GetPV(GetPVIndex(pv_type, elem_id));
  }
  const NucleotidePLVRef &operator()(const PVType pv_type,
                                     const DAGElementId elem_id) const {
    return GetPV(GetPVIndex(pv_type, elem_id));
  }
  // Get Spare PV by index from the vector of Partial Vectors.
  NucleotidePLVRef &GetSparePV(const PVId pv_id) {
    return GetPV(GetSparePVIndex(pv_id));
  }
  const NucleotidePLVRef &GetSparePV(const PVId pv_id) const {
    return GetPV(GetSparePVIndex(pv_id));
  }

  // Get total offset into PVs, indexed based on underlying DAG.
  static PVId GetPVIndex(const PVType pv_type, const DAGElementId elem_id,
                         const size_t elem_count) {
    return GetPVIndex(TypeEnum::GetIndex(pv_type), elem_id, elem_count);
  }
  PVId GetPVIndex(const PVType pv_type, const DAGElementId elem_id) const {
    Assert(elem_id.value_ < GetCount(), "Requested elem_id is out-of-range.");
    return GetPVIndex(pv_type, elem_id, GetCount());
  }

  // Get PVIndex as a pair of PVType and DAGElementId.
  std::pair<PVType, DAGElementId> GetReversePVIdIndex(const PVId pv_id) const {
    Assert(pv_id.value_ < GetPVCount(), "Requested pv_id is out-of-range.");
    DAGElementId elm_id = (pv_id.value_ % GetCount());
    PVType pv_type = PVType(pv_id.value_ / GetCount());
    return {pv_type, elm_id};
  }

  // Get total offset into temporary PVs, indexed based on underlying grafted DAG.
  PVId GetSparePVIndex(const PVId pv_id) const {
    const size_t pv_scratch_size = GetPaddedPVCount() - GetPVCount();
    Assert(pv_id < pv_scratch_size,
           "Requested temporary pv_id outside of allocated scratch space.");
    return PVId(pv_id.value_ + GetPVCount());
  }
  PVId GetSparePVIndex(const PVType pv_type, const DAGElementId elem_id) {
    Assert(elem_id.value_ < GetSpareCount(),
           "Requested spare elem_id is out-of-range.");
    PVId spare_pv_id = GetPVIndex(pv_type, elem_id, GetSpareCount());
    return GetSparePVIndex(spare_pv_id);
  }

  // PV Reindexer, which serves as the data map to sort PV data.
  const Reindexer &GetPVReindexer() const { return pv_reindexer_; }
  // Get array of all pv_ids for given node.
  PVIdArray GetPVIdArray(const DAGElementId elem_id) const {
    PVIdArray pv_array;
    for (auto pv_type : typename PVTypeEnum::Iterator()) {
      pv_array[pv_type] = GetPVIndex(pv_type, elem_id);
    }
    return pv_array;
  }

  void SetUseRemapping(bool use_remapping) { use_remapping_ = use_remapping; }
  bool GetUseRemapping() const { return use_remapping_; }

  // ** PV Operations

  std::pair<double, double> ValueRange(const PVId pvid) const {
    const auto &pv = GetPV(pvid);
    double max_value = -INFINITY;
    double min_value = INFINITY;
    for (int i = 0; i < pv.rows(); i++) {
      for (int j = 0; j < pv.cols(); j++) {
        double value = pv(i, j);
        if (max_value < value) {
          max_value = value;
        }
        if (min_value > value) {
          min_value = value;
        }
      }
    }
    return {min_value, max_value};
  }

  // Element-wise PV unary operations.
  using UnaryFunction = std::function<double(const double)>;
  void ApplyUnaryOperation(const PVId dest_pvid, const PVId src_pvid_a,
                           UnaryFunction una_fn) {
    auto &dest_pv = GetPV(dest_pvid);
    const auto &src_pv_a = GetPV(src_pvid_a);
    return ApplyUnaryOperation(dest_pv, src_pv_a, una_fn);
  }
  static void ApplyUnaryOperation(NucleotidePLVRef &dest_pv,
                                  const NucleotidePLVRef &src_pv_a,
                                  UnaryFunction una_fn) {
    for (int i = 0; i < src_pv_a.rows(); i++) {
      for (int j = 0; j < src_pv_a.cols(); j++) {
        dest_pv(i, j) = una_fn(src_pv_a(i, j));
      }
    }
  }

  // Element-wise unary operations.
  void AbsValue(const PVId dest_pvid, const PVId src_pvid_a) {
    auto AbsValueFunc = [](const double src) { return abs(src); };
    return ApplyUnaryOperation(dest_pvid, src_pvid_a, AbsValueFunc);
  }

  // Element-wise PV binary operations.
  using BinaryFunction = std::function<double(const double, const double)>;
  void ApplyBinaryOperation(const PVId dest_pvid, const PVId src_pvid_a,
                            const PVId src_pvid_b, BinaryFunction bin_fn) {
    auto &dest_pv = GetPV(dest_pvid);
    const auto &src_pv_a = GetPV(src_pvid_a);
    const auto &src_pv_b = GetPV(src_pvid_b);
    return ApplyBinaryOperation(dest_pv, src_pv_a, src_pv_b, bin_fn);
  }
  static void ApplyBinaryOperation(NucleotidePLVRef &dest_pv,
                                   const NucleotidePLVRef &src_pv_a,
                                   const NucleotidePLVRef &src_pv_b,
                                   BinaryFunction bin_fn) {
    for (int i = 0; i < src_pv_a.rows(); i++) {
      for (int j = 0; j < src_pv_b.cols(); j++) {
        dest_pv(i, j) = bin_fn(src_pv_a(i, j), src_pv_b(i, j));
      }
    }
  }

  // Element-wise binary operations.
  void Add(const PVId dest_pvid, const PVId src_pvid_a, const PVId src_pvid_b) {
    auto AddFunc = [](const double src_a, const double src_b) {
      return (src_a + src_b);
    };
    return ApplyBinaryOperation(dest_pvid, src_pvid_a, src_pvid_b, AddFunc);
  }
  void Subtract(const PVId dest_pvid, const PVId src_pvid_a, const PVId src_pvid_b) {
    auto SubtractFunc = [](const double src_a, const double src_b) {
      return (src_a - src_b);
    };
    return ApplyBinaryOperation(dest_pvid, src_pvid_a, src_pvid_b, SubtractFunc);
  }
  void AbsDiff(const PVId dest_pvid, const PVId src_pvid_a, const PVId src_pvid_b) {
    auto AbsDiffFunc = [](const double src_a, const double src_b) {
      return abs(src_a - src_b);
    };
    return ApplyBinaryOperation(dest_pvid, src_pvid_a, src_pvid_b, AbsDiffFunc);
  }
  void Multiply(const PVId dest_pvid, const PVId src_pvid_a, const PVId src_pvid_b) {
    auto MultiplyFunc = [](const double src_a, const double src_b) {
      return (src_a * src_b);
    };
    return ApplyBinaryOperation(dest_pvid, src_pvid_a, src_pvid_b, MultiplyFunc);
  }
  void Divide(const PVId dest_pvid, const PVId src_pvid_a, const PVId src_pvid_b) {
    auto DivideFunc = [](const double src_a, const double src_b) {
      return (src_a / src_b);
    };
    return ApplyBinaryOperation(dest_pvid, src_pvid_a, src_pvid_b, DivideFunc);
  }

  // Find the maximum element-wise absolute difference between two PVs.
  static double MaxDifference(const NucleotidePLVRef &pv_a,
                              const NucleotidePLVRef &pv_b) {
    double max_diff = 0;
    for (int i = 0; i < pv_a.rows(); i++) {
      for (int j = 0; j < pv_a.cols(); j++) {
        double diff = abs(pv_a(i, j) - pv_b(i, j));
        if (diff > max_diff) {
          max_diff = diff;
        }
      }
    }
    return max_diff;
  }
  double MaxDifference(const PVId pvid_a, const PVId pvid_b) const {
    const auto &pv_a = GetPV(pvid_a);
    const auto &pv_b = GetPV(pvid_b);
    return MaxDifference(pv_a, pv_b);
  }

  double Min(const PVId pvid) const { return GetPV(pvid).minCoeff(); }
  double Max(const PVId pvid) const { return GetPV(pvid).maxCoeff(); }

  // ** I/O

  // Output data to string.
  std::string ToString(const PVId pv_id, const bool show_labels = false) const {
    if (pv_id < GetPVCount()) {
      const auto &[pv_type, elem_id] = GetReversePVIdIndex(pv_id);
      return ToString(pv_type, elem_id, show_labels);
    }
    std::stringstream out;
    out << "PV[" << pv_id << "]: " << std::endl;
    for (auto &&row : GetPV(pv_id).rowwise()) {
      out << row << std::endl;
    }
    return out.str();
  }
  std::string ToString(const PVType pv_type, const DAGElementId elem_id,
                       const bool show_labels = false) const {
    std::stringstream out;
    out << "PV[" << PVTypeEnum::ToString(pv_type) << ", Element" << elem_id << ", PV"
        << GetPVIndex(pv_type, elem_id) << "]: " << std::endl;
    out << ToString(GetPV(pv_type, elem_id));
    return out.str();
  }
  std::string AllPVsToString(const bool show_labels = false) const {
    std::stringstream out;
    for (const auto pv_type : typename PVTypeEnum::Iterator()) {
      for (DAGElementId elem_id = 0; elem_id < GetCount(); elem_id++) {
        out << ToString(pv_type, elem_id, show_labels);
      }
    }
    return out.str();
  }
  size_t ToHash(const PVId pv_id) const { return ToHash(GetPV(pv_id)); }
  std::string ToHashString(const PVId pv_id, const size_t length = 16) const {
    return ToHashString(GetPV(pv_id), length);
  }
  DoubleVector ToDoubleVector(const PVId pv_id) const {
    return ToDoubleVector(GetPV(pv_id));
  }

  static std::string ToString(const NucleotidePLVRef &pv) {
    std::stringstream out;
    for (int i = 0; i < pv.rows(); i++) {
      out << "[";
      for (int j = 0; j < pv.cols(); j++) {
        out << pv(i, j) << ((j < (pv.cols() - 1)) ? ", " : "");
      }
      out << "]" << std::endl;
    }
    return out.str();
  }
  static size_t ToHash(const NucleotidePLVRef &pv) {
    size_t seed = pv.rows() * pv.cols();
    for (int i = 0; i < pv.rows(); i++) {
      for (int j = 0; j < pv.cols(); j++) {
        seed ^= std::hash<double>()(pv(i, j)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
    }
    return seed;
  }
  static std::string ToHashString(const NucleotidePLVRef &pv,
                                  const size_t length = 16) {
    return HashToString(ToHash(pv), length);
  }
  static DoubleVector ToDoubleVector(const NucleotidePLVRef &pv) {
    DoubleVector values;
    for (int i = 0; i < pv.rows(); i++) {
      for (int j = 0; j < pv.cols(); j++) {
        values.push_back(pv(i, j));
      }
    }
    return values;
  }

  // ** Miscellaneous

  std::unordered_map<PVId, std::pair<PVType, DAGElementId>> BuildPVIdMap() const {
    std::unordered_map<PVId, std::pair<PVType, DAGElementId>> pvid_map;
    for (DAGElementId elem_id = 0; elem_id < GetCount(); elem_id++) {
      for (PVType pv_type : typename PVTypeEnum::Iterator()) {
        const auto pv_id = GetPVIndex(pv_type, elem_id);
        pvid_map[pv_id] = {pv_type, elem_id};
      }
    }
    return pvid_map;
  }

  static int Compare(const PartialVectorHandler<PVTypeEnum, DAGElementId> &pv_lhs,
                     const PartialVectorHandler<PVTypeEnum, DAGElementId> &pv_rhs,
                     const bool is_quiet = true) {
    std::stringstream dev_null;
    std::ostream &os = (is_quiet ? dev_null : std::cerr);

    bool pv_count_match = (pv_lhs.GetPVCount() == pv_rhs.GetPVCount());
    if (!pv_count_match) {
      os << "PVHandler::Compare: PV Counts do not match." << std::endl;
      return false;
    }
    bool pvs_match = true;
    for (PVId pv_id(0); pv_id < pv_lhs.GetPVCount(); pv_id++) {
      bool pv_match = (pv_lhs.GetPV(pv_id) == pv_rhs.GetPV(pv_id));
      pvs_match &= pv_match;
      if (!pvs_match) {
        os << "PVHandler::Compare: PVs do not match at PV" << pv_id << std::endl;
        os << "PV_LHS:" << std::endl << pv_lhs.ToString(pv_id, true) << std::endl;
        os << "PV_RHS:" << std::endl << pv_rhs.ToString(pv_id, true) << std::endl;
        return false;
      }
    }
    return true;
  }

 protected:
  // Get total offset into PVs.
  static PVId GetPVIndex(const size_t pv_type_id, const DAGElementId elem_id,
                         const size_t elem_count) {
    return (pv_type_id * elem_count) + elem_id.value_;
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
  size_t element_spare_count_ = 16;
  // Number of nodes allocated for in PVHandler.
  size_t element_alloc_ = 0;
  // Size of Site Pattern.
  size_t pattern_count_ = 0;
  // Number of PVs for each node in DAG.
  size_t pv_count_per_element_ = PVTypeEnum::Count;
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

  // Reindex map for finding pv locations.
  Reindexer pv_reindexer_;
  size_t reindexer_init_size_ = 0;
  // Whether to use remapping to reindex PLVs, otherwise only use
  bool use_remapping_ = true;
  double reindex_ratio = 10;
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

  struct PVIdSet {
    DAGElementId id;
    PVId p_pvid;
    PVId phatright_pvid;
    PVId phatleft_pvid;
    PVId rhat_pvid;
    PVId rright_pvid;
    PVId rleft_pvid;
  };

  PLVHandler(const std::string &mmap_file_path, const size_t elem_count,
             const size_t pattern_count, const double resizing_factor = 2.0)
      : PartialVectorHandler<PLVTypeEnum, DAGElementId>(
            mmap_file_path, elem_count, pattern_count, resizing_factor) {}

  PVIdSet BuildPVIdSet(const DAGElementId elem_id) {
    PVIdSet result;
    result.id = elem_id;
    result.p_pvid = GetPVIndex(PLVType::P, elem_id);
    result.phatright_pvid = GetPVIndex(PLVType::PHatRight, elem_id);
    result.phatleft_pvid = GetPVIndex(PLVType::PHatLeft, elem_id);
    result.rhat_pvid = GetPVIndex(PLVType::RHat, elem_id);
    result.rright_pvid = GetPVIndex(PLVType::RRight, elem_id);
    result.rleft_pvid = GetPVIndex(PLVType::RLeft, elem_id);
    return result;
  }
};

using PLVNodeHandler = PLVHandler<NodeId>;
using PLVEdgeHandler = PLVHandler<EdgeId>;

// PSVHandler: Partial Sankoff Vector Handler
template <class DAGElementId>
class PSVHandler
    : public PartialVectorHandler<PartialVectorType::PSVTypeEnum, DAGElementId> {
 public:
  using PSVType = PartialVectorType::PSVType;
  using PSVTypeEnum = PartialVectorType::PSVTypeEnum;
  using PSVTypeIterator = PSVTypeEnum::Iterator;
  static const inline size_t psv_count_ = PartialVectorType::PSVTypeEnum::Count;

  struct PVIdSet {
    EdgeId edge_id;
    PVId pright_pvid;
    PVId pleft_pvid;
    PVId q_pvid;
  };

  PSVHandler(const std::string &mmap_file_path, const size_t elem_count,
             const size_t pattern_count, const double resizing_factor = 2.0)
      : PartialVectorHandler<PSVTypeEnum, DAGElementId>(
            mmap_file_path, elem_count, pattern_count, resizing_factor) {}

  PVIdSet BuildPVIdSet(const DAGElementId elem_id) {
    PVIdSet result;
    result.id = elem_id;
    result.pright_pvid = GetPVIndex(PSVType::PRight, elem_id);
    result.pleft_pvid = GetPVIndex(PSVType::PLeft, elem_id);
    result.q_pvid = GetPVIndex(PSVType::Q, elem_id);
    return result;
  }
};

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
