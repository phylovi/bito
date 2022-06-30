// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "plv_handler.hpp"

// template <class PVType, class PVTypeEnum>
// PVHandler::PVHandler(const std::string &mmap_file_path, const size_t node_count,
//                      const size_t pattern_count, const size_t plv_count_per_node,
//                      const double resizing_factor)
//     : pattern_count_(pattern_count),
//       node_count_(node_count),
//       resizing_factor_(resizing_factor),
//       mmap_file_path_(mmap_file_path),
//       mmapped_master_plvs_(mmap_file_path_,
//                            (node_count + node_padding_) * plv_count_per_node_ *
//                                size_t(resizing_factor_) * pattern_count) {}

// template <class PVType, class PVTypeEnum>
// double PVHandler<PVType, PVTypeEnum>::GetByteCount() const {
//   return mmapped_master_plvs_.ByteCount();
// }

//   // ** Counts

//   double GetByteCount() const { return mmapped_master_plvs_.ByteCount(); }
//   size_t GetPLVCountPerNode() const { return plv_count_per_node_; }
//   size_t GetSitePatternCount() const { return pattern_count_; }
//   size_t GetNodeCount() const { return node_count_; }
//   size_t GetTempNodeCount() const { return node_padding_; }
//   size_t GetPaddedNodeCount() const { return node_count_ + node_padding_; }
//   size_t GetAllocatedNodeCount() const { return node_alloc_; }
//   size_t GetPLVCount() const { return GetNodeCount() * GetPLVCountPerNode(); }
//   size_t GetTempPLVCount() const { return GetTempNodeCount() * GetPLVCountPerNode();
//   } size_t GetPaddedPLVCount() const {
//     return GetPaddedNodeCount() * GetPLVCountPerNode();
//   }
//   size_t GetAllocatedPLVCount() const {
//     return GetAllocatedNodeCount() * GetPLVCountPerNode();
//   }

//   // ** Resize

//   // Resize PVHandler to accomodate DAG with given number of nodes.
//   void Resize(const size_t new_node_count, const size_t new_node_alloc) {
//     const size_t old_plv_count = GetPLVCount();
//     node_count_ = new_node_count;
//     node_alloc_ = new_node_alloc;
//     // Allocate mmapped data block.
//     mmapped_master_plvs_.Resize(GetAllocatedPLVCount() * pattern_count_);
//     // Subdivide mmapped data in individual PLVs.
//     plvs_ = mmapped_master_plvs_.Subdivide(GetAllocatedPLVCount());
//     // Initialize new work space.
//     Assert((plvs_.back().rows() == MmappedNucleotidePLV::base_count_) &&
//                (plvs_.back().cols() == static_cast<Eigen::Index>(pattern_count_)) &&
//                (size_t(plvs_.size()) == GetAllocatedPLVCount()),
//            "Didn't get the right shape of PLVs out of Subdivide.");
//     for (size_t i = old_plv_count; i < GetPaddedPLVCount(); i++) {
//       plvs_.at(i).setZero();
//     }
//   }

//   // Reindex PV according to plv_reindexer.
//   void Reindex(const Reindexer plv_reindexer) {
//     Reindexer::ReindexInPlace(plvs_, plv_reindexer, GetPLVCount(),
//                               GetPLV(GetPLVCount()), GetPLV(GetPLVCount() + 1));
//   }

//   // Expand node_reindexer into plv_reindexer.
//   Reindexer BuildPLVReindexer(const Reindexer &node_reindexer,
//                               const size_t old_node_count,
//                               const size_t new_node_count) {
//     node_count_ = new_node_count;
//     Reindexer plv_reindexer(new_node_count * plv_count_per_node_);
//     size_t new_plvs_idx = old_node_count * plv_count_per_node_;
//     for (size_t i = 0; i < new_node_count; i++) {
//       const size_t old_node_idx = i;
//       const size_t new_node_idx = node_reindexer.GetNewIndexByOldIndex(old_node_idx);
//       for (const auto plv_type : typename PVTypeEnum::Iterator()) {
//         // Either get input plv_index from old plvs, or get new plv_index (new data
//         is
//         // irrelevant, so just get next available index).
//         size_t old_plv_idx;
//         if (old_node_idx < old_node_count) {
//           old_plv_idx = GetPLVIndex(plv_type, old_node_idx, old_node_count);
//         } else {
//           old_plv_idx = new_plvs_idx;
//           new_plvs_idx++;
//         }
//         const size_t new_plv_idx = GetPLVIndex(plv_type, new_node_idx,
//         new_node_count); plv_reindexer.SetReindex(old_plv_idx, new_plv_idx);
//       }
//     }
//     Assert(plv_reindexer.IsValid(GetPLVCount()), "PLV Reindexer is not valid.");
//     return plv_reindexer;
//   }

//   // ** Access

//   // Get vector of Partial Vectors.
//   NucleotidePLVRefVector &GetPLVs() { return plvs_; }
//   const NucleotidePLVRefVector &GetPLVs() const { return plvs_; }
//   // Get PLV by index from the vector of Partial Vectors.
//   NucleotidePLVRef &GetPLV(const size_t plv_id) { return plvs_.at(plv_id); };
//   const NucleotidePLVRef &GetPLV(const size_t plv_id) const {
//     return plvs_.at(plv_id);
//   };
//   // Get Temporary PLV by index from the vector of Partial Vectors.
//   NucleotidePLVRef &GetTempPLV(const size_t plv_id) {
//     return plvs_.at(GetTempPLVIndex(plv_id));
//   };
//   const NucleotidePLVRef &GetTempPLV(const size_t plv_id) const {
//     return plvs_.at(GetTempPLVIndex(plv_id));
//   };

//   // Get total offset into PLVs, indexed based on size of underlying DAG.
//   static size_t GetPLVIndex(const PVType plv_type, const size_t node_idx,
//                             const size_t node_count) {
//     // return GetPLVIndex(PVTypeEnum::GetIndex(plv_type), node_idx, node_count);
//     return GetPLVIndex(GetPLVTypeIndex(plv_type), node_idx, node_count);
//   };

//   size_t GetTempPLVIndex(const size_t plv_id) const {
//     const size_t plv_scratch_size = GetPaddedPLVCount() - GetPLVCount();
//     Assert(plv_id < plv_scratch_size,
//            "Requested temporary plv_id outside of allocated scratch space.");
//     return plv_id + GetPLVCount();
//   }

//   // Get vector of all node ids for given node.
//   static SizeVector GetPLVIndexVectorForNodeId(const size_t node_idx,
//                                                const size_t node_count) {
//     SizeVector plv_idxs;
//     for (const auto plv_type : PVTypeEnum::Iterator()) {
//       plv_idxs.push_back(GetPLVIndex(plv_type, node_idx, node_count));
//     }
//     return plv_idxs;
//   };

//  protected:
//   // Get total offset into PLVs.
//   static size_t GetPLVIndex(const size_t plv_type_idx, const size_t node_idx,
//                             const size_t node_count) {
//     return (plv_type_idx * node_count) + node_idx;
//   };
//   // Get index for given PLV enum.
//   static size_t GetPLVTypeIndex(const PLVType plv_type) {
//     return static_cast<typename std::underlying_type<PVType>::type>(plv_type);
//   };

//   // Size of Site Pattern.
//   size_t pattern_count_ = 0;
//   // Number of nodes in DAG.
//   size_t node_count_ = 0;
//   // Number of nodes allocated for in PLVHandler.
//   size_t node_alloc_ = 0;
//   // Number of nodes of additional padding for temporary graft nodes in DAG.
//   size_t node_padding_ = 2;
//   // Number of PLVs for each node in DAG.
//   const size_t plv_count_per_node_ = PVTypeEnum::Count;
//   // When size exceeds current allocation, ratio to grow new allocation.
//   double resizing_factor_ = 2.0;

//   // File path to data map.
//   std::string mmap_file_path_;
//   // Master PLV: Large data block of virtual memory for Partial Likelihood Vectors.
//   // Subdivided into sections for plvs_.
//   MmappedNucleotidePLV mmapped_master_plvs_;
//   // Partial Likelihood Vectors.
//   // Divides mmapped_master_plvs_.
//   // For example, GP PLVs are divided in the following:
//   // - [0, num_nodes): p(s).
//   // - [num_nodes, 2*num_nodes): phat(s_right).
//   // - [2*num_nodes, 3*num_nodes): phat(s_left).
//   // - [3*num_nodes, 4*num_nodes): rhat(s_right) = rhat(s_left).
//   // - [4*num_nodes, 5*num_nodes): r(s_right).
//   // - [5*num_nodes, 6*num_nodes): r(s_left).
//   NucleotidePLVRefVector plvs_;
// };

// // PLVHandler: Partial Likelihood Vector Handler
// class PLVHandler
//     : public PVHandler<PartialVectorType::PLVType, PartialVectorType::PLVTypeEnum> {
//  public:
//   using PLVTypeEnum = PartialVectorType::PLVTypeEnum;
//   using PLVType = PLVTypeEnum::Type;
//   using PLVTypeIterator = PLVTypeEnum::Iterator;
//   static const inline size_t plv_count_ = PLVTypeEnum::Count;

//   PLVHandler(const std::string &mmap_file_path, const size_t node_count,
//              const size_t pattern_count, const double resizing_factor)
//       : PVHandler<PartialVectorType::PLVType, PartialVectorType::PLVTypeEnum>(
//             mmap_file_path, node_count, pattern_count, PLVTypeEnum::Count,
//             resizing_factor){};

//   static Type RPLVType(const bool is_on_left) {
//     return is_on_left ? PLVType::RLeft : PLVType::RRight;
//   };

//   static Type PPLVType(const bool is_on_left) {
//     return is_on_left ? PLVType::PHatLeft : PLVType::PHatRight;
//   };
// };

// // PSVHandler: Partial Sankoff Vector Handler
// class PSVHandler
//     : public PVHandler<PartialVectorType::PSVType, PartialVectorType::PSVTypeEnum> {
//  public:
//   using PSVTypeEnum = PartialVectorType::PSVTypeEnum;
//   using PSVType = PSVTypeEnum::Type;
//   using PSVTypeIterator = PSVTypeEnum::Iterator;
//   static const inline size_t psv_count_ = PartialVectorType::PSVTypeEnum::Count;

//   static Type PPLVType(const bool is_on_left) {
//     return is_on_left ? PSVType::PLeft : PSVType::PRight;
//   }
// };

// explicit instatiation of
// template class PVHandler<PartialVectorType::PLVType, PartialVectorType::PLVTypeEnum>;
// template class PVHandler<PartialVectorType::PSVType, PartialVectorType::PSVTypeEnum>;
