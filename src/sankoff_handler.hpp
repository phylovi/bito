// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Methods to calculate Sankoff on a tree.

#pragma once

#include "eigen_sugar.hpp"
#include "sugar.hpp"
#include "sankoff_matrix.hpp"
#include "site_pattern.hpp"
#include "node.hpp"
#include "driver.hpp"
#include "pv_handler.hpp"
#include "gp_dag.hpp"

// Partial vector for one node across all sites
using SankoffPartial = NucleotidePLV;
// Each SankoffPartial represents calculations for one node
using SankoffPartialVec = std::vector<SankoffPartial>;
// references for SankoffPartials
using SankoffPartialRef = Eigen::Ref<SankoffPartial>;
using SankoffPartialRefVec = std::vector<SankoffPartialRef>;

class SankoffHandler {
 public:
  // DNA assumption
  static constexpr size_t state_count_ = 4;
  static constexpr double big_double_ = static_cast<double>(INT_MAX);

  // Constructors
  SankoffHandler(SitePattern site_pattern, const std::string &mmap_file_path,
                 double resizing_factor = 2.0)
      : mutation_costs_(SankoffMatrix()),
        site_pattern_(std::move(site_pattern)),
        resizing_factor_(resizing_factor),
        psv_handler_(mmap_file_path, 0, site_pattern_.PatternCount(),
                     resizing_factor_) {
    Assert(site_pattern_.SequenceCount() == site_pattern_.TaxonCount(),
           "Error in SankoffHandler constructor 1: Every sequence should be associated "
           "with a node.");
    psv_handler_.SetCount(site_pattern_.TaxonCount());
    psv_handler_.SetAllocatedCount(
        size_t(ceil(double(psv_handler_.GetPaddedCount()) * resizing_factor_)));
    psv_handler_.Resize(site_pattern_.TaxonCount(), psv_handler_.GetAllocatedCount());
  }

  SankoffHandler(CostMatrix cost_matrix, SitePattern site_pattern,
                 const std::string &mmap_file_path, double resizing_factor = 2.0)
      : mutation_costs_(SankoffMatrix(cost_matrix)),
        site_pattern_(std::move(site_pattern)),
        resizing_factor_(resizing_factor),
        psv_handler_(mmap_file_path, 0, site_pattern_.PatternCount(),
                     resizing_factor_) {
    Assert(site_pattern_.SequenceCount() == site_pattern_.TaxonCount(),
           "Error in SankoffHandler constructor 2: Every sequence should be associated "
           "with a node.");
    psv_handler_.SetCount(site_pattern_.TaxonCount());
    psv_handler_.SetAllocatedCount(
        size_t(ceil(double(psv_handler_.GetPaddedCount()) * resizing_factor_)));
    psv_handler_.Resize(site_pattern_.TaxonCount(), psv_handler_.GetAllocatedCount());
  }

  SankoffHandler(SankoffMatrix sankoff_matrix, SitePattern site_pattern,
                 const std::string &mmap_file_path, double resizing_factor = 2.0)
      : mutation_costs_(std::move(sankoff_matrix)),
        site_pattern_(std::move(site_pattern)),
        resizing_factor_(resizing_factor),
        psv_handler_(mmap_file_path, 0, site_pattern_.PatternCount(),
                     resizing_factor_) {
    Assert(site_pattern_.SequenceCount() == site_pattern_.TaxonCount(),
           "Error in SankoffHandler constructor 3: Every sequence should be associated "
           "with a node.");
    psv_handler_.SetCount(site_pattern_.TaxonCount());
    psv_handler_.SetAllocatedCount(
        size_t(ceil(double(psv_handler_.GetPaddedCount()) * resizing_factor_)));
    psv_handler_.Resize(site_pattern_.TaxonCount(), psv_handler_.GetAllocatedCount());
  }

  PSVNodeHandler &GetPSVHandler() { return psv_handler_; }
  SankoffMatrix &GetCostMatrix() { return mutation_costs_; }

  // Resize PVs to fit model.
  void Resize(const size_t new_node_count);

  // Partial Sankoff Vector Handler.
  SankoffPartialVec PartialsAtPattern(PSVType psv_type, size_t pattern_idx) {
    SankoffPartialVec partials_at_pattern(psv_handler_.GetCount());
    for (NodeId node = 0; node < psv_handler_.GetCount(); node++) {
      partials_at_pattern[node.value_] = psv_handler_(psv_type, node).col(pattern_idx);
    }
    return partials_at_pattern;
  }

  // Fill in leaf-values for P partials.
  void GenerateLeafPartials();

  // Sum p-partials for right and left children of node 'node_id'
  // In this case, we get the full p-partial of the given node after all p-partials
  // have been concatenated into one SankoffPartialVector
  EigenVectorXd TotalPPartial(NodeId node_id, size_t site_idx);

  // Calculate the partial for a given parent-child pair
  EigenVectorXd ParentPartial(EigenVectorXd child_partials);

  // Populate rootward parsimony PV for node.
  void PopulateRootwardParsimonyPVForNode(const NodeId parent_id,
                                          const NodeId left_child_id,
                                          const NodeId right_child_id);
  // Populate leafward parsimony PV for node.
  void PopulateLeafwardParsimonyPVForNode(const NodeId parent_id,
                                          const NodeId left_child_id,
                                          const NodeId right_child_id);

  // Calculates left p_partials, right p_partials, and q_partials for all nodes at all
  // sites in tree topology.
  void RunSankoff(Node::NodePtr topology);

  // Calculates parsimony score on given node across all sites.
  double ParsimonyScore(NodeId node_id = NodeId(0));

 private:
  SankoffMatrix mutation_costs_;
  SitePattern site_pattern_;
  double resizing_factor_;
  PSVNodeHandler psv_handler_;
};
