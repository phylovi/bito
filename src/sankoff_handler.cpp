// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "sankoff_handler.hpp"

void SankoffHandler::GenerateLeafPartials() {
  // first check that the psv_handler has been resized to deal with the leaf labels
  Assert(psv_handler_.GetNodeCount() >= site_pattern_.TaxonCount(),
         "Error in SankoffHandler::GenerateLeafPartials: "
         "psv_handler_ should be initialized to accomodate"
         "the number of leaf nodes in the site_pattern_.");

  // Iterate over all leaf nodes to instantiate each with P partial values
  for (size_t leaf_node = 0; leaf_node < site_pattern_.TaxonCount(); leaf_node++) {
    SankoffPartial node_partials(state_count_, site_pattern_.PatternCount());
    // set leaf node partial to have big_double_ infinity substitute
    node_partials.block(0, 0, state_count_, site_pattern_.PatternCount())
        .fill(big_double_);
    // now fill in appropriate entries of the leaf-partial where non-infinite
    for (size_t pattern_idx = 0; pattern_idx < site_pattern_.PatternCount();
         pattern_idx++) {
      auto site_val = site_pattern_.GetPatternSymbol(leaf_node, pattern_idx);
      if (site_val < state_count_) {
        node_partials(site_val, pattern_idx) = 0.;
      } else if (site_val == state_count_) {
        // Leaves with gaps in sequence and ambiguous nucleotides are assigned sankoff
        // partial vector [0, 0, 0, 0] at the corresponding site.
        node_partials.col(pattern_idx).fill(0);
      } else {
        Failwith(
            "Error in SankoffHandler::GenerateLeafPartials: Invalid nucleotide state "
            "in sequence alignment.");
      }
    }
    psv_handler_.GetPV(PSVType::PLeft, leaf_node) = node_partials;
    psv_handler_.GetPV(PSVType::PRight, leaf_node).fill(0);
  }
}

EigenVectorXd SankoffHandler::ParentPartial(EigenVectorXd child_partials) {
  Assert(child_partials.size() == state_count_,
         "child_partials in SankoffHandler::ParentPartial should have 4 states.");
  EigenVectorXd parent_partials(state_count_);
  parent_partials.setZero();
  for (size_t parent_state = 0; parent_state < state_count_; parent_state++) {
    EigenVectorXd partials_for_state(state_count_);
    for (size_t child_state = 0; child_state < state_count_; child_state++) {
      auto cost = mutation_costs_.GetCost(parent_state, child_state);
      partials_for_state[child_state] = cost + child_partials[child_state];
    }
    auto minimum_element =
        *std::min_element(partials_for_state.data(),
                          partials_for_state.data() + partials_for_state.size());
    parent_partials[parent_state] = minimum_element;
  }

  return parent_partials;
}

EigenVectorXd SankoffHandler::TotalPPartial(size_t node_id, size_t site_idx) {
  return psv_handler_.GetPV(PSVType::PLeft, node_id).col(site_idx) +
         psv_handler_.GetPV(PSVType::PRight, node_id).col(site_idx);
}

void SankoffHandler::RunSankoff(Node::NodePtr topology) {
  psv_handler_.SetNodeCount(topology->Id() + 1);
  psv_handler_.SetAllocatedNodeCount(
      size_t(ceil(double(psv_handler_.GetPaddedNodeCount()) * resizing_factor_)));
  psv_handler_.Resize(topology->Id() + 1, psv_handler_.GetAllocatedNodeCount());

  // fill in leaf node partials for PSV (stored in PLeft in PSVHandler instance)
  GenerateLeafPartials();

  // generating p_partials (right and left)
  topology->Postorder([this](const Node* node) {
    if (!node->IsLeaf()) {
      Assert(node->Children().size() == 2,
             "Error in SankoffHandler::RunSankoff: Tree should be bifurcating.");
      for (size_t pattern_idx = 0; pattern_idx < site_pattern_.PatternCount();
           pattern_idx++) {
        // Which child partial is in right or left doesn't actually matter because they
        // are summed when calculating q_partials. In DAG case, will want to ensure that
        // left children partials are in p_partials_left
        psv_handler_.GetPV(PSVType::PLeft, node->Id()).col(pattern_idx) =
            ParentPartial(TotalPPartial(node->Children()[0]->Id(), pattern_idx));

        psv_handler_.GetPV(PSVType::PRight, node->Id()).col(pattern_idx) =
            ParentPartial(TotalPPartial(node->Children()[1]->Id(), pattern_idx));
      }
    }
  });

  // generating q-partials
  topology->Preorder([this](const Node* node) {
    if (!node->IsLeaf()) {
      Assert(node->Children().size() == 2,
             "Error in SankoffHandler::RunSankoff: Tree should be bifurcating.");

      for (size_t pattern_idx = 0; pattern_idx < site_pattern_.PatternCount();
           pattern_idx++) {
        auto partials_from_parent =
            ParentPartial(psv_handler_.GetPV(PSVType::Q, node->Id()).col(pattern_idx));
        for (size_t child = 0; child < 2; child++) {
          size_t sister_idx = child == 0 ? 1 : 0;
          auto partials_from_sister = ParentPartial(
              TotalPPartial(node->Children()[sister_idx]->Id(), pattern_idx));
          psv_handler_.GetPV(PSVType::Q, node->Children()[child]->Id())
              .col(pattern_idx) = partials_from_sister + partials_from_parent;
        }
      }
    }
  });
}

double SankoffHandler::ParsimonyScore(size_t node_id) {
  auto weights = site_pattern_.GetWeights();
  double total_parsimony = 0.;
  for (size_t pattern = 0; pattern < site_pattern_.PatternCount(); pattern++) {
    // Note: doing ParentPartial first for the left and right p_partials and then adding
    // them together will give the same minimum parsimony score, but doesn't give
    // correct Sankoff Partial vector for the new rooting
    auto total_tree = ParentPartial(TotalPPartial(node_id, pattern));

    total_tree += ParentPartial(psv_handler_.GetPV(PSVType::Q, node_id).col(pattern));

    // If node_id is the root node, calculating the total_tree vector like so does not
    // yield the SankoffPartial of an actual rooting, but this will not change the
    // minimum value in the partial, so the root node can still be used to calculate the
    // parsimony score.
    total_parsimony +=
        *std::min_element(total_tree.begin(), total_tree.end()) * weights[pattern];
  }
  return total_parsimony;
}
