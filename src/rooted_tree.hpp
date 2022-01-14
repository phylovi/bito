// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// This is a rooted tree class that had the extra parameters required to do node height
// gradients. In fact, because RootedTree also has branch lengths (inherited from Tree)
// all of the extra state in this object is redundant other than the tip dates.
//
// Rooted trees can exist in 3 states:
// 1. No dates associated
// 2. Dates associated, which also initializes node_bounds_
// 3. As an initialized time tree, which means that all members are initialized.
//
// State 3 means that the branch lengths must be compatible with tip dates.
//
// In the terminology here, imagine that the tree is rooted at the top and hangs down.
// The "height" of a node is how far we have to go back in time to that divergence event
// from the present.
//
// The most important parameterization here is in terms of node height ratios, which are
// of the form n/d, where
//
// n = time difference between this node's height and that of its earliest descendant E
// d = time difference between the parent's height and that of E.

#pragma once

#include "eigen_sugar.hpp"
#include "tree.hpp"

class RootedTree : public Tree {
 public:
  using RootedTreeVector = std::vector<RootedTree>;

  RootedTree(const Node::NodePtr& topology, BranchLengthVector branch_lengths);
  explicit RootedTree(const Tree& tree);

  const std::vector<double>& GetNodeBounds() const {
    EnsureTipDatesHaveBeenSet();
    return node_bounds_;
  }
  const std::vector<double>& GetHeightRatios() const {
    EnsureTimeTreeHasBeenInitialized();
    return height_ratios_;
  }
  const std::vector<double>& GetNodeHeights() const {
    EnsureTimeTreeHasBeenInitialized();
    return node_heights_;
  }
  const std::vector<double>& GetRates() const {
    EnsureTimeTreeHasBeenInitialized();
    return rates_;
  }
  size_t RateCount() const { return rate_count_; }

  inline bool TipDatesHaveBeenSet() const { return !node_bounds_.empty(); }
  inline void EnsureTipDatesHaveBeenSet() const {
    if (!TipDatesHaveBeenSet()) {
      Failwith(
          "Attempted access of a time tree member that requires the tip dates to be "
          "set. Have you set dates for your time trees?");
    }
  }

  // Set the tip dates in the tree, and also set the node bounds.  Note that these dates
  // are the amount of time elapsed between the sampling date and the present. Thus,
  // older times have larger dates. This function requires the supplied branch lengths
  // to be clocklike.
  void SetTipDates(const TagDoubleMap& tag_date_map);

  inline bool TimeTreeHasBeenInitialized() const { return !height_ratios_.empty(); }
  inline void EnsureTimeTreeHasBeenInitialized() const {
    if (!TimeTreeHasBeenInitialized()) {
      Failwith(
          "Attempted access of a time tree member that requires the time tree to be "
          "initialized. Have you set dates for your time trees, and initialized the "
          "time trees?");
    }
  }

  // Use the branch lengths to set node heights and height ratios.
  void InitializeTimeTreeUsingBranchLengths();
  // Set node_heights_ so that they have the given height ratios.
  // This should become SetNodeHeights during #205, and actually set height ratios as
  // well.
  void InitializeTimeTreeUsingHeightRatios(EigenConstVectorXdRef height_ratios);

  TagDoubleMap TagDateMapOfDateVector(std::vector<double> leaf_date_vector);

  // The lower bound for the height of each node, which is the maximum of the tip dates
  // across all of the descendants of the node. See top of this file to read about how
  // this vector can be initialized even if the rest of the fields below are not.
  std::vector<double> node_bounds_;
  // This vector is of length equal to the number of internal nodes, and (except for the
  // last entry) has the node height ratios. The last entry is the root height.
  // The indexing is set up so that the ith entry has the node height ratio for the
  // (i+leaf_count)th node for all i except for the last.
  std::vector<double> height_ratios_;
  // The actual node heights for all nodes.
  std::vector<double> node_heights_;
  // The per-branch substitution rates.
  std::vector<double> rates_;
  // Number of substitution rates (e.g. 1 rate for strict clock)
  size_t rate_count_ = 0;

  bool operator==(const Tree& other) const = delete;
  bool operator==(const RootedTree& other) const;

  // The tree `(0:2,(1:1.5,(2:2,3:1):2.5):2.5):0;` as depicted in
  // https://github.com/phylovi/bito/issues/187#issuecomment-618421183
  static RootedTree Example();
  static RootedTree UnitBranchLengthTreeOf(Node::NodePtr topology);

 private:
  // As for SetTipDates, but only set the node bounds. No constraint on supplied
  // branch lengths.
  void SetNodeBoundsUsingDates(const TagDoubleMap& tag_date_map);
};

inline bool operator!=(const RootedTree& lhs, const RootedTree& rhs) {
  return !(lhs == rhs);
}

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("RootedTree") {
  // To understand this test, please see
  // https://github.com/phylovi/bito/issues/187#issuecomment-618421183
  auto tree = RootedTree::Example();
  std::vector<double> correct_height_ratios({1. / 3.5, 1.5 / 4., 7.});
  for (size_t i = 0; i < correct_height_ratios.size(); ++i) {
    CHECK_EQ(correct_height_ratios[i], tree.height_ratios_[i]);
  }
  std::vector<double> correct_node_heights({5., 3., 0., 1., 2., 4.5, 7.});
  std::vector<double> correct_node_bounds({5., 3., 0., 1., 1., 3., 5.});
  std::vector<double> correct_branch_lengths({2., 1.5, 2., 1., 2.5, 2.5});
  for (size_t i = 0; i < correct_node_heights.size(); ++i) {
    CHECK_EQ(correct_node_heights[i], tree.node_heights_[i]);
    CHECK_EQ(correct_node_bounds[i], tree.node_bounds_[i]);
  }
  for (size_t i = 0; i < correct_branch_lengths.size(); ++i) {
    CHECK_EQ(correct_branch_lengths[i], tree.branch_lengths_[i]);
  }
  // Test ratios to heights.
  const double arbitrary_dummy_number = -5.;
  std::fill(tree.LeafCount() + tree.node_heights_.begin(),  // First internal node.
            tree.node_heights_.end(), arbitrary_dummy_number);
  EigenVectorXd new_height_ratios(3);
  // Issue #205: eliminate this code duplication.
  // Root height is multiplied by 2.
  new_height_ratios << 1. / 3.5, 1.5 / 4., 14.;
  std::vector<double> new_correct_node_heights({5., 3., 0., 1., 2.75, 7.125, 14.});
  std::vector<double> new_correct_branch_lengths({9., 4.125, 2.75, 1.75, 4.375, 6.875});
  tree.InitializeTimeTreeUsingHeightRatios(new_height_ratios);
  for (size_t i = 0; i < correct_node_heights.size(); ++i) {
    CHECK_EQ(new_correct_node_heights[i], tree.node_heights_[i]);
  }
  for (size_t i = 0; i < correct_branch_lengths.size(); ++i) {
    CHECK_EQ(new_correct_branch_lengths[i], tree.branch_lengths_[i]);
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED
