#pragma once

#include "../src/rooted_tree.hpp"

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
