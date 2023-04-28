#pragma once

#include "../src/reindexer.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Reindexer: IdentityReindexer") {
  // Check that IdentityReindexer returns correctly.
  Reindexer correct_default({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
  CHECK_EQ(correct_default, Reindexer::IdentityReindexer(10));
}

TEST_CASE("IsValidReindexer") {
  // Index appears more than once.
  CHECK_FALSE(Reindexer({1, 3, 0, 0}).IsValid());
  // Missing an index and/or index is out of range.
  CHECK_FALSE(Reindexer({1, 3, 4, 2}).IsValid());
  // Valid reindexer.
  CHECK(Reindexer({1, 3, 0, 2}).IsValid());
}

TEST_CASE("Reindexer: Reindex") {
  // Check that Reindex throws if the given vector and reindexer have different
  // sizes.
  SizeVector old_size_vector{7, 8, 9};
  Reindexer reindexer({2, 0, 3, 1});
  CHECK_THROWS(Reindexer::Reindex(old_size_vector, reindexer));
  // Check that Reindex returns correctly.
  reindexer = Reindexer({2, 0, 1});
  SizeVector new_size_vector = Reindexer::Reindex(old_size_vector, reindexer);
  SizeVector correct_new_size_vector{8, 9, 7};
  CHECK_EQ(new_size_vector, correct_new_size_vector);
  // Check that Reindex also works with EigenVectorXd and additional values.
  EigenVectorXd old_eigen_vector(3);
  old_eigen_vector << 7, 8, 9;
  EigenVectorXd additional_values(2);
  additional_values << 10, 11;
  reindexer = Reindexer({2, 4, 0, 3, 1});
  EigenVectorXd new_eigen_vector =
      Reindexer::Reindex(old_eigen_vector, reindexer, additional_values);
  EigenVectorXd correct_new_eigen_vector(5);
  correct_new_eigen_vector << 9, 11, 7, 10, 8;
  CHECK_EQ(new_eigen_vector, correct_new_eigen_vector);
}

TEST_CASE("Reindexer: InvertReindexer") {
  // Check that inverting a vector twice results in the original vector.
  Reindexer reindexer({1, 3, 0, 2});
  Reindexer correct_inverted_reindexer({2, 0, 3, 1});
  Reindexer inverted_reindexer = reindexer.InvertReindexer();
  CHECK_EQ(inverted_reindexer, correct_inverted_reindexer);
  Reindexer correct_reindexer = Reindexer({1, 3, 0, 2});
  reindexer = inverted_reindexer.InvertReindexer();
  CHECK_EQ(reindexer, correct_reindexer);
}

TEST_CASE("Reindexer: RemapIdVector") {
  // Check that Reindex throws if the given vector has an index out of bounds of the
  // reindexer.
  SizeVector size_vector{3, 5};
  Reindexer reindexer({2, 0, 3, 1});
  CHECK_THROWS(Reindexer::RemapIdVector<size_t>(size_vector, reindexer));
  // Check that RemapIdVector returns correctly.
  size_vector = {3, 5};
  reindexer = Reindexer({2, 0, 3, 1, 6, 4, 5});
  Reindexer::RemapIdVector<size_t>(size_vector, reindexer);
  SizeVector correct_size_vector = {1, 4};
  CHECK_EQ(size_vector, correct_size_vector);
}

TEST_CASE("Reindexer: ReassignAndShift") {
  // Check that ReassignAndShift returns correctly when old_id > new_id.
  Reindexer reindexer({0, 1, 2, 3, 4, 5, 6});
  reindexer.ReassignAndShift(4, 1);
  Reindexer correct_reindexer({0, 2, 3, 4, 1, 5, 6});
  CHECK_EQ(reindexer, correct_reindexer);
  reindexer.ReassignAndShift(5, 2);
  correct_reindexer = Reindexer({0, 3, 4, 5, 1, 2, 6});
  CHECK_EQ(reindexer, correct_reindexer);
  reindexer.ReassignAndShift(1, 3);
  correct_reindexer = Reindexer({0, 2, 4, 5, 3, 1, 6});
  CHECK_EQ(reindexer, correct_reindexer);
  // Check that ReassignAndShift returns correctly when old_id = new_id.
  reindexer = Reindexer({1, 0, 4, 6, 5, 3, 2});
  reindexer.ReassignAndShift(4, 4);
  correct_reindexer = Reindexer({1, 0, 4, 6, 5, 3, 2});
  CHECK_EQ(reindexer, correct_reindexer);
  // Check that ReassignAndShift returns correctly when old_id < new_id.
  reindexer = Reindexer({6, 0, 4, 1, 5, 3, 2});
  reindexer.ReassignAndShift(1, 5);
  correct_reindexer = Reindexer({6, 0, 3, 5, 4, 2, 1});
  CHECK_EQ(reindexer, correct_reindexer);
}

TEST_CASE("Reindexer: ComposeWith") {
  // Check that identity reindexer composed with a second reindexer results in that
  // reindexer.
  Reindexer identity_reindexer, inverted_reindexer, pairswap_reindexer,
      composed_reindexer, correct_reindexer;
  identity_reindexer = Reindexer::IdentityReindexer(6);
  inverted_reindexer = Reindexer({5, 4, 3, 2, 1, 0});
  pairswap_reindexer = Reindexer({1, 0, 3, 2, 5, 4});
  composed_reindexer = identity_reindexer;
  composed_reindexer = composed_reindexer.ComposeWith(inverted_reindexer);
  composed_reindexer = composed_reindexer.ComposeWith(pairswap_reindexer);
  correct_reindexer = Reindexer({4, 5, 2, 3, 0, 1});
  CHECK_EQ(composed_reindexer, correct_reindexer);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
