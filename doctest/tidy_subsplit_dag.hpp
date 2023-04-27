#pragma once

#include "../src/tidy_subsplit_dag.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TidySubsplitDAG: slicing") {
  auto manual_dag = TidySubsplitDAG::ManualTrivialExample();

  // std::cout << manual_dag.AboveMatricesAsString() << std::endl;
  CHECK_EQ(GenericToString(manual_dag.AboveNode(NodeId(0))), "[1, 0, 0, 1, 1, 1]\n");
  CHECK_EQ(GenericToString(manual_dag.AboveNode(NodeId(1))), "[0, 1, 0, 1, 1, 1]\n");
  CHECK_EQ(GenericToString(manual_dag.AboveNode(NodeId(2))), "[0, 0, 1, 0, 1, 1]\n");
  CHECK_EQ(GenericToString(manual_dag.AboveNode(NodeId(3))), "[0, 0, 0, 1, 1, 1]\n");
  CHECK_EQ(GenericToString(manual_dag.AboveNode(NodeId(4))), "[0, 0, 0, 0, 1, 1]\n");
  CHECK_EQ(GenericToString(manual_dag.AboveNode(NodeId(5))), "[0, 0, 0, 0, 0, 1]\n");

  auto trivial_dag = TidySubsplitDAG::TrivialExample();
  CHECK_EQ(trivial_dag.AboveMatricesAsString(), manual_dag.AboveMatricesAsString());

  auto motivating_dag = TidySubsplitDAG::MotivatingExample();
  CHECK_EQ(GenericToString(motivating_dag.AboveNode(false, NodeId(4))),
           "[0, 0, 0, 0, 1, 1, 1, 1, 0, 0]\n");
  CHECK_EQ(GenericToString(motivating_dag.AboveNode(true, NodeId(4))),
           "[0, 0, 0, 0, 1, 0, 0, 0, 1, 1]\n");
  CHECK_EQ(GenericToString(motivating_dag.AboveNode(false, NodeId(7))),
           "[0, 0, 0, 0, 0, 0, 0, 1, 0, 0]\n");
  CHECK_EQ(GenericToString(motivating_dag.AboveNode(true, NodeId(7))),
           "[0, 0, 0, 0, 0, 0, 0, 1, 1, 1]\n");
  CHECK_EQ(GenericToString(motivating_dag.BelowNode(false, NodeId(7))),
           "[0, 0, 1, 1, 1, 0, 0, 1, 0, 0]\n");
  CHECK_EQ(GenericToString(motivating_dag.BelowNode(true, NodeId(7))),
           "[1, 0, 0, 0, 0, 0, 0, 1, 0, 0]\n");

  motivating_dag.SetDirtyStrictlyAbove(NodeId(4));
  CHECK_EQ(GenericToString(motivating_dag.DirtyVector(true)),
           "[0, 0, 0, 0, 0, 0, 0, 0, 1, 1]\n");
  CHECK_EQ(GenericToString(motivating_dag.DirtyVector(false)),
           "[0, 0, 0, 0, 0, 1, 1, 1, 0, 0]\n");

  motivating_dag.SetClean();
  // #321 Add test for Tidy traversal.
}
#endif  // DOCTEST_LIBRARY_INCLUDED
