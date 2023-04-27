#pragma once

#include "../src/subsplit_dag_node.hpp"

#ifdef DOCTEST_LIBRARY_INCLUDED

inline DAGVertex& GetStorage(const GenericSubsplitDAGNode<DAGVertex>& node) {
  return node.node_;
}
inline const DAGVertex& GetStorage(
    const GenericSubsplitDAGNode<const DAGVertex>& node) {
  return node.node_;
}

/* Create the following topology:
            [0]
            / \
           0   1
          /     \
        [1]     [2]
        / \
       2   3
      /     \
    [3]     [4]
 */

static SubsplitDAGStorage MakeStorage() {
  SubsplitDAGStorage storage;
  storage.AddLine({EdgeId(0), NodeId(0), NodeId(1), SubsplitClade::Left});
  storage.AddLine({EdgeId(1), NodeId(0), NodeId(2), SubsplitClade::Right});
  storage.AddLine({EdgeId(2), NodeId(1), NodeId(3), SubsplitClade::Left});
  storage.AddLine({EdgeId(3), NodeId(1), NodeId(4), SubsplitClade::Right});

  storage.AddVertex(DAGVertex{}.SetId(NodeId(0)))
      .AddNeighbor(Direction::Leafward, SubsplitClade::Left, NodeId(1), EdgeId(0))
      .AddNeighbor(Direction::Leafward, SubsplitClade::Right, NodeId(2), EdgeId(1));
  storage.AddVertex(DAGVertex{}.SetId(NodeId(1)))
      .AddNeighbor(Direction::Rootward, SubsplitClade::Left, NodeId(0), EdgeId(0))
      .AddNeighbor(Direction::Leafward, SubsplitClade::Left, NodeId(3), EdgeId(2))
      .AddNeighbor(Direction::Leafward, SubsplitClade::Right, NodeId(4), EdgeId(3));
  storage.AddVertex(DAGVertex{}.SetId(NodeId(2)))
      .AddNeighbor(Direction::Rootward, SubsplitClade::Right, NodeId(0), EdgeId(1));
  storage.AddVertex(DAGVertex{}.SetId(NodeId(3)))
      .AddNeighbor(Direction::Rootward, SubsplitClade::Left, NodeId(1), EdgeId(2));
  storage.AddVertex(DAGVertex{}.SetId(NodeId(4)))
      .AddNeighbor(Direction::Rootward, SubsplitClade::Right, NodeId(1), EdgeId(3));
  return storage;
}

TEST_CASE("SubsplitDAGStorage: LinesView structured binding") {
  auto storage = MakeStorage();

  size_t i = 0;
  for (auto [node_ids, line_id] : storage.GetLines()) {
    std::ignore = line_id;
    auto [parent_id, child_id] = node_ids;
    switch (i++) {
      case 0:
        CHECK_EQ(parent_id, 0);
        CHECK_EQ(child_id, 1);
        break;
      case 1:
        CHECK_EQ(parent_id, 0);
        CHECK_EQ(child_id, 2);
        break;
      case 2:
        CHECK_EQ(parent_id, 1);
        CHECK_EQ(child_id, 3);
        break;
      case 3:
        CHECK_EQ(parent_id, 1);
        CHECK_EQ(child_id, 4);
        break;
      default:
        Failwith("More lines than expected");
    }
  }
}

TEST_CASE("SubsplitDAGStorage: Neighbors iterator") {
  auto storage = MakeStorage();

  CHECK_EQ(*GetStorage(storage.GetVertices()[1])
                .GetNeighbors(Direction::Leafward, SubsplitClade::Left)
                .begin(),
           3);
  CHECK_EQ(GetStorage(storage.GetVertices()[1])
               .GetNeighbors(Direction::Leafward, SubsplitClade::Left)
               .begin()
               .GetEdge(),
           2);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
