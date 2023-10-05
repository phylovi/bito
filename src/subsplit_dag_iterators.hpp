// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//

class SubsplitDAG::LeafIterator {
 public:
  LeafIterator(const SubsplitDAG& dag, NodeId node_id = 0) : dag_(dag) {
    if (node_id > dag_.TaxonCount()) {
      Failwith("NodeId is out-of-range.");
    }
    node_id_ = node_id;
  }

  LeafIterator& operator++() {
    node_id_ = NodeId{node_id_.value_ + 1};
    return *this;
  }

  bool operator==(const NodeId node_id) const { return node_id_ == node_id; }

  bool operator!=(const NodeId node_id) const { return node_id_ != node_id; }

  SubsplitDAGNode operator*() { return dag_.GetDAGNode(node_id_); }

  NodeId begin() { return NodeId{0}; }

  NodeId end() { return NodeId{dag_.TaxonCount()}; }

 private:
  const SubsplitDAG& dag_;
  NodeId node_id_;
};
