// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "subsplit_dag_nni.hpp"

#include "bitset.hpp"
#include "subsplit_dag.hpp"

// ** NNIOperation Methods

int NNIOperation::Compare(const NNIOperation &nni_a, const NNIOperation &nni_b) {
  auto compare_parent = Bitset::Compare(nni_a.parent_, nni_b.parent_);
  if (compare_parent != 0) {
    return compare_parent;
  }
  auto compare_child = Bitset::Compare(nni_a.child_, nni_b.child_);
  return compare_child;
};

int NNIOperation::Compare(const NNIOperation &nni_b) const {
  const NNIOperation &nni_a = *this;
  return Compare(nni_a, nni_b);
}

bool operator<(const NNIOperation &lhs, const NNIOperation &rhs) {
  return NNIOperation::Compare(lhs, rhs) < 0;
}
bool operator>(const NNIOperation &lhs, const NNIOperation &rhs) {
  return NNIOperation::Compare(lhs, rhs) > 0;
}
bool operator==(const NNIOperation &lhs, const NNIOperation &rhs) {
  return NNIOperation::Compare(lhs, rhs) == 0;
}
bool operator!=(const NNIOperation &lhs, const NNIOperation &rhs) {
  return NNIOperation::Compare(lhs, rhs) != 0;
}

NNIOperation NNIOperation::NNIOperationFromNeighboringSubsplits(
    const Bitset parent_in, const Bitset child_in,
    const bool swap_which_child_clade_with_sister, const bool which_child_of_parent) {
  // Input: Parent(X,YZ) -> Child(Y,Z).
  Bitset X = parent_in.SubsplitGetClade(!which_child_of_parent);
  // "Y" clade can be chosen arbitrarily from (Y,Z), so "Y" is chosen based on which
  // we want to swap with "X".
  Bitset Y = child_in.SubsplitGetClade(swap_which_child_clade_with_sister);
  Bitset Z = child_in.SubsplitGetClade(!swap_which_child_clade_with_sister);
  // Output: Parent(Y,XZ) -> Child(X,Z).
  Bitset parent_out = Bitset::Subsplit(Y, X | Z);
  Bitset child_out = Bitset::Subsplit(X, Z);
  return NNIOperation(parent_out, child_out);
}

NNIOperation NNIOperation::NNIOperationFromNeighboringSubsplits(
    const Bitset parent_in, const Bitset child_in,
    const bool swap_which_child_clade_with_sister) {
  bool which_clade_of_parent = Bitset::SubsplitIsChildOfWhichParentClade(parent_in, child_in);
  return NNIOperationFromNeighboringSubsplits(
      parent_in, child_in, swap_which_child_clade_with_sister, which_clade_of_parent);
}

// ** SetOfNNIs Methods

bool operator==(const SetOfNNIs &lhs, const SetOfNNIs &rhs) {
  return lhs.set_ == rhs.set_;
}

bool operator!=(const SetOfNNIs &lhs, const SetOfNNIs &rhs) {
  return lhs.set_ != rhs.set_;
}

void SetOfNNIs::Insert(NNIOperation nni_op) { set_.insert(nni_op); };
void SetOfNNIs::Insert(Bitset parent, Bitset child) {
  Insert(NNIOperation(parent, child));
}

void SetOfNNIs::Erase(NNIOperation nni_op) { set_.erase(nni_op); };
void SetOfNNIs::Erase(Bitset parent, Bitset child) {
  Erase(NNIOperation(parent, child));
}

void SetOfNNIs::Clear() { set_.clear(); }

size_t SetOfNNIs::GetSize() const { return set_.size(); }

// ** NNIEvaluationEngine Methods

void SyncSetOfNNIsWithDAG(SetOfNNIs &set_of_nnis, const SubsplitDAG &dag) {
  set_of_nnis.Clear();
  // Only real node pairs are viable NNIs.
  dag.IterateOverRealNodes([&set_of_nnis, &dag](const SubsplitDAGNode *node) {
    dag.IterateOverParentAndChildAndLeafwardEdges(
        node, [&set_of_nnis, &dag](const size_t parent_id, const bool is_rotated,
                                   const size_t child_id, const size_t edge_idx) {
          // Only internal node pairs are viable NNIs.
          Bitset parent_bitset = dag.GetDAGNode(parent_id)->GetBitset();
          Bitset child_bitset = dag.GetDAGNode(child_id)->GetBitset();
          if (!(parent_bitset.SubsplitIsRoot() || child_bitset.SubsplitIsLeaf())) {
            SafeAddOutputNNIsToSetOfNNIs(set_of_nnis, dag, parent_bitset, child_bitset,
                                         is_rotated);
          }
        });
  });
}

void UpdateSetOfNNIsAfterDAGAddNodePair(SetOfNNIs &set_of_nnis, const SubsplitDAG &dag,
                                        const Bitset &parent_bitset,
                                        const Bitset &child_bitset) {
  size_t parent_id = dag.GetDAGNodeId(parent_bitset);
  size_t child_id = dag.GetDAGNodeId(child_bitset);
  // Every new edge added is a potential new NNI.
  // Iterate over the parent and child node of the new pair.
  for (const size_t &node_id : {parent_id, child_id}) {
    // Get nodes adjacent to current node from both sorted and rotated edges.
    for (const bool is_edge_leafward : {true, false}) {
      // Get nodes adjacent to current node from both leafward and rootward directions.
      for (const bool is_edge_rotated : {true, false}) {
        SizeVector adjacent_node_ids = dag.GetDAGNode(node_id)->GetLeafwardOrRootward(
            is_edge_leafward, is_edge_rotated);
        AddAllNNIsFromNodeVectorToSetOfNNIs(set_of_nnis, dag, node_id,
                                            adjacent_node_ids, is_edge_rotated,
                                            is_edge_leafward);
      }
    }
  }
  // Remove the pair that was just added to the DAG from NNI Set.
  NNIOperation new_nni = NNIOperation(parent_bitset, child_bitset);
  set_of_nnis.Erase(new_nni);
}

void AddAllNNIsFromNodeVectorToSetOfNNIs(SetOfNNIs &set_of_nnis, const SubsplitDAG &dag,
                                         const size_t &node_id,
                                         const SizeVector &adjacent_node_ids,
                                         const bool is_edge_rotated,
                                         const bool is_edge_leafward) {
  Bitset node_bitset = dag.GetDAGNode(node_id)->GetBitset();
  // Determine whether node_id corresponds to parent or child of the pair.
  // Add every edge's NNI to NNI Set.
  // If edges are leafward, node_id is the parent to all vector nodes.
  // If edges are rootward, node_id is the child to all vector nodes.
  if (is_edge_leafward) {
    const Bitset &parent_bitset = node_bitset;
    for (const auto &adjacent_node_id : adjacent_node_ids) {
      const Bitset child_bitset = dag.GetDAGNode(adjacent_node_id)->GetBitset();
      SafeAddOutputNNIsToSetOfNNIs(set_of_nnis, dag, parent_bitset, child_bitset,
                                   is_edge_rotated);
    }
  } else {
    const Bitset &child_bitset = node_bitset;
    for (const auto &adjacent_node_id : adjacent_node_ids) {
      const Bitset parent_bitset = dag.GetDAGNode(adjacent_node_id)->GetBitset();
      SafeAddOutputNNIsToSetOfNNIs(set_of_nnis, dag, parent_bitset, child_bitset,
                                   is_edge_rotated);
    }
  }
}

void SafeAddOutputNNIsToSetOfNNIs(SetOfNNIs &set_of_nnis, const SubsplitDAG &dag,
                                  const Bitset &parent_bitset,
                                  const Bitset &child_bitset,
                                  const bool is_edge_rotated) {
  // Soft assert that parent is not the root and child is not a leaf.
  if (parent_bitset.SubsplitIsRoot() || child_bitset.SubsplitIsLeaf()) {
    return;
  }
  // Input pair is in the DAG, so remove it from the Set if it exists.
  set_of_nnis.Erase(parent_bitset, child_bitset);
  // Add NNI for sorted clade swap and rotated clade swap.
  for (bool is_swap_with_sorted_child : {true, false}) {
    bool is_in_dag = false;
    const auto new_nni = NNIOperation::NNIOperationFromNeighboringSubsplits(
        parent_bitset, child_bitset, is_swap_with_sorted_child, !is_edge_rotated);
    if (dag.ContainsNode(new_nni.parent_) && dag.ContainsNode(new_nni.child_)) {
      const size_t parent_id = dag.GetDAGNodeId(new_nni.parent_);
      const size_t child_id = dag.GetDAGNodeId(new_nni.child_);
      is_in_dag = dag.ContainsEdge(parent_id, child_id);
    }
    if (is_in_dag == false) {
      set_of_nnis.Insert(new_nni);
    }
  }
}
