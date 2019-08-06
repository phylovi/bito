// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "node.hpp"
#include <limits.h>
#include <algorithm>
#include <cassert>
#include <deque>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

Node::Node(uint32_t leaf_id)
    : children_({}),
      index_(leaf_id),
      tag_(PackInts(leaf_id, 1)),
      hash_(SOHash(leaf_id)) {}

Node::Node(NodePtrVec children, size_t index)
    : children_(children), index_(index) {
  Assert(!children_.empty(),
         "Called internal Node constructor with no children.");
  // Order the children by their max leaf ids.
  std::sort(
      children_.begin(), children_.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs->MaxLeafID() == rhs->MaxLeafID()) {
          // Children should have non-overlapping leaf sets, so there
          // should not be ties.
          Failwith("Tie observed between " + lhs->Newick() + " and " +
                   rhs->Newick() + "\n" + "Do you have a taxon name repeated?");
        }
        return lhs->MaxLeafID() < rhs->MaxLeafID();
      });
  // Children are sorted by their max_leaf_id, so we can get the max by
  // looking at the last element.
  uint32_t max_leaf_id = children_.back()->MaxLeafID();
  uint32_t leaf_count = 0;
  hash_ = 0;
  for (const auto& child : children_) {
    leaf_count += child->LeafCount();
    hash_ ^= child->Hash();
  }
  tag_ = PackInts(max_leaf_id, leaf_count);
  // Bit rotation is necessary because if we only XOR then we can get
  // collisions when identical tips are in different
  // ordered subtrees (an example is in below doctest).
  hash_ = SORotate(hash_, 1);
}

bool Node::operator==(const Node& other) {
  if (this->Hash() != other.Hash()) {
    return false;
  }
  size_t child_count = this->Children().size();
  if (child_count != other.Children().size()) {
    return false;
  }
  for (size_t i = 0; i < child_count; i++) {
    if (!(*children_[i] == *other.Children()[i])) {
      return false;
    }
  }
  return true;
}

void Node::PreOrder(std::function<void(const Node*)> f) const {
  f(this);
  for (const auto& child : children_) {
    child->PreOrder(f);
  }
}

void Node::ConditionalPreOrder(std::function<bool(const Node*)> f) const {
  if (f(this)) {
    for (const auto& child : children_) {
      child->ConditionalPreOrder(f);
    }
  }
}

void Node::PostOrder(std::function<void(const Node*)> f) const {
  for (const auto& child : children_) {
    child->PostOrder(f);
  }
  f(this);
}

void Node::LevelOrder(std::function<void(const Node*)> f) const {
  std::deque<const Node*> to_visit = {this};
  while (to_visit.size()) {
    auto n = to_visit.front();
    f(n);
    to_visit.pop_front();

    for (const auto& child : n->children_) {
      to_visit.push_back(child.get());
    }
  }
}

void Node::TriplePreOrderBifurcating(
    std::function<void(const Node*, const Node*, const Node*)> f) const {
  if (!IsLeaf()) {
    Assert(children_.size() == 2,
           "TripleIndexPreOrderBifurcating expects a bifurcating tree.");
    f(this, children_[1].get(), children_[0].get());
    children_[0]->TriplePreOrderBifurcating(f);
    f(this, children_[0].get(), children_[1].get());
    children_[1]->TriplePreOrderBifurcating(f);
  }
}

static std::function<void(const Node*, const Node*, const Node*)> const
TripletIndexInfix(std::function<void(int, int, int)> f) {
  return [&f](const Node* node0, const Node* node1, const Node* node2) {
    f(static_cast<int>(node0->Index()), static_cast<int>(node1->Index()),
      static_cast<int>(node2->Index()));
  };
}

void Node::TripleIndexPreOrderBifurcating(
    std::function<void(int, int, int)> f) const {
  TriplePreOrderBifurcating(TripletIndexInfix(f));
}

static std::function<void(const Node*)> const BinaryIndexInfix(
    std::function<void(int, int, int)> f) {
  return [&f](const Node* node) {
    if (!node->IsLeaf()) {
      Assert(node->Children().size() == 2,
             "BinaryIndexInfix expects a bifurcating tree.");
      f(static_cast<int>(node->Index()),
        static_cast<int>(node->Children()[0]->Index()),
        static_cast<int>(node->Children()[1]->Index()));
    }
  };
}

void Node::BinaryIndexPreOrder(
    const std::function<void(int, int, int)> f) const {
  PreOrder(BinaryIndexInfix(f));
}

void Node::BinaryIndexPostOrder(
    const std::function<void(int, int, int)> f) const {
  PostOrder(BinaryIndexInfix(f));
}

void Node::TriplePreOrder(
    std::function<void(const Node*, const Node*, const Node*)> f_root,
    std::function<void(const Node*, const Node*, const Node*)> f_internal)
    const {
  Assert(children_.size() == 3,
         "TriplePreOrder expects a tree with a trifurcation at the root.");
  f_root(children_[0].get(), children_[1].get(), children_[2].get());
  f_root(children_[1].get(), children_[2].get(), children_[0].get());
  f_root(children_[2].get(), children_[0].get(), children_[1].get());
  for (const auto& child : children_) {
    child->TriplePreOrderInternal(f_internal);
  }
}

void Node::TriplePreOrderInternal(
    std::function<void(const Node*, const Node*, const Node*)> f) const {
  if (!IsLeaf()) {
    Assert(children_.size() == 2,
           "TriplePreOrderInternal expects a bifurcating tree.");
    f(this, children_[1].get(), children_[0].get());
    children_[0]->TriplePreOrderInternal(f);
    f(this, children_[0].get(), children_[1].get());
    children_[1]->TriplePreOrderInternal(f);
  }
}

// See the typedef of PCSSFun to understand the argument type to this
// function, and `doc/pcss.svg` for a diagram that will greatly help you
// understand the implementation.
void Node::PCSSPreOrder(PCSSFun f) const {
  this->TriplePreOrder(
      // f_root
      [&f](const Node* node0, const Node* node1, const Node* node2) {
        // Virtual root on node2's edge, with subsplit pointing up.
        f(node2, false, node2, true, node0, false, node1, false, nullptr);
        if (!node2->IsLeaf()) {
          Assert(node2->Children().size() == 2,
                 "PCSSPreOrder expects a bifurcating tree.");
          auto child0 = node2->Children()[0].get();
          auto child1 = node2->Children()[1].get();
          // Virtual root in node1.
          f(node0, false, node2, false, child0, false, child1, false, node1);
          // Virtual root in node0.
          f(node1, false, node2, false, child0, false, child1, false, node0);
          // Virtual root on node2's edge, with subsplit pointing down.
          f(node2, true, node2, false, child0, false, child1, false, nullptr);
          // Virtual root in child0.
          f(child1, false, node2, true, node0, false, node1, false, child0);
          // Virtual root in child1.
          f(child0, false, node2, true, node0, false, node1, false, child1);
        }
      },
      // f_internal
      [&f, this](const Node* parent, const Node* sister, const Node* node) {
        // Virtual root on node's edge, with subsplit pointing up.
        f(node, false, node, true, parent, true, sister, false, nullptr);
        if (!node->IsLeaf()) {
          Assert(node->Children().size() == 2,
                 "PCSSPreOrder expects a bifurcating tree.");
          auto child0 = node->Children()[0].get();
          auto child1 = node->Children()[1].get();
          // Virtual root up the tree.
          f(sister, false, node, false, child0, false, child1, false, this);
          // Virtual root in sister.
          f(parent, true, node, false, child0, false, child1, false, sister);
          // Virtual root on node's edge, with subsplit pointing down.
          f(node, true, node, false, child0, false, child1, false, nullptr);
          // Virtual root in child0.
          f(child1, false, node, true, sister, false, parent, true, child0);
          // Virtual root in child1.
          f(child0, false, node, true, sister, false, parent, true, child1);
        }
      });
}

// This function assigns indices to the nodes of the topology: the leaves get
// their indices (which are contiguously numbered from 0 through the leaf
// count -1) and the rest get ordered according to a postorder traversal. Thus
// the root always has index equal to the number of nodes in the tree.
//
// This function returns a map that maps the tags to their indices.
TagSizeMap Node::Reindex() {
  TagSizeMap tag_index_map;
  size_t next_index = 1 + MaxLeafID();
  MutablePostOrder([&tag_index_map, &next_index](Node* node) {
    if (node->IsLeaf()) {
      node->index_ = node->MaxLeafID();
    } else {
      node->index_ = next_index;
      next_index++;
    }
    SafeInsert(tag_index_map, node->Tag(), node->index_);
  });
  return tag_index_map;
}

std::string Node::Newick(std::function<std::string(const Node*)> node_labeler,
                         const DoubleVectorOption& branch_lengths) const {
  return NewickAux(node_labeler, branch_lengths) + ";";
}

std::string Node::NewickAux(
    std::function<std::string(const Node*)> node_labeler,
    const DoubleVectorOption& branch_lengths) const {
  std::string str;
  if (IsLeaf()) {
    str.assign(node_labeler(this));
  } else {
    str.assign("(");
    for (auto iter = children_.begin(); iter != children_.end(); iter++) {
      if (iter != children_.begin()) {
        str.append(",");
      }
      str.append((*iter)->NewickAux(node_labeler, branch_lengths));
    }
    str.append(")");
    str.append(node_labeler(this));
  }
  if (branch_lengths) {
    Assert(Index() < (*branch_lengths).size(),
           "branch_lengths vector is of insufficient length in NewickAux.");
    // ostringstream is the way to get scientific notation using the STL.
    std::ostringstream str_stream;
    str_stream << (*branch_lengths)[Index()];
    str.append(":" + str_stream.str());
  }
  return str;
}

std::string Node::Newick(const DoubleVectorOption& branch_lengths,
                         const TagStringMapOption& node_labels,
                         bool show_tags) const {
  return Newick(
      [&node_labels, &show_tags](const Node* node) {
        if (node->IsLeaf()) {
          if (node_labels) {
            return (*node_labels).at(node->Tag());
          } else if (show_tags) {
            return node->TagString();
          } else {
            return std::to_string(node->MaxLeafID());
          }
        } else {
          if (show_tags) {
            return node->TagString();
          }
        }
        return std::string("");
      },
      branch_lengths);
}

std::vector<size_t> Node::ParentIndexVector() {
  std::vector<size_t> indices(Index());
  PostOrder([&indices](const Node* node) {
    if (!node->IsLeaf()) {
      for (const auto& child : node->Children()) {
        Assert(child->Index() < indices.size(),
               "Problematic indices in ParentIndexVector.");
        indices[child->Index()] = node->Index();
      }
    }
  });
  return indices;
}

Node::NodePtr Node::Deroot() {
  Assert(LeafCount() >= 3, "Node::Deroot expects a tree with at least 3 tips.");
  Assert(Children().size() == 2, "Can't deroot a non-bifurcating tree.");
  auto deroot = [](const NodePtr other_child, const NodePtr has_descendants) {
    // Make a vector copy by passing a vector in.
    NodePtrVec children(has_descendants->Children());
    children.push_back(other_child);
    // has_descendants' index is now available.
    return Join(children, has_descendants->Index());
  };
  if (children_[1]->LeafCount() == 1) {
    return deroot(children_[1], children_[0]);
  }  // else
  return deroot(children_[0], children_[1]);
}

// Class methods
Node::NodePtr Node::Leaf(uint32_t id) { return std::make_shared<Node>(id); }
Node::NodePtr Node::Join(NodePtrVec children, size_t index) {
  return std::make_shared<Node>(children, index);
}
Node::NodePtr Node::Join(NodePtr left, NodePtr right, size_t index) {
  return Join(std::vector<NodePtr>({left, right}), index);
}

Node::NodePtr Node::OfParentIndexVector(std::vector<size_t> indices) {
  // We will fill this map with the indices of the descendants.
  std::unordered_map<size_t, std::vector<size_t>> downward_indices;
  for (size_t child_index = 0; child_index < indices.size(); child_index++) {
    const auto& parent_index = indices[child_index];
    auto search = downward_indices.find(parent_index);
    if (search == downward_indices.end()) {
      // The first time we have seen this parent.
      std::vector<size_t> child_indices({child_index});
      SafeInsert(downward_indices, parent_index, std::move(child_indices));
    } else {
      // We've seen the parent before, so append the child to the parent's
      // vector of descendants.
      search->second.push_back(child_index);
    }
  }
  std::function<NodePtr(size_t)> build_tree =
      [&build_tree, &downward_indices](size_t current_index) {
        auto search = downward_indices.find(current_index);
        if (search == downward_indices.end()) {
          // We assume that anything not in the map is a leaf, because leaves
          // don't have any children.
          return Leaf(static_cast<uint32_t>(current_index));
        } else {
          const auto& children_indices = search->second;
          std::vector<NodePtr> children;
          for (const auto& child_index : children_indices) {
            children.push_back(build_tree(child_index));
          }
          return Join(children, current_index);
        }
      };
  // We assume that the maximum index of the tree is the length of the input
  // index array. That makes sense because the root does not have a parent, so
  // is the first "missing" entry in the input index array.
  return build_tree(indices.size());
}

Node::NodePtrVec Node::ExampleTopologies() {
  NodePtrVec topologies = {
      // 0: (0,1,(2,3))
      Join(std::vector<NodePtr>({Leaf(0), Leaf(1), Join(Leaf(2), Leaf(3))})),
      // 1; (0,1,(2,3)) again
      Join(std::vector<NodePtr>({Leaf(1), Leaf(0), Join(Leaf(3), Leaf(2))})),
      // 2: (0,2,(1,3))
      Join(std::vector<NodePtr>({Leaf(0), Leaf(2), Join(Leaf(1), Leaf(3))})),
      // 3: (0,(1,(2,3)))
      Join(std::vector<NodePtr>(
          {Leaf(0), Join(Leaf(1), Join(Leaf(2), Leaf(3)))}))};
  for (auto& topology : topologies) {
    topology->Reindex();
  }
  return topologies;
}

inline uint32_t Node::SOHash(uint32_t x) {
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = (x >> 16) ^ x;
  return x;
}

inline size_t Node::SORotate(size_t n, uint32_t c) {
  const uint32_t mask =
      (CHAR_BIT * sizeof(n) - 1);  // assumes width is a power of 2.
  // assert ( (c<=mask) &&"rotate by type width or more");
  c &= mask;
  return (n << c) | (n >> ((-c) & mask));
}

void Node::MutablePostOrder(std::function<void(Node*)> f) {
  for (const auto& child : children_) {
    child->MutablePostOrder(f);
  }
  f(this);
}

void Node::PrePostOrder(std::function<void(const Node*)> pre,
                        std::function<void(const Node*)> post) const {
  pre(this);
  for (const auto& child : children_) {
    child->PrePostOrder(pre, post);
  }
  post(this);
}

SizeVectorVector Node::IndicesAbove() {
  SizeVectorVector indices_above(Index() + 1);
  SizeVector mutable_indices;
  PrePostOrder(
      [&indices_above, &mutable_indices](const Node* node) {
        // Store the current set of indices above.
        indices_above[node->Index()] = SizeVector(mutable_indices);
        // As we travel down the tree, the current node will be above.
        mutable_indices.push_back(node->Index());
      },
      // Going back up the tree, so remove the current node's index.
      [&mutable_indices](const Node*) { mutable_indices.pop_back(); });
  return indices_above;
}
