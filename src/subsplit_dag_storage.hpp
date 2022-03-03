// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

// The following classes are used internally by SubsplitDAG for storing the nodes and
// edges, and providing convenient views and lookups into the data. Terminology has been
// changed in order to distinguish from the public API of SubsplitDAG - Edge becomes
// DAGLine and Node becomes DAGVertex.
//
// Functionality interface of the individual components has been separated from the
// storage containers to allow for flexible representation of multi-element collections:
// class DAGLineStorage - owns the data representing an Edge (IDs of the two nodes, edge
// ID, clade) class DAGLineView - cheaply copyable wrapper for a reference to
// DAGLineStorage class DAGLine - a pure interface providing accessor methods to
// DAGLine{View,Storage} aliases LineView and ConstLineView - used within SubsplitDAG to
// access const and mutable edges class GenericLinesView - cheaply copyable view into a
// collection of edges
//
// Similar design is applied to nodes - they have storage, view, interface and
// collection classes. The node view class is called GenericSubsplitDAGNode and is
// defined in subsplit_dag_node.hpp. Additionally nodes exposes class
// GenericNeighborsView, which is a view that wraps a reference to the node neighbors
// collection - internally a std::map<VertexId, LineId>.
//
// The SubsplitDAGStorage class glues together all components, and is owned exclusively
// by SubsplitDAG.

#pragma once

#include <limits>
#include <vector>
#include <map>
#include <optional>
#include <functional>
#include <utility>
#include <type_traits>
#include <memory>

#include "bitset.hpp"

enum class Direction { Rootward, Leafward };
enum class Clade { Unspecified, Left, Right };

using VertexId = size_t;
using LineId = size_t;
static constexpr size_t NoId = std::numeric_limits<size_t>::max();

template <typename Derived>
class DAGLine {
 public:
  LineId GetId() const { return storage().id_; }
  VertexId GetParent() const { return storage().parent_; }
  VertexId GetChild() const { return storage().child_; }
  Clade GetClade() const { return storage().clade_; }

  Derived& SetId(LineId id) {
    storage().id_ = id;
    return derived();
  }
  Derived& SetParent(VertexId id) {
    storage().parent_ = id;
    return derived();
  }
  Derived& SetChild(VertexId id) {
    storage().child_ = id;
    return derived();
  }
  Derived& SetClade(Clade clade) {
    storage().clade_ = clade;
    return derived();
  }

  std::pair<VertexId, VertexId> GetVertexIds() const {
    return {GetParent(), GetChild()};
  }

 private:
  Derived& derived() { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<const Derived&>(*this); }
  auto& storage() { return derived().storage(); }
  const auto& storage() const { return derived().storage(); }
};

template <typename T>
class DAGLineView;

class DAGLineStorage : public DAGLine<DAGLineStorage> {
 public:
  DAGLineStorage() = default;
  DAGLineStorage(const DAGLineStorage&) = default;
  DAGLineStorage(LineId id, VertexId parent, VertexId child, Clade clade)
      : id_{id}, parent_{parent}, child_{child}, clade_{clade} {}

  template <typename T>
  DAGLineStorage& operator=(DAGLineView<T> other) {
    *this = other.line_;
    return *this;
  }

 private:
  // :: is workaround for GCC bug: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=52625
  template <typename>
  friend class ::DAGLine;
  DAGLineStorage& storage() { return *this; }
  const DAGLineStorage& storage() const { return *this; }

  LineId id_ = NoId;
  VertexId parent_ = NoId;
  VertexId child_ = NoId;
  Clade clade_ = Clade::Unspecified;
};

template <typename T>
class DAGLineView : public DAGLine<DAGLineView<T>> {
 public:
  DAGLineView(T& line) : line_{line} {
    Assert(line.GetClade() != Clade::Unspecified, "Uninitialized edge");
  }

  template <size_t I>
  auto get() const {
    static_assert(I < 2, "Index out of bounds");
    if constexpr (I == 0) return line_.GetVertexIds();
    if constexpr (I == 1) return line_.GetId();
  }

 private:
  template <typename>
  friend class DAGLine;
  friend class DAGLineStorage;
  DAGLineStorage& storage() { return line_; }
  const DAGLineStorage& storage() const { return line_; }

  T& line_;
};
using LineView = DAGLineView<DAGLineStorage>;
using ConstLineView = DAGLineView<const DAGLineStorage>;

namespace std {
template <>
struct tuple_size<::LineView> : public integral_constant<size_t, 2> {};
template <>
struct tuple_element<0, ::LineView> {
  using type = pair<::VertexId, ::VertexId>;
};
template <>
struct tuple_element<1, ::LineView> {
  using type = ::LineId;
};

template <>
struct tuple_size<::ConstLineView> : public integral_constant<size_t, 2> {};
template <>
struct tuple_element<0, ::ConstLineView> {
  using type = const pair<::VertexId, ::VertexId>;
};
template <>
struct tuple_element<1, ::ConstLineView> {
  using type = const ::LineId;
};
}  // namespace std

template <typename T>
class GenericNeighborsView {
  using iterator_type =
      std::conditional_t<std::is_const_v<T>, std::map<VertexId, LineId>::const_iterator,
                         std::map<VertexId, LineId>::iterator>;
  using map_type =
      std::conditional_t<std::is_const_v<T>, const std::map<VertexId, LineId>,
                         std::map<VertexId, LineId>>;

 public:
  class Iterator : public std::iterator<std::forward_iterator_tag, VertexId> {
   public:
    Iterator(const iterator_type& i, map_type& map) : i_{i}, map_{map} {}

    Iterator operator++() {
      ++i_;
      return *this;
    }
    bool operator!=(const Iterator& other) const { return i_ != other.i_; }
    bool operator==(const Iterator& other) const { return i_ == other.i_; }

    VertexId operator*() { return i_->first; }
    LineId GetEdge() const { return i_->second; }

   private:
    iterator_type i_;
    map_type& map_;
  };

  GenericNeighborsView(map_type& neighbors) : neighbors_{neighbors} {}
  template <typename U>
  GenericNeighborsView(const GenericNeighborsView<U>& other)
      : neighbors_{other.neighbors_} {}

  Iterator begin() const { return {neighbors_.begin(), neighbors_}; }
  Iterator end() const { return {neighbors_.end(), neighbors_}; }
  size_t size() const { return neighbors_.size(); }
  bool empty() const { return neighbors_.empty(); }

  void RemapIds(const SizeVector& reindexer) {
    for (auto [vertex_id, line_id] : neighbors_) {
      std::ignore = line_id;
      Assert(vertex_id < reindexer.size(),
            "Neighbors cannot contain an id out of bounds of the reindexer in "
            "GenericNeighborsView::RemapIds.");
    }
    T remaped{};
    for (auto [vertex_id, line_id] : neighbors_) {
      remaped[reindexer.at(vertex_id)] = line_id;
    }
    neighbors_ = remaped;
  }

  operator SizeVector() const {
    SizeVector result;
    for (auto i : neighbors_) result.push_back(i.first);
    return result;
  }

  void SetNeighbors(const T& neighbors) { neighbors_ = neighbors; }

 private:
  friend class Iterator;
  template <typename>
  friend class ::GenericNeighborsView;
  T& neighbors_;
};
using NeighborsView = GenericNeighborsView<std::map<VertexId, LineId>>;
using ConstNeighborsView = GenericNeighborsView<const std::map<VertexId, LineId>>;

class DAGVertex;
template <typename>
class GenericSubsplitDAGNode;
using MutableSubsplitDAGNode = GenericSubsplitDAGNode<DAGVertex>;
using SubsplitDAGNode = GenericSubsplitDAGNode<const DAGVertex>;

class DAGVertex {
 public:
  DAGVertex() = default;
  DAGVertex(const DAGVertex&) = default;
  inline DAGVertex(SubsplitDAGNode node);
  inline DAGVertex(MutableSubsplitDAGNode node);
  DAGVertex(VertexId id, Bitset subsplit) : id_{id}, subsplit_{std::move(subsplit)} {}

  VertexId GetId() const { return id_; }
  const Bitset& GetSubsplit() const { return subsplit_; }

  NeighborsView GetNeighbors(Direction direction, Clade clade) {
    return neighbors_.at({direction, clade});
  }

  ConstNeighborsView GetNeighbors(Direction direction, Clade clade) const {
    return neighbors_.at({direction, clade});
  }

  std::optional<std::tuple<LineId, Direction, Clade>> FindNeighbor(
      VertexId neighbor) const {
    for (auto i : neighbors_) {
      auto j = i.second.find(neighbor);
      if (j != i.second.end()) return {{j->second, i.first.first, i.first.second}};
    }
    return {};
  }

  DAGVertex& SetId(VertexId id) {
    id_ = id;
    return *this;
  }
  DAGVertex& SetSubsplit(Bitset subsplit) {
    subsplit_ = std::move(subsplit);
    return *this;
  }
  DAGVertex& AddNeighbor(Direction direction, Clade clade, VertexId neighbor,
                         LineId line) {
    neighbors_.at({direction, clade}).insert({neighbor, line});
    return *this;
  }

  void SetLineId(VertexId neighbor, LineId line) {
    for (auto& i : neighbors_) {
      auto j = i.second.find(neighbor);
      if (j != i.second.end()) {
        i.second.insert_or_assign(j, neighbor, line);
        return;
      }
    }
    Failwith("Neighbor not found");
  }

 private:
  VertexId id_ = NoId;
  Bitset subsplit_ = Bitset{{}};
  std::map<std::pair<Direction, Clade>, std::map<VertexId, LineId>> neighbors_ = {
      {{Direction::Rootward, Clade::Left}, {}},
      {{Direction::Rootward, Clade::Right}, {}},
      {{Direction::Leafward, Clade::Left}, {}},
      {{Direction::Leafward, Clade::Right}, {}},
  };
};

template <typename T>
class GenericLinesView {
  using view_type = std::conditional_t<std::is_const_v<T>, ConstLineView, LineView>;

 public:
  class Iterator : public std::iterator<std::forward_iterator_tag, view_type> {
   public:
    Iterator(const GenericLinesView& view, size_t index) : view_{view}, index_{index} {}

    Iterator operator++() {
      ++index_;
      return *this;
    }
    bool operator!=(const Iterator& other) const { return index_ != other.index_; }
    view_type operator*() const { return view_.storage_.lines_[index_]; }

   private:
    const GenericLinesView& view_;
    size_t index_;
  };

  GenericLinesView(T& storage) : storage_{storage} {}

  Iterator begin() const { return {*this, 0}; }
  Iterator end() const { return {*this, storage_.lines_.size()}; }
  size_t size() const { return storage_.lines_.size(); };
  view_type operator[](size_t i) const { return storage_.lines_[i]; };

 private:
  friend class Iterator;
  T& storage_;
};
class SubsplitDAGStorage;
using LinesView = GenericLinesView<SubsplitDAGStorage>;
using ConstLinesView = GenericLinesView<const SubsplitDAGStorage>;

template <typename T>
class GenericVerticesView {
  using view_type =
      std::conditional_t<std::is_const_v<T>, SubsplitDAGNode, MutableSubsplitDAGNode>;

 public:
  class Iterator : public std::iterator<std::forward_iterator_tag, view_type> {
   public:
    Iterator(T& storage, size_t index) : storage_{storage}, index_{index} {}

    Iterator& operator++() {
      ++index_;
      return *this;
    }
    Iterator operator++(int) { return {storage_, index_++}; }
    Iterator operator+(size_t i) { return {storage_, index_ + i}; }
    Iterator operator-(size_t i) { return {storage_, index_ - i}; }
    bool operator<(const Iterator& other) { return index_ < other.index_; }
    bool operator!=(const Iterator& other) const { return index_ != other.index_; }
    view_type operator*() const;

   private:
    T& storage_;
    size_t index_;
  };

  GenericVerticesView(T& storage) : storage_{storage} {}

  Iterator begin() const { return {storage_, 0}; }
  Iterator end() const { return {storage_, storage_.vertices_.size()}; }
  Iterator cbegin() const { return {storage_, 0}; }
  Iterator cend() const { return {storage_, storage_.vertices_.size()}; }
  size_t size() const { return storage_.vertices_.size(); };
  view_type operator[](size_t i) const;
  view_type at(size_t i) const;

 private:
  friend class Iterator;
  T& storage_;
};
using VerticesView = GenericVerticesView<SubsplitDAGStorage>;
using ConstVerticesView = GenericVerticesView<const SubsplitDAGStorage>;

class SubsplitDAGStorage {
 public:
  ConstVerticesView GetVertices() const { return *this; }
  VerticesView GetVertices() { return *this; }

  ConstLinesView GetLines() const { return *this; }
  LinesView GetLines() { return *this; }

  std::optional<ConstLineView> GetLine(VertexId parent, VertexId child, size_t vertex_offset = 0, size_t line_offset = 0) const {
    parent -= vertex_offset;
    child -= vertex_offset;
    if (parent >= vertices_.size() || child >= vertices_.size()) return {};
    auto result = vertices_.at(parent).FindNeighbor(child + vertex_offset);
    if (!result.has_value()) return {};
    return lines_.at(std::get<0>(result.value()) - line_offset);
  }

  std::optional<ConstLineView> GetLine(LineId id) const {
    if (id >= lines_.size()) return {};
    auto& line = lines_[id];
    if (line.GetId() == NoId) return {};
    return line;
  }

  bool ContainsVertex(VertexId id) const {
    if (id >= vertices_.size()) return false;
    return vertices_[id].GetId() != NoId;
  }

  DAGLineStorage& AddLine(const DAGLineStorage& newLine, size_t offset = 0) {
    if (newLine.GetId() == NoId) Failwith("Set line id before inserting to storage");
    if (newLine.GetClade() == Clade::Unspecified)
      Failwith("Set clade before inserting to storage");
    auto& line = GetOrInsert(lines_, newLine.GetId(), offset);
    line = newLine;
    return line;
  }

  DAGVertex& AddVertex(const DAGVertex& newVertex, size_t offset = 0) {
    if (newVertex.GetId() == NoId)
      Failwith("Set vertex id before inserting to storage");
    auto& vertex = GetOrInsert(vertices_, newVertex.GetId(), offset);
    vertex = newVertex;
    return vertex;
  }

  void ReindexLine(LineId line, VertexId parent, VertexId child) {
    lines_.at(line).SetParent(parent).SetChild(child);
  }

  void SetLines(const std::vector<DAGLineStorage>& lines) { lines_ = lines; }

  void SetVertices(const std::vector<DAGVertex>& vertices) { vertices_ = vertices; }

 private:
  template <typename>
  friend class GenericLinesView;
  template <typename>
  friend class GenericVerticesView;

  template <typename T, typename Id>
  static T& GetOrInsert(std::vector<T>& data, Id id, size_t offset = 0) {
    id -= offset;
    if (id >= data.size()) data.resize(id + 1);
    return data[id];
  }

  std::vector<DAGLineStorage> lines_;
  std::vector<DAGVertex> vertices_;
};
