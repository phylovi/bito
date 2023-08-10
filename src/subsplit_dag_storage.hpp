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
// collection - internally a std::map<NodeId, EdgeId>.
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
#include "reindexer.hpp"
#include "sugar.hpp"

using NodeId = GenericId<struct NodeIdTag>;
using EdgeId = GenericId<struct EdgeIdTag>;
using TaxonId = GenericId<struct TaxonIdTag>;

using StringTaxonIdMap = std::unordered_map<std::string, TaxonId>;
using BitsetNodeIdMap = std::unordered_map<Bitset, NodeId>;
using BitsetNodeIdVectorMap = std::unordered_map<Bitset, std::vector<NodeId>>;
using NodeIdBitsetMap = std::unordered_map<NodeId, Bitset>;
using EdgeIdPair = std::pair<EdgeId, EdgeId>;
using NodeIdPair = std::pair<NodeId, NodeId>;
using NodeIdEdgeIdPairMap = std::unordered_map<NodeId, EdgeIdPair>;
using NodeIdVector = std::vector<NodeId>;
using EdgeIdVector = std::vector<EdgeId>;
using TaxonIdVector = std::vector<TaxonId>;
using BitsetEdgeIdMap = std::unordered_map<Bitset, EdgeId>;
using EdgeIdBitsetMap = std::unordered_map<EdgeId, Bitset>;
using BitsetEdgeIdPairMap = std::unordered_map<Bitset, EdgeIdPair>;
using NodeIdVectorPair = std::pair<NodeIdVector, NodeIdVector>;

enum class Direction { Rootward, Leafward };
static const inline size_t DirectionCount = 2;
class DirectionEnum : public EnumWrapper<Direction, size_t, DirectionCount,
                                         Direction::Rootward, Direction::Leafward> {
 public:
  static inline const std::string Prefix = "Direction";
  static inline const Array<std::string> Labels = {{"Rootward", "Leafward"}};

  static std::string ToString(const Direction e) {
    std::stringstream ss;
    ss << Prefix << "::" << Labels[e];
    return ss.str();
  }
  friend std::ostream& operator<<(std::ostream& os, const Direction e) {
    os << ToString(e);
    return os;
  }
};

template <typename Derived>
class DAGLine {
 public:
  EdgeId GetId() const { return storage().id_; }
  NodeId GetParent() const { return storage().parent_; }
  NodeId GetChild() const { return storage().child_; }
  SubsplitClade GetSubsplitClade() const { return storage().clade_; }

  Derived& SetId(EdgeId id) {
    storage().id_ = id;
    return derived();
  }
  Derived& SetParent(NodeId id) {
    storage().parent_ = id;
    return derived();
  }
  Derived& SetChild(NodeId id) {
    storage().child_ = id;
    return derived();
  }
  Derived& SetSubsplitClade(SubsplitClade clade) {
    storage().clade_ = clade;
    return derived();
  }

  std::pair<NodeId, NodeId> GetNodeIds() const { return {GetParent(), GetChild()}; }

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
  DAGLineStorage(EdgeId id, NodeId parent, NodeId child, SubsplitClade clade)
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

  EdgeId id_ = EdgeId(NoId);
  NodeId parent_ = NodeId(NoId);
  NodeId child_ = NodeId(NoId);
  SubsplitClade clade_ = SubsplitClade::Unspecified;
};

template <typename T>
class DAGLineView : public DAGLine<DAGLineView<T>> {
 public:
  DAGLineView(T& line) : line_{line} {
    Assert(line.GetSubsplitClade() != SubsplitClade::Unspecified, "Uninitialized edge");
  }

  template <size_t I>
  auto get() const {
    static_assert(I < 2, "Index out of bounds");
    if constexpr (I == 0) return line_.GetNodeIds();
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
  using type = pair<::NodeId, ::NodeId>;
};
template <>
struct tuple_element<1, ::LineView> {
  using type = ::EdgeId;
};

template <>
struct tuple_size<::ConstLineView> : public integral_constant<size_t, 2> {};
template <>
struct tuple_element<0, ::ConstLineView> {
  using type = const pair<::NodeId, ::NodeId>;
};
template <>
struct tuple_element<1, ::ConstLineView> {
  using type = const ::EdgeId;
};
}  // namespace std

template <typename T>
class GenericNeighborsView {
  using iterator_type =
      std::conditional_t<std::is_const_v<T>, std::map<NodeId, EdgeId>::const_iterator,
                         std::map<NodeId, EdgeId>::iterator>;
  using map_type =
      std::conditional_t<std::is_const_v<T>, const std::map<NodeId, EdgeId>,
                         std::map<NodeId, EdgeId>>;

 public:
  class Iterator : public std::iterator<std::forward_iterator_tag, NodeId> {
   public:
    Iterator(const iterator_type& i, map_type& map) : iter_{i}, map_{map} {}

    Iterator operator++() {
      ++iter_;
      return *this;
    }
    bool operator!=(const Iterator& other) const { return iter_ != other.iter_; }
    bool operator==(const Iterator& other) const { return iter_ == other.iter_; }

    NodeId operator*() { return iter_->first; }

    NodeId GetNodeId() const { return iter_->first; }
    EdgeId GetEdge() const { return iter_->second; }

   private:
    iterator_type iter_;
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

  void RemapNodeIds(const Reindexer& reindexer) {
    for (auto [vertex_id, line_id] : neighbors_) {
      std::ignore = line_id;
      Assert(vertex_id < reindexer.size(),
             "Neighbors cannot contain an id out of bounds of the reindexer in "
             "GenericNeighborsView::RemapIds.");
    }
    T remapped{};
    for (auto [vertex_id, line_id] : neighbors_) {
      remapped[NodeId(reindexer.GetNewIndexByOldIndex(size_t(vertex_id)))] = line_id;
    }
    neighbors_ = remapped;
  }

  void RemapEdgeIdxs(const Reindexer& reindexer) {
    for (auto [vertex_id, line_id] : neighbors_) {
      std::ignore = vertex_id;
      Assert(line_id < reindexer.size(),
             "Neighbors cannot contain an id out of bounds of the reindexer in "
             "GenericNeighborsView::RemapIds.");
    }
    T remapped{};
    for (auto [vertex_id, line_id] : neighbors_) {
      remapped[vertex_id] = reindexer.GetNewIndexByOldIndex(size_t(line_id));
    }
    neighbors_ = remapped;
  }

  operator SizeVector() const {
    SizeVector result;
    for (auto i : neighbors_) {
      result.push_back(size_t(i.first));
    }
    return result;
  }

  operator NodeIdVector() const {
    NodeIdVector result;
    for (auto i : neighbors_) {
      result.push_back(i.first);
    }
    return result;
  }

  void SetNeighbors(const T& neighbors) { neighbors_ = neighbors; }

 private:
  friend class Iterator;
  template <typename>
  friend class ::GenericNeighborsView;
  T& neighbors_;
};
using NeighborsView = GenericNeighborsView<std::map<NodeId, EdgeId>>;
using ConstNeighborsView = GenericNeighborsView<const std::map<NodeId, EdgeId>>;

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
  DAGVertex(NodeId id, Bitset subsplit) : id_{id}, subsplit_{std::move(subsplit)} {}

  NodeId GetId() const { return id_; }
  const Bitset& GetSubsplit() const { return subsplit_; }

  NeighborsView GetNeighbors(Direction direction, SubsplitClade clade) {
    return neighbors_.at({direction, clade});
  }

  ConstNeighborsView GetNeighbors(Direction direction, SubsplitClade clade) const {
    return neighbors_.at({direction, clade});
  }

  bool IsRoot() const {
    return GetNeighbors(Direction::Rootward, SubsplitClade::Left).empty() &&
           GetNeighbors(Direction::Rootward, SubsplitClade::Right).empty();
  }

  bool IsLeaf() const {
    return GetNeighbors(Direction::Leafward, SubsplitClade::Left).empty() &&
           GetNeighbors(Direction::Leafward, SubsplitClade::Right).empty();
  }

  std::optional<std::tuple<EdgeId, Direction, SubsplitClade>> FindNeighbor(
      NodeId neighbor) const {
    for (auto i : neighbors_) {
      auto j = i.second.find(neighbor);
      if (j != i.second.end()) return {{j->second, i.first.first, i.first.second}};
    }
    return {};
  }

  DAGVertex& SetId(NodeId id) {
    id_ = id;
    return *this;
  }
  DAGVertex& SetSubsplit(Bitset subsplit) {
    subsplit_ = std::move(subsplit);
    return *this;
  }
  DAGVertex& AddNeighbor(Direction direction, SubsplitClade clade, NodeId neighbor,
                         EdgeId line) {
    neighbors_.at({direction, clade}).insert({neighbor, line});
    return *this;
  }
  void RemoveNeighbor(Direction direction, SubsplitClade clade, NodeId neighbor) {
    neighbors_.at({direction, clade}).erase(neighbor);
  }

  void SetEdgeId(NodeId neighbor, EdgeId line) {
    for (auto& i : neighbors_) {
      auto j = i.second.find(neighbor);
      if (j != i.second.end()) {
        i.second.insert_or_assign(j, neighbor, line);
        return;
      }
    }
    Failwith("Neighbor not found");
  }

  void ClearNeighbors() {
    neighbors_ = {
        {{Direction::Rootward, SubsplitClade::Left}, {}},
        {{Direction::Rootward, SubsplitClade::Right}, {}},
        {{Direction::Leafward, SubsplitClade::Left}, {}},
        {{Direction::Leafward, SubsplitClade::Right}, {}},
    };
  }

 private:
  NodeId id_ = NodeId(NoId);
  Bitset subsplit_ = Bitset{{}};
  std::map<std::pair<Direction, SubsplitClade>, std::map<NodeId, EdgeId>> neighbors_ = {
      {{Direction::Rootward, SubsplitClade::Left}, {}},
      {{Direction::Rootward, SubsplitClade::Right}, {}},
      {{Direction::Leafward, SubsplitClade::Left}, {}},
      {{Direction::Leafward, SubsplitClade::Right}, {}},
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

// A vector that can optionally be prepended with host data for Graft purposes.
// It has no ownership over the host data, so the lifetime of the host data object
// should be greater than the given HostableVector instance.
template <typename T>
class HostableVector {
 public:
  explicit HostableVector(HostableVector<T>* host = nullptr)
      : host_{host ? &host->data_ : nullptr} {}

  T& at(size_t i) {
    if (!host_) {
      return data_.at(i);
    }
    if (i < host_->size()) {
      return (*host_)[i];
    }
    return data_.at(i - host_->size());
  }

  const T& at(size_t i) const {
    if (!host_) {
      return data_.at(i);
    }
    if (i < host_->size()) {
      return (*host_)[i];
    }
    return data_.at(i - host_->size());
  }

  T& operator[](size_t i) {
    if (!host_) {
      return data_[i];
    }
    if (i < host_->size()) {
      return (*host_)[i];
    }
    return data_[i - host_->size()];
  }

  const T& operator[](size_t i) const {
    if (!host_) {
      return data_[i];
    }
    if (i < host_->size()) {
      return (*host_)[i];
    }
    return data_[i - host_->size()];
  }

  size_t size() const {
    if (!host_) {
      return data_.size();
    }
    return host_->size() + data_.size();
  }

  void resize(size_t new_size) {
    if (!host_) {
      data_.resize(new_size);
    } else {
      data_.resize(new_size - host_->size());
    }
  }

  HostableVector& operator=(const std::vector<T>& data) {
    data_ = data;
    return *this;
  }

  bool HaveHost() const { return host_; }

  size_t HostSize() const {
    if (!host_) {
      return 0;
    }
    return host_->size();
  }

  void ResetHost(HostableVector<T>* host) {
    host_ = host ? &host->data_ : nullptr;
    data_ = {};
  }

 private:
  std::vector<T> data_;
  std::vector<T>* host_;
};

// Tag dispatching type to avoid confusion with copy constructor.
// See https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Tag_Dispatching
struct HostDispatchTag {};

class SubsplitDAGStorage {
 public:
  SubsplitDAGStorage() = default;
  SubsplitDAGStorage(const SubsplitDAGStorage&) = default;
  SubsplitDAGStorage(SubsplitDAGStorage& host, HostDispatchTag)
      : lines_{&host.lines_}, vertices_{&host.vertices_} {}

  ConstVerticesView GetVertices() const { return *this; }
  VerticesView GetVertices() { return *this; }

  ConstLinesView GetLines() const { return *this; }
  LinesView GetLines() { return *this; }

  std::optional<ConstLineView> GetLine(NodeId parent, NodeId child) const {
    if (parent.value_ >= vertices_.size() || child.value_ >= vertices_.size()) {
      return {};
    }
    auto result = vertices_.at(parent.value_).FindNeighbor(child);
    if (!result.has_value()) {
      return {};
    }
    // if (std::get<0>(result.value()).value_ >= lines_.size()) {
    //   std::cerr << "ERROR: GetLine -- parent-child vertice pair reference value "
    //                "outside of line range."
    //             << std::endl;
    //   return {};
    // }
    return lines_.at(std::get<0>(result.value()).value_);
  }

  std::optional<ConstLineView> GetLine(EdgeId id) const {
    if (id.value_ >= lines_.size()) {
      return std::nullopt;
    }
    auto& line = lines_[id.value_];
    if (line.GetId() == EdgeId(NoId)) {
      return std::nullopt;
    }
    return line;
  }

  const DAGVertex& GetVertex(NodeId id) const { return vertices_.at(id.value_); }

  bool ContainsVertex(NodeId id) const {
    if (id.value_ >= vertices_.size()) return false;
    return vertices_[id.value_].GetId() != NodeId(NoId);
  }

  std::optional<std::reference_wrapper<DAGVertex>> FindVertex(const Bitset& subsplit) {
    for (NodeId id = 0; id < vertices_.size(); ++id.value_) {
      if (vertices_[id.value_].GetSubsplit() == subsplit) {
        return vertices_[id.value_];
      }
    }
    return std::nullopt;
  }

  std::optional<std::reference_wrapper<const DAGVertex>> FindVertex(
      const Bitset& subsplit) const {
    for (NodeId id = 0; id.value_ < vertices_.size(); ++id.value_) {
      if (vertices_[id.value_].GetSubsplit() == subsplit) {
        return vertices_[id.value_];
      }
    }
    return std::nullopt;
  }

  DAGLineStorage& AddLine(const DAGLineStorage& newLine) {
    if (newLine.GetId() == EdgeId(NoId))
      Failwith("Set line id before inserting to storage");
    if (newLine.GetSubsplitClade() == SubsplitClade::Unspecified)
      Failwith("Set clade before inserting to storage");
    auto& line = GetOrInsert(lines_, newLine.GetId());
    line = newLine;
    return line;
  }

  DAGVertex& AddVertex(const DAGVertex& newVertex) {
    if (newVertex.GetId() == NodeId(NoId))
      Failwith("Set vertex id before inserting to storage");
    auto& vertex = GetOrInsert(vertices_, newVertex.GetId());
    vertex = newVertex;
    return vertex;
  }

  void ReindexLine(EdgeId line, NodeId parent, NodeId child) {
    lines_.at(line.value_).SetParent(parent).SetChild(child);
  }

  void SetLines(const std::vector<DAGLineStorage>& lines) { lines_ = lines; }

  void SetVertices(const std::vector<DAGVertex>& vertices) { vertices_ = vertices; }

  bool HaveHost() const { return lines_.HaveHost(); }

  size_t HostLinesCount() const { return lines_.HostSize(); }

  size_t HostVerticesCount() const { return vertices_.HostSize(); }

  void ResetHost(SubsplitDAGStorage& host) {
    for (NodeId id = lines_.HostSize(); id < lines_.size(); ++id) {
      auto& line = lines_[id.value_];
      if (line.GetParent().value_ >= vertices_.HostSize()) {
        if (line.GetChild().value_ < vertices_.HostSize()) {
          vertices_[line.GetChild().value_].RemoveNeighbor(
              Direction::Rootward, line.GetSubsplitClade(), line.GetParent());
        }
      }
      if (line.GetChild().value_ >= vertices_.HostSize()) {
        if (line.GetParent().value_ < vertices_.HostSize()) {
          vertices_[line.GetParent().value_].RemoveNeighbor(
              Direction::Leafward, line.GetSubsplitClade(), line.GetChild());
        }
      }
    }
    lines_.ResetHost(&host.lines_);
    vertices_.ResetHost(&host.vertices_);
  }

  void ConnectVertices(EdgeId line_id) {
    auto& line = lines_.at(line_id.value_);
    auto& parent = vertices_.at(line.GetParent().value_);
    auto& child = vertices_.at(line.GetChild().value_);
    parent.AddNeighbor(Direction::Leafward, line.GetSubsplitClade(), child.GetId(),
                       line_id);
    child.AddNeighbor(Direction::Rootward, line.GetSubsplitClade(), parent.GetId(),
                      line_id);
  }

  void ConnectAllVertices() {
    for (size_t i = 0; i < vertices_.size(); ++i) {
      vertices_[i].ClearNeighbors();
    }
    for (size_t i = 0; i < lines_.size(); ++i) {
      if (lines_[i].GetId() == EdgeId(NoId)) {
        continue;
      }
      ConnectVertices(lines_[i].GetId());
    }
  }

  std::optional<std::reference_wrapper<const DAGVertex>> FindRoot() const {
    for (size_t i = 0; i < vertices_.size(); ++i) {
      auto& vertex = vertices_[i];
      if (vertex.GetId() == NoId) {
        continue;
      }
      if (vertex.GetNeighbors(Direction::Rootward, SubsplitClade::Left).empty() &&
          vertex.GetNeighbors(Direction::Rootward, SubsplitClade::Right).empty()) {
        return vertex;
      }
    }
    return std::nullopt;
  }

  void Print() const {
    std::cout << "== LINES ==" << std::endl;
    for (size_t id = 0; id < lines_.size(); id++) {
      auto& line = lines_[id];
      std::cout << "[" << id << "]: { id::" << line.GetId()
                << " parent::" << line.GetParent() << " child::" << line.GetChild()
                << " }" << std::endl;
    }
    std::cout << "== VERTICES ==" << std::endl;
    for (size_t id = 0; id < vertices_.size(); id++) {
      auto& node = vertices_[id];
      std::cout << "[" << id << " " << node.GetId() << "]: { ";
      for (const auto dir : DirectionEnum::Iterator()) {
        for (const auto clade : SubsplitCladeEnum::Iterator()) {
          std::cout << "[ ";
          for (const auto adj_node_id : node.GetNeighbors(dir, clade)) {
            std::cout << "adj::" << adj_node_id << " ";
            auto parent = (dir == Direction::Leafward) ? node.GetId() : adj_node_id;
            auto child = (dir == Direction::Rootward) ? adj_node_id : node.GetId();
            auto result = vertices_.at(parent.value_).FindNeighbor(child);
            if (result.has_value()) {
              auto line_id = std::get<0>(result.value()).value_;
              std::cout << "line::" << line_id << " ";
            }
          }
          std::cout << " ]";
        }
      }
      std::cout << "} " << std::endl;
    }
  }

 private:
  template <typename>
  friend class GenericLinesView;
  template <typename>
  friend class GenericVerticesView;

  template <typename T, typename Id>
  static T& GetOrInsert(HostableVector<T>& data, Id id) {
    if (id.value_ >= data.size()) {
      data.resize(id.value_ + 1);
    }
    return data[id.value_];
  }

  HostableVector<DAGLineStorage> lines_;
  HostableVector<DAGVertex> vertices_;
};
