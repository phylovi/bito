// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <vector>
#include <map>
#include <optional>
#include <functional>

#include "basic_graph.hpp"
#include "bitset.hpp"

enum class Direction { Rootward, Leafward };
enum class Clade { Left, Right };

template <typename T>
using OptRef = std::optional<std::reference_wrapper<T>>;

template <typename Traits>
class DAGLine {
  using Line = typename Traits::Line;
  
public:
  constexpr VertexId GetParent() const { return AsLine().GetFirst(); }
  constexpr VertexId GetChild() const { return AsLine().GetSecond(); }

  LineId GetId() const { return id_; }
  Clade GetClade() const { return clade_; }
  
  auto& SetClade(Clade clade) {
    clade_ = clade;
    return AsLine();
  }
  
  std::pair<VertexId, VertexId> GetVertexIds() const {
    return {GetParent(), GetChild()};
  }

protected:
  void OnLineChanged(const Line&, LineId id) {
    id_ = id;
  }

private:
  constexpr const Line& AsLine() const {
    return static_cast<const Line&>(*this);
  }

  LineId id_;
  Clade clade_;
};

template <typename Traits>
class DAGVertex {
  using Line = typename Traits::Line;
  using Vertex = typename Traits::Vertex;
  
public:
  VertexId GetId() const { return id_; }
  Vertex& SetId(VertexId id) { id_ = id; return AsVertex(); }
  const Bitset& GetSubsplit() const { return subsplit_; };
  Vertex& SetSubsplit(const Bitset& subsplit) { subsplit_ = subsplit; return AsVertex(); }

protected:
  void OnLineChanged(VertexId id, const Line&, LineId) {
    id_ = id;
  }

private:
  constexpr Vertex& AsVertex() {
    return static_cast<Vertex&>(*this);
  }

  VertexId id_ = NoId;
  Bitset subsplit_ = Bitset{0};
};

template <typename Traits>
class DAGNeighbors : public Traits {
  using typename Traits::Line;
  using typename Traits::Vertex;
  friend Vertex;
  
public:
  DAGNeighbors() {
    data_.insert({{Direction::Rootward, Clade::Left}, {}});
    data_.insert({{Direction::Rootward, Clade::Right}, {}});
    data_.insert({{Direction::Leafward, Clade::Left}, {}});
    data_.insert({{Direction::Leafward, Clade::Right}, {}});
  }

  const std::map<VertexId, LineId>& Get(Direction direction, Clade clade) const {
    return data_.at({direction, clade});
  }
  
  std::tuple<LineId, Direction, Clade> Find(VertexId id) const {
    for (auto i : data_) {
      auto j = i.second.find(id);
      if (j != i.second.end()) return {j->second, i.first.first, i.first.second};
    }
    return {NoId, {}, {}};
  }

private:
  
  void OnLineChanged(const Vertex&, VertexId vertexId, const Line& line, LineId lineId) {
    Direction direction;
    if (vertexId == line.GetParent()) direction = Direction::Leafward;
    else if (vertexId == line.GetChild()) direction = Direction::Rootward;
    else throw std::out_of_range("Vertex doesn't belong to the line");
    data_.at({direction, line.GetClade()}).insert({direction == Direction::Rootward ? line.GetParent() : line.GetChild(), lineId});
  }
  
  std::map<std::pair<Direction, Clade>, std::map<VertexId, LineId>> data_;
};

template <typename Traits>
class DAGLookup {
  using Line = typename Traits::Line;
  using Graph = typename Traits::Graph;
  
public:
  OptRef<const Line> GetLine(VertexId parent, VertexId child) const {
    auto id = std::get<0>(AsGraph().GetVertices()[parent].GetNeighbors().Find(child));
    if (id == NoId) return {};
    return AsGraph().GetLines()[id];
  }
  
  OptRef<const Line> GetLine(LineId id) const {
    if (id >= AsGraph().GetLines().size()) return {};
    auto& line = AsGraph().GetLines()[id];
    if (line.GetId() == NoId) return {};
    return line;
  }

  bool ContainsVertex(VertexId id) const {
      if (id >= AsGraph().GetVertices().size()) return false;
      return AsGraph().GetVertices()[id].GetId() != NoId;
  }
  
private:
  constexpr const Graph& AsGraph() const {
    return static_cast<const Graph&>(*this);
  }
};

struct DAGTraits {
  using Line = BasicLine<DAGTraits, DAGLine>;
  using Vertex = BasicVertex<DAGTraits, DAGVertex>;
  using Neighbors = DAGNeighbors<DAGTraits>;
  using Lines = std::vector<Line>;
  using Vertices = std::vector<Vertex>;
  using Graph = BasicGraph<DAGTraits, DAGLookup>;
};

using DAGStorage = DAGTraits::Graph;

template <typename T, typename Id>
T& GetOrInsert(std::vector<T>& data, Id id) {
  if (id >= data.size()) data.resize(id + 1);
  return data[id];
}
