// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <vector>
#include <map>
#include <optional>
#include <functional>

#include "basic_graph.hpp"

enum class Direction { Rootward, Leafward };
enum class Clade { Left, Right };

template <typename Traits>
class DAGLine {
public:

  using Line = typename Traits::Line;

  LineId GetId() const { return id_; }
  Clade GetClade() const { return clade_; }
  
  auto& SetClade(Clade clade) {
    clade_ = clade;
    return static_cast<typename Traits::Line&>(*this);
  }
  
  std::pair<VertexId, VertexId> GetVertexIds() const {
    return {GetLine().GetFirst(), GetLine().GetSecond()};
  }

  constexpr VertexId GetParent() const { return GetLine().GetFirst(); }
  constexpr VertexId GetChild() const { return GetLine().GetSecond(); }

private:

  template <typename, template <typename> typename ... > friend class BasicLine;

  constexpr const Line& GetLine() const {
    return static_cast<const Line&>(*this);
  }

  template <typename Line>
  void OnLineChanged(const Line&, LineId id) {
    id_ = id;
  }

  LineId id_;
  Clade clade_;
};

template <typename Traits>
class DAGVertex {
public:

  VertexId GetId() const { return id_; }

private:

  template <typename, template <typename> typename ... > friend class BasicVertex;

  template <typename Line>
  void OnLineChanged(VertexId id, const Line&, LineId) {
    id_ = id;
  }

  VertexId id_;
};

template <typename Traits>
class DAGNeighbors : public Traits {
public:

  using typename Traits::Line;
  using typename Traits::Vertex;

  const std::map<VertexId, LineId>& Get(Direction direction, Clade clade) const {
    return data_[{direction, clade}];
  }
  
  std::tuple<LineId, Direction, Clade> Find(VertexId id) const {
    for (auto i : data_) {
      auto j = i.second.find(id);
      if (j != i.second.end()) return {j->second, i.first.first, i.first.second};
    }
    return {NoId, {}, {}};
  }

private:
  friend Vertex;
  
  void OnLineChanged(const Vertex&, VertexId vertexId, const Line& line, LineId lineId) {
    Direction direction;
    if (vertexId == line.GetFirst()) direction = Direction::Leafward;
    else if (vertexId == line.GetSecond()) direction = Direction::Rootward;
    else throw std::out_of_range("Vertex doesn't belong to the line");
    data_[{direction, line.GetClade()}].insert({direction == Direction::Rootward ? line.GetFirst() : line.GetSecond(), lineId});
  }
  
  mutable std::map<std::pair<Direction, Clade>, std::map<VertexId, LineId>> data_;
};

template <typename T, typename Id>
T& GetOrInsert(std::vector<T>& data, Id id) {
  if (id >= data.size()) data.resize(id + 1);
  return data[id];
}

template <typename Traits>
class DAGLookup {
public:
  using Line = typename Traits::Line;
  using Graph = typename Traits::Graph;
  
  std::optional<std::reference_wrapper<const Line>> GetLine(VertexId parent, VertexId child) const {
    auto id = std::get<0>(GetGraph().GetVertices()[parent].GetNeighbors().Find(child));
    if (id == NoId) return {};
    return GetGraph().GetLines()[id];
  }
  
  std::optional<std::reference_wrapper<const Line>> GetLine(LineId id) const {
    if (id >= GetGraph().GetLines().size()) return {};
    auto& line = GetGraph().GetLines()[id];
    if (line.GetId() == NoId) return {};
    return line;
  }
  
private:
  constexpr const Graph& GetGraph() const {
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
